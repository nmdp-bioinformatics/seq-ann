# -*- coding: utf-8 -*-

#
#    seqann Sequence Annotation
#    Copyright (c) 2017 Be The Match operated by National Marrow Donor Program. All Rights Reserved.
#
#    This library is free software; you can redistribute it and/or modify it
#    under the terms of the GNU Lesser General Public License as published
#    by the Free Software Foundation; either version 3 of the License, or (at
#    your option) any later version.
#
#    This library is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; with out even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
#    License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License
#    along with this library;  if not, write to the Free Software Foundation,
#    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA.
#
#    > http://www.fsf.org/licensing/licenses/lgpl.html
#    > http://www.opensource.org/licenses/lgpl-license.php
#

from __future__ import absolute_import

import re
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import ExactPosition
from Bio.Align.Applications import ClustalOmegaCommandline

from seqann.util import cleanup
from seqann.util import randomid
from seqann.util import get_seqfeat
from seqann.seq_search import getblocks
from seqann.models.annotation import Annotation


def align_seqs(found_seqs, sequence, locus, verbose=False):

    randid = randomid()
    input_fasta = str(randid) + ".fasta"
    output_clu = str(randid) + ".clu"

    seqs = []
    seqs.append(found_seqs)
    seqs.append(sequence)
    SeqIO.write(seqs, input_fasta, "fasta")
    clustalomega_cline = ClustalOmegaCommandline(infile=input_fasta,
                                                 outfile=output_clu,
                                                 outfmt='clu', wrap=50000,
                                                 verbose=True, auto=True)
    stdout, stderr = clustalomega_cline()
    align = AlignIO.read(output_clu, "clustal")

    cleanup(randid)

    insers, dels = 0, 0
    all_features = []
    if len(align)-2 == 0:
        infeats = get_seqfeat(seqs[0])
        diffs = count_diffs(align, infeats, sequence, verbose)
        if isinstance(diffs, Annotation):
            return diffs, 0, 0
        else:
            insers, dels = diffs[0], diffs[1]
        f = find_features(infeats, align[0])
        all_features.append(f)
    else:
        for i in range(0, len(align)-2):
            infeats = get_seqfeat(seqs[i])
            f = find_features(infeats, align[i])
            all_features.append(f)

    annotation = resolve_feats(all_features, align[len(align)-1], verbose)
    return annotation, insers, dels


def find_features(feats, sequ):
    feats_a = list(feats.keys())

    j = 0
    s = 0
    en = 0
    start = 0
    for i in range(0, len(sequ)):
        if j <= len(feats_a)-1:
            if i > int(feats[feats_a[j]].location.end):
                j += 1
            if(sequ[i] == '-'):
                if i == 0:
                    start += 1
                    en += 1
                    s = 1
                else:
                    start += 1
                    en += 1
                    if s == 0:
                        feats[feats_a[j]] = SeqFeature(FeatureLocation(ExactPosition(feats[feats_a[j]].location.start), ExactPosition(int(feats[feats_a[j]].location.end + 1)), strand=1), type=feats[feats_a[j]].type)
                        if j != len(feats_a):
                            for l in range(j+1, len(feats_a)):
                                feats[feats_a[l]] = SeqFeature(FeatureLocation(ExactPosition(feats[feats_a[l]].location.start+1), ExactPosition(int(feats[feats_a[l]].location.end + 1)), strand=1), type=feats[feats_a[l]].type)
                      
            else:
                if s == 1:
                    st = feats[feats_a[j]].location.start + start
                    end = feats[feats_a[j]].location.end + en
                    feats[feats_a[j]] = SeqFeature(FeatureLocation(ExactPosition(st), ExactPosition(end), strand=1), type=feats[feats_a[j]].type)
                    if j != len(feats_a):
                        for l in range(j+1, len(feats_a)):
                            feats[feats_a[l]] = SeqFeature(FeatureLocation(ExactPosition(feats[feats_a[l]].location.start+st+1), ExactPosition(int(feats[feats_a[l]].location.end + end)), strand=1), type=feats[feats_a[l]].type)

                    s = 0
    return feats


def resolve_feats(feat_list, seq, verbose):

    # Look at consensus sequence
    #   gap_consensus
    seq_covered = len(seq.seq)
    coordinates = dict(map(lambda x: [x, 1],
                       [i for i in range(0, len(seq.seq)+1)]))

    mapping = dict(map(lambda x: [x, 1],
                       [i for i in range(0, len(seq.seq)+1)]))
    if len(feat_list) > 1:
        if verbose:
            print("****** resolve_feats error ******", file=sys.stderr)
        for i in range(0, len(feat_list)):
            for j in range(0, len(feat_list)):
                if i != j:
                    print(j, i)
    else:
        full_annotation = {}
        features = feat_list[0]
        for feat in features:
            f = features[feat]
            seqrec = f.extract(seq)
            seq_covered -= len(seqrec.seq)
            if re.search("-", str(seqrec.seq)):
                newseq = re.sub(r'-', '', str(seqrec.seq))
                seqrec.seq = Seq(newseq, IUPAC.unambiguous_dna)
            if seqrec.seq:
                full_annotation.update({feat: seqrec})
                for i in range(f.location.start, f.location.end):
                    if i in coordinates:
                        del coordinates[i]
                    mapping[i] = feat

        blocks = getblocks(coordinates)
        annotation = Annotation(annotation=full_annotation,
                                method="clustalo",
                                blocks=blocks)
        return annotation


def count_diffs(align, feats, inseq, verbose):

    nfeats = len(feats.keys())
    mm = 0
    insr = 0
    dels = 0
    gaps = 0
    match = 0
    lastb = ''
    l = len(align[0]) if len(align[0]) > len(align[1]) else len(align[1])

    for i in range(0, l):
        if align[0][i] == "-" or align[1][i] == "-":
            if align[0][i] == "-":
                insr += 1
                if lastb != '-':
                    gaps += 1
                lastb = "-"
            if align[1][i] == "-":
                dels += 1
                if lastb != '-':
                    gaps += 1
                lastb = "-"
        else:
            lastb = ''
            if align[0][i] != align[1][i]:
                mm += 1
            else:
                match += 1

    gper = gaps / nfeats
    delper = dels / l
    iper = insr / l
    mmper = mm / l
    mper = match / l
    mper2 = match / len(inseq)

    if verbose:
        print(list(feats.keys()), file=sys.stderr)
        print('{:<22}{:<6d}'.format("Number of feats: ", nfeats), file=sys.stderr)
        print('{:<22}{:<6d}{:<1.2f}'.format("Number of gaps: ", gaps, gper), file=sys.stderr)
        print('{:<22}{:<6d}{:<1.2f}'.format("Number of deletions: ", dels, delper), file=sys.stderr)
        print('{:<22}{:<6d}{:<1.2f}'.format("Number of insertions: ", insr, iper), file=sys.stderr)
        print('{:<22}{:<6d}{:<1.2f}'.format("Number of mismatches: ", mm, mmper), file=sys.stderr)
        print('{:<22}{:<6d}{:<1.2f}'.format("Number of matches: ", match, mper), file=sys.stderr)
        print('{:<22}{:<6d}{:<1.2f}'.format("Number of matches: ", match, mper2), file=sys.stderr)
        print("----", file=sys.stderr)
    indel = iper + delper
    if (indel > 0.5 or mmper > 0.05 or gper > .50) and mper2 < .9:
        return Annotation(complete_annotation=False)
    else:
        return insr, dels




