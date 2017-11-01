# -*- coding: utf-8 -*-

#
#    gene_feature_enumeration Gene Feature Enumeration.
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

from datetime import date, datetime
from typing import List, Dict
from seqann.util import deserialize_model, randomid, cleanup

from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from collections import OrderedDict
from Bio.Alphabet import IUPAC
from seqann.models.annotation import Annotation
import re


def align_seq(found_seq, sequence, locus, 
              match=2, mismatch=-1,
              startgap=-10,
              extendgap=-3):

    alignments = pairwise2.align.globalms(seq1, seq2, 2, -1, -10, -2)


def align_seqs(found_seqs, sequence, locus):

    randid = randomid()
    input_fasta = str(randid) + ".fasta"
    output_clu = str(randid) + ".clu"
    seqs = []
    seqs.append(found_seqs)
    print(found_seqs)
    seqs.append(sequence)
    SeqIO.write(seqs, input_fasta, "fasta")
    clustalomega_cline = ClustalOmegaCommandline(infile=input_fasta,
                                                 outfile=output_clu,
                                                 outfmt='clu', wrap=50000,
                                                 verbose=True, auto=True)
    stdout, stderr = clustalomega_cline()
    align = AlignIO.read(output_clu, "clustal")
    print("******************************")
    for a in align:
        print(str(a.seq))
    print("******************************")
    all_features = []
    if len(align)-2 == 0:
        feats = [feat for feat in seqs[0].features if feat.type != "source" and feat.type != "CDS" and isinstance(feat, SeqFeature)]
        l = 0
        feat_names = {}
        feat_nums = {}
        for f in feats:
            l += 1
            feat_name = ''
            if re.search("\d", f.type) or f.type == "five_prime_UTR" \
                    or f.type == "three_prime_UTR":
                        feat_names.update({f.type: f})
            else:
                if f.type not in feat_nums:
                    if(f.type == "UTR"):
                        feat_name = "five_prime_UTR"
                        feat_nums.update({f.type: 1})
                        feat_names.update({feat_name: f})
                    else:
                        feat_name = f.type + "_" + str(1)
                        feat_nums.update({f.type: 1})
                        feat_names.update({feat_name: f})
                else:
                    if(f.type == "UTR"):
                        feat_name = "three_prime_UTR"
                        feat_names.update({feat_name: f})
                    else:
                        num = feat_nums[f.type] + 1
                        feat_name = f.type + "_" + str(num)
                        feat_nums[f.type] = num
                        feat_names.update({feat_name: f})

        f = get_features(feat_names, align[0])
        all_features.append(f)
    else:
        for i in range(0, len(align)-2):
            feats = [feat for feat in seqs[i].features if feat.type != "source" and feat.type != "CDS" and isinstance(feat, SeqFeature)]
            l = 0
            feat_names = {}
            feat_nums = {}
            for f in feats:
                l += 1
                feat_name = ''
                if re.search("\d", f.type):
                    feat_names.update({f.type: f})
                else:
                    if f.type not in feat_nums:
                        if(f.type == "UTR"):
                            feat_name = "five_prime_UTR"
                            feat_nums.update({f.type: 1})
                            feat_names.update({feat_name: f})
                        else:
                            feat_name = f.type + "_" + str(1)
                            feat_nums.update({f.type: 1})
                            feat_names.update({feat_name: f})
                    else:
                        if(f.type == "UTR"):
                            feat_name = "three_prime_UTR"
                            feat_names.update({feat_name: f})
                        else:
                            num = feat_nums[f.type] + 1
                            feat_name = f.type + "_" + str(num)
                            feat_nums[f.type] = num
                            feat_names.update({feat_name: f})

            f = get_features(feat_names, align[i])
            all_features.append(f)

    cleanup(randid)

    return resolve_feats(all_features, align[len(align)-1])


def get_features(feats, sequ):
    feats_a = list(feats.keys())
    j = 0
    for i in range(0, len(sequ)):
        if j <= len(feats_a)-1:
            if i > int(feats[feats_a[j]].location.end):
                j += 1
            if(sequ[i] == '-'):
                feats[feats_a[j]] = SeqFeature(FeatureLocation(ExactPosition(feats[feats_a[j]].location.start), ExactPosition(int(feats[feats_a[j]].location.end + 1)), strand=1), type=feats[feats_a[j]].type)
                if j != len(feats_a):
                    for l in range(j+1, len(feats_a)):
                        feats[feats_a[l]] = SeqFeature(FeatureLocation(ExactPosition(feats[feats_a[l]].location.start+1), ExactPosition(int(feats[feats_a[l]].location.end + 1)), strand=1), type=feats[feats_a[l]].type)
    return feats


def resolve_feats(feat_list, seq):
    if len(feat_list) > 1:
        print("****** resolve_feats error ******")
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
            if re.search("-", str(seqrec.seq)):
                newseq = re.sub(r'-', '', str(seqrec.seq))
                seqrec.seq = Seq(newseq, IUPAC.unambiguous_dna)
            if seqrec.seq:
                print(seqrec.seq)
                full_annotation.update({feat: seqrec})
        annotation = Annotation(annotation=full_annotation, method="clustalo")
        return annotation



