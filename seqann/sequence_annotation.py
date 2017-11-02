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

from seqann.models.annotation import Annotation
from seqann.models.reference_data import ReferenceData
from seqann.blast_cmd import blastn
from seqann.seq_search import SeqSearch
from seqann.models.base_model_ import Model

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
from Bio.Alphabet import IUPAC
from Bio import SearchIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.SeqUtils import nt_search
import os
import glob
from Bio import AlignIO
from BioSQL import BioSeqDatabase
import re
from seqann.align import align_seqs
from seqann.util import randomid

isexon = lambda f: True if re.search("exon", f) else False


class BioSeqAnn(Model):
    '''
    import seqann
    seqanno = seqann.BioSeqAnn()
    annotations = [an.annotate(rec, loc) for rec in list(SeqIO.read(file,'fasta'))]

    '''
    def __init__(self, server=None, dbversion='3290', datfile='', rerun=False, rerun_n=3):
        self.server = server
        self.rerun = rerun
        self.rerun_n = rerun_n
        self.refdata = ReferenceData(server=server,
                                     dbversion=dbversion)
        self.seqsearch = SeqSearch(refdata=self.refdata)

    def annotate(self, sequence, locus, nseqs=8):
        matched_annotation = self.refdata.search_refdata(sequence, locus)
        if matched_annotation:
            return matched_annotation

        blast = blastn(sequence, locus, nseqs, refdata=self.refdata)

        if blast.failed:
            return self.broad_align(sequence, locus)

        partial_ann = None
        found = blast.match_seqs
        for i in range(0, len(found)):
            annotation = self.seqsearch.search_seqs(found[i],
                                                    sequence, locus,
                                                    partial_ann=partial_ann)
            if annotation.complete_annotation:
                return annotation
            else:
                aligned_ann = self.ref_align(found[i], sequence, locus,
                                             nseqs, annotation=annotation)
                if aligned_ann.complete_annotation:
                    return aligned_ann
                else:
                    print("Using partial annotation")
                    partial_ann = aligned_ann
        return self.ref_align(found, sequence, locus, nseqs,
                              partial_ann=partial_ann)

    def ref_align(self, found_seqs, sequence, locus, nseqs=3, annotation=None,
                  partial_ann=None):

        if annotation:
            # Extract the missing blocks and
            # only align those blocks to the known
            # missing features
            missing_blocks = []
            for b in annotation.blocks:
                    start = b[0]-1 if b[0] != 0 else 0
                    seq_feat = \
                        SeqFeature(
                            FeatureLocation(
                                ExactPosition(start),
                                ExactPosition(b[len(b)-1]),
                                strand=1),
                            type="unmapped")

                    feat = seq_feat.extract(annotation.seq)
                    combosrecs, exons, fullrec = self._refseqs(locus,
                                                               start,
                                                               annotation,
                                                               sequence.seq,
                                                               feat,
                                                               b)



                    an = align_seqs(refseqrec, feat, locus)
                    mapped_feat = list(an.annotation.keys())
                    if len(mapped_feat) >= 1:
                        print("Block found")
                    # refseqs = self._refseqs(locus, start, annotation,
                    #                          sequence.seq, feat,
                    #                           b)
                    # - what is size of input seq
                    # - whst is size of block
                    # - what is features have been found
                    # - Of the features that have not been found
                    #   and follow these rules:
                    #       ** ONLY align features that are expected
                    #           this block
                    #           - if 5'UTR and exon1 are ambigous
                    #             but intron1 is known, only align to
                    #             5'UTR and exon1
                    #       - LOOK for any features found above and
                    #         behind block. Only USE features within that
                    #         span.
                    #
                    #       - if it's at the start of the input seq
                    #         then skip anything that comes after the
                    #         block and is missing
                    #       - if it's at the end of the input seq
                    #       - then skip anything that comes before the
                    #         block and is missing
                    #       - Try first all combos of features where
                    #         len(block) <= ave len blocks - sd or
                    #         len(block) >= len blocks + sd
                    #
                    # print("FEATS")
                    # print(b[0]-1)
                    # print(b[len(b)-1])
                    # print(b)
                    #
                    #
                    an = align_seqs(refseqrec, feat, locus)
                    mapped_feat = list(an.annotation.keys())
                    if len(mapped_feat) >= 1:
                        print("Block found")
                    # print(b)
                    # TODO: Skip if input seq is
                    # too small to align
                        for f in an.annotation:
                            print("Found feat " + f)
                            print(str(an.annotation[f].seq))
                            print("")
                            annotation.annotation.update({f: an.annotation[f]})
                    else:
                        print("Block missing")
                        print(b)
                        missing_blocks.append(b)
            annotation.blocks = missing_blocks
            annotation.check_annotation()
            return annotation
        else:
            print("ref_align2")
            alignments = align_seqs(found_seqs, sequence, locus)
            if(alignments.complete_annotation):
                return alignments.annotation
            elif(self.rerun):
                nseqs_toblast = nseqs + self.rerun_n
                return self.get_features(sequence, locus, nseqs_toblast)
            else:
                print("failed")
                return
            #     return self._build_error('ref_align', found_seqs,
            #                              sequence, locus, nseqs)

    def broad_align(self, sequence, locus, N=8):
        #if self._valid(sequence):
            # refseqs = self.refdata.refseqs(locus, N)
            # aligned_ann = align_seqs(refseqs, sequence, locus)
            # if(aligned_ann.complete_annotation):
            #     return aligned_ann.annotation
            # else:
            #     print("Failed")
            #     return
        #else:
        #    print("Failed")
        print("Failed")
        return

    def _refseqs(self, locus, start_pos, annotation, inputseq, feat, block):

#                       # refseqs = self._refseqs(locus, start, annotation,
                        #                          sequence.seq, feat,
                        #                           b)
        # return blank if missing features
        # If exon only, then only extract exon sequences
        end_pos = start_pos + len(feat.seq)
        missing_feats = annotation.missing
        mapping = annotation.mapping
        next_feat = mapping[end_pos+1]
        prev_feat = mapping[start_pos-1]
        nxt_order = self.refdata.structure_order[next_feat]
        prv_order = self.refdata.structure_order[prev_feat]

        start = 0
        exstart = 0
        all_seqs = []
        all_feats = []
        exons_only = []
        exon_feats = []
        all_seqrecs = []
        for f in sorted(missing_feats,
                        key=lambda k: self.refdata.structures[locus][k]):
            mis_order = self.refdata.structures[locus][f]
            if mis_order > nxt_order or mis_order < prv_order:
                continue
            else:
                all_seqs.append(annotation.missing[f])
                seq_feat, start = self._make_seqfeat(start,
                                                     annotation.missing[f][1],
                                                     annotation.missing[f][2])
                all_feats.append(seq_feat)

                one_seqfeat, x = self._make_seqfeat(0,
                                                    annotation.missing[f][1],
                                                    annotation.missing[f][2])

                ref1 = Seq(str(annotation.missing[f]), IUPAC.unambiguous_dna)
                rid1 = "Ref1_" + str(f) + "_" + str(randomid(N=2))
                refrec1 = SeqRecord(seq=ref1, features=one_seqfeat, id=rid1)
                all_seqrecs.update({annotation.missing[f][2]: refrec1})

                if isexon(f):
                    ftx, exstart = self._make_seqfeat(exstart,
                                                      annotation.missing[f][1],
                                                      annotation.missing[f][2])
                    exons_only.append(annotation.missing[f][1])
                    exon_feats.append(ftx)

        # Exons only
        exononly_seq = "".join([str(ft) for ft in exons_only])
        refseq_exons = Seq(exononly_seq, IUPAC.unambiguous_dna)
        refid_exons = "RefExons_" + str(randomid())
        refrec_exons = SeqRecord(seq=refseq_exons, features=missing_feats,
                                 id=refid_exons)
        # All seqs
        seq = "".join([str(ft) for ft in all_seqs])
        refseq = Seq(seq, IUPAC.unambiguous_dna)
        refid = "RefAll_" + str(randomid())
        refrec = SeqRecord(seq=refseq, features=all_feats, id=refid)

        # ** Creat all combos of seqrecs and sort by
        #    which are most similar in length, gc content and molecular weight
        combos = []
        for f in sorted(all_seqrecs,
                        key=lambda k: self.refdata.structures[locus][k]):
            start_ord = self.refdata.structures[locus][f]

            if start_ord != nxt_order-1:
                nfeat = self.refdata.structure_order[locus][start_ord]
                rec = all_seqrecs[nfeat]

                length = self.refdata.structure_lengths[locus][f]['length']
                lengthsd = self.refdata.structure_lengths[locus][f]['sd']

                # min and max lengths expected
                max_length = length + lengthsd
                min_length = length - lengthsd
                if len(rec.seq) <= max_length \
                        and len(rec.seq) >= min_length \
                        and len(rec.seq) > 10:

                    recf, x = self._make_seqfeat(0, rec.seq, nfeat)
                    ctmpid1 = "Ref1_" + str(nfeat) + "_" + str(randomid(N=2))
                    tmprec1 = SeqRecord(seq=rec.seq, features=recf,
                                        id=ctmpid1)
                    combos.append(tmprec1)
            else:

                cstart = 0
                ref_feats = []
                combo_seq = []
                combo_feats = []
                for i in range(start_ord, nxt_order):
                    nfeat = self.refdata.structure_order[locus][i]
                    ref_feats.append(f)

                    length, lengthsd = 0, 0
                    for f in ref_feats:
                        length += self.refdata.structure_lengths[locus][f]['length']
                        lengthsd += self.refdata.structure_lengths[locus][f]['sd']

                    # min and max lengths expected
                    max_length = length + lengthsd
                    min_length = length - lengthsd

                    rec = all_seqrecs[nfeat]
                    combo_feats.append(recf)
                    combo_seq.append(rec.seq)

                    seqtmp = "".join([str(ft) for ft in combo_seq])
                    ctmpseq = Seq(seqtmp, IUPAC.unambiguous_dna)

                    if len(ctmpseq.seq) <= max_length \
                            and len(ctmpseq.seq) >= min_length \
                            and len(ctmpseq.seq) > 10:
                        recf, cstart = self._make_seqfeat(start,
                                                          rec.seq, nfeat)

                        ctmpid = "RefComboPart_" + str(randomid())
                        tmprec = SeqRecord(seq=ctmpseq, features=combo_feats,
                                           id=ctmpid)
                        combos.append(tmprec)

                combseq = "".join([str(ft) for ft in combo_seq])
                crefseq = Seq(combseq, IUPAC.unambiguous_dna)
                crefid = "RefCombo_" + str(randomid())
                crefrec = SeqRecord(seq=crefseq, features=combo_feats,
                                    id=crefid)
                combos.append(crefrec)

        return combos, refrec_exons, refrec

    def _make_seqfeat(self, start, sequence, featname):

        end = start + len(sequence)
        seq_feat = \
            SeqFeature(
                FeatureLocation(
                    ExactPosition(start),
                    ExactPosition(end),
                    strand=1),
                type=featname)
        start = end
        return seq_feat, start

        # start = 0
        # missing_feats = []
        # missing_seqs = []
        # exon_only = []
        # exon_seqs = []
        # for b in annotation.missing:


        #     missing_seqs.append(annotation.missing[b][1])
        #     end = start + len(annotation.missing[b][1])
        #     seq_feat = \
        #         SeqFeature(
        #             FeatureLocation(
        #                 ExactPosition(start),
        #                 ExactPosition(end),
        #                 strand=1),
        #             type=annotation.missing[b][2])
        #     start = end
        #     missing_feats.append(seq_feat)
        # seq = "".join([str(f) for f in missing_seqs])
        # refseq = Seq(seq, IUPAC.unambiguous_dna)
        # refid = "RefSeq_" + str(randomid())
        # refrec = SeqRecord(seq=refseq, features=missing_feats, id=refid)
        # print(refrec)
        # return SeqRecord(seq=refseq, features=missing_feats, id=refid)


    # def _refseqs(self, annotation):

    #     # return blank if missing features
    #     # If exon only, then only extract exon sequences
        

    #     start = 0
    #     missing_feats = []
    #     missing_seqs = []
    #     exon_only = []
    #     exon_seqs = []
    #     for b in annotation.missing:


    #         missing_seqs.append(annotation.missing[b][1])
    #         end = start + len(annotation.missing[b][1])
    #         seq_feat = \
    #             SeqFeature(
    #                 FeatureLocation(
    #                     ExactPosition(start),
    #                     ExactPosition(end),
    #                     strand=1),
    #                 type=annotation.missing[b][2])
    #         start = end
    #         missing_feats.append(seq_feat)
    #     seq = "".join([str(f) for f in missing_seqs])
    #     refseq = Seq(seq, IUPAC.unambiguous_dna)
    #     refid = "RefSeq_" + str(randomid())
    #     refrec = SeqRecord(seq=refseq, features=missing_feats, id=refid)
    #     print(refrec)
    #     return SeqRecord(seq=refseq, features=missing_feats, id=refid)

    # def _build_error(self, type, **kargs):












