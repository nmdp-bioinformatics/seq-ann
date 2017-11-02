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
            refseqrec = self._refseqs(annotation)
            if not refseqrec:
                return annotation
            else:
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

                        # refseqs = self._refseqs(locus, start, annotation,
                        #                          sequence.seq)
                        # - what is size of input seq
                        # - whst is size of block
                        # - what is features have been found
                        # - Of the features that have not been found
                        #   and follow these rules:
                        #       - if it's at the start of the input seq
                        #         then skip anything that comes after the
                        #         block and is missing
                        #       - if it's at the end of the input seq
                        #       - then skip anything that comes before the
                        #         block and is missing
                        #       - Try first all combos of features where
                        #         len(block) <= ave len blocks - sd or
                        #         len(block) >= len blocks + sd


                        # print("FEATS")
                        # print(b[0]-1)
                        # print(b[len(b)-1])
                        # print(b)

                        # TODO: Skip if input seq is too small to align
                        an = align_seqs(refseqrec, feat, locus)
                        mapped_feat = list(an.annotation.keys())
                        if len(mapped_feat) >= 1:
                            print("Block found")
                            #print(b)
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

    def _refseqs(self, annotation):

        # return blank if missing features
        # If exon only, then only extract exon sequences
        
        start = 0
        missing_feats = []
        missing_seqs = []
        exon_only = []
        exon_seqs = []
        for b in annotation.missing:


            missing_seqs.append(annotation.missing[b][1])
            end = start + len(annotation.missing[b][1])
            seq_feat = \
                SeqFeature(
                    FeatureLocation(
                        ExactPosition(start),
                        ExactPosition(end),
                        strand=1),
                    type=annotation.missing[b][2])
            start = end
            missing_feats.append(seq_feat)
        seq = "".join([str(f) for f in missing_seqs])
        refseq = Seq(seq, IUPAC.unambiguous_dna)
        refid = "RefSeq_" + str(randomid())
        refrec = SeqRecord(seq=refseq, features=missing_feats, id=refid)
        print(refrec)
        return SeqRecord(seq=refseq, features=missing_feats, id=refid)


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












