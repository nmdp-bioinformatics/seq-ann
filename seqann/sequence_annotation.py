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

# TODO: change file name to seq_annotation.py
#
import re
import sys

from Bio.Seq import Seq
from Bio import pairwise2
from Bio.Alphabet import IUPAC
from BioSQL.BioSeq import DBSeq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import ExactPosition
from Bio.SeqFeature import FeatureLocation

from seqann.models.reference_data import ReferenceData
from seqann.blast_cmd import blastn
from seqann.blast_cmd import get_locus
from seqann.seq_search import SeqSearch
from seqann.models.base_model_ import Model
from seqann.align import align_seqs
from seqann.util import randomid
from seqann.util import get_features
from seqann.util import get_seqs

from itertools import repeat

import logging

isexon = lambda f: True if re.search("exon", f) else False
isutr = lambda f: True if re.search("UTR", f) else False
isfive = lambda f: True if re.search("five", f) else False

is_classII = lambda x: True if re.search("HLA-D", x) else False

logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    datefmt='%m/%d/%Y %I:%M:%S %p',
                    level=logging.INFO)


class BioSeqAnn(Model):
    '''
    import seqann
    seqanno = seqann.BioSeqAnn()
    annotations = [an.annotate(rec, loc) for rec in list(SeqIO.read(file,'fasta'))]
    '''
    def __init__(self, server=None, dbversion='3310', datfile='',
                 rerun=False, rerun_n=3, verbose=False, verbosity=0,
                 pid='NA', kir=False, align=False):
        self.align = align
        self.kir = kir
        self.server = server
        self.rerun = rerun
        self.rerun_n = rerun_n
        self.verbose = verbose
        self.verbosity = verbosity
        self.refdata = ReferenceData(server=server,
                                     dbversion=dbversion,
                                     alignments=align,
                                     kir=kir)
        self.seqsearch = SeqSearch(refdata=self.refdata, verbose=self.verbose)
        self.logger = logging.getLogger("Logger." + __name__)
        self.logname = "ID {:<10} -".format(str(pid))

    def annotate(self, sequence, locus=None, nseqs=10):
        """
        annotate - method for annotating a Biopython sequence record
        :param sequence: The input sequence record.
        :type sequence: Seq
        :param locus: The gene locus associated with the sequence.
        :type locus: str
        :param nseqs: The number of blast sequences to use.
        :type nseqs: int

        Returns:
            Annotation: complete annotation is successful.

            The annotate function return an ``Annotation`` object that
            contains the sequence features and names associated with them.

            Example output::

                {
                     'annotation': {'exon_1': SeqRecord(seq=Seq('AGAGACTCTCCCG', SingleLetterAlphabet()), id='HLA:HLA00630', name='HLA:HLA00630', description='HLA:HLA00630 DQB1*03:04:01 597 bp', dbxrefs=[]),
                                    'exon_2': SeqRecord(seq=Seq('AGGATTTCGTGTACCAGTTTAAGGCCATGTGCTACTTCACCAACGGGACGGAGC...GAG', SingleLetterAlphabet()), id='HLA:HLA00630', name='HLA:HLA00630', description='HLA:HLA00630 DQB1*03:04:01 597 bp', dbxrefs=[]),
                                    'exon_3': SeqRecord(seq=Seq('TGGAGCCCACAGTGACCATCTCCCCATCCAGGACAGAGGCCCTCAACCACCACA...ATG', SingleLetterAlphabet()), id='HLA:HLA00630', name='<unknown name>', description='HLA:HLA00630', dbxrefs=[])},
                     'complete_annotation': True,
                     'features': {'exon_1': SeqFeature(FeatureLocation(ExactPosition(0), ExactPosition(13), strand=1), type='exon_1'),
                                  'exon_2': SeqFeature(FeatureLocation(ExactPosition(13), ExactPosition(283), strand=1), type='exon_2')},
                     'method': 'nt_search and clustalo',
                     'refmissing': ['five_prime_UTR',
                                    'intron_1',
                                    'intron_2',
                                    'exon_3',
                                    'intron_3',
                                    'exon_4',
                                    'intron_4',
                                    'exon_5',
                                    'intron_5',
                                    'exon_6',
                                    'three_prime_UTR'],
                     'seq': SeqRecord(seq=Seq('AGAGACTCTCCCGAGGATTTCGTGTACCAGTTTAAGGCCATGTGCTACTTCACC...ATG', SingleLetterAlphabet()), id='HLA:HLA00630', name='HLA:HLA00630', description='HLA:HLA00630 DQB1*03:04:01 597 bp', dbxrefs=[])
                }


        Examples:

            >>> from Bio.Seq import Seq
            >>> from Bio.SeqRecord import SeqRecord
            >>> import seqann
            >>> seqrec = SeqRecord(seq=Seq('AGAGACTCTCCCGAGGATTTCGTGTACCAGTTTAAGGCCATGTGCTACTTCACC...ATG', SingleLetterAlphabet()), id='HLA:HLA00630')
            >>> seqann = seqann.BioSeqAnn()
            >>> ann = seqann.annotate(seqrec)
            >>> for f in ann.annotation:
            ...    print(f, ann.method, ann.annotation[f], sep="\t")
            exon_2  nt_search and clustalo  AGGATTTCGTGTACCAGTTTAAGGCCATGTGCTACTTCACCAACGGGACGGAGCGCGTGCGTTATGTGACCAGATACATCTATAACCGAGAGGAGTACGCACGCTTCGACAGCGACGTGGAGGTGTACCGGGCGGTGACGCCGCTGGGGCCGCCTGCCGCCGAGTACTGGAACAGCCAGAAGGAAGTCCTGGAGAGGACCCGGGCGGAGTTGGACACGGTGTGCAGACACAACTACCAGTTGGAGCTCCGCACGACCTTGCAGCGGCGAG
            exon_3  nt_search and clustalo  TGGAGCCCACAGTGACCATCTCCCCATCCAGGACAGAGGCCCTCAACCACCACAACCTGCTGGTCTGCTCAGTGACAGATTTCTATCCAGCCCAGATCAAAGTCCGGTGGTTTCGGAATGACCAGGAGGAGACAACCGGCGTTGTGTCCACCCCCCTTATTAGGAACGGTGACTGGACCTTCCAGATCCTGGTGATGCTGGAAATGACTCCCCAGCATGGAGACGTCTACACCTGCCACGTGGAGCACCCCAGCCTCCAGAACCCCATCACCGTGGAGTGGC
            exon_1  nt_search and clustalo  AGAGACTCTCCCG
            exon_4  nt_search and clustalo  GGGCTCAGTCTGAATCTGCCCAGAGCAAGATG

        """
        # TODO: use logging
        if not locus:
            self.logger.info(self.logname + " No locus provided! ")
            locus = get_locus(sequence, kir=self.kir, refdata=self.refdata)
            self.logger.info(self.logname + " Locus prediction = " + locus)
            if not locus:
                self.logger.error(self.logname + " Locus count not be determined!")
                return ''

        # Exact match found
        matched_annotation = self.refdata.search_refdata(sequence, locus)
        if matched_annotation:
            if self.verbose:
                self.logger.info(self.logname + " exact match found")
            return matched_annotation

        blast = blastn(sequence, locus, nseqs,
                       kir=self.kir, verbose=self.verbose,
                       refdata=self.refdata)

        if blast.failed:
            # TODO: return error object
            # TODO: Add logging
            self.logger.error(self.logname + " blastn failed! " + locus + " " + str(sequence.seq))
            self.logger.error(self.logname + " loci invalid " + locus + " " + str(sequence.seq))
            locus = get_locus(sequence, kir=self.kir, refdata=self.refdata)
            self.logger.info(self.logname + " Locus prediction = " + locus)
            if not locus:
                self.logger.error(self.logname + " Locus couldn't be determined not be determined!")
                return
            return self.annotate(sequence, locus)

        # Do seq_search first on all blast sequences
        partial_ann = None
        found = blast.match_seqs
        for i in range(0, len(found)):
            if self.verbose:
                self.logger.info(self.logname + " running seq_search with " + found[i].name)
            annotation = self.seqsearch.search_seqs(found[i],
                                                    sequence, locus,
                                                    partial_ann=partial_ann)
            if annotation.complete_annotation:
                # TODO: Add options to clean
                # TODO: change clean to cleanup
                if self.verbose:
                    self.logger.info(self.logname + " Finished annotation with " + found[i].name)
                
                # Add alignment is specified
                if self.align:
                    if self.verbose:
                        self.logger.info(self.logname + " Adding alignment")
                    annotation = self.add_alignment(found[i], annotation)

                annotation.clean()
                return annotation
            else:
                partial_ann = annotation
                if self.verbose:
                    self.logger.info(self.logname + " Using partial annotation * run " + str(i) + " *")
                    self.logger.info(self.logname + " Features found = " + ",".join(list(annotation.features.keys())))
                    self.logger.info(self.logname + " Sequence unmapped = " + str(annotation.covered))

        # Now loop through doing alignment
        #for i in range(0, 2):
        #if not is_classII(locus):
        for i in range(0, 2):
            annotation = self.seqsearch.search_seqs(found[i],
                                                    sequence, locus,
                                                    partial_ann=partial_ann)
            aligned_ann = self.ref_align(found[i], sequence, locus,
                                         annotation=annotation)
            if aligned_ann.complete_annotation:
                if self.align:
                    if self.verbose:
                        self.logger.info(self.logname + " Adding alignment")
                    annotation = self.add_alignment(found[i], annotation)
                annotation.clean()
                return aligned_ann
            else:
                if self.verbose:
                    self.logger.info(self.logname + " Using partial annotation for alignment * run " + str(i) + " *")
                partial_ann = aligned_ann

        # TODO: make guess with full alignment
        # return self.ref_align(found, sequence, locus,
        #                       partial_ann=partial_ann)
        self.logger.error(self.logname + " No annotation produced!")
        return ''

    def ref_align(self, found_seqs, sequence, locus,
                  annotation=None, partial_ann=None):
        """
        ref_aling - method for annotating a Biopython sequence record
        :param found_seqs: The input sequence record.
        :type found_seqs: Seq
        :param sequence: The input sequence record.
        :type sequence: Seq
        :param locus: The gene locus associated with the sequence.
        :type locus: str
        :param annotation: The incomplete annotation from a previous iteration.
        :type annotation: Annotation
        :param partial_ann: The partial annotation after looping through all of the blast sequences.
        :type partial_ann: Annotation
        """
        if annotation:
            # Extract the missing blocks and
            # only align those blocks to the known
            # missing features
            # Start with all blocks missing
            # and then delete block if it is found
            missing_blocks = annotation.blocks
            for b in annotation.blocks:
                # **** Check if block equals full input sequence *** #
                # - If it does, then just align the ful
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
                                                           feat,
                                                           b)
                for combseqr in combosrecs:
                    mbtmp = []
                    an, ins, dels = align_seqs(combseqr, feat, locus, verbose=self.verbose)
                    mapped_feat = list(an.annotation.keys())
                    if len(mapped_feat) >= 1:
                        for f in an.annotation:
                            if f in annotation.missing:
                                length, lengthsd = 0, 0
                                length = float(self.refdata.feature_lengths[locus][f][0])
                                lengthsd = float(self.refdata.feature_lengths[locus][f][1])

                                # min and max lengths expected
                                max_length = length + (lengthsd*3) + ins
                                min_length = length - (lengthsd*3) - dels

                                if(len(an.annotation[f]) <= max_length and
                                        len(an.annotation[f]) >= min_length):
                                    annotation.annotation.update({f:
                                                                  an.annotation[f]
                                                                  })
                                    if an.blocks:
                                        mbtmp += an.blocks
                                    else:
                                        if self.verbose:
                                            self.logger.info(self.logname + " Mapped " + f + " with clustalo")
                                        if b in missing_blocks:
                                            del missing_blocks[missing_blocks.index(b)]
                                else:
                                    mbtmp.append(b)
                    else:
                        mbtmp.append(b)
                    annotation.blocks = mbtmp
                    annotation.check_annotation()
                    if annotation.complete_annotation:
                        return annotation

                if len(exons.seq) >= 4:
                    exonan, ins, dels = align_seqs(exons, feat, locus, verbose=self.verbose)
                    mapped_exons = list(exonan.annotation.keys())
                    if len(mapped_exons) >= 1:
                        if self.verbose:
                            self.logger.info(self.logname + " Mapped exons with clustalo")
                        for f in exonan.annotation:
                            annotation.annotation.update({f: exonan.annotation[f]})
                        del missing_blocks[missing_blocks.index(b)]

                        annotation.blocks = missing_blocks
                        annotation.check_annotation()
                        if annotation.complete_annotation:
                            return annotation

                # Run full sequence
                if len(fullrec.seq) >= 4:
                    fullref = align_seqs(fullrec, feat, locus, verbose=self.verbose)
                    if hasattr(fullref, 'annotation'):
                        mapped_full = list(fullref.annotation.keys())
                        if len(mapped_full) >= 1:
                            # If it wasn't found
                            del missing_blocks[missing_blocks.index(b)]

                            for f in fullref.annotation:
                                annotation.update({f: fullref.annotation[f]})

                        annotation.blocks = missing_blocks
                        annotation.check_annotation()
                        if annotation.complete_annotation:
                            if self.verbose:
                                self.logger.info(self.logname + " Mapped all features with clustalo")
                            return annotation
            return annotation
        elif partial_ann:
            self.logger.error(self.logname + " Partial alignment! SHOULDNT BE HERE!!")
            # Do full sequence alignments
            # any only extract out the part
            # that couldn't be explained from above
            return
        else:
            self.logger.error(self.logname + " Full alignment! SHOULDNT BE HERE!!")
            # Option for user that just want to
            # align a set of sequences
            return

    def add_alignment(self, ref_seq, annotation):
        """
        ref_aling - method for annotating a Biopython sequence record
        :param found_seqs: The input sequence record.
        :type found_seqs: Seq
        :param sequence: The input sequence record.
        :type sequence: Seq
        :param locus: The gene locus associated with the sequence.
        :type locus: str
        :param annotation: The incomplete annotation from a previous iteration.
        :type annotation: Annotation
        :param partial_ann: The partial annotation after looping through all of the blast sequences.
        :type partial_ann: Annotation
        """
        seq_features = get_seqs(ref_seq)
        annoated_align = {}
        allele = ref_seq.description.split(",")[0]
        locus = allele.split("*")[0].split("-")[1]
        for feat in seq_features:
            if feat in annotation.annotation:
                if isinstance(annotation.annotation[feat], DBSeq):
                    seq_len = len(str(annotation.annotation[feat]))
                    ref_len = len(seq_features[feat])
                else:
                    seq_len = len(str(annotation.annotation[feat].seq))
                    ref_len = len(seq_features[feat])
                if seq_len == ref_len:
                    seq = list(annotation.annotation[feat].seq)
                    gaps = self.refdata.annoated_alignments[locus][allele][feat]['Gaps']
                    if self.verbose and self.verbosity > 0:
                        self.logger.info(self.logname + " Gaps at " + feat)
                        self.logger.info(self.logname +
                            "-".join([",".join([str(s) for s in g]) for g in gaps]))
                    for i in range(0, len(gaps)):
                        for j in gaps[i]:
                            loc = j
                            seq.insert(loc, '-')
                    nseq = ''.join(seq)
                    annoated_align.update({feat: nseq})
                else:
                    in_seq = str(annotation.annotation[feat].seq)
                    ref_seq = self.refdata.annoated_alignments[locus][allele][feat]['Seq']
                    # TODO: add logging
                    # if in_seq == ref_len:
                    #     alignment = pairwise2.align.globalxx(self.refdata.annoated_alignments[locus][allele][feat]['Seq'], in_seq)
                    #     if self.verbose and self.verbosity > 0:
                    #         self.logger.info(self.logname + " Align1 -> in_seq == ref_len " + feat)
                    #         self.logger.info(self.logname + " " + str(len(in_seq)) + "==" + str(ref_len))
                    #         self.logger.info(alignment)
                    #     annoated_align.update({feat: alignment[0][0]})
                    # else:
                    alignment = pairwise2.align.globalxx(in_seq, ref_seq)
                    if self.verbose and self.verbosity > 0:
                        self.logger.info(self.logname + " Align2 -> in_seq != ref_len " + feat)
                        self.logger.info(self.logname + " " + str(len(in_seq)) + " == " + str(ref_len))
                    annoated_align.update({feat: alignment[0][0]})
            else:
                nseq = ''.join(list(repeat('-', len(seq_features[feat]))))
                annoated_align.update({feat: nseq})
        annotation.aligned = annoated_align
        return annotation

    def _refseqs(self, locus, start_pos, annotation, feat, block):
        # refseqs = _refseqs(locus, start, annotation,
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
        # return blank if missing features
        # If exon only, then only extract exon sequences
        end_pos = start_pos + len(feat.seq)
        missing_feats = annotation.missing
        mapping = annotation.mapping
        if start_pos == 0:
            prv_order = -1
        else:
            prev_feat = mapping[start_pos-1]
            prv_order = self.refdata.structures[locus][prev_feat]

        if len(mapping) == end_pos+1:
            nxt_order = len(self.refdata.structures[locus])+1
        else:
            next_feat = mapping[end_pos+1]
            nxt_order = self.refdata.structures[locus][next_feat]

        start = 0
        exstart = 0
        all_seqs = []
        all_feats = []
        exons_only = []
        exon_feats = []
        all_seqrecs = {}
        for f in sorted(missing_feats,
                        key=lambda k: self.refdata.structures[locus][k]):
            mis_order = self.refdata.structures[locus][f]
            if mis_order > nxt_order or mis_order < prv_order:
                continue
            else:
                all_seqs.append(annotation.missing[f])
                seq_feat, start = self._make_seqfeat(start,
                                                     annotation.missing[f],
                                                     f)
                all_feats.append(seq_feat)

                one_seqfeat, x = self._make_seqfeat(0,
                                                    annotation.missing[f],
                                                    f)

                ref1 = Seq(str(annotation.missing[f]), IUPAC.unambiguous_dna)
                rid1 = "Ref1_" + str(f) + "_" + str(randomid(N=2))
                refrec1 = SeqRecord(seq=ref1, features=[one_seqfeat], id=rid1)
                all_seqrecs.update({f: refrec1})

                if isexon(f):
                    ftx, exstart = self._make_seqfeat(exstart,
                                                      annotation.missing[f],
                                                      f)
                    exons_only.append(annotation.missing[f])
                    exon_feats.append(ftx)

        # Exons only
        exononly_seq = "".join([str(ft) for ft in exons_only])
        refseq_exons = Seq(exononly_seq, IUPAC.unambiguous_dna)
        refid_exons = "RefExons_" + str(randomid())
        refrec_exons = SeqRecord(seq=refseq_exons, features=exon_feats,
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

            if start_ord == nxt_order-1:
                nfeat = self.refdata.struct_order[locus][start_ord]
                if nfeat not in all_seqrecs:
                    continue
                rec = all_seqrecs[nfeat]

                length = float(self.refdata.feature_lengths[locus][f][0])
                lengthsd = float(self.refdata.feature_lengths[locus][f][1])

                max_length = length + lengthsd
                min_length = length - lengthsd
                if len(rec.seq) <= max_length \
                        and len(rec.seq) >= min_length:
                    recf, x = self._make_seqfeat(0, rec.seq, nfeat)
                    ctmpid1 = "Ref1_" + str(nfeat) + "_" + str(randomid(N=2))
                    tmprec1 = SeqRecord(seq=rec.seq, features=[recf],
                                        id=ctmpid1)
                    combos.append(tmprec1)
            else:
                cstart = 0
                ref_feats = []
                combo_seq = []
                feat_names = []
                combo_feats = []
                for i in range(start_ord, nxt_order):
                    nfeat = self.refdata.struct_order[locus][i]
                    if nfeat not in all_seqrecs:
                        continue
                    ref_feats.append(nfeat)

                    length, lengthsd = 0, 0

                    for fr in ref_feats:
                        length += float(self.refdata.feature_lengths[locus][fr][0])
                        lengthsd += float(self.refdata.feature_lengths[locus][fr][1])

                    # min and max lengths expected
                    max_length = length + lengthsd
                    min_length = length - lengthsd

                    rec = all_seqrecs[nfeat]
                    feat_names.append(nfeat)

                    cfeat, cstart = self._make_seqfeat(cstart, rec.seq, nfeat)
                    combo_feats.append(cfeat)
                    combo_seq.append(rec.seq)

                    seqtmp = "".join([str(ft) for ft in combo_seq])
                    ctmpseq = Seq(seqtmp, IUPAC.unambiguous_dna)
                    ctmpid = "|".join([str(ft) for ft in feat_names]) \
                             + "_" + str(randomid(N=2))

                    if len(ctmpseq) <= max_length \
                            and len(ctmpseq) >= min_length:
                        tmprec = SeqRecord(seq=ctmpseq, features=combo_feats,
                                           id=ctmpid)
                        combos.append(tmprec)
        if combos:
            # Sort combos by how close they
            # are in length to the input sequence block
            # * may need to reverse sort *
            sorted_combos = sorted(combos,
                                   key=lambda rec: abs(len(rec.seq) - len(feat)))
        else:
            sorted_combos = []

        return sorted_combos, refrec_exons, refrec

    def _make_seqfeat(self, start, sequence, featname):

        if isutr(featname):
            ftype = "UTR"
            end = start + len(sequence)
            seq_feat = \
                SeqFeature(
                    FeatureLocation(
                        ExactPosition(start),
                        ExactPosition(end),
                        strand=1),
                    type=ftype)
            start = end
            return seq_feat, start
        else:
            ftype, rank = featname.split("_")
            end = start + len(sequence)
            quals = {'number': [str(rank)]}
            seq_feat = \
                SeqFeature(
                    FeatureLocation(
                        ExactPosition(start),
                        ExactPosition(end),
                        strand=1),
                    type=ftype,
                    qualifiers=quals)
            start = end
            return seq_feat, start


class SeqAnnException(Exception):

    def __init__(self, status=None, reason=None, http_resp=None):
        if http_resp:
            self.status = http_resp.status
            self.reason = http_resp.reason
            self.body = http_resp.data
            self.headers = http_resp.getheaders()
        else:
            self.status = status
            self.reason = reason
            self.body = None
            self.headers = None

    def __str__(self):
        """
        Custom error messages for exception
        """
        error_message = "({0})\n"\
                        "Reason: {1}\n".format(self.status, self.reason)
        if self.headers:
            error_message += "HTTP response headers: {0}\n".format(self.headers)

        if self.body:
            error_message += "HTTP response body: {0}\n".format(self.body)

        return error_message



