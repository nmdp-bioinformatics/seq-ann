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
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.Alphabet import IUPAC
from BioSQL.BioSeq import DBSeq
from BioSQL import BioSeqDatabase
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import ExactPosition
from Bio.SeqFeature import FeatureLocation

from seqann.models.annotation import Annotation
from seqann.models.reference_data import ReferenceData
from seqann.blast_cmd import blastn
from seqann.blast_cmd import get_locus
from seqann.seq_search import SeqSearch
from seqann.models.base_model_ import Model
from seqann.align import align_seqs
from seqann.util import randomid
from seqann.util import get_seqs
from seqann.util import isexon
from seqann.util import isutr
from seqann.gfe import GFE

from itertools import repeat
from typing import Dict
from typing import List

import logging
import warnings

from Bio import BiopythonExperimentalWarning
warnings.simplefilter("ignore", BiopythonExperimentalWarning)


class BioSeqAnn(Model):
    '''
        ::

            from seqann import BioSeqAnn
            seqann = BioSeqAnn()

        :param server: A BioSQL database to use for retriving the sequence features. Using a BioSQL DB will speed up the annotations dramatically.
        :type server: BioSeqDatabase
        :param dbversion: The IPD-IMGT/HLA or KIR database release.
        :type dbversion: str
        :param datfile: The IPD-IMGT/HLA or KIR dat file to use in place of the server parameter.
        :type datfile: str
        :param pid: A process label that can be provided to help track the logging output.
        :type pid: str
        :param load_features: Flag for downloading all gene features and accessions from the feature service.
        :type load_features: bool
        :param store_features: Flag for caching all features and their corresponding accessions.
        :type store_features: bool
        :param cached_features: Dictionary containing all the features from the feature service.
        :type cached_features: Dict
        :param kir: Flag for indicating the input sequences are from the KIR gene system. 
        :type kir: bool
        :param align: Flag for producing the alignments along with the annotations.
        :type align: bool
        :param verbose: Flag for running in verbose mode.
        :type verbose: bool
        :param verbosity: Numerical value to indicate how verbose the output will be in verbose mode.
        :type verbosity: int
        :param debug: Dictionary containing names of steps that you want to debug.
        :type debug: Dict
    '''
    def __init__(self, server: BioSeqDatabase=None, dbversion: str='3310',
                 datfile: str='', verbose: bool=False, verbosity: int=0,
                 pid: str='NA', kir: bool=False, align: bool=False,
                 load_features: bool=False, store_features: bool=False,
                 cached_features: Dict=None,
                 debug: Dict=None):

        self.kir = kir
        self.align = align
        self.debug = debug
        self.server = server
        self.verbose = verbose
        self.verbosity = verbosity
        self.align_verbose = verbose
        self.align_verbosity = verbosity

        gfe_verbose = verbose
        gfe_verbosity = verbosity
        refdata_verbose = verbose
        seqsearch_verbose = verbose
        refdata_verbosity = verbosity
        seqsearch_verbosity = verbosity

        # Run with the most possible amount
        # of logging for the provided steps
        if self.debug:
            if 'seq_search' in self.debug:
                seqsearch_verbose = True
                seqsearch_verbosity = debug['seq_search']
            else:
                seqsearch_verbose = False
                seqsearch_verbosity = 0

            if 'refdata' in self.debug:
                refdata_verbose = True
                refdata_verbosity = debug['refdata']
            else:
                refdata_verbose = False
                refdata_verbosity = 0

            if 'seqann' in self.debug:
                self.verbose = True
                self.verbosity = debug['seqann']
            else:
                self.verbose = False
                self.verbosity = 0

            if 'align' in self.debug:
                self.align_verbose = True
                self.align_verbosity = debug['align']
            else:
                self.align_verbose = False
                self.align_verbosity = 0

            if 'gfe' in self.debug:
                gfe_verbose = True
                gfe_verbosity = debug['gfe']
            else:
                gfe_verbose = False
                gfe_verbosity = 0

        self.logger = logging.getLogger("Logger." + __name__)
        self.logname = "ID {:<10} -".format(str(pid))

        # Initalize GFE
        self.gfe = GFE(pid=pid, verbose=gfe_verbose,
                       verbosity=gfe_verbosity,
                       load_features=load_features,
                       store_features=store_features,
                       cached_features=cached_features)

        # Initalize reference data
        self.refdata = ReferenceData(server=server,
                                     dbversion=dbversion,
                                     alignments=align,
                                     verbose=refdata_verbose,
                                     verbosity=refdata_verbosity,
                                     kir=kir)

        # Initalize SeqSearch for matching features
        self.seqsearch = SeqSearch(refdata=self.refdata,
                                   verbose=seqsearch_verbose,
                                   verbosity=seqsearch_verbosity)

    def annotate(self, sequence: Seq=None, locus: str=None,
                 nseqs: int=10, alignseqs: int=10,
                 skip: List=[]) -> Annotation:
        """
        annotate - method for annotating a BioPython sequence

        :param sequence: The input consensus sequence.
        :type sequence: Seq
        :param locus: The gene locus associated with the sequence.
        :type locus: str
        :param nseqs: The number of blast sequences to use.
        :type nseqs: int
        :param alignseqs: The number of sequences to use for targeted alignments.
        :type alignseqs: int
        :param skip: A list of alleles to skip for using as a reference. This is used for validation and testing.
        :type skip: List
        :rtype: Annotation

        Returns:
            The annotate function return an ``Annotation`` object that
            contains the sequence features and names associated with them.

            Example output::

                {
                     'complete_annotation': True,
                     'annotation': {'exon_1': SeqRecord(seq=Seq('AGAGACTCTCCCG', SingleLetterAlphabet()), id='HLA:HLA00630', name='HLA:HLA00630', description='HLA:HLA00630 DQB1*03:04:01 597 bp', dbxrefs=[]),
                                    'exon_2': SeqRecord(seq=Seq('AGGATTTCGTGTACCAGTTTAAGGCCATGTGCTACTTCACCAACGGGACGGAGC...GAG', SingleLetterAlphabet()), id='HLA:HLA00630', name='HLA:HLA00630', description='HLA:HLA00630 DQB1*03:04:01 597 bp', dbxrefs=[]),
                                    'exon_3': SeqRecord(seq=Seq('TGGAGCCCACAGTGACCATCTCCCCATCCAGGACAGAGGCCCTCAACCACCACA...ATG', SingleLetterAlphabet()), id='HLA:HLA00630', name='<unknown name>', description='HLA:HLA00630', dbxrefs=[])},
                     'features': {'exon_1': SeqFeature(FeatureLocation(ExactPosition(0), ExactPosition(13), strand=1), type='exon_1'),
                                  'exon_2': SeqFeature(FeatureLocation(ExactPosition(13), ExactPosition(283), strand=1), type='exon_2')
                                  'exon_3': SeqFeature(FeatureLocation(ExactPosition(283), ExactPosition(503), strand=1), type='exon_3')},
                     'method': 'nt_search and clustalo',
                     'gfe': 'HLA-Aw2-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-4',
                     'seq': SeqRecord(seq=Seq('AGAGACTCTCCCGAGGATTTCGTGTACCAGTTTAAGGCCATGTGCTACTTCACC...ATG', SingleLetterAlphabet()), id='HLA:HLA00630', name='HLA:HLA00630', description='HLA:HLA00630 DQB1*03:04:01 597 bp', dbxrefs=[])
                }


        Example usage:

            >>> from Bio.Seq import Seq
            >>> from Bio.SeqRecord import SeqRecord
            >>> from seqann import BioSeqAnn
            >>> sequence = Seq('AGAGACTCTCCCGAGGATTTCGTGTACCAGTTTAAGGCCATGTGCTACTTCACC...ATG', SingleLetterAlphabet())
            >>> seqann = BioSeqAnn()
            >>> ann = seqann.annotate(sequence)
            >>> for f in ann.annotation:
            ...    print(f, ann.method, ann.annotation[f], sep="\t")
            exon_2  nt_search and clustalo  AGGATTTCGTGTACCAGTTTAAGGCCATGTGCTACTTCACCAACGGGACGGAGCGCGTGCGTTATGTGACCAGATACATCTATAACCGAGAGGAGTACGCACGCTTCGACAGCGACGTGGAGGTGTACCGGGCGGTGACGCCGCTGGGGCCGCCTGCCGCCGAGTACTGGAACAGCCAGAAGGAAGTCCTGGAGAGGACCCGGGCGGAGTTGGACACGGTGTGCAGACACAACTACCAGTTGGAGCTCCGCACGACCTTGCAGCGGCGAG
            exon_3  nt_search and clustalo  TGGAGCCCACAGTGACCATCTCCCCATCCAGGACAGAGGCCCTCAACCACCACAACCTGCTGGTCTGCTCAGTGACAGATTTCTATCCAGCCCAGATCAAAGTCCGGTGGTTTCGGAATGACCAGGAGGAGACAACCGGCGTTGTGTCCACCCCCCTTATTAGGAACGGTGACTGGACCTTCCAGATCCTGGTGATGCTGGAAATGACTCCCCAGCATGGAGACGTCTACACCTGCCACGTGGAGCACCCCAGCCTCCAGAACCCCATCACCGTGGAGTGGC
            exon_1  nt_search and clustalo  AGAGACTCTCCCG
            exon_4  nt_search and clustalo  GGGCTCAGTCTGAATCTGCCCAGAGCAAGATG

        """

        # Check it the locus exists
        if not locus:
            if self.verbose:
                self.logger.info(self.logname + " No locus provided! ")

            # Guessing locus with blastn
            locus = get_locus(sequence, kir=self.kir, refdata=self.refdata)

            if locus and self.verbose:
                self.logger.info(self.logname + " Locus prediction = " + locus)

            if not locus:
                self.logger.error(self.logname
                                  + " Locus could not be determined!")
                # TODO: Raise exception
                return ''

        # Exact match found
        matched_annotation = self.refdata.search_refdata(sequence, locus)
        if matched_annotation:
            # TODO: return name of allele
            if self.verbose:
                self.logger.info(self.logname + " exact match found")

            # Create GFE
            feats, gfe = self.gfe.get_gfe(matched_annotation, locus)
            matched_annotation.gfe = gfe
            matched_annotation.exact = True
            return matched_annotation

        # Run blast to get ref sequences
        blast = blastn(sequence, locus, nseqs,
                       kir=self.kir, verbose=self.verbose,
                       refdata=self.refdata)

        # If the blastn fails..
        if blast.failed:
            self.logger.error(self.logname + " blastn failed! "
                              + locus + " " + str(sequence.seq))

            # Try and determine the locus and rerun. This is
            # useful for cases when the sequences is associated
            # with the wrong locus.
            locus = get_locus(sequence, kir=self.kir, refdata=self.refdata)

            if locus and self.verbose:
                self.logger.info(self.logname + " Locus prediction = " + locus)

            # Check if the locus could not be found
            if not locus:
                self.logger.error(self.logname
                                  + " Locus could not be determined!")
                # TODO: Raise exception
                return
            return self.annotate(sequence, locus)

        # Do seq_search first on all blast sequences
        partial_ann = None
        found = blast.match_seqs
        for i in range(0, len(found)-1):

            # Skip a reference
            # * For validation *
            if found[i].name in skip:
                continue

            if self.verbose:
                self.logger.info(self.logname
                                 + " Running seq_search with "
                                 + found[i].name)

            # * Running sequence search *
            # This does a simple string search for the
            # reference features within the provided sequence
            annotation = self.seqsearch.search_seqs(found[i],
                                                    sequence, locus,
                                                    partial_ann=partial_ann)
            if annotation.complete_annotation:
                # TODO: change clean to cleanup
                if self.verbose:
                    self.logger.info(self.logname
                                     + " Finished annotation using "
                                     + found[i].name)

                # Add alignment flag is specified
                if self.align:
                    if self.verbose:
                        self.logger.info(self.logname + " Adding alignment")
                    annotation = self.add_alignment(found[i], annotation)

                if self.verbose and self.verbosity > 0:
                    self.logger.info(self.logname
                                     + " Features annotated = "
                                     + ",".join(list(annotation
                                                     .annotation.keys())))
                    if self.verbosity > 2:
                        for f in annotation.features:
                            self.logger.info(self.logname
                                             + " " + f + " = "
                                             + str(annotation
                                                   .annotation[f].seq))

                # Create GFE
                feats, gfe = self.gfe.get_gfe(annotation, locus)
                annotation.gfe = gfe
                annotation.clean()
                return annotation
            else:
                partial_ann = annotation
                if self.verbose:
                    self.logger.info(self.logname
                                     + " Using partial annotation * run "
                                     + str(i) + " *")
                    self.logger.info(self.logname
                                     + " Features found = "
                                     + ",".join(list(annotation
                                                     .features.keys())))
                    self.logger.info(self.logname
                                     + " Sequence unmapped = "
                                     + str(annotation.covered))

        # The number of sequences being used for alignment
        # can't be greater than the number of sequences
        # to be returned from the blast results
        if alignseqs > len(found):
            alignseqs = len(found)-1

        # Now loop through doing alignment
        # TODO: Add parameter for limiting this step
        for i in range(0, alignseqs):

            # Skip a reference
            # * For validation *
            if found[i].name in skip:
                continue

            if self.verbose:
                self.logger.info(self.logname
                                 + " running ref_align with "
                                 + found[i].name)

            aligned_ann = self.ref_align(found[i], sequence, locus,
                                         annotation=partial_ann)
            if aligned_ann.complete_annotation:
                if self.align:
                    if self.verbose:
                        self.logger.info(self.logname + " Adding alignment")
                    annotation = self.add_alignment(found[i], annotation)

                if self.verbose:
                    self.logger.info(self.logname
                                     + " Finished ref_align annotation using "
                                     + found[i].name)

                if self.verbose and self.verbosity > 0:
                    self.logger.info(self.logname
                                     + " Features annotated = "
                                     + ",".join(list(annotation
                                                     .annotation.keys())))
                    if self.verbosity > 2:
                        for f in annotation.features:
                            self.logger.info(self.logname
                                             + " " + f + " = "
                                             + str(annotation
                                                   .annotation[f].seq))

                # Create GFE
                feats, gfe = self.gfe.get_gfe(annotation, locus)
                annotation.gfe = gfe
                annotation.clean()
                return aligned_ann
            else:
                if self.verbose:
                    self.logger.info(self.logname
                                     + " Using partial annotation "
                                     + "for alignment * run "
                                     + str(i) + " *")
                    self.logger.info(self.logname
                                     + " Features found = "
                                     + ",".join(list(annotation
                                                     .features.keys())))

                partial_ann = aligned_ann

        # TODO: make guess with full alignment
        # return self.ref_align(found, sequence, locus,
        #                       partial_ann=partial_ann)
        # TODO: raise exception
        self.logger.error(self.logname + " No annotation produced!")
        return ''

    def ref_align(self, found_seqs, sequence: Seq=None,
                  locus: str=None, annotation: Annotation=None,
                  partial_ann: Annotation=None) -> Annotation:
        """
        ref_align - Method for doing targeted alignments on partial annotations

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
        :rtype: Annotation

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
                    an, ins, dels = align_seqs(combseqr, feat, locus,
                                               verbose=self.align_verbose,
                                               verbosity=self.align_verbosity)
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

                                if self.verbose and self.verbosity > 2:
                                    sl = str(len(an.annotation[f]))
                                    self.logger.info(self.logname + " " + locus
                                                     + " " + f
                                                     + " len = " + sl
                                                     + " | max = "
                                                     + str(max_length)
                                                     + " | min = "
                                                     + str(min_length))

                                    if len(an.annotation[f]) <= max_length:
                                        self.logger.info(self.logname
                                                         + " " + locus
                                                         + " " + f
                                                         + " " + sl + " <= "
                                                         + str(max_length))
                                    else:
                                        self.logger.info(self.logname
                                                         + " " + locus
                                                         + " " + f
                                                         + " " + sl + " ! <= !"
                                                         + str(max_length))

                                    if len(an.annotation[f]) >= min_length:
                                        self.logger.info(self.logname
                                                         + " " + locus
                                                         + " " + f
                                                         + " " + sl + " >= "
                                                         + str(min_length))
                                    else:
                                        self.logger.info(self.logname + " "
                                                         + locus
                                                         + " " + f
                                                         + " " + sl + " ! >= !"
                                                         + str(min_length))

                                if(len(an.annotation[f]) <= max_length and
                                        len(an.annotation[f]) >= min_length):
                                    if self.verbose and self.verbosity > 0:
                                        self.logger.info(self.logname
                                                         + " Annotated " + f
                                                         + " with clustalo")
                                    annotation.annotation.update({f:
                                                                  an.annotation[f]
                                                                  })
                                    if an.blocks:
                                        mbtmp += an.blocks
                                    else:
                                        if self.verbose and self.verbosity > 0:
                                            self.logger.info(self.logname
                                                             + " All blocks mapped")
                                        if b in missing_blocks:
                                            del missing_blocks[missing_blocks.index(b)]
                                else:
                                    mbtmp.append(b)
                    else:
                        mbtmp.append(b)
                    annotation.blocks = mbtmp
                    annotation.check_annotation()
                    if annotation.complete_annotation:
                        if self.verbose:
                            self.logger.info(self.logname
                                             + " Completed annotation"
                                             + " with targeted ref_align")
                        return annotation

                if len(exons.seq) >= 4:
                    exonan, ins, dels = align_seqs(exons, feat, locus,
                                                   verbose=self.align_verbose,
                                                   verbosity=self.align_verbosity)
                    mapped_exons = list(exonan.annotation.keys())
                    if len(mapped_exons) >= 1:
                        if self.verbose:
                            self.logger.info(self.logname
                                             + " Annotated exons with align")
                        for f in exonan.annotation:
                            if self.verbose and self.verbosity > 0:
                                self.logger.info(self.logname
                                                 + " Annotated "
                                                 + f + " len = "
                                                 + str(len(exonan
                                                           .annotation[f])))
                            annotation.annotation.update({f: exonan.annotation[f]})
                        del missing_blocks[missing_blocks.index(b)]

                        annotation.blocks = missing_blocks
                        annotation.check_annotation()
                        if annotation.complete_annotation:
                            if self.verbose:
                                self.logger.info(self.logname + " Completed annotation with targeted exons ref_align")
                            return annotation

                # Run full sequence
                if len(fullrec.seq) >= 4:
                    fullref = align_seqs(fullrec, feat, locus,
                                         verbose=self.align_verbose,
                                         verbosity=self.align_verbosity)
                    if hasattr(fullref, 'annotation'):
                        mapped_full = list(fullref.annotation.keys())
                        if len(mapped_full) >= 1:

                            if self.verbose:
                                self.logger.info(self.logname + " Annotated fullrec with clustalo")

                            # If it wasn't found
                            del missing_blocks[missing_blocks.index(b)]

                            for f in fullref.annotation:
                                if self.verbose and self.verbosity > 0:
                                    self.logger.info(self.logname + " Annotated " + f + " len = " + str(len(exonan.annotation[f])))
                                annotation.update({f: fullref.annotation[f]})

                        annotation.blocks = missing_blocks
                        annotation.check_annotation()
                        if annotation.complete_annotation:
                            if self.verbose:
                                self.logger.info(self.logname + " Annotated all features with clustalo")
                            return annotation
            if self.verbose:
                self.logger.info(self.logname
                                 + " Failed to annotate all features")
            return annotation
        elif partial_ann:
            self.logger.error(self.logname + " Partial alignment! ")
            # Do full sequence alignments
            # any only extract out the part
            # that couldn't be explained from above
            return
        else:
            self.logger.error(self.logname + " full alginment")
            # Option for user that just want to
            # align a set of sequences
            return

    def add_alignment(self, ref_seq, annotation) -> Annotation:
        """
        add_alignment - method for adding the alignment to an annotation

        :param ref_seq: List of reference sequences
        :type ref_seq: List
        :param annotation: The complete annotation
        :type annotation: Annotation

        :rtype: Annotation
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
                        self.logger.info(self.logname + " Lengths match for " + feat)
                        self.logger.info(self.logname + " Gaps at " + feat)
                        self.logger.info(self.logname +
                                         "-".join([",".join([str(s)
                                                             for s in g])
                                                   for g in gaps]))
                    for i in range(0, len(gaps)):
                        for j in gaps[i]:
                            loc = j
                            seq.insert(loc, '-')
                    nseq = ''.join(seq)
                    annoated_align.update({feat: nseq})
                else:
                    in_seq = str(annotation.annotation[feat].seq)
                    ref_seq = self.refdata.annoated_alignments[locus][allele][feat]['Seq']
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



