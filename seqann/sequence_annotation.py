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
import logging
import warnings

from Bio import BiopythonExperimentalWarning
warnings.simplefilter("ignore", BiopythonExperimentalWarning)

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
from seqann.util import checkseq
from seqann.util import isexon
from seqann.util import isfive
from seqann.util import isutr
from seqann.util import is_classII
from seqann.gfe import GFE

from itertools import repeat
from typing import Dict
from typing import List
import re


class BioSeqAnn(Model):
    '''
        ::

            from seqann import BioSeqAnn
            seqann1 = BioSeqAnn()
            seqann2 = BioSeqAnn(dbversion="3300", verbose=True, verbosity=3)
            seqann3 = BioSeqAnn(debug={"align":4}, safe)

        :param server: A BioSQL database to use for retriving the sequence features. Using a BioSQL DB will speed up the annotations dramatically.
        :type server: :ref:`bio`
        :param dbversion: The IPD-IMGT/HLA or KIR database release.
        :type dbversion: ``str``
        :param datfile: The IPD-IMGT/HLA or KIR dat file to use in place of the server parameter.
        :type datfile: ``str``
        :param pid: A process label that can be provided to help track the logging output.
        :type pid: ``str``
        :param load_features: Flag for downloading all gene features and accessions from the feature service.
        :type load_features: ``bool``
        :param store_features: Flag for caching all features and their corresponding accessions.
        :type store_features: ``bool``
        :param cached_features: Dictionary containing all the features from the feature service.
        :type cached_features: ``dict``
        :param kir: Flag for indicating the input sequences are from the KIR gene system.
        :type kir: ``bool``
        :param align: Flag for producing the alignments along with the annotations.
        :type align: ``bool``
        :param verbose: Flag for running in verbose mode.
        :type verbose: ``bool``
        :param verbosity: Numerical value to indicate how verbose the output will be in verbose mode.
        :type verbosity: ``int``
        :param debug: Dictionary containing names of steps that you want to debug.
        :type debug: ``dict``
        :param safemode: Flag for running the annotations in safemode. No alignments will be done if no feature matches were made. This can prevent the alignment step for running for too long on bad sequences.
        :type safemode: ``bool``
    '''
    def __init__(self, server: BioSeqDatabase=None, dbversion: str='3310',
                 datfile: str='', verbose: bool=False, verbosity: int=0,
                 pid: str='NA', kir: bool=False, align: bool=False,
                 load_features: bool=False, store_features: bool=False,
                 refdata: ReferenceData=None,
                 cached_features: Dict=None,
                 safemode: bool=False,
                 debug: Dict=None):

        self.kir = kir
        self.align = align
        self.debug = debug
        self.server = server
        self.verbose = verbose
        self.safemode = safemode
        self.verbosity = verbosity
        self.align_verbose = verbose
        self.align_verbosity = verbosity
        self.dbversion = dbversion

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

        if refdata:
            # Initialize with preloaded reference object
            self.refdata = refdata
        else:
            # Create new reference data object
            self.refdata = ReferenceData(server=server,
                                         dbversion=dbversion,
                                         alignments=align,
                                         verbose=refdata_verbose,
                                         verbosity=refdata_verbosity,
                                         kir=kir)

        # Initalize SeqSearch for matching features
        self.seqsearch = SeqSearch(verbose=seqsearch_verbose,
                                   verbosity=seqsearch_verbosity)

    def annotate(self, sequence: Seq=None, locus: str=None,
                 nseqs: int=20, alignseqs: int=10,
                 skip: List=[],
                 rerun: bool=True,
                 full: bool=True) -> Annotation:
        """
        annotate - method for annotating a BioPython sequence

        :param sequence: The input consensus sequence.
        :type sequence: Seq
        :param locus: The gene locus associated with the sequence.
        :type locus: ``str``
        :param nseqs: The number of blast sequences to use.
        :type nseqs: ``int``
        :param alignseqs: The number of sequences to use for targeted alignments.
        :type alignseqs: ``int``
        :param skip: A list of alleles to skip for using as a reference. This is used for validation and testing.
        :type skip: ``List``
        :rtype: :ref:`ann`

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
            >>> from seqann import BioSeqAnn
            >>> sequence = Seq('AGAGACTCTCCCGAGGATTTCGTGTACCAGTTTAAGGCCATGTGCTACTTCACC')
            >>> seqann = BioSeqAnn()
            >>> ann = seqann.annotate(sequence)
            >>> for f in ann.annotation:
            ...    print(f, ann.method, str(ann.annotation[f].seq), sep="\t")
            exon_2  nt_search and clustalo  AGGATTTCGTGTACCAGTTTAAGGCCATGTGCTACTTCACCAACGGGACGGAGCGCGTGCGTTATGTGACCAGATACATCTATAACCGAGAGGAGTACGCACGCTTCGACAGCGACGTGGAGGTGTACCGGGCGGTGACGCCGCTGGGGCCGCCTGCCGCCGAGTACTGGAACAGCCAGAAGGAAGTCCTGGAGAGGACCCGGGCGGAGTTGGACACGGTGTGCAGACACAACTACCAGTTGGAGCTCCGCACGACCTTGCAGCGGCGAG
            exon_3  nt_search and clustalo  TGGAGCCCACAGTGACCATCTCCCCATCCAGGACAGAGGCCCTCAACCACCACAACCTGCTGGTCTGCTCAGTGACAGATTTCTATCCAGCCCAGATCAAAGTCCGGTGGTTTCGGAATGACCAGGAGGAGACAACCGGCGTTGTGTCCACCCCCCTTATTAGGAACGGTGACTGGACCTTCCAGATCCTGGTGATGCTGGAAATGACTCCCCAGCATGGAGACGTCTACACCTGCCACGTGGAGCACCCCAGCCTCCAGAACCCCATCACCGTGGAGTGGC
            exon_1  nt_search and clustalo  AGAGACTCTCCCG
            exon_4  nt_search and clustalo  GGGCTCAGTCTGAATCTGCCCAGAGCAAGATG

        """

        # If sequence is now a biopython
        # sequence record convert it to one
        if isinstance(sequence, Seq):
            sequence = SeqRecord(seq=sequence,
                                 id="NO_ID")

        if isinstance(sequence, str):
            sequence = SeqRecord(seq=Seq(sequence),
                                 id="NO_ID")

        # If sequence contains any characters
        # other than ATCG then the GFE notation
        # can not be created
        create_gfe = checkseq(sequence.seq)

        if self.verbose and not create_gfe:
            self.logger.warning(self.logname + " Sequence alphabet "
                                + "contains non DNA")
            self.logger.warning(self.logname
                                + " No GFE string will be generated")

        # Check it the locus exists
        if not locus:
            if self.verbose:
                self.logger.info(self.logname + " No locus provided! ")

            # Guessing locus with blastn
            locus = get_locus(sequence, kir=self.kir, refdata=self.refdata)

            if locus and self.verbose:
                self.logger.info(self.logname + " Locus prediction = " + locus)

            if not locus:
                if self.verbose:
                    self.logger.error(self.logname
                                      + " Locus could not be determined!")
                # TODO: Raise exception
                #raise NoLocusException("")
                return

        # Exact match found
        matched_annotation = self.refdata.search_refdata(sequence, locus)
        if matched_annotation and not skip:

            matched_annotation.exact = True

            # TODO: return name of allele
            if self.verbose:
                self.logger.info(self.logname + " exact match found")

            # Create GFE
            if create_gfe:
                feats, gfe = self.gfe.get_gfe(matched_annotation, locus)
                matched_annotation.gfe = gfe
                matched_annotation.structure = feats
            return matched_annotation

        # Run blast to get ref sequences
        blast = blastn(sequence, locus, nseqs,
                       kir=self.kir, verbose=self.verbose,
                       refdata=self.refdata)

        # If the blastn fails..
        if blast.failed:

            if self.verbose:
                self.logger.info(self.logname + " Blastn failed!!")

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
            return self.annotate(sequence, locus,
                                 nseqs=nseqs,
                                 alignseqs=alignseqs)
            return

        # Do seq_search first on all blast sequences
        # retain what the largest ref seq is
        leastmissing = 100
        partial_ann = None
        leastmissing_feat = None
        found = blast.match_seqs
        for i in range(0, len(found)):

            run = 0
            if i == len(found)-1:
                run = 1

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
            ann = self.seqsearch.search_seqs(found[i],
                                             sequence, locus,
                                             partial_ann=partial_ann,
                                             run=run)

            if ann.complete_annotation:
                if self.verbose:
                    self.logger.info(self.logname
                                     + " Finished annotation using "
                                     + found[i].name)

                # Add alignment flag is specified
                if self.align:
                    if self.verbose:
                        self.logger.info(self.logname + " Adding alignment")
                    ann = self.add_alignment(found[i], ann)

                if self.verbose and self.verbosity > 0:
                    self.logger.info(self.logname
                                     + " Features annotated = "
                                     + ",".join(list(ann
                                                     .annotation.keys())))
                    if self.verbosity > 3:
                        for f in ann.features:
                            self.logger.info(self.logname
                                             + " " + f + " = "
                                             + str(ann
                                                   .annotation[f].seq))

                # Create GFE
                if create_gfe:
                    feats, gfe = self.gfe.get_gfe(ann, locus)
                    ann.gfe = gfe
                    ann.structure = feats
                ann.clean()
                return ann
            else:
                partial_ann = ann

                if hasattr(partial_ann, 'refmissing'):
                    if len(partial_ann.refmissing) < leastmissing:
                        leastmissing_feat = found[i]
                        leastmissing = len(partial_ann.refmissing)
                else:
                    leastmissing_feat = found[i]
                    leastmissing = 0

                if self.verbose and self.verbosity > 0:
                    self.logger.info(self.logname
                                     + " Using partial annotation * run "
                                     + str(i) + " *")
                    self.logger.info(self.logname
                                     + " Features found ("
                                     + str(len(ann.features.keys())) + ") = "
                                     + ",".join(list(ann
                                                     .features.keys())))
                    self.logger.info(self.logname
                                     + " Features missing ("
                                     + str(len(ann.missing.keys())) + ") = "
                                     + ",".join(list(ann
                                                     .missing.keys())))
                    self.logger.info(self.logname
                                     + " Sequence unmapped = "
                                     + str(ann.covered))
                    self.logger.info(self.logname + " ############" +
                                     "##################")

        # The number of sequences being used for alignment
        # can't be greater than the number of sequences
        # to be returned from the blast results
        if alignseqs > len(found):
            alignseqs = len(found)-1

        # * HARD CODED LOGIC * #
        # > After testing with multiple thresholds
        #   this value seemed to work best.
        #
        # Aligned % cutoff
        align_cutoff = .90
        if((not hasattr(partial_ann, 'features') or
           len(partial_ann.features) == 0)
           and len(sequence) > 700 and self.safemode):
            self.logger.error("No feature matches!")
            self.logger.error("Running in safe mode. " +
                              "No alignments will be done!")

            if rerun:
                # Check to see if reverse comp
                # TODO: Add note for reverse complement
                self.logger.info("Running with reverse complement.")
                sequence = sequence.reverse_complement()
                return self.annotate(sequence=sequence,
                                     locus=locus,
                                     rerun=False)
            return

        # Now loop through doing alignment
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
                                         annotation=partial_ann,
                                         run=i,
                                         cutoff=align_cutoff)

            # * HARD CODED LOGIC * #
            # > If sequences are very novel, then the alignment
            #   cutoff may be to stringent. Incrementally decreasing
            #   the cutoff improves the likelihood of these sequences
            #   being annotated.
            align_cutoff -= .01
            if aligned_ann and aligned_ann.complete_annotation:
                if self.align:
                    if self.verbose:
                        self.logger.info(self.logname + " Adding alignment")
                    aligned_ann = self.add_alignment(found[i], aligned_ann)

                if self.verbose:
                    self.logger.info(self.logname
                                     + " Finished ref_align annotation using "
                                     + found[i].name)

                if self.verbose and self.verbosity > 0:
                    self.logger.info(self.logname
                                     + " Features annotated = "
                                     + ",".join(list(aligned_ann
                                                     .annotation.keys())))
                    if self.verbosity > 2:
                        for f in aligned_ann.annotation:
                            self.logger.info(self.logname
                                             + " " + f + " = "
                                             + str(aligned_ann
                                                   .annotation[f].seq))

                # Create GFE
                if create_gfe:
                    feats, gfe = self.gfe.get_gfe(aligned_ann, locus)
                    aligned_ann.gfe = gfe
                    aligned_ann.structure = feats
                aligned_ann.clean()
                return aligned_ann
            elif(aligned_ann):

                if self.verbose and self.verbosity > 0:
                    self.logger.info(self.logname
                                     + " Using partial annotation "
                                     + "for alignment * run "
                                     + str(i) + " - cutoff = "
                                     + str(align_cutoff)
                                     + " *")
                    self.logger.info(self.logname
                                     + " Features found = "
                                     + ",".join(list(aligned_ann
                                                     .features.keys())))
                    self.logger.info(self.logname
                                     + " Features missing = "
                                     + ",".join(list(aligned_ann
                                                     .missing.keys())))
                    self.logger.info(self.logname + " ############" +
                                     "##################")
                partial_ann = aligned_ann

                if(hasattr(partial_ann, 'annotation')
                   and partial_ann.annotation):
                    exon_only = True
                    for f in partial_ann.annotation:
                        if re.search("intron", f) or re.search("UTR", f):
                            exon_only = False

                    if(is_classII(locus) and exon_only
                       and len(partial_ann.annotation.keys()) > 0
                       and align_cutoff < .9):
                        align_cutoff = .80

            if not is_classII(locus) and align_cutoff < .88:
                align_cutoff = .88

        # Don't run full
        # annotation if flag is passed
        if not full:
            return

        if self.verbose:
            self.logger.info(self.logname + " running full alignment")

        # Try doing full alignment
        full_align = self.ref_align(leastmissing_feat,
                                    sequence,
                                    locus,
                                    partial_ann=partial_ann,
                                    cutoff=.80)

        if self.verbose:
            self.logger.info(self.logname
                             + " Finished ref_align annotation using full "
                             + leastmissing_feat.name)

        # Check to see if an annotation was returned
        if(not isinstance(full_align, Annotation)
           or isinstance(full_align, str)):
            if(not rerun or len(sequence) > 4000):
                self.logger.info(self.logname + " Failed annotation!")
                return Annotation()
            else:
                if self.verbose and self.verbosity > 0:
                    self.logger.info(self.logname
                                     + " Reruning annotation!")
                return self.annotate(sequence=sequence,
                                     locus=locus,
                                     alignseqs=2,
                                     nseqs=nseqs+1,
                                     skip=[found[0].name],
                                     rerun=False)

        # Check if the annotation is complete
        if not full_align.complete_annotation and self.verbose:
            self.logger.info(self.logname + " Incomplete annotation!")

        # Add the alignment to the annotation
        if self.align and full_align.complete_annotation:
            if self.verbose:
                self.logger.info(self.logname + " Adding alignment")
            full_align = self.add_alignment(leastmissing_feat, full_align)

        if self.verbose and self.verbosity > 0:
            self.logger.info(self.logname
                             + " Features annotated = "
                             + ",".join(list(full_align
                                             .annotation.keys())))
            if self.verbosity > 2:
                for f in full_align.annotation:
                    self.logger.info(self.logname
                                     + " " + f + " = "
                                     + str(full_align
                                           .annotation[f].seq))
        # Create GFE
        if create_gfe and full_align.complete_annotation:
            feats, gfe = self.gfe.get_gfe(full_align, locus)
            full_align.gfe = gfe
            full_align.structure = feats

        full_align.clean()
        if(full_align.complete_annotation
           or not rerun or len(sequence) > 4000):
            return full_align
        else:
            if self.verbose and self.verbosity > 0:
                self.logger.info(self.logname
                                 + " Reruning annotation!")
            return self.annotate(sequence=sequence,
                                 locus=locus,
                                 nseqs=nseqs+1,
                                 alignseqs=2,
                                 skip=[found[0].name],
                                 rerun=False)

    def ref_align(self, found_seqs, sequence: Seq=None,
                  locus: str=None, annotation: Annotation=None,
                  partial_ann: Annotation=None,
                  run: int=0,
                  cutoff: float=.90) -> Annotation:
        """
        ref_align - Method for doing targeted alignments on partial annotations

        :param found_seqs: The input sequence record.
        :type found_seqs: Seq
        :param sequence: The input sequence record.
        :type sequence: Seq
        :param locus: The gene locus associated with the sequence.
        :type locus: ``str``
        :param annotation: The incomplete annotation from a previous iteration.
        :type annotation: :ref:`ann`
        :param partial_ann: The partial annotation after looping through all of the blast sequences.
        :type partial_ann: :ref:`ann`
        :rtype: :ref:`ann`

        """
        if annotation and isinstance(annotation, Annotation):
            if 0 in annotation.mapping \
                    and not isinstance(annotation.mapping[0], int):
                ft = annotation.mapping[0]
                start_order = self.refdata.structures[locus][ft]
            else:
                start_order = 0

            # Check whether it's exons only
            exon_only = True
            if hasattr(annotation, 'annotation') and annotation.annotation:
                for f in annotation.annotation:
                    if re.search("intron", f) or re.search("UTR", f):
                        exon_only = False
            elif(len(sequence.seq) > 900):
                exon_only = False

            annoated = []
            if hasattr(annotation, 'annotation') and annotation.annotation:
                annoated = list(annotation.annotation.keys())

            # Extract the missing blocks and
            # only align those blocks to the known
            # missing features
            # Start with all blocks missing
            # and then delete block if it is found
            tmp_missing = []
            missing_blocks = annotation.blocks
            for b in sorted(annotation.blocks):
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

                # Print out different sequence types being align
                if self.verbose and self.verbosity > 3:
                    for combseqr in combosrecs:
                        self.logger.info(self.logname + " " + combseqr.id)

                ic = 0

                # Looping through the generated combos
                # from the reference sequence. (ex. intron1-exon1)
                for combseqr in combosrecs:
                    if self.verbose:
                        self.logger.info(self.logname + " aligning ->" + combseqr.id)
                    
                    # Running the alignment
                    an, ins, dels = align_seqs(combseqr, feat, locus,
                                               start,
                                               annotation.missing,
                                               len(annoated),
                                               cutoff=cutoff,
                                               verbose=self.align_verbose,
                                               verbosity=self.align_verbosity)

                    # Checking if any features were mapped
                    mapped_feat = list(an.annotation.keys())
                    if len(mapped_feat) >= 1:

                        # loop through the annotated features
                        for f in an.annotation:
                            f_order = self.refdata.structures[locus][f]

                            # Only map features if they are in the correct order
                            if f in annotation.missing \
                                    and f_order >= start_order:

                                # * HARD CODED LOGIC *
                                # Only accept the alignments if they are
                                # somewhat similar to what's been observed already
                                length, lengthsd = 0, 0
                                length = float(self.refdata.feature_lengths[locus][f][0])
                                lengthsd = float(self.refdata.feature_lengths[locus][f][1])

                                incr = 3 if not is_classII(locus) else 4

                                max_length = length + (lengthsd*incr) + ins
                                min_length = length - (lengthsd*incr) - dels

                                if f == "exon_8" and not is_classII(f):
                                    max_length = 10

                                # Check ordering when it's only exon sequences
                                if exon_only:
                                    f_order = self.refdata.structures[locus][f]
                                    endp = an.features[f].location.end + 1

                                    #Make sure order of alignment make sense
                                    if an.features[f].location.start == 0 \
                                            and f != "five_prime_UTR" \
                                            and not isexon(f):
                                        del an.features[f]
                                        continue

                                    if endp in annotation.mapping and not isinstance(annotation.mapping[endp], int):
                                        mf = annotation.mapping[endp]
                                        expected_order = f_order + 1
                                        expected_order2 = f_order + 2
                                        if expected_order != self.refdata.structures[locus][mf] and expected_order2 != self.refdata.structures[locus][mf]:
                                            self.logger.info(self.logname + " out of order1 -> " + mf + " " + f)                                        
                                            del an.features[f]
                                            continue

                                    startp = an.features[f].location.start - 1
                                    if startp in annotation.mapping and not isinstance(annotation.mapping[startp], int):
                                        mf = annotation.mapping[startp]
                                        expected_order = f_order - 1
                                        expected_order2 = f_order - 2
                                        if expected_order != self.refdata.structures[locus][mf] and expected_order2 != self.refdata.structures[locus][mf]:
                                            self.logger.info(self.logname + " out of order2 -> " + mf + " " + f)
                                            del an.features[f]
                                            continue
                                else:

                                    ##Make sure order of alignment make sense
                                    if an.features[f].location.start == 0 \
                                            and f != "five_prime_UTR" \
                                            and 'three_prime_UTR' in annotation.annotation:
                                        del an.features[f]
                                        continue

                                    f_order = self.refdata.structures[locus][f]
                                    endp = an.features[f].location.end + 1
                                    if endp in annotation.mapping and not isinstance(annotation.mapping[endp], int):
                                        mf = annotation.mapping[endp]
                                        expected_order = f_order + 1
                                        if expected_order != self.refdata.structures[locus][mf]:
                                            self.logger.info(self.logname + " out of order12 -> " + mf + " " + f)
                                            del an.features[f]
                                            continue

                                    startp = an.features[f].location.start - 1
                                    if startp in annotation.mapping and not isinstance(annotation.mapping[startp], int):
                                        mf = annotation.mapping[startp]
                                        expected_order = f_order - 1
                                        if expected_order != self.refdata.structures[locus][mf]:
                                            self.logger.info(self.logname + " out of order22 -> " + mf + " " + f)
                                            del an.features[f]
                                            continue

                                if self.verbose and self.verbosity > 0:
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

                                # Update the original annotation
                                # with the returned annotation
                                if(len(an.annotation[f]) <= max_length and
                                        len(an.annotation[f]) >= min_length):
                                    if self.verbose and self.verbosity > 0:
                                        self.logger.info(self.logname
                                                         + " Annotated " + f
                                                         + " with clustalo using " +
                                                         combseqr.id)
                                        self.logger.info(self.logname
                                                         + " Coordinates for " + f
                                                         + str(an.features[f].location.start)
                                                         + " - " + str(an.features[f].location.end)
                                                         )
                                    if annotation.annotation:
                                        annotation.annotation.update({f:
                                                                    an.annotation[f]
                                                                  })
                                        annotation.features.update({f:
                                                                    an.features[f]
                                                                  })  
                                    else:
                                        annotation.annotation = {}
                                        annotation.annotation.update({f:
                                                                    an.annotation[f]
                                                                  })
                                        annotation.features.update({f:
                                                                    an.features[f]
                                                                  })
                                    if f in annotation.refmissing:
                                        i = annotation.refmissing.index(f)
                                        del annotation.refmissing[i]

                                    if f in annotation.missing:
                                        del annotation.missing[f]

                                    #if b in annotation.blocks:
                                    #    del annotation.blocks[annotation.blocks.index(b)]
                                    #     print(annotation.blocks)
                                    # if an.blocks:
                                    #     #print("PARTIAL BLOCK")
                                    #     #print(an.blocks)
                                    #     if b in missing_blocks:
                                    #         del missing_blocks[missing_blocks.index(b)]

                                    #     if self.verbose and self.verbosity > 0:
                                    #         self.logger.info(self.logname
                                    #                          + " Part of block mapped")

                                    # else:
                                    #     if self.verbose and self.verbosity > 0:
                                    #         self.logger.info(self.logname
                                    #                          + " All blocks mapped")
                                        
                                    #     del annotation.blocks[annotation.blocks.index(b)]
                                    #     if b in missing_blocks:
                                    #         del missing_blocks[missing_blocks.index(b)]
                                else:
                                    self.logger.info(self.logname + " FAILED HERE!!!!!!!!!")
                                #elif b not in mbtmp and b in missing_blocks:
                                #    mbtmp.append(b)
                            else:
                                self.logger.info(self.logname + " OUT OF ORDER !!!!!!!!!!!!!!!!!!!")

                    # Update the coordinates
                    coordinates = dict(map(lambda x: [x, 1], [i for i in range(0, len(sequence.seq)+1)]))
                    for f in annotation.features:
                        s = annotation.features[f].location.start
                        e = annotation.features[f].location.end
                        if s != 0:
                            s += 1
                            e += 1
                        else:
                            e += 1

                        for i in range(s, e):
                            annotation.mapping[i] = f
                            if i in coordinates:
                                del coordinates[i]

                    # Get any remaining blocks after updated
                    # annotation
                    blocks = getblocks(coordinates)
                    annotation.blocks = blocks

                    # Check to see if that annotation is complete
                    annotation.check_annotation()
                    if annotation.complete_annotation:
                        if self.verbose:
                            self.logger.info(self.logname
                                             + " Completed annotation"
                                             + " with targeted ref_align")
                        return annotation
                    else:

                        if an.features:
                            # for f in an.features:
                            #     f_order = self.refdata.structures[locus][f]
                            #     # Only add features if they are after the
                            #     # first feature mapped
                            #     if f_order >= start_order and f not in annotation.features \
                            #             and f in annotation.annotation:
                            #         annotation.features[f] = an.features[f]

                            # Rerunning seqsearch with
                            # new annotation from alignment
                            tmpann = self.seqsearch.search_seqs(found_seqs,
                                                                sequence,
                                                                locus,
                                                                partial_ann=annotation,
                                                                run=run)

                            if tmpann.complete_annotation:
                                for f in tmpann.annotation:
                                    if f not in annotation.annotation:
                                        annotation.annotation[f] = tmpann.annotation[f]
                                if self.verbose:
                                    self.logger.info(self.logname
                                                     + " Completed annotation"
                                                     + " with targeted ref_align and seqsearch!")

                                return tmpann
                            annotation = tmpann

                    ic += 1

                # Has to be missing exons
                exons_n = 0
                for f in annotation.missing:
                    if re.search("intron", f) or re.search("UTR", f):
                        exons_n += 1

                # Run exon only alignment
                if len(exons.seq) >= 4 and exons_n > 0:
                    exonan, ins, dels = align_seqs(exons, feat, locus, start,
                                                   annotation.missing,
                                                   len(annoated),
                                                   cutoff=cutoff,
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
                            annotation.features.update({f: exonan.features[f]})

                        coordinates = dict(map(lambda x: [x, 1], [i for i in range(0, len(sequence.seq)+1)]))
                        for f in annotation.features:
                            s = annotation.features[f].location.start
                            e = annotation.features[f].location.end
                            if s != 0:
                                s += 1
                                e += 1
                            else:
                                e += 1
                            for i in range(s, e):
                                annotation.mapping[i] = f
                                if i in coordinates:
                                    del coordinates[i]

                        blocks = getblocks(coordinates)
                        annotation.blocks = blocks
                        annotation.check_annotation()
                        if annotation.complete_annotation:
                            if self.verbose:
                                self.logger.info(self.logname + " Completed annotation with targeted exons ref_align")
                            return annotation

            return annotation
        elif partial_ann:

            annoated = []
            if hasattr(partial_ann, 'annotation') and partial_ann.annotation:
                annoated = list(partial_ann.annotation.keys())

            # Do full sequence alignments
            # any only extract out the part
            # that couldn't be explained from above
            if 0 in partial_ann.mapping \
                    and not isinstance(partial_ann.mapping[0], int):
                ft = partial_ann.mapping[0]
                start_order = self.refdata.structures[locus][ft]
            else:
                start_order = 0

            # Extract the missing blocks and
            # only align those blocks to the known
            # missing features
            # Start with all blocks missing
            # and then delete block if it is found
            tmp_missing = []
            missing_blocks = partial_ann.blocks
            for b in sorted(partial_ann.blocks):

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

                feat = seq_feat.extract(partial_ann.seq)
                combosrecs, exons, fullrec = self._refseqs(locus,
                                                           start,
                                                           partial_ann,
                                                           feat,
                                                           b)
                if len(fullrec.seq) >= 4:

                    fullref, ins, dels = align_seqs(fullrec, feat,
                                                    locus, start,
                                                    partial_ann.missing,
                                                    len(annoated),
                                                    cutoff=cutoff,
                                                    verbose=self.align_verbose,
                                                    verbosity=self.align_verbosity)

                    if hasattr(fullref, 'features') and fullref.features:
                        mapped_full = list(fullref.annotation.keys())
                        if len(mapped_full) >= 1:

                            if self.verbose:
                                self.logger.info(self.logname
                                                 + " Annotated fullrec"
                                                 + " with clustalo")

                            # If it wasn't found
                            del missing_blocks[missing_blocks.index(b)]

                            for f in fullref.annotation:
                                if self.verbose and self.verbosity > 0:
                                    self.logger.info(self.logname + " Annotated " + f + " len = " + str(len(fullref.annotation[f])))
                                partial_ann.annotation.update({f: fullref.annotation[f]})

                        if b in missing_blocks:
                            del missing_blocks[missing_blocks.index(b)]
                        else:
                            for bm in tmp_missing:
                                if bm in missing_blocks:
                                    del missing_blocks[missing_blocks.index(bm)]

                        for f in fullref.features:
                            f_order = self.refdata.structures[locus][f]
                            # Only add features if they are after the
                            # first feature mapped
                            if f_order >= start_order and f not in partial_ann.features \
                                    and f in partial_ann.annotation:
                                partial_ann.features[f] = fullref.features[f]

                        coordinates = dict(map(lambda x: [x, 1], [i for i in range(0, len(sequence.seq)+1)]))
                        for f in partial_ann.features:
                            s = partial_ann.features[f].location.start
                            e = partial_ann.features[f].location.end
                            if s != 0:
                                s += 1
                                e += 1
                            else:
                                e += 1
                            for i in range(s, e):
                                partial_ann.mapping[i] = f
                                if i in coordinates:
                                    del coordinates[i]

                        blocks = getblocks(coordinates)
                        partial_ann.check_annotation()
                        if partial_ann.complete_annotation:
                            if self.verbose:
                                self.logger.info(self.logname + " Annotated all features with clustalo")
                        
                    return partial_ann

            if self.verbose:
                self.logger.info(self.logname
                                 + " Failed to annotate features")
            return ''
        # else:
        #     self.logger.error(self.logname + " invalid use")
        #     # Option for user that just want to
        #     # align a set of sequences
        #     return

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

        # TODO: Do not use combos that
        # don't have missing features
        if start_pos == 0:
            prv_order = -1
        else:
            prev_feat = mapping[start_pos-1]
            if isinstance(prev_feat, int):
                prev_feat = mapping[start_pos]
            prv_order = self.refdata.structures[locus][prev_feat]

        if len(mapping) == end_pos+1:
            nxt_order = len(self.refdata.structures[locus])+1
        else:
            next_feat = mapping[end_pos+1]
            if isinstance(next_feat, int):
                nxt_order = prv_order + 2
            else:
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
            if mis_order > nxt_order \
                    or mis_order < prv_order:
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

                #length = float(self.refdata.feature_lengths[locus][f][0])
                #lengthsd = float(self.refdata.feature_lengths[locus][f][1])

                #max_length = length + lengthsd
                #min_length = length - lengthsd
                #if len(rec.seq) <= max_length \
                #       and len(rec.seq) >= min_length:
                recf, x = self._make_seqfeat(0, rec.seq, nfeat)
                ctmpid1 = "Ref1_" + str(nfeat) + "_" + str(randomid(N=2))
                tmprec1 = SeqRecord(seq=rec.seq, features=[recf],
                                    id=ctmpid1)
                combos.append(tmprec1)
            else:
                if isexon(f):

                    cstart = 0
                    ref_feats = []
                    combo_seq = []
                    feat_names = []
                    combo_feats = []
                    for i in range(start_ord, nxt_order):
                        nfeat = self.refdata.struct_order[locus][i]
                        if nfeat not in all_seqrecs or not isexon(nfeat):
                            continue

                        ref_feats.append(nfeat)
                        length, lengthsd = 0, 0

                        #for fr in ref_feats:
                        #    length += float(self.refdata.feature_lengths[locus][fr][0])
                        #    lengthsd += float(self.refdata.feature_lengths[locus][fr][1])

                        # min and max lengths expected
                        #max_length = length + lengthsd
                        #min_length = length - lengthsd

                        rec = all_seqrecs[nfeat]
                        feat_names.append(nfeat)

                        cfeat, cstart = self._make_seqfeat(cstart, rec.seq, nfeat)
                        combo_feats.append(cfeat)
                        combo_seq.append(rec.seq)

                        seqtmp = "".join([str(ft) for ft in combo_seq])
                        ctmpseq = Seq(seqtmp, IUPAC.unambiguous_dna)
                        ctmpid = "|".join([str(ft) for ft in feat_names]) \
                                 + "_" + str(randomid(N=2))

                        #if len(ctmpseq) <= max_length \
                        #      and len(ctmpseq) >= min_length:
                        tmprec = SeqRecord(seq=ctmpseq, features=combo_feats,
                                           id=ctmpid)
                        combos.append(tmprec)

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

                    #for fr in ref_feats:
                    #    length += float(self.refdata.feature_lengths[locus][fr][0])
                    #    lengthsd += float(self.refdata.feature_lengths[locus][fr][1])

                    # min and max lengths expected
                    #max_length = length + lengthsd
                    #min_length = length - lengthsd

                    rec = all_seqrecs[nfeat]
                    feat_names.append(nfeat)

                    cfeat, cstart = self._make_seqfeat(cstart, rec.seq, nfeat)
                    combo_feats.append(cfeat)
                    combo_seq.append(rec.seq)

                    seqtmp = "".join([str(ft) for ft in combo_seq])
                    ctmpseq = Seq(seqtmp, IUPAC.unambiguous_dna)
                    ctmpid = "|".join([str(ft) for ft in feat_names]) \
                             + "_" + str(randomid(N=2))

                    #if len(ctmpseq) <= max_length \
                    #      and len(ctmpseq) >= min_length:
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
            ftype = "5UTR" if isfive(featname) else "3UTR"
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


# TODO: Add to util
def getblocks(coords):
    block = []
    blocks = []
    sorted_i = sorted(coords.keys())
    if len(sorted_i) == 1:
        return [sorted_i]
    for i in range(0, len(sorted_i)-1):
        j = i+1
        if i == 0:
            block.append(sorted_i[i])
        if(j <= len(sorted_i)-1):
            if(sorted_i[i] == sorted_i[j]-1):
                block.append(sorted_i[j])
            else:
                blocks.append(block)
                block = []
                block.append(sorted_i[j])
        else:
            block.append(sorted_i[j])
    if len(block) > 1:
        blocks.append(block)
    return blocks

