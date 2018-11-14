# -*- coding: utf-8 -*-

#
#    seqann Sequence Annotation.
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
import logging
import warnings

from Bio import BiopythonExperimentalWarning
warnings.simplefilter("ignore", BiopythonExperimentalWarning)

from Bio.SeqUtils import nt_search
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import ExactPosition

from seqann.models.annotation import Annotation
from seqann.models.base_model_ import Model
from seqann.util import deserialize_model
from seqann.util import get_features
from seqann.util import get_structures
from seqann.util import get_structorder
from seqann.util import get_structmax
from seqann.util import is_classII


# NOTE: This doesn't need to be a class anymore and could be reduced to functions 

class SeqSearch(Model):
    '''
        This is a class for annotating a BioPython sequence without using alignment

        :param verbose: Flag for running in verbose mode.
        :type verbose: ``bool``
        :param verbosity: Numerical value to indicate how verbose the output will be in verbose mode.
        :type verbosity: ``int``

        Example usage:

            >>> from seqann.seq_search import SeqSearch
            >>> seqsrch = SeqSearch()


    '''
    def __init__(self, verbose: bool=False, verbosity: int=0):
        self.data_types = {
            'verbose': bool,
            'verbosity': int
        }

        self.attribute_map = {
            'verbose': 'verbose',
            'verbosity': 'verbosity'
        }

        self._verbose = verbose
        self._verbosity = verbosity
        self.logger = logging.getLogger("Logger." + __name__)

    @classmethod
    def from_dict(cls, dikt) -> 'SeqSearch':
        """
        Returns the dict as a model

        :param dikt: A dict.
        :type: dict
        :return: The SeqSearch of this SeqSearch.
        :rtype: SeqSearch
        """
        return deserialize_model(dikt, cls)

    def search_seqs(self, seqrec, in_seq, locus, run=0, partial_ann=None):
        """
        search_seqs - method for annotating a BioPython sequence without alignment

        :param seqrec: The reference sequence
        :type seqrec: SeqRecord
        :param locus: The gene locus associated with the sequence.
        :type locus: str
        :param in_seq: The input sequence
        :type in_seq: SeqRecord
        :param run: The number of runs that have been done
        :type run: int
        :param partial_ann: A partial annotation from a previous step
        :type partial_ann: :ref:`ann`
        :rtype: :ref:`ann`

        Example usage:

            >>> from Bio.Seq import Seq
            >>> from seqann.seq_search import SeqSearch
            >>> inseq = Seq('AGAGACTCTCCCGAGGATTTCGTGTACCAGTTTAAGGCCATGTGCTACTTCACC')
            >>> sqsrch = SeqSearch()
            >>> ann = sqsrch.search_seqs(refseqs, inseq)

        """
        # Extract out the sequences and feature names
        # from the reference sequences

        # The mapped features will be subtracted from seq_covered
        # so the final seq_covered number will reflect the remaining
        # number of base pairs that haven't been mapped.
        #
        # The coordinates and mapping will help determine what positions
        # in the sequence have been mapped and to what features. The
        # missing blocks variable will be generated using these.
        structures = get_structures()
        seq_covered = len(in_seq.seq)
        coordinates = dict(map(lambda x: [x, 1],
                           [i for i in range(0, len(in_seq.seq)+1)]))

        mapping = dict(map(lambda x: [x, 1],
                           [i for i in range(0, len(in_seq.seq)+1)]))

        ambig_map = {}
        found_feats = {}
        feat_missing = {}

        method = "nt_search" if not partial_ann else partial_ann.method

        # If the partial annotation is provided
        # then make the found_feats equal to
        # what has already been annotated
        feats = get_features(seqrec)
        if partial_ann:

            found_feats = partial_ann.features

            if self.verbose and self.verbosity > 4:
                self.logger.info("Found partial features:")
                for f in found_feats:
                    self.logger.info(f)

            # Skip references that only have features
            # that have already been annoated
            if len([f for f in feats if f in found_feats]) == len(feats):
                if self.verbose:
                    self.logger.info("Skipping incomplete refseq")
                return partial_ann

            if self.verbose and self.verbosity > 1:
                self.logger.info("Using partial annotation | "
                                 + locus + " "
                                 + str(len(partial_ann.features)))

            coordinates = dict(map(lambda l: [l, 1],
                                   [item for sublist
                                    in partial_ann.blocks
                                    for item in sublist]))
            seq_covered = partial_ann.covered
            mapping = partial_ann.mapping

            if self.verbose and self.verbosity > 2:
                self.logger.info("Partial sequence coverage = "
                                 + str(seq_covered))
                self.logger.info("Partial sequence metho = "
                                 + method)

        added_feat = {}
        deleted_coords = {}
        for feat_name in sorted(feats,
                                key=lambda k: structures[locus][k]):

            # skip if partial annotation is provided
            # and the feat name is not one of the
            # missing features
            if partial_ann and feat_name not in partial_ann.refmissing:
                if self.verbose and self.verbosity > 1:
                    self.logger.info("Skipping " + feat_name
                                     + " - Already annotated")
                continue

            if self.verbose and self.verbosity > 1:
                self.logger.info("Running seqsearch for " + feat_name)

            # Search for the reference feature sequence in the
            # input sequence. Record the coordinates if it's
            # found and if it's found in multiple spots. If it
            # is not found, then record that feature as missing.
            seq_search = nt_search(str(in_seq.seq), str(feats[feat_name]))

            if len(seq_search) == 2:

                if self.verbose and self.verbosity > 0:
                    self.logger.info("Found exact match for " + feat_name)

                seq_covered -= len(str(feats[feat_name]))
                end = int(len(str(feats[feat_name])) + seq_search[1])

                if feat_name == 'three_prime_UTR' \
                        and len(str(in_seq.seq)) > end:
                        end = len(str(in_seq.seq))

                # If the feature is found and it's a five_prime_UTR then
                # the start should always be 0, so insertions at the
                # beinging of the sequence will be found.
                start = seq_search[1] if feat_name != 'five_prime_UTR' else 0
                si = seq_search[1]+1 if seq_search[1] != 0 and \
                    feat_name != 'five_prime_UTR' else 0

                # check if this features has already been mapped
                mapcheck = set([0 if i in coordinates else 1
                               for i in range(si, end+1)])

                # Dont map features if they are out of order
                skip = False
                if found_feats and len(found_feats) > 0:
                    for f in found_feats:
                        o1 = structures[locus][feat_name]
                        o2 = structures[locus][f]
                        loctyp = loctype(found_feats[f].location.start,
                                         found_feats[f].location.end,
                                         start, end)

                        if o1 < o2 and loctyp:
                            skip = True
                            if self.verbose:
                                self.logger.info("Skipping map for "
                                                 + feat_name)
                        elif o2 < o1 and not loctyp:
                            skip = True
                            if self.verbose:
                                self.logger.info("Skipping map for "
                                                 + feat_name)

                if 1 not in mapcheck and not skip:
                    for i in range(si, end+1):
                        if i in coordinates:
                            if feat_name == "exon_8" or feat_name == 'three_prime_UTR':
                                deleted_coords.update({i: coordinates[i]})
                            del coordinates[i]
                        else:
                            if self.verbose:
                                self.logger.error("seqsearch - should't be here "
                                                  + locus + " - "
                                                  + " - " + feat_name)
                        mapping[i] = feat_name

                    found_feats.update({feat_name:
                                        SeqFeature(
                                            FeatureLocation(
                                                ExactPosition(start),
                                                ExactPosition(end), strand=1),
                                            type=feat_name)})

                    if feat_name == "exon_8" or feat_name == 'three_prime_UTR':
                        added_feat.update({feat_name: feats[feat_name]})
                    if self.verbose and self.verbosity > 3:
                        self.logger.info("Coordinates | Start = " + str(start) + " - End = " + str(end))

            elif(len(seq_search) > 2):
                if self.verbose and self.verbosity > 1:
                    self.logger.info("Found " + str(len(seq_search))
                                     + " matches for " + feat_name)

                new_seq = [seq_search[0]]
                for i in range(1, len(seq_search)):
                    tnp = seq_search[i]+1
                    if seq_search[i] in coordinates or tnp in coordinates:
                        new_seq.append(seq_search[i])

                seq_search = new_seq
                if(partial_ann and feat_name == "exon_8" and run > 0):
                    missing_feats = sorted(list(partial_ann.missing.keys()))

                    # * HARD CODED LOGIC * #
                    # > exon8 in class I maps to multiple spots in a sequence,
                    #   often in the 3' UTR. These features need to be mapped
                    #   last to make sure it's not mapping exon8 incorrectly.
                    if(missing_feats == ['exon_8', 'three_prime_UTR']
                       and len(seq_search) <= 3):
                        if self.verbose and self.verbosity > 0:
                            self.logger.info("Resolving exon_8")

                        seq_covered -= len(str(feats[feat_name]))
                        end = int(len(str(feats[feat_name])) + seq_search[1])

                        # If the feature is found and it's a five_prime_UTR then
                        # the start should always be 0, so insertions at the
                        # beinging of the sequence will be found.
                        start = seq_search[1]
                        si = seq_search[1]+1 if seq_search[1] != 0 else 0

                        # check if this features has already been mapped
                        mapcheck = set([0 if i in coordinates else 1
                                        for i in range(si, end+1)])

                        for i in range(si, end+1):
                            if i in coordinates:
                                del coordinates[i]
                            else:
                                if self.verbose:
                                    self.logger.error("seqsearch - should't be here "
                                                      + locus + " - "
                                                      + " - " + feat_name)
                            mapping[i] = feat_name

                        found_feats.update({feat_name:
                                            SeqFeature(
                                                FeatureLocation(
                                                    ExactPosition(start),
                                                    ExactPosition(end), strand=1),
                                                type=feat_name)})

                        if self.verbose and self.verbosity > 0:
                            self.logger.info("Coordinates | Start = " + str(start) + " - End = " + str(end))
                    else:
                        if self.verbose and self.verbosity > 0:
                            self.logger.info("Adding ambig feature " + feat_name)
                        feat_missing.update({feat_name: feats[feat_name]})
                        ambig_map.update({feat_name:
                                          seq_search[1:len(seq_search)]})
                else:
                    if self.verbose and self.verbosity > 0:
                        self.logger.info("Adding ambig feature " + feat_name)
                    feat_missing.update({feat_name: feats[feat_name]})
                    ambig_map.update({feat_name: seq_search[1:len(seq_search)]})
            else:
                if self.verbose and self.verbosity > 1:
                    self.logger.info("No match for " + feat_name)
                feat_missing.update({feat_name: feats[feat_name]})

        blocks = getblocks(coordinates)
        exact_matches = list(found_feats.keys())

        # * HARD CODED LOGIC * #
        # > 
        #
        #  HLA-DRB1 exon3 exact match - with intron1 and 3 missing
        if('exon_3' in exact_matches and run == 99 and locus == 'HLA-DRB1'
           and 'exon_2' in feat_missing and (len(blocks) == 1 or len(blocks) == 2)):

            for b in blocks:
                x = b[len(b)-1]
                if x == max(list(mapping.keys())):
                    featname = "intron_3"
                    found_feats.update({featname:
                                        SeqFeature(
                                            FeatureLocation(
                                                ExactPosition(b[0]-1),
                                                ExactPosition(b[len(b)-1]),
                                                strand=1),
                                            type=featname)})
                else:
                    featname = "exon_2"
                    found_feats.update({featname:
                                        SeqFeature(
                                            FeatureLocation(
                                                ExactPosition(b[0]),
                                                ExactPosition(b[len(b)-1]),
                                                strand=1),
                                            type=featname)})
                    seq_covered -= len(b)

                if self.verbose and self.verbosity > 1:
                    self.logger.info("Successfully annotated class DRB1 II sequence")

                return Annotation(features=found_feats,
                                  covered=seq_covered,
                                  seq=in_seq,
                                  missing=feat_missing,
                                  ambig=ambig_map,
                                  method=method,
                                  mapping=mapping,
                                  exact_match=exact_matches)

        # If it's a class II sequence and
        # exon_2 is an exact match
        # * HARD CODED LOGIC * #
        # > It's common for exon2 to be fully sequenced
        #   but intron_2 and intron_1 to be partially sequenced,
        #   which can make it hard to annotate those to features.
        #   If there are two missing blocks that is small enough
        #   and they are before and after exon2, then it's very
        #   very likely to be intron_2 and intron_1.
        if 'exon_2' in exact_matches and len(blocks) == 2 \
                and is_classII(locus) and seq_covered < 300:

            if self.verbose and self.verbosity > 1:
                self.logger.info("Running search for class II sequence")

            r = True
            for b in blocks:
                x = b[len(b)-1]
                if x == max(list(mapping.keys())):
                    x = b[0]-1
                else:
                    x += 1
                f = mapping[x]
                if f != 'exon_2':
                    r = False
            if r:
                for b in blocks:
                    x = b[len(b)-1]
                    if x == max(list(mapping.keys())):
                        featname = "intron_2"
                        found_feats.update({featname:
                                            SeqFeature(
                                                FeatureLocation(
                                                    ExactPosition(b[0]-1),
                                                    ExactPosition(b[len(b)-1]),
                                                    strand=1),
                                                type=featname)})
                    else:
                        featname = "intron_1"
                        found_feats.update({featname:
                                            SeqFeature(
                                                FeatureLocation(
                                                    ExactPosition(b[0]),
                                                    ExactPosition(b[len(b)-1]),
                                                    strand=1),
                                                type=featname)})
                    seq_covered -= len(b)

                if self.verbose and self.verbosity > 1:
                    self.logger.info("Successfully annotated class II sequence")

                return Annotation(features=found_feats,
                                  covered=seq_covered,
                                  seq=in_seq,
                                  missing=feat_missing,
                                  ambig=ambig_map,
                                  method=method,
                                  mapping=mapping,
                                  exact_match=exact_matches)

        annotated_feats, mb, mapping = self._resolve_unmapped(blocks,
                                                              feat_missing,
                                                              ambig_map,
                                                              mapping,
                                                              found_feats,
                                                              locus,
                                                              seq_covered
                                                              )

        # * HARD CODED LOGIC * #
        if(not mb and blocks and len(feat_missing.keys()) == 0
           and len(ambig_map.keys()) == 0):
            mb = blocks

        if mb:

            # Unmap exon 8
            if locus in ['HLA-C', 'HLA-A'] and len(in_seq.seq) < 3000 \
                    and 'exon_8' in exact_matches:
                for i in deleted_coords:
                    mapping[i] = 1
                coordinates.update(deleted_coords)
                mb = getblocks(coordinates)
                feat_missing.update(added_feat)

                # Delte from found features
                del exact_matches[exact_matches.index('exon_8')]
                del found_feats['exon_8']

                if 'exon_8' in annotated_feats:
                    del annotated_feats['exon_8']
                if 'three_prime_UTR' in found_feats:
                    del found_feats['three_prime_UTR']
                if 'three_prime_UTR' in annotated_feats:
                    del annotated_feats['three_prime_UTR']

            refmissing = [f for f in structures[locus]
                          if f not in annotated_feats]

            if self.verbose and self.verbosity > 1:
                self.logger.info("* Annotation not complete *")

            # Print out what features were missing by the ref
            if self.verbose and self.verbosity > 2:
                self.logger.info("Refseq was missing these features = " + ",".join(list(refmissing)))

            # Print out what features were ambig matches
            if self.verbose and self.verbosity > 1 and len(ambig_map) > 1:
                self.logger.info("Features with ambig matches = " + ",".join(list(ambig_map)))

            # Print out what features were exact matches
            if self.verbose and self.verbosity > 2 and len(exact_matches) > 1:
                self.logger.info("Features exact matches = " + ",".join(list(exact_matches)))

            # Print out what features have been annotated
            if self.verbose and self.verbosity > 1 and len(annotated_feats) > 1:
                self.logger.info("Features annotated = " + ",".join(list(annotated_feats)))

            # Print out what features are missing
            if self.verbose and self.verbosity > 1 and len(feat_missing) > 1:
                self.logger.info("Features missing = " + ",".join(list(feat_missing)))

            annotation = Annotation(features=annotated_feats,
                                    covered=seq_covered,
                                    seq=in_seq,
                                    missing=feat_missing,
                                    ambig=ambig_map,
                                    blocks=mb,
                                    method=method,
                                    refmissing=refmissing,
                                    mapping=mapping,
                                    exact_match=exact_matches,
                                    annotation=None)
        else:

            mb = None
            # Unmap exon 8
            if locus in ['HLA-C', 'HLA-A'] and len(in_seq.seq) < 600 \
                    and 'exon_8' in exact_matches \
                    and 'three_prime_UTR' in annotated_feats\
                    and 'three_prime_UTR' not in exact_matches:

                for i in deleted_coords:
                    mapping[i] = 1

                coordinates.update(deleted_coords)
                mb = getblocks(coordinates)
                feat_missing.update(added_feat)
                del exact_matches[exact_matches.index('exon_8')]
                del found_feats['exon_8']
                if 'exon_8' in annotated_feats:
                    del annotated_feats['exon_8']
                if 'three_prime_UTR' in found_feats:
                    del found_feats['three_prime_UTR']
                if 'three_prime_UTR' in annotated_feats:
                    del annotated_feats['three_prime_UTR']

            if self.verbose:
                self.logger.info("* No missing blocks after seq_search *")

            # Print out what features were ambig matches
            if self.verbose and self.verbosity > 0 and len(ambig_map) > 1:
                self.logger.info("Features with ambig matches = " + ",".join(list(ambig_map)))

            # Print out what features were exact matches
            if self.verbose and self.verbosity > 0 and len(exact_matches) > 1:
                self.logger.info("Features exact matches = " + ",".join(list(exact_matches)))

            # Print out what features have been annotated
            if self.verbose and self.verbosity > 0 and len(annotated_feats) > 1:
                self.logger.info("Features annotated = " + ",".join(list(annotated_feats)))

            # Print out what features are missing
            if self.verbose and self.verbosity > 0 and len(feat_missing) > 1:
                self.logger.info("Features missing = " + ",".join(list(feat_missing)))

            annotation = Annotation(features=annotated_feats,
                                    covered=seq_covered,
                                    seq=in_seq,
                                    missing=feat_missing,
                                    ambig=ambig_map,
                                    method=method,
                                    blocks=mb,
                                    mapping=mapping,
                                    exact_match=exact_matches,
                                    annotation=None)

        return annotation

    # TODO: This should be cleaned up and reduced..
    #       There is some repeating logic and overlap in
    #       this section.
    def _resolve_unmapped(self, blocks, feat_missing, ambig_map,
                          mapping, found_feats, loc, covered, rerun=False):

        structures = get_structures()
        struct_order = get_structorder()
        structure_max = get_structmax()
        exon_only = True
        found_exons = 0
        for f in found_feats:
            if re.search("intron", f) or re.search("UTR", f):
                exon_only = False

            if re.search("exon", f):
                found_exons += 1

        # Count the number of exons for the given loci
        num_exons = 0
        for f in structures[loc]:
            if re.search("exon", f):
                num_exons += 1

        # If all exons have been mapped
        # then it is not exon only data
        if found_exons == num_exons:
            exon_only = False

        # If it's exon only, then search two
        # features up rather than one
        add_num = 2 if exon_only and rerun and covered < 300 else 1

        block_mapped = []
        missing_blocks = []
        for b in blocks:
            for featname in ambig_map.keys():
                locats = ambig_map[featname]
                start_i = b[0]-1
                end_i = b[len(b)-1]+1

                # TODO: Catch ERROR
                #if not end_i in mapping:
                feat_num = structures[loc][featname]
                x = feat_num-add_num
                y = feat_num-add_num
                if feat_num+add_num <= structure_max[loc] \
                        and feat_num-add_num >= 0 and start_i >= 0 \
                        and end_i <= len(mapping) - 1 \
                        and x in struct_order[loc] \
                        and y in struct_order[loc]:
                        expected_p = struct_order[loc][feat_num-add_num]
                        expected_n = struct_order[loc][feat_num+add_num]
                        previous_feat = mapping[start_i]
                        next_feat = mapping[end_i]
                        if expected_p == previous_feat \
                            and expected_n == next_feat \
                            and expected_p != 1 \
                                and b[0]-1 in locats:
                                block_mapped.append(b)
                                found_feats.update({featname:
                                                    SeqFeature(
                                                        FeatureLocation(
                                                            ExactPosition(b[0]-1),
                                                            ExactPosition(b[len(b)-1]),
                                                            strand=1),
                                                        type=featname)})
                elif feat_num+add_num > structure_max[loc] \
                        and feat_num-add_num >= 0 and start_i >= 0 \
                        and end_i >= max(mapping) \
                        and y in struct_order[loc]:
                        expected_p = struct_order[loc][feat_num-add_num]
                        previous_feat = mapping[start_i]
                        if expected_p == previous_feat \
                            and expected_p != 1 \
                                and b[0]-1 in locats:
                                block_mapped.append(b)
                                found_feats.update({featname:
                                                    SeqFeature(
                                                        FeatureLocation(
                                                            ExactPosition(b[0]-1),
                                                            ExactPosition(b[len(b)-1]),
                                                            strand=1),
                                                        type=featname)})
                elif feat_num+add_num <= structure_max[loc] \
                        and feat_num-add_num < 0\
                        and x in struct_order[loc]:
                    expected_n = struct_order[loc][feat_num+add_num]
                    if not end_i in mapping:
                        next_feat = mapping[end_i-1]
                    else:
                        next_feat = mapping[end_i]
                    if expected_n == next_feat \
                        and expected_p != 1 \
                            and b[0]-1 in locats:
                            block_mapped.append(b)
                            found_feats.update({featname:
                                                SeqFeature(
                                                    FeatureLocation(
                                                        ExactPosition(b[0]-1),
                                                        ExactPosition(b[len(b)-1]),
                                                        strand=1),
                                                    type=featname)})
                else:
                    missing_blocks.append(b)

        for b in blocks:
            for featname in feat_missing.keys():
                if featname in ambig_map and b not in block_mapped:
                    if b not in missing_blocks:
                        missing_blocks.append(b)
                    continue

                if b not in block_mapped:
                    #featlen = feat_missing[featname]
                    start_i = b[0]-1
                    end_i = b[len(b)-1]+1
                    feat_num = structures[loc][featname]

                    if feat_num+add_num <= structure_max[loc] \
                        and feat_num-1 >= 1 \
                            and end_i <= max(mapping.keys()) \
                            and start_i >= 0 \
                            and feat_num-add_num > 0:

                        expected_p = struct_order[loc][feat_num-add_num]
                        expected_n = struct_order[loc][feat_num+add_num]
                        previous_feat = mapping[start_i]
                        next_feat = mapping[end_i]
                        if expected_p == previous_feat \
                            and expected_n == next_feat \
                                and expected_p != 1 \
                                and expected_n != 1:
                                if b in missing_blocks:
                                    del missing_blocks[missing_blocks.index(b)]
                                for i in b:
                                    mapping.update({i: featname})
                                block_mapped.append(b)
                                found_feats.update({featname:
                                                    SeqFeature(
                                                        FeatureLocation(
                                                            ExactPosition(b[0]-1),
                                                            ExactPosition(b[len(b)-1]),
                                                            strand=1),
                                                        type=featname)})
                        else:
                            if b not in missing_blocks:
                                missing_blocks.append(b)
                    elif feat_num+add_num > structure_max[loc] \
                            and feat_num-add_num >= 1 and start_i >= 0:
                        expected_p = struct_order[loc][feat_num-add_num]
                        previous_feat = mapping[start_i]
                        if expected_p == previous_feat \
                                and expected_p != 1:
                                if b in missing_blocks:
                                    del missing_blocks[missing_blocks.index(b)]
                                block_mapped.append(b)
                                for i in b:
                                    mapping.update({i: featname})
                                found_feats.update({featname:
                                                    SeqFeature(
                                                        FeatureLocation(
                                                            ExactPosition(b[0]-1),
                                                            ExactPosition(b[len(b)-1]),
                                                            strand=1),
                                                        type=featname)})
                        else:
                            if b not in missing_blocks:
                                missing_blocks.append(b)
                    elif(feat_num+add_num <= structure_max[loc]
                            and feat_num-add_num < 1
                            and end_i <= max(mapping)
                            and end_i in mapping):
                        expected_n = struct_order[loc][feat_num+add_num]
                        next_feat = mapping[end_i]
                        if expected_n == next_feat:
                            if b in missing_blocks:
                                del missing_blocks[missing_blocks.index(b)]
                            add = 0
                            block_mapped.append(b)
                            for i in b:
                                mapping.update({i: featname})
                            if add != 0:
                                for i in range(b[len(b)-1], b[len(b)-1]+add):
                                    mapping.update({i: featname})
                            found_feats.update({featname:
                                                SeqFeature(
                                                    FeatureLocation(
                                                        ExactPosition(b[0]),
                                                        ExactPosition(b[len(b)-1]+add),
                                                        strand=1),
                                                    type=featname)})
                        else:
                            if b not in missing_blocks:
                                missing_blocks.append(b)
                    else:
                        if b not in missing_blocks:
                            missing_blocks.append(b)

        # If it failed to map all features when only looking
        # at the exons, then try again and look at all features
        if exon_only and not rerun and missing_blocks:
            if self.verbose:
                self.logger.info("Rerunning seqsearch to look at all features")
            return self._resolve_unmapped(missing_blocks, feat_missing,
                                          ambig_map,
                                          mapping, found_feats,
                                          loc, covered, rerun=True)
        else:
            return found_feats, missing_blocks, mapping

    @property
    def verbose(self) -> bool:
        """
        Gets the verbose of this SeqSearch.

        :return: The verbose of this SeqSearch.
        :rtype: bool
        """
        return self._verbose

    @verbose.setter
    def verbose(self, verbose: bool):
        """
        Sets the verbose of this SeqSearch.

        :param verbose: The verbose of this SeqSearch.
        :type verbose: bool
        """
        self._verbose = verbose

    @property
    def verbosity(self) -> int:
        """
        Gets the verbosity of this SeqSearch.

        :return: The verbosity of this SeqSearch.
        :rtype: int
        """
        return self._verbosity

    @verbosity.setter
    def verbosity(self, verbosity: int):
        """
        Sets the verbosity of this SeqSearch.

        :param verbosity: The verbosity of this SeqSearch.
        :type verbosity: int
        """
        self._verbosity = verbosity


# TODO: Move to util.py
def loctype(s1, e1, s2, e2):
    if s1 < s2 and e1 < e2:
        return True
    else:
        return False


# TODO: Move to util.py
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
