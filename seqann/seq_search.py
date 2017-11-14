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
import sys

from Bio.SeqUtils import nt_search
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import ExactPosition

from seqann.models.reference_data import ReferenceData
from seqann.models.annotation import Annotation
from seqann.models.base_model_ import Model
from seqann.util import deserialize_model
from seqann.util import get_features


def getblocks(coords):
    block = []
    blocks = []
    sorted_i = sorted(coords.keys())
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


class SeqSearch(Model):
    '''
    classdocs
    '''
    def __init__(self, refdata: ReferenceData=None, verbose: bool=False):
        """
        ReferenceData
        :param refdata: The reference data of this SeqSearch.
        :type refdata: ReferenceData
        """
        self.data_types = {
            'refdata': ReferenceData,
            'verbose': bool
        }

        self.attribute_map = {
            'refdata': 'refdata',
            'verbose': 'verbose'
        }
        if not refdata:
            refdata = ReferenceData()

        self._refdata = refdata
        self._verbose = verbose

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

    def search_seqs(self, seqrec, in_seq, locus, partial_ann=None):

        # Extract out the sequences and feature names
        # from the reference sequences
        # feats = [[feat.type, feat.extract(seqrec.seq)]
        #          for feat in seqrec.features if feat.type != "source"
        #          and feat.type != "CDS" and isinstance(feat, SeqFeature)]

        # The mapped features will be subtracted from seq_covered
        # so the final seq_covered number will reflect the remaining
        # number of base pairs that haven't been mapped.
        #
        # The coordinates and mapping will help determine what positions
        # in the sequence have been mapped and to what features. The
        # missing blocks variable will be generated using these.
        seq_covered = len(in_seq.seq)
        coordinates = dict(map(lambda x: [x, 1],
                           [i for i in range(0, len(in_seq.seq)+1)]))

        mapping = dict(map(lambda x: [x, 1],
                           [i for i in range(0, len(in_seq.seq)+1)]))

        ambig_map = {}
        found_feats = {}
        feat_missing = {}

        # If the partial annotation is provided
        # then make the found_feats equal to
        # what has already been annotated
        if partial_ann:
            if self.verbose:
                print("seqsearch -> partial annotation", file=sys.stderr)
            found_feats = partial_ann.features
            coordinates = dict(map(lambda l: [l, 1],
                                   [item for sublist
                                    in partial_ann.blocks
                                    for item in sublist]))
            seq_covered = partial_ann.covered
            mapping = partial_ann.mapping

        feats = get_features(seqrec)
        for feat_name in feats:

            # skip if partial annotation is provided
            # and the feat name is not one of the
            # missing features
            if partial_ann and feat_name not in partial_ann.refmissing:
                continue

            # Search for the reference feature sequence in the
            # input sequence. Record the coordinates if it's
            # found and if it's found in multiple spots. If it
            # is not found, then record that feature as missing.
            seq_search = nt_search(str(in_seq.seq), str(feats[feat_name]))
            if len(seq_search) == 2:
                seq_covered -= len(str(feats[feat_name]))
                end = int(len(str(feats[feat_name])) + seq_search[1])

                # If the feature is found and it's a five_prime_UTR then
                # the start should always be 0, so insertions at the
                # beinging of the sequence will be found.
                start = seq_search[1] if feat_name != 'five_prime_UTR' else 0
                found_feats.update({feat_name:
                                    SeqFeature(
                                        FeatureLocation(
                                            ExactPosition(start),
                                            ExactPosition(end), strand=1),
                                        type=feat_name)})
                si = seq_search[1]+1 if seq_search[1] != 0 and \
                    feat_name != 'five_prime_UTR' else 0
                for i in range(si, end+1):
                    if i in coordinates:
                        del coordinates[i]
                    else:
                        if self.verbose:
                            print("seqsearch - shouldnt be here!!",
                                  file=sys.stderr)
                    mapping[i] = feat_name
            elif(len(seq_search) > 2):
                feat_missing.update({feat_name: feats[feat_name]})
                ambig_map.update({feat_name: seq_search[1:len(seq_search)]})
            else:
                feat_missing.update({feat_name: feats[feat_name]})

        blocks = getblocks(coordinates)

        # TODO: pass seq_covered and mapping, so the
        #       final annotation contains the updated values
        # print("FOUND")
        # print(found_feats)
        # print("MISSING")
        # print(list(feat_missing.keys()))
        # print("BLOCK")
        # print(blocks)
        #print(ambig_map)
        annotated_feats, mb = self._resolve_unmapped(blocks,
                                                     feat_missing,
                                                     ambig_map, mapping,
                                                     found_feats, locus)

        method = "nt_search" if not partial_ann else partial_ann.method

        #
        if mb:
            refmissing = [f for f in self.refdata.structures[locus]
                          if f not in annotated_feats]

            annotation = Annotation(features=annotated_feats,
                                    covered=seq_covered,
                                    seq=in_seq,
                                    missing=feat_missing,
                                    ambig=ambig_map,
                                    blocks=mb,
                                    method=method,
                                    refmissing=refmissing,
                                    mapping=mapping)
        else:
            annotation = Annotation(features=annotated_feats,
                                    covered=seq_covered,
                                    seq=in_seq,
                                    missing=feat_missing,
                                    ambig=ambig_map,
                                    method=method,
                                    mapping=mapping)

        return annotation

    def _resolve_unmapped(self, blocks, feat_missing, ambig_map,
                          mapping, found_feats, loc, rerun=False):

        exon_only = True
        found_exons = 0
        for f in found_feats:
            if re.search("intron", f) or re.search("UTR", f):
                exon_only = False
            if re.search("exon", f):
                found_exons += 1

        # Count the number of exons for the given loci
        num_exons = 0
        for f in self.refdata.structures[loc]:
            if re.search("exon", f):
                num_exons += 1

        # If all exons have been mapped
        # then it is not exon only data
        if found_exons == num_exons or rerun:
            exon_only = False

        # If it's exon only, then search two
        # features up rather than one
        add_num = 2 if exon_only else 1

        block_mapped = []
        missing_blocks = []
        for b in blocks:
            for featname in ambig_map.keys():
                locats = ambig_map[featname]
                start_i = b[0]-1
                end_i = b[len(b)-1]+1
                feat_num = self.refdata.structures[loc][featname]
                if feat_num+add_num <= self.refdata.structure_max[loc] \
                        and feat_num-add_num >= 0 and start_i >= 0 \
                        and end_i <= len(mapping) - 1:
                        x = feat_num-add_num
                        expected_p = self.refdata.struct_order[loc][feat_num-add_num]
                        expected_n = self.refdata.struct_order[loc][feat_num+add_num]
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
                elif feat_num+add_num > self.refdata.structure_max[loc] \
                        and feat_num-add_num >= 0 and start_i >= 0 \
                        and end_i > len(mapping) - 1:
                        expected_p = self.refdata.struct_order[loc][feat_num-add_num]
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
                elif feat_num+add_num <= self.refdata.structure_max[loc] and feat_num-add_num < 0:
                    expected_n = self.refdata.struct_order[loc][feat_num+1]
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
                if b not in block_mapped:
                    #featlen = feat_missing[featname]
                    start_i = b[0]-1
                    end_i = b[len(b)-1]+1
                    feat_num = self.refdata.structures[loc][featname]
                    if feat_num+1 <= self.refdata.structure_max[loc] \
                        and feat_num-1 >= 1 \
                            and end_i <= len(mapping) - 1 \
                            and start_i >= 0 \
                            and feat_num-add_num > 0:
                        expected_p = self.refdata.struct_order[loc][feat_num-add_num]
                        expected_n = self.refdata.struct_order[loc][feat_num+add_num]
                        previous_feat = mapping[start_i]
                        next_feat = mapping[end_i]
                        if expected_p == previous_feat \
                            and expected_n == next_feat \
                                and expected_p != 1 \
                                and expected_n != 1:
                                if b in missing_blocks:
                                    del missing_blocks[missing_blocks.index(b)]
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
                    elif feat_num+add_num > self.refdata.structure_max[loc] \
                            and feat_num-add_num >= 1 and start_i >= 0:
                        expected_p = self.refdata.struct_order[loc][feat_num-add_num]
                        previous_feat = mapping[start_i]
                        if expected_p == previous_feat \
                                and expected_p != 1:
                                if b in missing_blocks:
                                    del missing_blocks[missing_blocks.index(b)]
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
                    elif(feat_num+add_num <= self.refdata.structure_max[loc] and feat_num-add_num < 1 and end_i <= len(mapping) - 1):
                        expected_n = self.refdata.struct_order[loc][feat_num+add_num]
                        next_feat = mapping[end_i]
                        if expected_n == next_feat:
                            if b in missing_blocks:
                                del missing_blocks[missing_blocks.index(b)]
                            block_mapped.append(b)
                            found_feats.update({featname:
                                                SeqFeature(
                                                    FeatureLocation(
                                                        ExactPosition(b[0]),
                                                        ExactPosition(b[len(b)-1]),
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
                print("RERUNNING seqsearch", file=sys.stderr)
            return self._resolve_unmapped(missing_blocks, feat_missing,
                                          ambig_map,
                                          mapping, found_feats,
                                          loc, rerun=True)
        else:
            # print("_resolve_unmapped")
            # print(missing_blocks)
            return found_feats, missing_blocks

    @property
    def refdata(self) -> ReferenceData:
        """
        Gets the refdata of this SeqSearch.

        :return: The refdata of this SeqSearch.
        :rtype: ReferenceData
        """
        return self._refdata

    @refdata.setter
    def refdata(self, refdata: ReferenceData):
        """
        Sets the refdata of this SeqSearch.

        :param refdata: The refdata of this SeqSearch.
        :type refdata: ReferenceData
        """
        self._refdata = refdata

    @property
    def verbose(self) -> bool:
        """
        Gets the verbose of this SeqSearch.

        :return: The verbose of this SeqSearch.
        :rtype: ReferenceData
        """
        return self._verbose

    @verbose.setter
    def verbose(self, verbose: bool):
        """
        Sets the verbose of this SeqSearch.

        :param verbose: The verbose of this SeqSearch.
        :type verbose: ReferenceData
        """
        self._verbose = verbose
