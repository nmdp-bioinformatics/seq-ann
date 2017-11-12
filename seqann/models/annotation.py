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
from seqann.models.base_model_ import Model
from datetime import date, datetime
from typing import List, Dict
from BioSQL.BioSeq import DBSeqRecord
from ..util import deserialize_model
from Bio.SeqRecord import SeqRecord


class Annotation(Model):
    '''
    classdocs
    '''
    def __init__(self, missing: Dict=None, ambig: Dict=None,
                 seq: SeqRecord=None, features: Dict=None, covered: int=None,
                 annotation: Dict={}, blocks: List[List[int]]=None,
                 method: str=None, mapping: Dict=None,
                 refmissing: List[str]=None,
                 complete_annotation: bool=False):
        """
        Annotation

        :param features: The features of this Annotation.
        :type features: Dict
        :param covered: The sequence coverage of this Annotation.
        :type covered: str
        """
        self.data_types = {
            'seq': SeqRecord,
            'features': Dict,
            'covered': int,
            'missing': Dict,
            'ambig': Dict,
            'annotation': Dict,
            'complete_annotation': bool,
            'blocks': List[List[int]],
            'method': str,
            'mapping': Dict,
            'refmissing': List[str]
        }

        self.attribute_map = {
            'seq': 'seq',
            'features': 'features',
            'covered': 'covered',
            'missing': 'missing',
            'ambig': 'ambig',
            'annotation': 'annotation',
            'complete_annotation': 'complete_annotation',
            'blocks': 'blocks',
            'method': 'method',
            'mapping': 'mapping',
            'refmissing': 'refmissing'
        }

        self._refmissing = refmissing
        self._mapping = mapping
        self._method = method
        self._blocks = blocks
        self._seq = seq
        self._features = features
        self._covered = covered
        self._missing = missing
        self._ambig = ambig
        self._complete_annotation = complete_annotation
        self._annotation = annotation

        missing_blocks = {}
        if not annotation:
            self._complete_annotation = True
            if not features:
                self._complete_annotation = False
            else:
                for feat in missing:
                    if feat not in features:
                        self._complete_annotation = False
                        missing_blocks.update({feat: missing[feat]})

                for feat in ambig:
                    if feat not in features:
                        self._complete_annotation = False
                        missing_blocks.update({feat: missing[feat]})

                full_annotation = {}
                for feat in features:
                    f = features[feat]
                    full_annotation.update({feat: f.extract(seq)})
                self._annotation = full_annotation
            self._missing = missing_blocks

            if blocks:
                self._complete_annotation = False
            else:
                self._complete_annotation = True
        else:
            self._annotation = annotation

    @classmethod
    def from_dict(cls, dikt) -> 'Annotation':
        """
        Returns the dict as a model

        :param dikt: A dict.
        :type: dict
        :return: The Annotation of this Annotation.
        :rtype: Annotation
        """
        return deserialize_model(dikt, cls)

    @property
    def complete_annotation(self) -> bool:
        """
        Gets the complete_annotation of this Annotation.

        :return: The complete_annotation of this Annotation.
        :rtype: bool
        """
        return self._complete_annotation

    @complete_annotation.setter
    def complete_annotation(self, complete_annotation: bool):
        """
        Sets the complete_annotation of this Annotation.

        :param complete_annotation: The complete_annotation of this Annotation.
        :type complete_annotation: bool
        """
        self._complete_annotation = complete_annotation

    @property
    def features(self) -> Dict:
        """
        Gets the features of this Annotation.

        :return: The features of this Annotation.
        :rtype: Dict
        """
        return self._features

    @features.setter
    def features(self, features: Dict):
        """
        Sets the features of this Annotation.

        :param features: The features of this Annotation.
        :type features: Dict
        """
        self._features = features

    @property
    def covered(self) ->int:
        """
        Gets the coverage of this Annotation.

        :return: The covered of this Annotation.
        :rtype: int
        """
        return self._covered

    @covered.setter
    def covered(self, covered: int):
        """
        Sets the covered of this Annotation.

        :param covered: The covered of this Annotation.
        :type covered: int
        """
        self._covered = covered

    @property
    def seq(self) ->SeqRecord:
        """
        Gets the coverage of this Annotation.

        :return: The seq of this Annotation.
        :rtype: SeqRecord
        """
        return self._seq

    @seq.setter
    def seq(self, seq: SeqRecord):
        """
        Sets the seq of this Annotation.

        :param seq: The seq of this Annotation.
        :type seq: SeqRecord
        """
        self._seq = seq

    @property
    def ambig(self) ->Dict:
        """
        Gets the ambig of this Annotation.

        :return: The ambig of this Annotation.
        :rtype: Dict
        """
        return self._ambig

    @ambig.setter
    def ambig(self, ambig: Dict):
        """
        Sets the ambig of this Annotation.

        :param ambig: The ambig of this Annotation.
        :type ambig: Dict
        """
        self._ambig = ambig

    @property
    def method(self) ->str:
        """
        Gets the method of this Annotation.

        :return: The method of this Annotation.
        :rtype: str
        """
        return self._method

    @method.setter
    def method(self, method: str):
        """
        Sets the method of this Annotation.

        :param method: The method of this Annotation.
        :type method: str
        """
        self._method = method

    @property
    def missing(self) ->Dict:
        """
        Gets the missing of this Annotation.

        :return: The missing of this Annotation.
        :rtype: Dict
        """
        return self._missing

    @missing.setter
    def missing(self, missing: Dict):
        """
        Sets the missing of this Annotation.

        :param missing: The missing of this Annotation.
        :type missing: Dict
        """
        self._missing = missing

    @property
    def mapping(self) ->Dict:
        """
        Gets the coverage of this Annotation.

        :return: The seq of this Annotation.
        :rtype: Dict
        """
        return self._mapping

    @mapping.setter
    def mapping(self, mapping: Dict):
        """
        Sets the mapping of this Annotation.

        :param mapping: The mapping of this Annotation.
        :type mapping: Dict
        """
        self._mapping = mapping

    @property
    def refmissing(self) ->List[str]:
        """
        Gets the refmissing of this Annotation.

        :return: The refmissing of this Annotation.
        :rtype: List[str]
        """
        return self._refmissing

    @refmissing.setter
    def refmissing(self, refmissing: List[str]):
        """
        Sets the refmissing of this Annotation.

        :param refmissing: The refmissing of this Annotation.
        :type refmissing: List[str]
        """
        self._refmissing = refmissing

    @property
    def annotation(self) ->Dict:
        """
        Gets the coverage of this Annotation.

        :return: The seq of this Annotation.
        :rtype: Dict
        """
        return self._annotation

    @annotation.setter
    def annotation(self, annotation: Dict):
        """
        Sets the seq of this Annotation.

        :param annotation: The annotation of this Annotation.
        :type annotation: Dict
        """
        self._annotation = annotation

    @property
    def blocks(self) ->List[List[int]]:
        """
        Gets the coverage of this Annotation.

        :return: The blocks of this Annotation.
        :rtype: List[List[int]]
        """
        return self._blocks

    @blocks.setter
    def blocks(self, blocks: List[List[int]]):
        """
        Sets the blocks of this Annotation.

        :param blocks: The blocks of this Annotation.
        :type blocks: List[List[int]]
        """
        self._blocks = blocks

    def check_annotation(self):

        self.complete_annotation = True
        self.method = "nt_search and clustalo"
        # print("check_annotation")
        # print(self.missing)
        # print(self.ambig)
        # print("")
        if self.missing:
            for feat in self.missing:
                if feat not in self.annotation:
                    self.complete_annotation = False

        if self.ambig:
            for feat in self.ambig:
                if feat not in self.annotation:
                    self.complete_annotation = False

        if self.blocks:
            self.complete_annotation = False
        else:
            self.complete_annotation = True

    def clean(self):
        self.missing = ''
        self.ambig = ''
        self.mapping = ''


