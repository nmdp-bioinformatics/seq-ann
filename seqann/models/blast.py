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
from seqann.models.base_model_ import Model
from typing import List, Dict
from BioSQL.BioSeq import DBSeqRecord
from ..util import deserialize_model


# NOTE: This really doesn't need to be a class..

class Blast(Model):
    '''
    classdocs
    '''
    def __init__(self, failed: bool=None, match_seqs: List[DBSeqRecord]=None,
                 alleles: List[str]=None):
        """
        :param match_seqs: The match_seqs of this Blast.
        :type match_seqs: List[``SeqRecord``]
        :param alleles: The alleles of this Blast.
        :type alleles: List[``str``]
        """
        self.data_types = {
            'match_seqs': List[DBSeqRecord],
            'alleles': List[str],
            'failed': bool
        }

        self.attribute_map = {
            'match_seqs': 'match_seqs',
            'alleles': 'alleles',
            'failed': 'failed'
        }
        self._match_seqs = match_seqs
        self._alleles = alleles

        if not failed:
            if len(match_seqs) >= 1 and len(alleles) >= 1:
                self._failed = False
            else:
                self._failed = True
        else:
            self._failed = failed

    @classmethod
    def from_dict(cls, dikt) -> 'Blast':
        """
        Returns the dict as a model

        :param dikt: A dict.
        :type: dict
        :return: The Blast of this Blast.
        :rtype: Blast
        """
        return deserialize_model(dikt, cls)

    @property
    def match_seqs(self) -> List[DBSeqRecord]:
        """
        Gets the match_seqs of this Blast.

        :return: The match_seqs of this Blast.
        :rtype: List[str]
        """
        return self._match_seqs

    @match_seqs.setter
    def match_seqs(self, match_seqs: List[DBSeqRecord]):
        """
        Sets the match_seqs of this Blast.

        :param match_seqs: The match_seqs of this Blast.
        :type match_seqs: List[str]
        """
        self._match_seqs = match_seqs

    @property
    def alleles(self) -> List[str]:
        """
        Gets the alleles of this Blast.

        :return: The alleles of this Blast.
        :rtype: List[str]
        """
        return self._alleles

    @alleles.setter
    def alleles(self, alleles: List[str]):
        """
        Sets the alleles of this Blast.

        :param alleles: The alleles of this Blast.
        :type alleles: List[str]
        """
        self._alleles = alleles

    @property
    def failed(self) -> bool:
        """
        Gets the failed of this Blast.

        :return: The failed of this Blast.
        :rtype: bool
        """
        return self._failed

    @failed.setter
    def failed(self, failed: bool):
        """
        Sets the failed of this Blast.

        :param failed: The failed of this Blast.
        :type failed: List[str]
        """
        self._failed = failed
