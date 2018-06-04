#!/usr/bin/env python
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

"""
test_seqann
----------------------------------

Tests for `seqann` module.
"""
import os
import unittest

from Bio.Seq import Seq

from seqann.util import isutr
from seqann.util import is_kir
from seqann.util import checkseq
from seqann.util import randomid
from seqann.util import is_classII


class TestUtil(unittest.TestCase):

    def setUp(self):
        self.data_dir = os.path.dirname(__file__) + "/resources"
        pass

    def test_001_isutr(self):
        self.assertTrue(isutr('three_prime_UTR'))
        self.assertTrue(isutr('five_prime_UTR'))
        self.assertFalse(isutr('exon-3'))
        pass

    def test_002_iskir(self):
        self.assertTrue(is_kir('KIR*3DL1'))
        self.assertFalse(is_kir('exon-3'))
        pass

    def test_003_is_classII(self):
        self.assertTrue(is_classII('HLA-DRB1*15:01'))
        self.assertTrue(is_classII('HLA-DQB1*06:01'))
        self.assertFalse(is_classII('HLA-A*02:01'))
        pass

    def test_004_checkseq(self):
        self.assertTrue(checkseq(Seq('AAACTGATCG')))
        self.assertTrue(checkseq(Seq('AAACTGATCGGGGGAAACCCTTT')))
        self.assertFalse(checkseq(Seq('AAACTGATCGGGGGAAACCCTTTNN')))
        self.assertFalse(checkseq(Seq('NNNNAAACTGATCGGGGGAAACCCTTTNNNN')))
        self.assertFalse(checkseq(Seq('AAACTGATCGGGGGAAACCCTTTZ')))
        pass

    def test_005_randomid(self):
        self.assertEqual(len(randomid()), 12)
        self.assertEqual(len(randomid(N=6)), 6)
        self.assertEqual(len(randomid(N=13)), 13)
        pass
