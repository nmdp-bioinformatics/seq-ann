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
import json
import unittest
import warnings

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from seqann.gfe import GFE
from Bio.Alphabet import generic_dna
from seqann.models.annotation import Annotation
from seqann.feature_client.models.feature import Feature


def ignore_warnings(test_func):
    def do_test(self, *args, **kwargs):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ResourceWarning)
            test_func(self, *args, **kwargs)
    return do_test


class TestGfe(unittest.TestCase):

    def setUp(self):
        self.data_dir = os.path.dirname(__file__) + "/resources"
        self.dblist = [str(i) + str(0) for i in range(315, 331)]
        expected_json = self.data_dir + "/expected.json"
        with open(expected_json) as json_data:
            self.expected = json.load(json_data)
            json_data.close()
        pass

    def test_001_class(self):
        gfe = GFE()
        self.assertIsInstance(gfe, GFE)
        pass

    @ignore_warnings
    def test_002_gfe(self):
        gfe = GFE()
        for ex in self.expected['gfe']:
            loc = ex['locus']
            ann = ex['annotation']
            exp = ex['gfe']
            annotation = {}
            for f in ann:
                seqrec = SeqRecord(seq=Seq(ann[f], generic_dna), id="002_gfe")
                annotation.update({f: seqrec})
            a = Annotation(annotation=annotation)
            features, gfe = gfe.get_gfe(a, loc)
            for feat in features:
                self.assertIsInstance(feat, Feature)
            self.assertEqual(gfe, exp)
        pass

    def test_003_store(self):
        gfe = GFE(store_features=True)
        self.assertTrue(gfe.store_features)
        pass
