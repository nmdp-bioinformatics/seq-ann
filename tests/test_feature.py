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

from Bio import BiopythonExperimentalWarning
warnings.simplefilter("ignore", BiopythonExperimentalWarning)

from seqann.feature_client.api_client import ApiClient
from seqann.feature_client.rest import ApiException
from seqann.feature_client.models.feature import Feature
from seqann.feature_client.apis.features_api import FeaturesApi
from seqann.feature_client.models.feature_request import FeatureRequest


def ignore_warnings(test_func):
    def do_test(self, *args, **kwargs):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ResourceWarning)
            test_func(self, *args, **kwargs)
    return do_test


class TestFeature(unittest.TestCase):

    def setUp(self):
        self.data_dir = os.path.dirname(__file__) + "/resources"
        expected_json = self.data_dir + "/expected.json"
        with open(expected_json) as json_data:
            self.expected = json.load(json_data)
            json_data.close()
        pass

    def test_001_client(self):
        client = ApiClient(host="http://feature.nmdp-bioinformatics.org")
        self.assertIsInstance(client, ApiClient)
        self.assertEqual(client.host, "http://feature.nmdp-bioinformatics.org")
        pass

    def test_002_api(self):
        client = ApiClient(host="http://feature.nmdp-bioinformatics.org")
        api_instance = FeaturesApi(api_client=client)
        self.assertIsInstance(api_instance, FeaturesApi)
        pass

    def test_003_feature(self):
        feat = "exon"
        locus = "HLA-A"
        seq = "TGTGA"
        accession = 1
        feature = Feature(term=feat,
                          rank=8,
                          locus=locus,
                          sequence=seq,
                          accession=accession)
        self.assertIsInstance(feature, Feature)
        self.assertEqual(feature.term, feat)
        self.assertEqual(feature.rank, 8)
        self.assertEqual(feature.locus, locus)
        self.assertEqual(feature.sequence, seq)
        self.assertEqual(feature.accession, accession)

    def test_004_request(self):
        feat = "exon"
        locus = "HLA-A"
        seq = "TGTGA"
        request = FeatureRequest(locus=locus,
                                 term=feat,
                                 rank=8,
                                 sequence=seq)
        self.assertIsInstance(request, FeatureRequest)
        self.assertEqual(request.term, feat)
        self.assertEqual(request.rank, 8)
        self.assertEqual(request.locus, locus)
        self.assertEqual(request.sequence, seq)
        pass

    @ignore_warnings
    def test_005_create(self):
        feat = "exon"
        locus = "HLA-A"
        seq = "TGTGA"
        accession = 1
        client = ApiClient(host="http://feature.nmdp-bioinformatics.org")
        api = FeaturesApi(api_client=client)
        request = FeatureRequest(locus=locus,
                                 term=feat,
                                 rank=8,
                                 sequence=seq)
        feature = api.create_feature(body=request)
        self.assertIsInstance(feature, Feature)
        self.assertEqual(feature.accession, accession)
        pass

    @ignore_warnings
    def test_006_exception(self):
        client = ApiClient(host="http://feature.nmdp-bioinformatics.org")
        api = FeaturesApi(api_client=client)
        request = FeatureRequest(locus="HLA-Z")
        with self.assertRaises(ApiException) as context:
            api.create_feature(body=request)
        self.assertEqual(context.exception.reason, "Internal Server Error")
        pass

