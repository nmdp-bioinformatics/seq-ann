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
import time
import pymysql
import unittest
import warnings
from Bio import SeqIO
from BioSQL import BioSeqDatabase

import csv
import glob
import pickle
import collections

from seqann.util import get_features
from seqann.models.annotation import Annotation
from seqann.sequence_annotation import BioSeqAnn
from seqann.models.reference_data import ReferenceData
from seqann.feature_client.models.feature import Feature

neo4jpass = 'gfedb'
if os.getenv("NEO4JPASS"):
    neo4jpass = os.getenv("NEO4JPASS")

neo4juser = 'neo4j'
if os.getenv("NEO4JUSER"):
    neo4juser = os.getenv("NEO4JUSER")

neo4jurl = "http://neo4j.b12x.org:80"
if os.getenv("NEO4JURL"):
    neo4jurl = os.getenv("NEO4JURL")

biosqlpass = "my-secret-pw"
if os.getenv("BIOSQLPASS"):
    biosqlpass = os.getenv("BIOSQLPASS")

biosqluser = 'root'
if os.getenv("BIOSQLUSER"):
    biosqluser = os.getenv("BIOSQLUSER")

biosqlhost = "localhost"
if os.getenv("BIOSQLHOST"):
    biosqlhost = os.getenv("BIOSQLHOST")

biosqldb = "bioseqdb"
if os.getenv("BIOSQLDB"):
    biosqldb = os.getenv("BIOSQLDB")

biosqlport = 3307
if os.getenv("BIOSQLPORT"):
    biosqlport = int(os.getenv("BIOSQLPORT"))

# import logging
# logging.basicConfig(format='%(asctime)s - %(name)-35s - %(levelname)-5s - %(funcName)s %(lineno)d: - %(message)s',
#                     datefmt='%m/%d/%Y %I:%M:%S %p',
#                     level=logging.INFO)

verbose = False
# if os.getenv("VERBOSE"):
#     if os.getenv("VERBOSE") == "True" \
#             or int(os.getenv("VERBOSE")) == 1:
#         import logging
#         logging.basicConfig(format='%(asctime)s - %(name)-35s - %(levelname)-5s - %(funcName)s %(lineno)d: - %(message)s',
#                             datefmt='%m/%d/%Y %I:%M:%S %p',
#                             level=logging.INFO)
#         verbose = True

verbosity = 0
if os.getenv("VERBOSITY"):
    verbosity = int(os.getenv("VERBOSITY"))


def ignore_warnings(test_func):
    def do_test(self, *args, **kwargs):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ResourceWarning)
            test_func(self, *args, **kwargs)
    return do_test


def conn():
    try:
        conn = pymysql.connect(host=biosqlhost,
                               port=biosqlport, user=biosqluser,
                               passwd=biosqlpass, db=biosqldb)
        conn.close()
        return True
    except Exception as e:
        print("Exception while checking MYSQL Connection:" + str(e))
        return False


class TestBioSeqAnn(unittest.TestCase):

    def setUp(self):
        self.data_dir = os.path.dirname(__file__) + "/resources"
        self.dblist = [str(i) + str(0) for i in range(315, 331)]
        expected_json = self.data_dir + "/expected.json"
        with open(expected_json) as json_data:
            self.expected = json.load(json_data)
            json_data.close()
        pass

    @unittest.skipUnless(conn(), "TestBioSeqAnn 001 Requires MySQL connection")
    def test_001_seqann(self):
        server = BioSeqDatabase.open_database(driver="pymysql",
                                              user=biosqluser,
                                              passwd=biosqlpass,
                                              host=biosqlhost,
                                              db=biosqldb,
                                              port=biosqlport)
        seqann = BioSeqAnn(server=server,
                           verbose=False,
                           verbosity=verbosity,
                           pid="001_seqann")
        self.assertIsInstance(seqann, BioSeqAnn)
        self.assertIsInstance(seqann.refdata, ReferenceData)
        self.assertIsInstance(seqann.refdata, ReferenceData)
        self.assertGreater(len(seqann.refdata.hla_names), 10)
        self.assertEqual(seqann.refdata.structure_max['HLA-A'], 17)
        self.assertTrue(seqann.refdata.server_avail)
        server.close()
        pass

    def test_002_noserver(self):
        seqann = BioSeqAnn(verbose=False,
                           verbosity=verbosity,
                           pid="002_noserver")
        self.assertIsInstance(seqann, BioSeqAnn)
        self.assertIsInstance(seqann.refdata, ReferenceData)
        self.assertGreater(len(seqann.refdata.hla_names), 10)
        self.assertEqual(seqann.refdata.structure_max['HLA-A'], 17)
        self.assertFalse(seqann.refdata.server_avail)
        self.assertGreater(len(seqann.refdata.seqref), 0)
        self.assertGreater(len(seqann.refdata.hlaref), 0)
        pass

    @ignore_warnings
    @unittest.skipUnless(conn(), "TestBioSeqAnn 003 requires MySQL connection")
    def test_003_ambigserv(self):
        server = BioSeqDatabase.open_database(driver="pymysql",
                                              user=biosqluser,
                                              passwd=biosqlpass,
                                              host=biosqlhost,
                                              db=biosqldb,
                                              port=biosqlport)
        seqann = BioSeqAnn(server=server, verbose=False,
                           verbosity=verbosity, pid="003_ambig")
        input_seq = self.data_dir + '/ambig_seqs.fasta'

        for ex in self.expected['ambig']:
            i = int(ex['index'])
            locus = ex['locus']
            allele = ex['name']
            hla, loc = locus.split("-")
            in_seq = list(SeqIO.parse(input_seq, "fasta"))[i]
            ann = seqann.annotate(in_seq, locus)
            self.assertEqual(ann.method, "nt_search")
            self.assertFalse(ann.missing)
            self.assertFalse(ann.blocks)
            self.assertTrue(ann.structure)
            self.assertIsInstance(ann, Annotation)
            self.assertTrue(ann.complete_annotation)
            self.assertGreater(len(ann.annotation.keys()), 1)
            db = seqann.refdata.server[seqann.refdata.dbversion + "_" + loc]
            self.assertEqual(ann.gfe, ex['gfe'])
            n_diffs = 0
            expected = db.lookup(name=allele)
            expected_seqs = get_features(expected)
            self.assertGreater(len(expected_seqs.keys()), 1)
            self.assertGreater(len(ann.structure), 1)
            for feat in ann.structure:
                self.assertIsInstance(feat, Feature)
            for feat in expected_seqs:
                if feat not in ann.annotation:
                    self.assertEqual(feat, None)
                else:
                    if feat in ex['diff']:
                        n_diffs += 1
                        self.assertNotEqual(str(expected_seqs[feat]),
                                            str(ann.annotation[feat].seq))
                    else:
                        self.assertEqual(str(expected_seqs[feat]),
                                         str(ann.annotation[feat].seq))
            self.assertEqual(n_diffs, len(ex['diff']))
        server.close()
        pass

    def test_004_ambig(self):
        seqann = BioSeqAnn(verbose=False,
                           verbosity=verbosity,
                           pid="003_ambig")
        input_seq = self.data_dir + '/ambig_seqs.fasta'
        for ex in self.expected['ambig']:
            i = int(ex['index'])
            locus = ex['locus']
            allele = ex['name']
            hla, loc = locus.split("-")
            in_seq = list(SeqIO.parse(input_seq, "fasta"))[i]
            ann = seqann.annotate(in_seq, locus)
            self.assertEqual(ann.method, "nt_search")
            self.assertFalse(ann.missing)
            self.assertFalse(ann.blocks)
            self.assertTrue(ann.structure)
            self.assertIsInstance(ann, Annotation)
            self.assertTrue(ann.complete_annotation)
            self.assertGreater(len(ann.annotation.keys()), 1)
            self.assertEqual(ann.gfe, ex['gfe'])
            n_diffs = 0
            expected = seqann.refdata.hlaref[allele]
            expected_seqs = get_features(expected)
            self.assertGreater(len(expected_seqs.keys()), 1)
            self.assertGreater(len(ann.structure), 1)
            for feat in ann.structure:
                self.assertIsInstance(feat, Feature)
            for feat in expected_seqs:
                if feat not in ann.annotation:
                    self.assertEqual(feat, None)
                else:
                    if feat in ex['diff']:
                        n_diffs += 1
                        self.assertNotEqual(str(expected_seqs[feat]),
                                            str(ann.annotation[feat].seq))
                    else:
                        self.assertEqual(str(expected_seqs[feat]),
                                         str(ann.annotation[feat].seq))
            self.assertEqual(n_diffs, len(ex['diff']))
        pass

    @ignore_warnings
    @unittest.skipUnless(conn(), "TestBioSeqAnn 004 requires MySQL connection")
    def test_005_insertionserv(self):
        server = BioSeqDatabase.open_database(driver="pymysql",
                                              user=biosqluser,
                                              passwd=biosqlpass,
                                              host=biosqlhost,
                                              db=biosqldb,
                                              port=biosqlport)
        seqann = BioSeqAnn(server=server,
                           verbose=False,
                           verbosity=verbosity,
                           pid="004_insertion")
        input_seq = self.data_dir + '/insertion_seqs.fasta'
        for ex in self.expected['insertion']:
            i = int(ex['index'])
            locus = ex['locus']
            allele = ex['name']
            hla, loc = locus.split("-")
            in_seq = list(SeqIO.parse(input_seq, "fasta"))[i]
            ann = seqann.annotate(in_seq, locus)
            self.assertEqual(ann.method, "nt_search")
            self.assertFalse(ann.missing)
            self.assertFalse(ann.blocks)
            self.assertIsInstance(ann, Annotation)
            self.assertTrue(ann.complete_annotation)
            self.assertGreater(len(ann.annotation.keys()), 1)
            db = seqann.refdata.server[seqann.refdata.dbversion + "_" + loc]
            expected = db.lookup(name=allele)
            self.assertEqual(ann.gfe, ex['gfe'])
            self.assertGreater(len(ann.structure), 1)
            for feat in ann.structure:
                self.assertIsInstance(feat, Feature)
            n_diffs = 0
            expected_seqs = get_features(expected)
            self.assertGreater(len(expected_seqs.keys()), 1)
            for feat in expected_seqs:
                if feat not in ann.annotation:
                    self.assertEqual(feat, None)
                else:
                    if feat in ex['diff']:
                        n_diffs += 1
                        self.assertNotEqual(str(expected_seqs[feat]),
                                            str(ann.annotation[feat].seq))
                        diff_len = len(str(ann.annotation[feat].seq)) - \
                            len(str(expected_seqs[feat]))
                        self.assertEqual(diff_len, ex['lengths'][feat])
                    else:
                        self.assertEqual(str(expected_seqs[feat]),
                                         str(ann.annotation[feat].seq))
            self.assertEqual(n_diffs, len(ex['diff']))
        server.close()
        pass

    def test_006_insertion(self):
        seqann = BioSeqAnn(verbosity=verbosity, pid="004_insertion")
        input_seq = self.data_dir + '/insertion_seqs.fasta'
        for ex in self.expected['insertion']:
            i = int(ex['index'])
            locus = ex['locus']
            allele = ex['name']
            hla, loc = locus.split("-")
            in_seq = list(SeqIO.parse(input_seq, "fasta"))[i]
            ann = seqann.annotate(in_seq, locus)
            self.assertEqual(ann.method, "nt_search")
            self.assertFalse(ann.missing)
            self.assertFalse(ann.blocks)
            self.assertIsInstance(ann, Annotation)
            self.assertTrue(ann.complete_annotation)
            self.assertGreater(len(ann.annotation.keys()), 1)
            expected = seqann.refdata.hlaref[allele]
            self.assertEqual(ann.gfe, ex['gfe'])
            self.assertGreater(len(ann.structure), 1)
            for feat in ann.structure:
                self.assertIsInstance(feat, Feature)
            n_diffs = 0
            expected_seqs = get_features(expected)
            self.assertGreater(len(expected_seqs.keys()), 1)
            for feat in expected_seqs:
                if feat not in ann.annotation:
                    self.assertEqual(feat, None)
                else:
                    if feat in ex['diff']:
                        n_diffs += 1
                        self.assertNotEqual(str(expected_seqs[feat]),
                                            str(ann.annotation[feat].seq))
                        diff_len = len(str(ann.annotation[feat].seq)) - \
                            len(str(expected_seqs[feat]))
                        self.assertEqual(diff_len, ex['lengths'][feat])
                    else:
                        self.assertEqual(str(expected_seqs[feat]),
                                         str(ann.annotation[feat].seq))
            self.assertEqual(n_diffs, len(ex['diff']))
        pass

    @ignore_warnings
    @unittest.skipUnless(conn(), "TestBioSeqAnn 005 requires MySQL connection")
    def test_007_partialserv(self):
        server = BioSeqDatabase.open_database(driver="pymysql",
                                              user=biosqluser,
                                              passwd=biosqlpass,
                                              host=biosqlhost,
                                              db=biosqldb,
                                              port=biosqlport)
        seqann = BioSeqAnn(server=server,
                           verbose=False,
                           verbosity=verbosity,
                           pid="005_partial")
        input_seq = self.data_dir + '/partial_seqs.fasta'
        for ex in self.expected['partial']:
            i = int(ex['index'])
            locus = ex['locus']
            allele = ex['name']
            hla, loc = locus.split("-")
            in_seq = list(SeqIO.parse(input_seq, "fasta"))[i]
            ann = seqann.annotate(in_seq, locus)
            self.assertTrue(ann.complete_annotation)
            self.assertEqual(ann.method, "nt_search")
            self.assertFalse(ann.blocks)
            self.assertIsInstance(ann, Annotation)
            self.assertTrue(ann.complete_annotation)
            self.assertGreater(len(ann.annotation.keys()), 1)
            db = seqann.refdata.server[seqann.refdata.dbversion + "_" + loc]
            expected = db.lookup(name=allele)
            expected_seqs = get_features(expected)
            self.assertGreater(len(expected_seqs.keys()), 1)
            self.assertGreater(len(ann.annotation.keys()), 1)
            self.assertEqual(ann.gfe, ex['gfe'])
            self.assertGreater(len(ann.structure), 1)
            for feat in ann.structure:
                self.assertIsInstance(feat, Feature)
            # Make sure only mapped feats exist
            for mf in ex['missing_feats']:
                self.assertFalse(mf in ann.annotation)

            # Check that expected feats are mapped
            for feat in ex['feats']:
                if feat not in ann.annotation:
                    self.assertEqual(feat, None)
                else:
                    self.assertEqual(str(expected_seqs[feat]),
                                     str(ann.annotation[feat].seq))

        server.close()
        pass

    def test_008_partial(self):
        seqann = BioSeqAnn(verbose=False,
                           verbosity=verbosity,
                           pid="005_partial")
        input_seq = self.data_dir + '/partial_seqs.fasta'
        for ex in self.expected['partial']:
            i = int(ex['index'])
            locus = ex['locus']
            allele = ex['name']
            hla, loc = locus.split("-")
            in_seq = list(SeqIO.parse(input_seq, "fasta"))[i]
            ann = seqann.annotate(in_seq, locus)
            self.assertTrue(ann.complete_annotation)
            self.assertEqual(ann.method, "nt_search")
            self.assertFalse(ann.blocks)
            self.assertIsInstance(ann, Annotation)
            self.assertTrue(ann.complete_annotation)
            self.assertGreater(len(ann.annotation.keys()), 1)
            expected = seqann.refdata.hlaref[allele]
            expected_seqs = get_features(expected)
            self.assertGreater(len(expected_seqs.keys()), 1)
            self.assertGreater(len(ann.annotation.keys()), 1)
            self.assertEqual(ann.gfe, ex['gfe'])
            self.assertGreater(len(ann.structure), 1)
            for feat in ann.structure:
                self.assertIsInstance(feat, Feature)
            # Make sure only mapped feats exist
            for mf in ex['missing_feats']:
                self.assertFalse(mf in ann.annotation)

            # Check that expected feats are mapped
            for feat in ex['feats']:
                if feat not in ann.annotation:
                    self.assertEqual(feat, None)
                else:
                    self.assertEqual(str(expected_seqs[feat]),
                                     str(ann.annotation[feat].seq))

        pass

    @ignore_warnings
    @unittest.skipUnless(conn(), "TestBioSeqAnn 009 requires MySQL connection")
    def test_009_partialambigserv(self):
        server = BioSeqDatabase.open_database(driver="pymysql",
                                              user=biosqluser,
                                              passwd=biosqlpass,
                                              host=biosqlhost,
                                              db=biosqldb,
                                              port=biosqlport)
        seqann = BioSeqAnn(server=server,
                           verbose=False,
                           verbosity=verbosity,
                           pid="006_partialambig")
        input_seq = self.data_dir + '/partial_ambig.fasta'

        for ex in self.expected['partial_ambig']:
            i = int(ex['index'])
            locus = ex['locus']
            allele = ex['name']
            hla, loc = locus.split("-")
            print(str(i), allele)
            in_seq = list(SeqIO.parse(input_seq, "fasta"))[i]
            ann = seqann.annotate(in_seq, locus)
            self.assertTrue(ann.complete_annotation)
            self.assertEqual(ann.method, ex['method'])
            self.assertFalse(ann.blocks)
            self.assertIsInstance(ann, Annotation)
            self.assertTrue(ann.complete_annotation)
            self.assertGreater(len(ann.annotation.keys()), 1)
            db = seqann.refdata.server[seqann.refdata.dbversion + "_" + loc]
            expected = db.lookup(name=allele)
            expected_seqs = get_features(expected)
            self.assertGreater(len(expected_seqs.keys()), 1)
            self.assertGreater(len(ann.annotation.keys()), 1)
            self.assertEqual(ann.gfe, ex['gfe'])

            self.assertGreater(len(ann.structure), 1)
            for feat in ann.structure:
                self.assertIsInstance(feat, Feature)
            # Make sure only mapped feats exist
            for mf in ex['missing_feats']:
                self.assertFalse(mf in ann.annotation)

            for feat in ex['feats']:
                if feat in ex['diff']:
                    self.assertNotEqual(str(expected_seqs[feat]),
                                        str(ann.annotation[feat].seq))
                else:
                    self.assertEqual(str(expected_seqs[feat]),
                                     str(ann.annotation[feat].seq))

        server.close()
        pass

    def test_010_partialambig(self):
        seqann = BioSeqAnn(verbose=False,
                           verbosity=verbosity,
                           pid="006_partialambig")
        input_seq = self.data_dir + '/partial_ambig.fasta'
        for ex in self.expected['partial_ambig']:
            i = int(ex['index'])
            locus = ex['locus']
            allele = ex['name']
            hla, loc = locus.split("-")
            in_seq = list(SeqIO.parse(input_seq, "fasta"))[i]
            ann = seqann.annotate(in_seq, locus)
            self.assertTrue(ann.complete_annotation)
            self.assertEqual(ann.method, ex['method'])
            self.assertFalse(ann.blocks)
            self.assertIsInstance(ann, Annotation)
            self.assertTrue(ann.complete_annotation)
            self.assertGreater(len(ann.annotation.keys()), 1)
            expected = seqann.refdata.hlaref[allele]
            expected_seqs = get_features(expected)
            self.assertGreater(len(expected_seqs.keys()), 1)
            self.assertGreater(len(ann.annotation.keys()), 1)
            self.assertEqual(ann.gfe, ex['gfe'])
            self.assertGreater(len(ann.structure), 1)
            for feat in ann.structure:
                self.assertIsInstance(feat, Feature)
            # Make sure only mapped feats exist
            for mf in ex['missing_feats']:
                self.assertFalse(mf in ann.annotation)

            for feat in ex['feats']:
                if feat in ex['diff']:
                    self.assertNotEqual(str(expected_seqs[feat]),
                                        str(ann.annotation[feat].seq))
                else:
                    self.assertEqual(str(expected_seqs[feat]),
                                     str(ann.annotation[feat].seq))

        pass

    @ignore_warnings
    @unittest.skipUnless(conn(), "TestBioSeqAnn 011 requires MySQL connection")
    def test_011_exactserv(self):
        server = BioSeqDatabase.open_database(driver="pymysql",
                                              user=biosqluser,
                                              passwd=biosqlpass,
                                              host=biosqlhost,
                                              db=biosqldb,
                                              port=biosqlport)
        seqann = BioSeqAnn(server=server,
                           verbose=False,
                           verbosity=verbosity,
                           pid="007_exact")
        input_seq = self.data_dir + '/exact_seqs.fasta'

        for ex in self.expected['exact']:
            i = int(ex['index'])
            locus = ex['locus']
            allele = ex['name']
            hla, loc = locus.split("-")
            in_seq = list(SeqIO.parse(input_seq, "fasta"))[i]
            annotation = seqann.annotate(in_seq, locus)
            self.assertTrue(annotation.exact)
            self.assertIsNone(annotation.features)
            self.assertEqual(annotation.method, "match")
            self.assertIsInstance(annotation, Annotation)
            self.assertTrue(annotation.complete_annotation)
            self.assertGreater(len(annotation.annotation.keys()), 1)
            db = seqann.refdata.server[seqann.refdata.dbversion + "_" + loc]
            expected = db.lookup(name=allele)
            expected_seqs = get_features(expected)
            self.assertGreater(len(annotation.structure), 1)
            for feat in annotation.structure:
                self.assertIsInstance(feat, Feature)
            self.assertEqual(annotation.gfe, ex['gfe'])
            self.assertGreater(len(expected_seqs.keys()), 1)
            self.assertGreater(len(annotation.annotation.keys()), 1)
            for feat in expected_seqs:
                if feat not in annotation.annotation:
                    self.assertEqual(feat, None)
                else:
                    self.assertEqual(str(expected_seqs[feat]),
                                     str(annotation.annotation[feat]))
        server.close()
        pass

    def test_012_exact(self):
        seqann = BioSeqAnn(verbose=False,
                           verbosity=verbosity,
                           pid="007_exact")
        input_seq = self.data_dir + '/exact_seqs.fasta'
        for ex in self.expected['exact']:
            i = int(ex['index'])
            locus = ex['locus']
            allele = ex['name']
            hla, loc = locus.split("-")
            in_seq = list(SeqIO.parse(input_seq, "fasta"))[i]
            annotation = seqann.annotate(in_seq, locus)
            self.assertTrue(annotation.exact)
            self.assertIsNone(annotation.features)
            self.assertEqual(annotation.method, "match")
            self.assertIsInstance(annotation, Annotation)
            self.assertTrue(annotation.complete_annotation)
            self.assertGreater(len(annotation.annotation.keys()), 1)
            expected = seqann.refdata.hlaref[allele]
            expected_seqs = get_features(expected)
            self.assertGreater(len(annotation.structure), 1)
            for feat in annotation.structure:
                self.assertIsInstance(feat, Feature)
            self.assertEqual(annotation.gfe, ex['gfe'])
            self.assertGreater(len(expected_seqs.keys()), 1)
            self.assertGreater(len(annotation.annotation.keys()), 1)
            for feat in expected_seqs:
                if feat not in annotation.annotation:
                    self.assertEqual(feat, None)
                else:
                    self.assertEqual(str(expected_seqs[feat]),
                                     str(annotation.annotation[feat]))
        pass

    @ignore_warnings
    @unittest.skipUnless(conn(), "TestBioSeqAnn 013 Requires MySQL connection")
    def test_013_nomatch(self):
        server = BioSeqDatabase.open_database(driver="pymysql",
                                              user=biosqluser,
                                              passwd=biosqlpass,
                                              host=biosqlhost,
                                              db=biosqldb,
                                              port=biosqlport)
        seqann = BioSeqAnn(server=server,
                           verbose=False,
                           verbosity=verbosity,
                           pid="009_nomatch")
        self.assertIsInstance(seqann, BioSeqAnn)
        input_seq = self.data_dir + '/nomatch_seqs.fasta'
        in_seq = list(SeqIO.parse(input_seq, "fasta"))[0]
        annotation = seqann.annotate(in_seq, "HLA-A")
        self.assertIsInstance(annotation, Annotation)
        self.assertGreater(len(annotation.annotation.keys()), 1)
        self.assertTrue(annotation.complete_annotation)
        server.close()
        pass

    @ignore_warnings
    @unittest.skipUnless(conn(), "TestBioSeqAnn 014 Requires MySQL connection")
    def test_014_noloc(self):
        server = BioSeqDatabase.open_database(driver="pymysql",
                                              user=biosqluser,
                                              passwd=biosqlpass,
                                              host=biosqlhost,
                                              db=biosqldb,
                                              port=biosqlport)
        seqann = BioSeqAnn(server=server,
                           verbose=False,
                           verbosity=0,
                           pid="010_noloc")
        self.assertIsInstance(seqann, BioSeqAnn)
        input_seq = self.data_dir + '/nomatch_seqs.fasta'
        in_seq = list(SeqIO.parse(input_seq, "fasta"))[0]
        #print(in_seq)
        annotation = seqann.annotate(in_seq)
        self.assertIsInstance(annotation, Annotation)
        self.assertGreater(len(annotation.annotation.keys()), 1)
        self.assertTrue(annotation.complete_annotation)
        server.close()
        pass

    @unittest.skipUnless(conn(), "TestBioSeqAnn 011 Requires MySQL connection")
    def test_015_fail(self):
        input_seq = self.data_dir + '/failed_seqs.fasta'
        in_seq = list(SeqIO.parse(input_seq, "fasta"))[0]
        server = BioSeqDatabase.open_database(driver="pymysql",
                                              user=biosqluser,
                                              passwd=biosqlpass,
                                              host=biosqlhost,
                                              db=biosqldb,
                                              port=biosqlport)
        seqann = BioSeqAnn(server=server,
                           verbose=False,
                           verbosity=verbosity,
                           pid="011_fail")
        self.assertFalse(seqann.refdata.seqref)
        self.assertFalse(seqann.refdata.hlaref)
        annotation = seqann.annotate(in_seq)
        self.assertFalse(annotation)
        server.close()
        pass

    def test_016_debug(self):
        seqann = BioSeqAnn(debug={"seqann": 5, "align": 1,
                                  "seq_search": 3,
                                  "refdata": 2})
        self.assertTrue(seqann.debug)
        self.assertEqual(seqann.verbosity, 5)
        self.assertEqual(seqann.align_verbosity, 1)
        self.assertEqual(seqann.seqsearch.verbosity, 3)
        self.assertEqual(seqann.refdata.verbosity, 2)

        seqann = BioSeqAnn(debug={"seqann": 2,
                                  "seq_search": 5})
        self.assertTrue(seqann.debug)
        self.assertEqual(seqann.verbosity, 2)
        self.assertEqual(seqann.align_verbosity, 0)
        self.assertEqual(seqann.seqsearch.verbosity, 5)
        self.assertEqual(seqann.refdata.verbosity, 0)

        seqann = BioSeqAnn(debug={"gfe": 2,
                                  "seq_search": 5})
        self.assertTrue(seqann.debug)
        self.assertTrue(seqann.gfe.verbose)
        self.assertEqual(seqann.gfe.verbosity, 2)
        self.assertEqual(seqann.seqsearch.verbosity, 5)
        self.assertEqual(seqann.refdata.verbosity, 0)

        seqann = BioSeqAnn(verbose=True,
                           verbosity=3)
        self.assertFalse(seqann.debug)
        self.assertEqual(seqann.verbosity, 3)
        self.assertEqual(seqann.align_verbosity, 3)
        self.assertEqual(seqann.seqsearch.verbosity, 3)
        self.assertEqual(seqann.refdata.verbosity, 3)
        pass

    @ignore_warnings
    @unittest.skipUnless(conn(), "TestBioSeqAnn 013 Requires MySQL connection")
    def test_017_logging(self):
        server = BioSeqDatabase.open_database(driver="pymysql",
                                              user=biosqluser,
                                              passwd=biosqlpass,
                                              host=biosqlhost,
                                              db=biosqldb,
                                              port=biosqlport)

        with self.assertLogs(level='INFO') as cm:
            seqann = BioSeqAnn(server=server,
                               verbose=True)
            input_seq = self.data_dir + '/failed_seqs.fasta'
            in_seq = list(SeqIO.parse(input_seq, "fasta"))[0]
            annotation = seqann.annotate(in_seq)
            self.assertFalse(annotation)

        self.assertGreater(len(cm.output), 1)
        error = list(cm.output)[len(cm.output)-1].split(":")[0]
        error_msg = list(cm.output)[len(cm.output)-1].split("-")[1]
        self.assertEqual(error, "ERROR")
        self.assertEqual(error_msg, " Locus could not be determined!")
        server.close()
        pass

    @ignore_warnings
    @unittest.skipUnless(conn(), "TestBioSeqAnn 014 Requires MySQL connection")
    def test_018_nogfe(self):
        server = BioSeqDatabase.open_database(driver="pymysql",
                                              user=biosqluser,
                                              passwd=biosqlpass,
                                              host=biosqlhost,
                                              db=biosqldb,
                                              port=biosqlport)

        with self.assertLogs(level='INFO') as cm:
            seqann = BioSeqAnn(server=server, verbose=True)
            input_seq = self.data_dir + '/failed_seqs.fasta'
            in_seq = list(SeqIO.parse(input_seq, "fasta"))[1]
            annotation = seqann.annotate(in_seq)
            self.assertFalse(annotation.gfe)
            self.assertFalse(annotation.structure)
            self.assertTrue(annotation.annotation)

        self.assertGreater(len(cm.output), 2)
        error = list(cm.output)[0].split(":")[0]
        error_msg = list(cm.output)[0].split("-")[1]
        self.assertEqual(error, "WARNING")
        self.assertEqual(error_msg, " Sequence alphabet contains non DNA")
        server.close()
        pass

    # TODO: Break alleles into separate tests
    #       based on what they are testing
    # @ignore_warnings
    # @unittest.skipUnless(conn(), "TestBioSeqAnn 015 Requires MySQL connection")
    # def test_019_skipserv(self):
    #     # import logging
    #     # logging.basicConfig(format='%(asctime)s - %(name)-35s - %(levelname)-5s - %(funcName)s %(lineno)d: - %(message)s',
    #     #                    datefmt='%m/%d/%Y %I:%M:%S %p',
    #     #                   level=logging.INFO)
    #     server = BioSeqDatabase.open_database(driver="pymysql",
    #                                           user=biosqluser,
    #                                           passwd=biosqlpass,
    #                                           host=biosqlhost,
    #                                           db=biosqldb,
    #                                           port=biosqlport)

    #     seqann1 = BioSeqAnn(verbose=False, verbosity=0)
    #     seqann = BioSeqAnn(server=server, verbose=False, verbosity=0)
    #     refdata = seqann1.refdata

    #     # removed 'HLA-DRB1*04:04:01' because it's
    #     # # too large to test with travis
    #     # test_list = ['HLA-C*07:241', 'HLA-A*01:07', 'HLA-A*01:01:59',
    #     #              'HLA-A*01:09:01:01', 'HLA-A*02:545', 'HLA-B*40:02:05',
    #     #              'HLA-A*29:13', 'HLA-A*24:03:02', 'HLA-A*02:544',
    #     #              'HLA-DQA1*04:01:01:01', 'HLA-A*01:217', 'HLA-A*01:22N',
    #     #              'HLA-B*51:42', 'HLA-C*03:04:05', 'HLA-A*01:01:01:04',
    #     #              'HLA-A*01:09:01:01', 'HLA-B*82:01', 'HLA-A*03:04:01',
    #     #              'HLA-C*07:06:01:01', 'HLA-A*03:51', 'HLA-A*29:109',
    #     #              'HLA-A*02:01:130', 'HLA-B*07:271', "HLA-DRB1*13:247",
    #     #              "HLA-DRB4*01:03:05", "HLA-DRB4*01:03:06",
    #     #              "HLA-DQA1*05:01:02", "HLA-DRB1*13:02:02",
    #     #              "HLA-DRB4*01:03:04"]
    #     test_list = ['HLA-C*07:241', 'HLA-A*01:07', 'HLA-A*01:01:59',
    #                  'HLA-A*01:09:01:01', 'HLA-A*02:545', 'HLA-B*40:02:05',
    #                  'HLA-A*29:13', 'HLA-A*24:03:02', 'HLA-A*02:544',
    #                  'HLA-DQA1*04:01:01:01', 'HLA-A*01:217', 'HLA-A*01:22N',
    #                  'HLA-B*51:42', 'HLA-C*03:04:05', 'HLA-A*01:01:01:04',
    #                  'HLA-A*01:09:01:01', 'HLA-B*82:01', 'HLA-A*03:04:01',
    #                  'HLA-C*07:06:01:01', 'HLA-A*03:51', 'HLA-A*29:109',
    #                  'HLA-A*02:01:130', 'HLA-B*07:271', "HLA-DRB1*13:247",
    #                  "HLA-DRB4*01:03:05", "HLA-DRB4*01:03:06",
    #                  "HLA-DQA1*05:01:02", "HLA-DRB1*13:02:02",
    #                  "HLA-DRB4*01:03:04"]
    #     for seqname in refdata.hlaref:
    #         if seqname not in test_list:
    #             continue

    #         print(seqname)
    #         seqrec = refdata.hlaref[seqname]
    #         locus = seqrec.description.split("*")[0]
    #         ann1 = seqann.annotate(seqrec, locus=locus)
    #         ann2 = seqann.annotate(seqrec, locus=locus, skip=[seqname])
    #         self.assertTrue(ann1.exact)
    #         self.assertEqual(len(ann2.annotation), len(ann1.annotation))

    #         #for feat in ann2.structure:
    #         #    self.assertIsInstance(feat, Feature)

    #         for f in ann1.annotation:
    #             self.assertTrue(f in ann2.annotation)
    #             seq1 = str(ann1.annotation[f])
    #             # seq2 = '** NA **'
    #             # if f in ann2.annotation:
    #             # if f not in ann2.annotation:
    #             #     print(seqname, "MISSING", f)
    #             # else:
    #             # if f in ann2.annotation:
    #             seq2 = str(ann2.annotation[f].seq)
    #             #     if seq1 != seq2:
    #             #         print(seqname, "NOT EQUAL", f)
    #             #print(f, seq1, seq2)
    #             self.assertEqual(seq1, seq2)
    #             #print(f,seq2)

    #         self.assertEqual(ann1.gfe, ann2.gfe)

    #     server.close()
    #     pass

    def test_020_skip(self):
        # import logging
        # logging.basicConfig(format='%(asctime)s - %(name)-35s - %(levelname)-5s - %(funcName)s %(lineno)d: - %(message)s',
        #                     datefmt='%m/%d/%Y %I:%M:%S %p',
        #                     level=logging.INFO)
        seqann = BioSeqAnn(verbose=False)
        refdata = seqann.refdata
        test_list = ['HLA-C*07:241', 'HLA-A*01:07', 'HLA-A*01:01:59',
                     'HLA-A*01:09:01:01', 'HLA-A*02:545',
                     'HLA-A*29:13', 'HLA-A*24:03:02', 'HLA-A*02:544',
                     'HLA-DQA1*04:01:01:01', 'HLA-A*01:217', 'HLA-A*01:22N',
                     'HLA-B*51:42', 'HLA-C*03:04:05', 'HLA-A*01:01:01:04',
                     'HLA-A*01:09:01:01', 'HLA-B*82:01']

        for seqname in refdata.hlaref:
            if seqname not in test_list:
                continue
            print(seqname)
            seqrec = refdata.hlaref[seqname]
            locus = seqrec.description.split("*")[0]
            ann1 = seqann.annotate(seqrec, locus=locus)
            ann2 = seqann.annotate(seqrec, locus=locus, skip=[seqname])
            self.assertTrue(ann1.exact)
            self.assertEqual(len(ann2.annotation), len(ann1.annotation))
            self.assertEqual(ann1.gfe, ann2.gfe)
            self.assertGreater(len(ann2.structure), 1)
            for feat in ann2.structure:
                self.assertIsInstance(feat, Feature)
            for f in ann1.annotation:
                self.assertTrue(f in ann2.annotation)
                seq1 = str(ann1.annotation[f])
                seq2 = str(ann2.annotation[f].seq)
                self.assertEqual(seq1, seq2)
        pass

    @ignore_warnings
    @unittest.skipUnless(conn(), "TestBioSeqAnn 016 Requires MySQL connection")
    def test_021_stringseq(self):
        server = BioSeqDatabase.open_database(driver="pymysql",
                                              user=biosqluser,
                                              passwd=biosqlpass,
                                              host=biosqlhost,
                                              db=biosqldb,
                                              port=biosqlport)
        seqann = BioSeqAnn(server=server,
                           verbose=False,
                           verbosity=verbosity,
                           pid="015_stringseq")
        input_seq = self.data_dir + '/exact_seqs.fasta'
        ex = self.expected['exact'][0]
        locus = ex['locus']
        allele = ex['name']
        hla, loc = locus.split("-")
        in_seqrec = list(SeqIO.parse(input_seq, "fasta"))[0]
        in_str = str(in_seqrec.seq)
        in_seq = in_seqrec.seq
        ann_str = seqann.annotate(in_str, locus)
        ann_seq = seqann.annotate(in_seq, locus)
        for annotation in [ann_str, ann_seq]:
            self.assertTrue(annotation.exact)
            self.assertIsNone(annotation.features)
            self.assertEqual(annotation.method, "match")
            self.assertIsInstance(annotation, Annotation)
            self.assertTrue(annotation.complete_annotation)
            self.assertGreater(len(annotation.annotation.keys()), 1)
            db = seqann.refdata.server[seqann.refdata.dbversion + "_" + loc]
            expected = db.lookup(name=allele)
            expected_seqs = get_features(expected)
            self.assertEqual(annotation.gfe, ex['gfe'])
            self.assertGreater(len(expected_seqs.keys()), 1)
            self.assertGreater(len(annotation.annotation.keys()), 1)
            self.assertGreater(len(annotation.structure), 1)
            for feat in annotation.structure:
                self.assertIsInstance(feat, Feature)
            for feat in expected_seqs:
                if feat not in annotation.annotation:
                    self.assertEqual(feat, None)
                else:
                    self.assertEqual(str(expected_seqs[feat]),
                                     str(annotation.annotation[feat]))
        server.close()
        pass

