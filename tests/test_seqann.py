#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    gfe GFE.
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
import sys
import unittest
from Bio import SeqIO
from BioSQL import BioSeqDatabase
from seqann.models.reference_data import ReferenceData
from BioSQL.BioSeq import DBSeqRecord
from seqann.sequence_annotation import BioSeqAnn
from seqann.models.annotation import Annotation
from Bio.SeqFeature import SeqFeature
import os
import pymysql
import json


# TODO:
#   - add ability to detect insertions at beginning or end
#     of exons and introns.
#   - 

def conn():
    try:
        conn = pymysql.connect(host='localhost',
                               port=3306, user='root',
                               passwd='', db='bioseqdb')
        conn.close()
        return True
    except Exception as e:
        print("Exception while checking MYSQL Connection:" + str(e))
        return False


class TestBioSeqAnn(unittest.TestCase):

    def setUp(self):
        self.data_dir = os.path.dirname(__file__) + "/resources"
        self.dblist = [str(i) + str(0) for i in range(315, 329)]
        expected_json = self.data_dir + "/expected.json"
        with open(expected_json) as json_data:
            self.expected = json.load(json_data)
        pass

    @unittest.skipUnless(conn(), "TestBioSeqAnn 001 Requires MySQL connection")
    def test_001_seqann(self):
        server = BioSeqDatabase.open_database(driver="pymysql", user="root",
                                              passwd="", host="localhost",
                                              db="bioseqdb")
        seqann = BioSeqAnn(server=server)
        self.assertIsInstance(seqann, BioSeqAnn)
        self.assertIsInstance(seqann.refdata, ReferenceData)
        self.assertIsInstance(seqann.refdata, ReferenceData)
        self.assertGreater(len(seqann.refdata.hla_names), 10)
        self.assertEqual(seqann.refdata.structure_max['HLA-A'], 17)
        self.assertTrue(seqann.refdata.server_avail)
        self.assertEqual(seqann.refdata.dbversion, '3290')
        pass

    def test_002_noserver(self):
        seqann = BioSeqAnn()
        self.assertIsInstance(seqann, BioSeqAnn)
        self.assertIsInstance(seqann.refdata, ReferenceData)
        self.assertGreater(len(seqann.refdata.hla_names), 10)
        self.assertEqual(seqann.refdata.structure_max['HLA-A'], 17)
        self.assertFalse(seqann.refdata.server_avail)
        self.assertEqual(seqann.refdata.dbversion, '3290')
        self.assertGreater(len(seqann.refdata.imgtdat), 0)
        pass

    @unittest.skipUnless(conn(), "TestBioSeqAnn 003 requires MySQL connection")
    def test_003_ambig(self):
        server = BioSeqDatabase.open_database(driver="pymysql", user="root",
                                              passwd="", host="localhost",
                                              db="bioseqdb")
        seqann = BioSeqAnn(server=server)
        self.assertEqual(seqann.refdata.dbversion, '3290')
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
            self.assertIsInstance(ann, Annotation)
            self.assertTrue(ann.complete_annotation)
            self.assertGreater(len(ann.annotation.keys()), 1)
            db = seqann.refdata.server[seqann.refdata.dbversion + "_" + loc]
            expected = db.lookup(name=allele)
            feats = [[feat.type, feat.extract(expected.seq)]
                     for feat in expected.features if feat.type != "source"
                     and feat.type != "CDS" and isinstance(feat, SeqFeature)]

            n_diffs = 0
            feat_types = {}
            expected_seqs = {}
            for i in range(0, len(feats)):
                feat_name = ''
                if feats[i][0] not in feat_types:
                    if(feats[i][0] == "UTR"):
                        feat_name = "five_prime_UTR"
                        feat_types.update({feats[i][0]: 1})
                        expected_seqs.update({feat_name: feats[i][1]})
                    else:
                        feat_name = feats[i][0] + "_" + str(1)
                        feat_types.update({feats[i][0]: 1})
                        expected_seqs.update({feat_name: feats[i][1]})
                else:
                    if(feats[i][0] == "UTR"):
                        feat_name = "three_prime_UTR"
                        expected_seqs.update({feat_name: feats[i][1]})
                    else:
                        num = feat_types[feats[i][0]] + 1
                        feat_name = feats[i][0] + "_" + str(num)
                        feat_types[feats[i][0]] = num
                        expected_seqs.update({feat_name: feats[i][1]})

            self.assertGreater(len(expected_seqs.keys()), 1)
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

    @unittest.skipUnless(conn(), "TestBioSeqAnn 004 requires MySQL connection")
    def test_004_insertion(self):
        server = BioSeqDatabase.open_database(driver="pymysql", user="root",
                                              passwd="", host="localhost",
                                              db="bioseqdb")
        seqann = BioSeqAnn(server=server)
        self.assertEqual(seqann.refdata.dbversion, '3290')
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
            feats = [[feat.type, feat.extract(expected.seq)]
                     for feat in expected.features if feat.type != "source"
                     and feat.type != "CDS" and isinstance(feat, SeqFeature)]

            n_diffs = 0
            feat_types = {}
            expected_seqs = {}
            for i in range(0, len(feats)):
                feat_name = ''
                if feats[i][0] not in feat_types:
                    if(feats[i][0] == "UTR"):
                        feat_name = "five_prime_UTR"
                        feat_types.update({feats[i][0]: 1})
                        expected_seqs.update({feat_name: feats[i][1]})
                    else:
                        feat_name = feats[i][0] + "_" + str(1)
                        feat_types.update({feats[i][0]: 1})
                        expected_seqs.update({feat_name: feats[i][1]})
                else:
                    if(feats[i][0] == "UTR"):
                        feat_name = "three_prime_UTR"
                        expected_seqs.update({feat_name: feats[i][1]})
                    else:
                        num = feat_types[feats[i][0]] + 1
                        feat_name = feats[i][0] + "_" + str(num)
                        feat_types[feats[i][0]] = num
                        expected_seqs.update({feat_name: feats[i][1]})

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

    @unittest.skipUnless(conn(), "TestBioSeqAnn 004 requires MySQL connection")
    def test_005_partial(self):
        server = BioSeqDatabase.open_database(driver="pymysql", user="root",
                                              passwd="", host="localhost",
                                              db="bioseqdb")
        seqann = BioSeqAnn(server=server)
        self.assertEqual(seqann.refdata.dbversion, '3290')
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
            feats = [[feat.type, feat.extract(expected.seq)]
                     for feat in expected.features if feat.type != "source"
                     and feat.type != "CDS" and isinstance(feat, SeqFeature)]

            feat_types = {}
            expected_seqs = {}
            for i in range(0, len(feats)):
                feat_name = ''
                if feats[i][0] not in feat_types:
                    if(feats[i][0] == "UTR"):
                        feat_name = "five_prime_UTR"
                        feat_types.update({feats[i][0]: 1})
                        expected_seqs.update({feat_name: feats[i][1]})
                    else:
                        feat_name = feats[i][0] + "_" + str(1)
                        feat_types.update({feats[i][0]: 1})
                        expected_seqs.update({feat_name: feats[i][1]})
                else:
                    if(feats[i][0] == "UTR"):
                        feat_name = "three_prime_UTR"
                        expected_seqs.update({feat_name: feats[i][1]})
                    else:
                        num = feat_types[feats[i][0]] + 1
                        feat_name = feats[i][0] + "_" + str(num)
                        feat_types[feats[i][0]] = num
                        expected_seqs.update({feat_name: feats[i][1]})

            self.assertGreater(len(expected_seqs.keys()), 1)
            self.assertGreater(len(ann.annotation.keys()), 1)

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

    @unittest.skipUnless(conn(), "TestBioSeqAnn 006 requires MySQL connection")
    def test_006_partialambig(self):
        server = BioSeqDatabase.open_database(driver="pymysql", user="root",
                                              passwd="", host="localhost",
                                              db="bioseqdb")
        seqann = BioSeqAnn(server=server)
        self.assertEqual(seqann.refdata.dbversion, '3290')
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
            db = seqann.refdata.server[seqann.refdata.dbversion + "_" + loc]
            expected = db.lookup(name=allele)
            feats = [[feat.type, feat.extract(expected.seq)]
                     for feat in expected.features if feat.type != "source"
                     and feat.type != "CDS" and isinstance(feat, SeqFeature)]

            feat_types = {}
            expected_seqs = {}
            for i in range(0, len(feats)):
                feat_name = ''
                if feats[i][0] not in feat_types:
                    if(feats[i][0] == "UTR"):
                        feat_name = "five_prime_UTR"
                        feat_types.update({feats[i][0]: 1})
                        expected_seqs.update({feat_name: feats[i][1]})
                    else:
                        feat_name = feats[i][0] + "_" + str(1)
                        feat_types.update({feats[i][0]: 1})
                        expected_seqs.update({feat_name: feats[i][1]})
                else:
                    if(feats[i][0] == "UTR"):
                        feat_name = "three_prime_UTR"
                        expected_seqs.update({feat_name: feats[i][1]})
                    else:
                        num = feat_types[feats[i][0]] + 1
                        feat_name = feats[i][0] + "_" + str(num)
                        feat_types[feats[i][0]] = num
                        expected_seqs.update({feat_name: feats[i][1]})

            self.assertGreater(len(expected_seqs.keys()), 1)
            self.assertGreater(len(ann.annotation.keys()), 1)

            # Make sure only mapped feats exist
            for mf in ex['missing_feats']:
                self.assertFalse(mf in ann.annotation)

            for feat in ex['feats']:
                if feat in ex['diff']:
                    print(feat)
                    self.assertNotEqual(str(expected_seqs[feat]),
                                        str(ann.annotation[feat].seq))
                else:
                    self.assertEqual(str(expected_seqs[feat]),
                                     str(ann.annotation[feat].seq))

        server.close()
        pass

    @unittest.skipUnless(conn(), "TestBioSeqAnn 005 requires MySQL connection")
    def test_005_exact(self):
        server = BioSeqDatabase.open_database(driver="pymysql", user="root",
                                              passwd="", host="localhost",
                                              db="bioseqdb")
        seqann = BioSeqAnn(server=server, dbversion='3290')
        input_seq = self.data_dir + '/exact_seqs.fasta'

        for ex in self.expected['exact']:
            i = int(ex['index'])
            locus = ex['locus']
            allele = ex['name']
            hla, loc = locus.split("-")
            in_seq = list(SeqIO.parse(input_seq, "fasta"))[i]
            annotation = seqann.annotate(in_seq, locus)
            self.assertIsNone(annotation.features)
            self.assertEqual(annotation.method, "match")
            self.assertIsInstance(annotation, Annotation)
            self.assertTrue(annotation.complete_annotation)
            self.assertGreater(len(annotation.annotation.keys()), 1)
            db = seqann.refdata.server[seqann.refdata.dbversion + "_" + loc]
            expected = db.lookup(name=allele)
            feats = [[feat.type, feat.extract(expected.seq)]
                     for feat in expected.features if feat.type != "source"
                     and feat.type != "CDS" and isinstance(feat, SeqFeature)]

            feat_types = {}
            expected_seqs = {}
            for i in range(0, len(feats)):
                feat_name = ''
                if feats[i][0] not in feat_types:
                    if(feats[i][0] == "UTR"):
                        feat_name = "five_prime_UTR"
                        feat_types.update({feats[i][0]: 1})
                        expected_seqs.update({feat_name: feats[i][1]})
                    else:
                        feat_name = feats[i][0] + "_" + str(1)
                        feat_types.update({feats[i][0]: 1})
                        expected_seqs.update({feat_name: feats[i][1]})
                else:
                    if(feats[i][0] == "UTR"):
                        feat_name = "three_prime_UTR"
                        expected_seqs.update({feat_name: feats[i][1]})
                    else:
                        num = feat_types[feats[i][0]] + 1
                        feat_name = feats[i][0] + "_" + str(num)
                        feat_types[feats[i][0]] = num
                        expected_seqs.update({feat_name: feats[i][1]})

            self.assertGreater(len(expected_seqs.keys()), 1)
            self.assertGreater(len(annotation.annotation.keys()), 1)
            for feat in expected_seqs:
                if feat not in annotation.annotation:
                    self.assertEqual(feat, None)
                else:
                    self.assertEqual(str(expected_seqs[feat]),
                                     str(annotation.annotation[feat].seq))
        server.close()
        pass

    @unittest.skipUnless(not conn(), "TestBioSeqAnn 006 is run with no DB")
    def test_006_exactnoserv(self):
        seqann = BioSeqAnn(dbversion='3290')
        input_seq = self.data_dir + '/exact_seqs.fasta'

        for ex in self.expected['exact']:
            i = int(ex['index'])
            locus = ex['locus']
            allele = ex['name']
            hla, loc = locus.split("-")
            in_seq = list(SeqIO.parse(input_seq, "fasta"))[i]
            annotation = seqann.annotate(in_seq, locus)
            self.assertEqual(annotation.method, "nt_search")
            self.assertIsInstance(annotation, Annotation)
            self.assertTrue(annotation.complete_annotation)
            self.assertGreater(len(annotation.annotation.keys()), 1)

            expected = [a for a in seqann.refdata.imgtdat
                        if a.description.split(",")[0] == allele][0]

            feats = [[feat.type, feat.extract(expected.seq)]
                     for feat in expected.features if feat.type != "source"
                     and feat.type != "CDS" and isinstance(feat, SeqFeature)]

            feat_types = {}
            expected_seqs = {}
            for i in range(0, len(feats)):
                feat_name = ''
                if feats[i][0] not in feat_types:
                    if(feats[i][0] == "UTR"):
                        feat_name = "five_prime_UTR"
                        feat_types.update({feats[i][0]: 1})
                        expected_seqs.update({feat_name: feats[i][1]})
                    else:
                        feat_name = feats[i][0] + "_" + str(1)
                        feat_types.update({feats[i][0]: 1})
                        expected_seqs.update({feat_name: feats[i][1]})
                else:
                    if(feats[i][0] == "UTR"):
                        feat_name = "three_prime_UTR"
                        expected_seqs.update({feat_name: feats[i][1]})
                    else:
                        num = feat_types[feats[i][0]] + 1
                        feat_name = feats[i][0] + "_" + str(num)
                        feat_types[feats[i][0]] = num
                        expected_seqs.update({feat_name: feats[i][1]})

            self.assertGreater(len(expected_seqs.keys()), 1)
            self.assertGreater(len(annotation.annotation.keys()), 1)
            for feat in expected_seqs:
                if feat not in annotation.annotation:
                    self.assertEqual(feat, None)
                else:
                    self.assertEqual(str(expected_seqs[feat]),
                                     str(annotation.annotation[feat].seq))
        pass

    @unittest.skipUnless(conn(), "TestBioSeqAnn 007 Requires MySQL connection")
    def test_007_nomatch(self):
        server = BioSeqDatabase.open_database(driver="pymysql", user="root",
                                              passwd="", host="localhost",
                                              db="bioseqdb")
        seqann = BioSeqAnn(server=server, dbversion='3290')
        self.assertIsInstance(seqann, BioSeqAnn)
        input_seq = self.data_dir + '/classi_seq5-nomatch.fasta'
        in_seq = list(SeqIO.parse(input_seq, "fasta"))[0]
        annotation = seqann.annotate(in_seq, "HLA-A")
        self.assertIsInstance(annotation, Annotation)
        self.assertGreater(len(annotation.annotation.keys()), 1)
        self.assertTrue(annotation.complete_annotation)
        server.close()
        pass




