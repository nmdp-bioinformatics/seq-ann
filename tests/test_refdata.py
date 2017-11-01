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
test_refdata
----------------------------------

Tests for `seqann.models.reference_data` module.
"""


import sys
import unittest

from Bio import SeqIO
from BioSQL import BioSeqDatabase
from seqann.models.reference_data import ReferenceData
from BioSQL.BioSeq import DBSeqRecord
from seqann.models.blast import Blast
from seqann.models.annotation import Annotation
import os
import pymysql
from Bio.SeqFeature import SeqFeature
import json


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


class TestRefdata(unittest.TestCase):

    def setUp(self):
        self.data_dir = os.path.dirname(__file__) + "/resources"
        self.dblist = [str(i) + str(0) for i in range(315, 329)]
        expected_json = self.data_dir + "/expected.json"
        with open(expected_json) as json_data:
            self.expected = json.load(json_data)
        pass

    def test_001_noserver(self):
        refdata = ReferenceData()
        self.assertIsInstance(refdata, ReferenceData)
        self.assertFalse(refdata.server_avail)
        self.assertGreater(len(refdata.hla_names), 10)
        self.assertEqual(refdata.structure_max['HLA-A'], 17)
        self.assertFalse(refdata.server_avail)
        self.assertEqual(refdata.dbversion, '3290')
        self.assertGreater(len(refdata.imgtdat), 0)
        pass

    @unittest.skipUnless(conn(), "TestRefdata 002 Requires MySQL connection")
    def test_002_server(self):

        server = BioSeqDatabase.open_database(driver="pymysql", user="root",
                                              passwd="", host="localhost",
                                              db="bioseqdb")
        refdata = ReferenceData(server=server)
        self.assertIsInstance(refdata, ReferenceData)
        self.assertTrue(refdata.server_avail)
        self.assertFalse(refdata.imgtdat)
        server.close()
        pass

    @unittest.skipUnless(conn(), "TestRefdata 003 Requires MySQL connection")
    def test_003_dblist(self):

        server = BioSeqDatabase.open_database(driver="pymysql", user="root",
                                              passwd="", host="localhost",
                                              db="bioseqdb")
        for db in self.dblist:
            refdata = ReferenceData(server=server, dbversion=db)
            self.assertEqual(refdata.dbversion, db)
            self.assertTrue(refdata.server_avail)
            self.assertFalse(refdata.imgtdat)
        server.close()
        pass

    @unittest.skipUnless(conn(), "TestRefdata 004 requires MySQL connection")
    def test_004_select(self):
        server = BioSeqDatabase.open_database(driver="pymysql", user="root",
                                              passwd="", host="localhost",
                                              db="bioseqdb")
        refdata = ReferenceData(server=server, dbversion='3290')
        input_seq = self.data_dir + '/exact_seqs.fasta'
        self.assertFalse(refdata.imgtdat)

        for ex in self.expected['exact']:
            i = int(ex['index'])
            locus = ex['locus']
            allele = ex['name']
            hla, loc = locus.split("-")
            in_seq = list(SeqIO.parse(input_seq, "fasta"))[i]
            annotation = refdata.search_refdata(in_seq, locus)
            self.assertIsNone(annotation.features)
            self.assertEqual(annotation.method, "match")
            self.assertIsInstance(annotation, Annotation)
            self.assertTrue(annotation.complete_annotation)
            self.assertGreater(len(annotation.annotation.keys()), 1)
            db = refdata.server[refdata.dbversion + "_" + loc]
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

            for feat in expected_seqs:
                if feat not in annotation.annotation:
                    self.assertEqual(feat, None)
                else:
                    self.assertEqual(str(expected_seqs[feat]),
                                     str(annotation.annotation[feat].seq))
        server.close()
        pass

    @unittest.skipUnless(conn(), "TestRefdata 005 requires MySQL connection")
    def test_005_dbtodat(self):
        server = BioSeqDatabase.open_database(driver="pymysql", user="root",
                                              passwd="", host="localhost",
                                              db="bioseqdb")
        refdata1 = ReferenceData(server=server, dbversion='3290')
        refdata2 = ReferenceData(dbversion='3290')

        datseqs = [a for a in refdata2.imgtdat
                   if a.description.split(",")[0] == 'HLA-A*01:01:01:01'][0]

        feats1 = [[feat.type, feat.extract(datseqs.seq)]
                  for feat in datseqs.features if feat.type != "source"
                  and feat.type != "CDS" and isinstance(feat, SeqFeature)]

        db = refdata1.server[refdata1.dbversion + "_" + 'A']
        expected = db.lookup(name='HLA-A*01:01:01:01')
        feats2 = [[feat.type, feat.extract(expected.seq)]
                  for feat in expected.features if feat.type != "source"
                  and feat.type != "CDS" and isinstance(feat, SeqFeature)]

        feat_types1 = {}
        expected_seqs1 = {}
        for i in range(0, len(feats1)):
            feat_name = ''
            if feats1[i][0] not in feat_types1:
                if(feats1[i][0] == "UTR"):
                    feat_name = "five_prime_UTR"
                    feat_types1.update({feats1[i][0]: 1})
                    expected_seqs1.update({feat_name: feats1[i][1]})
                else:
                    feat_name = feats1[i][0] + "_" + str(1)
                    feat_types1.update({feats1[i][0]: 1})
                    expected_seqs1.update({feat_name: feats1[i][1]})
            else:
                if(feats1[i][0] == "UTR"):
                    feat_name = "three_prime_UTR"
                    expected_seqs1.update({feat_name: feats1[i][1]})
                else:
                    num = feat_types1[feats1[i][0]] + 1
                    feat_name = feats1[i][0] + "_" + str(num)
                    feat_types1[feats1[i][0]] = num
                    expected_seqs1.update({feat_name: feats1[i][1]})

        feat_types2 = {}
        expected_seqs2 = {}
        for i in range(0, len(feats2)):
            feat_name = ''
            if feats2[i][0] not in feat_types2:
                if(feats2[i][0] == "UTR"):
                    feat_name = "five_prime_UTR"
                    feat_types2.update({feats2[i][0]: 1})
                    expected_seqs2.update({feat_name: feats2[i][1]})
                else:
                    feat_name = feats2[i][0] + "_" + str(1)
                    feat_types2.update({feats2[i][0]: 1})
                    expected_seqs2.update({feat_name: feats2[i][1]})
            else:
                if(feats2[i][0] == "UTR"):
                    feat_name = "three_prime_UTR"
                    expected_seqs2.update({feat_name: feats2[i][1]})
                else:
                    num = feat_types2[feats2[i][0]] + 1
                    feat_name = feats2[i][0] + "_" + str(num)
                    feat_types2[feats2[i][0]] = num
                    expected_seqs2.update({feat_name: feats2[i][1]})

        for feat in expected_seqs1:
            if feat not in expected_seqs2:
                self.assertEqual(feat, None)
            else:
                self.assertEqual(str(expected_seqs1[feat]),
                                 str(expected_seqs2[feat]))
        server.close()
        pass







