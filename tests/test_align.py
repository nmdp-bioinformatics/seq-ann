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

Tests for `seqann.align` module.
"""


import sys
import unittest

from Bio import SeqIO
from BioSQL import BioSeqDatabase
from seqann.models.reference_data import ReferenceData
from seqann.blast_cmd import blastn
from BioSQL.BioSeq import DBSeqRecord
from seqann.sequence_annotation import BioSeqAnn
from Bio.SeqFeature import SeqFeature
from seqann.align import align_seqs
import os
import pymysql
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


class TestAlign(unittest.TestCase):

    def setUp(self):
        self.data_dir = os.path.dirname(__file__) + "/resources"
        expected_json = self.data_dir + "/expected.json"
        with open(expected_json) as json_data:
            self.expected = json.load(json_data)
        pass

    @unittest.skipUnless(conn(), "TestAlign 001 requires MySQL connection")
    def test_001_exact(self):
        input_seq = self.data_dir + '/exact_seqs.fasta'
        server = BioSeqDatabase.open_database(driver="pymysql", user="root",
                                              passwd="", host="localhost",
                                              db="bioseqdb")

        for ex in self.expected['exact']:
            i = int(ex['index'])
            locus = ex['locus']
            allele = ex['name']
            hla, loc = locus.split("-")
            in_seq = list(SeqIO.parse(input_seq, "fasta"))[i]
            db = server["3290_" + loc]
            refseq = db.lookup(name=allele)
            ann = align_seqs(refseq, in_seq, locus)
            feats = [[feat.type, feat.extract(refseq.seq)]
                     for feat in refseq.features if feat.type != "source"
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
            for feat in expected_seqs:
                if feat not in ann.annotation:
                    self.assertEqual(feat, None)
                else:
                    self.assertEqual(str(expected_seqs[feat]),
                                     str(ann.annotation[feat].seq))
        server.close()
        pass

    @unittest.skipUnless(conn(), "TestAlign 002 requires MySQL connection")
    def test_002_ambig(self):
        input_seq = self.data_dir + '/ambig_seqs.fasta'
        server = BioSeqDatabase.open_database(driver="pymysql", user="root",
                                              passwd="", host="localhost",
                                              db="bioseqdb")

        for ex in self.expected['ambig']:
            i = int(ex['index'])
            locus = ex['locus']
            allele = ex['name']
            hla, loc = locus.split("-")
            in_seq = list(SeqIO.parse(input_seq, "fasta"))[i]
            db = server["3290_" + loc]
            refseq = db.lookup(name=allele)
            ann = align_seqs(refseq, in_seq, locus)
            feats = [[feat.type, feat.extract(refseq.seq)]
                     for feat in refseq.features if feat.type != "source"
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






