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
test_blast
----------------------------------

Tests for `seqann.blast_cmd` module.
"""

import os
import pymysql
import unittest

from Bio import SeqIO
from BioSQL import BioSeqDatabase

from seqann.models.reference_data import ReferenceData
from seqann.models.blast import Blast
from seqann.blast_cmd import blastn


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


class TestBlast(unittest.TestCase):

    def setUp(self):
        self.data_dir = os.path.dirname(__file__) + "/resources"
        pass

    @unittest.skipUnless(conn(), "TestBlast 001 MySQL connection")
    def test_001_blast(self):
        input_seq = self.data_dir + '/ambig_seqs.fasta'
        in_seq = list(SeqIO.parse(input_seq, "fasta"))[0]
        server = BioSeqDatabase.open_database(driver="pymysql", user="root",
                                              passwd="", host="localhost",
                                              db="bioseqdb")
        refdata = ReferenceData(server=server)
        self.assertFalse(refdata.imgtdat)
        blast_o = blastn(in_seq, 'HLA-A', 3, refdata=refdata)
        self.assertIsInstance(blast_o, Blast)
        self.assertFalse(blast_o.failed)
        self.assertEqual(blast_o.alleles[0], "HLA-A*01:01:01:01")
        self.assertEqual(len(blast_o.alleles), 3)
        server.close()
        pass

    @unittest.skipUnless(conn(), "TestBlast 002 MySQL connection")
    def test_002_fail(self):
        input_seq = self.data_dir + '/failed_seqs.fasta'
        in_seq = list(SeqIO.parse(input_seq, "fasta"))[0]
        server = BioSeqDatabase.open_database(driver="pymysql", user="root",
                                              passwd="", host="localhost",
                                              db="bioseqdb")
        refdata = ReferenceData(server=server)
        self.assertFalse(refdata.imgtdat)
        blast_o = blastn(in_seq, 'HLA-A', 3, refdata=refdata)
        self.assertIsInstance(blast_o, Blast)
        self.assertTrue(blast_o.failed)
        self.assertFalse(blast_o.alleles)
        server.close()
        pass

    def test_003_noref(self):
        input_seq = self.data_dir + '/ambig_seqs.fasta'
        in_seq = list(SeqIO.parse(input_seq, "fasta"))[0]
        blast_o = blastn(in_seq, 'HLA-A', 3)
        self.assertIsInstance(blast_o, Blast)
        self.assertFalse(blast_o.failed)
        self.assertEqual(blast_o.alleles[0], "HLA-A*01:01:01:01")
        self.assertEqual(len(blast_o.alleles), 3)






