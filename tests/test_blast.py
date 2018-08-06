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
import logging
import unittest

from Bio import SeqIO
from BioSQL import BioSeqDatabase

from seqann.blast_cmd import blastn
from seqann.models.blast import Blast
from seqann.blast_cmd import get_locus
from seqann.models.reference_data import ReferenceData

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

verbose = False
# if os.getenv("VERBOSE"):
#     if os.getenv("VERBOSE") == "True" \
#             or int(os.getenv("VERBOSE")) == 1:
#         logging.basicConfig(format='%(asctime)s - %(name)-35s - %(levelname)-5s - %(funcName)s %(lineno)d: - %(message)s',
#                             datefmt='%m/%d/%Y %I:%M:%S %p',
#                             level=logging.INFO)
#         verbose = True

verbosity = 0
if os.getenv("VERBOSITY"):
    verbosity = int(os.getenv("VERBOSITY"))


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


class TestBlast(unittest.TestCase):

    def setUp(self):
        self.data_dir = os.path.dirname(__file__) + "/resources"
        pass

    @unittest.skipUnless(conn(), "TestBlast 001 MySQL connection")
    def test_001_blastserv(self):
        input_seq = self.data_dir + '/ambig_seqs.fasta'
        in_seq = list(SeqIO.parse(input_seq, "fasta"))[0]
        server = BioSeqDatabase.open_database(driver="pymysql",
                                              user=biosqluser,
                                              passwd=biosqlpass,
                                              host=biosqlhost,
                                              db=biosqldb,
                                              port=biosqlport)
        refdata = ReferenceData(server=server)
        self.assertFalse(refdata.seqref)
        self.assertFalse(refdata.hlaref)
        blast_o = blastn(in_seq, 'HLA-A', 3, refdata=refdata)
        self.assertIsInstance(blast_o, Blast)
        self.assertFalse(blast_o.failed)
        self.assertEqual(blast_o.alleles[0], "HLA-A*01:01:01:01")
        self.assertEqual(len(blast_o.match_seqs), 3)
        server.close()
        pass

    def test_002_blast(self):
        input_seq = self.data_dir + '/ambig_seqs.fasta'
        in_seq = list(SeqIO.parse(input_seq, "fasta"))[0]
        refdata = ReferenceData()
        self.assertTrue(refdata.seqref)
        self.assertTrue(refdata.hlaref)
        blast_o = blastn(in_seq, 'HLA-A', 3, refdata=refdata)
        self.assertIsInstance(blast_o, Blast)
        self.assertFalse(blast_o.failed)
        self.assertEqual(blast_o.alleles[0], "HLA-A*01:01:01:01")
        self.assertEqual(len(blast_o.match_seqs), 3)
        pass

    def test_002_blastnoloc(self):
        input_seq = self.data_dir + '/partial_seqs.fasta'
        in_seq = list(SeqIO.parse(input_seq, "fasta"))[0]
        refdata = ReferenceData()
        self.assertTrue(refdata.seqref)
        self.assertTrue(refdata.hlaref)
        locus = get_locus(in_seq, refdata=refdata, verbose=verbose)
        blast_o = blastn(in_seq, locus, 3, refdata=refdata, verbose=verbose)
        self.assertIsInstance(blast_o, Blast)
        self.assertFalse(blast_o.failed)
        self.assertEqual(blast_o.alleles[0], "HLA-A*01:01:01:12")
        self.assertEqual(len(blast_o.match_seqs), 3)
        pass

    @unittest.skipUnless(conn(), "TestBlast 002 MySQL connection")
    def test_003_fail(self):
        input_seq = self.data_dir + '/failed_seqs.fasta'
        in_seq = list(SeqIO.parse(input_seq, "fasta"))[0]
        server = BioSeqDatabase.open_database(driver="pymysql",
                                              user=biosqluser,
                                              passwd=biosqlpass,
                                              host=biosqlhost,
                                              db=biosqldb,
                                              port=biosqlport)
        refdata = ReferenceData(server=server)
        self.assertFalse(refdata.seqref)
        self.assertFalse(refdata.hlaref)
        blast_o = blastn(in_seq, 'HLA-A', 3, refdata=refdata)
        self.assertIsInstance(blast_o, Blast)
        self.assertTrue(blast_o.failed)
        self.assertFalse(blast_o.alleles)
        server.close()
        pass

    @unittest.skipUnless(conn(), "TestBlast 003 MySQL connection")
    def test_004_nolocserv(self):
        input_seq = self.data_dir + '/ambig_seqs.fasta'
        in_seq = list(SeqIO.parse(input_seq, "fasta"))[0]
        server = BioSeqDatabase.open_database(driver="pymysql",
                                              user=biosqluser,
                                              passwd=biosqlpass,
                                              host=biosqlhost,
                                              db=biosqldb,
                                              port=biosqlport)
        refdata = ReferenceData(server=server)
        self.assertFalse(refdata.seqref)
        self.assertFalse(refdata.hlaref)
        locus = get_locus(in_seq, refdata=refdata)
        self.assertIsInstance(locus, str)
        self.assertTrue(locus)
        self.assertEqual(locus, "HLA-A")
        server.close()
        pass

    def test_005_noloc(self):
        input_seq = self.data_dir + '/ambig_seqs.fasta'
        in_seq = list(SeqIO.parse(input_seq, "fasta"))[0]
        refdata = ReferenceData()
        self.assertTrue(refdata.seqref)
        self.assertTrue(refdata.hlaref)
        locus = get_locus(in_seq, refdata=refdata)
        self.assertIsInstance(locus, str)
        self.assertTrue(locus)
        self.assertEqual(locus, "HLA-A")
        pass

    def test_006_noref(self):
        input_seq = self.data_dir + '/ambig_seqs.fasta'
        in_seq = list(SeqIO.parse(input_seq, "fasta"))[0]
        blast_o = blastn(in_seq, 'HLA-A', 3)
        self.assertIsInstance(blast_o, Blast)
        self.assertFalse(blast_o.failed)
        self.assertEqual(blast_o.alleles[0], "HLA-A*01:01:01:01")
        self.assertEqual(len(blast_o.match_seqs), 3)


