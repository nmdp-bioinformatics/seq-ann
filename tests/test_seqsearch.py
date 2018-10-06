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
test_seqsearch
----------------------------------
pip install seqann

from seqann import SeqAnn

annotator = SeqAnn()

annotator.annotate()


Tests for `seqann.seqsearch` module.
"""

import os
import pymysql
import unittest

from BioSQL import BioSeqDatabase
from seqann.models.reference_data import ReferenceData
from seqann.seq_search import SeqSearch

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


class TestSeqSearch(unittest.TestCase):

    def setUp(self):
        self.data_dir = os.path.dirname(__file__) + "/resources"
        pass

    @unittest.skipUnless(conn(), "TestSeqAnn 001 Requires MySQL connection")
    def test_001_wrefdata(self):
        # server = BioSeqDatabase.open_database(driver="pymysql",
        #                                       user=biosqluser,
        #                                       passwd=biosqlpass,
        #                                       host=biosqlhost,
        #                                       db=biosqldb,
        #                                       port=biosqlport)
        #refdata = ReferenceData(server=server)
        seqsearch = SeqSearch()
        self.assertIsInstance(seqsearch, SeqSearch)
        #self.assertTrue(seqsearch.refdata.server_avail)
        #server.close()
        pass

    # @unittest.skipUnless(conn(), "TestSeqAnn 002 Requires MySQL connection")
    # def test_002_search(self):
    #     server = BioSeqDatabase.open_database(driver="pymysql",
    #                                           user=biosqluser,
    #                                           passwd=biosqlpass,
    #                                           host=biosqlhost,
    #                                           db=biosqldb,
    #                                           port=biosqlport)
    #     refdata = ReferenceData(server=server, dbversion='3290')
    #     seqsearch = SeqSearch(refdata=refdata)

    #     input_seq = self.data_dir + '/classi_seq1-ambig.fasta'
    #     in_seq = list(SeqIO.parse(input_seq, "fasta"))[0]
    #     bl = blastn(in_seq, 'HLA-A', 3, refdata=refdata)
    #     annotation = seqsearch.search_seqs(bl.match_seqs[0], in_seq, 'HLA-A')
    #     self.assertIsInstance(annotation, Annotation)
    #     self.assertTrue(annotation.complete_annotation)
    #     server.close()
    #     pass

    # @unittest.skipUnless(conn(), "TestSeqAnn 003 Requires MySQL connection")
    # def test_003_fail(self):
    #     server = BioSeqDatabase.open_database(driver="pymysql",
    #                                           user=biosqluser,
    #                                           passwd=biosqlpass,
    #                                           host=biosqlhost,
    #                                           db=biosqldb,
    #                                           port=biosqlport)
    #     refdata = ReferenceData(server=server, dbversion='3290')
    #     seqsearch = SeqSearch(refdata=refdata)

    #     input_seq = self.data_dir + '/classi_seq5-nomatch.fasta'
    #     in_seq = list(SeqIO.parse(input_seq, "fasta"))[0]
    #     bl = blastn(in_seq, 'HLA-A', 3, refdata=refdata)
    #     annotation = seqsearch.search_seqs(bl.match_seqs[0], in_seq, 'HLA-A')
    #     self.assertIsInstance(annotation, Annotation)
    #     server.close()
    #     pass

    # def test_004_init(self):
    #     seqsearch = SeqSearch()
    #     input_seq = self.data_dir + '/classi_seq1-ambig.fasta'
    #     in_seq = list(SeqIO.parse(input_seq, "fasta"))[0]
    #     bl = blastn(in_seq, 'HLA-A', 3)
    #     annotation = seqsearch.search_seqs(bl.match_seqs[0], in_seq, 'HLA-A')
    #     self.assertIsInstance(annotation, Annotation)
    #     self.assertTrue(annotation.complete_annotation)
    #     pass







