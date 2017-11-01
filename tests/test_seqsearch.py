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
test_seqsearch
----------------------------------
pip install seqann

from seqann import SeqAnn

annotator = SeqAnn()

annotator.annotate()


Tests for `seqann.seqsearch` module.
"""


import sys
import unittest
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from Bio import SeqIO
from BioSQL import BioSeqDatabase
from seqann.models.reference_data import ReferenceData
from seqann.blast_cmd import blastn
from seqann.seq_search import SeqSearch
from BioSQL.BioSeq import DBSeqRecord
from seqann.models.blast import Blast
from seqann.util import cleanup
from seqann.models.annotation import Annotation
import os
from seqann.align import align_seqs
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from seqann.util import randomid
import pprint
import pymysql


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


class TestSeqSearch(unittest.TestCase):

    def setUp(self):
        self.data_dir = os.path.dirname(__file__) + "/resources"
        pass

    @unittest.skipUnless(conn(), "TestSeqAnn 001 Requires MySQL connection")
    def test_001_wrefdata(self):
        # server = BioSeqDatabase.open_database(driver="pymysql", user="root",
        #                                       passwd="", host="localhost",
        #                                       db="bioseqdb")
        # refdata = ReferenceData(server=server, dbversion='3290')
        # seqsearch = SeqSearch(refdata=refdata)
        # self.assertIsInstance(seqsearch, SeqSearch)
        # self.assertTrue(seqsearch.refdata.server_avail)
        # server.close()
        pass

    # @unittest.skipUnless(conn(), "TestSeqAnn 002 Requires MySQL connection")
    # def test_002_search(self):
    #     server = BioSeqDatabase.open_database(driver="pymysql", user="root",
    #                                           passwd="", host="localhost",
    #                                           db="bioseqdb")
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
    #     server = BioSeqDatabase.open_database(driver="pymysql", user="root",
    #                                           passwd="", host="localhost",
    #                                           db="bioseqdb")
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







