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

Tests for `seqann.align` module.
"""

import os
import json
import pymysql
import unittest

from Bio import SeqIO
from BioSQL import BioSeqDatabase

from seqann.sequence_annotation import BioSeqAnn

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


class TestAlign(unittest.TestCase):

    def setUp(self):
        self.data_dir = os.path.dirname(__file__) + "/resources"
        expected_json = self.data_dir + "/expected.json"
        with open(expected_json) as json_data:
            self.expected = json.load(json_data)
        pass

    @unittest.skipUnless(conn(), "TestAlign 001 Requires MySQL connection")
    def test_001_align(self):
        server = BioSeqDatabase.open_database(driver="pymysql",
                                              user=biosqluser,
                                              passwd=biosqlpass,
                                              host=biosqlhost,
                                              db=biosqldb,
                                              port=biosqlport)

        # TODO *** NOT WORKING WITH NO LOCUS           ***
        # TODO *** FIX 3290 Alignments                 ***
        # TODO *** GET ALIGNMENTS WORKING WITH DB SEQS ***
        seqann = BioSeqAnn(server=server, align=True)
        input_seq = self.data_dir + '/align_tests.fasta'
        for ex in self.expected['align']:
            i = int(ex['index'])
            ex = self.expected['align'][i]
            locus = ex['locus']
            allele = ex['name']
            hla, loc = locus.split("-")
            in_seq = list(SeqIO.parse(input_seq, "fasta"))[i]
            annotation = seqann.annotate(in_seq, "HLA-A")
            align = "".join([annotation.aligned[s] for s in annotation.aligned])
            for i in range(0, len(align)-1):
                if str(i) in ex['diffs']:
                    self.assertEqual(list(align)[i],
                                     ex['diffs'][str(i)])
                else:
                    if list(align)[i] != list(ex['alignment'])[i]:
                        print("FAILED:", allele, i, list(align)[i], list(ex['alignment'])[i])
                    self.assertEqual(list(align)[i],
                                     list(ex['alignment'])[i])
        server.close()







