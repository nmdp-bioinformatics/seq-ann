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
from __future__ import absolute_import

import re
import os
import csv
import glob
import pymysql
import collections
import urllib.request
from typing import List, Dict
from datetime import date, datetime

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from BioSQL import BioSeqDatabase
from Bio.SeqRecord import SeqRecord

from seqann.util import deserialize_model
from seqann.models.base_model_ import Model
from seqann.models.annotation import Annotation
from seqann.util import get_features


is_kir = lambda x: True if re.search("KIR", x) else False


def download_dat(dat):
    url = 'https://raw.githubusercontent.com/ANHIG/IMGTHLA/3290/hla.dat'
    urllib.request.urlretrieve(url, dat)


class ReferenceData(Model):
    '''
    classdocs
    '''
    def __init__(self, server: BioSeqDatabase=None, datafile: str=None,
                 dbversion: str='3290', verbose: bool=False):
        """
        ReferenceData - a model defined in Swagger
        :param server: The server of this ReferenceData.
        :type server: BioSeqDatabase
        :param datafile: The datafile of this ReferenceData.
        :type datafile: str
        :param dbversion: The dbversion of this ReferenceData.
        :type dbversion: str
        """
        tree = lambda: collections.defaultdict(tree)
        self.data_types = {
            'server': BioSeqDatabase,
            'datafile': str,
            'dbversion': str,
            'hla_names': List[str],
            'feature_lengths': Dict,
            'structure_max': Dict,
            'struct_order': Dict,
            'structures': Dict,
            'blastdb': str,
            'server_avail': bool,
            'verbose': bool,
            'imgtdat': List[SeqRecord]
        }

        self.attribute_map = {
            'server': 'server',
            'datafile': 'datafile',
            'dbversion': 'dbversion',
            'hla_names': 'hla_names',
            'structure_max': 'structure_max',
            'feature_lengths': 'feature_lengths',
            'struct_order': 'struct_order',
            'structures': 'structures',
            'blastdb': 'blastdb',
            'hla_loci': 'hla_loci',
            'server_avail': 'server_avail',
            'imgtdat': 'imgtdat',
            'verbose': 'verbose'
        }
        self._verbose = verbose
        self._dbversion = dbversion
        self._server = server
        self._datafile = datafile
        self._server_avail = True if server else False

        hla_loci = ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-DRB1', 'HLA-DQB1',
                    'HLA-DPB1', 'HLA-DQA1', 'HLA-DPA1', 'HLA-DRB3',
                    'HLA-DRB4']
        hla_names = []
        data_dir = os.path.dirname(__file__)
        blastdb = data_dir + '/../data/blast/' + dbversion
        allele_list = data_dir + '/../data/allele_lists/Allelelist.' \
                               + dbversion + '.txt'

        # TODO: add try
        with open(allele_list, 'r') as f:
            for line in f:
                line = line.rstrip()
                accession, name = line.split(" ")
                hla_names.append("HLA-" + name)
            f.close()

        feature_lengths = tree()
        columns = ['mean', 'std', 'min', 'max']
        featurelength_file = data_dir + "/../data/feature_lengths.csv"

        # TODO: add try
        with open(featurelength_file, newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                #feat_loc = str(row['locus']) + "_" + str(row['feature'])
                ldata = [row[c] for c in columns]
                feature_lengths[row['locus']][row['feature']] = ldata

        self._feature_lengths = feature_lengths
        self._hla_names = hla_names
        struture_files = glob.glob(data_dir + '/../data/*.structure')
        structures = {}
        struct_order = {}
        structure_max = {}
        for inputfile in struture_files:
            file_path = inputfile.split("/")
            locus = file_path[len(file_path)-1].split(".")[0]
            # TODO: add try
            with open(inputfile, 'r') as f:
                features_order = {}
                features = {}
                n = 0
                for line in f:
                    line = line.rstrip()
                    [feature, rank] = line.split("\t")
                    feature_name = "_".join([feature, rank])
                    if feature == "three_prime_UTR" or feature == "five_prime_UTR":
                        feature_name = feature
                    n += 1
                    features.update({feature_name: n})
                    features_order.update({n: feature_name})
                    if is_kir(locus):
                        structures.update({locus: features})
                        struct_order.update({locus: features_order})
                    else:
                        structures.update({"HLA-" + locus: features})
                        struct_order.update({"HLA-" + locus: features_order})
                structure_max.update({"HLA-" + locus: n})
            f.close()
        self._structures = structures
        self._struct_order = struct_order
        self._structure_max = structure_max
        self._blastdb = blastdb
        self._hla_loci = hla_loci

        if not self._server_avail:
            hladat = data_dir + '/../data/hla.dat'
            if not os.path.isfile(hladat):
                download_dat(hladat)

            # TODO: add try
            hladata = list(SeqIO.parse(hladat, "imgt"))
            print("serverless load")
            self._imgtdat = hladata
        else:
            self._imgtdat = []

    @classmethod
    def from_dict(cls, dikt) -> 'ReferenceData':
        """
        Returns the dict as a model

        :param dikt: A dict.
        :type: dict
        :return: The ReferenceData of this ReferenceData.
        :rtype: ReferenceData
        """
        return deserialize_model(dikt, cls)

    @property
    def server(self) -> BioSeqDatabase:
        """
        Gets the server of this ReferenceData.

        :return: The server of this ReferenceData.
        :rtype: BioSeqDatabase
        """
        return self._server

    @server.setter
    def server(self, server: BioSeqDatabase):
        """
        Sets the server of this ReferenceData.

        :param server: The server of this ReferenceData.
        :type server: BioSeqDatabase
        """
        self._server = server

    @property
    def verbose(self) -> bool:
        """
        Gets the server of this ReferenceData.

        :return: The server of this ReferenceData.
        :rtype: BioSeqDatabase
        """
        return self._verbose

    @verbose.setter
    def verbose(self, verbose: bool):
        """
        Sets the verbose of this bool.

        :param verbose: The server of this ReferenceData.
        :type verbose: bool
        """
        self._verbose = verbose

    @property
    def datafile(self) -> str:
        """
        Gets the datafile of this ReferenceData.

        :return: The datafile of this ReferenceData.
        :rtype: str
        """
        return self._datafile

    @datafile.setter
    def datafile(self, datafile: str):
        """
        Sets the datafile of this ReferenceData.

        :param datafile: The datafile of this ReferenceData.
        :type datafile: str
        """
        self._datafile = datafile

    @property
    def dbversion(self) -> str:
        """
        Gets the dbversion of this ReferenceData.

        :return: The dbversion of this ReferenceData.
        :rtype: str
        """
        return self._dbversion

    @dbversion.setter
    def dbversion(self, dbversion: str):
        """
        Sets the dbversion of this ReferenceData.

        :param dbversion: The dbversion of this ReferenceData.
        :type dbversion: str
        """
        self._dbversion = dbversion

    @property
    def structures(self) -> Dict:
        """
        Gets the structures of this ReferenceData.

        :return: The structures of this ReferenceData.
        :rtype: Dict
        """
        return self._structures

    @property
    def structure_max(self) -> Dict:
        """
        Gets the structure_max of this ReferenceData.

        :return: The structure_max of this ReferenceData.
        :rtype: Dict
        """
        return self._structure_max

    @property
    def blastdb(self) -> str:
        """
        Gets the blastdb of this ReferenceData.

        :return: The blastdb of this ReferenceData.
        :rtype: str
        """
        return self._blastdb

    @property
    def struct_order(self) -> Dict:
        """
        Gets the struct_order of this ReferenceData.

        :return: The struct_order of this ReferenceData.
        :rtype: Dict
        """
        return self._struct_order

    @property
    def feature_lengths(self) -> Dict:
        """
        Gets the feature_lengths of this ReferenceData.

        :return: The feature_lengths of this ReferenceData.
        :rtype: Dict
        """
        return self._feature_lengths

    @property
    def hla_names(self) -> List[str]:
        """
        Gets the hla_names of this ReferenceData.

        :return: The hla_names of this ReferenceData.
        :rtype: Dict
        """
        return self._hla_names

    @property
    def hla_loci(self) -> List[str]:
        """
        Gets the hla_loci of this ReferenceData.

        :return: The hla_loci of this ReferenceData.
        :rtype: List[str]
        """
        return self._hla_loci

    @property
    def server_avail(self) -> bool:
        """
        Gets the server_avail of this ReferenceData.

        :return: The server_avail of this ReferenceData.
        :rtype: bool
        """
        return self._server_avail

    @property
    def imgtdat(self) -> List[SeqRecord]:
        """
        Gets the imgtdat of this ReferenceData.

        :return: The imgtdat of this ReferenceData.
        :rtype: List[SeqRecord]
        """
        return self._imgtdat

    def search_refdata(self, seq, locus):
        """
        Gets the Annotation of this ReferenceData. "select ent.name from bioentry limit 10"
        :return: The Annotation of this ReferenceData.
        :rtype: Annotation
        mysql --user=root bioseqdb -e "select * from biodatabase"
        """
        if self.server_avail:
            hla, loc = locus.split('-')
            p1 = "SELECT ent.name "
            p2 = "FROM bioentry ent,biosequence seq,biodatabase dbb "
            p3 = "WHERE dbb.biodatabase_id = ent.biodatabase_id AND seq.bioentry_id = ent.bioentry_id "
            p4 = "AND dbb.name = \"" + self.dbversion + "_" + loc + "\""
            p5 = "AND seq.seq = \"" + str(seq.seq) + "\""
            select_stm = p1 + p2 + p3 + p4 + p5

            # TODO: add try statement
            conn = pymysql.connect(host='localhost', port=3306, user='root', passwd='', db='bioseqdb')
            cur = conn.cursor()
            cur.execute(select_stm)

            typ = ''
            for row in cur:
                typ = row[0]

            cur.close()
            conn.close()

            if typ:
                return self.seqannotation(seq, typ, loc)
            else:
                return
        else:
            return

    def refseqs(self, locus, n):
        hla, loc = locus.split('-')
        if self.server_avail:
            select_stm = "SELECT ent.name " + \
                "FROM bioentry ent,biosequence seq,biodatabase dbb " + \
                "WHERE dbb.biodatabase_id = ent.biodatabase_id AND " + \
                "seq.bioentry_id = ent.bioentry_id " + \
                "AND dbb.name = \"" + self.dbversion + "_" + loc + "\" " + \
                "LIMIT " + n

            # TODO: add try statement
            conn = pymysql.connect(host='localhost', port=3306,
                                   user='root', passwd='', db='bioseqdb')
            cur = conn.cursor()
            cur.execute(select_stm)

            typing = []
            for row in cur:
                typing.append(self.seqrecord(row[0], loc))

            cur.close()
            conn.close()

            if typing:
                return typing
            else:
                return
        else:
            typings = [a for a in self.imgtdat
                       if a.description.split(",")[0].split("*")[0] == locus][0:n]
            return typings

    def seqrecord(self, allele, locus):
        """
        Gets the Annotation from the found sequence

        :return: The Annotation from the found sequence
        :rtype: Annotation
        """
        # TODO: add try statement
        db = self.server[self.dbversion + "_" + locus]
        seqrecord = db.lookup(name=allele)
        return seqrecord

    def seqannotation(self, seq, allele, locus):
        """
        Gets the Annotation from the found sequence

        :return: The Annotation from the found sequence
        :rtype: Annotation
        """
        seqrecord = self.seqrecord(allele, locus)
        complete_annotation = get_features(seqrecord)
        annotation = Annotation(annotation=complete_annotation,
                                method='match',
                                complete_annotation=True)
        return annotation






