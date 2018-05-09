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

from seqann.util import get_structures
from seqann.util import get_structorder

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from BioSQL import BioSeqDatabase
from Bio.SeqRecord import SeqRecord

from seqann.util import deserialize_model
from seqann.models.base_model_ import Model
from seqann.models.annotation import Annotation
from seqann.util import get_features
import sys
import pickle

is_kir = lambda x: True if re.search("KIR", x) else False


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


def download_dat(url, dat):
    urllib.request.urlretrieve(url, dat)


class ReferenceData(Model):
    '''
    classdocs
    '''
    def __init__(self, server: BioSeqDatabase=None, datafile: str=None,
                 dbversion: str='3310', verbose: bool=False,
                 kir: bool=False, alignments: bool=False):
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
            'alignments': bool,
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
            'kir': 'kir',
            'imgtdat': 'imgtdat',
            'alignments': 'alignments',
            'verbose': 'verbose'
        }
        self._alignments = alignments
        self._kir = kir
        self._verbose = verbose
        self._dbversion = dbversion
        self._server = server
        self._datafile = datafile
        self._server_avail = True if server else False

        hla_url = 'https://raw.githubusercontent.com/ANHIG/IMGTHLA/' + dbversion + '/hla.dat'
        kir_url = 'ftp://ftp.ebi.ac.uk/pub/databases/ipd/kir/KIR.dat'
        hla_loci = ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-DRB1', 'HLA-DQB1',
                    'HLA-DPB1', 'HLA-DQA1', 'HLA-DPA1', 'HLA-DRB3',
                    'HLA-DRB4', 'HLA-DRB5']

        # TODO: Download! Don't have in package!
        hla_names = []
        data_dir = os.path.dirname(__file__)
        if kir:
            blastdb = data_dir + '/../data/blast/KIR'
            allele_list = data_dir + '/../data/allele_lists/Allelelist.' \
                                   + 'KIR.txt'
        else:
            blastdb = data_dir + '/../data/blast/' + dbversion
            allele_list = data_dir + '/../data/allele_lists/Allelelist.' \
                                   + dbversion + '.txt'

        # TODO: add try
        with open(allele_list, 'r') as f:
            for line in f:
                line = line.rstrip()
                accession, name = line.split(" ")
                if not kir:
                    hla_names.append("HLA-" + name)
                else:
                    hla_names.append(name)
            f.close()

        feature_lengths = tree()
        columns = ['mean', 'std', 'min', 'max']

        featurelength_file = ''
        if kir:
            featurelength_file = data_dir + "/../data/kir-feature_lengths.csv"
        else:
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

        self._blastdb = blastdb
        self._hla_loci = hla_loci
        self._structures = get_structures()
        self._struct_order = get_structorder()
        self._structure_max = {'KIR2DP1': 20, 'KIR2DL5A': 20, 'KIR2DS4': 20,
                               'HLA-DPB1': 11, 'KIR2DS2': 20, 'KIR3DP1': 20,
                               'HLA-DRB4': 13, 'KIR2DL1': 20, 'KIR2DS5': 20,
                               'HLA-DRB3': 13, 'KIR2DS3': 20, 'KIR3DL1': 20,
                               'HLA-A': 17, 'HLA-DRB5': 13, 'KIR2DL4': 20,
                               'HLA-DQB1': 13, 'KIR3DL2': 20, 'HLA-B': 15,
                               'KIR3DS1': 20, 'KIR2DL5B': 20, 'HLA-DRB1': 13,
                               'KIR3DL3': 20, 'KIR2DS1': 20, 'HLA-C': 17}

        self.location = {"HLA-A": -300, "HLA-B": -284, "HLA-C": -283,
                         "HLA-DRB1": -599, "HLA-DRB3": -327, "HLA-DRB4": -313,
                         "HLA-DQB1": -525, "HLA-DPB1": -366, "HLA-DPA1": -523,
                         "HLA-DQA1": -746}

        self.align_coordinates = {}
        self.annoated_alignments = {}
        if alignments:
            # TODO: Use logging
            pickle_dir = data_dir + '/../data/alignments/' + dbversion
            pickle_files = glob.glob(pickle_dir + '/*.pickle')
            for pickle_file in pickle_files:
                locus = pickle_file.split("/")[len(pickle_file.split("/"))-1].split(".")[0].split("_")[0]
                with open(pickle_file, 'rb') as handle:
                    self.annoated_alignments.update({locus:
                                                     pickle.load(handle)})
                #print("Finished loading " + locus)
                allele = list(self.annoated_alignments[locus].keys())[0]
                if not locus in self.align_coordinates and "HLA-" + locus in self.struct_order:
                    start = 0
                    feat_order = list(self.struct_order["HLA-" + locus].keys())
                    feat_order.sort()
                    self.align_coordinates.update({locus: {}})
                    for i in feat_order:
                        feat = self.struct_order["HLA-" + locus][i]
                        seq = self.annoated_alignments[locus][allele][feat]['Seq']
                        end = start + len(seq)
                        # WHERE IS THIS USED??
                        for j in range(start, end):
                            self.align_coordinates[locus].update({j: feat})
                        start = end

        # TODO: ADD DB VERSION!
        # TODO: Be able to load KIR and HLA
        if not self._server_avail:
            datfile = ''

            if kir:
                datfile = data_dir + '/../data/KIR.dat'
            else:
                datfile = data_dir + '/../data/' + dbversion + '.hla.dat'

            if not os.path.isfile(datfile) and not kir:
                download_dat(hla_url, datfile)
            elif not os.path.isfile(datfile) and kir:
                download_dat(kir_url, datfile)

            # TODO: add try
            hladata = list(SeqIO.parse(datfile, "imgt"))
            #print("serverless load")
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
    def alignments(self) -> bool:
        """
        Gets the alignments of this ReferenceData.

        :return: The alignments of this ReferenceData.
        :rtype: BioSeqDatabase
        """
        return self._alignments

    @alignments.setter
    def alignments(self, alignments: bool):
        """
        Sets the alignments of this bool.

        :param alignments: The server of this ReferenceData.
        :type alignments: bool
        """
        self._alignments = alignments

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
    def kir(self) -> bool:
        """
        Gets the kir of this ReferenceData.

        :return: The kir of this ReferenceData.
        :rtype: bool
        """
        return self._kir

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
        # TODO: ONLY MAKE ONE CONNECTION
        # TODO: add try statement
        # TODO: take password from environment variable
        if self.server_avail:
            hla, loc = locus.split('-')
            p1 = "SELECT ent.name "
            p2 = "FROM bioentry ent,biosequence seq,biodatabase dbb "
            p3 = "WHERE dbb.biodatabase_id = ent.biodatabase_id AND seq.bioentry_id = ent.bioentry_id "
            p4 = " AND dbb.name = \"" + self.dbversion + "_" + loc + "\""
            p5 = " AND seq.seq = \"" + str(seq.seq) + "\""
            select_stm = p1 + p2 + p3 + p4 + p5

            # TODO: add try statement
            conn = pymysql.connect(host=biosqlhost, port=biosqlport,
                                   user=biosqluser, passwd=biosqlpass, db=biosqldb)
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

            # TODO: ONLY MAKE ONE CONNECTION
            # TODO: add try statement
            # TODO: take password from environment variable
            conn = pymysql.connect(host=biosqlhost, port=biosqlport,
                                   user=biosqluser, passwd=biosqlpass,
                                   db=biosqldb)
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
        # TODO: Add try statement
        db = self.server[self.dbversion + "_" + locus]
        seqrecord = db.lookup(name=allele)
        return seqrecord

    def seqannotation(self, seq, allele, loc):
        """
        Gets the Annotation from the found sequence

        :return: The Annotation from the found sequence
        :rtype: Annotation
        """

        seqrecord = self.seqrecord(allele, loc)
        complete_annotation = get_features(seqrecord)
        annotation = Annotation(annotation=complete_annotation,
                                method='match',
                                complete_annotation=True)

        if self.alignments:
            alignment = {f: self.annoated_alignments[loc][allele][f]['Seq']
                         for f in self.annoated_alignments[loc][allele].keys()}
            annotation.aligned = alignment

        return annotation






