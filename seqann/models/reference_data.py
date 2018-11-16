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
import sys
import csv
import glob
import pickle
import logging
import pymysql
import urllib.request
from typing import List, Dict

from seqann.util import get_structures
from seqann.util import get_structorder

from Bio import SeqIO
from BioSQL import BioSeqDatabase
from Bio.SeqRecord import SeqRecord

from seqann.util import get_features
from seqann.util import deserialize_model
from seqann.models.base_model_ import Model
from seqann.models.annotation import Annotation

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


# TODO: Use the AWS Lamba API for getting latest IMGT/DB

class ReferenceData(Model):
    '''
    classdocs
    '''
    def __init__(self, server: BioSeqDatabase=None,
                 datafile: str=None, dbversion: str='3310',
                 alleles: List=None, seqdata: Dict=None, hladata: Dict=None,
                 featuredata=None,
                 kir: bool=False, alignments: bool=False,
                 verbose: bool=False, verbosity: int=0):
        """
        ReferenceData - a model defined in Swagger
        :param server: The server of this ReferenceData.
        :type server: BioSeqDatabase
        :param datafile: The datafile of this ReferenceData.
        :type datafile: str
        :param dbversion: The dbversion of this ReferenceData.
        :type dbversion: str
        """
        self.data_types = {
            'server': BioSeqDatabase,
            'datafile': str,
            'dbversion': str,
            'hla_names': List[str],
            'feature_lengths': Dict,
            'hlaref': Dict,
            'seqref': Dict,
            'feature_lengths': Dict,
            'structure_max': Dict,
            'struct_order': Dict,
            'structures': Dict,
            'blastdb': str,
            'server_avail': bool,
            'verbose': bool,
            'verbosity': int,
            'alignments': bool
        }

        self.attribute_map = {
            'seqdata': 'seqdata',
            'hlaref': 'hlaref',
            'seqref': 'seqref',
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
            'alignments': 'alignments',
            'verbose': 'verbose',
            'verbosity': 'verbosity'
        }
        self._seqref = {}
        self._hlaref = {}
        self._kir = kir
        self._verbose = verbose
        self._verbosity = verbosity
        self._dbversion = dbversion
        self._server = server
        self._datafile = datafile
        self._alignments = alignments
        self._server_avail = True if server else False

        self.logger = logging.getLogger("Logger." + __name__)

        hla_url = 'https://raw.githubusercontent.com/ANHIG/IMGTHLA/' \
            + dbversion + '/hla.dat'
        kir_url = 'ftp://ftp.ebi.ac.uk/pub/databases/ipd/kir/KIR.dat'
        hla_loci = ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-DRB1', 'HLA-DQB1',
                    'HLA-DPB1', 'HLA-DQA1', 'HLA-DPA1', 'HLA-DRB3',
                    'HLA-DRB4', 'HLA-DRB5', 'HLA-DRA']

        if self.verbose and verbosity > 0:
            self.logger.info("IPD-IMGT/HLA release = " + str(dbversion))
            self.logger.info("HLA URL = " + hla_url)
            self.logger.info("KIR URL = " + kir_url)
            if self.server_avail:
                self.logger.info("Using BioSQL Server")
                self.logger.info("BIOSQLUSER = " + biosqluser)
                self.logger.info("BIOSQLHOST = " + biosqlhost)
                self.logger.info("BIOSQLDB = " + biosqldb)
                self.logger.info("BIOSQLPORT = " + str(biosqlport))

        # TODO: ** Have script seqann --setup (--latest|--release|--all)
        #           - downloads and creates all files
        #           - removes all data files except alignment files
        #           - Creates blast db
        #
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

        if alleles:
            self._hla_names = alleles
        else:
            splitter = " " if int(dbversion) < 3320 else ","
            # Open allele list file
            try:
                with open(allele_list, 'r') as f:
                    for line in f:
                        line = line.rstrip()
                        if re.search("#", line):
                            continue
                        accession, name = line.split(splitter)
                        if not kir:
                            hla_names.append("HLA-" + name)
                        else:
                            hla_names.append(name)
                    f.close()
                if self.verbose and verbosity > 0:
                    self.logger.info("Loaded " + str(len(hla_names)) + " allele names")
            except OSError as err:
                self.logger.error("OS error: {0}".format(err))
            except:
                self.logger.error("Unexpected error:", sys.exc_info()[0])
                raise
            self._hla_names = hla_names

        #if self.verbose:
        #    mem = "{:4.4f}".format(sys.getsizeof(self.all_feats) / 1000000)
        #    self.logger.info(self.logname + "Finished loading all features * all_feats = " + mem + " MB *")

        feature_lengths = {}
        columns = ['mean', 'std', 'min', 'max']

        featurelength_file = ''
        if kir:
            featurelength_file = data_dir + "/../data/kir-feature_lengths.csv"
        else:
            featurelength_file = data_dir + "/../data/feature_lengths.csv"

        if featuredata:
            self._feature_lengths = featuredata
        else:
            # TODO: use pandas
            try:
                columns = ['mean', 'std', 'min', 'max']
                with open(featurelength_file, newline='') as csvfile:
                    reader = csv.DictReader(csvfile)
                    for row in reader:
                        ldata = [row[c] for c in columns]
                        if row['locus'] in feature_lengths:
                            feature_lengths[row['locus']].update({row['feature']: ldata})
                        else:
                            feature_lengths.update({row['locus']: {row['feature']: ldata}})
                    csvfile.close()
            except OSError as err:
                self.logger.error("OS error: {0}".format(err))
            except:
                self.logger.error("Unexpected error:", sys.exc_info()[0])
                raise

            self._feature_lengths = feature_lengths

        self._blastdb = blastdb
        self._hla_loci = hla_loci

        self._structures = get_structures()
        self._struct_order = get_structorder()

        self._structure_max = {'KIR2DP1': 20, 'KIR2DL5A': 20, 'KIR2DS4': 20,
                               'HLA-DPA1': 9, 'HLA-DQA1': 9, 'KIR2DL2': 20,
                               'HLA-DPB1': 11, 'KIR2DS2': 20, 'KIR3DP1': 20,
                               'HLA-DRB4': 13, 'KIR2DL1': 20, 'KIR2DS5': 20,
                               'HLA-DRB3': 13, 'KIR2DS3': 20, 'KIR3DL1': 20,
                               'HLA-A': 17, 'HLA-DRB5': 13, 'KIR2DL4': 20,
                               'HLA-DQB1': 13, 'KIR3DL2': 20, 'HLA-B': 15,
                               'KIR3DS1': 20, 'KIR2DL5B': 20, 'HLA-DRB1': 13,
                               'KIR3DL3': 20, 'KIR2DS1': 20, 'HLA-C': 17,
                               'HLA-DRA': 9}

        # Starting location of sequence for IPD-IMGT/HLA alignments
        # ** DRA location is not right **
        self.location = {"HLA-A": -300, "HLA-B": -284, "HLA-C": -283,
                         "HLA-DRB1": -599, "HLA-DRB3": -327, "HLA-DRB4": -313,
                         "HLA-DQB1": -525, "HLA-DPB1": -366, "HLA-DPA1": -523,
                         "HLA-DQA1": -746,
                         "HLA-DRA": -400}

        self.align_coordinates = {}
        self.annoated_alignments = {}
        if alignments:
            pickle_dir = data_dir + '/../data/alignments/' + dbversion
            pickle_files = glob.glob(pickle_dir + '/*.pickle')
            for pickle_file in pickle_files:
                locus = pickle_file.split("/")[len(pickle_file.split("/"))-1].split(".")[0].split("_")[0]
                if self.verbose:
                    self.logger.info("Loading " + pickle_file)
                with open(pickle_file, 'rb') as handle:
                    self.annoated_alignments.update({locus:
                                                     pickle.load(handle)})
                    handle.close()
                allele = list(self.annoated_alignments[locus].keys())[0]
                if not locus in self.align_coordinates and "HLA-" + locus in self.struct_order:
                    start = 0
                    feat_order = list(self.struct_order["HLA-" + locus].keys())
                    feat_order.sort()
                    self.align_coordinates.update({locus: {}})
                    if self.verbose and self.verbosity > 2:
                        self.logger.info("* Alignment coordinates *")
                    for i in feat_order:
                        feat = self.struct_order["HLA-" + locus][i]
                        seq = self.annoated_alignments[locus][allele][feat]['Seq']
                        end = start + len(seq)
                        if self.verbose and self.verbosity > 2:
                            self.logger.info(feat + " start = " + str(start) + " | end = " + str(end))
                        for j in range(start, end):
                            self.align_coordinates[locus].update({j: feat})
                        start = end

        # If no server is provided
        # download the dat file
        if seqdata and hladata:
            self._hlaref = hladata
            self._seqref = seqdata
        elif not self._server_avail:
            if kir:
                datfile = data_dir + '/../data/KIR.dat'
            else:
                datfile = data_dir + '/../data/' + dbversion + '.hla.dat'

            if not os.path.isfile(datfile) and not kir:
                if self.verbose:
                    self.logger.info("Downloding KIR data file - " + datfile)
                download_dat(hla_url, datfile)
            elif not os.path.isfile(datfile) and kir:
                if self.verbose:
                    self.logger.info("Downloding HLA data file - " + datfile)
                download_dat(kir_url, datfile)

            # Load HLA dat file
            seqref_pickle = data_dir \
                + '/../data/seqref.' + dbversion + ".pickle"

            hlaref_pickle = data_dir \
                + '/../data/hlaref.' + dbversion + ".pickle"

            if not os.path.isfile(seqref_pickle) or \
                    not os.path.isfile(hlaref_pickle):

                hladata = SeqIO.parse(datfile, "imgt")
                for seqrec in hladata:
                    seqname = seqrec.description.split(",")[0]
                    locus = seqname.split("*")[0]
                    if locus in self.structure_max:
                        self._hlaref.update({seqname: seqrec})
                        self._seqref.update({str(seqrec.seq): seqname})

                if self.verbose:
                    self.logger.info("Finished loading dat file")
                    self.logger.info("Writing pickle of dat file")

                with open(seqref_pickle, 'wb') as handle:
                    pickle.dump(self._seqref, handle,
                                protocol=pickle.HIGHEST_PROTOCOL)
                    handle.close()
                with open(hlaref_pickle, 'wb') as handle:
                    pickle.dump(self._hlaref, handle,
                                protocol=pickle.HIGHEST_PROTOCOL)
                    handle.close()
            else:
                if self.verbose:
                    self.logger.info("Loading pickle dat file")
                with open(seqref_pickle, 'rb') as handle:
                    self._seqref = pickle.load(handle)
                    handle.close()
                with open(hlaref_pickle, 'rb') as handle:
                    self._hlaref = pickle.load(handle)
                    handle.close()

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
    def verbosity(self) -> int:
        """
        Gets the server of this ReferenceData.

        :return: The server of this ReferenceData.
        :rtype: BioSeqDatabase
        """
        return self._verbosity

    @verbosity.setter
    def verbosity(self, verbosity: int):
        """
        Sets the verbosity of this bool.

        :param verbosity: The server of this ReferenceData.
        :type verbosity: bool
        """
        self._verbosity = verbosity

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
    def hlaref(self) -> Dict:
        """
        Gets the hlaref of this ReferenceData.

        :return: The hlaref of this ReferenceData.
        :rtype: Dict
        """
        return self._hlaref

    @hlaref.setter
    def hlaref(self, hlaref: Dict):

        # Sets the seqref of this ReferenceData.

        # :param seqref: The seqref of this ReferenceData.
        # :type seqref: str

        self._hlaref = hlaref

    @property
    def seqref(self) -> Dict:
        """
        Gets the seqref of this ReferenceData.

        :return: The seqref of this ReferenceData.
        :rtype: Dict
        """
        return self._seqref

    @seqref.setter
    def seqref(self, seqref: Dict):

        # Sets the seqref of this ReferenceData.

        # :param seqref: The seqref of this ReferenceData.
        # :type seqref: str

        self._seqref = seqref

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

    def search_refdata(self, seq, locus):
        """
        This checks to see if a sequence already exists in the reference data. If it does, then it'll return the known annotation.

        :return: The Annotation of associated with the input sequence
        :rtype: :ref:`ann`

        Example:

            >>> from Bio.Seq import Seq
            >>> from seqann.models.reference_data import ReferenceData
            >>> sequence = Seq('AGAGACTCTCCCGAGGATTTCGTGTACCAGTTTAAGGCCATGTGCTACTTCACC')
            >>> refdata = ReferenceData()
            >>> matched_annotation = refdata.search_refdata(sequence, locus)

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
                                   user=biosqluser, passwd=biosqlpass,
                                   db=biosqldb)
            cur = conn.cursor()
            cur.execute(select_stm)

            typ = ''
            for row in cur:
                typ = row[0]

            cur.close()
            conn.close()

            if typ:
                if self.verbose:
                    self.logger.info("Exact typing found in BioSQL database")
                seqrecord = self.seqrecord(typ, loc)
                return self.seqannotation(seqrecord, typ, loc)
            else:
                return
        else:
            if str(seq.seq) in self.seqref:
                if self.verbose:
                    self.logger.info("Exact typing found in dat file")
                seqrec = self.hlaref[self.seqref[str(seq.seq)]]
                return self.seqannotation(seqrec,
                                          self.seqref[str(seq.seq)],
                                          locus)
            else:
                return

    # def refseqs(self, locus, n):
    #     hla, loc = locus.split('-')
    #     if self.server_avail:
    #         select_stm = "SELECT ent.name " + \
    #             "FROM bioentry ent,biosequence seq,biodatabase dbb " + \
    #             "WHERE dbb.biodatabase_id = ent.biodatabase_id AND " + \
    #             "seq.bioentry_id = ent.bioentry_id " + \
    #             "AND dbb.name = \"" + self.dbversion + "_" + loc + "\" " + \
    #             "LIMIT " + n

    #         # TODO: add try statement
    #         conn = pymysql.connect(host=biosqlhost, port=biosqlport,
    #                                user=biosqluser, passwd=biosqlpass,
    #                                db=biosqldb)
    #         cur = conn.cursor()
    #         cur.execute(select_stm)

    #         typing = []
    #         for row in cur:
    #             typing.append(self.seqrecord(row[0], loc))

    #         cur.close()
    #         conn.close()

    #         if typing:
    #             if self.verbose:
    #                 self.logger.info("Exact typing found in BioSQL database")
    #             return typing
    #         else:
    #             return
    #     else:
    #         typings = [a for a in self.imgtdat
    #                    if a.description.split(",")[0].split("*")[0] == locus][0:n]
    #         return typings

    def seqrecord(self, allele, locus):
        """
        Gets the Annotation from the found sequence

        :return: The Annotation from the found sequence
        :rtype: Annotation
        """
        try:
            db = self.server[self.dbversion + "_" + locus]
        except:
            self.logger.error("The database " + self.dbversion + "_"
                              + locus + " does not exist!")
            return ''

        seqrecord = db.lookup(name=allele)
        return seqrecord

    def seqannotation(self, seqrecord, allele, loc):
        """
        Gets the Annotation from the found sequence

        :return: The Annotation from the found sequence
        :rtype: Annotation
        """
        #seqrecord = self.seqrecord(allele, loc)
        complete_annotation = get_features(seqrecord)
        annotation = Annotation(annotation=complete_annotation,
                                method='match',
                                complete_annotation=True)

        if self.alignments:
            alignment = {f: self.annoated_alignments[loc][allele][f]['Seq']
                         for f in self.annoated_alignments[loc][allele].keys()}
            annotation.aligned = alignment

        return annotation






