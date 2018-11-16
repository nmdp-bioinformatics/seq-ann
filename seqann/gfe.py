# -*- coding: utf-8 -*-

#
#    pygfe pyGFE.
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

import sys
import logging

from Bio.Seq import Seq
from BioSQL.BioSeq import DBSeq

from seqann.util import get_structures
from seqann.util import get_structorder
from seqann.util import isutr

from seqann.feature_client.rest import ApiException
from seqann.feature_client.api_client import ApiClient
from seqann.feature_client.models.feature import Feature
from seqann.feature_client.apis.features_api import FeaturesApi
from seqann.feature_client.models.feature_request import FeatureRequest


class GFE(object):
    '''
    This class is used for converting annotations into GFE notations.

    Example:

        >>> from Bio import SeqIO
        >>> from BioSQL import BioSeqDatabase
        >>> from seqann.sequence_annotation import BioSeqAnn
        >>> from pygfe.pygfe import pyGFE
        >>> seq_file = 'test_dq.fasta'
        >>> gfe = pyGFE()
        >>> server = BioSeqDatabase.open_database(driver="pymysql", user="root",
        ...                                       passwd="", host="localhost",
        ...                                      db="bioseqdb")
        >>> seqann = BioSeqAnn(server=server)
        >>> seq_rec = list(SeqIO.parse(seq_file, 'fasta'))[0]
        >>> annotation = seqann.annotate(seq_rec, "HLA-DQB1")
        >>> features, gfe = gfe.get_gfe(annotation, "HLA-DQB1")
        >>> print(gfe)
        HLA-DQB1w0-4-0-141-0-12-0-4-0-0-0-0-0

    '''
    def __init__(self, url="http://feature.nmdp-bioinformatics.org",
                 loci=['KIR2DP1', 'KIR2DL5A', 'KIR2DS4', 'HLA-DRA', 'HLA-DPA1', 'HLA-DQA1', 'HLA-DPB1', 'KIR2DS2', 'KIR3DP1', 'HLA-DRB4', 'KIR2DL1', 'KIR2DS5', 'HLA-DRB3', 'KIR2DS3', 'KIR3DL1', 'HLA-A', 'HLA-DRB5', 'KIR2DL4', 'HLA-DQB1', 'KIR3DL2', 'HLA-B', 'KIR3DS1', 'KIR2DL5B', 'HLA-DRB1', 'KIR3DL3', 'KIR2DS1', 'HLA-C'],
                 load_features=False, store_features=False,
                 cached_features=None,
                 verbose=False,
                 pid="NA",
                 verbosity=0):

        self.loci = loci
        self.verbose = verbose
        self.verbosity = verbosity
        self.store_features = store_features
        self.logger = logging.getLogger("Logger." + __name__)
        self.logname = "ID {:<10} - ".format(str(pid))
        client = ApiClient(host=url)
        api_instance = FeaturesApi(api_client=client)
        self.api = api_instance
        self.all_feats = {loc: {} for loc in loci}
        self.structures = get_structures()
        self.struct_order = get_structorder()

        if cached_features:
            if verbose:
                self.logger.info(self.logname + "Using cached features")
            self.all_feats = cached_features

        # Load all features from feature service
        if load_features and not cached_features:
            if verbose:
                self.logger.info(self.logname + "Loading features...")

            # Calling load_features() to load
            # features at each locus
            self.load_features()

    def load_features(self):
        """
        Loads all the known features from the feature service
        """
        # Loading all loci that
        # are in self.loci variable defined
        # when the pyGFE object is created
        for loc in self.loci:
            if self.verbose:
                self.logger.info(self.logname + "Loading features for " + loc)

            # Loading all features for loc from feature service
            self.all_feats.update({loc: self.locus_features(loc)})

            if self.verbose:
                self.logger.info(self.logname + "Finished loading features for " + loc)

        if self.verbose:
            mem = "{:4.4f}".format(sys.getsizeof(self.all_feats) / 1000000)
            self.logger.info(self.logname + "Finished loading all features * all_feats = " + mem + " MB *")

    def locus_features(self, locus):
        """
        Returns all features associated with a locus

        :param locus: string containing HLA locus.
        :type locus: ``str``
        :rtype: ``dict``
        """
        features = self.api.list_features(locus=locus)
        feat_dict = {":".join([a.locus, str(a.rank), a.term, a.sequence]): a.accession for a in features}
        return feat_dict

    def get_gfe(self, annotation, locus):
        """
        creates GFE from a sequence annotation

        :param locus:  The gene locus
        :type locus: ``str``
        :param annotation: An sequence annotation object
        :type annotation: ``List``
        :rtype: ``List``

        Returns:
            The GFE notation and the associated features in an array

        """
        features = []
        accessions = {}
        for feat in annotation.annotation:
            if isinstance(annotation.annotation[feat], DBSeq) \
                    or isinstance(annotation.annotation[feat], Seq):
                seq = str(annotation.annotation[feat])
            else:
                seq = str(annotation.annotation[feat].seq)

            # TODO: Drop this if statement
            if isutr(feat):
                feat_str = ":".join([locus, str(1), feat, seq])

                # If the feature has been loaded or stored
                # then use that instead of making a feature request
                if self.verbose and self.verbosity > 2:
                    self.logger.info("Getting accession " + feat_str)

                if feat_str in self.all_feats[locus]:

                    if self.verbose and self.verbosity > 2:
                        self.logger.info("Feature found " + feat_str)

                    accession = self.all_feats[locus][feat_str]
                    feature = Feature(term=feat,
                                      rank=1,
                                      locus=locus,
                                      sequence=seq,
                                      accession=accession)
                    accessions.update({feat: accession})
                    features.append(feature)
                else:
                    if self.verbose and self.verbosity > 2:
                        self.logger.info(self.logname + "Making FeatureRequest " + feat_str)

                    # Create FeatureRequest object
                    request = FeatureRequest(locus=locus,
                                             term=feat,
                                             rank=1,
                                             sequence=seq)

                    # Attempt to make feature request
                    try:
                        feature = self.api.create_feature(body=request)
                        accessions.update({feat: feature.accession})
                        features.append(feature)
                    except ApiException as e:
                        self.logger.error(self.logname + "Exception when calling DefaultApi->create_feature" + e)
                        blank_feat = Feature(term=feat, rank=1, locus=locus,
                                             sequence=seq)
                        accessions.update({feat: 0})
                        features.append(blank_feat)

                    # Store new features for quick retrieval if flag passed
                    if self.store_features:

                        # Adding new feature to all_feats
                        self.all_feats[locus].update({feat_str: feature.accession})

                        # Calculating memory size of all_feats
                        if self.verbose and self.verbosity > 1:
                            self.logger.info(self.logname + "Storing new feature " + feat_str)
                            mem = "{:4.4f}".format(sys.getsizeof(self.all_feats) / 1000000)
                            self.logger.info(self.logname + "Updated * all_feats " + mem + " MB *")

            else:
                term, rank = feat.split("_")
                feat_str = ":".join([locus, str(rank), term, seq])

                # If the feature has been loaded or stored
                # then use that instead of making a feature request
                if feat_str in self.all_feats[locus]:

                    if self.verbose and self.verbosity > 2:
                        self.logger.info(self.logname + "Feature found " + feat_str)

                    accession = self.all_feats[locus][feat_str]
                    feature = Feature(term=term,
                                      rank=rank,
                                      locus=locus,
                                      sequence=seq,
                                      accession=accession)
                    accessions.update({feat: accession})
                    features.append(feature)
                else:

                    if self.verbose and self.verbosity > 2:
                        self.logger.info(self.logname + "Making FeatureRequest " + feat_str)

                    # Create FeatureRequest object
                    request = FeatureRequest(locus=locus,
                                             term=term,
                                             rank=rank,
                                             sequence=seq)

                    # Attempt to make feature request
                    try:
                        feature = self.api.create_feature(body=request)
                        accessions.update({feat: feature.accession})
                        features.append(feature)
                    except ApiException as e:
                        self.logger.error(self.logname + "Exception when calling DefaultApi->create_feature %e" + e)
                        blank_feat = Feature(term=term, rank=rank, locus=locus,
                                             sequence=seq)
                        accessions.update({feat: 0})
                        features.append(blank_feat)

                    # Store new features for quick retrieval if flag passed
                    if self.store_features:

                        # Adding new feature to all_feats
                        self.all_feats[locus].update({feat_str: feature.accession})

                        # Calculating memory size of all_feats
                        if self.verbose and self.verbosity > 1:
                            self.logger.info(self.logname + "Storing new feature " + feat_str)
                            mem = "{:4.4f}".format(sys.getsizeof(self.all_feats) / 1000000)
                            self.logger.info(self.logname + "Updated * all_feats " + mem + " MB *")

        # Creating GFE
        gfe = self._make_gfe(accessions, locus)

        if self.verbose:
            self.logger.info("GFE = " + gfe)

        return features, gfe

    def _seq(self, locus, term, rank, accession):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        try:
            feature = self.api.get_feature_by_path(locus,
                                                   term,
                                                   rank,
                                                   accession)
            return feature
        except ApiException as e:
            print("Exception when calling DefaultApi->get_feature_by_path: %s\n" % e)
            return ''

    def _breakup_gfe(self, gfe):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        [locus, feature_accessions] = gfe.split("w")
        accessions = feature_accessions.split("-")
        i = 0
        features = {}
        for feature_rank in self.structures[locus]:
            accession = accessions[i]
            features.update({feature_rank: accession})
            i += 1

        return(features)

    def _make_gfe(self, features, locus):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        gfe_list = []
        for feat in sorted(self.structures[locus],
                           key=lambda k: self.structures[locus][k]):
            acc = str(0)
            if feat in features:
                acc = str(features[feat])
            gfe_list.append(acc)

        gfea = '-'.join(gfe_list)
        return locus + "w" + gfea
