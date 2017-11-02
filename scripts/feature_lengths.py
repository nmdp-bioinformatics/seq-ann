# -*- coding: utf-8 -*-

#
#    gene_feature_enumeration Gene Feature Enumeration.
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

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
from Bio.Alphabet import IUPAC
import os
import glob
from BioSQL import BioSeqDatabase
import re
from datetime import date, datetime
from typing import List, Dict
import pymysql
import pandas as pd

is_kir = lambda x: True if re.search("KIR", x) else False

hladat = 'seqann/data/hla.dat'
hladata = SeqIO.parse(hladat, "imgt")

loci = ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-DRB1', 'HLA-DQB1',
        'HLA-DPB1', 'HLA-DQA1', 'HLA-DPA1', 'HLA-DRB3',
        'HLA-DRB4']

feature_cnts = []
for a in hladata:
    allele_name = a.description.split(",")[0]
    loc, allele = allele_name.split("*")
    if loc not in loci or not a.features:
        continue

    print(allele_name)
    fiveutr = [["five_prime_UTR", a.features[i].extract(a.seq)] for i in range(0, 3) if a.features[i].type != "source"
               and a.features[i].type != "CDS" and isinstance(a.features[i], SeqFeature)
               and not a.features[i].qualifiers]
    feats = [[str(feat.type + "_" + feat.qualifiers['number'][0]), feat.extract(a.seq)]
             for feat in a.features if feat.type != "source"
             and feat.type != "CDS" and isinstance(feat, SeqFeature)
             and 'number' in feat.qualifiers]
    threeutr = [["three_prime_UTR", a.features[i].extract(a.seq)] for i in range(len(a.features)-1, len(a.features)) if a.features[i].type != "source"
                and a.features[i].type != "CDS" and isinstance(a.features[i], SeqFeature)
                and not a.features[i].qualifiers]
    features = fiveutr + feats + threeutr
    for feature in features:
        feat_name = feature[0]
        feature_cnts.append([loc, feat_name, len(feature[1])])

counts_df = pd.DataFrame(feature_cnts, columns=['locus', 'feature', 'length'])
summary = counts_df.groupby(['locus', 'feature'])['length'].describe()
summary.to_csv("feature_lengths.csv")








