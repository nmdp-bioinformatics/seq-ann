# -*- coding: utf-8 -*-

#
#    seqann Sequence Annotation
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
import re
import os
import string
import random as r
from typing import GenericMeta
from datetime import datetime, date
from six import integer_types, iteritems
from Bio.SeqFeature import SeqFeature

# TODO: Add documentation


def checkseq(sequence: str=None,  code="ATGC") -> bool:
    """
    :param sequence: The input sequence.
    :type sequence: Seq
    :rtype: bool
    """
    for base in sequence:
        if base not in code:
            return False
    return True


def is_kir(feature: str=None) -> bool:
    """
    :param sequence: The input sequence record.
    :type sequence: Seq
    :param locus: The gene locus associated with the sequence.
    :type locus: str
    :param nseqs: The number of blast sequences to use.
    :type nseqs: int
    :rtype: bool
    """
    return True if re.search("KIR", feature) else False


def isexon(feature: str=None) -> bool:
    """
    :param sequence: The input sequence record.
    :type sequence: Seq
    :param locus: The gene locus associated with the sequence.
    :type locus: str
    :param nseqs: The number of blast sequences to use.
    :type nseqs: int
    :rtype: bool
    """
    return True if re.search("exon", feature) else False


def isutr(feature: str=None) -> bool:
    """
    :param sequence: The input sequence record.
    :type sequence: Seq
    :param locus: The gene locus associated with the sequence.
    :type locus: str
    :param nseqs: The number of blast sequences to use.
    :type nseqs: int
    :rtype: bool
    """
    return True if re.search("UTR", feature) else False


def isfive(feature: str=None) -> bool:
    """
    :param sequence: The input sequence record.
    :type sequence: Seq
    :param locus: The gene locus associated with the sequence.
    :type locus: str
    :param nseqs: The number of blast sequences to use.
    :type nseqs: int
    :rtype: bool
    """
    return True if re.search("five", feature) else False


def is_classII(feature: str=None) -> bool:
    """
    :param sequence: The input sequence record.
    :type sequence: Seq
    :param locus: The gene locus associated with the sequence.
    :type locus: str
    :param nseqs: The number of blast sequences to use.
    :type nseqs: int
    :rtype: bool
    """
    return True if re.search("HLA-D", feature) else False


def get_seqfeat(seqrecord):

    n = 3 if len(seqrecord.features) >= 3 else len(seqrecord.features)
    fiveutr = [["five_prime_UTR", seqrecord.features[i]] for i in range(0, n) if seqrecord.features[i].type != "source"
               and seqrecord.features[i].type != "CDS" and isinstance(seqrecord.features[i], SeqFeature)
               and seqrecord.features[i].type != "3UTR"
               and not seqrecord.features[i].qualifiers]
    feats = [[str(feat.type + "_" + feat.qualifiers['number'][0]), feat]
             for feat in seqrecord.features if feat.type != "source"
             and feat.type != "CDS" and isinstance(feat, SeqFeature)
             and 'number' in feat.qualifiers]
    threeutr = [["three_prime_UTR", seqrecord.features[i]] for i in range(len(seqrecord.features)-1, len(seqrecord.features)) if seqrecord.features[i].type != "source"
                and seqrecord.features[i].type != "CDS" and isinstance(seqrecord.features[i], SeqFeature)
                and seqrecord.features[i].type != "5UTR"
                and not seqrecord.features[i].qualifiers]

    feat_list = fiveutr + feats + threeutr
    annotation = {k[0]: k[1] for k in feat_list}
    return(annotation)


def get_seqs(seqrecord):
    n = 3 if len(seqrecord.features) >= 3 else len(seqrecord.features)
    fiveutr = [["five_prime_UTR", str(seqrecord.features[i].extract(seqrecord.seq))] for i in range(0, n) if seqrecord.features[i].type != "source"
               and seqrecord.features[i].type != "CDS"
               and seqrecord.features[i].type != "3UTR"
               and isinstance(seqrecord.features[i], SeqFeature)
               and not seqrecord.features[i].qualifiers]
    feats = [[str(feat.type + "_" + feat.qualifiers['number'][0]), str(feat.extract(seqrecord.seq))]
             for feat in seqrecord.features if feat.type != "source"
             and feat.type != "CDS" and isinstance(feat, SeqFeature)
             and 'number' in feat.qualifiers]
    threeutr = [["three_prime_UTR", str(seqrecord.features[i].extract(seqrecord.seq))] for i in range(len(seqrecord.features)-1, len(seqrecord.features)) if seqrecord.features[i].type != "source"
                and seqrecord.features[i].type != "CDS" and isinstance(seqrecord.features[i], SeqFeature)
                and seqrecord.features[i].type != "5UTR"
                and not seqrecord.features[i].qualifiers]
    feat_list = fiveutr + feats + threeutr
    annotation = {k[0]: str(k[1]) for k in feat_list}
    return(annotation)


# TODO: change name to get_featseq
def get_features(seqrecord):
    n = 3 if len(seqrecord.features) >= 3 else len(seqrecord.features)
    fiveutr = [["five_prime_UTR", seqrecord.features[i].extract(seqrecord.seq)] for i in range(0, n) if seqrecord.features[i].type != "source"
               and seqrecord.features[i].type != "CDS" and isinstance(seqrecord.features[i], SeqFeature)
               and not seqrecord.features[i].qualifiers]
    feats = [[str(feat.type + "_" + feat.qualifiers['number'][0]), feat.extract(seqrecord.seq)]
             for feat in seqrecord.features if feat.type != "source"
             and feat.type != "CDS" and isinstance(feat, SeqFeature)
             and 'number' in feat.qualifiers]
    threeutr = [["three_prime_UTR", seqrecord.features[i].extract(seqrecord.seq)] for i in range(len(seqrecord.features)-1, len(seqrecord.features)) if seqrecord.features[i].type != "source"
                and seqrecord.features[i].type != "CDS" and isinstance(seqrecord.features[i], SeqFeature)
                and not seqrecord.features[i].qualifiers]

    feat_list = fiveutr + feats + threeutr
    annotation = {k[0]: k[1] for k in feat_list}
    return(annotation)


def randomid(N=12):
    random_id = ''.join(r.choices(string.ascii_uppercase + string.digits,
                                  k=N))
    fastafile = str(random_id) + ".fasta"
    xmlfile = str(random_id) + ".xml"
    if os.path.isfile(fastafile) or os.path.isfile(xmlfile):
        N += 1
        return randomid(N=N)
    else:
        return random_id


def cleanup(randid):
    fastafile = str(randid) + ".fasta"
    xmlfile = str(randid) + ".xml"
    clufile = str(randid) + ".clu"

    if os.path.isfile(fastafile):
        os.remove(fastafile)

    if os.path.isfile(xmlfile):
        os.remove(xmlfile)

    if os.path.isfile(clufile):
        os.remove(clufile)


def _deserialize(data, klass):
    """
    Deserializes dict, list, str into an object.

    :param data: dict, list or str.
    :param klass: class literal, or string of class name.

    :return: object.
    """
    if data is None:
        return None

    if klass in integer_types or klass in (float, str, bool):
        return _deserialize_primitive(data, klass)
    elif klass == object:
        return _deserialize_object(data)
    elif klass == date:
        return deserialize_date(data)
    elif klass == datetime:
        return deserialize_datetime(data)
    elif type(klass) == GenericMeta:
        if klass.__extra__ == list:
            return _deserialize_list(data, klass.__args__[0])
        if klass.__extra__ == dict:
            return _deserialize_dict(data, klass.__args__[1])
    else:
        return deserialize_model(data, klass)


def _deserialize_primitive(data, klass):
    """
    Deserializes to primitive type.

    :param data: data to deserialize.
    :param klass: class literal.

    :return: int, long, float, str, bool.
    :rtype: int | long | float | str | bool
    """
    try:
        value = klass(data)
    except UnicodeEncodeError:
        value = unicode(data)
    except TypeError:
        value = data
    return value


def _deserialize_object(value):
    """
    Return a original value.

    :return: object.
    """
    return value


def deserialize_date(string):
    """
    Deserializes string to date.

    :param string: str.
    :type string: str
    :return: date.
    :rtype: date
    """
    try:
        from dateutil.parser import parse
        return parse(string).date()
    except ImportError:
        return string


def deserialize_datetime(string):
    """
    Deserializes string to datetime.

    The string should be in iso8601 datetime format.

    :param string: str.
    :type string: str
    :return: datetime.
    :rtype: datetime
    """
    try:
        from dateutil.parser import parse
        return parse(string)
    except ImportError:
        return string


def deserialize_model(data, klass):
    """
    Deserializes list or dict to model.

    :param data: dict, list.
    :type data: dict | list
    :param klass: class literal.
    :return: model object.
    """
    instance = klass()

    if not instance.swagger_types:
        return data

    for attr, attr_type in iteritems(instance.swagger_types):
        if data is not None \
                and instance.attribute_map[attr] in data \
                and isinstance(data, (list, dict)):
            value = data[instance.attribute_map[attr]]
            setattr(instance, attr, _deserialize(value, attr_type))

    return instance


def _deserialize_list(data, boxed_type):
    """
    Deserializes a list and its elements.

    :param data: list to deserialize.
    :type data: list
    :param boxed_type: class literal.

    :return: deserialized list.
    :rtype: list
    """
    return [_deserialize(sub_data, boxed_type)
            for sub_data in data]


def _deserialize_dict(data, boxed_type):
    """
    Deserializes a dict and its elements.

    :param data: dict to deserialize.
    :type data: dict
    :param boxed_type: class literal.

    :return: deserialized dict.
    :rtype: dict
    """
    return {k: _deserialize(v, boxed_type)
            for k, v in iteritems(data)}


def get_structures():
    return {
        'KIR2DL2':
        {
            'five_prime_UTR': 1,
            'exon_1': 2,
            'intron_1': 3,
            'exon_2': 4,
            'intron_2': 5,
            'exon_3': 6,
            'intron_3': 7,
            'exon_4': 8,
            'intron_4': 9,
            'exon_5': 10,
            'intron_5': 11,
            'exon_6': 12,
            'intron_6': 13,
            'exon_7': 14,
            'intron_7': 15,
            'exon_8': 16,
            'intron_8': 17,
            'exon_9': 18,
            'intron_9': 19,
            'three_prime_UTR': 20
        },
        'KIR2DL1':
        {
            'five_prime_UTR': 1,
            'exon_1': 2,
            'intron_1': 3,
            'exon_2': 4,
            'intron_2': 5,
            'exon_3': 6,
            'intron_3': 7,
            'exon_4': 8,
            'intron_4': 9,
            'exon_5': 10,
            'intron_5': 11,
            'exon_6': 12,
            'intron_6': 13,
            'exon_7': 14,
            'intron_7': 15,
            'exon_8': 16,
            'intron_8': 17,
            'exon_9': 18,
            'intron_9': 19,
            'three_prime_UTR': 20
        },
        'KIR2DP1':
        {
            'five_prime_UTR': 1,
            'exon_1': 2,
            'intron_1': 3,
            'exon_2': 4,
            'intron_2': 5,
            'exon_3': 6,
            'intron_3': 7,
            'exon_4': 8,
            'intron_4': 9,
            'exon_5': 10,
            'intron_5': 11,
            'exon_6': 12,
            'intron_6': 13,
            'exon_7': 14,
            'intron_7': 15,
            'exon_8': 16,
            'intron_8': 17,
            'exon_9': 18,
            'intron_9': 19,
            'three_prime_UTR': 20
        },
        'KIR2DL5A':
        {
            'five_prime_UTR': 1,
            'exon_1': 2,
            'intron_1': 3,
            'exon_2': 4,
            'intron_2': 5,
            'exon_3': 6,
            'intron_3': 7,
            'exon_4': 8,
            'intron_4': 9,
            'exon_5': 10,
            'intron_5': 11,
            'exon_6': 12,
            'intron_6': 13,
            'exon_7': 14,
            'intron_7': 15,
            'exon_8': 16,
            'intron_8': 17,
            'exon_9': 18,
            'intron_9': 19,
            'three_prime_UTR': 20
        },
        'KIR2DS4':
        {
            'five_prime_UTR': 1,
            'exon_1': 2,
            'intron_1': 3,
            'exon_2': 4,
            'intron_2': 5,
            'exon_3': 6,
            'intron_3': 7,
            'exon_4': 8,
            'intron_4': 9,
            'exon_5': 10,
            'intron_5': 11,
            'exon_6': 12,
            'intron_6': 13,
            'exon_7': 14,
            'intron_7': 15,
            'exon_8': 16,
            'intron_8': 17,
            'exon_9': 18,
            'intron_9': 19,
            'three_prime_UTR': 20
        },
        'HLA-DPB1':
        {
            'five_prime_UTR': 1,
            'exon_1': 2,
            'intron_1': 3,
            'exon_2': 4,
            'intron_2': 5,
            'exon_3': 6,
            'intron_3': 7,
            'exon_4': 8,
            'intron_4': 9,
            'exon_5': 10,
            'three_prime_UTR': 11
        },
        'KIR2DS2':
        {
            'five_prime_UTR': 1,
            'exon_1': 2,
            'intron_1': 3,
            'exon_2': 4,
            'intron_2': 5,
            'exon_3': 6,
            'intron_3': 7,
            'exon_4': 8,
            'intron_4': 9,
            'exon_5': 10,
            'intron_5': 11,
            'exon_6': 12,
            'intron_6': 13,
            'exon_7': 14,
            'intron_7': 15,
            'exon_8': 16,
            'intron_8': 17,
            'exon_9': 18,
            'intron_9': 19,
            'three_prime_UTR': 20
        },
        'KIR3DP1':
        {
            'five_prime_UTR': 1,
            'exon_1': 2,
            'intron_1': 3,
            'exon_2': 4,
            'intron_2': 5,
            'exon_3': 6,
            'intron_3': 7,
            'exon_4': 8,
            'intron_4': 9,
            'exon_5': 10,
            'intron_5': 11,
            'exon_6': 12,
            'intron_6': 13,
            'exon_7': 14,
            'intron_7': 15,
            'exon_8': 16,
            'intron_8': 17,
            'exon_9': 18,
            'intron_9': 19,
            'three_prime_UTR': 20
        },
        'HLA-DRB4':
        {
            'five_prime_UTR': 1,
            'exon_1': 2,
            'intron_1': 3,
            'exon_2': 4,
            'intron_2': 5,
            'exon_3': 6,
            'intron_3': 7,
            'exon_4': 8,
            'intron_4': 9,
            'exon_5': 10,
            'intron_5': 11,
            'exon_6': 12,
            'three_prime_UTR': 13
        },
        'KIR2DS5':
        {
            'five_prime_UTR': 1,
            'exon_1': 2,
            'intron_1': 3,
            'exon_2': 4,
            'intron_2': 5,
            'exon_3': 6,
            'intron_3': 7,
            'exon_4': 8,
            'intron_4': 9,
            'exon_5': 10,
            'intron_5': 11,
            'exon_6': 12,
            'intron_6': 13,
            'exon_7': 14,
            'intron_7': 15,
            'exon_8': 16,
            'intron_8': 17,
            'exon_9': 18,
            'intron_9': 19,
            'three_prime_UTR': 20
        },
        'HLA-DQA1':
        {
            'five_prime_UTR': 1,
            'exon_1': 2,
            'intron_1': 3,
            'exon_2': 4,
            'intron_2': 5,
            'exon_3': 6,
            'intron_3': 7,
            'exon_4': 8,
            'three_prime_UTR': 9
        },
        'HLA-DRB3':
        {
            'five_prime_UTR': 1,
            'exon_1': 2,
            'intron_1': 3,
            'exon_2': 4,
            'intron_2': 5,
            'exon_3': 6,
            'intron_3': 7,
            'exon_4': 8,
            'intron_4': 9,
            'exon_5': 10,
            'intron_5': 11,
            'exon_6': 12,
            'three_prime_UTR': 13
        },
        'KIR2DS3':
        {
            'five_prime_UTR': 1,
            'exon_1': 2,
            'intron_1': 3,
            'exon_2': 4,
            'intron_2': 5,
            'exon_3': 6,
            'intron_3': 7,
            'exon_4': 8,
            'intron_4': 9,
            'exon_5': 10,
            'intron_5': 11,
            'exon_6': 12,
            'intron_6': 13,
            'exon_7': 14,
            'intron_7': 15,
            'exon_8': 16,
            'intron_8': 17,
            'exon_9': 18,
            'intron_9': 19,
            'three_prime_UTR': 20
        },
        'KIR3DL1':
        {
            'five_prime_UTR': 1,
            'exon_1': 2,
            'intron_1': 3,
            'exon_2': 4,
            'intron_2': 5,
            'exon_3': 6,
            'intron_3': 7,
            'exon_4': 8,
            'intron_4': 9,
            'exon_5': 10,
            'intron_5': 11,
            'exon_6': 12,
            'intron_6': 13,
            'exon_7': 14,
            'intron_7': 15,
            'exon_8': 16,
            'intron_8': 17,
            'exon_9': 18,
            'intron_9': 19,
            'three_prime_UTR': 20
        },
        'HLA-A':
        {
            'five_prime_UTR': 1,
            'exon_1': 2,
            'intron_1': 3,
            'exon_2': 4,
            'intron_2': 5,
            'exon_3': 6,
            'intron_3': 7,
            'exon_4': 8,
            'intron_4': 9,
            'exon_5': 10,
            'intron_5': 11,
            'exon_6': 12,
            'intron_6': 13,
            'exon_7': 14,
            'intron_7': 15,
            'exon_8': 16,
            'three_prime_UTR': 17
        },
        'HLA-DRB5':
        {
            'five_prime_UTR': 1,
            'exon_1': 2,
            'intron_1': 3,
            'exon_2': 4,
            'intron_2': 5,
            'exon_3': 6,
            'intron_3': 7,
            'exon_4': 8,
            'intron_4': 9,
            'exon_5': 10,
            'intron_5': 11,
            'exon_6': 12,
            'three_prime_UTR': 13
        },
        'KIR2DL4':
        {
            'five_prime_UTR': 1,
            'exon_1': 2,
            'intron_1': 3,
            'exon_2': 4,
            'intron_2': 5,
            'exon_3': 6,
            'intron_3': 7,
            'exon_4': 8,
            'intron_4': 9,
            'exon_5': 10,
            'intron_5': 11,
            'exon_6': 12,
            'intron_6': 13,
            'exon_7': 14,
            'intron_7': 15,
            'exon_8': 16,
            'intron_8': 17,
            'exon_9': 18,
            'intron_9': 19,
            'three_prime_UTR': 20
        },
        'HLA-DQB1':
        {
            'five_prime_UTR': 1,
            'exon_1': 2,
            'intron_1': 3,
            'exon_2': 4,
            'intron_2': 5,
            'exon_3': 6,
            'intron_3': 7,
            'exon_4': 8,
            'intron_4': 9,
            'exon_5': 10,
            'intron_5': 11,
            'exon_6': 12,
            'three_prime_UTR': 13
        },
        'KIR3DL2':
        {
            'five_prime_UTR': 1,
            'exon_1': 2,
            'intron_1': 3,
            'exon_2': 4,
            'intron_2': 5,
            'exon_3': 6,
            'intron_3': 7,
            'exon_4': 8,
            'intron_4': 9,
            'exon_5': 10,
            'intron_5': 11,
            'exon_6': 12,
            'intron_6': 13,
            'exon_7': 14,
            'intron_7': 15,
            'exon_8': 16,
            'intron_8': 17,
            'exon_9': 18,
            'intron_9': 19,
            'three_prime_UTR': 20
        },
        'HLA-B':
        {
            'five_prime_UTR': 1,
            'exon_1': 2,
            'intron_1': 3,
            'exon_2': 4,
            'intron_2': 5,
            'exon_3': 6,
            'intron_3': 7,
            'exon_4': 8,
            'intron_4': 9,
            'exon_5': 10,
            'intron_5': 11,
            'exon_6': 12,
            'intron_6': 13,
            'exon_7': 14,
            'three_prime_UTR': 15
        },
        'KIR3DS1':
        {
            'five_prime_UTR': 1,
            'exon_1': 2,
            'intron_1': 3,
            'exon_2': 4,
            'intron_2': 5,
            'exon_3': 6,
            'intron_3': 7,
            'exon_4': 8,
            'intron_4': 9,
            'exon_5': 10,
            'intron_5': 11,
            'exon_6': 12,
            'intron_6': 13,
            'exon_7': 14,
            'intron_7': 15,
            'exon_8': 16,
            'intron_8': 17,
            'exon_9': 18,
            'intron_9': 19,
            'three_prime_UTR': 20
        },
        'KIR2DL5B':
        {
            'five_prime_UTR': 1,
            'exon_1': 2,
            'intron_1': 3,
            'exon_2': 4,
            'intron_2': 5,
            'exon_3': 6,
            'intron_3': 7,
            'exon_4': 8,
            'intron_4': 9,
            'exon_5': 10,
            'intron_5': 11,
            'exon_6': 12,
            'intron_6': 13,
            'exon_7': 14,
            'intron_7': 15,
            'exon_8': 16,
            'intron_8': 17,
            'exon_9': 18,
            'intron_9': 19,
            'three_prime_UTR': 20
        },
        'HLA-DRB1':
        {
            'five_prime_UTR': 1,
            'exon_1': 2,
            'intron_1': 3,
            'exon_2': 4,
            'intron_2': 5,
            'exon_3': 6,
            'intron_3': 7,
            'exon_4': 8,
            'intron_4': 9,
            'exon_5': 10,
            'intron_5': 11,
            'exon_6': 12,
            'three_prime_UTR': 13
        },
        'KIR3DL3':
        {
            'five_prime_UTR': 1,
            'exon_1': 2,
            'intron_1': 3,
            'exon_2': 4,
            'intron_2': 5,
            'exon_3': 6,
            'intron_3': 7,
            'exon_4': 8,
            'intron_4': 9,
            'exon_5': 10,
            'intron_5': 11,
            'exon_6': 12,
            'intron_6': 13,
            'exon_7': 14,
            'intron_7': 15,
            'exon_8': 16,
            'intron_8': 17,
            'exon_9': 18,
            'intron_9': 19,
            'three_prime_UTR': 20
        },
        'KIR2DS1':
        {
            'five_prime_UTR': 1,
            'exon_1': 2,
            'intron_1': 3,
            'exon_2': 4,
            'intron_2': 5,
            'exon_3': 6,
            'intron_3': 7,
            'exon_4': 8,
            'intron_4': 9,
            'exon_5': 10,
            'intron_5': 11,
            'exon_6': 12,
            'intron_6': 13,
            'exon_7': 14,
            'intron_7': 15,
            'exon_8': 16,
            'intron_8': 17,
            'exon_9': 18,
            'intron_9': 19,
            'three_prime_UTR': 20
        },
        'HLA-DPA1':
        {
            'five_prime_UTR': 1,
            'exon_1': 2,
            'intron_1': 3,
            'exon_2': 4,
            'intron_2': 5,
            'exon_3': 6,
            'intron_3': 7,
            'exon_4': 8,
            'three_prime_UTR': 9
        },
        'HLA-DRA':
        {
            'five_prime_UTR': 1,
            'exon_1': 2,
            'intron_1': 3,
            'exon_2': 4,
            'intron_2': 5,
            'exon_3': 6,
            'intron_3': 7,
            'exon_4': 8,
            'three_prime_UTR': 9
        },
        'HLA-C':
        {
            'five_prime_UTR': 1,
            'exon_1': 2,
            'intron_1': 3,
            'exon_2': 4,
            'intron_2': 5,
            'exon_3': 6,
            'intron_3': 7,
            'exon_4': 8,
            'intron_4': 9,
            'exon_5': 10,
            'intron_5': 11,
            'exon_6': 12,
            'intron_6': 13,
            'exon_7': 14,
            'intron_7': 15,
            'exon_8': 16,
            'three_prime_UTR': 17
        }
    }


def get_structmax():
    return {'HLA-DRA': 9, 'KIR2DL2': 20, 'KIR2DP1': 20, 'KIR2DL5A': 20, 'KIR2DS4': 20, 'HLA-DPA1': 9, 'HLA-DQA1': 9, 'HLA-DPB1': 11, 'KIR2DS2': 20, 'KIR3DP1': 20, 'HLA-DRB4': 13, 'KIR2DL1': 20, 'KIR2DS5': 20, 'HLA-DRB3': 13, 'KIR2DS3': 20, 'KIR3DL1': 20, 'HLA-A': 17, 'HLA-DRB5': 13, 'KIR2DL4': 20, 'HLA-DQB1': 13, 'KIR3DL2': 20, 'HLA-B': 15, 'KIR3DS1': 20, 'KIR2DL5B': 20, 'HLA-DRB1': 13, 'KIR3DL3': 20, 'KIR2DS1': 20, 'HLA-C': 17}


def get_structorder():
    return {
        'KIR2DL2':
        {
            1: 'five_prime_UTR',
            2: 'exon_1',
            3: 'intron_1',
            4: 'exon_2',
            5: 'intron_2',
            6: 'exon_3',
            7: 'intron_3',
            8: 'exon_4',
            9: 'intron_4',
            10: 'exon_5',
            11: 'intron_5',
            12: 'exon_6',
            13: 'intron_6',
            14: 'exon_7',
            15: 'intron_7',
            16: 'exon_8',
            17: 'intron_8',
            18: 'exon_9',
            19: 'intron_9',
            20: 'three_prime_UTR'
        },
        'KIR2DL1':
        {
            1: 'five_prime_UTR',
            2: 'exon_1',
            3: 'intron_1',
            4: 'exon_2',
            5: 'intron_2',
            6: 'exon_3',
            7: 'intron_3',
            8: 'exon_4',
            9: 'intron_4',
            10: 'exon_5',
            11: 'intron_5',
            12: 'exon_6',
            13: 'intron_6',
            14: 'exon_7',
            15: 'intron_7',
            16: 'exon_8',
            17: 'intron_8',
            18: 'exon_9',
            19: 'intron_9',
            20: 'three_prime_UTR'
        },
        'KIR2DP1':
        {
            1: 'five_prime_UTR',
            2: 'exon_1',
            3: 'intron_1',
            4: 'exon_2',
            5: 'intron_2',
            6: 'exon_3',
            7: 'intron_3',
            8: 'exon_4',
            9: 'intron_4',
            10: 'exon_5',
            11: 'intron_5',
            12: 'exon_6',
            13: 'intron_6',
            14: 'exon_7',
            15: 'intron_7',
            16: 'exon_8',
            17: 'intron_8',
            18: 'exon_9',
            19: 'intron_9',
            20: 'three_prime_UTR'
        },
        'KIR2DL5A':
        {
            1: 'five_prime_UTR',
            2: 'exon_1',
            3: 'intron_1',
            4: 'exon_2',
            5: 'intron_2',
            6: 'exon_3',
            7: 'intron_3',
            8: 'exon_4',
            9: 'intron_4',
            10: 'exon_5',
            11: 'intron_5',
            12: 'exon_6',
            13: 'intron_6',
            14: 'exon_7',
            15: 'intron_7',
            16: 'exon_8',
            17: 'intron_8',
            18: 'exon_9',
            19: 'intron_9',
            20: 'three_prime_UTR'
        },
        'KIR2DS4':
        {
            1: 'five_prime_UTR',
            2: 'exon_1',
            3: 'intron_1',
            4: 'exon_2',
            5: 'intron_2',
            6: 'exon_3',
            7: 'intron_3',
            8: 'exon_4',
            9: 'intron_4',
            10: 'exon_5',
            11: 'intron_5',
            12: 'exon_6',
            13: 'intron_6',
            14: 'exon_7',
            15: 'intron_7',
            16: 'exon_8',
            17: 'intron_8',
            18: 'exon_9',
            19: 'intron_9',
            20: 'three_prime_UTR'
        },
        'HLA-DPB1':
        {
            1: 'five_prime_UTR',
            2: 'exon_1',
            3: 'intron_1',
            4: 'exon_2',
            5: 'intron_2',
            6: 'exon_3',
            7: 'intron_3',
            8: 'exon_4',
            9: 'intron_4',
            10: 'exon_5',
            11: 'three_prime_UTR'
        },
        'HLA-DRA':
        {
            1: 'five_prime_UTR',
            2: 'exon_1',
            3: 'intron_1',
            4: 'exon_2',
            5: 'intron_2',
            6: 'exon_3',
            7: 'intron_3',
            8: 'exon_4',
            9: 'three_prime_UTR'
        },
        'KIR2DS2':
        {
            1: 'five_prime_UTR',
            2: 'exon_1',
            3: 'intron_1',
            4: 'exon_2',
            5: 'intron_2',
            6: 'exon_3',
            7: 'intron_3',
            8: 'exon_4',
            9: 'intron_4',
            10: 'exon_5',
            11: 'intron_5',
            12: 'exon_6',
            13: 'intron_6',
            14: 'exon_7',
            15: 'intron_7',
            16: 'exon_8',
            17: 'intron_8',
            18: 'exon_9',
            19: 'intron_9',
            20: 'three_prime_UTR'
        },
        'KIR3DP1':
        {
            1: 'five_prime_UTR',
            2: 'exon_1',
            3: 'intron_1',
            4: 'exon_2',
            5: 'intron_2',
            6: 'exon_3',
            7: 'intron_3',
            8: 'exon_4',
            9: 'intron_4',
            10: 'exon_5',
            11: 'intron_5',
            12: 'exon_6',
            13: 'intron_6',
            14: 'exon_7',
            15: 'intron_7',
            16: 'exon_8',
            17: 'intron_8',
            18: 'exon_9',
            19: 'intron_9',
            20: 'three_prime_UTR'
        },
        'HLA-DRB4':
        {
            1: 'five_prime_UTR',
            2: 'exon_1',
            3: 'intron_1',
            4: 'exon_2',
            5: 'intron_2',
            6: 'exon_3',
            7: 'intron_3',
            8: 'exon_4',
            9: 'intron_4',
            10: 'exon_5',
            11: 'intron_5',
            12: 'exon_6',
            13: 'three_prime_UTR'
        },
        'KIR2DS5':
        {
            1: 'five_prime_UTR',
            2: 'exon_1',
            3: 'intron_1',
            4: 'exon_2',
            5: 'intron_2',
            6: 'exon_3',
            7: 'intron_3',
            8: 'exon_4',
            9: 'intron_4',
            10: 'exon_5',
            11: 'intron_5',
            12: 'exon_6',
            13: 'intron_6',
            14: 'exon_7',
            15: 'intron_7',
            16: 'exon_8',
            17: 'intron_8',
            18: 'exon_9',
            19: 'intron_9',
            20: 'three_prime_UTR'
        },
        'HLA-DQA1':
        {
            1: 'five_prime_UTR',
            2: 'exon_1',
            3: 'intron_1',
            4: 'exon_2',
            5: 'intron_2',
            6: 'exon_3',
            7: 'intron_3',
            8: 'exon_4',
            9: 'three_prime_UTR'
        },
        'HLA-DRB3':
        {
            1: 'five_prime_UTR',
            2: 'exon_1',
            3: 'intron_1',
            4: 'exon_2',
            5: 'intron_2',
            6: 'exon_3',
            7: 'intron_3',
            8: 'exon_4',
            9: 'intron_4',
            10: 'exon_5',
            11: 'intron_5',
            12: 'exon_6',
            13: 'three_prime_UTR'
        },
        'KIR2DS3':
        {
            1: 'five_prime_UTR',
            2: 'exon_1',
            3: 'intron_1',
            4: 'exon_2',
            5: 'intron_2',
            6: 'exon_3',
            7: 'intron_3',
            8: 'exon_4',
            9: 'intron_4',
            10: 'exon_5',
            11: 'intron_5',
            12: 'exon_6',
            13: 'intron_6',
            14: 'exon_7',
            15: 'intron_7',
            16: 'exon_8',
            17: 'intron_8',
            18: 'exon_9',
            19: 'intron_9',
            20: 'three_prime_UTR'
        },
        'KIR3DL1':
        {
            1: 'five_prime_UTR',
            2: 'exon_1',
            3: 'intron_1',
            4: 'exon_2',
            5: 'intron_2',
            6: 'exon_3',
            7: 'intron_3',
            8: 'exon_4',
            9: 'intron_4',
            10: 'exon_5',
            11: 'intron_5',
            12: 'exon_6',
            13: 'intron_6',
            14: 'exon_7',
            15: 'intron_7',
            16: 'exon_8',
            17: 'intron_8',
            18: 'exon_9',
            19: 'intron_9',
            20: 'three_prime_UTR'
        },
        'HLA-A':
        {
            1: 'five_prime_UTR',
            2: 'exon_1',
            3: 'intron_1',
            4: 'exon_2',
            5: 'intron_2',
            6: 'exon_3',
            7: 'intron_3',
            8: 'exon_4',
            9: 'intron_4',
            10: 'exon_5',
            11: 'intron_5',
            12: 'exon_6',
            13: 'intron_6',
            14: 'exon_7',
            15: 'intron_7',
            16: 'exon_8',
            17: 'three_prime_UTR'
        },
        'HLA-DRB5':
        {
            1: 'five_prime_UTR',
            2: 'exon_1',
            3: 'intron_1',
            4: 'exon_2',
            5: 'intron_2',
            6: 'exon_3',
            7: 'intron_3',
            8: 'exon_4',
            9: 'intron_4',
            10: 'exon_5',
            11: 'intron_5',
            12: 'exon_6',
            13: 'three_prime_UTR'
        },
        'KIR2DL4':
        {
            1: 'five_prime_UTR',
            2: 'exon_1',
            3: 'intron_1',
            4: 'exon_2',
            5: 'intron_2',
            6: 'exon_3',
            7: 'intron_3',
            8: 'exon_4',
            9: 'intron_4',
            10: 'exon_5',
            11: 'intron_5',
            12: 'exon_6',
            13: 'intron_6',
            14: 'exon_7',
            15: 'intron_7',
            16: 'exon_8',
            17: 'intron_8',
            18: 'exon_9',
            19: 'intron_9',
            20: 'three_prime_UTR'
        },
        'HLA-DQB1':
        {
            1: 'five_prime_UTR',
            2: 'exon_1',
            3: 'intron_1',
            4: 'exon_2',
            5: 'intron_2',
            6: 'exon_3',
            7: 'intron_3',
            8: 'exon_4',
            9: 'intron_4',
            10: 'exon_5',
            11: 'intron_5',
            12: 'exon_6',
            13: 'three_prime_UTR'
        },
        'KIR3DL2':
        {
            1: 'five_prime_UTR',
            2: 'exon_1',
            3: 'intron_1',
            4: 'exon_2',
            5: 'intron_2',
            6: 'exon_3',
            7: 'intron_3',
            8: 'exon_4',
            9: 'intron_4',
            10: 'exon_5',
            11: 'intron_5',
            12: 'exon_6',
            13: 'intron_6',
            14: 'exon_7',
            15: 'intron_7',
            16: 'exon_8',
            17: 'intron_8',
            18: 'exon_9',
            19: 'intron_9',
            20: 'three_prime_UTR'
        },
        'HLA-B':
        {
            1: 'five_prime_UTR',
            2: 'exon_1',
            3: 'intron_1',
            4: 'exon_2',
            5: 'intron_2',
            6: 'exon_3',
            7: 'intron_3',
            8: 'exon_4',
            9: 'intron_4',
            10: 'exon_5',
            11: 'intron_5',
            12: 'exon_6',
            13: 'intron_6',
            14: 'exon_7',
            15: 'three_prime_UTR'
        },
        'KIR3DS1':
        {
            1: 'five_prime_UTR',
            2: 'exon_1',
            3: 'intron_1',
            4: 'exon_2',
            5: 'intron_2',
            6: 'exon_3',
            7: 'intron_3',
            8: 'exon_4',
            9: 'intron_4',
            10: 'exon_5',
            11: 'intron_5',
            12: 'exon_6',
            13: 'intron_6',
            14: 'exon_7',
            15: 'intron_7',
            16: 'exon_8',
            17: 'intron_8',
            18: 'exon_9',
            19: 'intron_9',
            20: 'three_prime_UTR'
        },
        'KIR2DL5B':
        {
            1: 'five_prime_UTR',
            2: 'exon_1',
            3: 'intron_1',
            4: 'exon_2',
            5: 'intron_2',
            6: 'exon_3',
            7: 'intron_3',
            8: 'exon_4',
            9: 'intron_4',
            10: 'exon_5',
            11: 'intron_5',
            12: 'exon_6',
            13: 'intron_6',
            14: 'exon_7',
            15: 'intron_7',
            16: 'exon_8',
            17: 'intron_8',
            18: 'exon_9',
            19: 'intron_9',
            20: 'three_prime_UTR'
        },
        'HLA-DRB1':
        {
            1: 'five_prime_UTR',
            2: 'exon_1',
            3: 'intron_1',
            4: 'exon_2',
            5: 'intron_2',
            6: 'exon_3',
            7: 'intron_3',
            8: 'exon_4',
            9: 'intron_4',
            10: 'exon_5',
            11: 'intron_5',
            12: 'exon_6',
            13: 'three_prime_UTR'
        },
        'KIR3DL3':
        {
            1: 'five_prime_UTR',
            2: 'exon_1',
            3: 'intron_1',
            4: 'exon_2',
            5: 'intron_2',
            6: 'exon_3',
            7: 'intron_3',
            8: 'exon_4',
            9: 'intron_4',
            10: 'exon_5',
            11: 'intron_5',
            12: 'exon_6',
            13: 'intron_6',
            14: 'exon_7',
            15: 'intron_7',
            16: 'exon_8',
            17: 'intron_8',
            18: 'exon_9',
            19: 'intron_9',
            20: 'three_prime_UTR'
        },
        'KIR2DS1':
        {
            1: 'five_prime_UTR',
            2: 'exon_1',
            3: 'intron_1',
            4: 'exon_2',
            5: 'intron_2',
            6: 'exon_3',
            7: 'intron_3',
            8: 'exon_4',
            9: 'intron_4',
            10: 'exon_5',
            11: 'intron_5',
            12: 'exon_6',
            13: 'intron_6',
            14: 'exon_7',
            15: 'intron_7',
            16: 'exon_8',
            17: 'intron_8',
            18: 'exon_9',
            19: 'intron_9',
            20: 'three_prime_UTR'
        },
        'HLA-DPA1':
        {
            1: 'five_prime_UTR',
            2: 'exon_1',
            3: 'intron_1',
            4: 'exon_2',
            5: 'intron_2',
            6: 'exon_3',
            7: 'intron_3',
            8: 'exon_4',
            9: 'three_prime_UTR'
        },
        'HLA-C':
        {
            1: 'five_prime_UTR',
            2: 'exon_1',
            3: 'intron_1',
            4: 'exon_2',
            5: 'intron_2',
            6: 'exon_3',
            7: 'intron_3',
            8: 'exon_4',
            9: 'intron_4',
            10: 'exon_5',
            11: 'intron_5',
            12: 'exon_6',
            13: 'intron_6',
            14: 'exon_7',
            15: 'intron_7',
            16: 'exon_8',
            17: 'three_prime_UTR'
        }
    }

