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
import os
import string
import random as rand
from typing import GenericMeta
from datetime import datetime, date
from six import integer_types, iteritems
from Bio.SeqFeature import SeqFeature


def get_seqfeat(seqrecord):

    n = 3 if len(seqrecord.features) >= 3 else len(seqrecord.features)
    fiveutr = [["five_prime_UTR", seqrecord.features[i]] for i in range(0, n) if seqrecord.features[i].type != "source"
               and seqrecord.features[i].type != "CDS" and isinstance(seqrecord.features[i], SeqFeature)
               and not seqrecord.features[i].qualifiers]
    feats = [[str(feat.type + "_" + feat.qualifiers['number'][0]), feat]
             for feat in seqrecord.features if feat.type != "source"
             and feat.type != "CDS" and isinstance(feat, SeqFeature)
             and 'number' in feat.qualifiers]
    threeutr = [["three_prime_UTR", seqrecord.features[i]] for i in range(len(seqrecord.features)-1, len(seqrecord.features)) if seqrecord.features[i].type != "source"
                and seqrecord.features[i].type != "CDS" and isinstance(seqrecord.features[i], SeqFeature)
                and not seqrecord.features[i].qualifiers]

    feat_list = fiveutr + feats + threeutr
    annotation = {k[0]: k[1] for k in feat_list}
    return(annotation)


# TODO: change name to get_featseq
def get_features(seqrecord):
    # print("^^^^^^^^^^^^^^^^^^^^^^^^^^")
    # print("get_features")
    # print(seqrecord)
    # print(seqrecord.features)
    # print("^^^^^^^^^^^^^^^^^^^^^^^^^^")
    # TODO: Make sure UTR's have type of UTR
    fiveutr = [["five_prime_UTR", seqrecord.features[i].extract(seqrecord.seq)] for i in range(0, 3) if seqrecord.features[i].type != "source"
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
    #print(annotation)
    return(annotation)


def randomid(N=6):
    random_id = ''.join(rand.choices(string.ascii_uppercase + string.digits,
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
