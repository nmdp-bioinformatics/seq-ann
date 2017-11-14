#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pymysql
from Bio import SeqIO
from BioSQL import BioSeqDatabase
from BioSQL.BioSeq import DBSeq

from seqann.sequence_annotation import BioSeqAnn
import time
import datetime
import argparse
import sys


def main():
    """This is run if file is directly executed, but not if imported as
    module. Having this in a separate function  allows importing the file
    into interactive python, and still able to execute the
    function for testing"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file",
                        required=True,
                        help="input file",
                        type=str)

    parser.add_argument("-l", "--locus",
                        required=True,
                        help="HLA locus",
                        type=str)

    parser.add_argument("-k", "--kir",
                        required=False,
                        help="Bool for KIR",
                        default=False,
                        type=bool)

    parser.add_argument("-s", "--server",
                        required=False,
                        help="BioSQL server",
                        default=False,
                        type=bool)

    parser.add_argument("-v", "--verbose",
                        help="Option for running in verbose",
                        default=False,
                        type=bool)

    args = parser.parse_args()
    fastafile = args.file
    loc = args.locus
    kir = args.kir
    serv = args.server

    server = None
    if serv:
        server = BioSeqDatabase.open_database(driver="pymysql", user="root",
                                              passwd="", host="localhost",
                                              db="bioseqdb")

    seqann = BioSeqAnn(verbose=True, server=server, kir=kir)
    for seq in SeqIO.parse(fastafile, "fasta"):
        ann = seqann.annotate(seq, loc)
        print('{:*^20} {:^20} {:*^20}'.format("", str(seq.description), ""))
        l = 0
        for f in ann.annotation:
            if isinstance(ann.annotation[f], DBSeq):
                print(f, ann.method, str(ann.annotation[f]), sep="\t")
                l += len(ann.annotation[f])
            else:
                print(f, ann.method, str(ann.annotation[f].seq), sep="\t")
                l += len(ann.annotation[f].seq)
        print("")

    if serv:
        server.close()

if __name__ == '__main__':
    """The following will be run if file is executed directly,
    but not if imported as a module"""
    main()


