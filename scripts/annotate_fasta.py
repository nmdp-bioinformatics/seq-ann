#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pymysql
from Bio import SeqIO
from BioSQL import BioSeqDatabase
from BioSQL.BioSeq import DBSeq

from seqann import BioSeqAnn
import time
import datetime
import argparse
import sys
import logging


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
                        help="Locus",
                        type=str)

    parser.add_argument("-k", "--kir",
                        help="Option for running with KIR",
                        action='store_true')

    parser.add_argument("-s", "--server",
                        help="Option for running with a server",
                        action='store_true')

    parser.add_argument("-v", "--verbose",
                        help="Option for running in verbose",
                        action='store_true')

    args = parser.parse_args()
    fastafile = args.file
    locus = args.locus

    verbose = False
    if args.verbose:
        verbose = True

    verbose = False
    if args.verbose:
        verbose = True

    kir = False
    if args.kir:
        kir = True

    serv = False
    if args.server:
        serv = True

    if verbose:
        logging.basicConfig(format='%(asctime)s - %(name)-35s - %(levelname)-5s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p',
                            level=logging.INFO)

    server = None
    if serv:
        server = BioSeqDatabase.open_database(driver="pymysql", user="root",
                                              passwd="", host="localhost",
                                              db="bioseqdb")

    seqann = BioSeqAnn(verbose=True, kir=kir)
    for seq in SeqIO.parse(fastafile, "fasta"):
        ann = seqann.annotate(seq, locus=locus)
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


