# -*- coding: utf-8 -*-

import os
from Bio import SeqIO
from BioSQL import BioSeqDatabase
import urllib.request
import sys


def download_dat(db):
    url = 'https://raw.githubusercontent.com/ANHIG/IMGTHLA/' + db + '/hla.dat'
    dat = ".".join([db, "hla", "dat"])
    urllib.request.urlretrieve(url, dat)
    return dat


def download_allelelist(db):
    url = 'https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Allelelist.' + db + '.txt'
    alist = ".".join([db, "Allelelist", "txt"])
    urllib.request.urlretrieve(url, alist)
    return alist


dblist = ["".join([str(i), str("0")]) for i in range(326, 330)]
server = BioSeqDatabase.open_database(driver="pymysql", user="root",
                                      passwd="", host="localhost",
                                      db="bioseqdb")

for dbv in dblist:

    hladat = download_dat(dbv)
    allele_list = download_allelelist(dbv)

    hla_names = {}
    try:
        s = "," if dbv == "3260" or dbv == "3270" else " "
        with open(allele_list, 'r') as f:
            for line in f:
                line = line.rstrip()
                accession, name = line.split(s)
                hla_names.update({accession: name})
        f.close()
        print("Loaded allele names ", allele_list)
    except ValueError as err:
        print("Allelelist error: {0}".format(err))
        server.close()
        os.remove(hladat)
        os.remove(allele_list)
        sys.exit()

    try:
        seq_list = SeqIO.parse(hladat, "imgt")
    except ValueError as err:
        print("Read dat error: {0}".format(err))
        server.close()
        os.remove(hladat)
        os.remove(allele_list)
        sys.exit()

    new_seqs = {"A": [], "B": [], "C": [], "DRB1": [],
                "DQB1": [], "DRB3": [], "DRB4": [], "DRB5": [],
                "DQA1": [], "DPA1": [], "DPB1": []}

    for seq in seq_list:
        if seq.name in hla_names:
            loc, allele = hla_names[seq.name].split("*")
            if loc in new_seqs:
                hla_name = "HLA-" + hla_names[seq.name]
                seq.name = hla_name
                new_seqs[loc].append(seq)

    dbsp = list(dbv)
    descr = ".".join([dbsp[0], dbsp[1]+dbsp[2], dbsp[3]])
    print("Loaded IMGT dat file ", descr)

    for locus in new_seqs:
        dbname = dbv + "_" + locus
        dbdescription = "IMGT/HLA " + descr + " " + locus
        db = server.new_database(dbname, description=dbdescription)
        try:
            count = db.load(new_seqs[locus])
        except:
            print("Faild to load", sys.exc_info()[0])
            server.close()
            os.remove(hladat)
            os.remove(allele_list)
            sys.exit()

        print("Loaded ", count, " for ", locus)
        server.commit()

    os.remove(hladat)
    os.remove(allele_list)
    print("Finished ", descr)

server.close()
