from Bio import SeqIO
from BioSQL import BioSeqDatabase

hla_names = {}
with open('Allelelist.3290.txt', 'r') as f:
    for line in f:
        line = line.rstrip()
        accession, name = line.split(" ")
        hla_names.update({accession: name})
f.close()

print("Loaded allele names")

seq_list = SeqIO.parse("hla.dat", "imgt")

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

print("Loaded IMGT dat file")
server = BioSeqDatabase.open_database(driver="pymysql", user="root",
                                      passwd="", host="localhost",
                                      db="bioseqdb")
for locus in new_seqs:
    dbname = "3290_" + locus
    dbdescription = "IMGT/HLA 3.29.0 " + locus
    db = server.new_database(dbname, description=dbdescription)
    count = db.load(new_seqs[locus])
    print("Loaded ", count, " for ", locus)
    server.commit()
