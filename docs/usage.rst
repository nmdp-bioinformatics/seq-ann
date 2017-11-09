=====
Usage
=====

To use SeqAnn in a project::

	import seqann
	from Bio import SeqIO
	from BioSQL import BioSeqDatabase

	server = BioSeqDatabase.open_database(driver="pymysql", user="root",
	                                      passwd="", host="localhost",
	                                      db="bioseqdb")
	seqann = seqann.BioSeqAnn(server=server)
	for seq in SeqIO.parse(input_seq, "fasta"):
		annotation = seqann.annotate(seq, "HLA-A")
		for feat in annotation.annotation:
			print(feat, annotation.annotation[feat], sep="\t")

