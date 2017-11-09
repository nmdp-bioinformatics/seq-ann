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

from Bio import SeqIO
from Bio import SearchIO
from Bio.Blast.Applications import NcbiblastnCommandline

from seqann.util import cleanup
from seqann.util import randomid
from seqann.models.blast import Blast
from seqann.models.reference_data import ReferenceData


def blastn(sequences, locus, nseqs, refdata=None, evalue=0.001):

    if not refdata:
        refdata = ReferenceData()

    file_id = str(randomid())
    input_fasta = file_id + ".fasta"
    output_xml = file_id + ".xml"
    SeqIO.write(sequences, input_fasta, "fasta")
    blastn_cline = NcbiblastnCommandline(query=input_fasta,
                                         db=refdata.blastdb,
                                         evalue=evalue, outfmt=5,
                                         out=output_xml)
    stdout, stderr = blastn_cline()
    loc = locus.split("-")[1]
    blast_qresult = SearchIO.read(output_xml, 'blast-xml')

    #   Delete files
    cleanup(file_id)

    if len(blast_qresult.hits) == 0:
        return Blast(failed=True)

    alleles = []
    full_sequences = []
    l = len(blast_qresult.hits) if nseqs > len(blast_qresult.hits) else nseqs

    if locus in refdata.hla_loci:
        alleles = ["HLA-" + blast_qresult[i].id.split("_")[0] for i in range(0, l)
                   if "HLA-" + blast_qresult[i].id.split("*")[0] == locus]

    # TODO: sort alleles by number of features they contain and evalue
    # Use biosql db if provided
    # otherwise use IMGT dat file
    if refdata.server_avail:
        db = refdata.server[refdata.dbversion + "_" + loc]
        full_sequences = [db.lookup(name=n) for n in alleles
                          if n in refdata.hla_names]
    else:
        full_sequences = [a for a in refdata.imgtdat
                          if a.description.split(",")[0] in alleles]
        full_sequences.reverse()

    #   Build Blast object
    blast_o = Blast(match_seqs=full_sequences, alleles=alleles)
    return blast_o

