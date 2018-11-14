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
import logging
import re

# TODO: Move to util.py
has_hla = lambda x: True if re.search("HLA", x) else False


def blastn(sequences, locus, nseqs, kir=False,
           verbose=False, refdata=None, evalue=10):
    """
    Gets the a list of alleles that are the most similar to the input sequence

    :param sequences: The input sequence record.
    :type sequences: SeqRecord
    :param locus: The gene locus associated with the sequence.
    :type locus: ``str``
    :param nseqs: The incomplete annotation from a previous iteration.
    :type nseqs: ``int``
    :param evalue: The evalue to use (default = 10)
    :type evalue: ``int``
    :param kir: Run with KIR or not
    :type kir: ``bool``
    :param verbose: Run in versboe
    :type verbose: ``bool``
    :param refdata: An object with reference data
    :type refdata: :ref:`ref`
    :rtype: :ref:`bl`

    Example usage:

        >>> from Bio.Seq import Seq
        >>> from seqann.blast_cmd import blastn
        >>> sequence = Seq('AGAGACTCTCCCGAGGATTTCGTGTACCAGTTTAAGGCCATGTGCTACTTCACC')
        >>> blast = blastn(sequence, locus, nseqs)

    """
    logger = logging.getLogger("Logger." + __name__)

    if not refdata:
        refdata = ReferenceData()

    file_id = str(randomid())
    input_fasta = file_id + ".fasta"
    output_xml = file_id + ".xml"
    SeqIO.write(sequences, input_fasta, "fasta")
    blastn_cline = NcbiblastnCommandline(query=input_fasta,
                                         db=refdata.blastdb,
                                         evalue=evalue,
                                         outfmt=5,
                                         reward=1,
                                         penalty=-3,
                                         gapopen=5,
                                         gapextend=2,
                                         dust='yes',
                                         out=output_xml)

    stdout, stderr = blastn_cline()
    loc = locus
    if not kir:
        loc = locus.split("-")[1]

    blast_qresult = SearchIO.read(output_xml, 'blast-xml')

    #   Delete files
    cleanup(file_id)

    # TODO: Use logging
    if len(blast_qresult.hits) == 0:
        if verbose:
            logger.error("Failed blast! No hits!")
            logger.error(stderr)
        return Blast(failed=True)

    alleles = []
    full_sequences = []
    load_blast = 70 if nseqs < 70 else nseqs
    l = len(blast_qresult.hits) if load_blast > len(blast_qresult.hits) else load_blast

    # TODO: update all blast files to have HLA-
    if locus in refdata.hla_loci and not kir:
        alleles = [blast_qresult[i].id.split("_")[0] for i in range(0, l)
                   if blast_qresult[i].id.split("*")[0] == locus or "HLA-" + blast_qresult[i].id.split("*")[0] == locus]
        alleles = ["HLA-" + a if not has_hla(a) else a for a in alleles]
    if kir:
        alleles = [blast_qresult[i].id.split("_")[0] for i in range(0, l)
                   if blast_qresult[i].id.split("*")[0] == locus]

    if verbose:
        logger.info("Blast alleles: " + ",".join(alleles))

    # TODO: sort alleles by number of features they contain and evalue
    # Use biosql db if provided
    # otherwise use IMGT dat file
    final_seqs = []
    rmax = refdata.structure_max[locus]
    if refdata.server_avail:
        db = refdata.server[refdata.dbversion + "_" + loc]
        full_sequences = []
        for n in alleles:
            if n in refdata.hla_names:
                try:
                    seq = db.lookup(name=n)
                    full_sequences.append(seq)
                except:
                    logger.error("Allele doesnt exist in IMGT BioSQL DB!! "
                                 + n)
    else:
        if verbose:
            logger.info("Getting sequences from HLA.dat file")

        full_sequences = [refdata.hlaref[a] for a in alleles
                          if a in refdata.hlaref]

        for s in full_sequences:
            s.name = s.description.split(",")[0]

    i = 1
    last_seq = []
    max_f = 0
    added_max = False
    full_feats = False
    for s in full_sequences:

        fs = len([f.type for f in s.features
                  if not f.type in ['source', 'CDS']])
        if i <= nseqs:
            final_seqs.append(s)
            max_f = fs if fs > max_f else max_f

        if i <= nseqs and max_f < rmax:
            full_feats = True

        if(i >= nseqs and fs == max_f and not added_max):
            if len(last_seq) >= 10:
                last_seq.insert(3, s)
            else:
                last_seq.append(s)
            added_max = True

        if(fs > max_f and len(last_seq) < 10
           and i >= nseqs and len(last_seq) < 10):
            last_seq.append(s)

        i += 1

    if full_feats:
        for s in last_seq:
            final_seqs.append(s)

    #   Build Blast object
    blast_o = Blast(match_seqs=final_seqs, alleles=alleles)
    return blast_o


def get_locus(sequences, kir=False, verbose=False, refdata=None, evalue=10):
    """
    Gets the locus of the sequence by running blastn

    :param sequences: sequenences to blast
    :param kir: bool whether the sequences are KIR or not
    :rtype: ``str``

    Example usage:

        >>> from Bio.Seq import Seq
        >>> from seqann.blast_cmd import get_locus
        >>> sequence = Seq('AGAGACTCTCCCGAGGATTTCGTGTACCAGTTTAAGGCCATGTGCTACTTCACC')
        >>> locus = get_locus(sequence)

    """
    if not refdata:
        refdata = ReferenceData()

    file_id = str(randomid())
    input_fasta = file_id + ".fasta"
    output_xml = file_id + ".xml"
    SeqIO.write(sequences, input_fasta, "fasta")
    blastn_cline = NcbiblastnCommandline(query=input_fasta,
                                         db=refdata.blastdb,
                                         evalue=evalue,
                                         outfmt=5,
                                         reward=1,
                                         penalty=-3,
                                         gapopen=5,
                                         gapextend=2,
                                         dust='yes',
                                         out=output_xml)

    stdout, stderr = blastn_cline()

    blast_qresult = SearchIO.read(output_xml, 'blast-xml')

    #   Delete files
    cleanup(file_id)

    if len(blast_qresult.hits) == 0:
        return ''

    loci = []
    for i in range(0, 3):
        if kir:
            loci.append(blast_qresult[i].id.split("*")[0])
        else:
            loci.append(blast_qresult[i].id.split("*")[0])

    locus = set(loci)
    if len(locus) == 1:
        if has_hla(loci[0]) or kir:
            return loci[0]
        else:
            return "HLA-" + loci[0]
    else:
        return ''

