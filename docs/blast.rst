.. highlight:: shell

======================
Creating blastn files
======================

.. note:: Make sure blastn is properly installed before running!

1) Download the allele list and the ``_gen`` and ``_nuc`` fasta files for each locus

2) Create the blast files by running the **ngs-imgt-db** perl script

.. code-block:: console

    $ ngs-imgt-db -i /path/to/imgt/files -o /output/dir

3) Add the new blast files to the seqann/data/blast directory and check them in.