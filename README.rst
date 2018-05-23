===============================
SeqAnn
===============================


.. image:: https://img.shields.io/pypi/v/seqann.svg
        :target: https://pypi.python.org/pypi/seqann

.. image:: https://img.shields.io/travis/nmdp-bioinformatics/seqann.svg
        :target: https://travis-ci.org/nmdp-bioinformatics/seqann

.. image:: https://readthedocs.org/projects/seqann/badge/?version=latest
        :target: https://seqann.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://pyup.io/repos/github/nmdp-bioinformatics/gfe/shield.svg
     :target: https://pyup.io/repos/github/nmdp-bioinformatics/seqann/
     :alt: Updates


Sequence Annotation


* Free software: LGPL 3.0
* Documentation: https://seqann.readthedocs.io.


Overview
---------

The ``BioSeqAnn`` class in the ``seqann`` package allows 
users to easily annotate the gene features in a consensus sequence.
No parameters are required when creating a ``BioSeqAnn`` object, but annotations can be
created significantly faster when using a ``BioSeqDatabase``. 
The lastest `hla.dat`_ file is downloaded and parsed when a ``BioSeqDatabase`` is not provided.
A BioSQL database containing all of IPD-IMGT/HLA is available on DockerHub_ and can be
run on any machine that has docker installed. Below are the list of parameters that
can be use to initalize a ``BioSeqAnn`` object.

.. table::
    :widths: 10 10 10 50

    +-------------+-------------------+---------+-------------------------------------------------------------------------------+
    | Parameter   | Type              | Default | Description                                                                   |
    +=============+===================+=========+===============================================================================+
    | server      | ``BioSeqDatabase``| None    | A BioSQL database containing all of the sequence data from IPD-IMGT/HLA.      |
    +-------------+-------------------+---------+-------------------------------------------------------------------------------+
    | dbversion   | ``str``           | Latest  | The IPD-IMGT/HLA or KIR database release.                                     |
    +-------------+-------------------+---------+-------------------------------------------------------------------------------+
    | datfile     | ``str``           | None    | The IPD-IMGT/HLA or KIR dat file to use in place of the **server** parameter. |
    +-------------+-------------------+---------+-------------------------------------------------------------------------------+
    | kir         | ``bool``          | False   | Flag for indicating the input sequences are from the KIR gene system.         |
    +-------------+-------------------+---------+-------------------------------------------------------------------------------+
    | align       | ``bool``          | False   | Flag for producing the alignments along with the annotations.                 |
    +-------------+-------------------+---------+-------------------------------------------------------------------------------+
    | verbose     | ``bool``          | False   | Flag for running in verbose mode.                                             |
    +-------------+-------------------+---------+-------------------------------------------------------------------------------+
    | verbosity   | ``int``           | None    | Numerical value to indicate how verbose the output will be in verbose mode.   |
    +-------------+-------------------+---------+-------------------------------------------------------------------------------+

Install
------------

.. code-block:: shell
    
    pip install seqann


Usage
---------

With default parameters:

.. code-block:: python3

    from seqann import BioSeqAnn
    from Bio import SeqIO

    # ** If you don't have a copy of the hla.dat
    # ** file it will download it
    seqann = BioSeqAnn()
    for seq in SeqIO.parse(input_seq, "fasta"):
        # Annotate sequence
        ann = seqann.annotate(seq, "HLA-A")
        # Loop through features in annotation
        for feat in ann.annotation:
            print(feat, ann.annotation[feat], sep="\t")



With BioSQL DB connection:

.. code-block:: python3

    from Bio import SeqIO
    from seqann import BioSeqAnn
    from BioSQL import BioSeqDatabase

    server = BioSeqDatabase.open_database(driver="pymysql", user="root",
                                          passwd="", host="localhost",
                                          db="bioseqdb")
    seqann = BioSeqAnn(server=server)
    for seq in SeqIO.parse(input_seq, "fasta"):
        # Annotate sequence
        ann = seqann.annotate(seq, "HLA-A")
        # Loop through features in annotation
        for feat in ann.annotation:
            print(feat, ann.annotation[feat], sep="\t")


Dependencies
------------
* `Clustal Omega`_ 1.2.0 or higher
* `Python 3.6`_
* blastn_

.. _DockerHub: http://google.com
.. _`GitHub page`: http://google.com
.. _`hla.dat`: http://google.com
.. _`Python 3.6`: https://www.python.org/downloads
.. _`Clustal Omega`: http://www.clustal.org/omega/
.. _blastn: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage

