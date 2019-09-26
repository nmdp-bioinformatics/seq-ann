===============================
SeqAnn
===============================


.. image:: https://img.shields.io/travis/nmdp-bioinformatics/SeqAnn.svg
        :target: https://travis-ci.org/nmdp-bioinformatics/SeqAnn

.. image:: https://readthedocs.org/projects/seqann/badge/?version=latest
        :target: https://seqann.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://pyup.io/repos/github/nmdp-bioinformatics/SeqAnn/shield.svg
     :target: https://pyup.io/repos/github/nmdp-bioinformatics/SeqAnn/
     :alt: Updates

.. image:: https://img.shields.io/pypi/v/seqann.svg
        :target: https://pypi.python.org/pypi/seqann

.. image:: https://coveralls.io/repos/github/nmdp-bioinformatics/SeqAnn/badge.svg?branch=master
        :target: https://coveralls.io/github/nmdp-bioinformatics/SeqAnn?branch=master


Python package for annotating gene features


* Free software: LGPL 3.0
* Documentation: https://seqann.readthedocs.io.
* `Jupyter Notebook`_


Overview
---------

The ``seqann`` package allows 
users to annotate gene features in consensus sequences. Annotations can be created by passing consensus sequences to the ``annotate`` method in the
``BioSeqAnn`` class. No parameters are required when initalizing a ``BioSeqAnn`` class. However, annotations can be
created significantly faster when using a BioSQL database. When a BioSQL database is not provided the lastest `hla.dat`_ file is downloaded and parsed.
A BioSQL database containing all of IPD-IMGT/HLA is available on DockerHub_ and can be
run on any machine that has docker installed. 

Install
------------

.. code-block:: shell
    
    pip install seq-ann


Parameters
----------

Below are the list of parameters and the default values used when initalizing a ``BioSeqAnn`` object.

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
    | debug       | ``Dict``          | None    | A dictionary containing a process names as the key and verbosity as the value |
    +-------------+-------------------+---------+-------------------------------------------------------------------------------+

Usage
---------

To annotated a sequence initialize a new ``BioSeqAnn`` object and then pass the sequence to the
``annotate`` method. The sequence must be a Biopython ``Seq``. The locus of the sequence is not required but it will improve the accuracy of the annotation.


  The packages `ncbi-blast+` and `clustalo` are required to be installed on your system.


  Set variables to BioSQL host/port if using BioSQL.

.. code-block:: shell

    export BIOSQLHOST="localhost"
    export BIOSQLPORT=3306

.. code-block:: python3

    from seqann import BioSeqAnn
    seqann = BioSeqAnn()
    ann = seqann.annotate(sequence, "HLA-A")


The annotation of sequence can be done with or without providing a ``BioSeqDatabase``. To use a BioSQL database 
initialize a ``BioSeqDatabase`` with the parameters that match the database you have running. If you are 
running the imgt_biosqldb from DockerHub_ then the following parameters we be the same. 

.. code-block:: python3

    from seqann import BioSeqAnn
    from BioSQL import BioSeqDatabase
    server = BioSeqDatabase.open_database(driver="pymysql", user="root",
                                          passwd="my-secret-pw", host="localhost",
                                          db="bioseqdb", port=3306)
    seqann = BioSeqAnn(server=server)
    ann = seqann.annotate(sequence, "HLA-A")


You may need to set environment variables: *BIOSQLHOST* (e.g. "localhost") and *BIOSQLPORT*  (e.g. 3306)  to your docker instance.


Annotations
------------

.. code-block:: shell

            {
                 'complete_annotation': True,
                 'annotation': {'exon_1': SeqRecord(seq=Seq('AGAGACTCTCCCG', SingleLetterAlphabet()), id='HLA:HLA00630', name='HLA:HLA00630', description='HLA:HLA00630 DQB1*03:04:01 597 bp', dbxrefs=[]),
                                'exon_2': SeqRecord(seq=Seq('AGGATTTCGTGTACCAGTTTAAGGCCATGTGCTACTTCACCAACGGGACGGAGC...GAG', SingleLetterAlphabet()), id='HLA:HLA00630', name='HLA:HLA00630', description='HLA:HLA00630 DQB1*03:04:01 597 bp', dbxrefs=[]),
                                'exon_3': SeqRecord(seq=Seq('TGGAGCCCACAGTGACCATCTCCCCATCCAGGACAGAGGCCCTCAACCACCACA...ATG', SingleLetterAlphabet()), id='HLA:HLA00630', name='<unknown name>', description='HLA:HLA00630', dbxrefs=[])},
                 'features': {'exon_1': SeqFeature(FeatureLocation(ExactPosition(0), ExactPosition(13), strand=1), type='exon_1'),
                              'exon_2': SeqFeature(FeatureLocation(ExactPosition(13), ExactPosition(283), strand=1), type='exon_2')
                              'exon_3': SeqFeature(FeatureLocation(ExactPosition(283), ExactPosition(503), strand=1), type='exon_3')},
                 'method': 'nt_search and clustalo',
                 'gfe': 'HLA-Aw2-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-4',
                 'seq': SeqRecord(seq=Seq('AGAGACTCTCCCGAGGATTTCGTGTACCAGTTTAAGGCCATGTGCTACTTCACC...ATG', SingleLetterAlphabet()), id='HLA:HLA00630', name='HLA:HLA00630', description='HLA:HLA00630 DQB1*03:04:01 597 bp', dbxrefs=[])
            }


Once a sequence has been annotated the gene features and their corresponding sequences are available in the returned ``Annotation`` object. If a full
annotation is not able to be produced then nothing will be returned. Below is an example showing how the features can be accessed and printed out.

.. code-block:: python3

    ann = seqann.annotate(sequence, "HLA-A")
    for feat in ann.annotation:
        print(feat, ann.gfe, str(ann.annotation[feat].seq), sep="\t")



Dependencies
------------
* `Clustal Omega`_ 1.2.0 or higher
* `Python 3.6`_
* blastn_

.. _DockerHub: https://hub.docker.com/r/nmdpbioinformatics/imgt_biosqldb/
.. _`GitHub page`: http://google.com
.. _`hla.dat`: https://github.com/ANHIG/IMGTHLA
.. _`Python 3.6`: https://www.python.org/downloads
.. _`Clustal Omega`: http://www.clustal.org/omega/
.. _blastn: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _`Jupyter Notebook`: https://github.com/nmdp-bioinformatics/SeqAnn/blob/master/notebook/Examples.ipynb
