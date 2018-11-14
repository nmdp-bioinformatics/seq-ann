.. highlight:: shell

======================
Testing
======================

.. warning:: Before running tests clustalo, blastn and all the required python packages must be installed!

To run all test simply execute the following command in the top directory of the SeqAnn repository.

.. code-block:: console

    $ python -m unittest tests


Different tests
---------------------

.. note:: If you don't have a imgt_biosql db running then not all of the test will run! 

* test_seqann
* test_align
* test_blast
* test_feature
* test_gfe
* test_refdata
* test_seqsearch
* test_util

Running specific tests
-----------------------

You can test a specific test by providing the full test path on the command line.

.. code-block:: console

    $ python -m unittest tests.test_seqann.TestBioSeqAnn.test_004_ambig

