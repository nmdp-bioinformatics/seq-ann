.. highlight:: shell

======================
Debugging
======================

In order to debug your code you must first create a logging instance in your script.

.. code-block:: python3

    import logging
    logging.basicConfig(format='%(asctime)s - %(name)-35s - %(levelname)-5s - %(message)s',
                datefmt='%m/%d/%Y %I:%M:%S %p',
                level=logging.INFO)

Once you have the logging set up, you can then pass a dictionary to the debug parameter. The
dictionary should have keys representing a main process and values for how verbose the output
should be for those proecesses. If a process is not present in the debug dictionary, then
no logging will be generated for that process.

.. code-block:: python3

    seqann = BioSeqAnn(debug={"seqann": 5, "align":1, "refdata":0, "seq_search": 5, "gfe": 4})

Instead of using the debug parameter, you can use the verbose and verbosity parameters to see the same
level logging for each process. Use these parameters when you want to see the most logging possible for each process.

.. code-block:: python3

    seqann = BioSeqAnn(verbose=True, verbosity=5)

