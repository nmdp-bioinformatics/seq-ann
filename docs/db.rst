.. highlight:: shell

.. _bio:

======================
BioSQL Database
======================


Running Container
-----------------

To get the IPD-IMGT/HLA BioSQL Docker image, run this command in your terminal:

.. code-block:: console

    $ docker pull nmdpbioinformatics/imgt_biosqldb

You can then run the datbase with the following command in your terminal:

.. code-block:: console

    $ docker run -d --name imgt_biosqldb -p 3306:3306 \
    	-e MYSQL_ROOT_PASSWORD=my-secret-pw nmdpbioinformatics/imgt_biosqldb

If you have a mysql database running locally already, then you can change the first number in the port mapping to something else. If
you change the port then remember to export the ``BIOSQLPORT`` environment variable to whatever you used. A password is required and the ``seqann`` package can access your password if you set the ``BIOSQLPASS`` environment variable.

Building Locally
----------------


If you want certain IPD-IMGT/HLA database verions that aren't loaded into the publicly available docker image, then you
can build the image locally with the database versions you want. The sources for IPD-IMGT/HLA BioSQL database can be downloaded from the `Github repo`_.

First clone the public repository:

.. code-block:: console

    $ git clone git://github.com/nmdp-bioinformatics/imgt_biosqldb

Then build the docker image and provide the ``RELEASES`` you want to use as an argument.

.. code-block:: console

    $ docker build -t imgt_biosqldb:3240-3250 --build-arg RELEASES="3240,3250" .


Once the image has finished building, you can run the database as described above.