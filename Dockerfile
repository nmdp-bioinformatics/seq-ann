FROM python:3.7
MAINTAINER NMDP Bioinformatics

RUN pip install biopython==1.74 PyMySQL==0.9.3 bson==0.5.8 requests==2.22.0

RUN apt-get update -q \
    && apt-get install python-mysqldb \
    && apt-get install clustalo -y \
	  && apt-get install ncbi-blast+ -y \
    && apt-get autoremove \
    && apt-get clean

RUN pip install seq-ann==1.0.5

