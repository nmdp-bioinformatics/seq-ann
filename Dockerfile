FROM python:3.7
MAINTAINER NMDP Bioinformatics

RUN apt-get update -q \
    && apt-get install python-mysqldb \
    && apt-get install clustalo -y \
	  && apt-get install ncbi-blast+ -y \
    && apt-get autoremove \
    && apt-get clean

RUN pip install -U pip seq-ann==1.1.0
