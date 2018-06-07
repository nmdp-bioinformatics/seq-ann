FROM ubuntu:17.10
MAINTAINER Mike Halagan <mhalagan@nmdp.org>

RUN apt-get update -q \
    && apt-get dist-upgrade -qy \
    && apt-get install -qyy wget curl build-essential cpp git \
    && apt-get -qyy install python3.6 python3-pip python3-dev python3-setuptools uwsgi-plugin-python3  python-mysqldb python3-mysql.connector

RUN apt-get install python3.6-dev -qy

RUN cd opt/ && git clone https://github.com/nmdp-bioinformatics/SeqAnn && cd SeqAnn \
    && curl https://bootstrap.pypa.io/get-pip.py | python3.6 \
    && pip install --upgrade pip

RUN cd opt/SeqAnn && pip install -r requirements.txt \
    && python3.6 setup.py install 

RUN apt-get install clustalo -y
RUN apt-get install ncbi-blast+ -y

