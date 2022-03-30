################## BASE IMAGE ######################

FROM python:2.7-slim

################## METADATA ######################

LABEL base_image="python:2.7-slim"
LABEL version="2"
LABEL software="BootstrapAnn"
LABEL software.version="e557dd3"
LABEL extra.binaries="BootstrapAnn.py"
LABEL maintainer="Clinical-Genomics/MIP"

RUN apt-get update && apt-get install -y git
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN pip install numpy scipy

## Clone git repository
RUN git clone https://github.com/J35P312/BootstrapAnn.git /usr/local/BootstrapAnn

RUN cd /usr/local/BootstrapAnn && git checkout e557dd3

RUN cd /usr/local/BootstrapAnn && \
    chmod a+x BootstrapAnn.py

RUN ln --symbolic --force /usr/local/BootstrapAnn/BootstrapAnn.py /usr/local/bin/BootstrapAnn.py

WORKDIR /data/
