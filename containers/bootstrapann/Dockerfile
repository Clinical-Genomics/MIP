################## BASE IMAGE ######################

FROM python:2.7-slim

################## METADATA ######################

LABEL "base_image"="python:2.7-slim"
LABEL version="2"
LABEL software="BootstrapAnn"
LABEL software.version="e557dd3"
LABEL extra.binaries="BootstrapAnn.py"
LABEL maintainer="Clinical-Genomics/MIP"

RUN apt-get update && apt-get install -y --no-install-recommends git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    pip install --no-cache-dir numpy scipy && \
    git clone https://github.com/J35P312/BootstrapAnn.git /usr/local/BootstrapAnn

WORKDIR  /usr/local/BootstrapAnn

RUN git checkout e557dd3 && \
    chmod a+x BootstrapAnn.py && \
    ln --symbolic --force /usr/local/BootstrapAnn/BootstrapAnn.py /usr/local/bin/BootstrapAnn.py

WORKDIR /data/
