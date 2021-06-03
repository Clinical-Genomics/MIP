################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="2"
LABEL software="upd"
LABEL software.version="0.1.1"
LABEL extra.binaries="upd"
LABEL maintainer="Clinical-Genomics/MIP"

RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends \
        ca-certificates && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN conda install pip python=3.7

## Clean up after conda
RUN /opt/conda/bin/conda clean -ya

## Download UPD
RUN wget --no-verbose https://github.com/bjhall/upd/archive/refs/tags/0.1.1.zip -O /opt/conda/share/UPD-0.1.1.zip

## Extract
RUN unzip -o /opt/conda/share/UPD-0.1.1.zip -d /opt/conda/share/

## Move to UPD directory
RUN cd /opt/conda/share/upd-0.1.1 && \
    python -m pip install --no-cache-dir --editable .

WORKDIR /data/
