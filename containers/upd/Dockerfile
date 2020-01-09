################## BASE IMAGE ######################

FROM clinicalgenomics/mip:latest

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:latest"
LABEL version="1"
LABEL software="upd"
LABEL software.version="0.1"
LABEL extra.binaries="upd"

################## MAINTAINER ######################
MAINTAINER Clinical-Genomics/MIP

RUN conda install pip python=3.7

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

## Download UPD
RUN wget --no-verbose https://github.com/bjhall/upd/archive/0.1.zip -O /opt/conda/share/UPD-0.1.zip

## Extract
RUN unzip -o /opt/conda/share/UPD-0.1.zip -d /opt/conda/share/

## Move to UPD directory
RUN cd /opt/conda/share/upd-0.1 && \
    python -m pip install --editable .

WORKDIR /data/
