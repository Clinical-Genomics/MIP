################## BASE IMAGE ######################

FROM clinicalgenomics/mip:latest

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:latest"
LABEL version="1"
LABEL software="SMNCopyNumberCaller"
LABEL software.version="4b2c1ad"
LABEL extra.binaries="smn_caller"

################## MAINTAINER ######################
MAINTAINER Clinical-Genomics/MIP

RUN conda install pip python=3.7 numpy scipy pysam statsmodels

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

## Clone git repository
RUN git clone https://github.com/Illumina/SMNCopyNumberCaller.git /opt/conda/share/SMNCopyNumberCaller

## Move to smn_caller directory
RUN cd /opt/conda/share/SMNCopyNumberCaller && \
    chmod a+x smn_caller.py

RUN ln --symbolic --force /opt/conda/share/SMNCopyNumberCaller/smn_caller.py /opt/conda/bin/smn_caller.py

WORKDIR /data/
