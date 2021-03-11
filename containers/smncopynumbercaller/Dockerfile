################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="3"
LABEL software="SMNCopyNumberCaller"
LABEL software.version="v1.1.1"
LABEL extra.binaries="smn_caller"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install pip python=3.7 numpy scipy pysam statsmodels

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

## Download
RUN wget --no-verbose --no-check-certificate https://github.com/Illumina/SMNCopyNumberCaller/archive/v1.1.1.zip -O /opt/conda/share/smncopynumbercaller-v.1.1.1.zip

## Extract
RUN unzip -o /opt/conda/share/smncopynumbercaller-v.1.1.1.zip -d /opt/conda/share/

## Move to smn_caller directory
RUN cd /opt/conda/share/SMNCopyNumberCaller-1.1.1 && \
    chmod a+x smn_caller.py

RUN ln --symbolic --force /opt/conda/share/SMNCopyNumberCaller-1.1.1/smn_caller.py /opt/conda/bin/smn_caller.py

RUN ln --symbolic --force /opt/conda/share/SMNCopyNumberCaller-1.1.1/data /opt/conda/bin/data

WORKDIR /data/
