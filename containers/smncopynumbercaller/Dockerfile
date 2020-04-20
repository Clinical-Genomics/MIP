################## BASE IMAGE ######################

FROM clinicalgenomics/mip:2.0

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:2.0"
LABEL version="1"
LABEL software="SMNCopyNumberCaller"
LABEL software.version="1.0"
LABEL extra.binaries="smn_caller"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install pip python=3.7 numpy scipy pysam statsmodels

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

## Download
RUN wget --no-verbose https://github.com/Illumina/SMNCopyNumberCaller/archive/v1.0.zip -O /opt/conda/share/smncopynumbercaller-v.1.0.zip

## Extract
RUN unzip -o /opt/conda/share/smncopynumbercaller-v.1.0.zip -d /opt/conda/share/

## Move to smn_caller directory
RUN cd /opt/conda/share/SMNCopyNumberCaller-1.0 && \
    chmod a+x smn_caller.py

RUN ln --symbolic --force /opt/conda/share/SMNCopyNumberCaller-1.0/smn_caller.py /opt/conda/bin/smn_caller.py

RUN ln --symbolic --force /opt/conda/share/SMNCopyNumberCaller-1.0/data /opt/conda/bin/data

WORKDIR /data/
