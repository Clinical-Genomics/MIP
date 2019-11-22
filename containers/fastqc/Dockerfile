################## BASE IMAGE ######################

FROM clinicalgenomics/mip:latest

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:latest"
LABEL version="1"
LABEL software="fastqc"
LABEL software.version="0.11.8"
LABEL extra.binaries="fastqc"

################## MAINTAINER ######################
MAINTAINER Clinical-Genomics/MIP

RUN conda install fastqc=0.11.8=1

WORKDIR /data/

