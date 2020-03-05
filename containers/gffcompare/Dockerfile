################## BASE IMAGE ######################

FROM clinicalgenomics/mip:latest

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:latest"
LABEL version="1"
LABEL software="gffcompare"
LABEL software.version="0.11.2"
LABEL extra.binaries="gffcompare"

################## MAINTAINER ######################
MAINTAINER Clinical-Genomics/MIP

RUN conda install --yes gffcompare=0.11.2

## Clean up after conda
RUN /opt/conda/bin/conda clean -ay

WORKDIR /data/
