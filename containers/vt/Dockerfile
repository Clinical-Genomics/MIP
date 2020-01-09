################## BASE IMAGE ######################

FROM clinicalgenomics/mip:latest

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:latest"
LABEL version="1"
LABEL software="vt"
LABEL software.version="2015.11.10"
LABEL extra.binaries="vt"

################## MAINTAINER ######################
MAINTAINER Clinical-Genomics/MIP

RUN conda install vt=2015.11.10

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

WORKDIR /data/
