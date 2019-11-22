################## BASE IMAGE ######################

FROM clinicalgenomics/mip:latest

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:latest"
LABEL version="1"
LABEL software="bwa"
LABEL software.version="0.7.17"
LABEL extra.binaries="bwa"

################## MAINTAINER ######################
MAINTAINER Clinical-Genomics/MIP

RUN conda install bwa=0.7.17=ha92aebf_3

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

WORKDIR /data/

