################## BASE IMAGE ######################

FROM clinicalgenomics/mip:latest

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:latest"
LABEL version="1"
LABEL software="salmon"
LABEL software.version="0.12.0"
LABEL extra.binaries="salmon"

################## MAINTAINER ######################
MAINTAINER Clinical-Genomics/MIP

RUN conda install salmon=0.12.0

## Clean up after conda
RUN /opt/conda/bin/conda clean -ya

WORKDIR /data/
