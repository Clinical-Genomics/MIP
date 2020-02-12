################## BASE IMAGE ######################

FROM clinicalgenomics/mip:latest

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:latest"
LABEL version="1"
LABEL software="utilities"
LABEL extra.binaries="pigz"

################## MAINTAINER ######################
MAINTAINER Clinical-Genomics/MIP

RUN conda install pigz

## Clean up after conda
RUN /opt/conda/bin/conda clean -ya

WORKDIR /data/
