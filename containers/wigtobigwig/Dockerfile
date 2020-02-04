################## BASE IMAGE ######################

FROM clinicalgenomics/mip:latest

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:latest"
LABEL version="1"
LABEL software="ucsc-wigtobigwig"
LABEL software.version="357"
LABEL extra.binaries="wigtobigwig"

################## MAINTAINER ######################
MAINTAINER Clinical-Genomics/MIP

RUN conda install ucsc-wigtobigwig=357=h35c10e6_3

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

WORKDIR /data/
