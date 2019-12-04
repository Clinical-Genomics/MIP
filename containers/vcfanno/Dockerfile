################## BASE IMAGE ######################

FROM clinicalgenomics/mip:latest

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:latest"
LABEL version="1"
LABEL software="vcfanno"
LABEL software.version="0.3.2"
LABEL extra.binaries="vcfanno"

################## MAINTAINER ######################
MAINTAINER Clinical-Genomics/MIP

RUN conda install vcfanno=0.3.2=0

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

WORKDIR /data/
