################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="2"
LABEL software="fastqc"
LABEL software.version="0.11.9"
LABEL extra.binaries="fastqc"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install fastqc=0.11.9=0

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

WORKDIR /data/
