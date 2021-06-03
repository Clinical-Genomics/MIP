################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="1"
LABEL software="gffcompare"
LABEL software.version="0.11.2"
LABEL extra.binaries="gffcompare"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install --yes gffcompare=0.11.2

## Clean up after conda
RUN /opt/conda/bin/conda clean -ay

WORKDIR /data/
