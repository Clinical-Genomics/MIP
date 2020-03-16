################## BASE IMAGE ######################

FROM clinicalgenomics/mip:2.0

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:2.0"
LABEL version="1"
LABEL software="vcfanno"
LABEL software.version="0.3.2"
LABEL extra.binaries="vcfanno"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install vcfanno=0.3.2=0

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

WORKDIR /data/
