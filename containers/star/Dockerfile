################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="2"
LABEL software="star"
LABEL software.version="2.7.8a"
LABEL extra.binaries="STAR"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install star=2.7.8a=h9ee0642_1

## Clean up after conda
RUN /opt/conda/bin/conda clean -ya

WORKDIR /data/
