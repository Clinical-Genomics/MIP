################## BASE IMAGE ######################

FROM clinicalgenomics/mip:2.0

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:2.0"
LABEL version="1"
LABEL software="star"
LABEL software.version="2.7.4a"
LABEL extra.binaries="STAR"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install star=2.7.4a

## Clean up after conda
RUN /opt/conda/bin/conda clean -ya

WORKDIR /data/
