################## BASE IMAGE ######################

FROM clinicalgenomics/mip:2.0

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:2.0"
LABEL version="1"
LABEL software="vt"
LABEL software.version="2015.11.10"
LABEL extra.binaries="vt"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install vt=2015.11.10

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

WORKDIR /data/
