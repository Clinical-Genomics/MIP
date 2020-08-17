################## BASE IMAGE ######################

FROM clinicalgenomics/mip:2.0

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:2.0"
LABEL version="1"
LABEL software="trim_galore"
LABEL software.version="0.6.4"
LABEL extra.binaries=""
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install trim-galore=0.6.4=1

## Clean up after conda
RUN /opt/conda/bin/conda clean -ya

WORKDIR /data/
