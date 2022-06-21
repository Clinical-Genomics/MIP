################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="3"
LABEL software="trim_galore"
LABEL software.version="0.6.7"
LABEL extra.binaries=""
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install trim-galore=0.6.7 fontconfig && \
    /opt/conda/bin/conda clean -ya

WORKDIR /data/
