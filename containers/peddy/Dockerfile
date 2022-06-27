################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="2"
LABEL software="peddy"
LABEL software.version="0.4.8"
LABEL extra.binaries="peddy"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install peddy=0.4.8 && \
    /opt/conda/bin/conda clean -ya

WORKDIR /data/
