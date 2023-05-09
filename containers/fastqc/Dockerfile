################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="3"
LABEL software="fastqc"
LABEL software.version="0.12.1"
LABEL extra.binaries="fastqc"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install fastqc=0.12.1 && \
    /opt/conda/bin/conda clean -ya

WORKDIR /data/
