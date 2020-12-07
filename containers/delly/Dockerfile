################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="1"
LABEL software="delly"
LABEL software.version="0.8.1"
LABEL extra.binaries="delly"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install delly=0.8.1 && \
    conda clean -ya

WORKDIR /data/
