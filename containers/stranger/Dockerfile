################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="6"
LABEL software="stranger"
LABEL software.version="0.8.1"
LABEL extra.binaries="stranger"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install pip python=3.7 && \
    /opt/conda/bin/conda clean -ya && \
    pip install --no-cache-dir stranger==0.8.1

WORKDIR /data/
