################## BASE IMAGE ######################

FROM condaforge/mambaforge:23.1.0-4

################## METADATA ######################

LABEL base_image="condaforge/mambaforge:23.1.0-4"
LABEL version="2"
LABEL software="genmod"
LABEL software.version="3.8.1"
LABEL extra.binaries="genmod"
LABEL maintainer="Clinical-Genomics/MIP"

RUN mamba install -y -c conda-forge -c bioconda pip pytabix python=3.8 && \
    /opt/conda/bin/conda clean -ya && \
    pip install --no-cache-dir genmod==3.8.2
