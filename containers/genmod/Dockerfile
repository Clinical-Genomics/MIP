################## BASE IMAGE ######################

FROM condaforge/mambaforge:23.1.0-4

################## METADATA ######################

LABEL base_image="condaforge/mambaforge:23.1.0-4"
LABEL version="2"
LABEL software="genmod"
LABEL software.version="3.8.0"
LABEL extra.binaries="genmod"
LABEL maintainer="Clinical-Genomics/MIP"

RUN mamba install -y -c conda-forge -c bioconda python=3.10 pip pytabix && \
    /opt/conda/bin/conda clean -ya && \
    pip install --no-cache-dir genmod==3.8.0

