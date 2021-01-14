################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="2"
LABEL software="utilities"
LABEL extra.binaries="gtf2bed, pigz"
LABEL maintainer="Clinical-Genomics/MIP"

## Install conda packages
RUN conda install pigz

## Clean up after conda
RUN /opt/conda/bin/conda clean -ya

## Install scripts
## gtf2bed
RUN wget --no-verbose -P /opt/conda/share/ https://github.com/ExpressionAnalysis/ea-utils/archive/1.04.807.zip && \
    unzip -o /opt/conda/share/1.04.807.zip -d /opt/conda/share/ && \
    ln --symbolic --force /opt/conda/share/ea-utils-1.04.807/clipper/gtf2bed /opt/conda/bin/gtf2bed && \
    rm /opt/conda/share/1.04.807.zip

WORKDIR /data/
