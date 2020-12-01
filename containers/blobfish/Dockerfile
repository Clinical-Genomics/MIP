################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="1"
LABEL software="BlobFish"
LABEL software.version="0.0.2"
LABEL extra.binaries="BlobFish.py, BlobFish.R"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install --channel bioconda --channel conda-forge bioconductor-tximport python=2.7 bioconductor-deseq2=1.22.1 r-readr r-optparse

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

## Download BlobFish
RUN wget --no-verbose https://github.com/J35P312/BlobFish/archive/0.0.2.zip -O /opt/conda/share/blobfish-0.0.2.zip

## Extract
RUN unzip -o /opt/conda/share/blobfish-0.0.2.zip -d /opt/conda/share/

## Move to blobfish dir, make executable and link
RUN cd /opt/conda/share/BlobFish-0.0.2 && \
    chmod a+x BlobFish.py && \
    chmod a+x BlobFish.R && \
    ln --symbolic --force /opt/conda/share/BlobFish-0.0.2/BlobFish.py /opt/conda/bin/BlobFish.py && \
    ln --symbolic --force /opt/conda/share/BlobFish-0.0.2/BlobFish.R /opt/conda/bin/BlobFish.R

WORKDIR /data/
