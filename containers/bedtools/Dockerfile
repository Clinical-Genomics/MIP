################## BASE IMAGE ######################

FROM clinicalgenomics/mip:2.0

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:2.0"
LABEL version="1"
LABEL software="bedtools"
LABEL software.version="2.29.0"
LABEL extra.binaries="bedtools"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install bedtools=2.29.0

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

WORKDIR /data/
