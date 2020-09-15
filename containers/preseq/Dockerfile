################## BASE IMAGE ######################

FROM clinicalgenomics/mip:2.0

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:2.0"
LABEL version="1"
LABEL software="preseq"
LABEL software.version="2.0.3"
LABEL extra.binaries="preseq"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install preseq=2.0.3

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

WORKDIR /data/
