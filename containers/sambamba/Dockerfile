################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="2"
LABEL software="sambamba"
LABEL software.version="0.8.0"
LABEL extra.binaries="sambamba"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install sambamba=0.8.0=h984e79f_0

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

WORKDIR /data/
