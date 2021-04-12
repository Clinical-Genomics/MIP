################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="2"
LABEL software="preseq"
LABEL software.version="3.1.2"
LABEL extra.binaries="preseq"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install preseq=3.1.2=h06ef8b0_1

## Clean up after conda
RUN /opt/conda/bin/conda clean -ya

WORKDIR /data/
