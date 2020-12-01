################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="2"
LABEL software="bwakit"
LABEL software.version="0.7.17"
LABEL extra.binaries="bwakit"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install bwakit=0.7.17

## Clean up after conda
RUN /opt/conda/bin/conda clean -ya

WORKDIR /data/
