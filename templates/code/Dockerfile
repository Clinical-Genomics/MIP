################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="1"
LABEL software="SOFTWARE"
LABEL software.version="VERSION"
LABEL extra.binaries=""
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install SOFTWARE=VERSION

## Clean up after conda
RUN /opt/conda/bin/conda clean -ya

WORKDIR /data/
