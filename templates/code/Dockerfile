################## BASE IMAGE ######################

FROM clinicalgenomics/mip:2.0

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:2.0"
LABEL version="1"
LABEL software="SOFTWARE"
LABEL software.version="VERSION"
LABEL extra.binaries=""
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install SOFTWARE=VERSION

## Clean up after conda
RUN /opt/conda/bin/conda clean -ya

WORKDIR /data/
