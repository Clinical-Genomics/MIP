################## BASE IMAGE ######################

FROM clinicalgenomics/mip:2.0

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:2.0"
LABEL version="1"
LABEL software="peddy"
LABEL software.version="0.4.3"
LABEL extra.binaries="peddy"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install python=2.7 peddy=0.4.3=py_0

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

WORKDIR /data/
