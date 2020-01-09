################## BASE IMAGE ######################

FROM clinicalgenomics/mip:latest

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:latest"
LABEL version="1"
LABEL software="peddy"
LABEL software.version="0.4.3"
LABEL extra.binaries="peddy"

################## MAINTAINER ######################
MAINTAINER Clinical-Genomics/MIP

RUN conda install python=2.7 peddy=0.4.3=py_0

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

WORKDIR /data/
