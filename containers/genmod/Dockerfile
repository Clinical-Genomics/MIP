################## BASE IMAGE ######################

FROM clinicalgenomics/mip:latest

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:latest"
LABEL version="1"
LABEL software="genmod"
LABEL software.version="3.7.3"
LABEL extra.binaries="genmod"

################## MAINTAINER ######################
MAINTAINER Clinical-Genomics/MIP

RUN conda install pip python=3.7

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

RUN pip install --no-cache-dir genmod==3.7.3

WORKDIR /data/
