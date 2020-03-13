################## BASE IMAGE ######################

FROM clinicalgenomics/mip:2.0

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:2.0"
LABEL version="1"
LABEL software="genmod"
LABEL software.version="3.7.3"
LABEL extra.binaries="genmod"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install pip python=3.7

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

RUN pip install --no-cache-dir genmod==3.7.3

WORKDIR /data/
