################## BASE IMAGE ######################

FROM clinicalgenomics/mip:latest

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:latest"
LABEL version="1"
LABEL software="sambamba"
LABEL software.version="0.6.8"
LABEL extra.binaries="sambamba"

################## MAINTAINER ######################
MAINTAINER Clinical-Genomics/MIP

RUN conda install sambamba=0.6.8=h682856c_0 ldc=1.11.0=hb2c9046_0

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

WORKDIR /data/
