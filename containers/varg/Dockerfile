################## BASE IMAGE ######################

FROM clinicalgenomics/mip:2.0

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:2.0"
LABEL version="1"
LABEL software="varg"
LABEL software.version="1.2.0"
LABEL extra.binaries="varg"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install pip python=3.7 cython=0.29.13 cyvcf2=0.11.5

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

## Move to chromograph directory
RUN pip install --no-cache-dir varg==1.2.0

WORKDIR /data/
