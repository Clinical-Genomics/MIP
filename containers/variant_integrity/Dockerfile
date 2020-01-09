################## BASE IMAGE ######################

FROM clinicalgenomics/mip:latest

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:latest"
LABEL version="1"
LABEL software="variant_integrity"
LABEL software.version="0.0.4"
LABEL extra.binaries="variant_integrity"

################## MAINTAINER ######################
MAINTAINER Clinical-Genomics/MIP

RUN conda install pip python=3.7

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

RUN pip install --no-cache-dir variant_integrity==0.0.4

WORKDIR /data/
