################## BASE IMAGE ######################

FROM clinicalgenomics/mip:2.0

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:2.0"
LABEL version="1"
LABEL software="plink"
LABEL software.version="1.90b3.35"
LABEL extra.binaries="plink"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install plink2=1.90b3.35

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

WORKDIR /data/
