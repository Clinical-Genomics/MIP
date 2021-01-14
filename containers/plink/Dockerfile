################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="1"
LABEL software="plink"
LABEL software.version="1.90b3.35"
LABEL extra.binaries="plink"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install plink2=1.90b3.35

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

WORKDIR /data/
