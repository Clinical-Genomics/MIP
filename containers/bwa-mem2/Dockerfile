################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="1"
LABEL software="bwa-mem2"
LABEL software.version="2.0"
LABEL extra.binaries="bwa-mem2"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install bwa-mem2=2.0=he513fc3_0

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

WORKDIR /data/
