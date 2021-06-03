################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="3"
LABEL software="expansionhunter"
LABEL software.version="4.0.2"
LABEL extra.binaries="expansionhunter"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install -c bioconda expansionhunter=4.0.2

## Clean up after conda
RUN /opt/conda/bin/conda clean -ya

WORKDIR /data/
