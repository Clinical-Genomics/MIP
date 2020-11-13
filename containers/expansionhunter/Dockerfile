################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="2"
LABEL software="expansionhunter"
LABEL software.version="3.2.2"
LABEL extra.binaries="expansionhunter"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install -c bioconda expansionhunter=3.2.2=hd1df264_0

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

WORKDIR /data/
