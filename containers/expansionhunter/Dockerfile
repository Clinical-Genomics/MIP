################## BASE IMAGE ######################

FROM clinicalgenomics/mip:latest

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:latest"
LABEL version="1"
LABEL software="expansionhunter"
LABEL software.version="3.1.2"
LABEL extra.binaries="expansionhunter"

################## MAINTAINER ######################
MAINTAINER Clinical-Genomics/MIP

RUN conda install -c bioconda expansionhunter=3.1.2=hca43064_0

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

WORKDIR /data/
