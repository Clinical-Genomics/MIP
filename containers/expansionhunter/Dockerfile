################## BASE IMAGE ######################

FROM clinicalgenomics/mip:2.0

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:2.0"
LABEL version="1"
LABEL software="expansionhunter"
LABEL software.version="3.1.2"
LABEL extra.binaries="expansionhunter"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install -c bioconda expansionhunter=3.1.2=hca43064_0

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

WORKDIR /data/
