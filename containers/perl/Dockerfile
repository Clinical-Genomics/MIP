################## BASE IMAGE ######################

FROM clinicalgenomics/mip:2.0

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:2.0"
LABEL version="1"
LABEL software="perl"
LABEL software.version="5.26"
LABEL extra.binaries="perl"
LABEL maintainer="Clinical-Genomics/MIP"

RUN locale-gen en_US.UTF-8  
ENV LANG en_US.UTF-8  
ENV LANGUAGE en_US:en  
ENV LC_ALL en_US.UTF-8

FROM perl:5.26

WORKDIR /data/
