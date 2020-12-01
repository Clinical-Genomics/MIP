################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="1"
LABEL software="rseqc"
LABEL software.version="3.0.1"
LABEL extra.binaries="rseqc"
LABEL maintainer="Clinical-Genomics/MIP"

RUN pip install --no-cache-dir rseqc==3.0.1

WORKDIR /data/
