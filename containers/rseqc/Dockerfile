################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="3"
LABEL software="rseqc"
LABEL software.version="4.0.0"
LABEL extra.binaries="rseqc"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install rseqc=4.0.0 ucsc-wigtobigwig
RUN conda clean -ya

WORKDIR /data/
