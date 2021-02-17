################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="2"
LABEL software="rseqc"
LABEL software.version="3.0.1"
LABEL extra.binaries="rseqc"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install rseqc=3.0.1 ucsc-wigtobigwig
RUN conda clean -ya

WORKDIR /data/
