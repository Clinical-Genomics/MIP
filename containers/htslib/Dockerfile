################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base-image="clinicalgenomics/mip_base:2.1"
LABEL version="8"
LABEL software="htslib"
LABEL software.version="1.16"
LABEL extra.binaries="bcftools, bgzip, samtools, tabix"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install bcftools=1.16 htslib=1.16 samtools=1.16.1 && \
    /opt/conda/bin/conda clean -ya

WORKDIR /data/
