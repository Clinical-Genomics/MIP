################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base-image="clinicalgenomics/mip_base:2.1"
LABEL version="7"
LABEL software="htslib"
LABEL software.version="1.15.1"
LABEL extra.binaries="bcftools, bgzip, samtools, tabix"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install bcftools=1.15.1 htslib=1.15.1 samtools=1.15.1 && \
    /opt/conda/bin/conda clean -ya

WORKDIR /data/
