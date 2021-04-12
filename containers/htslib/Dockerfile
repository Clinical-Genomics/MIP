################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="4"
LABEL software="htslib"
LABEL software.version="1.12"
LABEL extra.binaries="bcftools, bgzip, samtools, tabix"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install bcftools=1.12=h45bccc9_1 htslib=1.12=h9093b5e_1 samtools=1.12=h9aed4be_1

## Clean up after conda
RUN /opt/conda/bin/conda clean -ya

WORKDIR /data/
