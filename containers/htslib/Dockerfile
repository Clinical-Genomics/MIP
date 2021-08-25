################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="6"
LABEL software="htslib"
LABEL software.version="1.13"
LABEL extra.binaries="bcftools, bgzip, samtools, tabix"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install bcftools=1.13=h3a49de5_0 htslib=1.13=h9093b5e_0 samtools=1.13=h8c37831_0

## Clean up after conda
RUN /opt/conda/bin/conda clean -ya

WORKDIR /data/
