################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.0

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.0"
LABEL version="5"
LABEL software="htslib"
LABEL software.version="1.10.2"
LABEL extra.binaries="bcftools, bgzip, samtools, tabix"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install bcftools=1.10.2=hd2cd319_0 htslib=1.10.2=h78d89cc_0 samtools=1.10=h9402c20_2

## Clean up after conda
RUN /opt/conda/bin/conda clean -ya

WORKDIR /data/
