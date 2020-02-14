################## BASE IMAGE ######################

FROM clinicalgenomics/mip:latest

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:latest"
LABEL version="2"
LABEL software="htslib"
LABEL software.version="1.10"
LABEL extra.binaries="bcftools, bgzip, samtools, tabix"

################## MAINTAINER ######################
MAINTAINER Clinical-Genomics/MIP

RUN conda install bcftools=1.10.2=hd2cd319_0 htslib=1.10.2=h78d89cc_0 samtools=1.10=h9402c20_2

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

WORKDIR /data/
