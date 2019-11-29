################## BASE IMAGE ######################

FROM clinicalgenomics/mip:latest

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:latest"
LABEL version="1"
LABEL software="htslib"
LABEL software.version="1.9"
LABEL extra.binaries="bcftools, bgzip, samtools, tabix"

################## MAINTAINER ######################
MAINTAINER Clinical-Genomics/MIP

RUN conda install bcftools=1.9=ha228f0b_4 htslib=1.9=hc238db4_4 samtools=1.9=h8571acd_11

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

WORKDIR /data/
