################## BASE IMAGE ######################

FROM clinicalgenomics/mip:2.0

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:2.0"
LABEL version="2"
LABEL software="ucsc-bedtobigbed, ucsc-gtftogenepred, ucsc-wigtobigwig"
LABEL software.version="377"
LABEL extra.binaries="bedToBigBed, gtfToGenePredi, wigToBigWig"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install ucsc-wigtobigwig=377=h199ee4e_0 ucsc-bedtobigbed=377=h199ee4e_0 ucsc-gtftogenepred=377=h199ee4e_1

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

WORKDIR /data/
