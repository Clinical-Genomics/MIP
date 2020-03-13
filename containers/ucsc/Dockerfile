################## BASE IMAGE ######################

FROM clinicalgenomics/mip:2.0

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:2.0"
LABEL version="1"
LABEL software="ucsc-wigtobigwig, ucsc-bedtobigbed"
LABEL software.version="357"
LABEL extra.binaries="wigToBigWig, bedToBigBed"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install ucsc-wigtobigwig=357=h35c10e6_3 ucsc-bedtobigbed=357=h35c10e6_3

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

WORKDIR /data/
