################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="1"
LABEL software="tiddit"
LABEL software.version="2.12.1"
LABEL extra.binaries="TIDDIT.py"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install tiddit=2.12.1=py37h0498b6d_0

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

WORKDIR /data/
