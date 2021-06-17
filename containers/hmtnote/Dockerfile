################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="1"
LABEL software="HmtNote"
LABEL software.version="0.7.2"
LABEL extra.binaries="HmtNote"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install pip python=3.7

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

RUN pip install --no-cache-dir hmtnote==0.7.2

RUN hmtnote dump

WORKDIR /data/
