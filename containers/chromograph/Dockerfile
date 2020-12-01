################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="8"
LABEL software="chromograph"
LABEL software.version="0.3.3"
LABEL extra.binaries="chromograph"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install pip python=3.7 matplotlib

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

## Download release
RUN wget --no-verbose https://github.com/mikaell/chromograph/archive/v0.3.3.zip -O /opt/conda/share/chromograph-0.3.3.zip

## Extract
RUN unzip -o /opt/conda/share/chromograph-0.3.3.zip -d /opt/conda/share/ && \
    rm /opt/conda/share/chromograph-0.3.3.zip

## Move to chromograph directory
RUN cd /opt/conda/share/chromograph-0.3.3 && \
    python -m pip install --editable .

WORKDIR /data/
