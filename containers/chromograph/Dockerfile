################## BASE IMAGE ######################

FROM clinicalgenomics/mip:2.0

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:2.0"
LABEL version="4"
LABEL software="chromograph"
LABEL software.version="0.2.1"
LABEL extra.binaries="chromograph"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install pip python=3.7 matplotlib

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

## Clone git repository
RUN wget --no-verbose https://github.com/mikaell/chromograph/archive/v0.2.1.zip -O /opt/conda/share/chromograph-0.2.1.zip

## Extract
RUN unzip -o /opt/conda/share/chromograph-0.2.1.zip -d /opt/conda/share/ && \
    rm /opt/conda/share/chromograph-0.2.1.zip

## Move to chromograph directory
RUN cd /opt/conda/share/chromograph-0.2.1 && \
    python -m pip install --editable .

WORKDIR /data/
