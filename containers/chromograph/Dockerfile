################## BASE IMAGE ######################

FROM clinicalgenomics/mip:2.0

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:2.0"
LABEL version="10"
LABEL software="chromograph"
LABEL software.version="1.0.1"
LABEL extra.binaries="chromograph"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install pip python=3.7 matplotlib

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

## Download release
RUN wget --no-verbose https://github.com/mikaell/chromograph/archive/v1.0.1.zip -O /opt/conda/share/chromograph-1.0.1.zip

## Extract
RUN unzip -o /opt/conda/share/chromograph-1.0.1.zip -d /opt/conda/share/ && \
    rm /opt/conda/share/chromograph-1.0.1.zip

## Move to chromograph directory
RUN cd /opt/conda/share/chromograph-1.0.1 && \
    python -m pip install .

WORKDIR /data/
