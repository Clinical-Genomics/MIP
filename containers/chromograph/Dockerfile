################## BASE IMAGE ######################

FROM clinicalgenomics/mip:2.0

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:2.0"
LABEL version="5"
LABEL software="chromograph"
LABEL software.version="0.2.1"
LABEL extra.binaries="chromograph"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install pip python=3.7 matplotlib

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy



## Clone git repository
RUN git clone https://github.com/mikaell/chromograph.git /opt/conda/share/chromograph

## Move to chromograph directory
RUN cd /opt/conda/share/chromograph && \
    git checkout a3c5563 && \
    python -m pip install --editable .

WORKDIR /data/
