################## BASE IMAGE ######################

FROM clinicalgenomics/mip:2.0

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:2.0"
LABEL version="2"
LABEL software="chromograph"
LABEL software.version="91f0988"
LABEL extra.binaries="chromograph"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install pip python=3.7 matplotlib

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

## Clone git repository
RUN git clone https://github.com/mikaell/chromograph.git /opt/conda/share/chromograph

## Move to chromograph directory
RUN cd /opt/conda/share/chromograph && \
    git checkout 91f0988 && \
    python -m pip install --editable .

WORKDIR /data/
