################## BASE IMAGE ######################

FROM clinicalgenomics/mip:latest

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:latest"
LABEL version="1"
LABEL software="rhocall"
LABEL software.version="0.5.1"
LABEL extra.binaries="rhocall"

################## MAINTAINER ######################
MAINTAINER Clinical-Genomics/MIP

RUN conda install pip matplotlib

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

## Download rhocall
RUN wget --no-verbose https://github.com/dnil/rhocall/archive/0.5.1.zip -O /opt/conda/share/rhocall-0.5.1.zip

## Extract
RUN unzip -o /opt/conda/share/rhocall-0.5.1.zip -d /opt/conda/share/

##  Move to rhocall directory
RUN cd /opt/conda/share/rhocall-0.5.1 && \
    pip install --no-cache-dir numpy Cython && \
    pip install --requirement requirements.txt && \
    python -m pip install --editable .

WORKDIR /data/
