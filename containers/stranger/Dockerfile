################## BASE IMAGE ######################

FROM clinicalgenomics/mip:2.0

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:2.0"
LABEL version="2"
LABEL software="stranger"
LABEL software.version="0.7"
LABEL extra.binaries="stranger"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install pip python=3.7

## Clean up after conda
RUN /opt/conda/bin/conda clean -ya

RUN pip install --no-cache-dir stranger==0.7

WORKDIR /data/
