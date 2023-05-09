################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base-image="clinicalgenomics/mip_base:2.1"
LABEL version="2"
LABEL software="HmtNote"
LABEL software.version="0.7.2"
LABEL extra.binaries="HmtNote, hmtvar_oct2022"
LABEL maintainer="Clinical-Genomics/MIP"

## Conda env installation + clean up
RUN conda install pip python=3.7 && \
    /opt/conda/bin/conda clean -ya

# Install HmtNote
RUN pip install --no-cache-dir hmtnote==0.7.2 && \
    hmtnote dump

WORKDIR /data/
