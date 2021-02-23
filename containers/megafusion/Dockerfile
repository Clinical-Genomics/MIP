################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="1"
LABEL software="MegaFusion"
LABEL software.version="5828fbe"
LABEL extra.binaries="MegaFusion.py"
LABEL maintainer="Clinical-Genomics/MIP"

RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends ca-certificates && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN conda install python=3.7 && \
    /opt/conda/bin/conda clean -ya

## Clone git repository
RUN git clone https://github.com/J35P312/MegaFusion.git /opt/conda/share/MegaFusion

WORKDIR /opt/conda/share/MegaFusion

## Make sure we're on the right commit
RUN git reset --hard 5828fbe
