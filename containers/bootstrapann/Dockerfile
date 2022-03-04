################## BASE IMAGE ######################

FROM continuumio/miniconda3:4.7.12

################## METADATA ######################

LABEL base_image="continuumio/miniconda3:4.7.12"
LABEL version="2"
LABEL software="BootstrapAnn"
LABEL software.version="df02f35"
LABEL extra.binaries="BootstrapAnn.py"
LABEL maintainer="Clinical-Genomics/MIP"

RUN apt-get --allow-releaseinfo-change update --fix-missing && \
    apt-get install -y --no-install-recommends procps && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN conda install pip python=2.7 numpy scipy && \
    /opt/conda/bin/conda clean -ya

## Clone git repository
RUN git clone https://github.com/J35P312/BootstrapAnn.git /opt/conda/share/BootstrapAnn

WORKDIR /opt/conda/share/BootstrapAnn

RUN git checkout df02f35

RUN chmod a+x BootstrapAnn.py

RUN ln --symbolic --force /opt/conda/share/BootstrapAnn/BootstrapAnn.py /opt/conda/bin/BootstrapAnn.py

WORKDIR /data/
