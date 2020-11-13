################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="1"
LABEL software="BootstrapAnn"
LABEL software.version="df02f35"
LABEL extra.binaries="BootstrapAnn.py"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install pip python=2.7 numpy scipy

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

## Clone git repository
RUN git clone https://github.com/J35P312/BootstrapAnn.git /opt/conda/share/BootstrapAnn

RUN cd /opt/conda/share/BootstrapAnn && git checkout df02f35

## Move to smn_caller directory
RUN cd /opt/conda/share/BootstrapAnn && \
    chmod a+x BootstrapAnn.py

RUN ln --symbolic --force /opt/conda/share/BootstrapAnn/BootstrapAnn.py /opt/conda/bin/BootstrapAnn.py

WORKDIR /data/
