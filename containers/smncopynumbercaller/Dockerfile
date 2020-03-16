################## BASE IMAGE ######################

FROM clinicalgenomics/mip:2.0

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:2.0"
LABEL version="1"
LABEL software="SMNCopyNumberCaller"
LABEL software.version="4b2c1ad"
LABEL extra.binaries="smn_caller"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install pip python=3.7 numpy scipy pysam statsmodels

## Clean up after conda
RUN /opt/conda/bin/conda clean -tipsy

## Clone git repository
RUN git clone https://github.com/Illumina/SMNCopyNumberCaller.git /opt/conda/share/SMNCopyNumberCaller

RUN cd /opt/conda/share/SMNCopyNumberCaller && git checkout 4b2c1ad

## Move to smn_caller directory
RUN cd /opt/conda/share/SMNCopyNumberCaller && \
    chmod a+x smn_caller.py

RUN ln --symbolic --force /opt/conda/share/SMNCopyNumberCaller/smn_caller.py /opt/conda/bin/smn_caller.py

RUN ln --symbolic --force /opt/conda/share/SMNCopyNumberCaller/data /opt/conda/bin/data

WORKDIR /data/
