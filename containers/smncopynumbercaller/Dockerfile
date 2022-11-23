################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="4"
LABEL software="SMNCopyNumberCaller"
LABEL software.version="v1.1.2"
LABEL extra.binaries="smn_caller"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install pip python=3.7 numpy scipy pysam statsmodels && \
    /opt/conda/bin/conda clean -ya && \
    wget --no-verbose --no-check-certificate https://github.com/Illumina/SMNCopyNumberCaller/archive/v1.1.2.zip -O /opt/conda/share/smncopynumbercaller-v.1.1.2.zip && \
    unzip -o /opt/conda/share/smncopynumbercaller-v.1.1.2.zip -d /opt/conda/share/

WORKDIR /opt/conda/share/SMNCopyNumberCaller-1.1.2

RUN chmod a+x smn_caller.py && \
    ln --symbolic --force /opt/conda/share/SMNCopyNumberCaller-1.1.2/smn_caller.py /opt/conda/bin/smn_caller.py && \
    ln --symbolic --force /opt/conda/share/SMNCopyNumberCaller-1.1.2/data /opt/conda/bin/data

WORKDIR /data/
