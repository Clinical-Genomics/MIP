################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ########################


LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="1"
LABEL software="telomerecat"
LABEL software.version="4.0.0"
LABEL extra.binaries="telomerecat"
LABEL maintainer="Clinical-Genomics/MIP"

RUN conda install python=3.8 && \
    conda clean -ya

RUN pip3 install https://github.com/cancerit/parabam/releases/download/3.0.0/parabam-3.0.0.tar.gz
RUN pip3 install https://github.com/cancerit/telomerecat/archive/4.0.0.tar.gz

RUN pip cache purge