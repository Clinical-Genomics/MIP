################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base-image="clinicalgenomics/mip_base:2.1"
LABEL version="3"
LABEL software="Cyrius"
LABEL software.version="v1.1.1"
LABEL extra.binaries="star_caller.py"
LABEL maintainer="Clinical-Genomics/MIP"

WORKDIR /opt/conda/share

## Pysam errors when pip installing
RUN conda install pysam=0.16.0.1 && \
    conda clean -ya

## Download and extract
RUN wget --no-verbose --no-check-certificate https://github.com/Illumina/Cyrius/archive/v1.1.1.zip && \
    unzip v1.1.1.zip && \
    rm v1.1.1.zip

## Move to directory and install requirements
RUN cd Cyrius-1.1.1 && \
    python -m pip install --no-cache-dir -r requirements.txt

RUN chmod a+x /opt/conda/share/Cyrius-1.1.1/star_caller.py

ENV PATH ${PATH}:/opt/conda/share/Cyrius-1.1.1
