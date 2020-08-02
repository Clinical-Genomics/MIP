################## BASE IMAGE ######################

FROM clinicalgenomics/mip:2.0

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:2.0"
LABEL version="1"
LABEL software="Cyrius"
LABEL software.version="1.0"
LABEL extra.binaries="star_caller.py"
LABEL maintainer="Clinical-Genomics/MIP"

WORKDIR /opt/conda/share

## Pysam errors when pip installing
RUN conda install pysam=0.16.0 && \
    conda clean -ya

## Download and extract
RUN wget --no-verbose https://github.com/Illumina/Cyrius/archive/v1.0.zip && \
    unzip v1.0.zip && \
    rm v1.0.zip

## Move to directory and install requirements
RUN cd Cyrius-1.0 && \
   python -m pip install --no-cache-dir -r requirements.txt

RUN chmod a+x /opt/conda/share/Cyrius-1.0/star_caller.py 

ENV PATH ${PATH}:/opt/conda/share/Cyrius-1.0
