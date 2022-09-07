################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="13"
LABEL software="chromograph"
LABEL software.version="1.1.5"
LABEL extra.binaries="chromograph"
LABEL maintainer="Clinical-Genomics/MIP"

RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends \
        ca-certificates && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN conda install pip python=3.9 matplotlib && \
    /opt/conda/bin/conda clean -ya

WORKDIR /opt/conda/share

RUN wget --no-verbose https://github.com/mikaell/chromograph/archive/refs/tags/v1.1.5.zip && \
    unzip v1.1.5.zip && \
    cd chromograph-1.1.5 && \
    python -m pip install --no-cache-dir .

WORKDIR /data/
