################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="15"
LABEL software="chromograph"
LABEL software.version="1.3.1"
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

RUN wget --no-verbose https://github.com/Clinical-Genomics/chromograph/archive/refs/tags/v1.3.1.zip && \
    unzip v1.3.1.zip && \
    cd chromograph-1.3.1 && \
    python -m pip install --no-cache-dir .

WORKDIR /data/
