################## BASE IMAGE ######################

FROM ubuntu:bionic

################## METADATA ######################

LABEL base_image="ubuntu:bionic"
LABEL version="2"
LABEL software="bwa-mem2"
LABEL software.version="2.2"
LABEL extra.binaries="bwa-mem2"
LABEL maintainer="Clinical-Genomics/MIP"

## Install curl
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    libbz2-dev \
    libreadline-dev && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /app

## Install bwa-mem2
RUN curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2/bwa-mem2-2.2_x64-linux.tar.bz2 \
  | tar jxf -

ENV PATH=/app/bwa-mem2-2.2_x64-linux:${PATH}