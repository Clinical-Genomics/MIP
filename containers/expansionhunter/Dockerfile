################## BASE IMAGE ######################

FROM ubuntu:bionic

################## METADATA ######################

LABEL base_image="ubuntu:bionic"
LABEL version="4"
LABEL software="ExpanionHunter"
LABEL software.version="5.0.0"
LABEL extra.binaries="ExpanionHunter"
LABEL maintainer="Clinical-Genomics/MIP"

## Install wget
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    wget \
    libreadline-dev && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /app

RUN wget -nv https://github.com/Illumina/ExpansionHunter/releases/download/v5.0.0/ExpansionHunter-v5.0.0-linux_x86_64.tar.gz && \
    tar -xvf ExpansionHunter-v5.0.0-linux_x86_64.tar.gz && \
    rm ExpansionHunter-v5.0.0-linux_x86_64.tar.gz

ENV PATH=/app/ExpansionHunter-v5.0.0-linux_x86_64/bin:${PATH}
