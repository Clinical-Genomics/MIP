################## BASE IMAGE ######################

FROM perl:5.26

################## METADATA ######################

LABEL base_image="perl:5.26"
LABEL version="1"
LABEL software="perl"
LABEL software.version="5.26"
LABEL extra.binaries="perl"
LABEL maintainer="Clinical-Genomics/MIP"

RUN apt-get update && apt-get install -y --no-install-recommends locales locales-all \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8

WORKDIR /data/
