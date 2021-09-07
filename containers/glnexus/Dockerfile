################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base-image="clinicalgenomics/mip_base:2.1"
LABEL version="1"
LABEL software="glnexus"
LABEL software.version="1.4.1"
LABEL extra.binaries="glnexus"
LABEL maintainer="Clinical-Genomics/MIP"

WORKDIR /app

ENV PATH="/app:${PATH}"

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
            libc6 \
            ca-certificates && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* && \
    wget --no-verbose https://github.com/dnanexus-rnd/GLnexus/releases/download/v1.4.1/glnexus_cli && \
    chmod 700 /app/glnexus_cli
