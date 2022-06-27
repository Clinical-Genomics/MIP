################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="2"
LABEL software="stringtie"
LABEL software.version="2.2.2"
LABEL extra.binaries="stringtie"
LABEL maintainer="Clinical-Genomics/MIP"

WORKDIR /app

RUN wget --no-verbose http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.2.1.Linux_x86_64.tar.gz && \
    tar -xzf ./stringtie-2.2.1.Linux_x86_64.tar.gz && \
    rm ./stringtie-2.2.1.Linux_x86_64.tar.gz && \
    ln -s /app/stringtie-2.2.1.Linux_x86_64/stringtie /opt/conda/bin/stringtie
