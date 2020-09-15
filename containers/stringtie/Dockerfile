################## BASE IMAGE ######################

FROM clinicalgenomics/mip:2.0

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:2.0"
LABEL version="1"
LABEL software="stringtie"
LABEL software.version="2.1.3b"
LABEL extra.binaries="stringtie"
LABEL maintainer="Clinical-Genomics/MIP"

WORKDIR /app

RUN wget --no-verbose http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.1.3b.Linux_x86_64.tar.gz && \
    tar -xzf ./stringtie-2.1.3b.Linux_x86_64.tar.gz && \
    rm ./stringtie-2.1.3b.Linux_x86_64.tar.gz

RUN ln -s /app/stringtie-2.1.3b.Linux_x86_64/stringtie /opt/conda/bin/stringtie
