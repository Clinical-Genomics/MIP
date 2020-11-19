################## BASE IMAGE ######################

FROM perl:5.26

################## METADATA ######################

LABEL base_image="perl:5.26"
LABEL version="1"
LABEL software="MIP"
LABEL software.version="5.26"
LABEL extra.binaries="mip, perl, prove, cpanm"
LABEL maintainer="Clinical-Genomics/MIP"

RUN apt-get update && apt-get install -y --no-install-recommends locales locales-all \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8

# Add MIP into image
ENV MIP_INSTALL_DIR=/workspace/bin
COPY . "$MIP_INSTALL_DIR"
WORKDIR "$MIP_INSTALL_DIR"

# Install MIP dependencies for CPANM modules
RUN cpanm -f Net::DNS
RUN cpanm --installdeps .

# Make executable and add to binary to PATH
RUN chmod a+x "$MIP_INSTALL_DIR/mip"
ENV PATH "$PATH:$MIP_INSTALL_DIR"