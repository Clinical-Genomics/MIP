################## BASE IMAGE ######################

FROM perl:5.26

################## METADATA ######################

LABEL base_image="perl:5.26"
LABEL version="1"
LABEL software="MIP"
LABEL software.version="5.26"
LABEL extra.binaries="mip, perl, prove, cpanm, carton"
LABEL maintainer="Clinical-Genomics/MIP"

RUN apt-get update && apt-get install -y --no-install-recommends locales locales-all \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8

RUN cpanm install Carton

# Add MIP into image
ENV MIP_INSTALL_DIR=/workspace/bin
COPY . "$MIP_INSTALL_DIR"
WORKDIR "$MIP_INSTALL_DIR"
# Make executable and add to binary to PATH
RUN chmod a+x "$MIP_INSTALL_DIR/mip"
ENV PATH "$PATH:$MIP_INSTALL_DIR"

# Remove any outside local from previous carton installs
RUN rm -rf local

RUN carton install

ENTRYPOINT ["carton", "exec"]
