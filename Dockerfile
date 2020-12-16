################## BASE IMAGE ######################

FROM localhost/slurm-container-cluster

################## METADATA ######################

LABEL base_image="localhost/slurm-container-cluster"
LABEL version="1"
LABEL software="MIP"
LABEL software.version="5.26"
LABEL extra.binaries="mip, perl, prove, cpanm, carton"
LABEL maintainer="Clinical-Genomics/MIP"

RUN dnf -y install perl-App-cpanminus perl-CPAN perl-devel bcftools \
  && dnf clean all \
  && rm -rf /var/cache/yum

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

ENV PERL5LIB="$MIP_INSTALL_DIR"/local/lib/perl5