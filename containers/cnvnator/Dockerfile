################## BASE IMAGE ######################

FROM ubuntu:bionic

################## METADATA ######################

LABEL base_image="ubuntu:bionic"
LABEL version="1"
LABEL software="CNVnator"
LABEL software.version="0.4.1"
LABEL extra.binaries="cnvnator, cnvnator2VCF.pl, plotbaf.py, plotcircular.py, plotrdbaf.py, pytools.py"
LABEL maintainer="Clinical-Genomics/MIP"

## Install wget and unzip
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y \
    build-essential \
    cmake \
    language-pack-en-base \
    libbz2-dev \
    libcrypto++-dev \
    libcurl4-gnutls-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libreadline-dev \
    libssl-dev \
    libx11-xcb-dev \
    libxft-dev \
    unzip \
    wget \
    zlib1g-dev && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    update-locale LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8
    
WORKDIR /app

## Install ROOT
RUN wget --no-verbose https://root.cern.ch/download/root_v6.20.06.Linux-ubuntu18-x86_64-gcc7.5.tar.gz && \
    tar -xvf root_v6.20.06.Linux-ubuntu18-x86_64-gcc7.5.tar.gz && \
    rm root_v6.20.06.Linux-ubuntu18-x86_64-gcc7.5.tar.gz && \
    rm -rf root/tutorials

## Set paths for ROOT
ENV PATH=/usr/local/bin:${PATH} 
ENV ROOTSYS=/app/root 
ENV PATH=${PATH}:${ROOTSYS}/bin 
ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ROOTSYS}/lib

## Install samtools 
RUN wget --no-verbose https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 && \
    bunzip2 samtools-1.10.tar.bz2 && \
    tar -xf samtools-1.10.tar && \
    cd samtools-1.10 && \
    ./configure && \
    make all all-htslib && \
    make install install-htslib && \
    ln --symbolic --force /app/samtools-1.10/samtools /usr/local/bin && \
    rm /app/samtools-1.10.tar

## Install CNVnator
RUN wget --no-verbose https://github.com/abyzovlab/CNVnator/releases/download/v0.4.1/CNVnator_v0.4.1.zip && \
    unzip CNVnator_v0.4.1.zip && \
    cd CNVnator_v0.4.1/src/ && \
    ln -s /app/samtools-1.10 /app/CNVnator_v0.4.1/src/samtools && \
    make LIBS="-lcrypto" && \
    ln --symbolic --force /app/CNVnator_v0.4.1/src/cnvnator /usr/local/bin  && \
    ln --symbolic --force /app/CNVnator_v0.4.1/src/cnvnator2VCF.pl /usr/local/bin && \
    ln --symbolic --force /app/CNVnator_v0.4.1/src/plotbaf.py /usr/local/bin && \
    ln --symbolic --force /app/CNVnator_v0.4.1/src/plotcircular.py /usr/local/bin && \
    ln --symbolic --force /app/CNVnator_v0.4.1/src/plotrdbaf.py /usr/local/bin && \
    ln --symbolic --force /app/CNVnator_v0.4.1/src/plotrdbaf.py /usr/local/bin && \
    ln --symbolic --force /app/CNVnator_v0.4.1/src/pytools.py /usr/local/bin && \
    rm /app/CNVnator_v0.4.1.zip
