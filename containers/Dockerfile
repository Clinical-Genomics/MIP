# Base image
FROM ubuntu:16.04

################## METADATA ######################

LABEL base_image="ubuntu:16.04"
LABEL version="1"
LABEL software="mip"
LABEL software.version="7.1.4"
LABEL about.summary="Base image for MIP's dependency tools"
LABEL about.home="https://github.com/Clinical-Genomics/MIP"
LABEL about.documentation="https://clinical-genomics.gitbook.io/project-mip/"
LABEL about.license_file="https://github.com/Clinical-Genomics/MIP/blob/master/LICENSE"
LABEL about.license="MIT License (MIT)"
LABEL about.tags="Clinical,variants,analysis,pipeline"

################## MAINTAINER ######################
MAINTAINER Henrik Stranneheim <henrik.stranneheim@scilifelab.se>
MAINTAINER Anders Jemt <anders.jemt@scilifelab.se>

ENV DEBIAN_FRONTEND noninteractive

RUN mv /etc/apt/sources.list /etc/apt/sources.list.bkp && \
    bash -c 'echo -e "deb mirror://mirrors.ubuntu.com/mirrors.txt xenial main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-updates main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-backports main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-security main restricted universe multiverse\n\n" > /etc/apt/sources.list' && \
    cat /etc/apt/sources.list.bkp >> /etc/apt/sources.list && \
    cat /etc/apt/sources.list

RUN apt-get update --fix-missing && \
    apt-get install -y \
        build-essential \
        coreutils \
        dialog \
        git \
        language-pack-en-base \
        libbz2-dev \
        libcurl4-openssl-dev \
        liblzma-dev \
        libncurses5-dev \
        libncursesw5-dev \
        libreadline-dev \
        libssl-dev \
        unzip \
        wget \
        zlib1g-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install conda and give write permissions to conda folder
RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda3-4.7.10-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

RUN mkdir /data /config

ENV PATH=$PATH:/opt/conda/bin

RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

# Clean up after conda
RUN  /opt/conda/bin/conda clean -tipsy

VOLUME ["/data", "/config"]

WORKDIR /data
