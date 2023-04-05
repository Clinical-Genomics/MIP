################## BASE IMAGE ######################

FROM clinicalgenomics/mip_base:2.1

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip_base:2.1"
LABEL version="1"
LABEL software="retroseq"
LABEL software.version="1.5_9d4f3b5"
LABEL extra.binaries="retroseq"
LABEL maintainer="Clinical-Genomics/MIP"

RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends \
        ca-certificates && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN conda install samtools exonerate bedtools bcftools && \
    /opt/conda/bin/conda clean -ya

RUN git clone https://github.com/tk2/RetroSeq.git /opt/conda/share/RetroSeq

WORKDIR /opt/conda/share/RetroSeq

## Remove samtool check
## Make sure we're on the right commit and remove samtools check
RUN git reset --hard 9d4f3b5 && \
    sed -i '/RetroSeq::Utilities::checkBinary( q\[samtools\].*/d' ./bin/retroseq.pl
