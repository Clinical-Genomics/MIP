################## BASE IMAGE ######################

FROM clinicalgenomics/mip:2.0

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:2.0"
LABEL version="1"
LABEL software="CADD"
LABEL software.version="1.6"
LABEL extra.binaries="CADD.sh"
LABEL maintainer="Clinical-Genomics/MIP"

## Get snakemake
RUN conda install -c conda-forge -c bioconda snakemake=5.17.0

## Download CADD
RUN wget --no-verbose https://github.com/kircherlab/CADD-scripts/archive/CADD1.6.zip -O /opt/conda/share/CADD1.6.zip

## Extract
RUN unzip -o /opt/conda/share/CADD1.6.zip -d /opt/conda/share/ && \
    rm /opt/conda/share/CADD1.6.zip

## Install 
RUN cd /opt/conda/share/CADD-scripts-CADD1.6 && \
    snakemake test/input.tsv.gz --use-conda --conda-create-envs-only \
    --conda-prefix /opt/conda/share/CADD-scripts-CADD1.6/envs \
    --cores 1 --configfile config/config_GRCh38_v1.6.yml --snakefile Snakefile

## Link CADD.sh
RUN  ln -s /opt/conda/share/CADD-scripts-CADD1.6/CADD.sh /opt/conda/bin/CADD.sh

## Clean up after conda
RUN /opt/conda/bin/conda clean -ya

WORKDIR /data/
