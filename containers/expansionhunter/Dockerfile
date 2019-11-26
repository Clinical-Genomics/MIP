################## BASE IMAGE ######################

FROM clinicalgenomics/mip:latest

################## METADATA ######################

LABEL base_image="clinicalgenomics/mip:latest"
LABEL version="1"
LABEL software="expansionhunter"
LABEL software.version="3.1.2"
LABEL extra.binaries="expansionhunter"

################## MAINTAINER ######################
MAINTAINER Clinical-Genomics/MIP

RUN wget --no-verbose https://github.com/Illumina/ExpansionHunter/releases/download/v3.1.2/ExpansionHunter-v3.1.2-linux_x86_64.tar.gz -O /opt/conda/share/ExpansionHunter-v3.1.2-linux_x86_64.tar.gz && \
    tar --extract --file=/opt/conda/share/ExpansionHunter-v3.1.2-linux_x86_64.tar.gz --directory=/opt/conda/share/ && \
    rm /opt/conda/share/ExpansionHunter-v3.1.2-linux_x86_64.tar.gz

RUN chmod -R 755 /opt/conda/share/ExpansionHunter-v3.1.2-linux_x86_64 && \
    chmod a+x /opt/conda/share/ExpansionHunter-v3.1.2-linux_x86_64/bin/ExpansionHunter && \
    ln --symbolic --force /opt/conda/share/ExpansionHunter-v3.1.2-linux_x86_64/bin/ExpansionHunter /opt/conda/bin/

WORKDIR /data/
