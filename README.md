#MIP - Mutation Identification Pipeline


MIP enables identification of potential disease causing variants from sequencing data. 

##Overview

MIP analyses paired end reads from the Illumina plattform in fastq(.gz) format to annotated 
ranked potential disease causing variants. 

MIP performs QC, aligns reads using Mosaik or BWA, performs variant discovery and 
annotation as well as ranking the found variants according to disease potential.

##Example Usage
```
perl mip.pl -pMosaikBuild 0 -configFile 1_config.yaml
```

##Getting Started

###Installation

MIP is written in perl and therfore requires that perl is installed on your OS. 

####Prerequisites

- Perl YAML.pm module from CPAN (if you want to supply config files to MIP)
- Simple Linux Utility for Resource Management (SLURM)

Dependening on what programs you run you also need to add these programs to your bashrc

- FastQC
- Mosaik
- BWA
- SamTools
- BedTools
- PicardTools
