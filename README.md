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

#####Programs/Modules
- Perl YAML.pm module from CPAN since this is not included in the perl standard distribution
(if you want to supply config files to MIP)
- Simple Linux Utility for Resource Management (SLURM)
- FastQC
- Mosaik
- BWA
- SamTools
- BedTools
- PicardTools
- GATK
- ANNOVAR

Dependening on what programs you include in the MIP analysis you also need to add these programs to your bashrc

- FastQC
- Mosaik
- BWA
- SamTools
- BedTools
- PicardTools

#####Meta-Data
- Pedigree file (PLINK-format)

#####Infiles
Needs to be on the format: "lane_date_flow-cellID_sampleID_indexX_direction.fastq".  For instance:

```
6_120313_AC0GTUACXX_9-1-1A_indexACTTGA_1.fastq.gz
6_120313_AC0GTUACXX_9-1-1A_indexACTTGA_2.fastq.gz
```

###Usage
MIP is called from the command line and takes input from the command line (precedence) or a config file (yaml-format).
Lists are supplied as comma separated input or repeated flag entries. Only flags that will actually be used needs to 
be specified and MIP will check that all required parameters are set before submitting to SLURM. 
Program parameters always begins with "p". 
You can always supply ```perl mip.pl -h``` to list all availaible parameters and defaults.  

