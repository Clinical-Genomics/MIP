#MIP - Mutation Identification Pipeline


MIP enables identification of potential disease causing variants from sequencing data. 

##Overview

MIP performs whole genome or target region analysis of sequenced paired end reads from the Illumina plattform in 
fastq(.gz) format to generate annotated ranked potential disease causing variants. 

MIP performs QC, aligns reads using Mosaik or BWA, variant discovery and annotation as well as ranking the found 
variants according to disease potential with a minimum of manual intervention.

##Example Usage
```
perl mip.pl -pMosaikBuild 0 -configFile 1_config.yaml
```

##Getting Started

###Installation

MIP is written in perl and therfore requires that perl is installed on your OS. 

####Prerequisites

#####Programs/Modules
- Perl YAML.pm module from CPAN, since this is not included in the perl standard distribution (if you want to 
  supply config files to MIP)
- Simple Linux Utility for Resource Management (SLURM)
- FastQC
- Mosaik
- BWA
- SamTools
- BedTools
- PicardTools
- GATK
- ANNOVAR
- intersectCollect.pl
- add_depth.pl
- rank_list_filter.pl
- VcfTools
- PLINK

Depending on what programs you include in the MIP analysis you also need to add these programs to your `bashrc`:

- FastQC
- Mosaik
- BWA
- SamTools
- BedTools
- PicardTools
- VcfTools
- PLINK

#####Meta-Data
- Pedigree file (PLINK-format)
- Template files for intersectCollect.pl

#####Infiles
Needs to be on the format: `lane_date_flow-cellID_sampleID_indexX_direction.fastq`. For instance:

```
6_120313_AC0GTUACXX_9-1-1A_indexACTTGA_1.fastq.gz
6_120313_AC0GTUACXX_9-1-1A_indexACTTGA_2.fastq.gz
```

###Usage
MIP is called from the command line and takes input from the command line (precedence) or a config file (yaml-format).
Lists are supplied as comma separated input or repeated flag entries. Only flags that will actually be used needs to 
be specified and MIP will check that all required parameters are set before submitting to SLURM. 

Program parameters always begins with "p". Program parameters can be set to "0" (=do not run), "1" (=run) and "2" 
(=dry run mode). Any progam can be set to dry run mode and MIP will create sbatch scripts, but not submitted them to 
SLURM. MIP can be restarted from any module, but you need to supply previous dependent programs in dry run mode to 
ensure proper file handling. 

MIP allows individual target file calculations if supplied with a pedigree file containing the supported capture kits.

You can always supply ```perl mip.pl -h``` to list all availaible parameters and defaults.  

Example usage:
```
perl mip.pl -f 3 -sampleid 3-1-1A,3-2-1U -sampleid 3-2-2U -pFQC 0 -pMosaikBuild 2 -pMosaikAlign 2 -c 3_config.yaml
```
This will analyse family 3 using 3 individuals from that family and begin the analysis with programs after 
MosaikAlign and use all parameter values as specified in the config file except those supplied on the command line, 
which has precedence.

####Input

All MIP scripts (including mip.pl) should be placed in the script directory specified by `-inScriptDir`.

All references and template files should be placed directly in the reference directory specified by `-referencesDir`,
except for ANNOVAR db files, which should be located in annovar/humandb.

####Output

Analyses done per individual is found under respective sampleID and analyses done including all samples can be found
under the family directory

#####Sbatch Scripts
MIP will create sbatch scripts (.sh) and submit them in proper order with attached dependencies to SLURM. These sbatch 
script are placed in the output script directory specified by `-outScriptDir`. The sbatch scripts are versioned and will
not be overwritten if you begin a new analysis.

#####Data
MIP will place any generated datafiles in the output data directory specified by `-outDataDir`. All datatfiles are 
regenerated for each analysis. STDOUT and STDERR for each program is written in the program/info directory prior to 
alignment and in the aligner/info directory post alignment.
