#MIP - Mutation Identification Pipeline

MIP enables identification of potential disease causing variants from sequencing
data. 

![DOI](https://zenodo.org/badge/4091/henrikstranneheim/MIP.png)
=======

##Citing MIP
```
Rapid pulsed whole genome sequencing for comprehensive acute diagnostics of inborn errors of metabolism
Stranneheim H, Engvall M, Naess K, Lesko N, Larsson P, Dahlberg M, Andeer R, Wredenberg A, Freyer C, Barbaro M, Bruhn H, Emahazion T, Magnusson M, Wibom R, Zetterström RH, Wirta V, von Döbeln U, Wedell A.
BMC Genomics. 2014 Dec 11;15(1):1090. doi: 10.1186/1471-2164-15-1090.
PMID:25495354
```

##Overview

MIP performs whole genome or target region analysis of sequenced single-end and/or paired-end
reads from the Illumina plattform in fastq(.gz) format to generate annotated
ranked potential disease causing variants. 

MIP performs QC, alignment, coverage analysis, variant discovery and
annotation, sample checks as well as ranking the found variants according to disease potential
with a minimum of manual intervention. MIP is compatible with Scout for visualization of
identified variants.

##Example Usage
```
perl mip.pl -pMosaikBuild 0 -configFile 1_config.yaml
```

##Features
 - Autonomous
	* Automated install of dependencies via install script
 	* Checks that all dependencies are fulfilled before launching
 	* Builds/downloads references and/or files missing before launching
 	* Splits and merges files for samples and families when relevant
 - Automatic
	* A minimal amount of hands-on time
 	* Tracks and executes all module without manual intervention
 	* Creates internal queues at nodes to optimize processing
 	* Minimal IO between nodes and login node
 - Flexible:
 	* Design your own workflow by turning on/off relevant modules 
 	* Restart an analysis from anywhere in your workflow
 	* Process one, or multiple samples using the module(s) of your choice
 	* Supply parameters on the command line, in a pedigree file or via config files
 	* Simulate your analysis before performing it
 	* Redirect each modules analysis process to a temporary directory (@nodes or @login)
 	* Limit a run to a specific set of genomic intervals
 - Fast
 	* Analyses an exome trio in approximately 6 h
 	* Analyses a genome in approximately 35 h
 	* Rapid mode analyzes a WGS sample in approximately 4 h using a data reduction and parallelization scheme
 - Traceability
 	* Recreate your analysis from the MIP log or generated config files
 	* Logs sample meta-data and sequence meta-data
 	* Logs version numbers of softwares and databases
 	* Checks sample integrity (sex and relationship)
 - Annotation
 	* Gene annotation
 		* Summarise over all transcript and output on gene level
 	* Transcript level annotation
 		* Separate pathogenic transcripts for correct downstream annotation
 	* Annotate all alleles for a position
 		* Split multi-allelic records into single records to ease annotation
 - Standardized
 	* Use standard formats whenever possible
 - Visualization
 	* Output is directly compatibel with Scout

##Getting Started

###Installation

MIP is written in perl and therfore requires that perl is installed on your OS. 

####Prerequisites

#####Programs/Modules
- Perl modules: YAML.pm, Log4perl.pm, List::MoreUtils, DateTime, DateTime::Format::ISO8601, DateTime::Format::HTTP, DateTime::Format::Mail, Set::IntervalTree from CPAN, since these are not included in the perl standard distribution
- Simple Linux Utility for Resource Management (SLURM)
- FastQC
- Mosaik
- BWA
- Sambamba
- SamTools
- BedTools
- PicardTools
- Chanjo
- GATK
- vt
- VEP
- vcfParser.pl (Supplied with MIP)
- SnpEff
- Annovar
- Genmod
- VcfTools
- PLINK

Depending on what programs you include in the MIP analysis you also need to add
these programs to your `bashrc`:

- FastQC
- Mosaik
- BWA
- SamTools
- Tabix
- BedTools
- VcfTools
- PLINK

and these to your python ``virtualenvironment``:

- Chanjo
- GENMOD
- Cosmid (version: 0.4.9.1) for automatic download

#####Meta-Data
- Pedigree file (PLINK-format)
- Configuration file (YAML-format)

###Usage
IP is called from the command line and takes input from the command line
(precedence), a config file (yaml-format) or falls back on defaults where applicable.

Lists are supplied as comma separated input, repeated flag entries on the command line or 
in the config using the yaml format for arrays.
Only flags that will actually be used needs to be specified and MIP will check that all
required parameters are set before submitting to SLURM. 

Program parameters always begins with "p" followed by a capital letter. Program parameters can be set to "0"
(=off), "1" (=on) and "2" (=dry run mode). Any progam can be set to dry
run mode and MIP will create sbatch scripts, but not submit them to SLURM. MIP
can be restarted from any module, but you need to supply previous dependent
programs in dry run mode to ensure proper file handling. 

MIP will overwrite data files when reanalyzing, but keeps all "versioned" sbatch scripts for traceability.

MIP allows individual target file calculations if supplied with a pedigree file or config file
containing the supported capture kits for each sample.

You can always supply ```perl mip.pl -h``` to list all availaible parameters and
defaults.

Example usage:
```
perl mip.pl -f 3 -sampleid 3-1-1A,3-2-1U -sampleid 3-2-2U -pFQC 0 -pMosaikBuild 2 -pMosaikAlign 2 -c 3_config.yaml
```
This will analyse family 3 using 3 individuals from that family and begin the
analysis with programs after MosaikAlign and use all parameter values as
specified in the config file except those supplied on the command line, which
has precedence.

####Input

MIP scripts locations are specified by `-inScriptDir`.

All references and template files should be placed directly in the reference
directory specified by `-referencesDir`, except for ANNOVAR db files, which
should be located in annovar/humandb.

####Output

Analyses done per individual is found under respective sampleID and analyses
done including all samples can be found under the family directory

#####Sbatch Scripts
MIP will create sbatch scripts (.sh) and submit them in proper order with
attached dependencies to SLURM. These sbatch script are placed in the output
script directory specified by `-outScriptDir`. The sbatch scripts are versioned
and will not be overwritten if you begin a new analysis. Versioned "xargs" scripts will also
be created where possible to maximize the use of the cores processing power.

#####Data
MIP will place any generated datafiles in the output data directory specified by
`-outDataDir`. All datatfiles are regenerated for each analysis. STDOUT and
STDERR for each program is written in the program/info directory prior to
alignment and in the aligner/info directory post alignment.
