# MIP - Mutation Identification Pipeline

[![Build Status](https://travis-ci.org/henrikstranneheim/MIP.svg?branch=develop)](https://travis-ci.org/henrikstranneheim/MIP)

MIP enables identification of potential disease causing variants from sequencing data.

# [![DOI](https://zenodo.org/badge/7667877.svg)](https://zenodo.org/badge/latestdoi/7667877)

## Citing MIP

```
Rapid pulsed whole genome sequencing for comprehensive acute diagnostics of inborn errors of metabolism
Stranneheim H, Engvall M, Naess K, Lesko N, Larsson P, Dahlberg M, Andeer R, Wredenberg A, Freyer C, Barbaro M, Bruhn H, Emahazion T, Magnusson M, Wibom R, Zetterström RH, Wirta V, von Döbeln U, Wedell A.
BMC Genomics. 2014 Dec 11;15(1):1090. doi: 10.1186/1471-2164-15-1090.
PMID:25495354
```

## Overview

MIP performs whole genome or target region analysis of sequenced single-end and/or paired-end reads from the Illumina plattform in fastq\(.gz\) format to generate annotated ranked potential disease causing variants.

MIP performs QC, alignment, coverage analysis, variant discovery and annotation, sample checks as well as ranking the found variants according to disease potential with a minimum of manual intervention. MIP is compatible with Scout for visualization of identified variants. MIP analyses snv, indels and SV.

MIP has been in use in the clinical production at the Clinical Genomics facility at Science for Life Laboratory since 2014.

## Example Usage

```
$ mip --family_id [family_id] --pbwa_mem 1 --config_file [mip_config.yaml] --pedigree_file [family_id_pedigree.yaml]
```

## Features

* Installation
  * Simple automated install of all programs using conda/SHELL via supplied install script
  * Downloads and prepares references in the installation process
* Autonomous
  * Checks that all dependencies are fulfilled before launching
  * Builds and prepares references and/or files missing before launching
  * Decompose and normalise reference\(s\) and variant vcf\(s\)
  * Splits and merges files/contigs for samples and families when relevant
* Automatic
  * A minimal amount of hands-on time
  * Tracks and executes all module without manual intervention
  * Creates internal queues at nodes to optimize processing
  * Minimal IO between nodes and login node
* Flexible:
  * Design your own workflow by turning on/off relevant modules
  * Restart an analysis from anywhere in your workflow
  * Process one, or multiple samples using the module\(s\) of your choice
  * Supply parameters on the command line, in a pedigree.yaml file or via config files
  * Simulate your analysis before performing it
  * Redirect each modules analysis process to a temporary directory \(@nodes or @login\)
  * Limit a run to a specific set of genomic intervals
  * Use multiple variant callers for both snv, indels and SV
  * Use multiple annotation programs
  * Optionally split data into clinical variants and research variants
* Fast
  * Analyses an exome trio in approximately 4 h
  * Analyses a genome in approximately 21 h
  * Rapid mode analyzes a WGS sample in approximately 4 h using a data reduction and parallelization scheme
* Traceability
  * Track the status of each modules through dynamically updated status logs
  * Recreate your analysis from the MIP log or generated config files
  * Log sample meta-data and sequence meta-data
  * Log version numbers of softwares and databases
  * Checks sample integrity \(sex, contamination, duplications, ancestry, inbreeding and relationship\)
  * Test data output existens and integrity using automated tests
* Annotation
  * Gene annotation
    * Summarise over all transcript and output on gene level
  * Transcript level annotation
    * Separate pathogenic transcripts for correct downstream annotation
  * Annotate all alleles for a position
    * Split multi-allelic records into single records to ease annotation
    * Left align and trim variants to normalise them prior to annotation
  * Extracts QC-metrics and stores them in YAML format
  * Annotate coverage across genetic regions via Sambamba and Chanjo
* Standardized
  * Use standard formats whenever possible
* Visualization
  * Ranks variants according to pathogenic potential
  * Output is directly compatible with [Scout](https://github.com/Clinical-Genomics/scout)

## Getting Started

### Installation

MIP is written in perl and therefore requires that perl is installed on your OS.

#### Prerequisites
* [Perl], version 5.18.2 or above
* [Cpanm](http://search.cpan.org/~miyagawa/App-cpanminus-1.7043/lib/App/cpanminus.pm)
* [Miniconda]

#### Automated Installation \(Linux x86\_64\)
1. Clone the official git repository
```Bash
$ git clone https://github.com/Clinical-Genomics/MIP.git
$ cd MIP
```
2. Install required modules from cpan
```Bash
$ cd definitions
$ cpanm --installdeps .
$ cd -
```
3. Test conda and mip_install (optional)
```Bash
$ cd t; prove mip_install.t
$ cd -
```

2. Create the install instructions for MIP
```
$ perl mip_install.pl
```
This will generate a batch script "mip.sh" for the install in your working directory.

  ###### *Note:*  
  The batch script will install the MIP dependencies in Conda's root environment. Often it is beneficial to create a separate environment for each of your applications. In order to create a separate environment for MIP supply the ``-env [env_name]`` flag when running *mip_install.pl*.  

  There is currrently an issue with installing GATK via Conda causing the installation process to freeze. The solution at the moment is to turn of the GATK installation by running the install script together with the ``-skip gatk`` flag and install [GATK] manually.

  For a full list of available options and parameters, run: ``$ perl mip_install.pl --help``.  
  For a full list of parameter defaults, run: ``$ perl mip_install.pl -ppd``

3. Run the bash script
```Bash
$ bash mip.sh
```
This will install all the dependencies of MIP and other modules included in MIP into a conda environment.

  ###### *Note:*  
  Some references are quite large and will take time to download. You might want to run this using screen or tmux.

4. Test your MIP installation (optional)  
  Make sure to activate your conda environment if that option was used above.  

  ```Bash
  $ cd t; prove run_tests.t
  $ cd -
  ```

5. Installing tools with conflicting dependencies
Create seperate environment for each tool that has to have a seperate dependecy demands. For instance a tools that require python v3 or higher and are not compatible with older versions of pyton. 

```conda create --name python_v3.6_tools --quiet --yes python=3.6 pip
$ source activate python_v3.6_tools
$ pip install genmod==3.7.2 chanjo==4.2.0
```

In your ~/.bashrc add a variable and and source of a file containing aliases of the tools you want to activate:
```# source aliases
CONDA_BIN="[MINICONDA_PATH]/envs/python_v3.6_tools/bin"
source ~/dotfiles/aliases.sh 
```

Add the aliases of the binaries in the aliases.sh file:
alias chanjo="${CONDA_BIN}/chanjo"
alias genmod="${CONDA_BIN}/genmod"

### Usage

MIP is called from the command line and takes input from the command line \(precedence\ or falls back on defaults where applicable.

Lists are supplied as repeated flag entries on the command line or in the config using the yaml format for arrays.  
Only flags that will actually be used needs to be specified and MIP will check that all required parameters are set before submitting to SLURM.

Program parameters always begins with "p" followed by a capital letter. Program parameters can be set to "0" \(=off\), "1" \(=on\) and "2" \(=dry run mode\). Any progam can be set to dry run mode and MIP will create sbatch scripts, but not submit them to SLURM. MIP can be restarted from any module, but you need to supply previous dependent programs in dry run mode to ensure proper file handling.

MIP will overwrite data files when reanalyzing, but keeps all "versioned" sbatch scripts for traceability.

You can always supply `perl mip.pl -h` to list all availaible parameters and defaults.

Example usage:
```
$ mip -f 3 -sampleid 3-1-1A,3-2-1U -sampleid 3-2-2U -pfqc 0 --pbwa_mem 2 -c 3_config.yaml
```

This will analyse family 3 using 3 individuals from that family and begin the analysis with programs after Bwa mem and use all parameter values as specified in the config file except those supplied on the command line, which has precedence.

#### Input

All references and template files should be placed directly in the reference directory specified by `--reference_dir`, except for ANNOVAR db files, which should be located in annovar/humandb.

##### Meta-Data

* [Pedigree file] \(YAML-format\)
* [Configuration file] \(YAML-format\)


#### Output

Analyses done per individual is found in each sample_id directory and analyses done including all samples can be found in the family directory.

##### Sbatch Scripts

MIP will create sbatch scripts \(.sh\) and submit them in proper order with attached dependencies to SLURM. These sbatch script are placed in the output script directory specified by `--outscript_dir`. The sbatch scripts are versioned and will not be overwritten if you begin a new analysis. Versioned "xargs" scripts will also be created where possible to maximize the use of the cores processing power.

##### Data

MIP will place any generated datafiles in the output data directory specified by `--outdata_dir`. All data files are regenerated for each analysis. STDOUT and STDERR for each program is written in the program/info directory prior to alignment and in the aligner/info directory post alignment.

[Miniconda]: http://conda.pydata.org/miniconda.html
[CPAN]: https://www.cpan.org/
[Pedigree file]: https://github.com/Clinical-GenomicsMIP/tree/master/templates/643594-miptest_pedigree.yaml
[Configuration file]: https://github.com/Clinical-Genomics/MIP/blob/master/templates/mip_config.yaml
[GATK]:https://software.broadinstitute.org/gatk/
[Perl]:https://www.perl.org/
