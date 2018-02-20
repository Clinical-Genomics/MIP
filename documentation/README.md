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

MIP performs whole genome or target region analysis of sequenced single-end and/or paired-end reads from the Illumina platform in fastq\(.gz\) format to generate annotated ranked potential disease causing variants.

MIP performs QC, alignment, coverage analysis, variant discovery and annotation, sample checks as well as ranking the found variants according to disease potential with a minimum of manual intervention. MIP is compatible with Scout for visualization of identified variants. MIP analyses snv, indels and SV.

MIP has been in use in the clinical production at the Clinical Genomics facility at Science for Life Laboratory since 2014.

## Example Usage

```Bash
$ mip --family_id [family_id] --pbwa_mem 1 --config_file [mip_config.yaml] --pedigree_file [family_id_pedigree.yaml]
```

## Features

* Installation
  * Simple automated install of all programs using conda/SHELL via supplied install script
  * Downloads and prepares references in the installation process
  * Handle conflicting tool dependencies
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
    * Summarize over all transcript and output on gene level
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
* [Perl], version 5.26.0 or above
* [Cpanm](http://search.cpan.org/~miyagawa/App-cpanminus-1.7043/lib/App/cpanminus.pm)
* [Miniconda]

We recommend perlbrew for installing and managing perl and cpanm libraries. Installation instructions and setting up specific cpanm libraries can be found [here](https://github.com/Clinical-Genomics/development/blob/master/perl/installation/installation.md).

#### Automated Installation \(Linux x86\_64\)  
##### 1.Clone the official git repository

```Bash
$ git clone https://github.com/Clinical-Genomics/MIP.git
$ cd MIP
```
##### 2.Install required modules from cpan

```Bash
$ cd definitions
$ cpanm --installdeps .
$ cd -
```  

##### 3.Test conda and mip_installation (optional)

```Bash
$ cd t; prove mip_install.t
$ cd -
```

##### 4.Create the install instructions for MIP  
```Bash
$ perl mip_install.pl
```
This will generate a batch script "mip.sh" for the install in your working directory.

  ###### *Note:*  
  - The batch script will install the MIP dependencies in Conda's root environment. Often it is beneficial to create a separate environment for each of your applications. In order to create a separate environment for MIP supply the ``-env [env_name]`` flag when running *mip_install.pl*.  

  - For a full list of available options and parameters, run: ``$ perl mip_install.pl --help``

  - For a full list of parameter defaults, run: ``$ perl mip_install.pl -ppd``

##### 5.Run the bash script

```Bash
$ bash mip.sh
```
This will install MIP and most of its dependencies into a conda environment.

  ###### *Note:*  
  - Some references are quite large and will take time to download. You might want to run this using screen or tmux.

##### 6.Test your MIP installation (optional)

Make sure to activate your conda environment if that option was used above.  

```Bash
$ cd t; prove -r
$ cd -
```
##### 7.Install tools with conflicting dependencies  
Tools that have conflicting dependencies needs to be installed in separate conda environments. Currently these programs requires separate environments:  

  * Genmod, Chanjo, Multiqc and Variant_integrity
    - requires python 3
  * Peddy
  * CNVnator  
    - Requires access to ROOT which disturbs the normal linking of C libraries  
  * SVDB  
  * VEP


  ```bash
  ## Python 3 tools
  $ perl mip_install.pl -env mip_pyv3.6 --python_version 3.6 --select_program genmod --select_program chanjo --select_program variant_integrity --select_program multiqc
  $ bash mip.sh

  ## Peddy
  $ perl mip_install.pl -env mip_peddy --select_program peddy
  $ bash mip.sh

  ## CNVnator
  $ perl mip_install.pl -env mip_cnvnator --select_program cnvnator
  $ bash mip.sh

  ## SVDB
  $ perl mip_install.pl -env mip_svdb --select_program svdb
  $ bash mip.sh

  ## VEP
  $ perl mip_install.pl -env mip_vep --select_program vep
  $ bash mip.sh
  ```

  In your config yaml file or on the command line you will have to supply the ``module_source_environment_command`` parameter to activate the conda environment specific for the tool. Here is an example with three Python 3 tools in their own environment and Peddy, CNVnator, SVDB and VEP in each own, with some extra initialization:

  ```Yml
  program_source_environment_command:
    pgenmod: "source activate mip_pyv3.6"
  module_source_environment_command:
    pchanjo_sexcheck: "source activate mip_pyv3.6"
    pcnvnator: "LD_LIBRARY_PATH=[CONDA_PATH]/lib/:$LD_LIBRARY_PATH; export LD_LIBRARY_PATH; source [CONDA_PATH]/envs/mip_cnvnator/root/bin/thisroot.sh; source activate mip_cnvnator"
    pmultiqc: "source activate mip_pyv3.6"
    ppeddy: "source activate mip_peddy"
    prankvariant: "source activate mip_pyv3.6"
    psv_rankvariant: "source activate mip_pyv3.6"
    psv_combinevariantcallsets: "source activate mip_svdb"
    psv_varianteffectpredictor: "LD_LIBRARY_PATH=[CONDA_PATH]/envs/mip_vep/lib/:$LD_LIBRARY_PATH; export LD_LIBRARY_PATH; source activate mip_vep"
    pvarianteffectpredictor: "LD_LIBRARY_PATH=[CONDA_PATH]/envs/mip_vep/lib/:$LD_LIBRARY_PATH; export LD_LIBRARY_PATH; source activate mip_vep"
    pvariant_integrity: "source activate mip_pyv3.6"
  source_main_environment_commands:
    - source
    - activate
    - mip
  ```

  MIP will execute this on the node before executing the program and then revert to the ``--source_main_environment_command`` if set. Otherwise ``source deactivate`` is used to return to the conda root environment.

### Usage

MIP is called from the command line and takes input from the command line \(precedence\) or falls back on defaults where applicable.

Lists are supplied as repeated flag entries on the command line or in the config using the yaml format for arrays.  
Only flags that will actually be used needs to be specified and MIP will check that all required parameters are set before submitting to SLURM.

Program parameters always begins with "p" followed by a capital letter. Program parameters can be set to "0" \(=off\), "1" \(=on\) and "2" \(=dry run mode\). Any progam can be set to dry run mode and MIP will create sbatch scripts, but not submit them to SLURM. MIP can be restarted from any module, but you need to supply previous dependent programs in dry run mode to ensure proper file handling.

MIP will overwrite data files when reanalyzing, but keeps all "versioned" sbatch scripts for traceability.

You can always supply `perl mip.pl --help` to list all available parameters and defaults.

Example usage:
```Bash
$ mip -f 3 --sample_ids 3-1-1A --sample_ids 3-2-1U --sample_ids 3-2-2U -pfqc 0 --pbwa_mem 2 -c 3_config.yaml
```

This will analyse family 3 using 3 individuals from that family and begin the analysis with programs after Bwa mem and use all parameter values as specified in the config file except those supplied on the command line, which has precedence.

#### Input

All references and template files should be placed directly in the reference directory specified by `--reference_dir`.

##### Meta-Data

* [Configuration file] \(YAML-format\)
* [Gene panel file]
* [Pedigree file] \(YAML-format\)
* [Rank model file] \(Ini-format; Snv/indel\)
* [SV rank model file] \(Ini-format; SV\)
* [Qc regexp file] \(YAML-format\)

#### Output

Analyses done per individual is found in each sample_id directory and analyses done including all samples can be found in the family directory.

##### Sbatch Scripts

MIP will create sbatch scripts \(.sh\) and submit them in proper order with attached dependencies to SLURM. These sbatch script are placed in the output script directory specified by `--outscript_dir`. The sbatch scripts are versioned and will not be overwritten if you begin a new analysis. Versioned "xargs" scripts will also be created where possible to maximize the use of the cores processing power.

##### Data

MIP will place any generated datafiles in the output data directory specified by `--outdata_dir`. All data files are regenerated for each analysis. STDOUT and STDERR for each program is written in the program/info directory prior to alignment and in the aligner/info directory post alignment.

[Configuration file]: https://github.com/Clinical-Genomics/MIP/blob/master/templates/mip_config.yaml
[CPAN]: https://www.cpan.org/
[GATK]:https://software.broadinstitute.org/gatk/
[Gene panel file]: https://github.com/Clinical-Genomics/MIP/blob/master/templates/aggregated_master.txt
[Miniconda]: http://conda.pydata.org/miniconda.html
[Pedigree file]: https://github.com/Clinical-Genomics/MIP/tree/master/templates/643594-miptest_pedigree.yaml
[Perl]:https://www.perl.org/
[Rank model file]: https://github.com/Clinical-Genomics/MIP/blob/master/templates/rank_model_cmms_-v1.20-.ini
[SV rank model file]: https://github.com/Clinical-Genomics/MIP/blob/master/templates/svrank_model_cmms_-v1.2-.ini
[Qc regexp file]: https://github.com/Clinical-Genomics/MIP/blob/master/templates/qc_regexp_-v1.17-.yaml
