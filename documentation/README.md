# MIP - Mutation Identification Pipeline

![MIP CI conda production install](https://github.com/Clinical-Genomics/MIP/workflows/MIP%20CI%20conda%20production%20install/badge.svg?branch=master)
[![Coverage Status](https://coveralls.io/repos/github/Clinical-Genomics/MIP/badge.svg?branch=master)](https://coveralls.io/github/Clinical-Genomics/MIP?branch=master)
[![GitHub license](https://img.shields.io/badge/License-MIT-blue.svg)](https://raw.githubusercontent.com/Clinical-Genomics/MIP/master/LICENSE)
[![GitHub Releases](https://img.shields.io/github/release/Clinical-Genomics/MIP.svg)](https://github.com/Clinical-Genomics/MIP/releases)
[![GitHub Issues](https://img.shields.io/github/issues/Clinical-Genomics/MIP.svg)](https://github.com/Clinical-Genomics/MIP/issues)
[![CodeFactor](https://www.codefactor.io/repository/github/clinical-genomics/mip/badge)](https://www.codefactor.io/repository/github/clinical-genomics/mip)

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

MIP performs QC, alignment, coverage analysis, variant discovery and annotation, sample checks as well as ranking the found variants according to disease potential with a minimum of manual intervention. MIP is compatible with [Scout](https://github.com/Clinical-Genomics/scout) for visualization of identified variants.

MIP rare disease DNA analyses single nucleotide variants (SNVs), insertions and deletions (INDELs) and structural variants (SVs).

MIP rare disease RNA analyses mono allelic expression, fusion transcripts, transcript expression and alternative splicing.

MIP rare disease DNA vcf rerun performs re-runs starting from BCFs or VCFs.

MIP has been in use in the clinical production at the Clinical Genomics facility at Science for Life Laboratory since 2014.

## Example Usage
### MIP analyse rare disease DNA
```Bash
$ mip analyse rd_dna [case_id] --config_file [mip_config_dna.yaml] --pedigree_file [case_id_pedigree.yaml]
```

### MIP analyse rare disease DNA VCF rerun
```Bash
mip analyse rd_dna_vcf_rerun [case_id] --config_file [mip_config_dna_vcf_rerun.yaml] --vcf_rerun_file vcf.bcf  --sv_vcf_rerun_file sv_vcf.bcf --pedigree [case_id_pedigree_vcf_rerun.yaml]
```
### MIP analyse rare disease RNA
```Bash
$ mip analyse rd_rna [case_id] --config_file [mip_config_rna.yaml] --pedigree_file [case_id_pedigree_rna.yaml]
```
## Features

* Installation
  * Simple automated install of all programs using conda/docker/singularity via supplied install application
  * Downloads and prepares references in the installation process
* Autonomous
  * Checks that all dependencies are fulfilled before launching
  * Builds and prepares references and/or files missing before launching
  * Decompose and normalise reference\(s\) and variant VCF\(s\)
* Automatic
  * A minimal amount of hands-on time
  * Tracks and executes all recipes without manual intervention
  * Creates internal queues at nodes to optimize processing
* Flexible:
  * Design your own workflow by turning on/off relevant recipes in predefined pipelines
  * Restart an analysis from anywhere in your workflow
  * Process one, or multiple samples
  * Supply parameters on the command line, in a pedigree.yaml file or via config files
  * Simulate your analysis before performing it
  * Limit a run to a specific set of genomic intervals or chromosomes
  * Use multiple variant callers for both SNV, INDELs and SV
  * Use multiple annotation programs
  * Optionally split data into clinical variants and research variants
* Fast
  * Analyses an exome trio in approximately 4 h
  * Analyses a genome in approximately 21 h
* Traceability
  * Track the status of each recipe through dynamically updated status logs
  * Recreate your analysis from the MIP log or generated config files
  * Log sample meta-data and sequence meta-data
  * Log version numbers of softwares and databases
  * Checks sample integrity \(sex, contamination, duplications, ancestry, inbreeding and relationship\)
  * Test data output file creation and integrity using automated tests
* Annotation
  * Gene annotation
    * Summarize over all transcript and output on gene level
  * Transcript level annotation
    * Separate pathogenic transcripts for correct downstream annotation
  * Annotate all alleles for a position
    * Split multi-allelic records into single records to facilitate annotation
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
* [Miniconda] version 4.5.11
* [Singularity] version 3.2.1

We recommend miniconda for installing perl and cpanm libraries. However, perlbrew can also be used for installing and managing perl and cpanm libraries together with MIP.
Installation instructions and setting up specific cpanm libraries using perlbrew can be found [here](https://github.com/Clinical-Genomics/development/blob/master/docs/perl/installation/perlbrew.md).

#### Automated Installation \(Linux x86\_64\)
Below are instructions for installing the Mutation Identification Pipeline (MIP).

##### 1. Clone the official git repository

```Bash
$ git clone https://github.com/Clinical-Genomics/MIP.git
$ cd MIP
```
##### 2. Install required perl modules from cpan to a specified conda environment

```Bash
$ bash mip_install_perl.sh -e [mip] -p [$HOME/miniconda3]
```  

##### 3. Test conda and mip installation files (optional, but recommended)

```Bash
$ perl t/mip_install.test
```
A conda environment will be created where MIP with all dependencies will be installed.

##### 4. Install MIP
```Bash
$ perl mip install --environment_name [mip] --reference_dir [$HOME/mip_references]
```
This will cache the containers that are used by MIP.

###### *Note:*
  - For a full list of available options and parameters, run: ``$ perl mip install --help``

##### 6. Test your MIP installation (optional, but recommended)

Make sure to activate your MIP conda environment before executing prove.

```Bash
$ prove t -r
$ perl t/mip_analyse_rd_dna.test
```

###### When setting up your analysis config file
  A starting point for the config is provided in MIP's template directory. You will have to modify the load_env keys to whatever you named the environment. If you are using the default environment name the load_env part of the config should look like this:

  ```Yml
  load_env:
    mip:
      mip:
      method: conda
  ```

### Usage

MIP is called from the command line and takes input from the command line \(precedence\) or falls back on defaults where applicable.

Lists are supplied as repeated flag entries on the command line or in the config using the yaml format for arrays.  
Only flags that will actually be used needs to be specified and MIP will check that all required parameters are set before submitting to SLURM.

Recipe parameters can be set to "0" \(=off\), "1" \(=on\) and "2" \(=dry run mode\). Any recipe can be set to dry run mode and MIP will create the sbatch scripts, but not submit them to SLURM. MIP can be restarted from any recipe using the ``--start_with_recipe`` flag.

MIP will overwrite data files when reanalyzing, but keeps all "versioned" sbatch scripts for traceability.

You can always supply `mip [process] [pipeline] --help` to list all available parameters and defaults.

Example usage:
```Bash
$ mip analyse rd_dna case_3 --sample_ids 3-1-1A --sample_ids 3-2-1U --sample_ids 3-2-2U --start_with_recipe samtools_merge --config 3_config.yaml
```

This will analyse case 3 using 3 individuals from that case and begin the analysis with recipes after Bwa mem and use all parameter values as specified in the config file except those supplied on the command line, which has precedence.

###### Running programs in containers
Aside from a conda environment, MIP uses containers to run programs. You can use either singularity or docker as your container manager. Containers that are downloaded using MIP's automated installer will need no extra setup. By default MIP will make the reference-, outdata- and temp directory available to the container. Extra directories can be made available to each recipe by adding the key `recipe_bind_path` in the config.

In the example below the config has been modified to include the infile directories for the bwa_mem recipe:
  ```Yml
  recipe_bind_path:
    bwa_mem:
      - <path_to_directory_with_fastq_files>
  ```

#### Input

* Fastq file directories can be supplied with `--infile_dirs [PATH_TO_FASTQ_DIR=SAMPLE_ID]`
* All references and template files should be placed directly in the reference directory specified by `--reference_dir`.

##### Meta-Data

* [Configuration file] \(YAML-format\)
* [Gene panel file]
* [Pedigree file] \(YAML-format\)
* [Rank model file] \(Ini-format; SNV/INDEL\)
* [SV rank model file] \(Ini-format; SV\)
* [Qc regexp file] \(YAML-format\)

#### Output

Analyses done per individual is found in each sample_id directory and analyses done including all samples can be found in the case directory.

##### Sbatch Scripts

MIP will create sbatch scripts \(.sh\) and submit them in proper order with attached dependencies to SLURM. These sbatch script are placed in the output script directory specified by `--outscript_dir`. The sbatch scripts are versioned and will not be overwritten if you begin a new analysis. Versioned "xargs" scripts will also be created where possible to maximize the use of the cores processing power.

##### Data

MIP will place any generated data files in the output data directory specified by `--outdata_dir`. All data files are regenerated for each analysis. STDOUT and STDERR for each recipe is written in the recipe/info directory.

[Configuration file]: https://github.com/Clinical-Genomics/MIP/blob/master/templates/mip_rd_dna_config.yaml
[Gene panel file]: https://github.com/Clinical-Genomics/MIP/blob/master/templates/gene_panels.bed
[Miniconda]: http://conda.pydata.org/miniconda.html
[Pedigree file]: https://github.com/Clinical-Genomics/MIP/tree/master/templates/643594-miptest_pedigree.yaml
[Perl]:https://www.perl.org/
[Rank model file]: https://github.com/Clinical-Genomics/MIP/blob/master/templates/rank_model_-v1.31-.ini
[SV rank model file]: https://github.com/Clinical-Genomics/MIP/blob/master/templates/svrank_model_-v1.8-.ini
[Qc regexp file]: https://github.com/Clinical-Genomics/MIP/blob/master/templates/qc_regexp_-v1.26-.yaml
