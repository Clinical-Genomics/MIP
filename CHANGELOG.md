# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [8.2.6]
- Updates multiqc 

**Tools**
- multiqc: 1.7 -> 1.9

## [8.2.5]
- Adds output files to store for gatk_combinevariants, sambamba depth, chromograph recipes
- Use MIPs bcftools singularity image in the conda env when checking references 

## [8.2.4]
- Chromograph patch

## [8.2.3]
- Chromograph now creates ideograms over the chromosomes
- Proband chromograph files are named with sample id instead of case id
- Increases java memory for picardtools markduplicates

## [8.2.2]
- Fixes a crash in MIP caused by not looping over the Y chromosomes for females in the GATK SplitNCigarReads recipe

## [8.2.1]
- Removes multiple identical entries in the store file. Keeps the most recent

## [8.2.0]
- Minor changes to the install processes
- Added yamllint
- New store format file in new file: <case>_deliverables.yaml
- Fixed memory error in version_collect due to GATK version command
- Internal code refactoring
- Updated docs for code style and best practise

**Tools**
- ucsc-wigToBigWig: 357
- ucsc-bedToBigBed: 357

## [8.1.0]
- STAR aligns sample fastq files with the same sequnece mode (eg. single-end or paired-end) within one job
- Dropped option to run sambamba markduplicates for RNA
- Added SMNCopyNumberCaller to MIP for SMN calling with WGS data
- Added downlod of vcf2cytosure blacklist file for grch37

**CLI**
- New CLI option for picard markduplicates optical duplicate distance. Default: 2500
- Removed sambamba markduplicates CLI options from RNA
- Removed CLI optono markduplicates_no_bam_to_cram
- New analysis recipe option for smncopynumbercaller
- New option for vcf2cytosure "--vcf2cytosure_blacklist"

**Tools**
- Arriba: 1.1.0
- SMNCopyNumberCaller: 4b2c1ad
- vcf2cytosure: 0.4.3 -> 0.5.1

## [8.0.0]

- Added BAM contigs set to allow for single ladies to get a proper check sex (no pun intended:). Actually it is for all female only analyses.
- Add new recipe to get executable version dynamically and in separate file
- Remove Snpeff and snpsift from MIP
- Add clinvar annotation as vep custom file: CLINVAR_CLNSIG,CLINVAR_CLNVID,CLINVAR_CLNREVSTAT now in CSQ field
- Move frequency annotations from Snpsift to frequency_annotation recipe
- Moves frequency annotation to separate recipe
- Adds upd for trios
- Adds rhocall viz
- Adds chromograph for chromosome visualization
- Moves CNVnator to singularity container
- Moves Manta to singularity container
- Moves VEP to singularity container
- Moves Svdb to singularity container and update version to 2.2.0
- Moves delly to singularity container
- Moved dbNSFP processing from snpsift to VEP as plugin: GERP++_NR,GERP++_RS,REVEL_rankscore,phastCons100way_vertebrate,phyloP100way_vertebrate is now part of VEP CSQ Schema instead of separate keys in the VCF INFO field 
- Install CADD via MIPs installer
- Moves STAR to singularity container
- Moves STAR-Fusion to singularity container
- Moves RSeQC to singularity container
- Moves Trim Galore to singularity container
- Removes the py3 and perl5 conda environment for the RNA pipeline
- Moves stringtie to singularity container
- Removed cramtools
- Removes cutadapt 
- Supply wgs variantcalls for ASE analysis

**New references**
- grch37_sv_frequency_vcfanno_filter_config_-v1.2-.toml
- grch37_frequency_vcfanno_filter_config_-v1.3-.toml

**Reference**
- clinvar: 20191013
- gnomad: r2.0.1 -> r2.1.1
- loqusdb: 2018-12-18 -> 2019-11-04
- expansionhunter: 3.0.0 -> 3.1.2

**Tools**
- bedtools: 2.27.1-he941832_2 -> 2.29.0=h6ed99ea_1
- chromograph:
- expansionhunter: 3.0.0 -> 3.1.2
- GATK: 4.1.2.0-1 -> 4.1.3.0-0 
- gffcompare: 0.10.6 -> 0.11.2
- manta: 1.5.0-py27_0 -> 1.6.0-py27_0 
- multiqc: 1.6 -> 1.7
- picard: 2.18.14-0 -> 2.20.7-0
- rseqc: 3.0.0 -> 3.0.1
- rhocall: 0.4 -> 0.5.1
- rtg-tools: 3.9.1-1 -> 3.10.1-0
- star: 2.6.1d -> 2.7.3a
- star-fusion: 1.5.0 -> 1.7.0
- stranger: 0.5.4 -> 0.5.5
- stringtie: 1.3.4 -> 2.0.3
- svdb: 2.0.0 -> 2.2.0
- trim-galore: 0.5.0 -> 0.6.4
- upd: 0.1
- vcfanno: 0.3.1-0 -> 0.3.2-0
- VEP: 95 -> 97

## [7.1.12]
- Remove upper case in reference file name from test data

## [7.1.11]
- Fixed bug that will casue select vcf files for snv/indel to not be produced if you turn off all SV programs
 
## [7.1.10]
- Added BAM contigs set to allow for single ladies to get a proper check sex (no pun intended:). Actually it is for all female only analyses.

## [7.1.9]
- Remove race condition between expansionhunter and sambamba depth when reading index files

## [7.1.8]
- Remove MT contig also for analysis run with analysis_type "mixed"

## [7.1.7]
- Use MQ annotation when running GATK VariantRecalibration for SNV

## [7.1.6]
- Use correct recipe name for qccollect path when adding to qc_sample_info.yaml

## [7.1.5]
- Increased sv_varianteffectpredictor memory parameter 9 -> 18 Gb

## [7.1.4]
- Fix bug in outfile_path when mitochondria contig is not part of gene panel

## [7.1.3]
- Increased sv_varianteffectpredictor memory parameter 8 -> 9 Gb

## [7.1.2]
- Update samtools_subsample_mt to fix bug in downsampling of MT bam

## [7.1.1]
- Fixed bug when skipping evaluation in QC_collect

## [7.1.0]
- Updated TIDDIT to enable faster processing
- Updated GATK for faster haplotypecalling

**Tools**
- TIDDIT: 2.5.0 -> 2.7.1
- bcftools: 1.9-h4da6232_0 -> 1.9=ha228f0b_4
- bioconductor-deseq2: 1.18.1=r3.4.1_1 -> 1.22.1=r351hf484d3e_0
- bioconductor-tximport: 1.8.0=r341_0 -> 1.12.0=r351_0 
- GATK: 4.1.0.0-0 -> 4.1.2.0-1
- samtools: 1.9-h8ee4bcc_1 -> 1.9=h8571acd_11 

## [7.0.10]
- option: tiddit_bin_size -> tiddit_coverage_bin_size
- Added generation of wig coverage data file using tiddit_coverage
- Added set_recipe_resource options for setting core_number, time, memory on CLI per recipe(s)

**New recipes**
*Rd_dna*
- tiddit_coverage

## [7.0.9]
- Removed plink memory allocation from rd_dna_vcf_rerun

## [7.0.8]
- Increased memory allocation for vep (snv/indel) again

## [7.0.7]
- Increased memory allocation for vep (snv/indel)

## [7.0.6]
- Allow "unknown" sex when using expansionhunter by then not using gender in expansionhunter CLI

## [7.0.5]
- Updated stranger to version 0.5.4 to avoid repeat id warnings

**Tools**
- stranger: 0.5.1 -> 0.5.4

## [7.0.4]
- Increased recipe memory for plink and vcf2cytosure

## [7.0.3]
- Added removal of genomicsDB dir from potential previous analysis as it causes gatk genotyping to crash

## [7.0.2]
- New framework structure with sub commands - for analysis, install and download
- New pipelines: rd_dna (previous MIP), rd_dna_vcf_rerun (light rerun from rd_dna data) and rd_rna
- Install now has sbatch features
- Download is now only sbatch
- Added initiation_maps for pipeline engines
- Changed family to case
- Changed output data dir structure to flat for each ID
- Removed call type value in file names
- Rename module time and cores to recipe time and core
- Renamed option `start_with_program` to `start_with_recipe`
- Removed use of "p" before recipe names
- Add check for recipe when using start_with_flag
- Modify parsing of pedigree to allow new RNA DE keys Fix https://github.com/Clinical-Genomics/MIP/issues/554
- Two-step model for reruns. Fix https://github.com/Clinical-Genomics/MIP/issues/546
- Add input SV vcf for rd_dna_vcf_rerun to qc_sample_info. Fix https://github.com/Clinical-Genomics/MIP/issues/548
- Added io to all recipes
- Updated GATK to version 4.1.0 for most GATK recipes
- Removed bed_cov and corresponding R scripts from rare disease analysis 
- Removed bamcalibrationblock and variantannotation block
- Removed "--rio" option
- Refactored and updated Delly to "0.7.8". Removed small-indel calling for better speed.
- Use "--use-best-n-alleles" in freebayes and added default of "4"
- Add Expansion Hunter Fix https://github.com/Clinical-Genomics/MIP/issues/442
- One case one Multiqc report Fix https://github.com/Clinical-Genomics/MIP/issues/515
- Added exclude contig option. Fix https://github.com/Clinical-Genomics/MIP/issues/509.
- Add UCSC genomicsSuperDups to annotation and rank model. Fix https://github.com/Clinical-Genomics/MIP/issues/574
- Switched to using conda instead of source with conda e.g. "conda activate [ENV]" instead of "source activate [ENV]"
- Changed default for gatk_calculategenotypeposteriors to 0 (=no). 
- Switched 1000G phase3_v4_2013-05-02 to gnomad r2.0.1 as default for gatk_calculategenotypeposteriors_support_set option
- Added switch to add all research variants to clinical file for MT. Required to display all MT variants in Scout clinical view as they are all deemed clinically relevant.
- Added gatk GatherBqsrReports to gather bqsr reports after parallelization
- Renamed flag "gatk_calculategenotypeposteriors_support_set" to "gatk_calculate_genotype_call_set"
- Added recipe_memory parameter to parameter definitions and all recipes

**New Pipeline**
- rd_dna
- rd_dna_vcf_rerun
- rd_rna

**New recipes**
*Rd_dna*
- cadd_ar
- sv_annotate

*rd_dna_vcf_rerun*
- sv_vcf_rerun_reformat
- vcf_rerun_reformat

**Annotations**
- genomic_superdups_frac_match
- REVEL_rankscore
- CADD for indels
- MaxEntScan

**References**
- clinvar: 20180429 -> 20190305
- dbnsfp: v2.9 -> v3.5a (GRCh37)
- GRCh37_gatk_merged_reference_samples.txt
- GRCh37_mip_sv_svdb_export_-2018-10-09-.vcf
- VEPs cache: 91 -> 94
- GRCh37_loqusdb_-2017-05-22-.vcf.gz -> GRCh37_loqusdb_snv_indel_-2018-12-18-.vcf.gz
- rank_model: 1.21 -> 1.24
- svrank_model: 1.3 -> 1.8
- qc_regexp: 1.19 -> 1.22 

**Tools**
- bcftools: 1.6 -> 1.9-h4da6232_0
- bedtools: 2.26.0 -> 2.27.1-he941832_2
- bwa: 0.7.15 -> 0.7.17-ha92aebf_3
- cadd: 1.4
- bwakit: 0.7.12 -> 0.7.15-1
- cutadapt: 1.14 -> 1.18-py27_0
- cramtools: 3.0.b47 -> 3.0.b127-2
- expansionhunter: 3.0.0
- delly: 0.7.7 ->0.8.1-h4037b6b_00
- fastqc: 0.11.4 -> 0.11.8-0
- freebayes: 1.1.0 -> 1.2.0-py27h82df9c4_3
- genmod: 3.7.2 -> 3.7.3
- GATK: 3.8 -> 4.1.0.0-0 (and some modules still use 3.8)
- htslib: 1.6 -> 1.9-hc238db4_4
- manta: 1.1.0 -> 1.5.0-py27_0
- multiqc: 1.4 -> 1.6
- peddy: 0.3.1 -> 0.4.2-py_0
- picardtools: 2.14.1 -> 2.18.14-0
- rtg-tools: 3.8.4-0 -> 3.9.1-1
- sambamba: 0.6.6 -> 0.6.8-h682856c_0
- samtools: 1.6 -> 1.9-h8ee4bcc_1
- stranger: 0.5.1
- svdb: 1.1.2 -> 1.3.0-py27h7eb728f_0
- tiddit: 2.2.5 -> 2.5.0
- vcf2cytosure: 0.3.0 -> 0.4.3
- vcfanno: 0.1.0 -> 0.3.1-0
- vep: 91 -> 94.4

## [6.0.7]
- Increased ploidy when calling variants in the mitochondria. Fix https://github.com/Clinical-Genomics/MIP/issues/507
- New option: 'start_with_program'
- Fixed chdir inconsistency. Fix #402 
- Added boolean option to qccollect to evaluate plink_gender check with `-epg`.
- Moved vcfcytosure processing to the last sv module

## [6.0.0]
- Version 6.0.0. https://github.com/Clinical-Genomics/MIP/issues/186
