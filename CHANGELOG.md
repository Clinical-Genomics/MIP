# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [7.0.0]
- New framework structure with sub commands - for analysis, install and download
- New pipelines: rd_dna (previous MIP), rd_dna_vcf_rerun (light rerun from rd_dna data) and rd_rna
- Install now has sbatch features
- Added initiation_maps for pipeline engine
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
- Updated GATK to version 4.0.10 for most GATK recipes
- Removed bed_cov and corresponding R scripts from rare disease analysis 
- Removed variantannotation block - "--rio" now only operates on BAM block
- Refactored and updated Delly to "0.7.8". Removed small-indel calling for better speed.
- Use "--use-best-n-alleles" in freebayes and added default of "4"
- Add Expansion Hunter Fix https://github.com/Clinical-Genomics/MIP/issues/442
- One case one Multiqc report Fix https://github.com/Clinical-Genomics/MIP/issues/515
- Added exclude contig option. Fix https://github.com/Clinical-Genomics/MIP/issues/509.

**New Pipeline**
- rna
- rd_dna_vcf_rerun


**New recipes**
*Rd_dna*
- sv_annotate

*rd_dna_vcf_rerun*
- sv_vcf_rerun_reformat
- vcf_rerun_reformat

**Tools**
- bcftools: 1.6 -> 1.9-h4da6232_0
- bwa: 0.7.15 -> 0.7.1-ha92aebf_3
- bwakit: 0.7.12 -> 0.7.15-1
- delly: 0.7.7 -> 0.7.8-h278814d_3
- fastqc: 0.11.4 -> 0.11.8-0
- freebayes: 1.1.0 -> 1.2.0-py27h82df9c4_3
- GATK: 3.8 -> 4.0.10.0-0
- htslib: 1.6 -> 1.9-hc238db4_4
- manta: 1.1.0 -> 1.4.0-py27_1
- peddy: 0.3.1 -> 0.4.2-py_0
- picardtools: 2.14.1 -> 2.18.14-0
- rtg-tools: 3.8.4-0 -> 3.9.1-1
- sambamba: 0.6.6 -> 0.6.8-h682856c_0
- samtools: 1.6 -> 1.9-h8ee4bcc_1
- svdb: 1.1.2 -> 1.3.0-py27h7eb728f_0
- tiddit: 2.2.5 -> 2.3.1

## [6.0.7]
- Increased ploidy when calling variants in the mitochondria. Fix https://github.com/Clinical-Genomics/MIP/issues/507
- New option: 'start_with_program'
- Fixed chdir inconsistency. Fix #402 
- Added boolean option to qccollect to evaluate plink_gender check with `-epg`.
- Moved vcfcytosure processing to the last sv module

## [6.0.0]
- Version 6.0.0. https://github.com/Clinical-Genomics/MIP/issues/186
