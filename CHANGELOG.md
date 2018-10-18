# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [7.0.0]
- New framework structure with sub commands - for analysis, install and download
- New pipelines: rd_dna (previous MIP), rd_dna_vcf_rerun (light rerun from rd_dna data) and rd_rna
- Install now has sbatch features
- Added initiation_maps for pipeline engine
- Changed output data dir structure to flat for each ID
- Removed call type value in file names
- Removed use of "p" before program names
- Add check for program when using start_with_flag
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
- bwa: 0.7.15 -> 0.7.1-ha92aebf_3
- bwakit: 0.7.12 -> 0.7.15-1
- delly: 0.7.7 -> 0.7.8-h278814d_3
- fastqc: 0.11.4 -> 0.11.8-0
- GATK: 3.8 -> 4.0.10.0-0
- picardtools: 2.14.1 -> 2.18.14-0
- samtools: 1.6 -> 1.9-h8ee4bcc_1


## [6.0.7]
- Increased ploidy when calling variants in the mitochondria. Fix https://github.com/Clinical-Genomics/MIP/issues/507
- New option: 'start_with_program'
- Fixed chdir inconsistency. Fix #402 
- Added boolean option to qccollect to evaluate plink_gender check with `-epg`.
- Moved vcfcytosure processing to the last sv module

## [6.0.0]
- Version 6.0.0. https://github.com/Clinical-Genomics/MIP/issues/186
