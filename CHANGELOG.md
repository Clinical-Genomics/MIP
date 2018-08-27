# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [7.0.0]
- Add input SV vcf for vcf_rerun to qc_sample_info. Fix https://github.com/Clinical-Genomics/MIP/issues/548
- Add Expansion Hunter Fix https://github.com/Clinical-Genomics/MIP/issues/442
- One case one Multiqc report Fix https://github.com/Clinical-Genomics/MIP/issues/515
- Added exclude contig option. Fix https://github.com/Clinical-Genomics/MIP/issues/509.

**New Pipeline**
- rna
- vcf_rerun


**New recipes**
*Rare_disease*
- sv_annotate

*vcf_rerun*
- sv_vcf_rerun_reformat
- vcf_rerun_reformat

## [6.0.7]
- Increased ploidy when calling variants in the mitochondria. Fix https://github.com/Clinical-Genomics/MIP/issues/507
- New option: 'start_with_program'
- Fixed chdir inconsistency. Fix #402 
- Added boolean option to qccollect to evaluate plink_gender check with `-epg`.
- Moved vcfcytosure processing to the last sv module

## [6.0.0]
- Version 6.0.0. https://github.com/Clinical-Genomics/MIP/issues/186
