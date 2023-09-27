# Change Log

All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## Develop

- Store PCT_PF_READS_IMPROPER_PAIRS in qc file
- Evaluate and fail analysis where PCT_PF_READS_IMPROPER_PAIRS is above 5 %

## [12.0.0]

- Set overlap requirement for merging two SVs to 0.8, down from the 0.95 default
- Set overlap requirement for SV annotation to 0.6, up from 0.5
- Update Tiddit to improve SV positioning
- Increased memory allocation for salmon and picardtools_mergersamfiles (RNA)
- New version of MegaFusion. A bug in the previous version prevented SVDB from writing the format and sample field in the vcf.
- Remove the RSeQC read duplication analysis as it often fails.
- Increased run time allocation for gatk_asereadcounter.
- Increased the default required length for a trimmed rna read to be retained from 20 bp to 40 bp. Configurable via CLI or config.
- Fixed a bug where gnomad SV version 2.0 instead of version 2.1 was used to annotate SVs
- One-pass instead of two-pass mapping with STAR-Fusion, as recommended for STAR-Fusion 1.12
- Bump max run time for retroseq to 15 hours.

### Tools

- Arriba 2.3.0 -> 2.4.0
- DeepVariant 1.4.0 -> 1.5.0
- FastQC: 0.11.9 -> 0.12.1
- GATK: 4.2.6.1 -> 4.4.0.0
- Gffcompare 0.11.2 -> 0.12.6
- Htslib 1.15.1 -> 1.16
- MegaFusion 66a3a80 -> d3feacf
- Picard 2.27.2 -> 2.27.5
- STAR-Fusion 1.11.0 -> 1.12.0
- SVDB: 2.7.0 -> 2.8.2
- Tiddit 3.3.2 -> 3.6.0

### Databases

clinvar: 20220829 -> 20230508
loqusdb snapshot: 20230208 -> 20230512
hmtvar: oct2022

## [11.2.2]

- New patch of gens pre processing container.

### Tools

- gens_preproc 1.0.8 -> 1.0.11

## [11.2.1]

- Patching of gens pre processing container to solve an issue with incomplete bed files.

### Tools

- gens_preproc 1.0.2 -> 1.0.8

## [11.2.0]

- Adds retroseq for mobile element detection

### Databases

- expansionhunter variant catalog: v4.0.2 -> v5.0.0

### Tools

- RetroSeq: 9d4f3b5

## [11.1.3]

- Adds Gens' bed index file to deliverables

## [11.1.2]

- Fixed a bug in the mitochondrial deletion recipe affecting samples with 0 intermediate discordant MT reads.

## [11.1.1]

- Updates chromgraph to version 1.3.1
- Added vcfanno config version 1.18 to the download config

## [11.1.0]

- Save raw files from ExpansionHunter
- Run UPD and subsequently chromograph on unaffected children
- Annotate SV variants with the caller that reported the variant
- Produce files for CNV analysis in Gens
- Updated SO terms for new version of VEP
- ExACpLI -> pLI, see [vep issue 108](https://github.com/Ensembl/VEP_plugins/issues/108)
- Use REVEL_score rather than REVEL_rankscore for the ranking algorithm
- Use BWA-mem2 instead of BWA mem for mapping
- Set default annotation overlap for structural variants to 0.5 (previously 0.8), due to change in TIDDIT
- Turn on Stringtie and gffcompare by default
- Run varg on research vcf
- Increase max for coverage calculation to 500x
- Separate list of ranked SO terms for structural variants to ensure that the right SO term gets picked as the most severe for SVs
- Adds option to use bedpe files with svdb query

### Tools

- Arriba: 2.1.0 -> 2.3.0
- Chromograph 1.1.4 -> 1.3.0
- DeepVariant: 1.1.0 -> 1.4.0
- ExpansionHunter: 4.0.2 -> 5.0.0
- GATK: 4.2.2.0 -> 4.2.6.1
- HTSlib: 1.13 -> 1.15.1
- MultiQC: 1.11 -> 1.12
- Peddy: 0.4.3 -> 0.4.8
- Picard: 2.25.0 -> 2.27.2
- SMNCopyNumberCaller 1.1.1 -> 1.1.2
- Star Fusion: 1.10.1 -> 1.11.0
- Stranger: 0.8.0 -> 0.8.1
- Stringtie: 2.1.3b -> 2.2.1
- Tiddit: 2.12.1 -> 3.3.2
- Trimgalore: 0.6.4 -> 0.6.7
- VEP: 104.3 -> 107.0
- svdb: 2.4.0 -> 2.7.0
- vcf2cytosure v0.5.1 -> v0.8

### Databases

- clinvar: 20211010 -> 20220829
- dbnsfp: 4.1a -> 4.3a (grch38 only)
- gnomad: r3.1.1 -> r3.1.2 (grch38 only)
- giab: 3.3.2 -> 4.2.1
- loqusdb dump: 20210921 -> c
- nist: v3.3.2 -> v4.2.1
- vcf2cytosure blacklist: 200520

## [11.0.3]

- Initiate conda prior to activation

## [11.0.2]

Updates chromograph

### Tools

chromograph 1.1.4 -> 1.1.5

## [11.0.1]

- When running Deepvariant, set tmpdir to analysis folder and use `intermediate_results_dir`.
- When running DeepVariant via singularity, use the options `--no-home` and `--cleanenv`.

## [11.0.0]

- HmtNote: annotate mitochondrial variants in VCF file
- Updating to latest and greatest versions
- Mitochondrial deletion analysis
- GATK Haplotypecaller has been turned off in favour of Deepvariant
- Introduces possibility to store singularity images locally as a .sif file
- Increased allele frequency cut off for when a variant is filtered out to 0.7
- Turned off Star_caller and Telomerecat by default

### Tools

- cyrius v1.1 -> v1.1.1
- deepvariant 1.1.0 -> 1.2.0
- deeptrio 1.1.0 -> 1.2.0
- gatk 4.2.0.0 -> 4.2.2.0
- glnexus v1.3.1 -> v1.4.1
- HmtNote: 0.7.2
- htslib: 1.10.2 -> 1.13
- multiqc 1.10.1 -> v1.11
- star-fusion 1.10.0 -> 1.10.1
- vep release_103.1 -> release_104.3

### References

- gnomad: r3.0 -> r3.1.1
- [NEW] gnomad mt: r3.1
- clinvar: 20210415 -> 20211010

## [10.2.5]

- Allow slurm quality of service flag to be set to 'express'

## [10.2.4]

- Split Star-Fusion alignment and detection into two recipes
- Use temp directory with Star-Fusion
- Resource bump for RNA
- Limit memory for glnexus
- Use non-gpu version of Deepvariant by default

## [10.2.3]

- Updates Chromograph to version 1.1.4

## [10.2.5]

- Allow slurm quality of service flag to be set to 'express'

## [10.2.4]

- Split Star-Fusion alignment and detection into two recipes
- Use temp directory with Star-Fusion
- Resource bump for RNA
- Limit memory for glnexus
- Use non-gpu version of Deepvariant by default

## [10.2.3]

- Updates Chromograph to version 1.1.4

## [10.2.2]

- Adds missing median coverage metrics to metrics deliverable file

## [10.2.1]

- Removed frequency filtering for mitochondrial sites

## [10.2.0]

- Introduced the option `--start_after_recipe <recipe_name>` to start the pipeline after a given recipe

## [10.1.0]

- Only store qc_metrics_deliverables path in file store for downstream parsing

## [10.0.3]

- Remove duplicates from Glnexus output

## [10.0.2]

- Glnexus are used to genotype the gvcf regardless of how many samples that are analysed.
- Resource bump to the MIP RNA recipe fusion_report

## [10.0.1]

- Fix to gene panel regexp
- Use automated build of MIP
- Increased memory allocation for version_collect

### Tools

stranger: 0.7.1 -> 0.8.0

## [10.0.0]

- Remove unused recipe split_fastq_file
- Align with the same bwa mem options as used by Broad
  - Fixed number of bases in each batch
  - Use bwa mem instead of run-bwamem for alignment to grch38
- Use Chanjo repos Docker file instead of MIPs
- Removed support to run bcftools_mpileup as a variant caller
- Added perl and MIP docker file and use it in recipes
- Removed PATH check for proxy bins
- Switched to using container manager instead of proxy bins in recipes
- Removed support to run variant_integrity
- Added deepvariant as variant caller and glnexus to merge samples to case vcf

### Tools

arriba: v1.2.0 -> v2.1.0
bedtools: 2.29.0 -> 2.30.0
bwa-mem2 2.0 -> 2.2
chanjo: 4.2.0 -> 4.6
chromograph: 1.0.1 -> 1.1
cyrius: v1.0 -> v1.1
deepvariant: 1.1.0
delly: 0.8.1 -> 0.8.7
gatk: v4.1.8.1 -> v4.2.0.0
glnexus: v1.3.1
megafusion: 66a3a80
multiqc: 1.9 -> 1.10.1
pdfmerger: v1.0
picardtools: 2.23.4 -> 2.25.0
preseq: 2.0.3 -> 3.1.2
rseqc: 3.0.1 -> 4.0.0
rtg-tools: 3.10.1 -> 3.12
salmon: 0.12.0 -> 1.4.0
smncopynumbercaller: 4b2c1ad -> v1.1.1
star: 2.7.4a -> 2.7.8a
stranger: 0.7 -> 0.7.1
svdb: 2.2.0 -> 2.4.0--py37h77a2a36_4
tiddit: 2.8.1 -> 2.12.1
upd: 0.1 -> 0.1.1

### References

- clinvar_20200905 -> clinvar_20210415
- VEP cache: 100 -> 103.1
- grch37_gencode_v19_ctat_lib_plug-n-play_-apr032020-.tar.gz -> grch37_gencode_v19_CTAT_lib_Mar012021.plug-n-play.tar.gz
- grch38_gencode_v31_ctat_lib_plug-n-play_-apr062020-.tar.gz -> grch38_gencode_v37_CTAT_lib_Mar012021.plug-n-play.tar.gz
- gencode annotation: v34 -> v37

## [9.1.3]
- Fix memory allocation for mip-rna markduplicates.
- Update to repeat expansion calling
- Adds possibility to rename arriba fusion report from with sample display name
- Updates documentation

**Tools**
expansionhunter 3.1.2 -> 4.0.2
stranger 0.5.5 -> 0.7

## [9.1.2]
- Increase markduplicates java memory allocation for chromosome 2.
- Turn off chromograph_viz for wes analysis
- Chromograph exits gracefully on empty infile
- Use median coverage instead of expected coverage when evaluating whether expected coverage has been reached or not [#1719](https://github.com/Clinical-Genomics/MIP/issues/1719)
- Fixes infile error in rtg_vcfeval recipe

**Tools**
chromograph 1.0.0 -> 1.0.1

## [9.1.1]
- Fix MIP's gender estaimation for wgs samples with gender set to unknown.

## [9.1.0]
- PNGs generated by chromograph are now uniform in dpi and image size
- Adds chromograph recipe to generate images of rhocall viz output (regions of autozygosity and fraction of homozygous snps)

**Tools**
chromograph 0.3.3 -> 1.0.0

## [9.0.6]
- Use "PAN" key for slurm_jobs_ids file instead of "ALL" as "ALL" has a size constraint

## [9.0.5]
- Predicted gender from wgs samples are now used in the generated fam files.
- Restrict plink analysis to intersected target capture kits for mixed wgs/wes cases
- Update to chromograph in order to fix the renderering coverage images

**Tools**
chromograph: 0.3.1 -> 0.3.3

## [9.0.4]
- Increased memory allocation for samtools_subsample_mt
- Check that vep plugin paths exists prior to executing mip
- Cd into cadd temp directory before executing cadd in order to escape snakemake lock errors

## [9.0.3]
- Changed path and name of slurm job ids file to facilitate analysis monitoring
- SpliceAI annotation with VEP instead of vcfanno
- Files from Chromograph are no longer compressed into a tarball
- Sample specific naming of outfiles from rhocall viz
- Use temporary contig directory for CADD indel calling in order to avoid race condition

**Tools**
chromograph: 0.1.3 -> 0.3.1

## [9.0.2]
**References**
- clinvar_20200728 -> clinvar_20200905

## [9.0.1]
- Use sample_id in smncopynumber caller instead of file_name_prefix

## [9.0.0]
- Moved annotationof CADD and SPIDEX to vcfanno´s toml config
- Removed CADD and SPIDEX annotations from Rankvariants recipe, CLI and parameters
- Turned off bcftools_mpileup by default
- Replaced sambamba sort with samtools sort after alignment
- Replaced recipe picartools_mergesamfiles with samtools_merge
- Replaced sambamba flagstat with samtools flagstat in markduplicates recipe
- Rename frequency_annotation to variant_annotation
- Removed option to run sambamba markduplicates for markduplicates recipe
- Added SpliceAI annotation
- Collect and evaluate QC metrics generated in the RNA pipeline
- Per default MIP now installs all programs needed for the different pipelines into one conda environment
- Add picardtools CollectRnaSeqMetrics to the RNA pipeline
- Call CYP2D6 alleles with star_caller from the Cyrius package
- Added bwa_mem2 as an alignemnt option instead of bwa_mem
- Added option "genomicsdb-shared-posixfs-optimizations" to gatk_genomicsDB to turn off file lock
- Moved smncopynumbercaller from sample level to case level
- Added telomerecat analysis for estimating telomere length from wgs

**Tools**
- Arriba: 1.1.0 -> 1.2.0
- bcftools: 1.9=ha228f0b_4 -> 1.10.2-hd2cd319_0 (DNA)
- bwa-mem2: 2.0-he513fc3_0
- CADD: v1.5 -> v1.6
- Cyrius: 1.0
- expansionhunter: 3.1.2 -> 3.2.2
- fastqc: 0.11.8-0 -> 0.11.9
- gatk: 4.1.3.0 -> 4.1.8.1
- htslib: 1.9-hc238db4_4 -> 1.10.2=h78d89cc_0 (DNA)
- picard: 2.20.7 -> 2.22.4
- samtools: 1.9=h8571acd_11 -> 1.10-h9402c20_2 (DNA)
- SMNCopyNumberCaller: 4b2c1ad -> 1.0
- STAR 2.7.3a -> 2.7.4a
- STAR-Fusion v1.8.0 -> v1.9.0
- stringtie 2.0.3 -> 2.1.3b
- VEP: 97 -> 100

**References**
- clinvar_20191013 -> clinvar_20200728
- dbNSFP4.0b2a.zip -> dbNSFP4.1a.zip
- delly_exclude grch37 20150227 -> 20200310
- grch37_frequency_vcfanno_filter_config_-v1.3-.toml -> grch37_vcfanno_config_-v1.10-.toml
- grch37_gencode_annotation_-v31-.gtf.gz -> grch37_gencode_annotation_-v34-.gtf.gz
- grch37_gencode_v19_ctat_lib_plug-n-play_-oct012019-.tar.gz -> grch37_gencode_v19_ctat_lib_plug-n-play_-apr032020-.tar.gz
- grch37_loqusdb_snv_indel_-2019-11-04-.vcf.gz -> grch37_loqusdb_snv_indel_-2020-03-24-.vcf.gz
- grch37_loqusdb_sv_-2020-04-20.vcf
- grch38_frequency_vcfanno_filter_config_-v1.2-.toml -> grch38_frequency_vcfanno_filter_config_-v1.3-.toml
- grch38_gencode_annotation_-v31-.gtf.gz -> grch38_gencode_annotation_-v34-.gtf.gz
- grch38_gencode_v31_ctat_lib_plug-n-play_-oct012019-.tar.gz -> grch38_gencode_v31_ctat_lib_plug-n-play_-apr062020-.tar.gz
- VEP cache: 97 -> 100

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
