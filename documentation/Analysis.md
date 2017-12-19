# MIP Analysis

You can modify all parameters to MIP in order of precedence using:

1. Command line
2. Config file
3. Pedigree file
4. Definitions file (typically not done by user)

## Start standard analysis
```Bash
$ mip -f 0 -c GRCh37_config_-v1.4-.yaml -rio 1
```

``-rio 1`` means that two blocks will be performed at the nodes without transfer of files between HDS and SLURM nodes within the block. These two blocks are:

*bamcalibrationblock*

- picardtool_mergesamfiles
- markduplicates
- gatk_realignertargetcreator/indelrealigner
- gatk_baserecalibrator/printreads


*variantannotationblock*

- prepareforvariantannotationblock
- rhocall
- vt
- frequency_filter
- varianteffectpredictor
- vcfparser
- snpeff/snpsift
- rankvariants
- endvariantannotationblock

``-rio 0`` means that MIP will copy in and out files from HDS and SLURM nodes between each module. Thus increasing the network traffic.

### Excluding a program from the analysis

```Bash
$ mip -f 0 -c GRCh37_config_-v1.4-.yaml -rio 1 --pMarkduplicates 0
```

### Skipping a already processed module i.e expect that the ouput has already been generated

```Bash
$ mip -f 0 -c GRCh37_config_-v1.4-.yaml -rio 1 --pMarkduplicates 2
```

### Simulate standard analysis

```Bash
$ mip -f 0 -c GRCh37_config_-v1.4-.yaml -rio 1 -dra 1
```

``-dra`` means "dry run mode" i.e simulation mode. If enabled MIP will execute everything except the final sbatch submission to SLURM and updates to qc_sample_info.yaml.

``-dra 1`` will simulate analysis.

``-dra 0`` means no simulation

One can use ``-dra 1`` to generate sbatch scripts which then can be submitted manually by the user individually or sequentially using ``sbatch --dependency=[type]:[jobid]``. Note that this will not update qc_sampleInfo.yaml as this is done at MIP runtime.

### Rerun analysis using exactly the same parameters as last analysis run for family 0

```Bash
$ mip -c 0/analysis/0_config.yaml
```

### Rerun analysis using exactly the same parameters as last analysis run, but in simulation mode

```Bash
$ mip -c 0/analysis/0_config.yaml -dra 2
```
### Generate all supported standard programs

```Bash
$ mip -f 0 -c GRCh37_config_-v1.4-.yaml -rio 1 -pp
```

This will print a string with programs in mode 2 (expect ouput) in chronological order (as far as possible, some things are processed in parallel):

```Bash
$ --psplit_fastq_file 2 --pgzip_fastq 2 --pfastqc 2 --pbwa_mem 2 --ppicardtools_mergesamfiles 2 --pmarkduplicates 2 --pgatk_realigner 2 --pgatk_baserecalibration 2 --pchanjo_sexcheck 2 --psambamba_depth 2 --pbedtools_genomecov 2 --ppicardtools_collectmultiplemetrics 2 --ppicardtools_collecthsmetrics 2 --prcovplots 2 --pcnvnator 2 --pdelly_call 2 --pdelly_reformat 2 --pmanta 2 --ptiddit 2 --psv_combinevariantcallsets 2 --psv_varianteffectpredictor 2 --psv_vcfparser 2 --psv_rankvariant 2 --psv_reformat 2 --psamtools_mpileup 2 --pfreebayes 2 --pgatk_haplotypecaller 2 --pgatk_genotypegvcfs 2 --pgatk_variantrecalibration 2 --pgatk_combinevariantcallsets 2 --pprepareforvariantannotationblock 2 --prhocall 2 --pvt 2 --pfrequency_filter 2 --pgatk_variantevalall 2 --pgatk_variantevalexome 2 --pvarianteffectpredictor 2 --pvcfparser 2 --psnpeff 2 --ppeddy 2 --pplink 2 --pvariant_integrity 2 --pevaluation 2 --prankvariant 2 --pendvariantannotationblock 2 --pqccollect 2 --pmultiqc 2 --panalysisrunstatus 2 --psacct 2
```

Thus you will always have the actual program names that are supported facilitating starting from any step in the analysis for instance updating qc_sampleInfo.yaml and rerunning module in bamcalibrationblock skipping markduplicates:

```Bash
$ mip -f 0 -c GRCh37_config_-v1.4-.yaml -rio 1 --psplit_fastq_file 2 --pgzip_fastq 2 --pfastqc 2 --pbwa_mem 2 --ppicardtools_mergesamfiles 2 --pmarkduplicates 0 --pgatk_realigner 2 --pgatk_baserecalibration 2 --pchanjo_sexcheck 2 --psambamba_depth 2 --pbedtools_genomecov 2 --ppicardtools_collectmultiplemetrics 2 --ppicardtools_collecthsmetrics 2 --prcovplots 2 --pcnvnator 2 --pdelly_call 2 --pdelly_reformat 2 --pmanta 2 --ptiddit 2 --psv_combinevariantcallsets 2 --psv_varianteffectpredictor 2 --psv_vcfparser 2 --psv_rankvariant 2 --psv_reformat 2 --psamtools_mpileup 2 --pfreebayes 2 --pgatk_haplotypecaller 2 --pgatk_genotypegvcfs 2 --pgatk_variantrecalibration 2 --pgatk_combinevariantcallsets 2 --pprepareforvariantannotationblock 2 --prhocall 2 --pvt 2 --pfrequency_filter 2 --pgatk_variantevalall 2 --pgatk_variantevalexome 2 --pvarianteffectpredictor 2 --pvcfparser 2 --psnpeff 2 --ppeddy 2 --pplink 2 --pvariant_integrity 2 --pevaluation 2 --prankvariant 2 --pendvariantannotationblock 2 --pqccollect 2 --pmultiqc 2 --panalysisrunstatus 2 --psacct 2
```

You can of course start or skip any number of modules as long as it is sane to do so (MIP will not check this but just execute)

### You can also modulate the mode of '-pp' using -ppm:
```	  
$ mip -f 0 -c GRCh37_config_-v1.4-.yaml -rio 1 -pp -ppm 1	
$ --psplit_fastq_file 1 --pgzip_fastq 1 --pfastqc 1 --pbwa_mem 1 --ppicardtools_mergesamfiles 1 --pmarkduplicates 1 --pgatk_realigner 1 --pgatk_baserecalibration 1 --pchanjo_sexcheck 1 --psambamba_depth 1 --pbedtools_genomecov 1 --ppicardtools_collectmultiplemetrics 1 --ppicardtools_collecthsmetrics 1 --prcovplots 1 --pcnvnator 1 --pdelly_call 1 --pdelly_reformat 1 --pmanta 1 --ptiddit 1 --psv_combinevariantcallsets 1 --psv_varianteffectpredictor 1 --psv_vcfparser 1 --psv_rankvariant 1 --psv_reformat 1 --psamtools_mpileup 1 --pfreebayes 1 --pgatk_haplotypecaller 1 --pgatk_genotypegvcfs 1 --pgatk_variantrecalibration 1 --pgatk_combinevariantcallsets 1 --pprepareforvariantannotationblock 1 --prhocall 1 --pvt 1 --pfrequency_filter 1 --pgatk_variantevalall 1 --pgatk_variantevalexome 1 --pvarianteffectpredictor 1 --pvcfparser 1 --psnpeff 1 --ppeddy 1 --pplink 1 --pvariant_integrity 1 --pevaluation 1 --prankvariant 1 --pendvariantannotationblock 1 --pqccollect 1 --pmultiqc 1 --panalysisrunstatus 1 --psacct 1
```
