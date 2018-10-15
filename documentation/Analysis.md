# MIP Analysis

You can modify all parameters to MIP in order of precedence using:

1. Command line
2. Config file
3. Pedigree file
4. Definitions file (typically not done by user)

## Start standard analysis
```Bash
$ mip analyse rare_disease --fam 0 --config GRCh37_config_-v1.4-.yaml --rio
```

``--rio`` means that one block will be performed at the nodes without transfer of files between HDS and SLURM nodes within the block. These block is:

*bamcalibrationblock*

- picardtool_mergesamfiles
- markduplicates
- gatk_baserecalibrator/gatk_applybqsr

When not supplying the ``-rio`` flag MIP will copy in and out files from HDS and SLURM nodes between each module. Thus increasing the network traffic.

## Start reannotation and reranking from vcf files
```Bash
$ mip analyse vcf_rerun --config GRCh37_config_vcf_rerun-v1.0-.yaml --vcf_rerun_file vcf_BOTH.bcf --sv_vcf_rerun_file vcf_SV.bcf
```

### Excluding a program from the analysis

```Bash
$ mip analyse rare_disease --fam 0 --config GRCh37_config_-v1.4-.yaml --rio --markduplicates 0
```

### Skipping a already processed module i.e expect that the ouput has already been generated

```Bash
$ mip analyse rare_disease --fam 0 --config GRCh37_config_-v1.4-.yaml --rio --markduplicates 2
```

### Simulate standard analysis

```Bash
$ mip analyse rare_disease --fam 0 --config GRCh37_config_-v1.4-.yaml --rio --dra
```

``-dra`` means "dry run mode" i.e simulation mode. If enabled MIP will execute everything except the final sbatch submission to SLURM and updates to qc_sample_info.yaml.

When not supplying the ``--dra`` flag MIP will launch sbatch submission to slurm.

One can use ``--dra`` to generate sbatch scripts which then can be submitted manually by the user individually or sequentially using ``sbatch --dependency=[type]:[jobid]``. Note that this will not update qc_sampleInfo.yaml as this is done at MIP run time.

### Rerun analysis using exactly the same parameters as last analysis run for family 0

```Bash
$ mip analyse rare_disease --config 0/analysis/0_config.yaml
```

### Rerun analysis using exactly the same parameters as last analysis run, but in simulation mode

```Bash
$ mip analyse rare_disease --config 0/analysis/0_config.yaml --dra
```

### Turn on individual or consecutive programs 
After performing a dry run all programs are set to simulation mode in the config. If you want run one or several programs set them to "1" in the config or supply them on the command line.
```Bash
$ mip analyse rare_disease --config 0/analysis/0_config.yaml --bwa_mem 1 --peddy 1
```

### Run all downstream dependencies starting from a program
```Bash
$ mip analyse disease --config 0/analysis/0_config.yaml --start_with_program gatk_variantrecalibration
```
This will swith the mode for all downstream dependencies to run and all programs upstream of the program to simulation mode.

### Generate all supported standard programs

```Bash
$ mip analyse rare_disease --fam 0 --config GRCh37_config_-v1.4-.yaml --rio --pp
```

This will print a string with programs in mode 2 (expect ouput) in chronological order (as far as possible, some things are processed in parallel):

```Bash
$ --split_fastq_file 2 --gzip_fastq 2 --fastqc 2 --bwa_mem 2 --picardtools_mergesamfiles 2 --markduplicates 2 --gatk_baserecalibration 2 --chanjo_sexcheck 2 --sambamba_depth 2 --picardtools_collectmultiplemetrics 2 --picardtools_collecthsmetrics 2 --cnvnator 2 --delly_call 2 --delly_reformat 2 --manta 2 --tiddit 2 --sv_combinevariantcallsets 2 --sv_varianteffectpredictor 2 --sv_vcfparser 2 --sv_rankvariant 2 --sv_reformat 2 --bcftools_mpileup 2 --freebayes 2 --gatk_haplotypecaller 2 --gatk_genotypegvcfs 2 --gatk_variantrecalibration 2 --gatk_combinevariantcallsets 2 --prepareforvariantannotationblock 2 --rhocall 2 --vt 2 --frequency_filter 2 --gatk_variantevalall 2 --gatk_variantevalexome 2 --varianteffectpredictor 2 --vcfparser 2 --snpeff 2 --peddy 2 --plink 2 --variant_integrity 2 --evaluation 2 --rankvariant 2 --endvariantannotationblock 2 --qccollect 2 --multiqc 2 --analysisrunstatus 2 --sacct 2
```

Thus you will always have the actual program names that are supported facilitating starting from any step in the analysis for instance updating qc_sampleInfo.yaml and rerunning module in bamcalibrationblock skipping markduplicates:

```Bash
$ mip analyse rare_disease --fam 0 --config GRCh37_config_-v1.4-.yaml --rio --split_fastq_file 2 --gzip_fastq 2 --fastqc 2 --bwa_mem 2 --picardtools_mergesamfiles 2 --markduplicates 0 --gatk_baserecalibration 2 --chanjo_sexcheck 2 --sambamba_depth 2 --picardtools_collectmultiplemetrics 2 --picardtools_collecthsmetrics 2 --cnvnator 2 --delly_call 2 --delly_reformat 2 --manta 2 --tiddit 2 --sv_combinevariantcallsets 2 --sv_varianteffectpredictor 2 --sv_vcfparser 2 --sv_rankvariant 2 --sv_reformat 2 --bcftools_mpileup 2 --freebayes 2 --gatk_haplotypecaller 2 --gatk_genotypegvcfs 2 --gatk_variantrecalibration 2 --gatk_combinevariantcallsets 2 --prepareforvariantannotationblock 2 --rhocall 2 --vt 2 --frequency_filter 2 --gatk_variantevalall 2 --gatk_variantevalexome 2 --varianteffectpredictor 2 --vcfparser 2 --snpeff 2 --peddy 2 --plink 2 --variant_integrity 2 --evaluation 2 --rankvariant 2 --endvariantannotationblock 2 --qccollect 2 --multiqc 2 --analysisrunstatus 2 --sacct 2
```

You can of course start or skip any number of modules as long as it is sane to do so (MIP will not check this but just execute)

### You can also modulate the mode of '--pp' using --ppm:
```	  
$ mip analyse rare_disease --fam 0 --config GRCh37_config_-v1.4-.yaml --rio --pp --ppm 1	
$ --split_fastq_file 1 --gzip_fastq 1 --fastqc 1 --bwa_mem 1 --picardtools_mergesamfiles 1 --markduplicates 1 --gatk_baserecalibration 1 --chanjo_sexcheck 1 --sambamba_depth 1 --picardtools_collectmultiplemetrics 1 --picardtools_collecthsmetrics 1 --cnvnator 1 --delly_call 1 --delly_reformat 1 --manta 1 --tiddit 1 --sv_combinevariantcallsets 1 --sv_varianteffectpredictor 1 --sv_vcfparser 1 --sv_rankvariant 1 --sv_reformat 1 --bcftools_mpileup 1 --freebayes 1 --gatk_haplotypecaller 1 --gatk_genotypegvcfs 1 --gatk_variantrecalibration 1 --gatk_combinevariantcallsets 1 --prepareforvariantannotationblock 1 --rhocall 1 --vt 1 --frequency_filter 1 --gatk_variantevalall 1 --gatk_variantevalexome 1 --varianteffectpredictor 1 --vcfparser 1 --snpeff 1 --peddy 1 --plink 1 --variant_integrity 1 --evaluation 1 --rankvariant 1 --endvariantannotationblock 1 --qccollect 1 --multiqc 1 --analysisrunstatus 1 --sacct 1
```
