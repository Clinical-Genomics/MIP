# MIP Analysis

You can modify all parameters to MIP in order of precedence using:

1. Command line
2. Config file
3. Pedigree file
4. Definitions file (typically not done by user)

## Start standard analysis
```Bash
$ mip analyse rd_dna [case_id] --config_file [mip_config_rd_dna]
```

## Start reannotation and reranking from vcf files
```Bash
$ mip analyse rd_dna_vcf_rerun [case_id] --config_file [mip_config_rd_dna_vcf_rerun.yaml] --pedigree [case_id_pedigree_vcf_rerun.yaml] --vcf_rerun_file vcf_snv_indel.bcf --sv_vcf_rerun_file vcf_SV.bcf
```

### Excluding a recipe from the analysis

```Bash
$ mip analyse rd_dna [case_id] --config_file [mip_config_rd_dna] --markduplicates 0
```

### Skipping a already processed recipe i.e expect that the output has already been generated

```Bash
$ mip analyse rd_dna [case_id] --config_file [mip_config_rd_dna] --markduplicates 2
```

### Simulate an analysis

```Bash
$ mip analyse rd_dna [case_id] --config_file [mip_config_rd_dna] --dra
```

``-dra`` means "dry run all" i.e simulation mode. If enabled MIP will execute everything except the final sbatch submission to SLURM and updates to qc_sample_info.yaml.

When not supplying the ``--dra`` flag MIP will launch sbatch submission to SLURM.

One can use ``--dra`` to generate sbatch scripts which then can be submitted manually by the user individually or sequentially using ``sbatch --dependency=[type]:[jobid]``. Note that this will not update qc_sample_info.yaml as this is done at MIP run time.

### Rerun analysis using exactly the same parameters as last analysis run for [case_id]

```Bash
$ mip analyse rd_dna [case_id] --config_file [case_id]/analysis/[case_id]_config.yaml
```

### Rerun analysis using exactly the same parameters as last analysis run, but in simulation mode

```Bash
$ mip analyse rd_dna [case_id]--config_file [case_id]/analysis/[case_id]_config.yaml --dra
```

### Turn on individual or consecutive recipes
After performing a dry run all recipes are set to simulation mode in the config. If you want run one or several recipes set them to "1" in the config or supply them on the command line.
```Bash
$ mip analyse rd_dna [case_id] --config_file [case_id]/analysis/[case_id]_config.yaml --bwa_mem 1 --peddy_ar 1
```

### Run all downstream dependencies starting from a recipe
```Bash
$ mip analyse rd_dna [case_id] --config_file [case_id]/analysis/[case_id]_config.yaml --start_with_recipe gatk_variantrecalibration
```
This will switch the mode for all downstream dependencies to run and all recipes upstream of the recipe to simulation mode.

### Generate all supported standard recipes

```Bash
$ mip analyse rd_dna [case_id] --config_file [mip_config_rd_dna] --print_recipe
```

This will print a string with recipes in mode 2 (expect output) in chronological order (as far as possible, as some things are processed in parallel):

```Bash
$ --gzip_fastq 2 --fastqc_ar 2 --bwa_mem 2 --samtools_merge 2 --markduplicates 2 --gatk_baserecalibration 2 --chanjo_sexcheck 2 --sambamba_depth 2 --picardtools_collectmultiplemetrics 2 --picardtools_collecthsmetrics 2 --cnvnator_ar 2 --delly_call 2 --delly_reformat 2 --manta 2 --tiddit 2 --sv_combinevariantcallsets 2 --sv_varianteffectpredictor 2 --sv_vcfparser 2 --sv_rankvariant 2 --sv_reformat 2 --gatk_haplotypecaller 2 --gatk_genotypegvcfs 2 --gatk_variantrecalibration 2 --gatk_combinevariantcallsets 2 --prepareforvariantannotationblock 2 --rhocall_ar 2 --vt_ar 2 --frequency_filter 2 --gatk_variantevalall 2 --gatk_variantevalexome 2 --varianteffectpredictor 2 --vcfparser 2 --peddy_ar 2 --plink 2 --rankvariant 2 --endvariantannotationblock 2 --qccollect_ar 2 --multiqc_ar 2 --analysisrunstatus 2 --sacct 2
```

Thus you will always have the actual recipe names that are supported facilitating starting from any step in the analysis for instance updating qc_sample_info.yaml and rerunning some recipes, but skipping markduplicates:

```Bash
$ mip analyse rd_dna [case_id] --config_file [mip_config_rd_dna] --gzip_fastq 2 --fastqc_ar 2 --bwa_mem 2 --samtools_merge 2 --markduplicates 0 --gatk_baserecalibration 1 --chanjo_sexcheck 2 --sambamba_depth 2 --picardtools_collectmultiplemetrics 2 --picardtools_collecthsmetrics 2 --cnvnator_ar 2 --delly_call 2 --delly_reformat 2 --manta 2 --tiddit 2 --sv_combinevariantcallsets 2 --sv_varianteffectpredictor 2 --sv_vcfparser 2 --sv_rankvariant 2 --sv_reformat 2 --gatk_haplotypecaller 2 --gatk_genotypegvcfs 2 --gatk_variantrecalibration 2 --gatk_combinevariantcallsets 2 --prepareforvariantannotationblock 2 --rhocall_ar 2 --vt_ar 2 --frequency_filter 2 --gatk_variantevalall 2 --gatk_variantevalexome 2 --varianteffectpredictor 2 --vcfparser 2 --peddy_ar 2 --plink 2 --rankvariant 2 --endvariantannotationblock 2 --qccollect_ar 2 --multiqc 2 --analysisrunstatus 2 --sacct 2
```

You can of course start or skip any number of recipes as long as it is sane to do so (MIP will not check this but just execute)

### You can also modulate the mode of '--print_recipe' using --print_recipe_mode:
```  
$ mip analyse rd_dna [case_id] --config_file [mip_config_rd_dna] --print_recipe --print_recipe_mode 1
$ --gzip_fastq 1 --fastqc_ar 1 --bwa_mem 1 --samtools_merge 1 --markduplicates 1 --gatk_baserecalibration 1 --chanjo_sexcheck 1 --sambamba_depth 1 --picardtools_collectmultiplemetrics 1 --picardtools_collecthsmetrics 1 --cnvnator_ar 1 --delly_call 1 --delly_reformat 1 --manta 1 --tiddit 1 --sv_combinevariantcallsets 1 --sv_varianteffectpredictor 1 --sv_vcfparser 1 --sv_rankvariant 1 --sv_reformat 1 --gatk_haplotypecaller 1 --gatk_genotypegvcfs 1 --gatk_variantrecalibration 1 --gatk_combinevariantcallsets 1 --prepareforvariantannotationblock 1 --rhocall_ar 1 --vt_ar 1 --frequency_filter 1 --gatk_variantevalall 1 --gatk_variantevalexome 1 --varianteffectpredictor 1 --vcfparser 1 --peddy_ar 1 --plink 1 --rankvariant 1 --endvariantannotationblock 1 --qccollect_ar 1 --multiqc 1 --analysisrunstatus 1 --sacct 1
```
