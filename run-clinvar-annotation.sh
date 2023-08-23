# . init-cg-conda
# us
# conda activate D_mip_torb
perl mip analyse rd_dna_vcf_rerun clinvar --container_config /home/proj/stage/bin/miniconda3/envs/D_mip_torb/bin/templates/mip_container_config.yaml -c templates/clinvar-config.yaml --vcf_rerun_file cases/clinvar/clinvar_20230617-cleaned.vcf.gz --sv_vcf_rerun_file cases/clinvar/clinvar_20230617-cleaned.vcf.gz
