package MIP::Recipes::Pipeline::Analyse_rd_dna_panel;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use List::MoreUtils qw { any };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME $SPACE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ parse_rd_dna_panel pipeline_analyse_rd_dna_panel };
}

sub parse_rd_dna_panel {

## Function : Rare disease panel DNA pipeline specific checks and parsing
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $broadcasts_ref        => Holds the parameters info for broadcasting later {REF}
##          : $file_info_href        => File info hash {REF}
##          : $order_parameters_ref  => Order of parameters (for structured output) {REF}
##          : $parameter_href        => Parameter hash {REF}
##          : $sample_info_href      => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $broadcasts_ref;
    my $file_info_href;
    my $order_parameters_ref;
    my $parameter_href;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        broadcasts_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$broadcasts_ref,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        order_parameters_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$order_parameters_ref,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Active_parameter qw{
      check_sample_id_in_hash_parameter
      check_sample_id_in_hash_parameter_path
      parse_infiles
      parse_vep_plugin
      set_vcfparser_outfile_counter
      write_references
    };
    use MIP::Analysis qw{ broadcast_parameters parse_prioritize_variant_callers };
    use MIP::Constants qw{ set_container_constants };
    use MIP::Config qw{ write_mip_config };
    use MIP::Environment::Container qw{ parse_containers };
    use MIP::Fastq qw{ parse_fastq_infiles };
    use MIP::File_info qw{ check_parameter_metafiles };
    use MIP::Gatk qw{ check_gatk_sample_map_paths };
    use MIP::Reference qw{ parse_exome_target_bed parse_nist_parameters };
    use MIP::Sample_info qw{ set_parameter_in_sample_info };
    use MIP::Vep qw{
      check_vep_api_cache_versions
      check_vep_custom_annotation
    };
    use MIP::Vcfanno qw{ parse_toml_config_parameters };

    ## Constants
    Readonly my @REMOVE_CONFIG_KEYS => qw{ associated_recipe };

    ## Set analysis constants
    set_container_constants( { active_parameter_href => $active_parameter_href, } );

    parse_containers(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
        }
    );

    ## Update exome_target_bed files with human_genome_reference_source and human_genome_reference_version
    parse_exome_target_bed(
        {
            exome_target_bed_file_href     => $active_parameter_href->{exome_target_bed},
            human_genome_reference_source  => $file_info_href->{human_genome_reference_source},
            human_genome_reference_version => $file_info_href->{human_genome_reference_version},
        }
    );

    ## Checks parameter metafile exists and set build_file parameter
    check_parameter_metafiles(
        {
            active_parameter_href => $active_parameter_href,
            file_info_href        => $file_info_href,
            parameter_href        => $parameter_href,
        }
    );

    ## Update the expected number of outfiles after vcfparser
    set_vcfparser_outfile_counter( { active_parameter_href => $active_parameter_href, } );

    ## Check that VEP directory and VEP cache match
    check_vep_api_cache_versions(
        {
            vep_directory_cache => $active_parameter_href->{vep_directory_cache},
        }
    );

    ## Check VEP custom annotations options
    check_vep_custom_annotation(
        {
            vep_custom_ann_href => \%{ $active_parameter_href->{vep_custom_annotation} },
        }
    );

    parse_vep_plugin(
        {
            active_parameter_href => $active_parameter_href,
            mip_vep_plugins_ref   => [qw{ vep_plugin }],
        }
    );

    ## Check sample_id provided in hash parameter is included in the analysis
    check_sample_id_in_hash_parameter(
        {
            active_parameter_href => $active_parameter_href,
            parameter_names_ref   => [qw{ analysis_type expected_coverage }],
            parameter_href        => $parameter_href,
            sample_ids_ref        => \@{ $active_parameter_href->{sample_ids} },
        }
    );

    ## Check sample_id provided in hash path parameter is included in the analysis and only represented once
    check_sample_id_in_hash_parameter_path(
        {
            active_parameter_href => $active_parameter_href,
            parameter_names_ref   => [qw{ exome_target_bed infile_dirs }],
            sample_ids_ref        => \@{ $active_parameter_href->{sample_ids} },
        }
    );

    ## Check that the supplied gatk sample map file paths exists
    check_gatk_sample_map_paths(
        {
            gatk_genotypegvcfs_mode => $active_parameter_href->{gatk_genotypegvcfs},
            sample_map_path         => $active_parameter_href->{gatk_genotypegvcfs_ref_gvcf},
        }
    );

    ## Parse parameters with TOML config files
    parse_toml_config_parameters(
        {
            active_parameter_href => $active_parameter_href,
        }
    );

    parse_nist_parameters(
        {
            active_parameter_href => $active_parameter_href,
        }
    );

    broadcast_parameters(
        {
            active_parameter_href => $active_parameter_href,
            broadcasts_ref        => $broadcasts_ref,
            order_parameters_ref  => $order_parameters_ref,
        }
    );

    ## Write references for this analysis to yaml
    write_references(
        {
            active_parameter_href => $active_parameter_href,
            outfile_path          => $active_parameter_href->{reference_info_file},
            parameter_href        => $parameter_href,
        }
    );

    ## Check that all active variant callers have a prioritization order and that the prioritization elements match a supported variant caller
    ## parse_prioritize_variant_callers(
    ##    {
    ##        active_parameter_href => $active_parameter_href,
    ##        parameter_href        => $parameter_href,
    ##     }
    ## );

    ## Write config file for case
    write_mip_config(
        {
            active_parameter_href => $active_parameter_href,
            remove_keys_ref       => \@REMOVE_CONFIG_KEYS,
            sample_info_href      => $sample_info_href,
        }
    );

    ## Get the ".fastq(.gz)" files from the supplied infiles directory. Checks if the files exist
    parse_infiles(
        {
            active_parameter_href => $active_parameter_href,
            file_info_href        => $file_info_href,
        }
    );

    ## Reformat file names to MIP format, get file name info and add info to sample_info
    parse_fastq_infiles(
        {
            active_parameter_href => $active_parameter_href,
            file_info_href        => $file_info_href,
            sample_info_href      => $sample_info_href,
        }
    );

    ## Add to sample info
    set_parameter_in_sample_info(
        {
            active_parameter_href => $active_parameter_href,
            file_info_href        => $file_info_href,
            sample_info_href      => $sample_info_href,
        }
    );

    return;
}

sub pipeline_analyse_rd_dna_panel {

## Function : Pipeline recipe for dna panel data analysis.
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $broadcasts_ref        => Holds the parameters info for broadcasting later {REF}
##          : $file_info_href        => File info hash {REF}
##          : $job_id_href           => Job id hash {REF}
##          : $log                   => Log object to write to
##          : $order_parameters_ref  => Order of parameters (for structured output) {REF}
##          : $order_recipes_ref     => Order of recipes
##          : $parameter_href        => Parameter hash {REF}
##          : $sample_info_href      => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $broadcasts_ref;
    my $file_info_href;
    my $job_id_href;
    my $log;
    my $order_parameters_ref;
    my $order_recipes_ref;
    my $parameter_href;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        broadcasts_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$broadcasts_ref,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        order_parameters_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$order_parameters_ref,
            strict_type => 1,
        },
        order_recipes_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$order_recipes_ref,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Analysis qw{ set_recipe_bwa_mem };
    use MIP::Log::MIP_log4perl qw{ log_display_recipe_for_user };
    use MIP::Parse::Reference qw{ parse_references };

    ## Recipes
    use MIP::Recipes::Analysis::Analysisrunstatus qw{ analysis_analysisrunstatus };
    use MIP::Recipes::Analysis::Bcftools_norm qw{ analysis_bcftools_norm_panel };
    use MIP::Recipes::Analysis::Bwa_mem qw{ analysis_bwa_mem2 };
    use MIP::Recipes::Analysis::Cadd qw{ analysis_cadd_panel };
    use MIP::Recipes::Analysis::Endvariantannotationblock
      qw{ analysis_endvariantannotationblock_panel };
    use MIP::Recipes::Analysis::Fastqc qw{ analysis_fastqc };
    use MIP::Recipes::Analysis::Frequency_filter qw{ analysis_frequency_filter_panel };
    use MIP::Recipes::Analysis::Gatk_baserecalibration qw{ analysis_gatk_baserecalibration_panel };
    use MIP::Recipes::Analysis::Gatk_combinevariantcallsets
      qw{ analysis_gatk_combinevariantcallsets };
    use MIP::Recipes::Analysis::Gatk_haplotypecaller qw{ analysis_gatk_haplotypecaller_panel};
    use MIP::Recipes::Analysis::Gzip_fastq qw{ analysis_gzip_fastq };
    use MIP::Recipes::Analysis::Gatk_gathervcfs qw{ analysis_gatk_gathervcfs };
    use MIP::Recipes::Analysis::Gatk_genotypegvcfs qw{ analysis_gatk_genotypegvcfs };
    use MIP::Recipes::Analysis::Gatk_variantevalall qw{ analysis_gatk_variantevalall };
    use MIP::Recipes::Analysis::Gatk_variantrecalibration
      qw{ analysis_gatk_variantrecalibration_wes };
    use MIP::Recipes::Analysis::Markduplicates qw{ analysis_markduplicates_panel };
    use MIP::Recipes::Analysis::Mip_qccollect qw{ analysis_mip_qccollect };
    use MIP::Recipes::Analysis::Mip_vcfparser qw{ analysis_mip_vcfparser_panel };
    use MIP::Recipes::Analysis::Mip_vercollect qw{ analysis_mip_vercollect };
    use MIP::Recipes::Analysis::Multiqc qw{ analysis_multiqc };
    use MIP::Recipes::Analysis::Picardtools_collecthsmetrics
      qw{ analysis_picardtools_collecthsmetrics };
    use MIP::Recipes::Analysis::Picardtools_collectmultiplemetrics
      qw{ analysis_picardtools_collectmultiplemetrics };
    use MIP::Recipes::Analysis::Rankvariant qw{ analysis_rankvariant };
    use MIP::Recipes::Analysis::Rtg_vcfeval qw{ analysis_rtg_vcfeval  };
    use MIP::Recipes::Analysis::Sambamba_depth qw{ analysis_sambamba_depth };
    use MIP::Recipes::Analysis::Samtools_merge qw{ analysis_samtools_merge_panel };
    use MIP::Recipes::Analysis::Variant_annotation qw{ analysis_variant_annotation_panel };
    use MIP::Recipes::Analysis::Vep qw{ analysis_vep };
    use MIP::Recipes::Build::Rd_dna qw{ build_rd_dna_meta_files };

    ### Pipeline specific checks
    parse_rd_dna_panel(
        {
            active_parameter_href => $active_parameter_href,
            broadcasts_ref        => $broadcasts_ref,
            file_info_href        => $file_info_href,
            order_parameters_ref  => $order_parameters_ref,
            parameter_href        => $parameter_href,
            sample_info_href      => $sample_info_href,
        }
    );

    ### Build recipes
    $log->info(q{[Reference check - Reference prerequisites]});

    build_rd_dna_meta_files(
        {
            active_parameter_href => $active_parameter_href,
            file_info_href        => $file_info_href,
            job_id_href           => $job_id_href,
            log                   => $log,
            parameter_href        => $parameter_href,
            sample_info_href      => $sample_info_href,
        }
    );

    ## Check if references needs preprocessing
    ## If not try to reprocesses them before launching recipes
    parse_references(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            parameter_href        => $parameter_href,
        }
    );

    ### Analysis recipes
    ## Create code reference table for pipeline analysis recipes
    my %analysis_recipe = (
        analysisrunstatus                  => \&analysis_analysisrunstatus,
        bcftools_norm                      => \&analysis_bcftools_norm_panel,
        bwa_mem                            => undef,
        bwa_mem2                           => \&analysis_bwa_mem2,
        cadd_ar                            => \&analysis_cadd_panel,
        endvariantannotationblock          => \&analysis_endvariantannotationblock_panel,
        fastqc_ar                          => \&analysis_fastqc,
        frequency_filter                   => \&analysis_frequency_filter_panel,
        gatk_baserecalibration             => \&analysis_gatk_baserecalibration_panel,
        gatk_combinevariantcallsets        => \&analysis_gatk_combinevariantcallsets,
        gatk_haplotypecaller               => \&analysis_gatk_haplotypecaller_panel,
        gatk_gathervcfs                    => \&analysis_gatk_gathervcfs,
        gatk_genotypegvcfs                 => \&analysis_gatk_genotypegvcfs,
        gatk_variantevalall                => \&analysis_gatk_variantevalall,
        gatk_variantrecalibration          => \&analysis_gatk_variantrecalibration_wes,
        gzip_fastq                         => \&analysis_gzip_fastq,
        markduplicates                     => \&analysis_markduplicates_panel,
        multiqc_ar                         => \&analysis_multiqc,
        picardtools_collecthsmetrics       => \&analysis_picardtools_collecthsmetrics,
        picardtools_collectmultiplemetrics => \&analysis_picardtools_collectmultiplemetrics,
        qccollect_ar                       => \&analysis_mip_qccollect,
        rankvariant                        => \&analysis_rankvariant,
        rtg_vcfeval                        => \&analysis_rtg_vcfeval,
        sambamba_depth                     => \&analysis_sambamba_depth,
        samtools_merge                     => \&analysis_samtools_merge_panel,
        variant_annotation                 => \&analysis_variant_annotation_panel,
        varianteffectpredictor             => \&analysis_vep,
        vcfparser_ar                       => \&analysis_mip_vcfparser_panel,
        version_collect_ar                 => \&analysis_mip_vercollect,
    );

    ## Set correct bwa_mem recipe depending on version and source of the human_genome_reference: Source (hg19 or grch)
    set_recipe_bwa_mem(
        {
            analysis_recipe_href           => \%analysis_recipe,
            human_genome_reference_version => $file_info_href->{human_genome_reference_version},
            run_bwakit                     => $active_parameter_href->{bwa_mem_run_bwakit},
        }
    );

  RECIPE:
    foreach my $recipe ( @{$order_recipes_ref} ) {

        ## Skip not active recipes
        next RECIPE if ( not $active_parameter_href->{$recipe} );

        ## Skip recipe if not part of dispatch table (such as gzip_fastq)
        next RECIPE if ( not $analysis_recipe{$recipe} );

        ### Analysis recipes
        ## For displaying
        log_display_recipe_for_user(
            {
                log    => $log,
                recipe => $recipe,
            }
        );
        ## Sample mode
        if ( $parameter_href->{$recipe}{analysis_mode} eq q{sample} ) {

          SAMPLE_ID:
            foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

                $analysis_recipe{$recipe}->(
                    {
                        active_parameter_href => $active_parameter_href,
                        file_info_href        => $file_info_href,
                        job_id_href           => $job_id_href,
                        parameter_href        => $parameter_href,
                        recipe_name           => $recipe,
                        sample_id             => $sample_id,
                        sample_info_href      => $sample_info_href,
                    }
                );
            }
        }

        ## Family mode
        elsif ( $parameter_href->{$recipe}{analysis_mode} eq q{case} ) {

            $analysis_recipe{$recipe}->(
                {
                    active_parameter_href => $active_parameter_href,
                    file_info_href        => $file_info_href,
                    job_id_href           => $job_id_href,
                    parameter_href        => $parameter_href,
                    recipe_name           => $recipe,
                    sample_info_href      => $sample_info_href,
                }
            );
        }

    }
    return;
}

1;
