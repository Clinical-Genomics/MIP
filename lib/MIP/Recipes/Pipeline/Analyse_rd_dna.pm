package MIP::Recipes::Pipeline::Analyse_rd_dna;

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
use MIP::Constants qw{ $SPACE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.39;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ parse_rd_dna pipeline_analyse_rd_dna };
}

sub parse_rd_dna {

## Function : Rare disease DNA pipeline specific checks and parsing
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
    use MIP::Analysis qw{
      broadcast_parameters
      parse_prioritize_variant_callers
      update_prioritize_flag
      update_recipe_mode_for_wes
    };
    use MIP::Config qw{ write_mip_config };
    use MIP::Constants qw{ set_container_constants };
    use MIP::Contigs qw{ update_contigs_for_run };
    use MIP::Environment::Container qw{ parse_containers };
    use MIP::Fastq qw{ parse_fastq_infiles };
    use MIP::File_info qw{ check_parameter_metafiles parse_select_file_contigs };
    use MIP::Gatk qw{ check_gatk_sample_map_paths };
    use MIP::Parameter qw{ get_cache };
    use MIP::Parse::Gender qw{ parse_fastq_for_gender };
    use MIP::Reference
      qw{ get_select_file_contigs parse_exome_target_bed parse_nist_parameters };
    use MIP::Sample_info qw{ set_parameter_in_sample_info };
    use MIP::Vep qw{
      check_vep_api_cache_versions
      check_vep_custom_annotation
    };
    use MIP::Vcfanno qw{ parse_toml_config_parameters };

    ## Constants
    Readonly my @MIP_VEP_PLUGINS => qw{ sv_vep_plugin vep_plugin };
    Readonly my @ONLY_WGS_VARIANT_CALLER_RECIPES =>
      qw{ cnvnator_ar delly_reformat tiddit };
    Readonly my @ONLY_WGS_RECIPIES =>
      qw{ cnvnator_ar delly_call delly_reformat expansionhunter
      samtools_subsample_mt smncopynumbercaller star_caller telomerecat_ar tiddit };
    Readonly my @REMOVE_CONFIG_KEYS => qw{ associated_recipe };

    my $consensus_analysis_type = get_cache(
        {
            parameter_href => $parameter_href,
            parameter_name => q{consensus_analysis_type},
        }
    );

    ## Update exome_target_bed files with human_genome_reference_source and human_genome_reference_version
    parse_exome_target_bed(
        {
            exome_target_bed_file_href => $active_parameter_href->{exome_target_bed},
            human_genome_reference_source =>
              $file_info_href->{human_genome_reference_source},
            human_genome_reference_version =>
              $file_info_href->{human_genome_reference_version},
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

    ## Collect select file contigs to loop over downstream
    parse_select_file_contigs(
        {
            consensus_analysis_type => $consensus_analysis_type,
            file_info_href          => $file_info_href,
            select_file_path        => $active_parameter_href->{vcfparser_select_file},
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
            sample_map_path => $active_parameter_href->{gatk_genotypegvcfs_ref_gvcf},
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

    ## Check that all active variant callers have a prioritization order and that the prioritization elements match a
    ## supported variant caller
    parse_prioritize_variant_callers(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
        }
    );

    ## Update prioritize flag depending on analysis run value as some recipes are not applicable for e.g. wes
    $active_parameter_href->{sv_svdb_merge_prioritize} = update_prioritize_flag(
        {
            consensus_analysis_type => $consensus_analysis_type,
            parameter_href          => $parameter_href,
            prioritize_key          => $active_parameter_href->{sv_svdb_merge_prioritize},
            recipes_ref             => \@ONLY_WGS_VARIANT_CALLER_RECIPES,
        }
    );

    ## Update recipe mode depending on analysis run value as some recipes are not applicable for e.g. wes
    update_recipe_mode_for_wes(
        {
            active_parameter_href   => $active_parameter_href,
            consensus_analysis_type => $consensus_analysis_type,
            recipes_ref             => \@ONLY_WGS_RECIPIES,
        }
    );

    ## Write config file for case
    write_mip_config(
        {
            active_parameter_href => $active_parameter_href,
            remove_keys_ref       => \@REMOVE_CONFIG_KEYS,
            sample_info_href      => $sample_info_href,
        }
    );

    ## Update contigs depending on settings in run (wes or if only male samples)
    update_contigs_for_run(
        {
            consensus_analysis_type => $consensus_analysis_type,
            exclude_contigs_ref     => \@{ $active_parameter_href->{exclude_contigs} },
            file_info_href          => $file_info_href,
            include_y               => $active_parameter_href->{include_y},
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

    parse_fastq_for_gender(
        {
            active_parameter_href   => $active_parameter_href,
            consensus_analysis_type => $consensus_analysis_type,
            file_info_href          => $file_info_href,
            sample_info_href        => $sample_info_href,
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

    ## Set analysis constants
    set_container_constants( { active_parameter_href => $active_parameter_href, } );

    parse_containers(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
        }
    );

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
            mip_vep_plugins_ref   => \@MIP_VEP_PLUGINS,
        }
    );

    return;
}

sub pipeline_analyse_rd_dna {

## Function : Pipeline recipe for wes and or wgs data analysis.
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

    use MIP::Log::MIP_log4perl qw{ log_display_recipe_for_user };
    use MIP::Parse::Reference qw{ parse_references };
    use MIP::Set::Analysis
      qw{ set_recipe_bwa_mem set_recipe_gatk_variantrecalibration set_recipe_on_analysis_type set_rankvariants_ar };

    ## Recipes
    use MIP::Recipes::Analysis::Analysisrunstatus qw{ analysis_analysisrunstatus };
    use MIP::Recipes::Analysis::Cadd qw{ analysis_cadd };
    use MIP::Recipes::Analysis::Chanjo_sex_check qw{ analysis_chanjo_sex_check };
    use MIP::Recipes::Analysis::Chromograph
      qw{ analysis_chromograph_cov analysis_chromograph_upd };
    use MIP::Recipes::Analysis::Cnvnator qw{ analysis_cnvnator };
    use MIP::Recipes::Analysis::Deepvariant qw { analysis_deepvariant };
    use MIP::Recipes::Analysis::Delly_call qw{ analysis_delly_call };
    use MIP::Recipes::Analysis::Delly_reformat qw{ analysis_delly_reformat };
    use MIP::Recipes::Analysis::Endvariantannotationblock
      qw{ analysis_endvariantannotationblock };
    use MIP::Recipes::Analysis::Expansionhunter qw{ analysis_expansionhunter };
    use MIP::Recipes::Analysis::Fastqc qw{ analysis_fastqc };
    use MIP::Recipes::Analysis::Frequency_filter qw{ analysis_frequency_filter };
    use MIP::Recipes::Analysis::Gatk_baserecalibration
      qw{ analysis_gatk_baserecalibration };
    use MIP::Recipes::Analysis::Gatk_combinevariantcallsets
      qw{ analysis_gatk_combinevariantcallsets };
    use MIP::Recipes::Analysis::Gatk_gathervcfs qw{ analysis_gatk_gathervcfs };
    use MIP::Recipes::Analysis::Gatk_genotypegvcfs qw{ analysis_gatk_genotypegvcfs };
    use MIP::Recipes::Analysis::Gatk_haplotypecaller qw{ analysis_gatk_haplotypecaller };
    use MIP::Recipes::Analysis::Gatk_variantevalall qw{ analysis_gatk_variantevalall };
    use MIP::Recipes::Analysis::Gatk_variantevalexome
      qw{ analysis_gatk_variantevalexome };
    use MIP::Recipes::Analysis::Glnexus qw{ analysis_glnexus };
    use MIP::Recipes::Analysis::Gzip_fastq qw{ analysis_gzip_fastq };
    use MIP::Recipes::Analysis::Manta qw{ analysis_manta };
    use MIP::Recipes::Analysis::Markduplicates qw{ analysis_markduplicates };
    use MIP::Recipes::Analysis::Mip_qccollect qw{ analysis_mip_qccollect };
    use MIP::Recipes::Analysis::Mip_vcfparser qw{ analysis_mip_vcfparser };
    use MIP::Recipes::Analysis::Mip_vercollect qw{ analysis_mip_vercollect };
    use MIP::Recipes::Analysis::Multiqc qw{ analysis_multiqc };
    use MIP::Recipes::Analysis::Peddy qw{ analysis_peddy };
    use MIP::Recipes::Analysis::Picardtools_collecthsmetrics
      qw{ analysis_picardtools_collecthsmetrics };
    use MIP::Recipes::Analysis::Picardtools_collectmultiplemetrics
      qw{ analysis_picardtools_collectmultiplemetrics };
    use MIP::Recipes::Analysis::Plink qw{ analysis_plink };
    use MIP::Recipes::Analysis::Prepareforvariantannotationblock
      qw{ analysis_prepareforvariantannotationblock };
    use MIP::Recipes::Analysis::Rankvariant
      qw{ analysis_rankvariant analysis_rankvariant_unaffected analysis_rankvariant_sv analysis_rankvariant_sv_unaffected };
    use MIP::Recipes::Analysis::Rhocall
      qw{ analysis_rhocall_annotate analysis_rhocall_viz };
    use MIP::Recipes::Analysis::Rtg_vcfeval qw{ analysis_rtg_vcfeval  };
    use MIP::Recipes::Analysis::Sacct qw{ analysis_sacct };
    use MIP::Recipes::Analysis::Sambamba_depth qw{ analysis_sambamba_depth };
    use MIP::Recipes::Analysis::Samtools_merge qw{ analysis_samtools_merge };
    use MIP::Recipes::Analysis::Samtools_subsample_mt
      qw{ analysis_samtools_subsample_mt };
    use MIP::Recipes::Analysis::Smncopynumbercaller qw{ analysis_smncopynumbercaller };
    use MIP::Recipes::Analysis::Split_fastq_file qw{ analysis_split_fastq_file };
    use MIP::Recipes::Analysis::Star_caller qw{ analysis_star_caller };
    use MIP::Recipes::Analysis::Sv_annotate qw{ analysis_sv_annotate };
    use MIP::Recipes::Analysis::Sv_reformat qw{ analysis_reformat_sv };
    use MIP::Recipes::Analysis::Sv_combinevariantcallsets
      qw{ analysis_sv_combinevariantcallsets };
    use MIP::Recipes::Analysis::Split_fastq_file qw{ analysis_split_fastq_file };
    use MIP::Recipes::Analysis::Telomerecat qw{ analysis_telomerecat };
    use MIP::Recipes::Analysis::Tiddit qw{ analysis_tiddit };
    use MIP::Recipes::Analysis::Tiddit_coverage qw{ analysis_tiddit_coverage };
    use MIP::Recipes::Analysis::Upd qw{ analysis_upd };
    use MIP::Recipes::Analysis::Varg qw{ analysis_varg };
    use MIP::Recipes::Analysis::Variant_annotation qw{ analysis_variant_annotation };
    use MIP::Recipes::Analysis::Vcf2cytosure qw{ analysis_vcf2cytosure };
    use MIP::Recipes::Analysis::Vep qw{ analysis_vep_wgs };
    use MIP::Recipes::Analysis::Vt qw{ analysis_vt };
    use MIP::Recipes::Build::Rd_dna qw{build_rd_dna_meta_files};

    ### Pipeline specific checks
    parse_rd_dna(
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
        analysisrunstatus => \&analysis_analysisrunstatus,
        bwa_mem           => undef,                           # Depends on genome build
        bwa_mem2          => undef,
        cadd_ar           => \&analysis_cadd,
        chanjo_sexcheck   => \&analysis_chanjo_sex_check,
        chromograph_cov   => \&analysis_chromograph_cov,
        chromograph_upd   => \$sample_info_href->{has_trio}
        ? \&analysis_chromograph_upd
        : undef,                                              # Depends on pedigree
        cnvnator_ar                 => \&analysis_cnvnator,
        deepvariant                 => \&analysis_deepvariant,
        delly_call                  => \&analysis_delly_call,
        delly_reformat              => \&analysis_delly_reformat,
        endvariantannotationblock   => \&analysis_endvariantannotationblock,
        expansionhunter             => \&analysis_expansionhunter,
        fastqc_ar                   => \&analysis_fastqc,
        frequency_filter            => \&analysis_frequency_filter,
        gatk_baserecalibration      => \&analysis_gatk_baserecalibration,
        gatk_gathervcfs             => \&analysis_gatk_gathervcfs,
        gatk_combinevariantcallsets => \&analysis_gatk_combinevariantcallsets,
        gatk_genotypegvcfs          => \&analysis_gatk_genotypegvcfs,
        gatk_haplotypecaller        => \&analysis_gatk_haplotypecaller,
        gatk_variantevalall         => \&analysis_gatk_variantevalall,
        gatk_variantevalexome       => \&analysis_gatk_variantevalexome,
        gatk_variantrecalibration =>
          undef,    # Depends on analysis type and/or number of samples
        glnexus_merge                => \&analysis_glnexus,
        gzip_fastq                   => \&analysis_gzip_fastq,
        manta                        => \&analysis_manta,
        markduplicates               => \&analysis_markduplicates,
        multiqc_ar                   => \&analysis_multiqc,
        peddy_ar                     => \&analysis_peddy,
        picardtools_collecthsmetrics => \&analysis_picardtools_collecthsmetrics,
        picardtools_collectmultiplemetrics =>
          \&analysis_picardtools_collectmultiplemetrics,
        plink                            => \&analysis_plink,
        prepareforvariantannotationblock => \&analysis_prepareforvariantannotationblock,
        qccollect_ar                     => \&analysis_mip_qccollect,
        rankvariant    => undef,                         # Depends on sample features
        rhocall_ar     => \&analysis_rhocall_annotate,
        rhocall_viz    => \&analysis_rhocall_viz,
        rtg_vcfeval    => \&analysis_rtg_vcfeval,
        sacct          => \&analysis_sacct,
        sambamba_depth => \&analysis_sambamba_depth,
        samtools_merge => \&analysis_samtools_merge,
        samtools_subsample_mt     => \&analysis_samtools_subsample_mt,
        smncopynumbercaller       => \&analysis_smncopynumbercaller,
        split_fastq_file          => \&analysis_split_fastq_file,
        star_caller               => \&analysis_star_caller,
        sv_annotate               => \&analysis_sv_annotate,
        sv_combinevariantcallsets => \&analysis_sv_combinevariantcallsets,
        sv_rankvariant            => undef,                   # Depends on sample features
        sv_reformat               => \&analysis_reformat_sv,
        sv_varianteffectpredictor => undef,                   # Depends on analysis type
        sv_vcfparser              => undef,                   # Depends on analysis type
        telomerecat_ar            => \&analysis_telomerecat,
        tiddit                    => \&analysis_tiddit,
        tiddit_coverage        => \&analysis_tiddit_coverage,
        upd_ar                 => $sample_info_href->{has_trio} ? \&analysis_upd : undef,
        varg_ar                => \&analysis_varg,
        varianteffectpredictor => \&analysis_vep_wgs,
        variant_annotation     => \&analysis_variant_annotation,
        version_collect_ar     => \&analysis_mip_vercollect,
        vcfparser_ar           => \&analysis_mip_vcfparser,
        vcf2cytosure_ar        => \&analysis_vcf2cytosure,
        vt_ar                  => \&analysis_vt,
    );

    ## Special case for rankvariants recipe
    set_rankvariants_ar(
        {
            analysis_recipe_href => \%analysis_recipe,
            log                  => $log,
            parameter_href       => $parameter_href,
            sample_ids_ref       => $active_parameter_href->{sample_ids},
        }
    );

    ## Update which recipe to use depending on consensus analysis type
    set_recipe_on_analysis_type(
        {
            analysis_recipe_href    => \%analysis_recipe,
            consensus_analysis_type => $parameter_href->{cache}{consensus_analysis_type},
        }
    );

    ## Set correct bwa_mem recipe depending on version and source of the human_genome_reference: Source (hg19 or grch)
    set_recipe_bwa_mem(
        {
            analysis_recipe_href => \%analysis_recipe,
            human_genome_reference_version =>
              $file_info_href->{human_genome_reference_version},
        }
    );

    ## Update which recipe to use depending on number of samples
    set_recipe_gatk_variantrecalibration(
        {
            analysis_recipe_href => \%analysis_recipe,
            log                  => $log,
            sample_ids_ref       => $active_parameter_href->{sample_ids},
            use_cnnscorevariants => $active_parameter_href->{gatk_cnnscorevariants},
        }
    );

  RECIPE:
    foreach my $recipe ( @{$order_recipes_ref} ) {

        ## Skip not active recipes
        next RECIPE if ( not $active_parameter_href->{$recipe} );

        ## Skip recipe if not part of dispatch table (such as gzip_fastq)
        next RECIPE if ( not $analysis_recipe{$recipe} );

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

        ## Special case
        exit if ( $recipe eq q{split_fastq_file} );
    }
    return;
}

1;
