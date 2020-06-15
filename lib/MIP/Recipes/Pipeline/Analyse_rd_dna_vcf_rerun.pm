package MIP::Recipes::Pipeline::Analyse_rd_dna_vcf_rerun;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use List::MoreUtils qw { any };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $CLOSE_BRACKET $OPEN_BRACKET $SPACE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.17;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ parse_rd_dna_vcf_rerun pipeline_analyse_rd_dna_vcf_rerun };
}

sub parse_rd_dna_vcf_rerun {

## Function : Rare disease DNA vcf rerun pipeline specific checks and parsing
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $broadcasts_ref          => Holds the parameters info for broadcasting later {REF}
##          : $file_info_href          => File info hash {REF}
##          : $order_parameters_ref    => Order of parameters (for structured output) {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $sample_info_href        => Info on samples and case hash {REF}

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
      parse_vep_plugin
      set_vcfparser_outfile_counter
      write_references
    };
    use MIP::Analysis qw{ broadcast_parameters };
    use MIP::Config qw{ write_mip_config };
    use MIP::Contigs qw{ update_contigs_for_run };
    use MIP::File_info qw{ check_parameter_metafiles parse_select_file_contigs };
    use MIP::Parameter qw{ get_cache };
    use MIP::Reference qw{ get_select_file_contigs };
    use MIP::Sample_info qw{ set_parameter_in_sample_info };
    use MIP::Vep qw{
      check_vep_api_cache_versions
      check_vep_custom_annotation
    };

    ## Constants
    Readonly my @MIP_VEP_PLUGINS    => qw{ sv_vep_plugin vep_plugin };
    Readonly my @REMOVE_CONFIG_KEYS => qw{ associated_recipe };

    my $consensus_analysis_type = get_cache(
        {
            parameter_href => $parameter_href,
            parameter_name => q{consensus_analysis_type},
        }
    );

    ## Check sample_id provided in hash parameter is included in the analysis
    check_sample_id_in_hash_parameter(
        {
            active_parameter_href => $active_parameter_href,
            parameter_names_ref   => [qw{ analysis_type }],
            parameter_href        => $parameter_href,
            sample_ids_ref        => \@{ $active_parameter_href->{sample_ids} },
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

sub pipeline_analyse_rd_dna_vcf_rerun {

## Function : Pipeline recipe for wes and or wgs data analysis rerun
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

    use MIP::Constants qw{ set_singularity_constants };
    use MIP::Parse::Reference qw{ parse_references };
    use MIP::Set::Analysis qw{ set_recipe_on_analysis_type set_rankvariants_ar };

    ## Recipes
    use MIP::Log::MIP_log4perl qw{ log_display_recipe_for_user };
    use MIP::Recipes::Analysis::Analysisrunstatus qw{ analysis_analysisrunstatus };
    use MIP::Recipes::Analysis::Cadd qw{ analysis_cadd };
    use MIP::Recipes::Analysis::Endvariantannotationblock
      qw{ analysis_endvariantannotationblock };
    use MIP::Recipes::Analysis::Frequency_filter qw{ analysis_frequency_filter };
    use MIP::Recipes::Analysis::Mip_vcfparser
      qw{ analysis_mip_vcfparser analysis_mip_vcfparser_sv_wes analysis_mip_vcfparser_sv_wgs };
    use MIP::Recipes::Analysis::Mip_vercollect qw{ analysis_mip_vercollect };
    use MIP::Recipes::Analysis::Prepareforvariantannotationblock
      qw{ analysis_prepareforvariantannotationblock };
    use MIP::Recipes::Analysis::Rankvariant
      qw{ analysis_rankvariant analysis_rankvariant_unaffected analysis_rankvariant_sv analysis_rankvariant_sv_unaffected };
    use MIP::Recipes::Analysis::Rhocall qw{ analysis_rhocall_annotate };
    use MIP::Recipes::Analysis::Sacct qw{ analysis_sacct };
    use MIP::Recipes::Analysis::Sv_annotate qw{ analysis_sv_annotate };
    use MIP::Recipes::Analysis::Sv_reformat qw{ analysis_reformat_sv };
    use MIP::Recipes::Analysis::Variant_annotation qw{ analysis_variant_annotation };
    use MIP::Recipes::Analysis::Vcf_rerun_reformat
      qw{ analysis_vcf_rerun_reformat_sv analysis_vcf_rerun_reformat };
    use MIP::Recipes::Analysis::Vep
      qw{ analysis_vep_wgs analysis_vep_sv_wes analysis_vep_sv_wgs };
    use MIP::Recipes::Analysis::Vt qw{ analysis_vt };
    use MIP::Recipes::Build::Human_genome_prerequisites
      qw{ build_human_genome_prerequisites };
    use MIP::Recipes::Build::Rd_dna_vcf_rerun qw{build_rd_dna_vcf_rerun_meta_files};

    ### Pipeline specific checks
    parse_rd_dna_vcf_rerun(
        {
            active_parameter_href => $active_parameter_href,
            broadcasts_ref        => $broadcasts_ref,
            file_info_href        => $file_info_href,
            order_parameters_ref  => $order_parameters_ref,
            parameter_href        => $parameter_href,
            sample_info_href      => $sample_info_href,
        }
    );

    ## Set analysis constants
    set_singularity_constants( { active_parameter_href => $active_parameter_href, } );

    ### Build recipes
    $log->info(q{[Reference check - Reference prerequisites]});

    build_rd_dna_vcf_rerun_meta_files(
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
        analysisrunstatus                => \&analysis_analysisrunstatus,
        cadd_ar                          => \&analysis_cadd,
        endvariantannotationblock        => \&analysis_endvariantannotationblock,
        frequency_filter                 => \&analysis_frequency_filter,
        prepareforvariantannotationblock => \&analysis_prepareforvariantannotationblock,
        rankvariant    => undef,                         # Depends on sample features
        rhocall_ar     => \&analysis_rhocall_annotate,
        sacct          => \&analysis_sacct,
        sv_annotate    => \&analysis_sv_annotate,
        sv_rankvariant => undef,                         # Depends on sample features
        sv_reformat    => \&analysis_reformat_sv,
        sv_vcf_rerun_reformat => \&analysis_vcf_rerun_reformat_sv,
        sv_varianteffectpredictor => undef,                # Depends on analysis type,
        sv_vcfparser              => undef,                # Depends on analysis type
        varianteffectpredictor    => \&analysis_vep_wgs,
        variant_annotation => \&analysis_variant_annotation,
        vcfparser_ar       => \&analysis_mip_vcfparser,
        vcf_rerun_reformat => \&analysis_vcf_rerun_reformat,
        version_collect_ar => \&analysis_mip_vercollect,
        vt_ar              => \&analysis_vt,
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

        if ( $parameter_href->{$recipe}{analysis_mode} eq q{case} ) {

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
