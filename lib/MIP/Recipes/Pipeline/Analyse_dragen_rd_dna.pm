package MIP::Recipes::Pipeline::Analyse_dragen_rd_dna;

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

    our $VERSION = 1.09;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ pipeline_analyse_dragen_rd_dna };
}

sub pipeline_analyse_dragen_rd_dna {

## Function : Pipeline recipe for Dragen wgs data analysis
## Returns  :
## Arguments: $active_parameter_href           => Active parameters for this analysis hash {REF}
##          : $broadcasts_ref                  => Holds the parameters info for broadcasting later {REF}
##          : $file_info_href                  => File info hash {REF}
##          : $infile_both_strands_prefix_href => The infile(s) without the ".ending" and strand info {REF}
##          : $infile_lane_prefix_href         => Infile(s) without the ".ending" {REF}
##          : $job_id_href                     => Job id hash {REF}
##          : $log                             => Log object to write to
##          : $order_parameters_ref            => Order of parameters (for structured output) {REF}
##          : $order_recipes_ref               => Order of recipes
##          : $parameter_href                  => Parameter hash {REF}
##          : $sample_info_href                => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $broadcasts_ref;
    my $file_info_href;
    my $infile_both_strands_prefix_href;
    my $infile_lane_prefix_href;
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
        infile_both_strands_prefix_href => {
            default     => {},
            store       => \$infile_both_strands_prefix_href,
            strict_type => 1,
        },
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_lane_prefix_href,
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

    use MIP::Check::Pipeline qw{ check_dragen_rd_dna };
    use MIP::Constants qw{ set_singularity_constants };
    use MIP::Parse::Reference qw{ parse_reference_for_vt };
    use MIP::Set::Analysis qw{ set_recipe_on_analysis_type set_rankvariants_ar };

    ## Recipes
    use MIP::Log::MIP_log4perl qw{ log_display_recipe_for_user };
    use MIP::Recipes::Analysis::Analysisrunstatus qw{ analysis_analysisrunstatus };
    use MIP::Recipes::Analysis::Cadd qw{ analysis_cadd };
    use MIP::Recipes::Analysis::Dragen_dna
      qw{ analysis_dragen_dna_align_vc analysis_dragen_dna_joint_calling };
    use MIP::Recipes::Analysis::Endvariantannotationblock
      qw{ analysis_endvariantannotationblock };
    use MIP::Recipes::Analysis::Frequency_annotation qw{ analysis_frequency_annotation };
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
    use MIP::Recipes::Analysis::Vcf_rerun_reformat
      qw{ analysis_vcf_rerun_reformat_sv analysis_vcf_rerun_reformat };
    use MIP::Recipes::Analysis::Vep
      qw{ analysis_vep analysis_vep_sv_wes analysis_vep_sv_wgs };
    use MIP::Recipes::Analysis::Vt qw{ analysis_vt };
    use MIP::Recipes::Build::Human_genome_prerequisites
      qw{ build_human_genome_prerequisites };
    use MIP::Recipes::Build::Dragen_rd_dna qw{build_dragen_rd_dna_meta_files};

    ### Pipeline specific checks
    check_dragen_rd_dna(
        {
            active_parameter_href           => $active_parameter_href,
            broadcasts_ref                  => $broadcasts_ref,
            file_info_href                  => $file_info_href,
            infile_both_strands_prefix_href => $infile_both_strands_prefix_href,
            infile_lane_prefix_href         => $infile_lane_prefix_href,
            log                             => $log,
            order_parameters_ref            => $order_parameters_ref,
            parameter_href                  => $parameter_href,
            sample_info_href                => $sample_info_href,
        }
    );

    ## Set analysis constants
    set_singularity_constants( { active_parameter_href => $active_parameter_href, } );

    ### Build recipes
    $log->info(q{[Reference check - Reference prerequisites]});

    build_dragen_rd_dna_meta_files(
        {
            active_parameter_href   => $active_parameter_href,
            file_info_href          => $file_info_href,
            infile_lane_prefix_href => $infile_lane_prefix_href,
            job_id_href             => $job_id_href,
            log                     => $log,
            parameter_href          => $parameter_href,
            sample_info_href        => $sample_info_href,
        }
    );

    ## Check if vt has processed references
    ## If not try to reprocesses them before launching recipes
    $log->info(q{[Reference check - Reference processed by VT]});
    parse_reference_for_vt(
        {
            active_parameter_href   => $active_parameter_href,
            infile_lane_prefix_href => $infile_lane_prefix_href,
            job_id_href             => $job_id_href,
            log                     => $log,
            parameter_href          => $parameter_href,
        }
    );

    ### Analysis recipes
    ## Create code reference table for pipeline analysis recipes
    my %analysis_recipe = (
        analysisrunstatus                => \&analysis_analysisrunstatus,
        cadd_ar                          => \&analysis_cadd,
        dragen_dna_align_vc              => \&analysis_dragen_dna_align_vc,
        dragen_dna_joint_calling         => \&analysis_dragen_dna_joint_calling,
        endvariantannotationblock        => \&analysis_endvariantannotationblock,
        frequency_annotation             => \&analysis_frequency_annotation,
        frequency_filter                 => \&analysis_frequency_filter,
        prepareforvariantannotationblock => \&analysis_prepareforvariantannotationblock,
        rankvariant    => undef,                         # Depends on sample features
        rhocall_ar     => \&analysis_rhocall_annotate,
        sacct          => \&analysis_sacct,
        sv_annotate    => \&analysis_sv_annotate,
        sv_rankvariant => undef,                         # Depends on sample features
        sv_reformat    => \&analysis_reformat_sv,
        sv_vcf_rerun_reformat => \&analysis_vcf_rerun_reformat_sv,
        sv_varianteffectpredictor => undef,                    # Depends on analysis type,
        sv_vcfparser              => undef,                    # Depends on analysis type
        varianteffectpredictor    => \&analysis_vep,
        vcfparser_ar              => \&analysis_mip_vcfparser,
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

        ## Sample mode
        if ( $parameter_href->{$recipe}{analysis_mode} eq q{sample} ) {

          SAMPLE_ID:
            foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

                $analysis_recipe{$recipe}->(
                    {
                        active_parameter_href   => $active_parameter_href,
                        file_info_href          => $file_info_href,
                        infile_lane_prefix_href => $infile_lane_prefix_href,
                        job_id_href             => $job_id_href,
                        parameter_href          => $parameter_href,
                        recipe_name             => $recipe,
                        sample_id               => $sample_id,
                        sample_info_href        => $sample_info_href,
                    }
                );
            }
        }

        ## Family mode
        elsif ( $parameter_href->{$recipe}{analysis_mode} eq q{case} ) {

            $analysis_recipe{$recipe}->(
                {
                    active_parameter_href   => $active_parameter_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    parameter_href          => $parameter_href,
                    recipe_name             => $recipe,
                    sample_info_href        => $sample_info_href,
                }
            );
        }
    }
    return;
}

1;
