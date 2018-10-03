package MIP::Recipes::Pipeline::Vcf_rerun;

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

##MIPs lib/
use MIP::Delete::List qw{ delete_male_contig };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ pipeline_vcf_rerun };
}

## Constants
Readonly my $CLOSE_BRACKET => q{]};
Readonly my $OPEN_BRACKET  => q{[};
Readonly my $SPACE         => q{ };

sub pipeline_vcf_rerun {

## Function : Pipeline recipe for wes and or wgs data analysis.
## Returns  :
## Arguments: $active_parameter_href           => Active parameters for this analysis hash {REF}
##          : $broadcasts_ref                  => Holds the parameters info for broadcasting later {REF}
##          : $file_info_href                  => File info hash {REF}
##          : $infile_lane_prefix_href         => Infile(s) without the ".ending" {REF}
##          : $job_id_href                     => Job id hash {REF}
##          : $log                             => Log object to write to
##          : $order_parameters_ref            => Order of parameters (for structured output) {REF}
##          : $order_programs_ref              => Order of programs
##          : $outaligner_dir                  => Outaligner dir used in the analysis
##          : $parameter_href                  => Parameter hash {REF}
##          : $sample_info_href                => Info on samples and family hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $broadcasts_ref;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $log;
    my $order_parameters_ref;
    my $order_programs_ref;
    my $parameter_href;
    my $sample_info_href;

    ## Default(s)
    my $outaligner_dir;

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
        order_programs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$order_programs_ref,
            strict_type => 1,
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            store       => \$outaligner_dir,
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

    use MIP::Check::Pipeline qw{ check_vcf_rerun };

    ## Recipes
    use MIP::Log::MIP_log4perl qw{ log_display_program_for_user };
    use MIP::Recipes::Analysis::Analysisrunstatus
      qw{ analysis_analysisrunstatus };
    use MIP::Recipes::Analysis::Endvariantannotationblock
      qw{ analysis_endvariantannotationblock analysis_endvariantannotationblock_rio };
    use MIP::Recipes::Analysis::Frequency_filter
      qw{ analysis_frequency_filter analysis_frequency_filter_rio };
    use MIP::Recipes::Analysis::Mip_vcfparser
      qw{ analysis_mip_vcfparser analysis_mip_vcfparser_rio analysis_vcfparser_sv_wes analysis_vcfparser_sv_wgs };
    use MIP::Recipes::Analysis::Prepareforvariantannotationblock
      qw{ analysis_prepareforvariantannotationblock analysis_prepareforvariantannotationblock_rio };
    use MIP::Recipes::Analysis::Rankvariant
      qw{ analysis_rankvariant analysis_rankvariant_rio analysis_rankvariant_rio_unaffected analysis_rankvariant_unaffected analysis_rankvariant_sv analysis_rankvariant_sv_unaffected };
    use MIP::Recipes::Analysis::Rhocall
      qw{ analysis_rhocall_annotate analysis_rhocall_annotate_rio };
    use MIP::Recipes::Analysis::Sacct qw{ analysis_sacct };
    use MIP::Recipes::Analysis::Snpeff
      qw{ analysis_snpeff analysis_snpeff_rio };
    use MIP::Recipes::Analysis::Sv_annotate qw{ analysis_sv_annotate };
    use MIP::Recipes::Analysis::Sv_reformat qw{ analysis_reformat_sv };
    use MIP::Recipes::Analysis::Variantannotationblock
      qw{ analysis_variantannotationblock };
    use MIP::Recipes::Analysis::Vcf_rerun_reformat
      qw{ analysis_sv_vcf_rerun_reformat analysis_vcf_rerun_reformat };
    use MIP::Recipes::Analysis::Vep
      qw{ analysis_vep analysis_vep_rio analysis_vep_sv_wes analysis_vep_sv_wgs };
    use MIP::Recipes::Analysis::Vt qw{ analysis_vt analysis_vt_rio };
    use MIP::Recipes::Build::Human_genome_prerequisites
      qw{ build_human_genome_prerequisites };
    use MIP::Recipes::Build::Vcf_rerun qw{build_vcf_rerun_meta_files};
    use MIP::Set::Analysis qw{ set_recipe_on_analysis_type };

    ### Pipeline specific checks
    check_vcf_rerun(
        {
            active_parameter_href   => $active_parameter_href,
            broadcasts_ref          => $broadcasts_ref,
            file_info_href          => $file_info_href,
            infile_lane_prefix_href => $infile_lane_prefix_href,
            log                     => $log,
            order_parameters_ref    => $order_parameters_ref,
            parameter_href          => $parameter_href,
            sample_info_href        => $sample_info_href,
        }
    );

    ### Build recipes
    $log->info(q{[Reference check - Reference prerequisites]});

    build_vcf_rerun_meta_files(
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

    ### Analysis recipes
    ## Create code reference table for pipeline analysis recipes
    my %analysis_recipe = (
        analysisrunstatus         => \&analysis_analysisrunstatus,
        endvariantannotationblock => \&analysis_endvariantannotationblock,
        frequency_filter          => \&analysis_frequency_filter,
        prepareforvariantannotationblock =>
          \&analysis_prepareforvariantannotationblock,
        rankvariant => undef,                       # Depends on sample features
        rhocall     => \&analysis_rhocall_annotate,
        sacct       => \&analysis_sacct,
        snpeff      => \&analysis_snpeff,
        sv_annotate => \&analysis_sv_annotate,
        sv_rankvariant => undef,                    # Depends on sample features
        sv_reformat    => \&analysis_reformat_sv,
        sv_vcf_rerun_reformat => \&analysis_sv_vcf_rerun_reformat,
        sv_varianteffectpredictor => undef,          # Depends on analysis type,
        sv_vcfparser              => undef,          # Depends on analysis type
        varianteffectpredictor    => \&analysis_vep,
        vcfparser          => \&analysis_mip_vcfparser,
        vcf_rerun_reformat => \&analysis_vcf_rerun_reformat,
        vt                 => \&analysis_vt,
    );

    ### Special case for '--rio' capable analysis recipes
    ## Define rio blocks programs and order
    my $is_variantannotationblock_done;
    my @order_varann_programs;
    my %varann_ar;
    _define_variantannotationblock_ar(
        {
            active_parameter_href     => $active_parameter_href,
            order_varann_programs_ref => \@order_varann_programs,
            varann_ar_href            => \%varann_ar,
        }
    );

    ## Special case for rankvariants recipe
    _update_rankvariants_ar(
        {
            active_parameter_href => $active_parameter_href,
            analysis_recipe_href  => \%analysis_recipe,
            log                   => $log,
            parameter_href        => $parameter_href,
            varann_ar_href        => \%varann_ar,
        }
    );

    ## Update which recipe to use depending on consensus analysis type
    set_recipe_on_analysis_type(
        {
            analysis_recipe_href => \%analysis_recipe,
            consensus_analysis_type =>
              $parameter_href->{dynamic_parameter}{consensus_analysis_type},
        }
    );

  PROGRAM:
    foreach my $program ( @{$order_programs_ref} ) {

        ## Skip not active programs
        next PROGRAM if ( not $active_parameter_href->{$program} );

        ## Skip program if not part of dispatch table (such as gzip_fastq)
        next PROGRAM if ( not $analysis_recipe{$program} );

        ## Skip program if variant annotation block is done
        ## and program is part of variantannotation  block
        next PROGRAM
          if ( $is_variantannotationblock_done
            and any { $_ eq $program } @order_varann_programs );

        if ( $active_parameter_href->{reduce_io}
            and any { $_ eq $program } @order_varann_programs )
        {
            ## rio enabled and variantannotation block analysis recipe

            ## For displaying
            log_display_program_for_user(
                {
                    log     => $log,
                    program => q{variantannotationblock},
                }
            );

            analysis_variantannotationblock(
                {
                    active_parameter_href   => $active_parameter_href,
                    call_type               => q{BOTH},
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    outaligner_dir => $active_parameter_href->{outaligner_dir},
                    log            => $log,
                    order_programs_ref => \@order_varann_programs,
                    parameter_href     => $parameter_href,
                    program_name       => q{variantannotationblock},
                    sample_info_href   => $sample_info_href,
                    varann_ar_href     => \%varann_ar,
                }
            );

            ## Done with variantannotationblock block
            $is_variantannotationblock_done = 1;
        }
        else {

            ## For displaying
            log_display_program_for_user(
                {
                    log     => $log,
                    program => $program,
                }
            );

            ## Sample mode
            if ( $parameter_href->{$program}{analysis_mode} eq q{sample} ) {

              SAMPLE_ID:
                foreach
                  my $sample_id ( @{ $active_parameter_href->{sample_ids} } )
                {

                    $analysis_recipe{$program}->(
                        {
                            active_parameter_href   => $active_parameter_href,
                            file_info_href          => $file_info_href,
                            infile_lane_prefix_href => $infile_lane_prefix_href,
                            job_id_href             => $job_id_href,
                            parameter_href          => $parameter_href,
                            program_name            => $program,
                            sample_id               => $sample_id,
                            sample_info_href        => $sample_info_href,
                        }
                    );
                }
            }

            ## Family mode
            elsif ( $parameter_href->{$program}{analysis_mode} eq q{family} ) {

                $analysis_recipe{$program}->(
                    {
                        active_parameter_href   => $active_parameter_href,
                        file_info_href          => $file_info_href,
                        infile_lane_prefix_href => $infile_lane_prefix_href,
                        job_id_href             => $job_id_href,
                        parameter_href          => $parameter_href,
                        program_name            => $program,
                        sample_info_href        => $sample_info_href,
                    }
                );
            }
        }
    }
    return;
}

sub _define_variantannotationblock_ar {

## Function : Define variantannotationblock recipes, order, coderefs and activate
## Returns  :
## Arguments: $active_parameter_href     => Active parameters for this analysis hash {REF}
##          : $order_varann_programs_ref => Order of programs in variant annotation block {REF}
##          : $varann_ar_href            => Variant annotation analysis recipe hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $order_varann_programs_ref;
    my $varann_ar_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        order_varann_programs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$order_varann_programs_ref,
            strict_type => 1,
        },
        varann_ar_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$varann_ar_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Define rio blocks programs and order
    @{$order_varann_programs_ref} =
      qw{ prepareforvariantannotationblock rhocall vt frequency_filter varianteffectpredictor vcfparser snpeff rankvariant endvariantannotationblock };

    %{$varann_ar_href} = (
        endvariantannotationblock => \&analysis_endvariantannotationblock_rio,
        frequency_filter          => \&analysis_frequency_filter_rio,
        prepareforvariantannotationblock =>
          \&analysis_prepareforvariantannotationblock_rio,
        rankvariant => undef,    # Depends on sample features
        rhocall                => \&analysis_rhocall_annotate_rio,
        snpeff                 => \&analysis_snpeff_rio,
        varianteffectpredictor => \&analysis_vep_rio,
        vcfparser              => \&analysis_mip_vcfparser_rio,
        vt                     => \&analysis_vt_rio,
    );

    ## Enable varann as analysis recipe
    $active_parameter_href->{variantannotationblock} = 1;

    if ( $active_parameter_href->{dry_run_all} ) {

        ## Dry run
        $active_parameter_href->{variantannotationblock} = 2;
    }
    return;
}

sub _update_rankvariants_ar {

## Function : Update which rankvariants recipe to use
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $analysis_recipe_href    => Analysis recipe hash {REF}
##          : $log                     => Log object to write to
##          : $parameter_href          => Parameter hash {REF}
##          : $varann_ar_href            => Variant annotation analysis recipe hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $analysis_recipe_href;
    my $log;
    my $parameter_href;
    my $varann_ar_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        analysis_recipe_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$analysis_recipe_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        varann_ar_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$varann_ar_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    if ( defined $parameter_href->{dynamic_parameter}{unaffected}
        && @{ $parameter_href->{dynamic_parameter}{unaffected} } eq
        @{ $active_parameter_href->{sample_ids} } )
    {

        $log->warn(
q{Only unaffected sample(s) in pedigree - skipping genmod 'models', 'score' and 'compound'}
        );

        $analysis_recipe_href->{sv_rankvariant} =
          \&analysis_rankvariant_sv_unaffected;

        ## Rio recipe
        if ( $active_parameter_href->{reduce_io} ) {

            $varann_ar_href->{rankvariant} =
              \&analysis_rankvariant_rio_unaffected;
            return;
        }
        $analysis_recipe_href->{rankvariant} =
          \&analysis_rankvariant_unaffected;
    }
    else {

        $analysis_recipe_href->{sv_rankvariant} = \&analysis_rankvariant_sv;

        ## Rio recipe
        if ( $active_parameter_href->{reduce_io} ) {

            $varann_ar_href->{rankvariant} = \&analysis_rankvariant_rio;
            return;
        }
        $analysis_recipe_href->{rankvariant} = \&analysis_rankvariant;
    }
    return;
}

1;
