package MIP::Recipes::Pipeline::Rna;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use Readonly;

## MIPs lib/

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.10;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ pipeline_rna };
}

## Constants
Readonly my $CLOSE_BRACKET => q{]};
Readonly my $OPEN_BRACKET  => q{[};
Readonly my $SPACE         => q{ };

sub pipeline_rna {

## Function : Pipeline recipe for rna data analysis.
## Returns  :

## Arguments: $active_parameter_href           => Active parameters for this analysis hash {REF}
##          : $broadcasts_ref                  => Holds the parameters info for broadcasting later {REF}
##          : $file_info_href                  => File info hash {REF}
##          : $infile_both_strands_prefix_href => The infile(s) without the ".ending" and strand info {REF}
##          : $infile_lane_prefix_href         => Infile(s) without the ".ending" {REF}
##          : $job_id_href                     => Job id hash {REF}
##          : $log                             => Log object to write to
##          : $order_parameters_ref            => Order of parameters (for structured output) {REF}
##          : $order_programs_ref              => Execution order of programs {REF}
##          : $parameter_href                  => Parameter hash {REF}
##          : $sample_info_href                => Info on samples and family hash {REF}

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
    my $order_programs_ref;
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
            defined     => 1,
            required    => 1,
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
        order_programs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$order_programs_ref,
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

    use MIP::Check::Pipeline qw{ check_rna };

    ## Recipes
    use MIP::Recipes::Analysis::BootstrapAnn qw{ analysis_bootstrapann };
    use MIP::Recipes::Analysis::Fastqc qw{ analysis_fastqc };
    use MIP::Recipes::Analysis::Gatk_asereadcounter
      qw{ analysis_gatk_asereadcounter };
    use MIP::Recipes::Analysis::Gatk_baserecalibration
      qw{ analysis_gatk_baserecalibration };
    use MIP::Recipes::Analysis::Gatk_haplotypecaller
      qw{ analysis_gatk_haplotypecaller_rna };
    use MIP::Recipes::Analysis::Gatk_splitncigarreads
      qw{ analysis_gatk_splitncigarreads };
    use MIP::Recipes::Analysis::Gatk_variantfiltration
      qw{ analysis_gatk_variantfiltration };
    use MIP::Recipes::Analysis::Gzip_fastq qw{ analysis_gzip_fastq };
    use MIP::Recipes::Analysis::Multiqc qw{ analysis_multiqc };
    use MIP::Recipes::Analysis::Markduplicates qw{ analysis_markduplicates };
    use MIP::Recipes::Analysis::Picardtools_mergesamfiles
      qw{ analysis_picardtools_mergesamfiles };
    use MIP::Recipes::Analysis::Rseqc qw{ analysis_rseqc };
    use MIP::Recipes::Analysis::Sacct qw{ analysis_sacct };
    use MIP::Recipes::Analysis::Salmon_quant qw{ analysis_salmon_quant };
    use MIP::Recipes::Analysis::Star_aln qw{ analysis_star_aln };
    use MIP::Recipes::Analysis::Star_fusion qw{ analysis_star_fusion };
    use MIP::Recipes::Build::Rna qw{build_rna_meta_files};

    ### Pipeline specific checks
    check_rna(
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

    ### Build recipes
    $log->info(q{[Reference check - Reference prerequisites]});

    build_rna_meta_files(
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

    ## Dispatch table
    my %analysis_recipe = (
        bootstrapann              => \&analysis_bootstrapann,
        fastqc                    => \&analysis_fastqc,
        gatk_asereadcounter       => \&analysis_gatk_asereadcounter,
        gatk_baserecalibration    => \&analysis_gatk_baserecalibration,
        gatk_haplotypecaller      => \&analysis_gatk_haplotypecaller_rna,
        gatk_splitncigarreads     => \&analysis_gatk_splitncigarreads,
        gatk_variantfiltration    => \&analysis_gatk_variantfiltration,
        markduplicates            => \&analysis_markduplicates,
        multiqc                   => \&analysis_multiqc,
        picardtools_mergesamfiles => \&analysis_picardtools_mergesamfiles,
        rseqc                     => \&analysis_rseqc,
        sacct                     => \&analysis_sacct,
        salmon_quant              => \&analysis_salmon_quant,
        star_aln                  => \&analysis_star_aln,
        star_fusion               => \&analysis_star_fusion,
    );

    ## Program names for the log
    my %program_name = (
        bootstrapann              => q{BootstrapAnn},
        fastqc                    => q{FastQC},
        gatk_asereadcounter       => q{GATK ASEReadCounter},
        gatk_baserecalibration    => q{GATK BaseRecalibrator/ApplyBQSR},
        gatk_haplotypecaller      => q{GATK Haplotypecaller},
        gatk_splitncigarreads     => q{GATK SplitNCigarReads},
        gatk_variantfiltration    => q{GATK VariantFiltration},
        multiqc                   => q{Multiqc},
        markduplicates            => q{Markduplicates},
        picardtools_mergesamfiles => q{Picardtools MergeSamFiles},
        rseqc                     => q{Rseqc},
        sacct                     => q{Sacct},
        salmon_quant              => q{Salmon Quant},
        star_aln                  => q{STAR},
        star_fusion               => q{STAR-Fusion},
    );

  PROGRAM:
    foreach my $program ( @{$order_programs_ref} ) {

        ## Skip not active programs
        next PROGRAM if ( not $active_parameter_href->{$program} );

        ## Skip program if not part of dispatch table (such as gzip fastq)
        next PROGRAM if ( not $analysis_recipe{$program} );

        $log->info( $OPEN_BRACKET . $program_name{$program} . $CLOSE_BRACKET );

        ## Sample mode
        if ( $parameter_href->{$program}{analysis_mode} eq q{sample} ) {

          SAMPLE:
            foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } )
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

    return;
}

1;
