package MIP::Recipes::Pipeline::Cancer;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use File::Spec::Functions qw{ catdir catfile };

## CPANM
use Readonly;

##MIPs lib/
use MIP::Delete::List qw{ delete_male_contig };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ pipeline_cancer };
}

## Constants
Readonly my $SPACE => q{ };

sub pipeline_cancer {

## Function : Pipeline recipe for cancer data analysis.
## Returns  :

## Arguments: $parameter_href          => Parameter hash {REF}
##          : $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $file_info_href          => File info hash {REF}
##          : $indir_path_href         => Indirectory hash {REF}
##          : $infile_href             => Infile hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $lane_href               => The lane info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner dir used in the analysis
##          : $log                     => Log object to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $indir_path_href;
    my $infile_href;
    my $infile_lane_prefix_href;
    my $lane_href;
    my $job_id_href;
    my $log;

    ## Default(s)
    my $outaligner_dir;

    my $tmpl = {
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        indir_path_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$indir_path_href,
            strict_type => 1,
        },
        infile_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_href,
            strict_type => 1,
        },
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_lane_prefix_href,
            strict_type => 1,
        },
        lane_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$lane_href,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            store       => \$outaligner_dir,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Recipes
    use MIP::Recipes::Analysis::Analysisrunstatus
      qw{ analysis_analysisrunstatus };
    use MIP::Recipes::Analysis::Bwa_mem qw{ analysis_bwa_mem };
    use MIP::Recipes::Analysis::Chanjo_sex_check
      qw{ analysis_chanjo_sex_check };
    use MIP::Recipes::Analysis::Fastqc qw{ analysis_fastqc };
    use MIP::Recipes::Analysis::Gatk_baserecalibration
      qw{ analysis_gatk_baserecalibration };
    use MIP::Recipes::Analysis::Markduplicates qw{ analysis_markduplicates };
    use MIP::Recipes::Analysis::Multiqc qw{ analysis_multiqc };
    use MIP::Recipes::Analysis::Picardtools_collecthsmetrics
      qw{ analysis_picardtools_collecthsmetrics };
    use MIP::Recipes::Analysis::Picardtools_collectmultiplemetrics
      qw{ analysis_picardtools_collectmultiplemetrics };
    use MIP::Recipes::Analysis::Picardtools_mergesamfiles
      qw{ analysis_picardtools_mergesamfiles };
    use MIP::Recipes::Analysis::Qccollect qw{ analysis_qccollect };
    use MIP::Recipes::Analysis::Sacct qw{ analysis_sacct };
    use MIP::Recipes::Analysis::Sambamba_depth qw{ analysis_sambamba_depth };

## FastQC oer sample_id
    if ( $active_parameter_href->{pfastqc} ) {

        $log->info(q{[Fastqc]});

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            my $fastqc_program_name = q{fastqc};

            my $outsample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $fastqc_program_name );
            analysis_fastqc(
                {
                    parameter_href        => $parameter_href,
                    active_parameter_href => $active_parameter_href,
                    sample_info_href      => $sample_info_href,
                    infiles_ref           => \@{ $infile_href->{$sample_id} },

                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    insample_directory      => $indir_path_href->{$sample_id},
                    outsample_directory     => $outsample_directory,
                    sample_id               => $sample_id,
                    program_name            => $fastqc_program_name,
                }
            );
        }
    }

## Aligning fastq files based on sample_id
    if ( $active_parameter_href->{pbwa_mem} ) {

        $log->info(q{[BWA Mem]});

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            my $outsample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $active_parameter_href->{outaligner_dir} );

            analysis_bwa_mem(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    file_info_href          => $file_info_href,
                    infiles_ref             => \@{ $infile_href->{$sample_id} },
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    insample_directory      => $indir_path_href->{$sample_id},
                    outsample_directory     => $outsample_directory,
                    sample_id               => $sample_id,
                    program_name            => q{bwa_mem},
                }
            );
        }
    }

## Merge bam files
    if ( $active_parameter_href->{ppicardtools_mergesamfiles} ) {

        $log->info(q{[Picardtools mergesamfiles]});

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            my $insample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $active_parameter_href->{outaligner_dir} );
            my $outsample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $active_parameter_href->{outaligner_dir} );

            analysis_picardtools_mergesamfiles(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    lane_href               => $lane_href,
                    job_id_href             => $job_id_href,
                    insample_directory      => $insample_directory,
                    outsample_directory     => $outsample_directory,
                    sample_id               => $sample_id,
                    program_name            => q{picardtools_mergesamfiles},
                }
            );
        }
    }

## MarkDuplicates
    if ( $active_parameter_href->{pmarkduplicates} ) {

        $log->info(q{[Markduplicates]});

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            ## Assign directories
            my $insample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $active_parameter_href->{outaligner_dir} );
            my $outsample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $active_parameter_href->{outaligner_dir} );

            analysis_markduplicates(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    insample_directory      => $insample_directory,
                    outsample_directory     => $outsample_directory,
                    sample_id               => $sample_id,
                    program_name            => q{markduplicates},
                }
            );
        }
    }

## Base recalibrator
    if ( $active_parameter_href->{pgatk_baserecalibration} ) {

        $log->info(q{[GATK baserecalibrator/printreads]});

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            ## Assign directories
            my $insample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $active_parameter_href->{outaligner_dir} );
            my $outsample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $active_parameter_href->{outaligner_dir} );

            analysis_gatk_baserecalibration(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    sample_id               => $sample_id,
                    insample_directory      => $insample_directory,
                    outsample_directory     => $outsample_directory,
                    program_name            => q{gatk_baserecalibration},
                }
            );
        }
    }

## Chanjo sex check
    if ( $active_parameter_href->{pchanjo_sexcheck} ) {

        $log->info(q{[Chanjo sexcheck]});

      SAMPLE_IDS:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            my $insample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $active_parameter_href->{outaligner_dir} );
            my $outsample_directory = catdir(
                $active_parameter_href->{outdata_dir},    $sample_id,
                $active_parameter_href->{outaligner_dir}, q{coveragereport}
            );

            analysis_chanjo_sex_check(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    sample_id               => $sample_id,
                    insample_directory      => $insample_directory,
                    outsample_directory     => $outsample_directory,
                    program_name            => q{chanjo_sexcheck},
                }
            );
        }
    }

## Sambamba depth
    if ( $active_parameter_href->{psambamba_depth} ) {

        $log->info(q{[Sambamba depth]});

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            my $insample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $active_parameter_href->{outaligner_dir} );
            my $outsample_directory = catdir(
                $active_parameter_href->{outdata_dir},    $sample_id,
                $active_parameter_href->{outaligner_dir}, q{coveragereport}
            );

            analysis_sambamba_depth(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    sample_id               => $sample_id,
                    insample_directory      => $insample_directory,
                    outsample_directory     => $outsample_directory,
                    program_name            => q{sambamba_depth},
                }
            );
        }
    }

## Picard collect multiple metrics
    if ( $active_parameter_href->{ppicardtools_collectmultiplemetrics} ) {

        $log->info(q{[Picardtools collectmultiplemetrics]});

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            ## Assign directories
            my $insample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $active_parameter_href->{outaligner_dir} );
            my $outsample_directory = catdir(
                $active_parameter_href->{outdata_dir},    $sample_id,
                $active_parameter_href->{outaligner_dir}, q{coveragereport}
            );

            analysis_picardtools_collectmultiplemetrics(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    sample_id               => $sample_id,
                    insample_directory      => $insample_directory,
                    outsample_directory     => $outsample_directory,
                    program_name => q{picardtools_collectmultiplemetrics},
                }
            );
        }
    }

## Picard Collect HS metrics
    if ( $active_parameter_href->{ppicardtools_collecthsmetrics} ) {

        $log->info(q{[Picardtools collecthsmetrics]});

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            ## Assign directories
            my $insample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $active_parameter_href->{outaligner_dir} );
            my $outsample_directory = catdir(
                $active_parameter_href->{outdata_dir},    $sample_id,
                $active_parameter_href->{outaligner_dir}, q{coveragereport}
            );

            analysis_picardtools_collecthsmetrics(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    sample_id               => $sample_id,
                    insample_directory      => $insample_directory,
                    outsample_directory     => $outsample_directory,
                    program_name            => q{picardtools_collecthsmetrics},
                }
            );
        }
    }

## QC collect
    if ( $active_parameter_href->{pqccollect} ) {

        $log->info(q{[Qccollect]});

        my $outfamily_directory = catdir(
            $active_parameter_href->{outdata_dir},
            $active_parameter_href->{family_id}
        );

        analysis_qccollect(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                sample_info_href        => $sample_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                program_name            => q{qccollect},
                outfamily_directory     => $outfamily_directory,
            }
        );
    }

## MultiQC
    if ( $active_parameter_href->{pmultiqc} ) {

        $log->info(q{[Multiqc]});

        analysis_multiqc(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                sample_info_href        => $sample_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                program_name            => q{multiqc},
            }
        );
    }

    if ( $active_parameter_href->{panalysisrunstatus} ) {

        $log->info(q{[Analysis run status]});

        analysis_analysisrunstatus(
            {
                active_parameter_href   => $active_parameter_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                parameter_href          => $parameter_href,
                program_name            => q{analysisrunstatus},
                sample_info_href        => $sample_info_href,
            }
        );
    }

    if ( $active_parameter_href->{psacct} ) {

        $log->info(q{[Sacct]});

        analysis_sacct(
            {
                active_parameter_href   => $active_parameter_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                parameter_href          => $parameter_href,
                program_name            => q{sacct},
                sample_info_href        => $sample_info_href,
            }
        );
    }
    return;
}

1;
