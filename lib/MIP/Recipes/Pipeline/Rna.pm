package MIP::Recipes::Pipeline::Rna;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use File::Spec::Functions qw{ catdir };

## CPANM
use Readonly;

## MIPs lib/

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ pipeline_rna };
}

## Constants
Readonly my $SPACE => q{ };

sub pipeline_rna {

## Function : Pipeline recipe for rna data analysis.
## Returns  :

## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $file_info_href          => File info hash {REF}
##          : $indir_path_href         => Indirectory hash {REF}
##          : $infile_href             => Infile hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $lane_href               => The lane info hash {REF}
##          : $log                     => Log object to write to
##          : $outaligner_dir          => Outaligner dir used in the analysis
##          : $parameter_href          => Parameter hash {REF}
##          : $sample_info_href        => Info on samples and family hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $indir_path_href;
    my $infile_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $lane_href;
    my $log;
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
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        lane_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$lane_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
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

    use MIP::Recipes::Analysis::Fastqc qw{ analysis_fastqc };
    use MIP::Recipes::Analysis::Star_aln qw{ analysis_star_aln };
    use MIP::Recipes::Analysis::Star_fusion qw{ analysis_star_fusion };
    use MIP::Recipes::Build::Rna qw{build_rna_meta_files};

    ### Build recipes
    $log->info(q{[Reference check - Reference prerequisites]});

    build_rna_meta_files(
        {
            parameter_href          => $parameter_href,
            active_parameter_href   => $active_parameter_href,
            sample_info_href        => $sample_info_href,
            file_info_href          => $file_info_href,
            infile_lane_prefix_href => $infile_lane_prefix_href,
            job_id_href             => $job_id_href,
            log                     => $log,
        }
    );

    ### Analysis recipes
    ## Run FastQC
    if ( $active_parameter_href->{pfastqc} ) {

        $log->info(q{[Fastqc]});

      SAMPLE:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            my $fastqc_program_name = q{fastqc};
            my $fastq_outsample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $fastqc_program_name );
            analysis_fastqc(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    infiles_ref             => \@{ $infile_href->{$sample_id} },
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    insample_directory      => $indir_path_href->{$sample_id},
                    outsample_directory     => $fastq_outsample_directory,
                    sample_id               => $sample_id,
                    program_name            => $fastqc_program_name,
                }
            );
        }
    }

    ## Star aln
    if ( $active_parameter_href->{pstar_aln} ) {

        $log->info(q{[Star Aln]});

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            my $star_outsample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $active_parameter_href->{outaligner_dir} );

            analysis_star_aln(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    file_info_href          => $file_info_href,
                    infiles_ref             => \@{ $infile_href->{$sample_id} },
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    insample_directory      => $indir_path_href->{$sample_id},
                    outsample_directory     => $star_outsample_directory,
                    sample_id               => $sample_id,
                    program_name            => q{star_aln},
                }
            );
        }
    }
    ## Star fusion
    if ( $active_parameter_href->{pstar_fusion} ) {

        $log->info(q{[Star Fusion]});

      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            my $star_insample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $active_parameter_href->{outaligner_dir} );
            my $star_outsample_directory =
              catdir( $active_parameter_href->{outdata_dir},
                $sample_id, $active_parameter_href->{outaligner_dir} );

            analysis_star_fusion(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    insample_directory      => $star_insample_directory,
                    outsample_directory     => $star_outsample_directory,
                    sample_id               => $sample_id,
                    program_name            => q{star_fusion},
                }
            );
        }
    }
    return;
}

1;
