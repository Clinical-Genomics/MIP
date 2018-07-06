package MIP::Recipes::Analysis::Chanjo_sex_check;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile devnull };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use List::MoreUtils qw { any };
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.06;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_chanjo_sex_check };

}

##Constants
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub analysis_chanjo_sex_check {

## Function : Predicts gender from BAM files.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $file_info_href          => The file_info hash {REF}
##          : indir_path_href          => Indirectories path(s) hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $sample_id               => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $indir_path_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $program_name;
    my $sample_info_href;
    my $sample_id;

    ## Default(s)
    my $family_id;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            store       => \$family_id,
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
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        program_name => {
            defined     => 1,
            required    => 1,
            store       => \$program_name,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
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

    use MIP::Get::File qw{ get_file_suffix get_merged_infile_prefix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_dead_end };
    use MIP::Program::Alignment::Chanjo qw{ chanjo_sex };
    use MIP::QC::Record
      qw{ add_program_metafile_to_sample_info add_program_outfile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $program_mode = $active_parameter_href->{$program_name};

    ## Alias
    my $job_id_chain = $parameter_href->{$program_name}{chain};
    my ( $core_number, $time, @source_environment_cmds ) =
      get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
      );

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Assign directories
    my $insample_directory  = $indir_path_href->{$sample_id};
    my $outsample_directory = catdir(
        $active_parameter_href->{outdata_dir},    $sample_id,
        $active_parameter_href->{outaligner_dir}, q{coveragereport}
    );

    ## Add merged infile name prefix after merging all BAM files per sample_id
    my $merged_infile_prefix = get_merged_infile_prefix(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
    );

    ## Assign file_tags
    my $infile_tag =
      $file_info_href->{$sample_id}{gatk_baserecalibration}{file_tag};
    my $outfile_tag =
      $file_info_href->{$sample_id}{$program_name}{file_tag};

    my $infile_prefix  = $merged_infile_prefix . $infile_tag;
    my $outfile_prefix = $merged_infile_prefix . $outfile_tag;

    ## Get infile_suffix from baserecalibration jobid chain
    my $infile_suffix = get_file_suffix(
        {
            jobid_chain    => $parameter_href->{gatk_baserecalibration}{chain},
            parameter_href => $parameter_href,
            suffix_key     => q{alignment_file_suffix},
        }
    );
    my $outfile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            program_name   => $program_name,
            suffix_key     => q{outfile_suffix},
        }
    );

    ## Files
    my $infile_name  = $infile_prefix . $infile_suffix;
    my $outfile_name = $outfile_prefix . $outfile_suffix;

    ## Paths
    my $infile_path  = catfile( $insample_directory,  $infile_name );
    my $outfile_path = catfile( $outsample_directory, $outfile_name );
    my $log_file_path = catfile( $outsample_directory,
        $infile_prefix . $UNDERSCORE . q{chanjo_sexcheck.log} );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_name, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $core_number,
            directory_id          => $sample_id,
            FILEHANDLE            => $FILEHANDLE,
            log                   => $log,
            job_id_href           => $job_id_href,
            process_time          => $time,
            program_directory     => catfile(
                lc( $active_parameter_href->{outaligner_dir} ),
                q{coveragereport}
            ),
            program_name                    => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
        }
    );

    ## chanjo_sexcheck
    say {$FILEHANDLE} q{## Predicting sex from alignment};

    ## Get parameters

    my $chr_prefix;

    # If element is part of array
    if ( any { $_ eq q{chrX} } @{ $file_info_href->{contigs_size_ordered} } ) {
        $chr_prefix = q{chr};
    }
    chanjo_sex(
        {
            chr_prefix  => $chr_prefix,
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $infile_path,
            log_level   => $active_parameter_href->{chanjo_sexcheck_log_level},
            log_file_path => $log_file_path,
            outfile_path  => $outfile_path,
        }
    );
    say {$FILEHANDLE} $NEWLINE;
    close $FILEHANDLE;

    if ( $program_mode == 1 ) {

        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                infile       => $merged_infile_prefix,
                path         => catfile( $outsample_directory, $outfile_name ),
                program_name => q{chanjo_sexcheck},
                sample_id    => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );
        add_program_metafile_to_sample_info(
            {
                infile       => $merged_infile_prefix,
                metafile_tag => q{log},
                path         => catfile(
                    $outsample_directory,
                    $infile_prefix . $UNDERSCORE . q{chanjo_sexcheck.log}
                ),
                program_name     => q{chanjo_sexcheck},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );
        slurm_submit_job_sample_id_dependency_dead_end(
            {
                family_id               => $family_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                path                    => $job_id_chain,
                sample_id               => $sample_id,
                sbatch_file_name        => $file_name,
            }
        );
    }
    return;
}

1;
