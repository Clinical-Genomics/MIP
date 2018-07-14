package MIP::Recipes::Analysis::Qccollect;

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
use autodie qw{ :all };
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_qccollect };
}

## Constants
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub analysis_qccollect {

## Function : Collect qc metrics for this analysis run.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $infile_path             => Infile path
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $sample_info_href        => Info on samples and family hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $infile_path;
    my $job_id_href;
    my $parameter_href;
    my $program_name;
    my $sample_info_href;

    ## Default(s)
    my $family_id;
    my $outaligner_dir;

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
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_lane_prefix_href,
            strict_type => 1,
        },
        infile_path => {
            default =>
              $arg_href->{active_parameter_href}{qccollect_sampleinfo_file},
            store       => \$infile_path,
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
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            store       => \$outaligner_dir,
            strict_type => 1,
        },
        program_name => {
            defined     => 1,
            required    => 1,
            store       => \$program_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_file_suffix get_merged_infile_prefix };
    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_chain_job_ids_dependency_add_to_path };
    use MIP::Program::Qc::Qccollect qw{ qccollect };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set program mode
    my $program_mode = $active_parameter_href->{$program_name};

    ## Unpack parameters
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

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($file_path) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $family_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            log                             => $log,
            program_directory               => $program_name,
            program_name                    => $program_name,
            process_time                    => $time,
            source_environment_commands_ref => \@source_environment_cmds,
        }
    );

    ## Assign directories
    my $outfamily_directory =
      catdir( $active_parameter_href->{outdata_dir}, $family_id );

    ## Paths
    my $outfile_path = catfile( $outfamily_directory,
        $family_id . $UNDERSCORE . q{qc_metrics.yaml} );
    my $log_file_path = catfile( $outfamily_directory,
        $family_id . $UNDERSCORE . q{qccollect.log} );

    qccollect(
        {
            infile_path      => $infile_path,
            outfile_path     => $outfile_path,
            log_file_path    => $log_file_path,
            regexp_file_path => $active_parameter_href->{qccollect_regexp_file},
            skip_evaluation =>
              $active_parameter_href->{qccollect_skip_evaluation},
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

    if ( $program_mode == 1 ) {

        ## Collect QC metadata info for later use
        my $qc_metric_outfile = $family_id . $UNDERSCORE . q{qc_metrics.yaml};
        my $path = catfile( $outfamily_directory, $qc_metric_outfile );

        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => q{qccollect},
                path             => $path,
            }
        );

        slurm_submit_chain_job_ids_dependency_add_to_path(
            {
                job_id_href      => $job_id_href,
                path             => $job_id_chain,
                log              => $log,
                sbatch_file_name => $file_path,
            }
        );
    }
    return;
}

1;
