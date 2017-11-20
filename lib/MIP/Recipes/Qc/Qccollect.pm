package MIP::Recipes::QC::Qccollect;

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
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_qccollect };
}

## Constants
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub analysis_qccollect {

    ## Function : Collect qc metrics for this analysis run.
    ## Returns  :
    ## Arguments: $parameter_href          => Parameter hash {REF}
    ##          : $active_parameter_href   => Active parameters for this analysis hash {REF}
    ##          : $sample_info_href        => Info on samples and family hash {REF}
    ##          : $file_info_href          => File_info hash {REF}
    ##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
    ##          : $job_id_href             => Job id hash {REF}
    ##          : $infile_path             => Infile path
    ##          : $outfamily_directory     => out Family Directory
    ##          : $program_name            => Program name
    ##          : $family_id               => Family id
    ##          : $outaligner_dir          => Outaligner_dir used in the analysis

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $infile_path;
    my $outfamily_directory;
    my $program_name;

    ## Default(s)
    my $family_id;
    my $outaligner_dir;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        infile_path => {
            default =>
              $arg_href->{active_parameter_href}{qccollect_sampleinfo_file},
            strict_type => 1,
            store       => \$infile_path,
        },
        outfamily_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfamily_directory,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id,
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ### SET IMPORT IN ALPHABETIC ORDER
    use MIP::Get::File qw{ get_file_suffix get_merged_infile_prefix };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_chain_job_ids_dependency_add_to_path };
    use MIP::Program::Qc::Qccollect qw{ qccollect };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Unpack parameters
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
    my $core_number =
      $active_parameter_href->{module_core_number}{$mip_program_name};
    my $time = $active_parameter_href->{module_time}{$mip_program_name};
    my $reduce_io_ref = $active_parameter_href->{reduce_io};

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($file_path) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $family_id,
            program_name          => $program_name,
            program_directory     => $program_name,
            core_number           => $core_number,
            process_time          => $time,
        }
    );

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

    if ( $mip_program_mode == 1 ) {

        ## Collect QC metadata info for later use
        my $qc_metric_outfile = $family_id . $UNDERSCORE . q{qc_metrics.yaml};
        my $path = catfile( $outfamily_directory, $qc_metric_outfile );

        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => q{qccollect},
                outdirectory     => $outfamily_directory,
                outfile          => $qc_metric_outfile,
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
