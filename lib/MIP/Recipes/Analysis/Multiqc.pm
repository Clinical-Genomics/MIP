package MIP::Recipes::Analysis::Multiqc;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use Params::Check qw{ check allow last_error };
use open qw{ :encoding(UTF-8) :std };
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
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_multiqc };

}

## Constants
Readonly my $NEWLINE => qq{\n};

sub analysis_multiqc {

## Function : Aggregate bioinforamtics reports per case
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $sample_info_href        => Info on samples and family hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $program_name;
    my $sample_info_href;

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
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_chain_job_ids_dependency_add_to_path };
    use MIP::Program::Qc::Multiqc qw{ multiqc };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::QC::Record qw{ add_program_metafile_to_sample_info };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set program mode
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

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($file_path) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $family_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            log                             => $log,
            process_time                    => $time,
            program_directory               => $program_name,
            program_name                    => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
        }
    );

    ## Assign directories
    my $program_outdirectory_name =
      $parameter_href->{$program_name}{outdir_name};

    ## Always analyse case
    my @report_ids = ($family_id);

    ## Generate report per sample id
    if ( $active_parameter_href->{multiqc_per_sample} ) {

        ## Add samples to analysis
        push @report_ids, @{ $active_parameter_href->{sample_ids} };
    }

  REPORT_ID:
    foreach my $report_id (@report_ids) {

        ## Assign directories
        my $indirectory  = catdir( $active_parameter_href->{outdata_dir} );
        my $outdirectory = catdir( $active_parameter_href->{outdata_dir},
            $report_id, $program_outdirectory_name );

        ## Analyse sample id only for this report
        if ( $report_id ne $family_id ) {

            $indirectory =
              catdir( $active_parameter_href->{outdata_dir}, $report_id );
        }

        multiqc(
            {
                FILEHANDLE  => $FILEHANDLE,
                force       => 1,
                indir_path  => $indirectory,
                outdir_path => $outdirectory,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        if ( $program_mode == 1 ) {

            ## Collect QC metadata info for later use
            add_program_metafile_to_sample_info(
                {
                    metafile_tag => $report_id,
                    path => catfile( $outdirectory, q{multiqc_report.html} ),
                    program_name     => q{multiqc},
                    sample_info_href => $sample_info_href,
                }
            );
        }

    }

    close $FILEHANDLE;

    if ( $program_mode == 1 ) {

        slurm_submit_chain_job_ids_dependency_add_to_path(
            {
                job_id_href      => $job_id_href,
                log              => $log,
                path             => $job_id_chain,
                sbatch_file_name => $file_path,
            }
        );
    }
    return;
}

1;
