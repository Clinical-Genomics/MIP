package MIP::Recipes::Qc::Multiqc;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use autodie qw{ :all };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use File::Spec::Functions qw{ catdir catfile };

## CPANM
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_multiqc };

}

## Constants
Readonly my $NEWLINE => qq{\n};

sub analysis_multiqc {

## Function : Aggregate bioinforamtics reports per case
## Returns  :
## Arguments: $parameter_href          => Parameter hash {REF}
##          : $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $program_name            => Program name
##          : $family_id               => Family id {REF}

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $program_name;

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
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_chain_job_ids_dependency_add_to_path };
    use MIP::Program::Qc::Multiqc qw{ multiqc };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Alias
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
    my $core_number =
      $active_parameter_href->{module_core_number}{$mip_program_name};
    my $time = $active_parameter_href->{module_time}{$mip_program_name};

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
            program_directory     => lc($program_name),
            core_number           => $core_number,
            process_time          => $time,
        }
    );

    ## Assign directories
    my $program_outdirectory_name =
      $parameter_href->{$mip_program_name}{outdir_name};

  SAMPLE:
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        ## Assign directories
        my $insample_directory =
          catdir( $active_parameter_href->{outdata_dir}, $sample_id );
        my $outsample_directory = catdir( $active_parameter_href->{outdata_dir},
            $sample_id, $program_outdirectory_name );

        multiqc(
            {
                indir_path  => $insample_directory,
                outdir_path => $outsample_directory,
                force       => 1,
                FILEHANDLE  => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    close $FILEHANDLE;

    if ( $mip_program_mode == 1 ) {

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
