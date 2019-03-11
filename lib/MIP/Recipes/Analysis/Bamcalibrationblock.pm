package MIP::Recipes::Analysis::Bamcalibrationblock;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
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
    our @EXPORT_OK = qw{ analysis_bamcalibrationblock };

}

## Constants
Readonly my $TAB        => qq{\t};
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};
Readonly my $SPACE => q{ };

sub analysis_bamcalibrationblock {

## Function : Run consecutive module
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $lane_href               => The lane info hash {REF}
##          : $log                     => Log object to write to
##          : $outaligner_dir          => Outaligner dir used
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $lane_href;
    my $log;
    my $parameter_href;
    my $program_name;
    my $sample_info_href;

    ## Default(s)
    my $outaligner_dir;
    my $temp_directory;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
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
        lane_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$lane_href
        },
        log => {
            required => 1,
            defined  => 1,
            store    => \$log
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir,
        },
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Recipes::Analysis::Gatk_realigner
      qw{ analysis_gatk_realigner_rio };
    use MIP::Recipes::Analysis::Gatk_baserecalibration
      qw{ analysis_gatk_baserecalibration_rio };
    use MIP::Recipes::Analysis::Markduplicates
      qw{ analysis_markduplicates_rio };
    use MIP::Recipes::Analysis::Picardtools_mergesamfiles
      qw{ analysis_picardtools_mergesamfiles_rio };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Constants
    Readonly my $PROCESS_TIME => 80;

    my $core_number = $active_parameter_href->{max_cores_per_node};
    my $source_environment_cmd = join $SPACE, @{ $active_parameter_href->{source_main_environment_commands} };

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    # Set order of supplying user info
    my @rio_program_order =
      qw{ ppicardtools_mergesamfiles pmarkduplicates pgatk_realigner pgatk_baserecalibration };

    # Store what to supply to user
    my %rio_program = (
        ppicardtools_mergesamfiles => q{[Picardtools mergesamfiles]},
        pmarkduplicates            => q{[Markduplicates]},
        pgatk_realigner => q{[GATK realignertargetcreator/indelrealigner]},
        pgatk_baserecalibration => q{[GATK baserecalibrator/printreads]},
    );

  RIO_PROGRAM:
    foreach my $mip_program_name (@rio_program_order) {

        if ( $active_parameter_href->{$mip_program_name} ) {

            my $program_header = $rio_program{$mip_program_name};

            $log->info( $TAB . $program_header );
        }
    }

  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        my $xargs_file_counter = 0;

        ## Assign directories
        my $insample_directory = catdir( $active_parameter_href->{outdata_dir},
            $sample_id, $outaligner_dir );
        my $outsample_directory = catdir( $active_parameter_href->{outdata_dir},
            $sample_id, $outaligner_dir );

        ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
        my ( $file_path, $program_info_path ) = setup_script(
            {
                active_parameter_href => $active_parameter_href,
                core_number           => $core_number,
                directory_id          => $sample_id,
                FILEHANDLE            => $FILEHANDLE,
                job_id_href           => $job_id_href,
                process_time          => $PROCESS_TIME,
                program_directory     => $outaligner_dir,
                program_name          => $program_name,
                source_environment_commands_ref => [$source_environment_cmd],
                temp_directory        => $temp_directory,
            }
        );

        ## Always run Picardtools mergesamfiles even for single samples to rename them correctly for standardised downstream processing.
        ## Will also split alignment per contig and copy to temporary directory for -rio 1 block to enable selective removal of block submodules.
        if ( $active_parameter_href->{ppicardtools_mergesamfiles} ) {

            ($xargs_file_counter) = analysis_picardtools_mergesamfiles_rio(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    lane_href               => $lane_href,
                    job_id_href             => $job_id_href,
                    insample_directory      => $insample_directory,
                    sample_id               => $sample_id,
                    program_name            => q{picardtools_mergesamfiles},
                    file_path               => $file_path,
                    program_info_path       => $program_info_path,
                    FILEHANDLE              => $FILEHANDLE,
                }
            );
        }

        # Markduplicates
        if ( $active_parameter_href->{pmarkduplicates} ) {

            ($xargs_file_counter) = analysis_markduplicates_rio(
                {
                    parameter_href        => $parameter_href,
                    active_parameter_href => $active_parameter_href,
                    sample_info_href      => $sample_info_href,
                    file_info_href        => $file_info_href,
                    sample_id             => $sample_id,
                    program_name          => q{markduplicates},
                    file_path             => $file_path,
                    program_info_path     => $program_info_path,
                    FILEHANDLE            => $FILEHANDLE,
                    xargs_file_counter    => $xargs_file_counter,
                }
            );
        }

        ## Run GATK realignertargetcreator/indelrealigner
        if ( $active_parameter_href->{pgatk_realigner} ) {

            ($xargs_file_counter) = analysis_gatk_realigner_rio(
                {
                    parameter_href        => $parameter_href,
                    active_parameter_href => $active_parameter_href,
                    file_info_href        => $file_info_href,
                    sample_id             => $sample_id,
                    program_name          => q{gatk_realigner},
                    file_path             => $file_path,
                    program_info_path     => $program_info_path,
                    FILEHANDLE            => $FILEHANDLE,
                    xargs_file_counter    => $xargs_file_counter,
                }
            );
        }

        ## Run GATK baserecalibrator/printreads
        if ( $active_parameter_href->{pgatk_baserecalibration} ) {

            ($xargs_file_counter) = analysis_gatk_baserecalibration_rio(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    sample_info_href        => $sample_info_href,
                    file_info_href          => $file_info_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    sample_id               => $sample_id,
                    outsample_directory     => $outsample_directory,
                    program_name            => q{gatk_baserecalibration},
                    file_path               => $file_path,
                    program_info_path       => $program_info_path,
                    FILEHANDLE              => $FILEHANDLE,
                    xargs_file_counter      => $xargs_file_counter,
                }
            );
        }
    }
    return;
}

1;
