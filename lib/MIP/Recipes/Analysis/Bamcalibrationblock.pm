package MIP::Recipes::Analysis::Bamcalibrationblock;

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
    our @EXPORT_OK = qw{ analysis_bamcalibrationblock };

}

##Constants
Readonly my $TAB        => qq{\t};
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub analysis_bamcalibrationblock {

## Function : Run consecutive module
## Returns  :
## Arguments: $parameter_href          => Parameter hash {REF}
##          : $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $lane_href               => The lane info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => The outaligner_dir used
##          : $program_name            => Program name
##          : $temp_directory          => Temporary directory
##          : $outaligner_dir          => Outaligner dir used in the analysis
##          : $log                     => Log object to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $lane_href;
    my $job_id_href;
    my $program_name;
    my $log;

    ## Default(s)
    my $temp_directory;
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
        lane_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$lane_href
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
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory,
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir,
        },
        log => {
            required => 1,
            defined  => 1,
            store    => \$log
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Recipes::Analysis::Markduplicates qw{ analysis_markduplicates_rio };
    use MIP::Recipes::Analysis::Picardtools_mergesamfiles
      qw{ analysis_picardtools_mergesamfiles_rio };
    use MIP::Recipes::Analysis::Gatk_realigner qw{ analysis_gatk_realigner_rio };
    use MIP::Recipes::Analysis::Gatk_baserecalibration
      qw{ analysis_gatk_baserecalibration_rio };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Constants
    Readonly my $PROCESS_TIME => 80;

    my $core_number = $active_parameter_href->{max_cores_per_node};

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ### Always run Picardtools mergesamfiles even for single samples to rename them correctly for standardised downstream processing.
    ## Will also split alignment per contig and copy to temporary directory for '-rio 1' block to enable selective removal of block submodules.
    if ( $active_parameter_href->{ppicardtools_mergesamfiles} > 0 ) {

        $log->info( $TAB . q{[Picardtools mergesamfiles]} );
    }
    ## Markduplicates
    if ( $active_parameter_href->{pmarkduplicates} > 0 ) {

        $log->info( $TAB . q{[Markduplicates]} );
    }

    ## Run GATK realignertargetcreator/indelrealigner
    if ( $active_parameter_href->{pgatk_realigner} > 0 ) {

        $log->info( $TAB . q{[GATK realignertargetcreator/indelrealigner]} );
    }
    ## Run GATK baserecalibrator/printreads
    if ( $active_parameter_href->{pgatk_baserecalibration} > 0 ) {

        $log->info( $TAB . q{[GATK baserecalibrator/printreads]} );
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
                job_id_href           => $job_id_href,
                FILEHANDLE            => $FILEHANDLE,
                directory_id          => $sample_id,
                program_name          => $program_name,
                program_directory     => $outaligner_dir,
                core_number           => $core_number,
                process_time          => $PROCESS_TIME,
                temp_directory        => $temp_directory,
            }
        );

        ## Always run Picardtools mergesamfiles even for single samples to rename them correctly for standardised downstream processing.
        ## Will also split alignment per contig and copy to temporary directory for -rio 1 block to enable selective removal of block submodules.
        if ( $active_parameter_href->{ppicardtools_mergesamfiles} > 0 ) {

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
        if ( $active_parameter_href->{pmarkduplicates} > 0 ) {

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
        if ( $active_parameter_href->{pgatk_realigner} > 0 ) {

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
        if ( $active_parameter_href->{pgatk_baserecalibration} > 0 ) {

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
