package MIP::Script::Setup_script;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Path qw{ make_path };
use File::Spec::Functions qw{ catdir catfile devnull };
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

    use base qw{ Exporter };
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ setup_script write_return_to_conda_environment write_source_environment_command };
}

## Constant
Readonly my $DOT        => q{.};
Readonly my $EMPTY_STR  => q{};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub setup_script {

## Function : Creates program directories (info & program data & program script), program script filenames and writes sbatch header.
## Returns  : $file_path, $file_info_path . $file_name_version
## Arguments: $active_parameter_href           => The active parameters for this analysis hash {REF}
##          : $core_number                     => Number of cores to allocate {Optional}
##          : $directory_id                    => $samplID|$family_id
##          : $email_types_ref                 => Email type
##          : $error_trap                      => Error trap switch {Optional}
##          : $FILEHANDLE                      => FILEHANDLE to write to
##          : $job_id_href                     => The job_id hash {REF}
##          : $log                             => Log object
##          : $program_directory               => Builds from $directory_id/$outaligner_dir
##          : $program_name                    => Assigns filename to sbatch script
##          : $outdata_dir                     => MIP outdata directory {Optional}
##          : $outscript_dir                   => MIP outscript directory {Optional}
##          : $process_time                    => Allowed process time (Hours) {Optional}
##          : $set_errexit                     => Bash set -e {Optional}
##          : $set_nounset                     => Bash set -u {Optional}
##          : $set_pipefail                    => Pipe fail switch {Optional}
##          : $sleep                           => Sleep for X seconds {Optional}
##          : $slurm_quality_of_service        => SLURM quality of service priority {Optional}
##          : $source_environment_commands_ref => Source environment command {REF}
##          : $temp_directory                  => Temporary directory for program {Optional}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $directory_id;
    my $FILEHANDLE;
    my $job_id_href;
    my $log;
    my $program_directory;
    my $program_name;
    my $source_environment_commands_ref;

    ## Default(s)
    my $core_number;
    my $email_types_ref;
    my $error_trap;
    my $outdata_dir;
    my $outscript_dir;
    my $process_time;
    my $set_errexit;
    my $set_nounset;
    my $set_pipefail;
    my $sleep;
    my $slurm_quality_of_service;
    my $temp_directory;

    use MIP::Check::Parameter qw{ check_allowed_array_values };

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        core_number => {
            allow       => qr{ \A\d+\z }xsm,
            default     => 1,
            store       => \$core_number,
            strict_type => 1,
        },
        directory_id => {
            defined     => 1,
            required    => 1,
            store       => \$directory_id,
            strict_type => 1,
        },
        email_types_ref => {
            allow => [
                sub {
                    check_allowed_array_values(
                        {
                            allowed_values_ref => [qw{ NONE BEGIN END FAIL REQUEUE ALL }],
                            values_ref         => $arg_href->{email_types_ref},
                        }
                    );
                }
            ],
            default     => \@{ $arg_href->{active_parameter_href}{email_types} },
            store       => \$email_types_ref,
            strict_type => 1,
        },
        error_trap => {
            allow       => [ 0, 1 ],
            default     => 1,
            store       => \$error_trap,
            strict_type => 1,
        },
        FILEHANDLE  => { store => \$FILEHANDLE, },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        log         => { defined => 1, required => 1, store => \$log, },
        outdata_dir => {
            default     => $arg_href->{active_parameter_href}{outdata_dir},
            store       => \$outdata_dir,
            strict_type => 1,
        },
        outscript_dir => {
            default     => $arg_href->{active_parameter_href}{outscript_dir},
            store       => \$outscript_dir,
            strict_type => 1,
        },
        process_time => {
            allow       => qr{ \A\d+\z }xsm,
            default     => 1,
            store       => \$process_time,
            strict_type => 1,
        },
        program_directory => {
            defined     => 1,
            required    => 1,
            store       => \$program_directory,
            strict_type => 1,
        },
        program_name => {
            defined     => 1,
            required    => 1,
            store       => \$program_name,
            strict_type => 1,
        },
        set_errexit => {
            allow       => [ 0, 1 ],
            default     => $arg_href->{active_parameter_href}{bash_set_errexit},
            store       => \$set_errexit,
            strict_type => 1,
        },
        set_nounset => {
            allow       => [ 0, 1 ],
            default     => $arg_href->{active_parameter_href}{bash_set_nounset},
            store       => \$set_nounset,
            strict_type => 1,
        },
        set_pipefail => {
            allow       => [ 0, 1 ],
            default     => $arg_href->{active_parameter_href}{bash_set_pipefail},
            store       => \$set_pipefail,
            strict_type => 1,
        },
        sleep => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$sleep,
            strict_type => 1,
        },
        slurm_quality_of_service => {
            allow       => [qw{ low high normal }],
            default     => $arg_href->{active_parameter_href}{slurm_quality_of_service},
            store       => \$slurm_quality_of_service,
            strict_type => 1,
        },
        source_environment_commands_ref => {
            default     => [],
            store       => \$source_environment_commands_ref,
            strict_type => 1,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Path qw{ check_file_version_exist };
    use MIP::Gnu::Bash qw{ gnu_set };
    use MIP::Gnu::Coreutils qw{ gnu_echo gnu_mkdir gnu_sleep };
    use MIP::Language::Shell
      qw{ build_shebang create_housekeeping_function create_error_trap_function enable_trap quote_bash_variable };
    use MIP::Workloadmanager::Slurm qw{ slurm_build_sbatch_header };

    ## Constants
    Readonly my $MAX_SECONDS_TO_SLEEP => 240;

    ## Unpack parameters
    my $program_mode       = $active_parameter_href->{$program_name};
    my $submission_profile = $active_parameter_href->{submission_profile};

    my %submission_method = ( slurm => q{sbatch}, );
    my $script_type = $submission_method{$submission_profile};

    ### Script names and directory creation
    ## File
    my $file_name_prefix = $program_name . $UNDERSCORE . $directory_id . $DOT;
    my $file_name_version;
    my $file_name_suffix = $DOT . q{sh};

    # Path
    my $file_path;

    ## Directory paths
    my $program_data_directory_path =
      catdir( $outdata_dir, $directory_id, $program_directory );
    my $program_info_directory_path = catdir( $program_data_directory_path, q{info} );
    my $program_script_directory_path =
      catdir( $outscript_dir, $directory_id, $program_directory );

    ## Create directories
    make_path( $program_info_directory_path, $program_data_directory_path,
        $program_script_directory_path );

    ## File paths
    my $file_path_prefix = catfile( $program_script_directory_path, $file_name_prefix );
    my $file_info_path   = catfile( $program_info_directory_path,   $file_name_prefix );
    my $dry_run_file_path_prefix = catfile( $program_script_directory_path,
        q{dry_run} . $UNDERSCORE . $file_name_prefix );
    my $dry_run_file_info_path =
      catfile( $program_info_directory_path,
        q{dry_run} . $UNDERSCORE . $file_name_prefix );

    ## Dry run program - update file paths
    if ( $program_mode == 2 ) {

        $file_path_prefix = $dry_run_file_path_prefix;
        $file_info_path   = $dry_run_file_info_path;
        $log->info( q{Dry run:} . $NEWLINE );
    }

    ## Check if a file with with a filename consisting of
    ## $file_path_prefix.$file_name_version.$file_path_suffix exist
    ( $file_path, $file_name_version ) = check_file_version_exist(
        {
            file_path_prefix => $file_path_prefix,
            file_path_suffix => $file_name_suffix,
        }
    );

    ### Info and Log
    $log->info( q{Creating }
          . $script_type
          . q{ script for }
          . $program_name
          . q{ and writing script file(s) to: }
          . $file_path
          . $NEWLINE );
    $log->info(
            ucfirst $script_type
          . q{ script }
          . $program_name
          . q{ data files will be written to: }
          . $program_data_directory_path
          . $NEWLINE );

    ## Script file
    open $FILEHANDLE, q{>}, $file_path
      or
      $log->logdie( q{Cannot write to '} . $file_path . q{' :} . $OS_ERROR . $NEWLINE );

    # Build bash shebang line
    build_shebang(
        {
            bash_bin_path => catfile( dirname( dirname( devnull() ) ), qw{ bin bash } ),
            FILEHANDLE    => $FILEHANDLE,
            invoke_login_shell => 1,
        }
    );

    ## Set shell attributes
    gnu_set(
        {
            FILEHANDLE   => $FILEHANDLE,
            set_errexit  => $set_errexit,
            set_nounset  => $set_nounset,
            set_pipefail => $set_pipefail,
        }
    );

    ### Sbatch header
    ## Get parameters
    my $job_name        = $program_name . $UNDERSCORE . $directory_id;
    my $stderrfile_path = $file_info_path . $file_name_version . $DOT . q{stderr.txt};
    my $stdoutfile_path = $file_info_path . $file_name_version . $DOT . q{stdout.txt};

    ## SLURM specific headers and parameters
    my @sbatch_headers;
    my @sacct_format_fields;
    if ( $submission_profile eq q{slurm} ) {

        @sacct_format_fields = @{ $active_parameter_href->{sacct_format_fields} };
        @sbatch_headers      = slurm_build_sbatch_header(
            {
                core_number              => $core_number,
                email                    => $active_parameter_href->{email},
                email_types_ref          => $email_types_ref,
                FILEHANDLE               => $FILEHANDLE,
                job_name                 => $job_name,
                process_time             => $process_time . q{:00:00},
                project_id               => $active_parameter_href->{project_id},
                slurm_quality_of_service => $slurm_quality_of_service,
                stderrfile_path          => $stderrfile_path,
                stdoutfile_path          => $stdoutfile_path,
            }
        );
    }
    say {$FILEHANDLE} q{readonly PROGNAME=$(basename "$0")}, $NEWLINE;

    gnu_echo(
        {
            FILEHANDLE  => $FILEHANDLE,
            strings_ref => [q{Running on: $(hostname)}],
        }
    );
    say {$FILEHANDLE} $NEWLINE;

# Let the process sleep for a random couple of seconds (0-60) to avoid race conditions in mainly conda sourcing activate
    if ($sleep) {

        gnu_sleep(
            {
                FILEHANDLE       => $FILEHANDLE,
                seconds_to_sleep => int rand $MAX_SECONDS_TO_SLEEP,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }
    if ( @{$source_environment_commands_ref}
        && $source_environment_commands_ref->[0] )
    {

        write_source_environment_command(
            {
                FILEHANDLE                      => $FILEHANDLE,
                source_environment_commands_ref => $source_environment_commands_ref,
            }
        );
    }

    # Not all programs need a temporary directory
    if ( defined $temp_directory ) {

        say {$FILEHANDLE} q{## Create temporary directory};

        ## Double quote incoming variables in string
        my $temp_directory_quoted =
          quote_bash_variable( { string_with_variable_to_quote => $temp_directory, } );

        # Assign batch variable
        say {$FILEHANDLE} q{readonly TEMP_DIRECTORY=} . $temp_directory_quoted;

        # Update perl scalar to bash variable
        my $temp_directory_bash = q{"$TEMP_DIRECTORY"};

        gnu_mkdir(
            {
                FILEHANDLE       => $FILEHANDLE,
                indirectory_path => $temp_directory_bash,
                parents          => 1,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        create_housekeeping_function(
            {
                FILEHANDLE              => $FILEHANDLE,
                job_ids_ref             => \@{ $job_id_href->{PAN}{PAN} },
                log_file_path           => $active_parameter_href->{log_file},
                remove_dir              => $temp_directory_bash,
                sacct_format_fields_ref => \@sacct_format_fields,
                trap_function_call      => q{$(finish } . $temp_directory_bash . q{)},
                trap_function_name      => q{finish},
                trap_signals_ref        => [qw{ EXIT TERM INT }],
            }
        );
    }

    if ($error_trap) {

        ## Create debug trap
        enable_trap(
            {
                FILEHANDLE         => $FILEHANDLE,
                trap_function_call => q{previous_command="$BASH_COMMAND"},
                trap_signals_ref   => [qw{ DEBUG }],
            }
        );

        ## Create error handling function and trap
        create_error_trap_function(
            {
                FILEHANDLE              => $FILEHANDLE,
                job_ids_ref             => \@{ $job_id_href->{PAN}{PAN} },
                log_file_path           => $active_parameter_href->{log_file},
                sacct_format_fields_ref => \@sacct_format_fields,
                trap_signals_ref        => [qw{ ERR }],
                trap_function_call      => q{$(error "$previous_command" "$?")},
                trap_function_name      => q{error},
            }
        );
    }

    # Return filen path, file path for stdout/stderr for QC check later
    return ( $file_path, $file_info_path . $file_name_version );
}

sub write_return_to_conda_environment {

## Function : Return to main or default environment using conda
## Returns  :
## Arguments: $FILEHANDLE                           => Filehandle to write to
##          : $source_main_environment_commands_ref => Source main environment command {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $source_main_environment_commands_ref;

    my $tmpl = {
        FILEHANDLE                           => { required => 1, store => \$FILEHANDLE, },
        source_main_environment_commands_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$source_main_environment_commands_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Package_manager::Conda qw{ conda_source_deactivate };

    ## Return to main environment
    if ( @{$source_main_environment_commands_ref}
        && $source_main_environment_commands_ref->[0] )
    {

        write_source_environment_command(
            {
                FILEHANDLE => $FILEHANDLE,
                source_environment_commands_ref =>
                  \@{$source_main_environment_commands_ref},
            }
        );
    }
    else {
        ## Return to login shell environment

        say {$FILEHANDLE} q{## Deactivate environment};
        conda_source_deactivate(
            {
                FILEHANDLE => $FILEHANDLE,
            }
        );
        print {$FILEHANDLE} $NEWLINE;
    }
    return;
}

sub write_source_environment_command {

## Function : Write source environment commmands to filehandle
## Returns  :
## Arguments: $FILEHANDLE                      => Filehandle to write to
##          : $source_environment_commands_ref => Source environment command {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $source_environment_commands_ref;

    my $tmpl = {
        FILEHANDLE                      => { required => 1, store => \$FILEHANDLE, },
        source_environment_commands_ref => {
            default     => [],
            store       => \$source_environment_commands_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Unix::Write_to_file qw{ unix_write_to_file };

    if ( @{$source_environment_commands_ref} ) {

        say {$FILEHANDLE} q{## Activate environment};

        unix_write_to_file(
            {
                commands_ref => $source_environment_commands_ref,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }
    return;
}

1;
