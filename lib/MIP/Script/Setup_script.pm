package MIP::Script::Setup_script;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Basename qw{ dirname fileparse };
use File::Path qw{ make_path };
use File::Spec::Functions qw{ catdir catfile devnull };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use Time::Piece;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants
  qw{ $DOT $EMPTY_STR $LOG_NAME $NEWLINE $SINGLE_QUOTE $SPACE $UNDERSCORE };

BEGIN {

    use base qw{ Exporter };
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.14;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      build_script_directories_and_paths
      check_script_file_path_exist
      setup_install_script
      setup_script
      write_return_to_environment
      write_return_to_conda_environment
      write_source_environment_command };
}

sub build_script_directories_and_paths {

## Function : Builds and makes recipe directories (info & data & script) and recipe script paths
## Returns  : $file_info_path, $file_path_prefix, $recipe_data_directory_path
## Arguments: $directory_id               => $sample id | $case_id
##          : $outdata_dir                => MIP outdata directory
##          : $outscript_dir              => MIP outscript directory
##          : $recipe_data_directory_path => Set recipe data directory path
##          : $recipe_directory           => Builds from $directory_id
##          : $recipe_mode                => Recipe mode
##          : $recipe_name                => Assigns filename to sbatch script

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $directory_id;
    my $outdata_dir;
    my $outscript_dir;
    my $recipe_data_directory_path;
    my $recipe_directory;
    my $recipe_mode;
    my $recipe_name;

    my $tmpl = {
        directory_id => {
            defined     => 1,
            required    => 1,
            store       => \$directory_id,
            strict_type => 1,
        },
        outdata_dir => {
            defined     => 1,
            required    => 1,
            store       => \$outdata_dir,
            strict_type => 1,
        },
        outscript_dir => {
            defined     => 1,
            required    => 1,
            store       => \$outscript_dir,
            strict_type => 1,
        },
        recipe_data_directory_path => {
            required    => 1,
            store       => \$recipe_data_directory_path,
            strict_type => 1,
        },
        recipe_directory => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_directory,
            strict_type => 1,
        },
        recipe_mode => => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_mode,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my $file_name_prefix = $recipe_name . $UNDERSCORE . $directory_id . $DOT;

    ## Directory paths
    if ( not defined $recipe_data_directory_path ) {

        $recipe_data_directory_path =
          catdir( $outdata_dir, $directory_id, $recipe_directory );
    }
    my $recipe_info_directory_path = catdir( $recipe_data_directory_path, q{info} );
    my $recipe_script_directory_path =
      catdir( $outscript_dir, $directory_id, $recipe_directory );

    ## Create directories
    make_path( $recipe_info_directory_path, $recipe_data_directory_path,
        $recipe_script_directory_path );

    ## File paths
    my $file_path_prefix = catfile( $recipe_script_directory_path, $file_name_prefix );
    my $file_info_path   = catfile( $recipe_info_directory_path,   $file_name_prefix );
    my $dry_run_file_path_prefix = catfile( $recipe_script_directory_path,
        q{dry_run} . $UNDERSCORE . $file_name_prefix );
    my $dry_run_file_info_path =
      catfile( $recipe_info_directory_path,
        q{dry_run} . $UNDERSCORE . $file_name_prefix );

    ## Dry run recipe - update file paths
    if ( $recipe_mode == 2 ) {

        $file_path_prefix = $dry_run_file_path_prefix;
        $file_info_path   = $dry_run_file_info_path;
        $log->info( q{Dry run:} . $NEWLINE );
    }
    return ( $file_info_path, $file_path_prefix, $recipe_data_directory_path );
}

sub check_script_file_path_exist {

## Function : Check if a file with with a filename consisting of $file_path_prefix.$file_counter.$file_path_suffix exist.
##          : If so bumps the version number and return new file path and version number
## Returns  : $file_path, $file_name_version
## Arguments: $file_path_prefix => File path prefix
##          : $file_path_suffix => File path suffix

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path_prefix;
    my $file_path_suffix;

    my $tmpl = {
        file_path_prefix => {
            defined     => 1,
            required    => 1,
            store       => \$file_path_prefix,
            strict_type => 1,
        },
        file_path_suffix => {
            defined     => 1,
            required    => 1,
            store       => \$file_path_suffix,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Nr of scripts with identical file names i.e. version number
    my $file_name_version = 0;

    my $file_path = $file_path_prefix . $file_name_version . $file_path_suffix;

  FILE_PATHS:
    while ( -e $file_path ) {

        $file_name_version++;

        ## New file_path to test for existence
        $file_path = $file_path_prefix . $file_name_version . $file_path_suffix;
    }
    return ( $file_path, $file_name_version );
}

sub setup_install_script {

## Function : Build bash file with header
## Returns  :
## Arguments: $active_parameter_href => Master hash, used when running in sbatch mode {REF}
##          : $file_name             => File name
##          : $filehandle            => Filehandle to write to
##          : $remove_dir            => Directory to remove when caught by trap function
##          : $log                   => Log object to write to
##          : $invoke_login_shell    => Invoked as a login shell. Reinitilize bashrc and bash_profile
##          : $sbatch_mode           => Create headers for sbatch submission;
##          : $set_errexit           => Halt script if command has non-zero exit code (-e)
##          : $set_nounset           => Halt script if variable is uninitialised (-u)
##          : $set_pipefail          => Detect errors within pipes (-o pipefail)

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_name;
    my $filehandle;
    my $log;
    my $remove_dir;

    ## Default(s)
    my $invoke_login_shell;
    my $sbatch_mode;
    my $set_errexit;
    my $set_nounset;
    my $set_pipefail;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        filehandle => {
            required => 1,
            store    => \$filehandle,
        },
        file_name => {
            defined     => 1,
            required    => 1,
            store       => \$file_name,
            strict_type => 1,
        },
        invoke_login_shell => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$invoke_login_shell,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        remove_dir => {
            allow       => qr/ ^\S+$ /xsm,
            store       => \$remove_dir,
            strict_type => 1,
        },
        sbatch_mode => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$sbatch_mode,
            strict_type => 1,
        },
        set_errexit => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$set_errexit,
            strict_type => 1,
        },
        set_nounset => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$set_nounset,
            strict_type => 1,
        },
        set_pipefail => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$set_pipefail,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Gnu::Bash qw{ gnu_set };
    use MIP::Workloadmanager::Slurm qw{ slurm_build_sbatch_header };

    ## Set $bash_bin_path default
    my $bash_bin_path =
      catfile( dirname( dirname( devnull() ) ), qw{ usr bin env bash } );

    if ($sbatch_mode) {
        $bash_bin_path =
          catfile( dirname( dirname( devnull() ) ), qw{ bin bash } );
    }

    ## Build bash shebang line
    build_shebang(
        {
            bash_bin_path      => $bash_bin_path,
            filehandle         => $filehandle,
            invoke_login_shell => $invoke_login_shell,
        }
    );

    ## Set shell attributes
    gnu_set(
        {
            filehandle   => $filehandle,
            set_errexit  => $set_errexit,
            set_nounset  => $set_nounset,
            set_pipefail => $set_pipefail,
        }
    );

    if ($sbatch_mode) {

        ## Get local time
        my $date_time       = localtime;
        my $date_time_stamp = $date_time->datetime;

        ## Get bash_file_name minus suffix and add time stamp.
        my $job_name =
          fileparse( $file_name, qr/\.[^.]*/xms ) . $UNDERSCORE . $date_time_stamp;

        ## Set STDERR/STDOUT paths
        my $stderrfile_path = catfile( cwd(), $job_name . $DOT . q{stderr.txt} );
        my $stdoutfile_path = catfile( cwd(), $job_name . $DOT . q{stdout.txt} );

        slurm_build_sbatch_header(
            {
                core_number     => $active_parameter_href->{core_number},
                email           => $active_parameter_href->{email},
                email_types_ref => $active_parameter_href->{email_types},
                filehandle      => $filehandle,
                job_name        => $job_name,
                process_time    => $active_parameter_href->{process_time},
                project_id      => $active_parameter_href->{project_id},
                slurm_quality_of_service =>
                  $active_parameter_href->{slurm_quality_of_service},
                stderrfile_path => $stderrfile_path,
                stdoutfile_path => $stdoutfile_path,
            }
        );
    }

    ## Create housekeeping function which removes entire directory when finished
    create_housekeeping_function(
        {
            filehandle         => $filehandle,
            remove_dir         => $remove_dir,
            trap_function_name => q{finish},
        }
    );

    ## Create debug trap
    enable_trap(
        {
            filehandle         => $filehandle,
            trap_function_call => q{previous_command="$BASH_COMMAND"},
            trap_signals_ref   => [qw{ DEBUG }],
        }
    );

    ## Create error handling function and trap
    create_error_trap_function( { filehandle => $filehandle, } );

    $log->info( q{Created bash file: '} . catfile($file_name), $SINGLE_QUOTE );

    return;
}

sub setup_script {

## Function : Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header.
## Returns  : $file_path, $file_info_path . $file_name_version
## Arguments: $active_parameter_href           => The active parameters for this analysis hash {REF}
##          : $core_number                     => Number of cores to allocate {Optional}
##          : $directory_id                    => $sample id | $case_id
##          : $email_types_ref                 => Email type
##          : $error_trap                      => Error trap switch {Optional}
##          : $filehandle                      => filehandle to write to
##          : $job_id_href                     => The job_id hash {REF}
##          : $log                             => Log object
##          : $memory_allocation               => Memory allocation
##          : $outdata_dir                     => MIP outdata directory {Optional}
##          : $outscript_dir                   => MIP outscript directory {Optional}
##          : $process_time                    => Allowed process time (Hours) {Optional}
##          : $recipe_data_directory_path      => Set recipe data directory path
##          : $recipe_directory                => Builds from $directory_id
##          : $recipe_name                     => Assigns filename to sbatch script
##          : $set_errexit                     => Bash set -e {Optional}
##          : $set_nounset                     => Bash set -u {Optional}
##          : $set_pipefail                    => Pipe fail switch {Optional}
##          : $sleep                           => Sleep for X seconds {Optional}
##          : $slurm_quality_of_service        => SLURM quality of service priority {Optional}
##          : $source_environment_commands_ref => Source environment command {REF}
##          : $temp_directory                  => Temporary directory for recipe {Optional}
##          : $ulimit_n                        => Set ulimit -n for recipe {Optional}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $directory_id;
    my $filehandle;
    my $job_id_href;
    my $memory_allocation;
    my $log;
    my $recipe_data_directory_path;
    my $recipe_directory;
    my $recipe_name;
    my $source_environment_commands_ref;
    my $ulimit_n;

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
        filehandle  => { store => \$filehandle, },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        log               => { defined => 1, required => 1, store => \$log, },
        memory_allocation => {
            allow       => [ undef, qr{ \A\d+\z }sxm ],
            store       => \$memory_allocation,
            strict_type => 1,
        },
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
        recipe_data_directory_path => {
            store       => \$recipe_data_directory_path,
            strict_type => 1,
        },
        recipe_directory => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_directory,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
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
        ulimit_n => {
            allow       => [ undef, qr/ \A \d+ \z /xms ],
            store       => \$ulimit_n,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::List qw{ check_allowed_array_values };
    use MIP::Program::Gnu::Bash qw{ gnu_set gnu_ulimit };
    use MIP::Program::Gnu::Coreutils qw{ gnu_echo gnu_mkdir gnu_sleep };
    use MIP::Language::Shell
      qw{ build_shebang create_housekeeping_function create_error_trap_function enable_trap quote_bash_variable };
    use MIP::Workloadmanager::Slurm qw{ slurm_build_sbatch_header };

    ## Constants
    Readonly my $MAX_SECONDS_TO_SLEEP => 240;
    Readonly my %SUBMISSION_METHOD    => ( slurm => q{sbatch}, );

    ## Unpack parameters
    my $submission_profile = $active_parameter_href->{submission_profile};

    my $script_type = $SUBMISSION_METHOD{$submission_profile};

    ### Script names and directory creation
    ## File
    my $file_name_suffix = $DOT . q{sh};

    ( my ( $file_info_path, $file_path_prefix ), $recipe_data_directory_path ) =
      build_script_directories_and_paths(
        {
            directory_id               => $directory_id,
            outdata_dir                => $outdata_dir,
            outscript_dir              => $outscript_dir,
            recipe_data_directory_path => $recipe_data_directory_path,
            recipe_directory           => $recipe_directory,
            recipe_mode                => $active_parameter_href->{$recipe_name},
            recipe_name                => $recipe_name,
        }
      );

    ## Check if a file with with a filename consisting of
    ## $file_path_prefix.$file_name_version.$file_path_suffix exist
    my ( $file_path, $file_name_version ) = check_script_file_path_exist(
        {
            file_path_prefix => $file_path_prefix,
            file_path_suffix => $file_name_suffix,
        }
    );

    ### Info and Log
    $log->info( q{Creating }
          . $script_type
          . q{ script for }
          . $recipe_name
          . q{ and writing script file(s) to: }
          . $file_path
          . $NEWLINE );
    $log->info(
            ucfirst $script_type
          . q{ script }
          . $recipe_name
          . q{ data files will be written to: }
          . $recipe_data_directory_path
          . $NEWLINE );

    ## Script file
    open $filehandle, q{>}, $file_path
      or
      $log->logdie( q{Cannot write to '} . $file_path . q{' :} . $OS_ERROR . $NEWLINE );

    # Build bash shebang line
    build_shebang(
        {
            bash_bin_path => catfile( dirname( dirname( devnull() ) ), qw{ bin bash } ),
            filehandle    => $filehandle,
            invoke_login_shell => 1,
        }
    );

    ### Sbatch header
    ## Get parameters
    my $job_name        = $recipe_name . $UNDERSCORE . $directory_id;
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
                filehandle               => $filehandle,
                job_name                 => $job_name,
                memory_allocation        => $memory_allocation,
                process_time             => $process_time . q{:00:00},
                project_id               => $active_parameter_href->{project_id},
                slurm_quality_of_service => $slurm_quality_of_service,
                stderrfile_path          => $stderrfile_path,
                stdoutfile_path          => $stdoutfile_path,
            }
        );
    }

    ## Set shell attributes
    gnu_set(
        {
            filehandle   => $filehandle,
            set_errexit  => $set_errexit,
            set_nounset  => $set_nounset,
            set_pipefail => $set_pipefail,
        }
    );
    if ($ulimit_n) {
        gnu_ulimit(
            {
                filehandle     => $filehandle,
                max_open_files => $ulimit_n,
            }
        );
        say {$filehandle} $NEWLINE;
    }

    say {$filehandle} q{readonly PROGNAME=$(basename "$0")}, $NEWLINE;

    gnu_echo(
        {
            filehandle  => $filehandle,
            strings_ref => [q{Running on: $(hostname)}],
        }
    );
    say {$filehandle} $NEWLINE;

# Let the process sleep for a random couple of seconds (0-60) to avoid race conditions in mainly conda sourcing activate
    if ($sleep) {

        gnu_sleep(
            {
                filehandle       => $filehandle,
                seconds_to_sleep => int rand $MAX_SECONDS_TO_SLEEP,
            }
        );
        say {$filehandle} $NEWLINE;
    }
    if ( @{$source_environment_commands_ref}
        && $source_environment_commands_ref->[0] )
    {

        write_source_environment_command(
            {
                filehandle                      => $filehandle,
                source_environment_commands_ref => $source_environment_commands_ref,
            }
        );
    }

    # Not all recipes need a temporary directory
    if ( defined $temp_directory ) {

        say {$filehandle} q{## Create temporary directory};

        ## Double quote incoming variables in string
        my $temp_directory_quoted =
          quote_bash_variable( { string_with_variable_to_quote => $temp_directory, } );

        # Assign batch variable
        say {$filehandle} q{readonly TEMP_DIRECTORY=} . $temp_directory_quoted;

        # Update perl scalar to bash variable
        my $temp_directory_bash = q{"$TEMP_DIRECTORY"};

        gnu_mkdir(
            {
                filehandle       => $filehandle,
                indirectory_path => $temp_directory_bash,
                parents          => 1,
            }
        );
        say {$filehandle} $NEWLINE;

        create_housekeeping_function(
            {
                filehandle              => $filehandle,
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
                filehandle         => $filehandle,
                trap_function_call => q{previous_command="$BASH_COMMAND"},
                trap_signals_ref   => [qw{ DEBUG }],
            }
        );

        ## Create error handling function and trap
        create_error_trap_function(
            {
                filehandle              => $filehandle,
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

sub write_return_to_environment {

## Function : Return to MIP MAIN conda environment or default environment
## Returns  :
## Arguments: $active_parameter_href => The active parameters for this analysis hash {REF}
##          : $filehandle            => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $filehandle;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        filehandle => { required => 1, store => \$filehandle, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Active_parameter qw{ get_package_env_attributes };
    use MIP::Environment::Manager qw{ get_env_method_cmds };

    my @env_method_cmds;

    ## Get MIPs MAIN env
    my ( $env_name, $env_method ) = get_package_env_attributes(
        {
            load_env_href => $active_parameter_href->{load_env},
            package_name  => q{mip},
        }
    );

    ## Get env load command
    @env_method_cmds = get_env_method_cmds(
        {
            action     => q{load},
            env_name   => $env_name,
            env_method => $env_method,
        }
    );

    write_source_environment_command(
        {
            filehandle                      => $filehandle,
            source_environment_commands_ref => \@env_method_cmds,
        }
    );
    return;
}

sub write_return_to_conda_environment {

## Function : Return to main or default environment using conda
## Returns  :
## Arguments: $filehandle                           => Filehandle to write to
##          : $source_main_environment_commands_ref => Source main environment command {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $source_main_environment_commands_ref;

    my $tmpl = {
        filehandle                           => { required => 1, store => \$filehandle, },
        source_main_environment_commands_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$source_main_environment_commands_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Conda qw{ conda_deactivate };

    ## Return to main environment
    if ( @{$source_main_environment_commands_ref}
        && $source_main_environment_commands_ref->[0] )
    {

        write_source_environment_command(
            {
                filehandle => $filehandle,
                source_environment_commands_ref =>
                  \@{$source_main_environment_commands_ref},
            }
        );
    }
    else {
        ## Return to login shell environment

        say {$filehandle} q{## Deactivate environment};
        conda_deactivate(
            {
                filehandle => $filehandle,
            }
        );
        print {$filehandle} $NEWLINE;
    }
    return;
}

sub write_source_environment_command {

## Function : Write source environment commmands to filehandle
## Returns  :
## Arguments: $filehandle                      => Filehandle to write to
##          : $source_environment_commands_ref => Source environment command {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $source_environment_commands_ref;

    my $tmpl = {
        filehandle                      => { required => 1, store => \$filehandle, },
        source_environment_commands_ref => {
            default     => [],
            store       => \$source_environment_commands_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Unix::Write_to_file qw{ unix_write_to_file };

    if ( @{$source_environment_commands_ref} ) {

        say {$filehandle} q{## Activate environment};

        unix_write_to_file(
            {
                commands_ref => $source_environment_commands_ref,
                filehandle   => $filehandle,
            }
        );
        say {$filehandle} $NEWLINE;
    }
    return;
}

1;
