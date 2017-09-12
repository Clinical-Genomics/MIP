package MIP::Script::Setup_script;

#### Copyright 2017 Henrik Stranneheim

use strict;
use warnings;
use warnings qw(FATAL utf8);
use utf8;    #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use autodie;
use English qw(-no_match_vars);
use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;    #Do not convert to lower case
use File::Basename qw(dirname);
use File::Spec::Functions qw(catdir catfile devnull);
use File::Path qw(make_path);

##CPANM
use Readonly;

BEGIN {

    use base qw(Exporter);
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(setup_script);
}

##Constant
Readonly my $EMPTY_STR => q{};
Readonly my $SPACE     => q{ };

sub setup_script {

##setup_script

##Function : Creates program directories (info & programData & programScript), program script filenames and writes sbatch header.
##Returns  : Path to stdout
##Arguments: $active_parameter_href, $job_id_href, $FILEHANDLE, $directory_id, $program_directory, $program_name, $call_type, $outdata_dir, $outscript_dir, $temp_directory, $email_types_ref, $source_environment_commands_ref, $slurm_quality_of_service, $core_number, $process_time, $error_trap, $set_errexit, $set_nounset, $set_pipefail, $sleep
##         : $active_parameter_href           => The active parameters for this analysis hash {REF}
##         : $job_id_href                     => The job_id hash {REF}
##         : $FILEHANDLE                      => FILEHANDLE to write to
##         : $directory_id                    => $samplID|$family_id
##         : $program_directory               => Builds from $directory_id/$outaligner_dir
##         : $program_name                    => Assigns filename to sbatch script
##         : $call_type                       => SNV,INDEL or BOTH
##         : $source_environment_commands_ref => Source environment command {REF}
##         : $outdata_dir                     => The MIP out data directory {Optional}
##         : $outscript_dir                   => The MIP out script directory {Optional}
##         : $temp_directory                  => Temporary directory for program {Optional}
##         : $email_types_ref                 => The email type
##         : $slurm_quality_of_service        => SLURM quality of service priority {Optional}
##         : $core_number                     => The number of cores to allocate {Optional}
##         : $process_time                    => Allowed process time (Hours) {Optional}
##         : $error_trap                      => Error trap switch {Optional}
##         : $set_errexit                     => Bash set -e {Optional}
##         : $set_nounset                     => BAsh set -u {Optional}
##         : $set_pipefail                    => Pipe fail switch {Optional}
##         : $sleep                           => Sleep for X seconds {Optional}

    my ($arg_href) = @_;

    ## Default(s)
    my $outdata_dir;
    my $outscript_dir;
    my $temp_directory;
    my $email_types_ref;
    my $source_environment_commands_ref;
    my $slurm_quality_of_service;
    my $core_number;
    my $process_time;
    my $set_errexit;
    my $set_nounset;
    my $set_pipefail;
    my $error_trap;
    my $sleep;

    if ( defined( $arg_href->{call_type} ) ) {

        $arg_href->{call_type} = q{_} . $arg_href->{call_type};
    }
    $arg_href->{call_type} //= $EMPTY_STR;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $job_id_href;
    my $FILEHANDLE;
    my $directory_id;
    my $program_directory;
    my $program_name;
    my $call_type;

    use MIP::Check::Parameter qw(check_allowed_array_values);

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href
        },
        FILEHANDLE   => { store => \$FILEHANDLE },
        directory_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$directory_id
        },
        program_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_directory
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name
        },
        call_type   => { strict_type => 1, store => \$call_type },
        outdata_dir => {
            default     => $arg_href->{active_parameter_href}{outdata_dir},
            strict_type => 1,
            store       => \$outdata_dir
        },
        outscript_dir => {
            default     => $arg_href->{active_parameter_href}{outscript_dir},
            strict_type => 1,
            store       => \$outscript_dir
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory
        },
        email_types_ref => {
            default => $arg_href->{active_parameter_href}{email_types},
            allow   => [
                sub {
                    check_allowed_array_values(
                        {
                            allowed_values_ref =>
                              [qw(NONE BEGIN END FAIL REQUEUE ALL)],
                            values_ref => $arg_href->{email_types_ref},
                        }
                    );
                }
            ],
            strict_type => 1,
            store       => \$email_types_ref
        },
        source_environment_commands_ref => {
            default =>
              $arg_href->{active_parameter_href}{source_environment_commands},
            strict_type => 1,
            store       => \$source_environment_commands_ref
        },
        core_number => {
            default     => 1,
            allow       => qr/^\d+$/,
            strict_type => 1,
            store       => \$core_number
        },
        process_time => {
            default     => 1,
            allow       => qr/^\d+$/,
            strict_type => 1,
            store       => \$process_time
        },
        slurm_quality_of_service => {
            default =>
              $arg_href->{active_parameter_href}{slurm_quality_of_service},
            allow       => [qw(low high normal)],
            strict_type => 1,
            store       => \$slurm_quality_of_service
        },
        set_nounset => {
            default     => $arg_href->{active_parameter_href}{bash_set_nounset},
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$set_nounset
        },
        set_errexit => {
            default     => $arg_href->{active_parameter_href}{bash_set_errexit},
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$set_errexit
        },
        set_pipefail => {
            default => $arg_href->{active_parameter_href}{bash_set_pipefail},
            allow   => [ 0, 1 ],
            strict_type => 1,
            store       => \$set_pipefail
        },
        error_trap => {
            default     => 1,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$error_trap
        },
        sleep => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$sleep
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    use MIP::Language::Shell
      qw(build_shebang create_housekeeping_function create_error_trap_function enable_trap);
    use MIP::Workloadmanager::Slurm qw(slurm_build_sbatch_header);
    use MIP::Gnu::Bash qw(gnu_set);
    use MIP::Gnu::Coreutils qw(gnu_echo gnu_mkdir gnu_sleep);
    use MIP::Check::Path qw(check_file_version_exist);
    use MIP::Language::Shell qw{quote_bash_variable};

    ##Constants
    Readonly my $MAX_SECONDS_TO_SLEEP => 60;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger('MIP');

    ### Sbatch script names and directory creation
    my $file_name_end = q{.sh};

    # The sbatch script - to be created filename
    my $file_name;
    my $file_name_tracker;

    ## Directories
    my $program_data_directory =
      catdir( $outdata_dir, $directory_id, $program_directory );
    my $program_info_directory = catdir( $program_data_directory, 'info' );
    my $program_script_directory =
      catdir( $outscript_dir, $directory_id, $program_directory );

    ## Create directories
    make_path( $program_info_directory, $program_data_directory,
        $program_script_directory );

    ## File
    my $file_name_suffix =
      $program_name . q{_} . $directory_id . $call_type . q{.};
    my $file_name_path =
      catfile( $program_script_directory, $file_name_suffix );
    my $dry_run_file_name_path =
      catfile( $program_script_directory, q{dry_run_} . $file_name_suffix );
    my $file_info_path = catfile( $program_info_directory, $file_name_suffix );
    my $dry_run_file_info_path =
      catfile( $program_info_directory, q{dry_run_} . $file_name_suffix );

    ## Set paths depending on dry run or not
    if (   ( $active_parameter_href->{ 'p' . $program_name } == 1 )
        && ( !$active_parameter_href->{dry_run_all} ) )
    {

        $file_name = $file_name_path;
    }
    elsif ( $active_parameter_href->{ 'p' . $program_name } == 2 )
    {    #Dry run single program

        $file_name      = $dry_run_file_name_path;
        $file_info_path = $dry_run_file_info_path;
        $log->info( q{Dry run:}, "\n" );
    }
    else {    #Dry run

        $file_name      = $dry_run_file_name_path;
        $file_info_path = $dry_run_file_info_path;
        $log->info( q{Dry run:}, "\n" );
    }

    ## Check if a file with with a filename consisting of
    ## $file_path_prefix_ref.$file_counter.$file_path_suffix_ref exist
    ( $file_name, $file_name_tracker ) = check_file_version_exist(
        {
            file_path_prefix_ref => \$file_name,
            file_path_suffix_ref => \$file_name_end,
        }
    );

###Info and Log
    $log->info( q{Creating sbatch script for }
          . $program_name
          . q{ and writing script file(s) to: }
          . $file_name
          . "\n" );
    $log->info( q{Sbatch script }
          . $program_name
          . q{ data files will be written to: }
          . $program_data_directory
          . "\n" );

    ## Script file
    open $FILEHANDLE, q{>}, $file_name
      or $log->logdie(
        q{Can't write to '} . $file_name . q{' :} . $OS_ERROR . "\n" );

    # Build bash shebang line
    build_shebang(
        {
            FILEHANDLE => $FILEHANDLE,
            bash_bin_path =>
              catfile( dirname( dirname( devnull() ) ), qw(bin bash) ),
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
    my $job_name        = $program_name . q{_} . $directory_id . $call_type;
    my $stderrfile_path = $file_info_path . $file_name_tracker . q{.stderr.txt};
    my $stdoutfile_path = $file_info_path . $file_name_tracker . q{.stdout.txt};

    my @sbatch_headers = slurm_build_sbatch_header(
        {
            project_id               => $active_parameter_href->{project_id},
            core_number              => $core_number,
            process_time             => $process_time . q{:00:00},
            slurm_quality_of_service => $slurm_quality_of_service,
            job_name                 => $job_name,
            stderrfile_path          => $stderrfile_path,
            stdoutfile_path          => $stdoutfile_path,
            email                    => $active_parameter_href->{email},
            email_types_ref          => $email_types_ref,
            FILEHANDLE               => $FILEHANDLE,
        }
    );

    say {$FILEHANDLE} q{readonly PROGNAME=$(basename "$0")}, "\n";

    gnu_echo(
        {
            strings_ref => [q{Running on: $(hostname)}],
            FILEHANDLE  => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} "\n";

# Let the process sleep for a random couple of seconds (0-60) to avoid race conditions in mainly conda sourcing activate
    if ($sleep) {

        gnu_sleep(
            {
                seconds_to_sleep => int rand $MAX_SECONDS_TO_SLEEP,
                FILEHANDLE       => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} "\n";
    }
    if (   ($source_environment_commands_ref)
        && ( @{$source_environment_commands_ref} ) )
    {

        say {$FILEHANDLE} q{##Activate environment};
        say {$FILEHANDLE} join( $SPACE, @{$source_environment_commands_ref} ),
          "\n";
    }

    # Not all programs need a temporary directory
    if ( defined $temp_directory ) {

        say {$FILEHANDLE} q{## Create temporary directory};

        ## Double quote incoming variables in string
        my $temp_directory_quoted =
          quote_bash_variable(
            { string_with_variable_to_quote => $temp_directory, } );

        # Assign batch variable
        say {$FILEHANDLE} q{readonly TEMP_DIRECTORY=} . $temp_directory_quoted;

        # Update perl scalar to bash variable
        my $temp_directory_bash = q{"$TEMP_DIRECTORY"};

        gnu_mkdir(
            {
                indirectory_path => $temp_directory_bash,
                parents          => 1,
                FILEHANDLE       => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} "\n";

        create_housekeeping_function(
            {
                job_ids_ref => \@{ $job_id_href->{PAN}{PAN} },
                sacct_format_fields_ref =>
                  \@{ $active_parameter_href->{sacct_format_fields} },
                log_file_path      => $active_parameter_href->{log_file},
                FILEHANDLE         => $FILEHANDLE,
                remove_dir         => $temp_directory_bash,
                trap_signals_ref   => [qw(EXIT TERM INT)],
                trap_function_name => 'finish',
                trap_function_call => q{$(finish }
                  . $temp_directory_bash . q{)},
            }
        );
    }

    if ($error_trap) {

        ## Create debug trap
        enable_trap(
            {
                FILEHANDLE         => $FILEHANDLE,
                trap_signals_ref   => ['DEBUG'],
                trap_function_call => q?previous_command="$BASH_COMMAND"?,
            }
        );

        ## Create error handling function and trap
        create_error_trap_function(
            {
                job_ids_ref => \@{ $job_id_href->{PAN}{PAN} },
                sacct_format_fields_ref =>
                  \@{ $active_parameter_href->{sacct_format_fields} },
                log_file_path      => $active_parameter_href->{log_file},
                FILEHANDLE         => $FILEHANDLE,
                trap_signals_ref   => ['ERR'],
                trap_function_name => 'error',
                trap_function_call => q{$(error "$previous_command" "$?")},
            }
        );
    }

    # Return filen name, file path for stdout/stderr for QC check later
    return ( $file_name, $file_info_path . $file_name_tracker );
}

1;
