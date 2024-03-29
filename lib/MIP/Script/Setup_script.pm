package MIP::Script::Setup_script;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Basename qw{ dirname fileparse };
use File::Path qw{ make_path };
use File::Spec::Functions qw{ catdir catfile rootdir };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use Time::Piece;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants
  qw{ $CLOSE_PARENTHESIS $DOLLAR_SIGN $DOT $DOUBLE_QUOTE $EMPTY_STR $EQUALS $LOG_NAME $NEWLINE $OPEN_PARENTHESIS $SINGLE_QUOTE $SPACE $UNDERSCORE };
use MIP::Validate::Data qw{ %constraint };

BEGIN {

    use base qw{ Exporter };
    require Exporter;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      build_script_directories_and_paths
      check_script_file_path_exist
      create_script_error_trap
      create_script_temp_dir
      set_script_env_variables
      setup_script
      set_script_shell_attributes
    };
}

sub build_script_directories_and_paths {

## Function : Builds and makes recipe directories (info & data & script) and recipe script paths
## Returns  : $file_info_path, $file_path_prefix, $recipe_data_directory_path
## Arguments: $directory_id               => $sample id | $case_id
##          : $info_file_id               => Extra id for stderrr/stdout file
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

    ## Default(s)
    my $info_file_id;

    my $tmpl = {
        directory_id => {
            defined     => 1,
            required    => 1,
            store       => \$directory_id,
            strict_type => 1,
        },
        info_file_id => {
            default     => $arg_href->{info_file_id} ||= $arg_href->{directory_id},
            defined     => 1,
            required    => 1,
            store       => \$info_file_id,
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

    my $file_name_prefix = $recipe_name . $UNDERSCORE . $info_file_id . $DOT;

    ## Directory paths
    if ( not defined $recipe_data_directory_path ) {

        $recipe_data_directory_path = catdir( $outdata_dir, $directory_id, $recipe_directory );
    }
    my $recipe_info_directory_path   = catdir( $recipe_data_directory_path, q{info} );
    my $recipe_script_directory_path = catdir( $outscript_dir, $directory_id, $recipe_directory );

    ## Create directories
    make_path( $recipe_info_directory_path, $recipe_data_directory_path,
        $recipe_script_directory_path );

    ## File paths
    my $file_path_prefix = catfile( $recipe_script_directory_path, $file_name_prefix );
    my $file_info_path   = catfile( $recipe_info_directory_path,   $file_name_prefix );
    my $dry_run_file_path_prefix =
      catfile( $recipe_script_directory_path, q{dry_run} . $UNDERSCORE . $file_name_prefix );
    my $dry_run_file_info_path =
      catfile( $recipe_info_directory_path, q{dry_run} . $UNDERSCORE . $file_name_prefix );

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
    while ( $constraint{file_exists}->($file_path) ) {

        $file_name_version++;

        ## New file_path to test for existence
        $file_path = $file_path_prefix . $file_name_version . $file_path_suffix;
    }
    return ( $file_path, $file_name_version );
}

sub create_script_error_trap {

## Function : Create script error trap
## Returns  : 0 | undef
## Arguments: $error_trap              => Create error trap
##          : $filehandle              => Filehandle to write to
##          : $log_file_path           => Log file to write job_id progress to {REF}
##          : $job_ids_ref             => Job_ids to update status on {REF}
##          : $sacct_format_fields_ref => Format and fields of sacct output
##          : $submission_profile      => Process manager

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $error_trap;
    my $filehandle;
    my $log_file_path;
    my $job_ids_ref;
    my $sacct_format_fields_ref;
    my $submission_profile;

    my $tmpl = {
        error_trap => {
            allow       => [ 0, 1 ],
            default     => 1,
            store       => \$error_trap,
            strict_type => 1,
        },
        filehandle    => { required => 1,               store       => \$filehandle, },
        log_file_path => { store    => \$log_file_path, strict_type => 1, },
        job_ids_ref   => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$job_ids_ref,
            strict_type => 1,
        },
        sacct_format_fields_ref => {
            default     => [],
            store       => \$sacct_format_fields_ref,
            strict_type => 1,
        },
        submission_profile => {
            default     => q{slurm},
            store       => \$submission_profile,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Language::Shell qw{ create_error_trap_function enable_trap quote_bash_variable };

    return 0 if ( not $error_trap );

    my $enable_trap_function_call = quote_bash_variable(
        { string_with_variable_to_quote => q{previous_command=$BASH_COMMAND}, } );

    ## Create debug trap
    enable_trap(
        {
            filehandle         => $filehandle,
            trap_function_call => $enable_trap_function_call,
            trap_signals_ref   => [qw{ DEBUG }],
        }
    );

    ## Create error handling function and trap
    my $trap_function_name = q{error};
    my $previous_cmd =
      quote_bash_variable( { string_with_variable_to_quote => q{$previous_command}, } );
    my $previous_cmd_value = $SPACE . $DOUBLE_QUOTE . q{$?} . $DOUBLE_QUOTE . $CLOSE_PARENTHESIS;

    create_error_trap_function(
        {
            filehandle              => $filehandle,
            job_ids_ref             => $job_ids_ref,
            log_file_path           => $log_file_path,
            sacct_format_fields_ref => $sacct_format_fields_ref,
            submission_profile      => $submission_profile,
            trap_signals_ref        => [qw{ ERR }],
            trap_function_call      => $DOLLAR_SIGN
              . $OPEN_PARENTHESIS
              . $trap_function_name
              . $previous_cmd
              . $previous_cmd_value,
            trap_function_name => $trap_function_name,
        }
    );
    return;
}

sub create_script_temp_dir {

## Function : Create script temporary directory to use and set trap to remove it
## Returns  : $temp_directory_bash
## Arguments: $filehandle              => Filehandle to write to
##          : $log_file_path           => Log file to write job_id progress to {REF}
##          : $job_ids_ref             => Job_ids to update status on {REF}
##          : $sacct_format_fields_ref => Format and fields of sacct output
##          : $submission_profile      => Process manager
##          : $temp_directory          => Temporary directory for recipe

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $log_file_path;
    my $job_ids_ref;
    my $sacct_format_fields_ref;
    my $submission_profile;
    my $temp_directory;

    my $tmpl = {
        filehandle    => { required => 1,               store       => \$filehandle, },
        log_file_path => { store    => \$log_file_path, strict_type => 1, },
        job_ids_ref   => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$job_ids_ref,
            strict_type => 1,
        },
        sacct_format_fields_ref => {
            default     => [],
            store       => \$sacct_format_fields_ref,
            strict_type => 1,
        },
        submission_profile => {
            default     => q{slurm},
            store       => \$submission_profile,
            strict_type => 1,
        },
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Language::Shell qw{ create_housekeeping_function quote_bash_variable };
    use MIP::Program::Gnu::Coreutils qw{ gnu_mkdir };

    # Not all recipes need a temporary directory
    return 0 if ( not defined $temp_directory );

    say {$filehandle} q{## Create temporary directory};

    ## Double quote incoming variables in string
    my $temp_directory_quoted =
      quote_bash_variable( { string_with_variable_to_quote => $temp_directory, } );

    # Assign batch variable
    say {$filehandle} q{readonly TEMP_DIRECTORY} . $EQUALS . $temp_directory_quoted;

    # Update perl scalar to bash variable
    my $temp_directory_bash =
      quote_bash_variable( { string_with_variable_to_quote => $DOLLAR_SIGN . q{TEMP_DIRECTORY}, } );

    gnu_mkdir(
        {
            filehandle       => $filehandle,
            indirectory_path => $temp_directory_bash,
            parents          => 1,
        }
    );
    say {$filehandle} $NEWLINE;

    my $trap_function_name = q{finish};
    create_housekeeping_function(
        {
            filehandle              => $filehandle,
            job_ids_ref             => $job_ids_ref,
            log_file_path           => $log_file_path,
            remove_dir              => $temp_directory_bash,
            sacct_format_fields_ref => $sacct_format_fields_ref,
            submission_profile      => $submission_profile,
            trap_function_call      => $DOLLAR_SIGN
              . $OPEN_PARENTHESIS
              . $trap_function_name
              . $SPACE
              . $temp_directory_bash
              . $CLOSE_PARENTHESIS,
            trap_function_name => $trap_function_name,
            trap_signals_ref   => [qw{ EXIT TERM INT }],
        }
    );
    return $temp_directory_bash;
}

sub set_script_env_variables {

    ## Function : Set environment variables
    ## Returns  :
    ## Arguments: $filehandle          => Filehandle to write to
    ##          : $temp_directory_bash => Bash prepared temp directory
    ##          : $xdg_runtime_dir     => XDG_RUNTIME_DIR environment variable

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $temp_directory_bash;
    my $xdg_runtime_dir;

    my $tmpl = {
        filehandle          => { required => 1, store => \$filehandle, },
        temp_directory_bash => {
            store       => \$temp_directory_bash,
            strict_type => 1,
        },
        xdg_runtime_dir => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$xdg_runtime_dir,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    if ( $temp_directory_bash and $xdg_runtime_dir ) {

        say {$filehandle} q{XDG_RUNTIME_DIR} . $EQUALS . $temp_directory_bash;
    }
    return;
}

sub setup_script {

## Function : Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header.
## Returns  : $file_path, $file_info_path . $file_name_version
## Arguments: $active_parameter_href           => The active parameters for this analysis hash {REF}
##          : $core_number                     => Number of cores to allocate {Optional}
##          : $directory_id                    => $sample id | $case_id
##          : $error_trap                      => Error trap switch {Optional}
##          : $filehandle                      => filehandle to write to
##          : $gpu_number                      => Number of GPUs to use
##          : $info_file_id                    => Extra id for stderrr/stdout file
##          : $job_id_href                     => The job_id hash {REF}
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
##          : $source_environment_commands_ref => Source environment command {REF}
##          : $temp_directory                  => Temporary directory for recipe {Optional}
##          : $ulimit_n                        => Set ulimit -n for recipe {Optional}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $directory_id;
    my $filehandle;
    my $info_file_id;
    my $job_id_href;
    my $memory_allocation;
    my $recipe_data_directory_path;
    my $recipe_directory;
    my $recipe_name;
    my $source_environment_commands_ref;
    my $ulimit_n;

    ## Default(s)
    my $core_number;
    my $error_trap;
    my $gpu_number;
    my $outdata_dir;
    my $outscript_dir;
    my $process_time;
    my $set_errexit;
    my $set_nounset;
    my $set_pipefail;
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
            allow       => sub { $constraint{is_digit}->( $_[0] ) },
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
        error_trap => {
            allow       => [ 0, 1 ],
            default     => 1,
            store       => \$error_trap,
            strict_type => 1,
        },
        filehandle => { store => \$filehandle, },
        gpu_number => {
            allow       => [ undef, sub { $constraint{is_digit}->( $_[0] ) }, ],
            store       => \$gpu_number,
            strict_type => 1,
        },
        info_file_id => {
            store       => \$info_file_id,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        memory_allocation => {
            allow       => [ undef, sub { $constraint{is_digit}->( $_[0] ) }, ],
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
            allow       => sub { $constraint{is_digit}->( $_[0] ) },
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
            allow       => [ undef, sub { $constraint{is_digit}->( $_[0] ) }, ],
            store       => \$ulimit_n,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Environment::Manager qw{ write_source_environment_command };
    use MIP::Language::Shell qw{ build_shebang create_housekeeping_function };
    use MIP::Program::Gnu::Bash qw{ gnu_set gnu_ulimit };
    use MIP::Program::Gnu::Coreutils qw{ gnu_echo };
    use MIP::Program::Slurm qw{ slurm_build_sbatch_header };

    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Constants
    Readonly my %SUBMISSION_METHOD => ( slurm => q{sbatch}, );

    ## Unpack parameters
    my $submission_profile  = $active_parameter_href->{submission_profile};
    my @sacct_format_fields = @{ $active_parameter_href->{sacct_format_fields} };

    my $script_type = $SUBMISSION_METHOD{$submission_profile};

    ### Script names and directory creation
    ## File
    my $file_name_suffix = $DOT . q{sh};

    ( my ( $file_info_path, $file_path_prefix ), $recipe_data_directory_path ) =
      build_script_directories_and_paths(
        {
            directory_id               => $directory_id,
            info_file_id               => $info_file_id,
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

    _broadcast_script_paths(
        {
            file_path                  => $file_path,
            recipe_data_directory_path => $recipe_data_directory_path,
            recipe_name                => $recipe_name,
            script_type                => $script_type,
        }
    );

    ## Script file
    open $filehandle, q{>}, $file_path
      or $log->logdie( q{Cannot write to '} . $file_path . q{' :} . $OS_ERROR . $NEWLINE );

    # Build bash shebang line
    build_shebang(
        {
            bash_bin_path      => catfile( rootdir(), qw{ bin bash } ),
            filehandle         => $filehandle,
            invoke_login_shell => 0,
        }
    );
    print {$filehandle} $NEWLINE;

    ## SLURM specific headers and parameters
    if ( $submission_profile eq q{slurm} ) {

        ## Get parameters
        my $job_name        = $recipe_name . $UNDERSCORE . $directory_id;
        my $stderrfile_path = $file_info_path . $file_name_version . $DOT . q{stderr.txt};
        my $stdoutfile_path = $file_info_path . $file_name_version . $DOT . q{stdout.txt};

        slurm_build_sbatch_header(
            {
                core_number              => $core_number,
                email                    => $active_parameter_href->{email},
                email_types_ref          => $active_parameter_href->{email_types},
                filehandle               => $filehandle,
                gpu_number               => $gpu_number,
                job_name                 => $job_name,
                memory_allocation        => $memory_allocation,
                process_time             => $process_time . q{:00:00},
                project_id               => $active_parameter_href->{project_id},
                slurm_quality_of_service => $active_parameter_href->{slurm_quality_of_service},
                stderrfile_path          => $stderrfile_path,
                stdoutfile_path          => $stdoutfile_path,
            }
        );
    }

    ## Set shell attributes
    set_script_shell_attributes(
        {
            filehandle   => $filehandle,
            set_errexit  => $set_errexit,
            set_nounset  => $set_nounset,
            set_pipefail => $set_pipefail,
            ulimit_n     => $ulimit_n,
        }
    );

    write_source_environment_command(
        {
            filehandle                      => $filehandle,
            source_environment_commands_ref => $source_environment_commands_ref,
        }
    );

    my $temp_directory_bash = create_script_temp_dir(
        {
            filehandle              => $filehandle,
            job_ids_ref             => \@{ $job_id_href->{ALL}{ALL} },
            log_file_path           => $active_parameter_href->{log_file},
            sacct_format_fields_ref => \@sacct_format_fields,
            temp_directory          => $temp_directory,
        }
    );

    set_script_env_variables(
        {
            filehandle          => $filehandle,
            temp_directory_bash => $temp_directory_bash,
            xdg_runtime_dir     => 1,
        }
    );

    create_script_error_trap(
        {
            filehandle              => $filehandle,
            job_ids_ref             => \@{ $job_id_href->{ALL}{ALL} },
            log_file_path           => $active_parameter_href->{log_file},
            sacct_format_fields_ref => \@sacct_format_fields,
        }
    );

    # Return filen path, file path for stdout/stderr for QC check later
    return ( $file_path, $file_info_path . $file_name_version );
}

sub set_script_shell_attributes {

## Function : Set script shell attributes, such as "set -e", "ulimit" etc
## Returns  :
## Arguments: $filehandle   => filehandle to write to
##          : $set_errexit  => Bash set -e {Optional}
##          : $set_nounset  => Bash set -u {Optional}
##          : $set_pipefail => Pipe fail switch {Optional}
##          : $ulimit_n     => Set ulimit -n for recipe {Optional}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $ulimit_n;

    ## Default(s)
    my $set_errexit;
    my $set_nounset;
    my $set_pipefail;

    my $tmpl = {
        filehandle  => { store => \$filehandle, },
        set_errexit => {
            allow       => [ 0, 1 ],
            default     => 1,
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
            default     => 1,
            store       => \$set_pipefail,
            strict_type => 1,
        },
        ulimit_n => {
            allow       => [ undef, qr/ \A \d+ \z /xms ],
            store       => \$ulimit_n,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Language::Shell qw{ log_host_name };
    use MIP::Program::Gnu::Bash qw{ gnu_set gnu_ulimit };

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

    log_host_name( { filehandle => $filehandle, } );
    say {$filehandle} $NEWLINE;

    return;
}

sub _broadcast_script_paths {

## Function : Broadcast script and data directory paths
## Returns  :
## Arguments: $file_path                  => Script file path
##          : $recipe_data_directory_path => Recipe data file output directory
##          : $recipe_name                => Name of recipe
##          : $script_type                => Dry run or sharp

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path;
    my $recipe_data_directory_path;
    my $recipe_name;
    my $script_type;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        recipe_data_directory_path => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_data_directory_path,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        script_type => {
            defined     => 1,
            required    => 1,
            store       => \$script_type,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $log = Log::Log4perl->get_logger($LOG_NAME);

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
    return;
}

1;
