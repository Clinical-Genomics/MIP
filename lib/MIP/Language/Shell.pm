package MIP::Language::Shell;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use FindBin qw{ $Bin };
use File::Basename qw{ dirname fileparse };
use File::Spec::Functions qw{ catfile catdir devnull };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use Time::Piece;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.09;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      build_shebang
      check_exist_and_move_file
      clear_trap
      create_error_trap_function
      create_housekeeping_function
      enable_trap
      quote_bash_variable
      track_progress
    };
}

## Constants
Readonly my $AMPERSAND    => q{&};
Readonly my $COMMA        => q{,};
Readonly my $DOT          => q{.};
Readonly my $NEWLINE      => qq{\n};
Readonly my $PIPE         => q{|};
Readonly my $SINGLE_QUOTE => q{'};
Readonly my $SPACE        => q{ };
Readonly my $TAB          => qq{\t};
Readonly my $UNDERSCORE   => q{_};

sub build_shebang {

## Function : Build bash shebang line. Returns @commands or writes to already opened filehandle.
## Returns  : @commands
## Arguments: $FILEHANDLE         => Filehandle to write to
##          : $bash_bin_path      => Location of bash bin
##          : $invoke_login_shell => Invoked as a login shell (-l). Reinitilize bashrc and bash_profile
##          : $separator          => Separator to use when writing

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;

    ## Default(s)
    my $bash_bin_path;
    my $invoke_login_shell;
    my $separator;

    my $tmpl = {
        bash_bin_path => {
            allow   => qr/ ^\S+$ /sxm,
            default => catfile( dirname( dirname( devnull() ) ), qw{ usr bin env bash } ),
            store   => \$bash_bin_path,
            strict_type => 1,
        },
        FILEHANDLE         => { store => \$FILEHANDLE, },
        invoke_login_shell => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$invoke_login_shell,
            strict_type => 1,
        },
        separator => {
            default     => $NEWLINE,
            store       => \$separator,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Unix::Write_to_file qw{ unix_write_to_file };

    ## Build shebang
    # Stores commands depending on input parameters
    my @commands = ( q{#!} . $bash_bin_path );

    ##Invoke as login shell
    if ($invoke_login_shell) {

        $commands[0] .= $SPACE . q{--login};
    }

    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
            separator    => $separator,
        }
    );
    return @commands;
}

sub create_housekeeping_function {

## Function : Create housekeeping function which removes entire directory when finished
## Returns  :
## Arguments: $FILEHANDLE              => Filehandle to write to
##          : $job_ids_ref             => Job ids
##          : $log_file_path           => Log file to write job_id progress to {REF}
##          : $remove_dir              => Directory to remove when caught by trap function
##          : $sacct_format_fields_ref => Format and fields of sacct output
##          : $trap_function_call      => Trap function call
##          : $trap_function_name      => The trap function argument
##          : $trap_signals_ref        => Array with signals to enable trap for {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $job_ids_ref;
    my $log_file_path;
    my $remove_dir;
    my $sacct_format_fields_ref;

    ## Default(s)
    my $trap_function_call;
    my $trap_function_name;
    my $trap_signals_ref;

    my $tmpl = {
        FILEHANDLE => { required => 1, store => \$FILEHANDLE, },
        job_ids_ref => { default => [], store => \$job_ids_ref, strict_type => 1, },
        log_file_path           => { store => \$log_file_path, strict_type => 1, },
        remove_dir              => { store => \$remove_dir,    strict_type => 1, },
        sacct_format_fields_ref => {
            default     => [],
            store       => \$sacct_format_fields_ref,
            strict_type => 1,
        },
        trap_function_call => {
            default => q{$(}
              . $arg_href->{trap_function_name}
              . $SPACE
              . $arg_href->{remove_dir} . q{)},
            store       => \$trap_function_call,
            strict_type => 1,
        },
        trap_function_name => {
            default     => q{finish},
            strict_type => 1,
            store       => \$trap_function_name,
        },
        trap_signals_ref => {
            default     => [qw{ EXIT TERM INT }],
            strict_type => 1,
            store       => \$trap_signals_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Gnu::Coreutils qw{ gnu_rm };

    ## Create housekeeping function and trap
    say {$FILEHANDLE} $trap_function_name . q?() {?, $NEWLINE;

    if ( defined $remove_dir && $remove_dir ) {

        say   {$FILEHANDLE} $TAB . q{local directory="$1"};
        say   {$FILEHANDLE} $TAB . q{## Perform exit housekeeping};
        print {$FILEHANDLE} $TAB;

        gnu_rm(
            {
                FILEHANDLE  => $FILEHANDLE,
                force       => 1,
                infile_path => q{"$directory"},
                recursive   => 1,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }
    if (   defined $job_ids_ref
        && @{$job_ids_ref}
        && defined $log_file_path
        && $log_file_path )
    {

        ## Output SLURM info on each job via sacct command
        ## and write to log file(.status)
        track_progress(
            {
                FILEHANDLE              => $FILEHANDLE,
                job_ids_ref             => \@{$job_ids_ref},
                log_file_path           => $log_file_path,
                sacct_format_fields_ref => \@{$sacct_format_fields_ref},
            }
        );
    }

    ## End of trap function
    say {$FILEHANDLE} q?}?;

    ## Enable trap function with trap signal(s)
    enable_trap(
        {
            FILEHANDLE         => $FILEHANDLE,
            trap_signals_ref   => \@{$trap_signals_ref},
            trap_function_call => $trap_function_call,
        }
    );
    return;
}

sub create_error_trap_function {

## Function : Create error handling function and trap
## Returns  :
## Arguments: $FILEHANDLE              => Filehandle to write to
##          : $job_ids_ref             => Job ids
##          : $log_file_path           => Log file to write job_id progress to {REF}
##          : $sacct_format_fields_ref => Format and fields of sacct output
##          : $trap_function_call      => Trap function call
##          : $trap_function_name      => The trap function argument
##          : $trap_signals_ref        => Array with signals to enable trap for {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $job_ids_ref;
    my $log_file_path;
    my $sacct_format_fields_ref;
    my $trap_function_call;

    ## Default(s)
    my $trap_function_name;
    my $trap_signals_ref;

    my $tmpl = {
        FILEHANDLE => { required => 1, store => \$FILEHANDLE, },
        job_ids_ref   => { default     => [], store => \$job_ids_ref, strict_type => 1, },
        log_file_path => { strict_type => 1,  store => \$log_file_path, },
        sacct_format_fields_ref => {
            default     => [],
            store       => \$sacct_format_fields_ref,
            strict_type => 1,
        },
        trap_function_call => {
            default     => q{$(error "$previous_command" "$?")},
            store       => \$trap_function_call,
            strict_type => 1,
        },
        trap_function_name => {
            default     => q{error},
            store       => \$trap_function_name,
            strict_type => 1,
        },
        trap_signals_ref => {
            default     => [qw{ ERR }],
            store       => \$trap_signals_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Create error handling function and trap
    say {$FILEHANDLE} $trap_function_name . q?() {?,    $NEWLINE;
    say {$FILEHANDLE} $TAB . q{local program="$1"},     $NEWLINE;
    say {$FILEHANDLE} $TAB . q{local return_code="$2"}, $NEWLINE;

    if (   defined $job_ids_ref
        && @{$job_ids_ref}
        && defined $log_file_path
        && $log_file_path )
    {

        ## Output SLURM info on each job via sacct command
        ## and write to log file(.status)
        track_progress(
            {
                FILEHANDLE              => $FILEHANDLE,
                job_ids_ref             => \@{$job_ids_ref},
                log_file_path           => $log_file_path,
                sacct_format_fields_ref => \@{$sacct_format_fields_ref},
            }
        );
    }

    say {$FILEHANDLE} $TAB . q{## Display error message and exit};
    say {$FILEHANDLE} $TAB
      . q?echo "${program}: ${return_code}: Unknown Error - ExitCode=$return_code" 1>&2?;
    say {$FILEHANDLE} $TAB . q{exit 1};
    say {$FILEHANDLE} q?}?;

    ## Enable trap function with trap signal(s)
    enable_trap(
        {
            FILEHANDLE         => $FILEHANDLE,
            trap_function_call => $trap_function_call,
            trap_signals_ref   => \@{$trap_signals_ref},
        }
    );
    return;
}

sub clear_trap {

## Function : Clear trap for signal(s), e.g. in exome analysis since the might be no variants in MT or Y contigs. This will cause premature exit from sbatch
## Returns  :
## Arguments: $FILEHANDLE       => The FILEHANDLE to write to
##          : $trap_signals_ref => Array with signals to clear trap for {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;

    ## Default(s)
    my $trap_signals_ref;

    my $tmpl = {
        FILEHANDLE       => { required => 1, store => \$FILEHANDLE, },
        trap_signals_ref => {
            default     => [qw{ ERR }],
            store       => \$trap_signals_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Gnu::Bash qw{ gnu_trap };

    ## Clear trap for signal ERR
    say {$FILEHANDLE} $NEWLINE . q{## Clear trap for signal(s) } . join $SPACE,
      @{$trap_signals_ref};

    gnu_trap(
        {
            FILEHANDLE         => $FILEHANDLE,
            trap_function_call => q{-},
            trap_signals_ref   => $trap_signals_ref,
        }
    );
    gnu_trap( { FILEHANDLE => $FILEHANDLE, } );
    say {$FILEHANDLE} $NEWLINE;
    return;
}

sub enable_trap {

## Function : Enable trap function with trap signal(s).
## Returns  :
## Arguments: $FILEHANDLE         => The FILEHANDLE to write to
##          : $trap_function_call => The trap function argument
##          : $trap_signals_ref   => Array with signals to enable trap for {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;

    ## Default(s)
    my $trap_function_call;
    my $trap_signals_ref;

    my $tmpl = {
        FILEHANDLE         => { required => 1, store => \$FILEHANDLE, },
        trap_function_call => {
            default     => q{error},
            store       => \$trap_function_call,
            strict_type => 1,
        },
        trap_signals_ref => {
            default     => [qw{ ERR }],
            store       => \$trap_signals_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Gnu::Bash qw{ gnu_trap };

    say {$FILEHANDLE} $NEWLINE . q{## Enable trap for signal(s) } . join $SPACE,
      @{$trap_signals_ref};

    gnu_trap(
        {
            trap_signals_ref   => $trap_signals_ref,
            trap_function_call => $trap_function_call,
            FILEHANDLE         => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;
    return;
}

sub track_progress {

## Function : Output SLURM info on each job via sacct command and write to log file(.status)
## Returns  :
## Arguments: $FILEHANDLE              => Sbatch filehandle to write to
##          : $job_ids_ref             => Job ids
##          : $log_file_path           => The log file {REF}
##          : $sacct_format_fields_ref => Format and fields of sacct output

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $job_ids_ref;
    my $log_file_path;

    ## Default(s)
    my $sacct_format_fields_ref;

    my $tmpl = {
        FILEHANDLE    => { store   => \$FILEHANDLE, },
        job_ids_ref   => { default => [], store => \$job_ids_ref, strict_type => 1, },
        log_file_path => { store   => \$log_file_path, strict_type => 1, },
        sacct_format_fields_ref => {
            default => [
                qw{
                  jobid jobname%50 account partition alloccpus TotalCPU elapsed start end state exitcode }
            ],
            store       => \$sacct_format_fields_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Workloadmanager::Slurm qw{ slurm_sacct slurm_reformat_sacct_output };

    if ( @{$job_ids_ref} ) {

        ## Copy array
        my @reformat_sacct_headers = @{$sacct_format_fields_ref};

        ## Remove "%digits" from headers
      HEADER_ELEMENT:
        foreach my $element (@reformat_sacct_headers) {

            $element =~ s/%\d+//gsxm;
        }
        my @commands = slurm_sacct(
            {
                fields_format_ref => \@{$sacct_format_fields_ref},
                job_ids_ref       => \@{$job_ids_ref},
            }
        );

        slurm_reformat_sacct_output(
            {
                commands_ref               => \@commands,
                FILEHANDLE                 => $FILEHANDLE,
                log_file_path              => $log_file_path,
                reformat_sacct_headers_ref => \@reformat_sacct_headers,
            }
        );
    }
    return;
}

sub quote_bash_variable {

## Function : Double quote incoming variables in string
## Returns  : "String with double quoted variables"
## Arguments: $string_with_variable_to_quote => String to find variables in

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $string_with_variable_to_quote;

    my $tmpl = {
        string_with_variable_to_quote => {
            defined     => 1,
            required    => 1,
            store       => \$string_with_variable_to_quote,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $string_with_double_quoted_variable = $string_with_variable_to_quote;

    ## Find and double quote variables in string
    $string_with_double_quoted_variable =~ s/(\$\w+)/"$1"/gsxm;

    return $string_with_double_quoted_variable;
}

sub check_exist_and_move_file {

## Function : Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
## Returns  :
## Arguments: $FILEHANDLE          => FILEHANDLE to write to
##          : $intended_file_path  => Path to file to check for existence {REF}
##          : $temporary_file_path => File that has been created {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $intended_file_path;
    my $temporary_file_path;

    my $tmpl = {
        FILEHANDLE         => { required => 1, store => \$FILEHANDLE, },
        intended_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$intended_file_path,
            strict_type => 1,
        },
        temporary_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$temporary_file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Gnu::Coreutils qw{ gnu_rm gnu_mv};

    ## Check file exists and is larger than 0
    print {$FILEHANDLE} q{[ -s } . $intended_file_path . q{ ]} . $SPACE;
    print {$FILEHANDLE} $AMPERSAND x 2 . $SPACE;

    ## If other processes already has created file, remove temp file
    gnu_rm(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $temporary_file_path,
        }
    );
    ## File has not been created by other processes
    print {$FILEHANDLE} $PIPE x 2 . $SPACE;

    gnu_mv(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $temporary_file_path,
            outfile_path => $intended_file_path,
        }
    );
    say {$FILEHANDLE} $NEWLINE;
    return;
}

1;
