package MIP::Recipes::Install::Post_installation;

use 5.026;
use Carp;
use Cwd qw{ abs_path };
use English qw{ -no_match_vars };
use File::Basename qw{ fileparse };
use File::Spec::Functions qw{ catfile catdir };
use FindBin qw{ $Bin };
use List::Util qw{ any };
use Params::Check qw{ allow check last_error };
use charnames qw{ :full :short };
use open qw{ :encoding(UTF-8) :std };
use strict;
use Time::Piece;
use utf8;
use warnings qw{ FATAL utf8 };
use warnings;

## CPANM
use autodie qw{ :all };
use List::MoreUtils qw{ natatime };
use Readonly;

## MIPs lib/
use MIP::Constants
  qw{ $DOUBLE_QUOTE $NEWLINE $LOG_NAME $SEMICOLON $SINGLE_QUOTE $SPACE $TAB };
use MIP::Get::Parameter qw{ get_env_method_cmds };
use MIP::Gnu::Bash qw{ gnu_set };
use MIP::Gnu::Coreutils qw{ gnu_cp gnu_echo gnu_printf gnu_rm };
use MIP::Log::MIP_log4perl qw{ retrieve_log };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.09;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      build_perl_program_check_command
      check_mip_installation
      check_program_installations
    };
}

sub check_mip_installation {

## Function : Write installation check oneliner to open filehandle
## Returns  :
## Arguments: $active_parameter_href => Active parameter hash {REF}
##          : $filehandle            => Open filehandle

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
        filehandle => {
            required => 1,
            store    => \$filehandle,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Io::Read qw{ read_from_file };

    my %program_test_cmds = read_from_file(
        {
            format => q{yaml},
            path   => $active_parameter_href->{program_test_file},
        }
    );

    ## Get the programs that mip has tried to install
    my @programs_to_test = (
        keys %{ $active_parameter_href->{conda} },
        keys %{ $active_parameter_href->{pip} },
        keys %{ $active_parameter_href->{shell} },
        keys %{ $active_parameter_href->{singularity} },
    );

    check_program_installations(
        {
            env_name                  => $active_parameter_href->{environment_name},
            filehandle                => $filehandle,
            programs_ref              => \@programs_to_test,
            program_test_command_href => \%program_test_cmds,
        }
    );
    return;
}

sub check_program_installations {

## Function : Write installation check oneliner to open filehandle
## Returns  :
## Arguments: $env_name                   => Program environment name
##          : $filehandle                 => open filehandle
##          : $programs_ref               => Programs to check {REF}
##          : $program_test_command_href  => Hash with test commands {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $env_name;
    my $filehandle;
    my $programs_ref;
    my $program_test_command_href;

    my $tmpl = {
        env_name => {
            defined     => 1,
            required    => 1,
            store       => \$env_name,
            strict_type => 1,
        },
        filehandle => {
            required => 1,
            store    => \$filehandle,
        },
        programs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$programs_ref,
            strict_type => 1,
        },
        program_test_command_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$program_test_command_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    Readonly my $COLUMNS      => 4;
    Readonly my $COLUMN_WIDTH => 25;
    Readonly my $HASH_SIGN    => q{#};

    ## Retrieve logger object
    my $log = retrieve_log( { log_name => $LOG_NAME, } );

    $log->info(qq{Writing tests for programs installed in environment: $env_name});

    ## Sort program list for readability
    my @sorted_programs = sort { lc $a cmp lc $b } @{$programs_ref};

    ## Get environment commands
    my @env_method_load_cmds = get_env_method_cmds(
        {
            action     => q{load},
            env_method => q{conda},
            env_name   => $env_name,
        }
    );

    my @env_method_unload_cmds = get_env_method_cmds(
        {
            action     => q{unload},
            env_method => q{conda},
            env_name   => $env_name,
        }
    );

    ## Construct perl test for programs
    my @perl_commands = build_perl_program_check_command(
        {
            programs_ref              => \@sorted_programs,
            program_test_command_href => $program_test_command_href,
        },
    );

    ## Return if there are no tests to be written;
    if ( not @perl_commands ) {
        $log->info(qq{No tests available for programs in environment: $env_name});
        return;
    }

    say {$filehandle} qq{## Testing programs installed in $env_name};

    ## Build header string to echo
    my @header = ( q{\n} . ( $HASH_SIGN x $COLUMN_WIDTH x $COLUMNS ) . q{\n\n} );
    push @header, (q{\tMIP has attempted to install the following programs\n});
    push @header, ( q{\tin environment: } . $env_name . q{\n\n} );

    gnu_echo(
        {
            enable_interpretation => 1,
            filehandle            => $filehandle,
            strings_ref           => \@header,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Print programs selected for installation in sets
    my $program_iterator = natatime $COLUMNS, @sorted_programs;

    ## Create field layout for printf
    my $field = q{"\t} . ( ( q{%-} . $COLUMN_WIDTH . q{s} ) x $COLUMNS ) . q{\n"};

  PRINT_SET:
    while ( my @print_set = $program_iterator->() ) {
        ## Print programs according to field
        my $format_string =
            $field
          . $SPACE
          . $DOUBLE_QUOTE
          . join( $DOUBLE_QUOTE . $SPACE . $DOUBLE_QUOTE, @print_set )
          . $DOUBLE_QUOTE;
        gnu_printf(
            {
                filehandle    => $filehandle,
                format_string => $format_string,
            }
        );
        print {$filehandle} $NEWLINE;
    }

    gnu_echo(
        {
            enable_interpretation => 1,
            filehandle            => $filehandle,
            strings_ref           => [q{\n\tTesting installation\n}],
        }
    );
    say {$filehandle} $NEWLINE;

    ## Load env
    say   {$filehandle} qq{## Load environment: $env_name};
    say   {$filehandle} join $SPACE, @env_method_load_cmds;
    print {$filehandle} $NEWLINE;

    ## Create success and fail case
    my $success_message =
      q{\n\tAll programs were succesfully installed in: } . $env_name . q{\n};
    my $fail_message =
      q{\n\tMIP failed to install some programs in: } . $env_name . q{\n};

    my $success_echo = join $SPACE,
      gnu_echo(
        {
            enable_interpretation => 1,
            strings_ref           => [$success_message],
        }
      );
    my $fail_echo = join $SPACE,
      gnu_echo(
        {
            enable_interpretation => 1,
            strings_ref           => [$fail_message],
        }
      );

    my $success_case = qq?&& { $success_echo; }?;
    my $fail_case    = qq?|| { $fail_echo; }?;

    ## Enabling querying of $?
    gnu_set(
        {
            filehandle    => $filehandle,
            unset_errexit => 1,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Write test oneliner
    say {$filehandle} q{## Test programs and capture outcome in bash variable};
    print {$filehandle} join $SPACE, @perl_commands;
    say {$filehandle} $SPACE . $success_case . $SPACE . $fail_case . $NEWLINE;

    ## Restore errexit
    gnu_set(
        {
            filehandle  => $filehandle,
            set_errexit => 1,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Unload env
    say   {$filehandle} qq{## Unload environment: $env_name};
    say   {$filehandle} join $SPACE, @env_method_unload_cmds;
    print {$filehandle} $NEWLINE;

    return 1;
}

sub build_perl_program_check_command {

## Function : Build perl oneliner for testing program installation
## Returns  : $perl_commands
## Arguments: $programs_ref                => Programs to check {REF}
##          : $program_test_command_href  => Hash with test commands {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $programs_ref;
    my $program_test_command_href;

    my $tmpl = {
        programs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$programs_ref,
            strict_type => 1,
        },
        program_test_command_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$program_test_command_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Constants
    Readonly my $OPEN_STRING  => q/q{/;
    Readonly my $CLOSE_STRING => q/}/;
    Readonly my $TIMEOUT      => 20;

    ## Array for storing test commands
    my @program_test_commands;

  PROGRAM:
    foreach my $program ( @{$programs_ref} ) {

        ## Skip programs that lacks a test_command
        next PROGRAM if ( not $program_test_command_href->{$program} );

        ## Add path test
        if ( $program_test_command_href->{$program}{path} ) {
            my $path_test =
              $OPEN_STRING . $program_test_command_href->{$program}{path} . $CLOSE_STRING;
            my $path_test_name =
              $OPEN_STRING . q{Program in path: } . $program . $CLOSE_STRING;
            my $path_test_command =
              qq{ok(can_run( $path_test ), $path_test_name)} . $SEMICOLON;
            push @program_test_commands, $path_test_command;
        }
        ## Add execution test
        if ( $program_test_command_href->{$program}{execution} ) {
            my $execution_test =
                $OPEN_STRING
              . $program_test_command_href->{$program}{execution}
              . $CLOSE_STRING;
            my $execution_test_name =
              $OPEN_STRING . q{Can execute: } . $program . $CLOSE_STRING;
            my $execution_test_command =
qq{ok(run(command => $execution_test, timeout => $TIMEOUT), $execution_test_name)}
              . $SEMICOLON;
            push @program_test_commands, $execution_test_command;
        }
    }

    ## Return nothing if no tests are available for the programs
    return if ( not scalar @program_test_commands );

    ## Start perl oneliner
    my @perl_commands = q{perl -e};

    ## Add opening quote
    push @perl_commands, $SINGLE_QUOTE;

    ## Add required IPC module
    push @perl_commands, q{use IPC::Cmd qw{ can_run run };};

    ## Add required IPC module
    push @perl_commands, q{use Test::More;};

    ## Add program tests
    @perl_commands = ( @perl_commands, @program_test_commands );

    ## Done testing
    push @perl_commands, q{done_testing;};

    ## Add closing single quote
    push @perl_commands, $SINGLE_QUOTE;

    return @perl_commands;
}

1;
