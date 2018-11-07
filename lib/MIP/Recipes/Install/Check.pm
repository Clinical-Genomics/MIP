package MIP::Recipes::Install::Check;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use List::MoreUtils qw{ natatime };
use Readonly;

## MIPs lib/
use MIP::Unix::Write_to_file qw{ unix_write_to_file };
use MIP::Get::Parameter qw{ get_env_method_cmds };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_program_installations build_perl_program_check_command };
}

## Constants
Readonly my $DOUBLE_QUOTE => q{"};
Readonly my $NEWLINE      => qq{\n};
Readonly my $SINGLE_QUOTE => q{'};
Readonly my $SPACE        => q{ };

sub check_program_installations {

## Function : Write installation check oneliner to open filehandle
## Returns  :
## Arguments: $env                        => Program environment
##          : $FILEHANDLE                 => open filehandle
##          : $log                        => Log
##  		: $programs_ref                => Programs to check {REF}
##          : $program_test_command_href  => Hash with test commands {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $env;
    my $FILEHANDLE;
    my $log;
    my $programs_ref;
    my $program_test_command_href;

    my $tmpl = {
        env => {
            defined     => 1,
            required    => 1,
            store       => \$env,
            strict_type => 1,
        },
        FILEHANDLE => {
            required => 1,
            store    => \$FILEHANDLE,
        },
        log => {
            required => 1,
            store    => \$log,
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

    $log->info(qq{Writing tests for programs installed in environment: $env});

    ## Sort program list for readability
    my @sorted_programs = sort { lc $a cmp lc $b } @{$programs_ref};

    ## Getting environment commands
    my @env_method_load_cmds = get_env_method_cmds(
        {
            action     => q{load},
            env_method => q{conda},
            env_name   => $env,
        }
    );

    my @env_method_unload_cmds = get_env_method_cmds(
        {
            action     => q{unload},
            env_method => q{conda},
            env_name   => $env,
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
        $log->info(qq{No tests available for programs in environment: $env});
        return;
    }

    say {$FILEHANDLE} qq{## Testing programs installed in $env};

    ## Make bash output
    say {$FILEHANDLE}
q{echo -e '\n########################################################################\n'};
    say {$FILEHANDLE} q{echo -e "\tMIP has attempted to install the following programs"};
    say {$FILEHANDLE} q{echo -e "\tin in environment: } . $env . q{\n"};

    ## Print programs selected for installation in set of 5.
    my $program_iterator = natatime $COLUMNS, @sorted_programs;

    ## Create field layout for printf
    my $field = q{"\t} . ( ( q{%-} . $COLUMN_WIDTH . q{s} ) x $COLUMNS ) . q{\n"};

    while ( my @print_set = $program_iterator->() ) {
        ## Print programs according to field
        say {$FILEHANDLE} qq{printf $field } . q{"}
          . join( $DOUBLE_QUOTE . $SPACE . $DOUBLE_QUOTE, @print_set ) . q{"};
    }

    say {$FILEHANDLE} q{echo -e "\n\tTesting installation\n"};

    ## Load env
    say   {$FILEHANDLE} qq{## Load environment: $env};
    say   {$FILEHANDLE} join $SPACE, @env_method_load_cmds;
    print {$FILEHANDLE} $NEWLINE;

    ## Write test oneliner
    say   {$FILEHANDLE} q{## Test programs};
    say   {$FILEHANDLE} join $SPACE, @perl_commands;
    print {$FILEHANDLE} $NEWLINE;

    ## Unload env
    say   {$FILEHANDLE} qq{## Unload environment: $env};
    say   {$FILEHANDLE} join $SPACE, @env_method_unload_cmds;
    print {$FILEHANDLE} $NEWLINE;

    return;
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
    Readonly my $SEMICOLON    => q{;};
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
    return if scalar @program_test_commands == 0;

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
