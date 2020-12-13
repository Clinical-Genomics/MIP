package MIP::Environment::Child_process;

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
use MIP::Constants qw{ $SPACE $TEST_MODE %TEST_PROCESS_RETURN };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ child_process };
}

sub child_process {

## Function : Wrapper to open a child process. Return the output from the child process in hash
## Returns  : %child_process
## Arguments: $commands_ref => Command to process {REF}
##          : $process_type => Type of child process to use
##          : $verbose      => Verbosity level

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $commands_ref;
    my $process_type;
    my $verbose;

    my $tmpl = {
        commands_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$commands_ref,
            strict_type => 1,
        },
        process_type => {
            allow       => [qw{ ipc_cmd_run open3 }],
            defined     => 1,
            required    => 1,
            store       => \$process_type,
            strict_type => 1,
        },
        verbose => {
            default     => 0,
            store       => \$verbose,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::System_call qw{ ipc_cmd_run ipc_open3 };

    ## Return mock results if in test mode
    return %TEST_PROCESS_RETURN if ($TEST_MODE);

    my %process_api = (
        open3 => {
            method   => \&ipc_open3,
            arg_href => { command_string => join $SPACE, @{$commands_ref}, },
        },
        ipc_cmd_run => {
            method   => \&ipc_cmd_run,
            arg_href => {
                commands_ref => $commands_ref,
                verbose      => $verbose,
            },
        },
    );

    my %process_return = $process_api{$process_type}{method}
      ->( { %{ $process_api{$process_type}{arg_href} } } );

    return %process_return;
}

1;
