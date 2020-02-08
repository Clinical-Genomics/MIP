package MIP::System_call;

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

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ ipc_cmd_run ipc_open3 };
}

sub ipc_cmd_run {

## Function : Wrapper for ipc cmd run sub
## Returns  :
## Arguments: $commands_ref => Commands to run {REF}
##          : $verbose       => Verbosity level

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $commands_ref;

    ## Defaults(s)
    my $verbose;

    my $tmpl = {
        commands_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$commands_ref,
            strict_type => 1,
        },
        verbose => {
            default     => 0,
            store       => \$verbose,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use IPC::Cmd qw{ run };

    # System call
    my ( $success, $error_message_ref, $full_buf_ref, $stdout_buf_ref, $stderr_buf_ref )
      = run( command => $commands_ref, verbose => $verbose );

    my %process_return = (
        error_messages_ref => $error_message_ref,
        buffers_ref        => $full_buf_ref,
        stdouts_ref        => $stdout_buf_ref,
        stderrs_ref        => $stderr_buf_ref,
        success            => $success,
    );
    return %process_return;
}

sub ipc_open3 {

## Function : Open a process for reading, writing, and error handling using open3(). Return the output from the child process in %chld_handlers{output|error}.
## Returns  : %chld_handlers
## Arguments: $command_string => Command to system

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $command_string;

    my $tmpl = {
        command_string => {
            defined     => 1,
            required    => 1,
            store       => \$command_string,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use IPC::Open3 qw{ open3 };
    use Symbol qw{ gensym };

    # Store handler of process
    my %chld_handlers;

    # Initilize
    my ( $writer, $reader, $error );

    ## Creates an anonymous glob and returns a reference to it

    ## Required to capture error
    $error = gensym();

    # System call
    my $pid = open3( $writer, $reader, $error, qq{$command_string} );

    # Terminate process
    waitpid $pid, 0 or croak(qq{Child process died: $OS_ERROR});

    # Capture output
    if ($reader) {

        @{ $chld_handlers{stdouts_ref} } = <$reader>;
    }

    # Capture error
    if ($error) {

        @{ $chld_handlers{stderrs_ref} } = <$error>;
    }

    return %chld_handlers;
}

1;
