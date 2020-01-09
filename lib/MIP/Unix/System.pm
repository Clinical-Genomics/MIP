package MIP::Unix::System;

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
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ system_cmd_call };
}

sub system_cmd_call {

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
    waitpid $pid, 0;

    # Capture output
    if ($reader) {

        @{ $chld_handlers{output} } = <$reader>;
    }

    # Capture error
    if ($error) {

        @{ $chld_handlers{error} } = <$error>;
    }

    return %chld_handlers;
}

1;
