package MIP::Check::Unix;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use Params::Check qw{ check allow last_error };
$Params::Check::PRESERVE_CASE = 1;    #Do not convert to lower case

use IPC::Cmd qw{ can_run };

## MIPs lib/

## Constants

BEGIN {

    use base qw{Exporter};
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{check_binary_in_path};
}

sub check_binary_in_path {

## binary_in_path

## Function   : Scans through PATH for supplied binary
## Returns    :
## Arguments  : $binary
##            : $binary => $binary to search for

    my ($arg_href) = @_;

    ## Flatten argument(s)

    my $binary;

    my $tmpl = {
        binary => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$binary,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Search for conda in PATH and exit if not present
    if ( can_run($binary) ) {
        say STDERR q{Program check: } . $binary . q{ in PATH};
    }
    else {
        say STDERR q{Could not detect } . $binary . q{in your PATH};
        exit 1;
    }

    return;
}

1;
