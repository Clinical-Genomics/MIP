package MIP::Script::Utils;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Basename qw{ basename };
use File::Spec::Functions qw{ catdir };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings qw{ FATAL utf8 };
use warnings;

## CPANM
use Readonly;
use autodie;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.10;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ help };
}

sub help {

## Function : Print help text and exit with supplied exit code
## Returns  :
## Arguments: $exit_code => Exit code
##          : $usage     => Help text

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $usage;

    ## Default(s)
    my $exit_code;

    my $tmpl = {
        exit_code => {
            allow       => qr/ ^\d+$ /xsm,
            default     => 0,
            store       => \$exit_code,
            strict_type => 1,
        },
        usage => {
            defined     => 1,
            required    => 1,
            store       => \$usage,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    say {*STDOUT} $usage;
    exit $exit_code;
}

1;
