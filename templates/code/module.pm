package MIP::PATH::TO::MODULE;

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

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ space separated subroutines };
}

## Constants
Readonly my $SPACE => q{ };

sub name_of_subroutine {

## Function :
## Returns  :
## Arguments: $arrays_ref => Array ref description {REF}
##          : $hash_href  => Hash ref description {REF}
##          : $scalar     => Scalar description

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $hash_href;
    my $arrays_ref;

    ## Default(s)
    my $scalar;

    my $tmpl = {
        arrays_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$arrays_ref,
        },
        hash_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$hash_href,
        },
        scalar => {
            default     => 1,
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$scalar,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return;
}

1;
