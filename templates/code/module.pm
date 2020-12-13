package MIP::PATH::TO::MODULE;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ space separated subroutines };
}

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
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$arrays_ref,
            strict_type => 1,
        },
        hash_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$hash_href,
            strict_type => 1,
        },
        scalar => {
            allow       => qr/ \A \d+ \z /sxm,
            default     => 1,
            store       => \$scalar,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return;
}

1;
