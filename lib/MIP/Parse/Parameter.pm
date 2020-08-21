package MIP::Parse::Parameter;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile splitpath };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.19;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      parse_download_reference_parameter
    };

}

sub parse_download_reference_parameter {

## Function : Remodel depending on if "--reference" was used or not as the user info is stored as a scalar per reference_id while yaml is stored as arrays per reference_id
## Returns  :
## Arguments: $reference_href => Reference hash {REF}

    my ($arg_href) = @_;

## Flatten argument(s)
    my $reference_href;

    my $tmpl = {
        reference_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$reference_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  VERSION_REF:
    foreach my $versions_ref ( values %{$reference_href} ) {

        if ( ref $versions_ref ne q{ARRAY} ) {

            ## Make scalar from CLI '--ref key=value' option into array
            $versions_ref = [$versions_ref];
        }
    }

    return;
}

1;
