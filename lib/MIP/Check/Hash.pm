package MIP::Check::Hash;

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
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_element_exist_hash_of_array };
}

sub check_element_exist_hash_of_array {

## Function : Test element for being part of hash of array at supplied key.
## Returns  : Return "1" if element is not part of array
## Arguments: $element  => Element to look for in hash of array
##          : $hash_ref => Hash {REF}
##          : $key      => The key pointing to the array in the $hash_ref

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $element;
    my $hash_ref;
    my $key;

    my $tmpl = {
        element => {
            defined     => 1,
            required    => 1,
            store       => \$element,
            strict_type => 1,
        },
        hash_ref => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$hash_ref,
            strict_type => 1,
        },
        key => { required => 1, defined => 1, strict_type => 1, store => \$key },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use List::Util qw{ any };

    # Information on entry present
    if ( exists ${$hash_ref}{$key} ) {

        # If element is part of array
        if ( any { $_ eq $element } @{ $hash_ref->{$key} } ) {
            return 1;
        }
    }
    return;
}

1;
