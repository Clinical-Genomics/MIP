package MIP::Check::Hash;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use List::MoreUtils qw { any };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_element_exist_hash_of_array };
}

## Constants
Readonly my $SPACE => q{ };

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
        element =>
          { required => 1, defined => 1, strict_type => 1, store => \$element },
        hash_ref => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$hash_ref
        },
        key =>
          { required => 1, defined => 1, strict_type => 1, store => \$key },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Information on entry present
    if ( exists ${$hash_ref}{$key} ) {

        # If element is not part of array
        if ( not any { $_ eq $element } @{ $hash_ref->{$key} } ) {
            return 1;
        }
    }
    return;
}

1;
