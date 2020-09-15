package MIP::List;

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
    our @EXPORT_OK = qw{ check_allowed_array_values get_splitted_lists };
}

sub check_allowed_array_values {

## Function : Check that the array values are allowed
## Returns  :
## Arguments: $allowed_values_ref => Allowed values for parameter {REF}
##          : $values_ref         => Values for parameter {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $allowed_values_ref;
    my $values_ref;

    my $tmpl = {
        allowed_values_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$allowed_values_ref,
            strict_type => 1,
        },
        values_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$values_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Array::Utils qw{ array_minus };

    # Check if arrays contains allowed values
    return 0 if ( array_minus( @{$values_ref}, @{$allowed_values_ref} ) );

    # All ok
    return 1;
}

sub get_splitted_lists {

## Function : Split list on regexp
## Returns  : \@match_list, \@no_match_list
## Arguments: $list_ref => List to split
##          : $regexp   => Regexp to split list on

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $list_ref;
    my $regexp;

    my $tmpl = {
        list_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$list_ref,
            strict_type => 1,
        },
        regexp => {
            default     => qr//,
            defined     => 1,
            required    => 1,
            store       => \$regexp,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @match_list    = grep { m/$regexp/xms } @{$list_ref};
    my @no_match_list = grep { not m/$regexp/xms } @{$list_ref};

    return \@match_list, \@no_match_list;
}

1;
