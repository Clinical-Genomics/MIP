package MIP::List;

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

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ check_allowed_array_values check_element_exist_hash_of_array get_splitted_lists };
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

    return 0 if ( not exists ${$hash_ref}{$key} );

    return 1 if ( any { $_ eq $element } @{ $hash_ref->{$key} } );

    return;
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
