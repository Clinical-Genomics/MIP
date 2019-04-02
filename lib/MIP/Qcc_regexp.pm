package MIP::Qcc_regexp;

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
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ get_qcc_regexp_recipe_attribute };
}

sub get_qcc_regexp_recipe_attribute {

## Function : Get recipe attributes from qccollect reg exp hash
## Returns  : "$attribute" or "$attribute_href"
## Arguments: $attribute       => Attribute key
##          : $recipe_name     => Recipe to get attributes from
##          : $qcc_regexp_href => qccollect regexp hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $attribute;
    my $recipe_name;
    my $qcc_regexp_href;

    my $tmpl = {
        attribute => {
            store       => \$attribute,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        qcc_regexp_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qcc_regexp_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Get and return attribute value
    if ( defined $attribute && $attribute ) {

        return $qcc_regexp_href->{$recipe_name}{$attribute};
    }

    ## Get recipe attribute hash
    return %{ $qcc_regexp_href->{$recipe_name} };
}

1;
