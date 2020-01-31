package MIP::Recipes::Parse;

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
    our @EXPORT_OK = qw{ parse_recipes };
}

sub parse_recipes {

## Function : Set default analysis type to active parameters
## Returns  :
## Arguments: $active_parameter_href  => Holds all set parameter for analysis {REF}
##          : $parameter_href         => Hash with parameters from yaml file {REF}
##          : parameter_to_check_href => Map options to check depending on type {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_href;
    my $parameter_to_check_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        parameter_to_check_href => {
            default     => {},
            required    => 1,
            store       => \$parameter_to_check_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Active_parameter qw{ set_recipe_resource };
    use MIP::Parameter qw{ parse_parameter_recipe_names };
    use MIP::Recipes::Check qw{ check_recipe_exists_in_hash };
    ## Unpack
    my @parameter_keys_to_check     = @{ $parameter_to_check_href->{keys} };
    my @parameter_elements_to_check = @{ $parameter_to_check_href->{elements} };

  PARAMETER_NAME:
    foreach my $parameter_name (@parameter_keys_to_check) {

        ## Test if key from query hash exists truth hash
        check_recipe_exists_in_hash(
            {
                parameter_name => $parameter_name,
                query_ref      => \%{ $active_parameter_href->{$parameter_name} },
                truth_href     => $parameter_href,
            }
        );
    }

    ## Set recipe resource allocation for specific recipe(s)
    set_recipe_resource( { active_parameter_href => $active_parameter_href, } );

    parse_parameter_recipe_names( { parameter_href => $parameter_href, } );

    ## Parameters that have elements as MIP recipe names
    foreach my $parameter_name (@parameter_elements_to_check) {

        ## Test if element from query array exists truth hash
        check_recipe_exists_in_hash(
            {
                parameter_name => $parameter_name,
                query_ref      => \@{ $active_parameter_href->{$parameter_name} },
                truth_href     => $parameter_href,
            }
        );
    }
    return 1;
}

1;
