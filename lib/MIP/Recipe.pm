package MIP::Recipe;

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
    our @EXPORT_OK = qw{ parse_recipe_prerequisites };
}

sub parse_recipe_prerequisites {

## Function : Parse recipe prerequisites and return return them
## Returns  : %recipe
## Arguments: $active_parameter_href => The active parameters for this analysis hash {REF}
##          : $parameter_href        => Holds all parameters
##          : $recipe_name           => Recipe name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_href;
    my $recipe_name;

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
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            strict_type => 1,
            store       => \$recipe_name,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Active_parameter qw{ get_recipe_resources };
    use MIP::Parameter qw{ get_recipe_attributes };

    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );
    $recipe_resource{job_id_chain} = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            attribute      => q{chain},
        }
    );

    return %recipe_resource;
}

1;
