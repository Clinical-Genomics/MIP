package MIP::Environment::Cluster;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use List::Util qw{ min };
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
    our @EXPORT_OK = qw{ check_max_core_number };
}

sub check_max_core_number {

## Function : Limit number of cores requested to the maximum number of cores available
## Returns  : $core_number
## Arguments: $core_number_requested => Number of cores requested to allocate
##          : $max_cores_per_node    => Max number of cores per node
##          : $recipe_core_number    => Number of cores allocated by the recipe

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $core_number_requested;
    my $max_cores_per_node;
    my $recipe_core_number;

    my $tmpl = {
        core_number_requested => {
            required    => 1,
            defined     => 1,
            store       => \$core_number_requested,
            strict_type => 1,
        },
        max_cores_per_node => {
            defined     => 1,
            required    => 1,
            store       => \$max_cores_per_node,
            strict_type => 1,
        },
        recipe_core_number => {
            defined     => 1,
            store       => \$recipe_core_number,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Set default core constraint
    my $core_constraint = $max_cores_per_node;

    ## Update $core_constraint
    if ($recipe_core_number) {

        $core_constraint = min( $max_cores_per_node, $recipe_core_number );
    }

    ## Check number of cores available
    if ( $core_number_requested > $core_constraint ) {

        ## Return maximum available cores
        return $core_constraint;
    }

        ## Core number requested is lower than the number of available cores
        return $core_number_requested;
}

1;
