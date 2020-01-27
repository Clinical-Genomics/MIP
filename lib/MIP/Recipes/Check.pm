package MIP::Recipes::Check;

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
use MIP::Constants qw{ $LOG_NAME $SINGLE_QUOTE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_recipe_exists_in_hash };
}

sub check_recipe_exists_in_hash {

## Function : Test if parameter "recipe name" from query parameter exists in truth hash
## Returns  :
## Arguments: $parameter_name => Parameter name
##          : $query_ref      => Query (ARRAY|HASH|SCALAR) {REF}
##          : $truth_href     => Truth hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_name;
    my $query_ref;
    my $truth_href;

    my $tmpl = {
        parameter_name => { defined => 1, required => 1, store => \$parameter_name, },
        query_ref => {
            defined  => 1,
            required => 1,
            store    => \$query_ref,
        },
        truth_href     => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$truth_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my $error_msg = qq{$SINGLE_QUOTE - Does not exist as recipe parameter in MIP};

    if ( ref $query_ref eq q{HASH} ) {

      RECIPE_NAME:
        foreach my $recipe_name ( keys %{$query_ref} ) {

            next RECIPE_NAME if ( exists $truth_href->{$recipe_name} );

            $log->fatal(
                $parameter_name . qq{ key $SINGLE_QUOTE} . $recipe_name . $error_msg );
            exit 1;
        }
    }
    if ( ref $query_ref eq q{ARRAY} ) {

      RECIPE_NAME:
        foreach my $recipe_name ( @{$query_ref} ) {

            next RECIPE_NAME if ( exists $truth_href->{$recipe_name} );

            $log->fatal( $parameter_name
                  . qq{ element $SINGLE_QUOTE}
                  . $recipe_name
                  . $error_msg );
            exit 1;
        }
    }
    if ( ref $query_ref eq q{SCALAR} ) {

        return if ( exists $truth_href->{$parameter_name} );

        $log->fatal(
            $parameter_name . qq{ element $SINGLE_QUOTE} . $parameter_name . $error_msg );
        exit 1;
    }
    return;
}

1;
