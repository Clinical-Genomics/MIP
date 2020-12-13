#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Environment::Cluster} => [qw{ check_max_core_number }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Environment::Cluster qw{ check_max_core_number };

diag(   q{Test check_max_core_number from Cluster.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Set constants
Readonly my $BIG_REQUEST  => 10;
Readonly my $MAX_CORES    => 8;
Readonly my $OK_REQUEST   => 6;
Readonly my $RECIPE_CORES => 4;

my $core_number;

## Given a request within node limit
$core_number = check_max_core_number(
    {
        core_number_requested => $OK_REQUEST,
        max_cores_per_node    => $MAX_CORES,
    }
);
## Then return requested cores
is( $core_number, $OK_REQUEST, q{Ok request} );

## Given a request that exceeds the cores on the nodes
$core_number = check_max_core_number(
    {
        core_number_requested => $BIG_REQUEST,
        max_cores_per_node    => $MAX_CORES,
    }
);

## Then return max cores on node
is( $core_number, $MAX_CORES, q{Node constrained request} );

## Given a request that exceeds the number of cores allocated
$core_number = check_max_core_number(
    {
        core_number_requested => $OK_REQUEST,
        max_cores_per_node    => $MAX_CORES,
        recipe_core_number    => $RECIPE_CORES,
    }
);

## Then return max cores allocated by recipe
is( $core_number, $RECIPE_CORES, q{Recipe constrained request} );

done_testing();
