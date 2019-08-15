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
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Cluster}        => [qw{ get_core_number }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Cluster qw{ get_core_number };

diag(   q{Test get_core_number from Cluster.pm v}
      . $MIP::Cluster::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given recipe core number and modifier core number
Readonly my $RECIPE_CORE_NUMBER   => 1;
Readonly my $MODIFIER_CORE_NUMBER => 2;
Readonly my $MAX_CORES_PER_NODE   => 2;

my $returned_recipe_core_number = get_core_number(
    {
        recipe_core_number   => $RECIPE_CORE_NUMBER,
        modifier_core_number => $MODIFIER_CORE_NUMBER,
        max_cores_per_node   => $MAX_CORES_PER_NODE,
    }
);

## Then return recipe core number and ignore modifier
is( $returned_recipe_core_number, 1, q{Got module core number} );

## Given modifier core number, and no recipe core number
my $returned_modified_core_number = get_core_number(
    {
        recipe_core_number   => undef,
        modifier_core_number => $MODIFIER_CORE_NUMBER,
        max_cores_per_node   => $MAX_CORES_PER_NODE,
    }
);

## Then return modifier core number
is( $returned_modified_core_number, 2, q{Got modifier core number} );

done_testing();
