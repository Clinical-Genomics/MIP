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
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.01;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Dependency_tree} => [qw{ get_dependency_tree_order }],
        q{MIP::Test::Fixtures}  => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Dependency_tree qw{ get_dependency_tree_order };

diag(   q{Test get_dependency_tree_order from Dependency_tree.pm v}
      . $MIP::Dependency_tree::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my %dependency_tree = (
    ALL_CHAINS => [
        q{recipe_0},
        q{recipe_1},
        {
            CHAIN_0 => {
                PARALLEL => [
                    q{parallel_recipe_0},
                    q{parallel_recipe_1},
                    {
                        parallel_recipe_2 =>
                          [ q{parallel_recipe_2}, q{parallel_recipe_3}, ],
                    },
                ],
            },
        },
        q{recipe_2},
        {
            CHAIN_1 => {
                PARALLEL => [
                    q{parallel_recipe_4},
                    q{parallel_recipe_5},
                    {
                        parallel_recipe_6 =>
                          [ q{parallel_recipe_6}, q{parallel_recipe_7}, ],
                    },
                ],
            },
        },
        q{recipe_3},
    ],
);

## Expected results from traversing the dependency tree
my @expected_recipes =
  qw{ recipe_0 recipe_1 parallel_recipe_0 parallel_recipe_1 parallel_recipe_2 parallel_recipe_3 recipe_2 parallel_recipe_4 parallel_recipe_5 parallel_recipe_6 parallel_recipe_7 recipe_3 };

my @recipes;

get_dependency_tree_order(
    {
        dependency_tree_href => \%dependency_tree,
        recipes_ref          => \@recipes,
    }
);

is_deeply( \@recipes, \@expected_recipes, q{Got order of recipes } );

done_testing();
