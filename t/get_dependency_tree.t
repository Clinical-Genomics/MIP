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
        q{MIP::Dependency_tree} => [qw{ get_dependency_tree }],
        q{MIP::Test::Fixtures}  => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Dependency_tree qw{ get_dependency_tree };

diag(   q{Test get_dependency_tree from Dependency_tree.pm v}
      . $MIP::Dependency_tree::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given an initiation map
my %dependency_tree = (
    CHAIN_ALL => [
        q{recipe_12},
        {
            CHAIN_MAIN => [ qw{ recipe_0 }, ],
        },
        {
            CHAIN_0 => [ qw{ recipe_1 }, ],
        },
        q{recipe_2},
        {
            CHAIN_MAIN => [
                qw{ recipe_3 recipe_4 },
                {
                    CHAIN_1 => [qw{ recipe_5 }],
                },
                {
                    CHAIN_2 => [
                        {
                            PARALLEL => [
                                qw{ parallel_recipe_0 parallel_recipe_1 },
                                {
                                    PARALLEL_RECIPE_2 =>
                                      [ qw{ parallel_recipe_2 parallel_recipe_3 }, ],
                                },
                            ],
                        },
                        qw{ recipe_6 recipe_7 },
                    ],
                },
            ],
        },
        {
            CHAIN_MAIN => [
                {
                    PARALLEL => [
                        qw{ parallel_recipe_4 parallel_recipe_5 },
                        {
                            PARALLEL_RECIPE_6 =>
                              [ qw{ parallel_recipe_6 parallel_recipe_7 }, ],
                        },
                    ],
                },
                qw{ recipe_8 recipe_9 },
                {
                    CHAIN_3 => [ qw{ recipe_10 }, ],
                },
            ],
        },
        q{recipe_11},
    ],
);

## Expected results from traversing the dependency tree with different starting points
my @expected_recipes_12 = qw{ recipe_12 recipe_2 recipe_11 };
my @expected_recipes_0 =
  qw{ recipe_0 recipe_1 recipe_2 recipe_3 recipe_4 recipe_5 parallel_recipe_0 parallel_recipe_1 parallel_recipe_2 parallel_recipe_3 recipe_6 recipe_7 parallel_recipe_4 parallel_recipe_5 parallel_recipe_6 parallel_recipe_7 recipe_8 recipe_9 recipe_10 recipe_11 };
my @expected_recipes_1          = qw{ recipe_1 recipe_2 recipe_11 };
my @expected_parallel_recipes_0 = qw{ parallel_recipe_0 recipe_6 recipe_7 recipe_11 };
my @expected_parallel_recipes_2 =
  qw{ parallel_recipe_2 parallel_recipe_3 recipe_6 recipe_7 recipe_11 };
my @expected_parallel_recipes_3 = qw{ parallel_recipe_3 recipe_6 recipe_7 recipe_11 };
my @expected_recipe_2           = qw{ recipe_2 recipe_11 };
my @expected_parallel_recipes_5 =
  qw{ parallel_recipe_5 recipe_8 recipe_9 recipe_10 recipe_11 };
my @expected_recipe_11 = qw{ recipe_11 };

## Define tests
my %test = (
    recipe_12         => \@expected_recipes_12,
    recipe_0          => \@expected_recipes_0,
    recipe_1          => \@expected_recipes_1,
    parallel_recipe_0 => \@expected_parallel_recipes_0,
    parallel_recipe_2 => \@expected_parallel_recipes_2,
    parallel_recipe_3 => \@expected_parallel_recipes_3,
    recipe_2          => \@expected_recipe_2,
    parallel_recipe_5 => \@expected_parallel_recipes_5,
    recipe_11         => \@expected_recipe_11,
);

## Run tests
while ( my ( $test_key, $expected_recipes_ref ) = each %test ) {

    my @start_with_recipes;
    my $is_recipe_found   = 0;
    my $is_chain_found    = 0;
    my $start_with_recipe = $test_key;

    get_dependency_tree(
        {
            dependency_tree_href   => \%dependency_tree,
            is_recipe_found_ref    => \$is_recipe_found,
            is_chain_found_ref     => \$is_chain_found,
            recipe                 => $start_with_recipe,
            start_with_recipes_ref => \@start_with_recipes,
        }
    );

## Then match the expected result
    is_deeply(
        \@start_with_recipes,
        \@{$expected_recipes_ref},
        q{Start with } . $test_key
    );
}

done_testing();
