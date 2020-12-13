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
        q{MIP::Dependency_tree} => [qw{ get_recipe_dependency_tree_chain }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Dependency_tree qw{ get_recipe_dependency_tree_chain };

diag(   q{Test get_recipe_dependency_tree_chain from Dependency_tree.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a dependency tree and a recipe name when not on a PARALLEL chain
my $recipe_name     = q{program_2};
my %dependency_tree = (
    CHAIN_ALL => [
        q{program_0},
        q{program_1},
        {
            CHAIN_0 => {
                PARALLEL => [
                    q{parallel_program_0},
                    q{parallel_program_1},
                    {
                        PARALLEL_PROGRAM_2 =>
                          [ q{parallel_program_2}, q{parallel_program_3}, ],
                    },
                ],
            },
        },
        q{program_2},
        {
            CHAIN_1 => {
                PARALLEL => [
                    q{parallel_program_4},
                    q{parallel_program_5},
                    {
                        CHAIN_MAIN => [ q{parallel_program_6}, q{parallel_program_7}, ],
                    },
                ],
            },
        },
        q{program_3},
    ],
);

my $recipe_chain;
get_recipe_dependency_tree_chain(
    {
        recipe               => $recipe_name,
        dependency_tree_href => \%dependency_tree,
        chain_id_ref         => \$recipe_chain,
    }
);
my $expected_ordinary_chain = q{CHAIN_ALL};

## Then recipe chain should be set
is( $recipe_chain, $expected_ordinary_chain, q{Got chain which recipe belongs to} );

## Given a recipe name within a parallel chain, but in a chain part
$recipe_name  = q{parallel_program_7};
$recipe_chain = undef;

get_recipe_dependency_tree_chain(
    {
        recipe               => $recipe_name,
        dependency_tree_href => \%dependency_tree,
        chain_id_ref         => \$recipe_chain,
    }
);
my $expected_within_parallel_block = q{CHAIN_MAIN};

## Then recipe chain should be set
is(
    $recipe_chain,
    $expected_within_parallel_block,
    q{Got chain which recipe belongs to within parallel block}
);

## Given a recipe name in a parallel chain
$recipe_name  = q{parallel_program_1};
$recipe_chain = undef;

get_recipe_dependency_tree_chain(
    {
        recipe               => $recipe_name,
        dependency_tree_href => \%dependency_tree,
        chain_id_ref         => \$recipe_chain,
    }
);
my $expected_in_parallel_block = q{CHAIN_0};

## Then recipe chain should upstream of parallel block should be set
is( $recipe_chain, $expected_in_parallel_block,
    q{Got upstream chain for recipe in parallel block} );

done_testing();
