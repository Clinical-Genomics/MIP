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
our $VERSION = 1.00;

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
        q{MIP::Dependency_tree} => [qw{ get_recipes_for_dependency_tree_chain }],
        q{MIP::Test::Fixtures}  => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Dependency_tree qw{ get_recipes_for_dependency_tree_chain };
use MIP::Test::Fixtures qw{ test_mip_hashes };

diag(   q{Test get_recipes_for_dependency_tree_chain from Dependency_tree.pm v}
      . $MIP::Dependency_tree::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my %dependency_tree = test_mip_hashes(
    {
        mip_hash_name => q{dependency_tree_dna},
    }
);
my $chain_initiation_point  = q{DELLY_CALL};
my $dependency_subtree_href = {};
my $recipe_initiation_point = q{delly_reformat};

## Given request to get recipes in chain DELLY_CALL starting with delly_reformat
my @recipes = get_recipes_for_dependency_tree_chain(
    {
        chain_initiation_point  => $chain_initiation_point,
        dependency_tree_href    => \%dependency_tree,
        recipe_initiation_point => $recipe_initiation_point,
    }
);

## Then get it
my @expected_recipes = qw{ delly_reformat };
is_deeply( \@recipes, \@expected_recipes, q{Get recipes} );

done_testing();
