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
use MIP::Constants qw{ $COLON $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_mip_hashes };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Active_parameter} => [qw{ get_recipe_resources }],
        q{MIP::Test::Fixtures}   => [qw{ test_log test_mip_hashes }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Active_parameter qw{ get_recipe_resources };

diag(   q{Test get_recipe_resources from Active_parameter.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $CORE_MEMORY        => 5;
Readonly my $DEFAULT_MEMORY     => 175;
Readonly my $RECIPE_CORE_MEMORY => 1;

test_log( { no_screen => 1, } );

## Given a recipe name and active parameter hash
my %active_parameter =
  test_mip_hashes( { mip_hash_name => q{active_parameter}, recipe_name => q{deepvariant}, } );
my $recipe_name = q{deepvariant};

my %recipe_resource = get_recipe_resources(
    {
        active_parameter_href => \%active_parameter,
        recipe_name           => $recipe_name,
    }
);

## Then return recipe resource hash
my %expected = (
    core_number  => 35,
    gpu_number   => 1,
    load_env_ref => [qw{ conda activate test }],
    memory       => 175,
    mode         => 2,
    time         => 10,
);
is_deeply( \%recipe_resource, \%expected, q{Got recipe resource hash} );

## Given a request for a specific resource
RESOURCE:
foreach my $resource ( keys %expected ) {

    my $recipe_resource = get_recipe_resources(
        {
            active_parameter_href => \%active_parameter,
            recipe_name           => $recipe_name,
            resource              => $resource,
        }
    );

    ## Then return recipe resource
    is_deeply( $recipe_resource, $expected{$resource}, q{Got} . $SPACE . $resource );
}

## Given a recipe that lacks memory specification
$active_parameter{recipe_memory}{deepvariant} = undef;

my $recipe_memory = get_recipe_resources(
    {
        active_parameter_href => \%active_parameter,
        recipe_name           => $recipe_name,
        resource              => q{memory},
    }
);

## Then return 5 gigs times the number of cores
is( $recipe_memory, $DEFAULT_MEMORY, q{Got default memory} );

## Given a recipe that lacks memory and core specification
$active_parameter{recipe_core_number}{deepvariant} = undef;

$recipe_memory = get_recipe_resources(
    {
        active_parameter_href => \%active_parameter,
        recipe_name           => $recipe_name,
        resource              => q{memory},
    }
);

## Then return the core ram memory
is( $recipe_memory, $CORE_MEMORY, q{Got core memory} );

## Given a recipe that lacks core specification but has a memory specification
$active_parameter{recipe_memory}{deepvariant} = 1;

$recipe_memory = get_recipe_resources(
    {
        active_parameter_href => \%active_parameter,
        recipe_name           => $recipe_name,
        resource              => q{memory},
    }
);

## Then return the recipe ram memory
is( $recipe_memory, $RECIPE_CORE_MEMORY, q{Got core memory} );

done_testing();
