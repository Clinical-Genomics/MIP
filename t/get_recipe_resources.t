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
use Modern::Perl qw{ 2014 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COLON $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_mip_hashes test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.02;

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
        q{MIP::Get::Parameter} => [qw{ get_recipe_resources }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Get::Parameter qw{ get_recipe_resources };

diag(   q{Test get_recipe_resources from Parameter.pm v}
      . $MIP::Get::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

test_log( {} );

## Given a recipe name and active parameter hash
my %active_parameter = test_mip_hashes( { mip_hash_name => q{active_parameter}, } );
my $recipe_name      = q{bwa_mem};

my %recipe_resource = get_recipe_resources(
    {
        active_parameter_href => \%active_parameter,
        recipe_name           => $recipe_name,
    }
);

## Then return recipe resource hash
my %expected = (
    core_number  => 30,
    memory       => 120,
    time         => 30,
    load_env_ref => [qw{conda activate test }],
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
$active_parameter{recipe_memory}{bwa_mem} = undef;

my $recipe_memory = get_recipe_resources(
    {
        active_parameter_href => \%active_parameter,
        recipe_name           => $recipe_name,
        resource              => q{memory},
    }
);

## Then return 5 gigs times the number of cores
Readonly my $DEFAULT_MEMORY => 150;
is( $recipe_memory, $DEFAULT_MEMORY, q{Got default memory} );

## Given a recipe that lacks memory and core specification
$active_parameter{recipe_core_number}{bwa_mem} = undef;

$recipe_memory = get_recipe_resources(
    {
        active_parameter_href => \%active_parameter,
        recipe_name           => $recipe_name,
        resource              => q{memory},
    }
);

## Then return the core ram memory
Readonly my $CORE_MEMORY => 5;
is( $recipe_memory, $CORE_MEMORY, q{Got core memory} );

## Given a recipe that lacks core specification but has a memory specification
$active_parameter{recipe_memory}{bwa_mem} = 1;

$recipe_memory = get_recipe_resources(
    {
        active_parameter_href => \%active_parameter,
        recipe_name           => $recipe_name,
        resource              => q{memory},
    }
);

## Then return the recipe ram memory
Readonly my $CORE_MEMORY_1 => 1;
is( $recipe_memory, $CORE_MEMORY_1, q{Got core memory} );
done_testing();
