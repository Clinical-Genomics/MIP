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
use Test::Trap;

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
        q{MIP::Cluster}        => [qw{ get_memory_constrained_core_number }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Cluster qw{ get_memory_constrained_core_number };

diag(   q{Test get_memory_constrained_core_number from Cluster.pm v}
      . $MIP::Cluster::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Set constants
Readonly my $EXCESSIVE_MEMORY_REQUEST => 100;
Readonly my $BIG_MEMORY_REQUEST       => 20;
Readonly my $MAX_CORES                => 16;
Readonly my $OK_MEMORY_REQUEST        => 4;
Readonly my $RECIPE_CORES             => 16;
Readonly my $NODE_RAM_MEMORY          => 80;
Readonly my $CONSTRAINED_CORES        => 4;

my $core_number;

## Given an OK request
$core_number = get_memory_constrained_core_number(
    {
        memory_allocation  => $OK_MEMORY_REQUEST,
        node_ram_memory    => $NODE_RAM_MEMORY,
        max_cores_per_node => $MAX_CORES,
        recipe_core_number => $RECIPE_CORES,
    }
);

## Then return number of cores allocated by recipe
is( $core_number, $RECIPE_CORES, q{OK request} );

## Given a memory allocation greater than the core ram memory
$core_number = get_memory_constrained_core_number(
    {
        memory_allocation  => $BIG_MEMORY_REQUEST,
        node_ram_memory    => $NODE_RAM_MEMORY,
        max_cores_per_node => $MAX_CORES,
        recipe_core_number => $RECIPE_CORES,
    }
);

## Then limit the number of cores
is( $core_number, $CONSTRAINED_CORES, q{Limit program cores} );

## Given a memory allocation greater than the node ram memory
trap {
    $core_number = get_memory_constrained_core_number(
        {
            memory_allocation  => $EXCESSIVE_MEMORY_REQUEST,
            node_ram_memory    => $NODE_RAM_MEMORY,
            max_cores_per_node => $MAX_CORES,
            recipe_core_number => $RECIPE_CORES,
        }
      )
};
## Then croak
is( $trap->leaveby, q{die}, q{Croak on exessive memory request} );

done_testing();
