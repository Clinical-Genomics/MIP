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
use MIP::Constants qw{ $COMMA $LOG $SPACE };
use MIP::Log::MIP_log4perl qw{ initiate_logger };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

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
        q{MIP::Cluster}        => [qw{ check_recipe_memory_allocation }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Cluster qw{ check_recipe_memory_allocation };

diag(   q{Test check_recipe_memory_allocation from Cluster.pm v}
      . $MIP::Cluster::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

test_log( {} );

## Given a memory allocation that stays within the node ram memory
Readonly my $NODE_RAM_MEMORY   => 180;
Readonly my $MEMORY_ALLOCATION => 30;

my $memory_allocation = check_recipe_memory_allocation(
    {
        node_ram_memory          => $NODE_RAM_MEMORY,
        recipe_memory_allocation => $MEMORY_ALLOCATION,
    }
);

## Then return the memory allocation
is( $memory_allocation, $MEMORY_ALLOCATION, q{Got memory allocation} );

## Given a memory allocation that exceeds the available memory
Readonly my $MEMORY_ALLOCATION_2 => 200;

trap {
    $memory_allocation = check_recipe_memory_allocation(
        {
            node_ram_memory          => $NODE_RAM_MEMORY,
            recipe_memory_allocation => $MEMORY_ALLOCATION_2,
        }
    )
};

## Then limit to available memory and warn
is( $memory_allocation, $NODE_RAM_MEMORY, q{Got max memory} );
like( $trap->stderr, qr{WARN}xms, q{Throw warning} );

done_testing();
