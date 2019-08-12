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
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

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
        q{MIP::Cluster}        => [qw{ update_memory_allocation }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Cluster qw{ update_memory_allocation };

diag(   q{Test update_memory_allocation from Cluster.pm v}
      . $MIP::Cluster::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given the number of parallel processes and the amount of memory to reserve for each process
Readonly my $NODE_RAM_MEMORY           => 180;
Readonly my $PROCESS_MEMORY_ALLOCATION => 3;
Readonly my $PARALLEL_PROCESSES        => 5;

test_log( { no_screen => 1, } );

my $memory_allocation = update_memory_allocation(
    {
        node_ram_memory           => $NODE_RAM_MEMORY,
        parallel_processes        => $PARALLEL_PROCESSES,
        process_memory_allocation => $PROCESS_MEMORY_ALLOCATION,
    }
);

## Then return the memory allocation
Readonly my $MEMORY_ALLOCATION => 15;
is( $memory_allocation, $MEMORY_ALLOCATION, q{Got memory allocation} );

done_testing();
