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
        q{MIP::Cluster}        => [qw{ get_parallel_processes }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Cluster qw{ get_parallel_processes };

diag(   q{Test get_parallel_processes from Cluster.pm v}
      . $MIP::Cluster::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given input that puts a limitation on memory
Readonly my $CORE_NUMBER               => 5;
Readonly my $PROCESS_MEMORY_ALLOCATION => 5;
Readonly my $RECIPE_MEMORY_ALLOCATION  => 12;

my $parallel_processes = get_parallel_processes(
    {
        core_number               => $CORE_NUMBER,
        process_memory_allocation => $PROCESS_MEMORY_ALLOCATION,
        recipe_memory_allocation  => $RECIPE_MEMORY_ALLOCATION,
    }
);

## Then constrain the number of parallel processes
Readonly my $OK_PARALLEL_PROCESSES => 2;
is( $parallel_processes, $OK_PARALLEL_PROCESSES,
    q{Constrain parallel processes to available memory} );

## Given input that puts a limitation on cores
Readonly my $RECIPE_MEMORY_ALLOCATION_2 => 40;

$parallel_processes = get_parallel_processes(
    {
        core_number               => $CORE_NUMBER,
        process_memory_allocation => $PROCESS_MEMORY_ALLOCATION,
        recipe_memory_allocation  => $RECIPE_MEMORY_ALLOCATION_2,
    }
);

## Then constrain the number of parallel processes
Readonly my $OK_PARALLEL_PROCESSES_2 => 5;
is( $parallel_processes, $OK_PARALLEL_PROCESSES_2,
    q{Constrain parallel processes to available cores} );

done_testing();
