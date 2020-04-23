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
use Test::Trap qw{ :stderr:output(systemsafe) };

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.05;

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
        q{MIP::Vep}            => [qw{ check_vep_api_cache_versions }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Vep qw{ check_vep_api_cache_versions };

diag(   q{Test check_vep_api_cache_versions from Vep.pm v}
      . $MIP::Vep::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $test_dir = File::Temp->newdir();

## Creates log object
my $log = test_log( {} );

## Given matching vep API and cache version
my $vep_binary_path =
  catdir( $Bin, qw{ data modules miniconda envs mip_travis bin vep } );
my $vep_directory_cache =
  catdir( $Bin, qw{ data references ensembl-tools-data-99 cache } );

## When comparing API and cache version
my $match = check_vep_api_cache_versions(
    {
        vep_directory_cache => $vep_directory_cache,
        vep_binary_path     => $vep_binary_path,
    }
);

## Then return true
ok( $match, q{Return on matching versions} );

## Given non matching API and cache
$vep_directory_cache = catdir( $Bin, qw{ data references ensembl-tools-data-99 cache2 } );

## When comparing API and cache version
trap {
    check_vep_api_cache_versions(
        {
            log                 => $log,
            vep_directory_cache => $vep_directory_cache,
            vep_binary_path     => $vep_binary_path,
        }
    )
};

## Then exit and print fatal error message
like( $trap->stderr, qr/FATAL/xms, q{Throw FATAL log message} );
ok( $trap->exit, q{Exit on non matching versions} );

## Given a cache folder lacking the homo_sapiens subdirectory
$vep_directory_cache =
  catdir( $Bin, qw{ data modules miniconda envs test_env ensembl-vep } );

## When trying to retireve the cache versions
trap {
    check_vep_api_cache_versions(
        {
            log                 => $log,
            vep_directory_cache => $vep_directory_cache,
            vep_binary_path     => $vep_binary_path,
        }
    )
};

## Then return and print warning message
ok( $trap->return, q{Return on unknown cache version} );
like( $trap->stderr, qr/WARN/xms, q{Warn for unknown VEP cache version} );

## Given a direcory that lacks a working vep bin
$vep_binary_path = catdir( $Bin, qw{ data modules miniconda envs test_env bin vep} );
$vep_directory_cache =
  catdir( $Bin, qw{ data modules miniconda envs test_env ensembl-tools-91 cache } );

## When trying to retrieve API version
trap {
    check_vep_api_cache_versions(
        {
            log                 => $log,
            vep_directory_cache => $vep_directory_cache,
            vep_binary_path     => $vep_binary_path,
        }
    )
};

## Then return and print warning message
ok( $trap->return, q{Return on unknown API version} );
like( $trap->stderr, qr/WARN/xms, q{Warn for unknown VEP api version} );

done_testing();
