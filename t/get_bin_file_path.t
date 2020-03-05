#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2018 };

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_mip_hashes test_standard_cli };

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
        q{MIP::Environment::Path} => [qw{ get_bin_file_path }],
        q{MIP::Test::Fixtures}    => [qw{ test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Environment::Path qw{ get_bin_file_path };

diag(   q{Test get_bin_file_path from Path.pm v}
      . $MIP::Environment::Path::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a conda bin and bwa bin file path, when exists
my %active_parameter = test_mip_hashes( { mip_hash_name => q{active_parameter}, } );
my $conda_path       = catfile( $Bin, qw{ data modules miniconda } );
my $recipe_name      = q{bwa_mem};

my %environment = ( $recipe_name => [qw{ source activate test }] );
my $bin_file    = q{bwa};

my ( $bin_file_path, $conda_environment ) = get_bin_file_path(
    {
        bin_file         => $bin_file,
        conda_path       => $conda_path,
        environment_href => \%environment,
        environment_key  => $recipe_name,
    }
);

my $expected_bin_file_path =
  catfile( $Bin, qw{ data modules miniconda envs test bin bwa } );
my $expected_conda_env = q{test};

## Then return absolute path for binary file and the corresponding environment
is( $bin_file_path, $expected_bin_file_path, q{Got absolute path for binary file} );

is( $conda_environment, $expected_conda_env, q{Got conda env for binary file} );

done_testing();
