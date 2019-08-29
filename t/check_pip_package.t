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
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::File::Format::Yaml qw{ load_yaml };
use MIP::Test::Fixtures qw{ test_standard_cli };

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
        q{MIP::Package_manager::Pip} => [qw{ check_pip_package }],
        q{MIP::Test::Fixtures}       => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Package_manager::Pip qw{ check_pip_package };

diag(   q{Test check_pip_package from Pip.pm v}
      . $MIP::Package_manager::Pip::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

Readonly my $TRAVIS_TESTS => 4;

my $conda_environment = q{mip_travis_svdb};
my $conda_prefix_path = catfile( $Bin, qw{ data modules miniconda } );


my %installation = load_yaml(
    {
        yaml_file => catfile($Bin, qw{ .. templates mip_install_rd_dna_config_-1.0-.yaml }),
    }
);
my $rhocall_version = $installation{emip}{shell}{rhocall}{version};
my $svdb_version = $installation{esvdb}{shell}{svdb}{version};

## Only run on travis
SKIP: {
    skip q{No control of local environment names}, $TRAVIS_TESTS, if not $ENV{q{TRAVIS}};

    ## Given a conda environment and pip package

    my $found = check_pip_package(
        {
            conda_environment => $conda_environment,
            conda_prefix_path => $conda_prefix_path,
            package           => q{svdb},
            version           => $svdb_version,
        }
    );
    ## Then return 1
    ok( $found, q{Found correct package and version in conda env} );

    ## Given an existing package but the wrong version
    my $not_found = check_pip_package(
        {
            conda_environment => $conda_environment,
            conda_prefix_path => $conda_prefix_path,
            package           => q{svdb},
            version           => q{1.2.3.4},
        }
    );
    ## Then return undef
    is( $not_found, undef, q{Not found in env} );

    ## Given a non existing package
    $not_found = check_pip_package(
        {
            conda_environment => $conda_environment,
            conda_prefix_path => $conda_prefix_path,
            package           => q{testing_mip},
            version           => $rhocall_version,
        }
    );
    ## Then return undef
    is( $not_found, undef, q{No package found in env} );

    ## Given a pip package without a conda env, when package exist
    my $is_ok = check_pip_package(
        {
            package => q{rhocall},
            version => $rhocall_version,
        }
    );
    ## Then return true
    ok( $is_ok, q{Checked for pip package in current environment} );
}

## Given a non existing conda environment
my $not_a_conda_env = check_pip_package(
    {
        conda_environment => q{testing_mip},
        conda_prefix_path => $conda_prefix_path,
        package           => q{svdb},
        version           => $svdb_version,
    }
);

## Then return undef
is( $not_a_conda_env, undef, q{Not a conda env} );

## Given an non existing package in the current environment
my $not_found = check_pip_package(
    {
        package => q{non_existing_package},
    }
);

## Then return undef
is( $not_found, undef, q{Not found in current environment} );

done_testing();
