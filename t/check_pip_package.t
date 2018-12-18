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
use Modern::Perl qw{ 2014 };
use Readonly;

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

## Given a conda environment and pip package, but not installed via conda
my $conda_environment = q{mip_travis_svdb};
my $svdb_version      = 1.0;
my $conda_prefix_path = catfile($Bin, qw{ data modules miniconda });

my $is_not_ok_conda = check_pip_package(
    {
        conda_environment => $conda_environment,
     conda_prefix_path => $conda_prefix_path,
        package           => q{svdb},
        version           => $svdb_version,
    }
);

## Then return undef
like( $is_not_ok_conda, qr/EnvironmentLocationNotFound/, q{Checked pip package in conda env} );

## Given a conda environment and pip package, when package exist
my $is_ok = check_pip_package(
    {
        package => q{conda},
    }
);

## Then return true
ok( $is_ok, q{Checked pip package in base env} );

done_testing();
