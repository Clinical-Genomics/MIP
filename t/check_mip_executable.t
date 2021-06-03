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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import test_log };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Install}        => [qw{ check_mip_executable }],
        q{MIP::Test::Fixtures} => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Install qw{ check_mip_executable };

diag(   q{Test check_mip_executable from Install.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

test_log( {} );

## Given no existing mip binary
my $conda_environment_path = q{does_not_exists};

## When checking for executable
my $is_not_found = check_mip_executable(
    {
        conda_environment_path => $conda_environment_path,
    }
);

## Then return true
ok( $is_not_found, q{Found no existing executable} );

## Given an existing mip binary
$conda_environment_path = catfile( $Bin, qw{ data modules miniconda envs mip_ci } );

## When checking for executable
trap {
    check_mip_executable(
        {
            conda_environment_path => $conda_environment_path,
        }
    )
};

## Then throw warning
like( $trap->stderr, qr/This\s+will\s+overwrite/xms, q{Throw overwrite warning} );

done_testing();
