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
        q{MIP::Check::Installation} => [qw{ check_mip_executable }],
        q{MIP::Test::Fixtures}      => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Check::Installation qw{ check_mip_executable };

diag(   q{Test check_mip_executable from Installation.pm v}
      . $MIP::Check::Installation::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Given no existing mip binary
my $conda_prefix_path = q{does_not_exists};

my $is_ok = check_mip_executable(
    {
        conda_prefix_path => $conda_prefix_path,
        log               => $log,
    }
);

## Then
ok( $is_ok, q{Found no existing executable} );

## Given an existing mip binary
$conda_prefix_path = catfile( $Bin, qw{ data modules miniconda envs mip_travis } );

trap {
    check_mip_executable(
        {
            conda_prefix_path => $conda_prefix_path,
            log               => $log,
        }
    )
};

## Then throw warning
like( $trap->stderr, qr/This\s+will\s+overwrite/xms, q{Throw overwrite warning} );

done_testing();
