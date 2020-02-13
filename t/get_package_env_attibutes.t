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
use MIP::Test::Fixtures qw{ test_mip_hashes test_standard_cli };

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
        q{MIP::Active_parameter} => [qw{ get_package_env_attributes }],
        q{MIP::Test::Fixtures}   => [qw{ test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Active_parameter qw{ get_package_env_attributes };

diag(   q{Test get_package_env_attributes from Active_parameter.pm v}
      . $MIP::Active_parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given package attributes and a package name
my %active_parameter = test_mip_hashes( { mip_hash_name => q{active_parameter}, } );

my ( $env_name, $env_method ) = get_package_env_attributes(
    {
        load_env_href => $active_parameter{load_env},
        package_name  => q{mip},
    }
);

my $expected_env_name   = q{test};
my $expected_env_method = q{conda};

## Then return environment name and method for package
is( $expected_env_name,   $env_name,   q{Got environment name for package} );
is( $expected_env_method, $env_method, q{Got environment method for package} );

## Given a not existing package name
my $is_ok = get_package_env_attributes(
    {
        load_env_href => $active_parameter{load_env},
        package_name  => q{package_does_not_exist},
    }
);

## Then return false
is( $is_ok, undef, q{Return if no existing package name} );

done_testing();
