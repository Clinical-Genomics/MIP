#! /usr/bin/env perl

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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import test_log };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Active_parameter} => [qw{ set_conda_path }],
        q{MIP::Test::Fixtures}   => [qw{ test_log }],

    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Environment::Path qw{ get_conda_path };
use MIP::Active_parameter qw{ set_conda_path };

diag(   q{Test set_conda_path from Active_parameter.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

test_log( {} );

## Given a set conda path in active parameters
my %active_parameter;
my $conda_path = q{a_conda_path};
$active_parameter{conda_path} = $conda_path;

## Given a set environment_name in active parameters
my $environment_name = q{an_env_name};

## When the conda path is already set active_parameter
set_conda_path(
    {
        active_parameter_href => \%active_parameter,
        environment_name      => $environment_name,
    }
);

## Then the conda path supplied should be set
is( $active_parameter{conda_path}, $conda_path, q{Set supplied conda path} );

## Then the conda_environment_path should be set
my $expected_conda_environment_path = catdir( $conda_path, q{envs}, $environment_name );

is(
    $active_parameter{conda_environment_path},
    $expected_conda_environment_path,
    q{Set supplied conda_environment_path}
);

## Given no set conda_path
delete $active_parameter{conda_path};

## When the conda path is not already set active_parameter
set_conda_path(
    {
        active_parameter_href => \%active_parameter,
        environment_name      => $environment_name,
    }
);

## Then the conda path should be set dynamically
my $expected_conda_path = get_conda_path( {} );

is( $active_parameter{conda_path}, $expected_conda_path, q{Set conda path dynamically} );

done_testing();
