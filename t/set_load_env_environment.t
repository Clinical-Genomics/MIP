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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Active_parameter} => [qw{ set_load_env_environment }],
        q{MIP::Test::Fixtures}   => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Active_parameter qw{ set_load_env_environment };

diag(   q{Test set_load_env_environment from Active_parameter.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Given a single load_env environment
my $user_env_name    = q{user_env_name};
my %active_parameter = (
    load_env => {
        mip_rd_dna => {
            method => q{conda},
            mip    => undef,
        },
    },
    environment_name => $user_env_name,
);

set_load_env_environment( { active_parameter_href => \%active_parameter, } );

my %expected_load_env = (
    load_env => {
        $user_env_name => {
            method => q{conda},
            mip    => undef,
        },
    },
);

## Then user supplied environment name should be set in active_parameters load_env
is_deeply(
    $active_parameter{load_env},
    $expected_load_env{load_env},
    q{Set user environment name in load_env}
);

## When multiple load_env environments
$active_parameter{load_env}{another_env} = {
    method => q{conda},
    mip    => undef,
};

trap {
    set_load_env_environment( { active_parameter_href => \%active_parameter, } )
};

like( $trap->stderr, qr/Could\s+not\s+use/xms, q{Throw warning log message} );

## When environment_name is not defined
$active_parameter{environment_name} = undef;

my $return = set_load_env_environment( { active_parameter_href => \%active_parameter, } );

## Then skip sub and return zero
is( $return, 0, q{Skip - no environment name defined} );

done_testing();
