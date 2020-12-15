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

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Active_parameter} => [qw{ check_recipe_mode }],
        q{MIP::Test::Fixtures}   => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Active_parameter qw{ check_recipe_mode };

diag(   q{Test check_recipe_mode from Active_parameter.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Creates log object
my $log = test_log( {} );

## Given recipes when correct recipe modes
my %parameter =
  ( cache => { recipe => [ qw{ bwa_mem fastqc_ar genmod}, ], }, );
my %active_parameter = (
    bwa_mem   => 1,
    fastqc_ar => 0,
    genmod    => 2,
);

my $is_ok = check_recipe_mode(
    {
        active_parameter_href => \%active_parameter,
        parameter_href        => \%parameter,
    }
);

## Then all is ok
ok( $is_ok, q{All recipe modes passed} );

## Given recipes when a recipe has a not allowed value
$active_parameter{bwa_mem} = q{not allowed value};

trap {
    check_recipe_mode(
        {
            active_parameter_href => \%active_parameter,
            parameter_href        => \%parameter,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if recipe mode is not allowed} );
like( $trap->stderr, qr/FATAL/xms,
    q{Throw fatal log message if recipe mode is not allowed} );

## Given a incorrect recipe name
# Reinstate correct mode
$active_parameter{bwa_mem} = 1;
push @{ $parameter{cache}{recipe} }, q{not_a_recipe};

trap {
    check_recipe_mode(
        {
            active_parameter_href => \%active_parameter,
            parameter_href        => \%parameter,
        }
    )
};

## Then exit and throw FATAL log message
is( $trap->leaveby, q{die}, q{Exit if recipe key is not found} );

done_testing();
