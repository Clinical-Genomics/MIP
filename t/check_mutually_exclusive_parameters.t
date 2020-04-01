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
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.02;

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
        q{MIP::Active_parameter} => [qw{ check_mutually_exclusive_parameters }],
        q{MIP::Test::Fixtures}   => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Active_parameter qw{ check_mutually_exclusive_parameters };

diag(   q{Test check_mutually_exclusive_parameters from Active_parameter.pm v}
      . $MIP::Active_parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Given no mutually exclusive parameter
my %active_parameter = (
    bwa_mem                                   => 1,
    markduplicates_picardtools_markduplicates => 1,
);

my $is_ok = check_mutually_exclusive_parameters(
    {
        active_parameter_href => \%active_parameter,
    }
);

## Then return TRUE
ok( $is_ok, q{No mutually exclusive options} );

## Given mutually exclusive parameter
$active_parameter{markduplicates_sambamba_markdup} = 1;

trap {
    check_mutually_exclusive_parameters(
        {
            active_parameter_href => \%active_parameter,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if mutually exclusive parameters are found} );
like( $trap->stderr, qr/FATAL/xms, q{Throw fatal mutually exclusive log message} );

done_testing();
