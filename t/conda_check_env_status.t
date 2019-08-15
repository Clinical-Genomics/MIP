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
use Test::Trap qw{ :stderr:output(systemsafe) };    # Define layer to trap

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

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
        q{MIP::Package_manager::Conda} => [qw{ conda_check_env_status }],
        q{MIP::Test::Fixtures}         => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Package_manager::Conda qw{ conda_check_env_status };

diag(   q{Test conda_check_env_status from Conda.pm v}
      . $MIP::Package_manager::Conda::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Given an activate conda env, when running test suite (=TRUE)
# This test will fail otherwise
trap {
    conda_check_env_status( { log => $log, } )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if in conda env } );
like(
    $trap->stderr,
    qr/conda \s+ deactivate' \s+ prior \s+ to/xms,
    q{Throw fatal log message if in conda env}
);

done_testing();
