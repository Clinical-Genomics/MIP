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
use Test::Trap;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.03;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

BEGIN {

    use MIP::Test::Fixtures qw{ test_import test_log };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Environment::Path} => [qw{ get_conda_path }],
        q{MIP::Test::Fixtures}    => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Environment::Path qw{ get_conda_path };

diag(   q{Test get_conda_path from Parameter.pm v}
      . $MIP::Environment::Path::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Given a sub call when conda should be installed
my $is_ok = get_conda_path( {} );

## Then return true
ok( $is_ok, q{Found conda path} );

# Given a sub call when conda should not be installed
trap {
    get_conda_path( { bin_file => q{not_a_conda_bin} } )
};
## Then print fatal log message and exit
ok( $trap->exit, q{Exit when no conda is found } );
like( $trap->stderr, qr/Failed \s+ to \s+ find \s+ path/xms, q{Throw fatal log message} );

done_testing();
