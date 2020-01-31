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
use MIP::Test::Fixtures qw{ test_log test_standard_cli test_mip_hashes };

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
        q{MIP::Environment}    => [qw{ is_binary_in_path }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli test_mip_hashes}],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Environment qw{ is_binary_in_path };

diag(   q{Test is_binary_in_path from Environment.pm v}
      . $MIP::Environment::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Given existing binary
my $binary       = q{ls};

trap {
    is_binary_in_path(
        {
            binary                => $binary,
        }
    )
};

## Then return true
ok( $trap->return, q{Binary is found} );

## Given no existing binary
my $no_binary = q{Nothing to see here};

trap {
    is_binary_in_path(
        {
            binary                => $no_binary,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if binary cannot be found} );

done_testing();
