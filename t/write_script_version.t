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
use Modern::Perl qw{ 2014 };
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = '1.0.0';

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
        q{MIP::Script::Utils}  => [qw{ write_script_version }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Script::Utils qw{ write_script_version };

diag(   q{Test write_script_version from Utils.pm v}
      . $MIP::Script::Utils::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given version, when not supposed to print and exit
my $write_version;

my $is_ok = write_script_version(
    {
        write_version => $write_version,
        version       => $VERSION,
    }
);

## Then return true
ok( $is_ok, q{Did not print version} );

## Given version, when supposed to print and exit
$write_version = 1;

trap {
    write_script_version(
        {
            write_version => $write_version,
            version       => $VERSION,
        }
      )
};

## Then exit and print version
is( $trap->exit, 0, q{Exit and print version } );

done_testing();
