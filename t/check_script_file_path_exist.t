#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Path qw{ make_path };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp qw{ tempdir };
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
use MIP::Test::Fixtures qw{ test_standard_cli };

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
        q{MIP::Script::Setup_script} => [qw{ check_script_file_path_exist }],
        q{MIP::Test::Fixtures}       => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Script::Setup_script qw{ check_script_file_path_exist };

diag(   q{Test check_script_file_path_exist from Setup_script.pm v}
      . $MIP::Script::Setup_script::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $ZERO => 0;
Readonly my $ONE  => 1;

my $test_dir = tempdir( CLEANUP => 1, );

## Given a file prefix and suffix
my $prefix = catfile( $test_dir, q{before} );
my $suffix = q{after};

## When no file path exists
my ( $file_path, $file_name_version ) = check_script_file_path_exist(
    {
        file_path_prefix => $prefix,
        file_path_suffix => $suffix,
    }
);

my $expected_file_path = $prefix . $ZERO . $suffix;

## Then return file path with version zero
is( $file_path, $expected_file_path, q{Initilize file path} );
is( $file_name_version, $ZERO, q{Set version to zero} );

## When file path exists
make_path($file_path);

( $file_path, $file_name_version ) = check_script_file_path_exist(
    {
        file_path_prefix => $prefix,
        file_path_suffix => $suffix,
    }
);
$expected_file_path = $prefix . $ONE . $suffix;

## Then return file path with version one
is( $file_path, $expected_file_path, q{File path version 1} );
is( $file_name_version, $ONE, q{Set version to one} );

done_testing();
