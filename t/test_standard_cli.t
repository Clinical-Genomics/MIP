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
    my %perl_module = ( q{MIP::Test::Fixtures} => [qw{ test_standard_cli }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Test::Fixtures qw{ test_standard_cli };
use MIP::Unix::System qw{ system_cmd_call };

diag(   q{Test test_standard_cli from Fixtures.pm v}
      . $MIP::Test::Fixtures::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a cli
my $is_ok = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Then return true
ok( $is_ok, q{Generated standard cli} );

## Given verbose and version (this scripts)
my $command_version_string = qq{perl $PROGRAM_NAME -v};
my %return = system_cmd_call( { command_string => $command_version_string, } );

my $expected_version_return = join $SPACE, @{ $return{output} };

## Then show version
like( $expected_version_return, qr/test_standard_cli.t \s+ 1/xms, q{Show version} );

## Given verbose and help
my $command_help_string = qq{perl $PROGRAM_NAME -h};
%return = system_cmd_call( { command_string => $command_help_string, } );

my $expected_help_return = join $SPACE, @{ $return{output} };

## Then show help text
like( $expected_help_return, qr/test_standard_cli.t \s+ [[]options[]]/xms, q{Show help} );

done_testing();
