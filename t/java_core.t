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
use autodie qw{ :all };
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $EQUALS $SPACE };
use MIP::Test::Commands qw{ test_function };
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
        q{MIP::Language::Java} => [qw{ java_core }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Language::Java qw{ java_core };

diag(   q{Test java_core from Java.pm v}
      . $MIP::Language::Java::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ java };

my %base_argument = (
    FILEHANDLE => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

my %specific_argument = (
    FILEHANDLE => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
    java_jar => {
        input           => q{test.jar},
        expected_output => q{-jar test.jar},
    },
    java_use_large_pages => {
        input           => 1,
        expected_output => q{-XX:-UseLargePages},
    },
    memory_allocation => {
        input           => q{Xmx2g},
        expected_output => q{-Xmx2g},
    },
    picard_use_barclay_parser => {
        input           => 1,
        expected_output => q{-Dpicard.useLegacyParser} . $EQUALS . q{false},
    },
    temp_directory => {
        input           => catdir(qw{ path to temp dir }),
        expected_output => q{-Djava.io.tmpdir=} . catdir(qw{ path to temp dir }),
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&java_core;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

HASHES_OF_ARGUMENTS:
foreach my $argument_href (@arguments) {
    my @commands = test_function(
        {
            argument_href              => $argument_href,
            do_test_base_command       => 1,
            function_base_commands_ref => \@function_base_commands,
            module_function_cref       => $module_function_cref,
        }
    );
}

done_testing();
