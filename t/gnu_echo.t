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
use autodie qw{ :all };
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $DOUBLE_QUOTE $SPACE };
use MIP::Test::Commands qw{ test_function };
use MIP::Test::Fixtures qw{ test_standard_cli };

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
    my %perl_module = (
        q{MIP::Program::Gnu::Coreutils} => [qw{ gnu_echo }],
        q{MIP::Test::Fixtures}          => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Gnu::Coreutils qw{ gnu_echo };

diag(   q{Test gnu_echo from Coreutils.pm v}
      . $MIP::Program::Gnu::Coreutils::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ echo };

my %base_argument = (
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
    stderrfile_path => {
        input           => q{stderrfile.test},
        expected_output => q{2> stderrfile.test},
    },
    stderrfile_path_append => {
        input           => q{stderrfile.test},
        expected_output => q{2>> stderrfile.test},
    },
    stdoutfile_path_append => {
        input           => q{stdoutfile.test},
        expected_output => q{1>> stdoutfile.test},
    },
);

## Can be duplicated with %base and/or %specific to enable testing of each individual argument
my %required_argument = (
    strings_ref => {
        inputs_ref      => [q{This is my test string}],
        expected_output => $DOUBLE_QUOTE . q{This is my test string} . $DOUBLE_QUOTE,
    },
);

## Specific arguments
my %specific_argument = (
    strings_ref => {
        inputs_ref      => [q{This is my test string}],
        expected_output => $DOUBLE_QUOTE . q{This is my test string} . $DOUBLE_QUOTE,
    },
    outfile_path => {
        input           => q{outfile.test},
        expected_output => q{> outfile.test},
    },
    enable_interpretation => {
        input           => 1,
        expected_output => q{-e},
    },
    no_trailing_newline => {
        input           => 1,
        expected_output => q{-n},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&gnu_echo;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

foreach my $argument_href (@arguments) {
    my @commands = test_function(
        {
            argument_href              => $argument_href,
            required_argument_href     => \%required_argument,
            module_function_cref       => $module_function_cref,
            function_base_commands_ref => \@function_base_commands,
        }
    );
}
done_testing();
