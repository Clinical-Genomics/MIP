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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Commands qw{ test_function };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Program::MODULE} => [qw{ SUB_ROUTINE }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::MODULE qw{ SUB_ROUTINE };

diag(   q{Test SUB_ROUTINE from MODULE.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ BASE_COMMAND };

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
    stdoutfile_path => {
        input           => q{stdoutfile.test},
        expected_output => q{1> stdoutfile.test},
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    ARRAY => {
        inputs_ref      => [qw{ TEST_STRING_1 TEST_STRING_2 }],
        expected_output => q{PROGRAM OUTPUT},
    },
    HASH => {
        input_href => {
            key_1 => q{value_1},
            key_2 => q{value_2},
        },

        # Always sorted to an alphabetical order according to ASCII table
        expected_output => q{--hash_arg key_1=value_1 --hash_arg key_2=value_2},
    },
    SCALAR => {
        input           => q{TEST_STRING},
        expected_output => q{PROGRAM_OUTPUT},
    },
);

my %specific_argument = (
    ARRAY => {
        inputs_ref      => [qw{ TEST_STRING_1 TEST_STRING_2 }],
        expected_output => q{PROGRAM OUTPUT},
    },
    HASH => {
        input_href => {
            key_1 => q{value_1},
            key_2 => q{value_2},
        },

        # Always sorted to an alphabetical order according to ASCII table
        expected_output => q{--hash_arg key_1=value_1 --hash_arg key_2=value_2},
    },
    SCALAR => {
        input           => q{TEST_STRING},
        expected_output => q{PROGRAM_OUTPUT},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&SUB_ROUTINE;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

ARGUMENT_HASH_REF:
foreach my $argument_href (@arguments) {

    my @commands = test_function(
        {
            argument_href              => $argument_href,
            do_test_base_command       => 1,
            function_base_commands_ref => \@function_base_commands,
            module_function_cref       => $module_function_cref,
            required_argument_href     => \%required_argument,
        }
    );
}

done_testing();
