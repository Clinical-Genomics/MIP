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
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Commands qw{ test_function };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Program::Gnu::Coreutils} => [qw{ gnu_split }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Gnu::Coreutils qw{ gnu_split };

diag(   q{Test gnu_split from Coreutils.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

Readonly my $NUMBER_OF_FILES => 12;

## Base arguments
my @function_base_commands = qw{ split };

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
);

## Can be duplicated with %base and/or %specific to enable testing of each individual argument
my %required_argument = (
    infile_path => {
        input           => q{infile.test},
        expected_output => q{infile.test},
    },
);

## Specific arguments
my %specific_argument = (
    files => {
        input           => $NUMBER_OF_FILES,
        expected_output => q{--number=l/} . $NUMBER_OF_FILES,
    },
    infile_path => {
        input           => q{infile.test},
        expected_output => q{infile.test},
    },
    lines => {
        input           => 10,
        expected_output => q{--lines=10},
    },
    numeric_suffixes => {
        input           => 1,
        expected_output => q{--numeric-suffixes},
    },
    prefix => {
        input           => q{test},
        expected_output => q{test},
    },
    quiet => {
        input           => 1,
        expected_output => q{--quiet},
    },
    suffix_length => {
        input           => 2,
        expected_output => q{--suffix-length=2},
    },
    verbose => {
        input           => 1,
        expected_output => q{--verbose},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&gnu_split;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

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
