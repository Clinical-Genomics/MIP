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
        q{MIP::Program::Variant_integrity} => [qw{ variant_integrity_mendel }],
        q{MIP::Test::Fixtures}             => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Variant_integrity qw{ variant_integrity_mendel };

diag(   q{Test variant_integrity_mendel from Variant_integrity v}
      . $MIP::Program::Variant_integrity::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ variant_integrity };

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

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    case_file => {
        input           => q{casefile},
        expected_output => q{--family_file casefile},
    },
    infile_path => {
        input           => q{infile.test},
        expected_output => q{infile.test},
    },
);

my %specific_argument = (
    case_type => {
        input           => q{ped},
        expected_output => q{--family_type ped},
    },
    outfile_path => {
        input           => q{outfile_mendel.txt},
        expected_output => q{--outfile outfile_mendel.txt},
    },
    verbosity => {
        input           => 1,
        expected_output => q{-1},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&variant_integrity_mendel;

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
            required_argument_href     => \%required_argument,
        }
    );
}

done_testing();
