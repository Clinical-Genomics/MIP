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

## Constants
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Program::Variantcalling::Genmod} => [qw{ genmod_models }],
        q{MIP::Test::Fixtures}                  => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Variantcalling::Genmod qw{ genmod_models };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test genmod_models from Genmod.pm v}
      . $MIP::Program::Variantcalling::Genmod::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ genmod };

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
    infile_path => {
        input           => catfile(qw{ a test infile }),
        expected_output => catfile(qw{ a test infile }),
    },
    case_file => {
        input           => catfile(qw{ a test case_file }),
        expected_output => q{--family_file} . $SPACE . catfile(qw{ a test case_file }),
    },
);

my %specific_argument = (
    case_file => {
        input           => catfile(qw{ a test case_file }),
        expected_output => q{--family_file} . $SPACE . catfile(qw{ a test case_file }),
    },
    case_type => {
        input           => q{mip},
        expected_output => q{--family_type} . $SPACE . q{mip},
    },
    infile_path => {
        input           => catfile(qw{ a test infile }),
        expected_output => catfile(qw{ a test infile }),
    },
    outfile_path => {
        input           => catfile(qw{ a test outfile }),
        expected_output => q{--outfile} . $SPACE . catfile(qw{ a test outfile }),
    },
    reduced_penetrance_file_path => {
        input           => catfile(qw{ a test reduce_pen_file }),
        expected_output => q{--reduced_penetrance}
          . $SPACE
          . catfile(qw{ a test reduce_pen_file}),
    },
    temp_directory_path => {
        input           => catfile(qw{ a test directory }),
        expected_output => q{--temp_dir} . $SPACE . catfile(qw{ a test directory }),
    },
    thread_number => {
        input           => q{2},
        expected_output => q{--processes} . $SPACE . q{2},
    },
    verbosity => {
        input           => q{v},
        expected_output => q{-v},
    },
    vep => {
        input           => q{1},
        expected_output => q{--vep},
    },
    whole_gene => {
        input           => q{1},
        expected_output => q{--whole_gene},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&genmod_models;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

ARGUMENT_HASH_REF:
foreach my $argument_href (@arguments) {
    my @commands = test_function(
        {
            argument_href              => $argument_href,
            required_argument_href     => \%required_argument,
            module_function_cref       => $module_function_cref,
            function_base_commands_ref => \@function_base_commands,
            do_test_base_command       => 1,
        }
    );
}

done_testing();
