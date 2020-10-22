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
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Commands qw{ test_function };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.02;

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
        q{MIP::Program::Genmod} => [qw{ genmod_annotate }],
        q{MIP::Test::Fixtures}  => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Genmod qw{ genmod_annotate };

diag(   q{Test genmod_annotate from Genmod.pm v}
      . $MIP::Program::Genmod::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

Readonly my $GENOME_BUILD_VERSION => 38;

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
);

my %specific_argument = (
    annotate_region => {
        input           => 1,
        expected_output => q{--annotate_regions},
    },
    cadd_file_paths_ref => {
        inputs_ref =>
          [ catfile(qw{ a test cadd_file_1 }), catfile(qw{ a test cadd_file_2 }) ],
        expected_output => q{--cadd_file}
          . $SPACE
          . catfile(qw{ a test cadd_file_1 })
          . $SPACE
          . q{--cadd_file}
          . $SPACE
          . catfile(qw{ a test cadd_file_2 }),
    },
    genome_build => {
        input           => $THIRTYEIGHT,
        expected_output => q{--genome-build} . $SPACE . $THIRTYEIGHT,
    },
    infile_path => {
        input           => catfile(qw{ a test infile }),
        expected_output => catfile(qw{ a test infile }),
    },
    max_af => {
        input           => q{1},
        expected_output => q{--max_af},
    },
    outfile_path => {
        input           => catfile(qw{ a test outfile }),
        expected_output => q{--outfile} . $SPACE . catfile(qw{ a test outfile }),
    },
    temp_directory_path => {
        input           => catfile(qw{ a test directory }),
        expected_output => q{--temp_dir} . $SPACE . catfile(qw{ a test directory }),
    },
    spidex_file_path => {
        input           => catfile(qw{ a test spidex_file }),
        expected_output => q{--spidex} . $SPACE . catfile(qw{ a test spidex_file }),
    },
    thousand_g_file_path => {
        input           => catfile(qw{ a test thousand_g_file }),
        expected_output => q{--thousand-g}
          . $SPACE
          . catfile(qw{ a test thousand_g_file }),
    },
    verbosity => {
        input           => q{v},
        expected_output => q{-v},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&genmod_annotate;

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
