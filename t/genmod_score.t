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


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Program::Genmod} => [qw{ genmod_score }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Genmod qw{ genmod_score };

diag(   q{Test genmod_score from Genmod.pm}
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
    case_file => {
        input           => catfile(qw{ a test case_file }),
        expected_output => q{--family_file} . $SPACE . catfile(qw{ a test case_file }),
    },
    infile_path => {
        input           => catfile(qw{ a test infile }),
        expected_output => catfile(qw{ a test infile }),
    },
    rank_model_file_path => {
        input           => catfile(qw{ a test rank_file }),
        expected_output => q{--score_config} . $SPACE . catfile(qw{ a test rank_file }),
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
    rank_model_file_path => {
        input           => catfile(qw{ a test rank_file }),
        expected_output => q{--score_config} . $SPACE . catfile(qw{ a test rank_file }),
    },
    rank_result => {
        input           => q{1},
        expected_output => q{--rank_results},
    },
    temp_directory_path => {
        input           => catfile(qw{ a test directory }),
        expected_output => q{--temp_dir} . $SPACE . catfile(qw{ a test directory }),
    },
    verbosity => {
        input           => q{v},
        expected_output => q{-v},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&genmod_score;

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
