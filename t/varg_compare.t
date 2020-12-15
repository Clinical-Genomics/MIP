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
        q{MIP::Program::Varg} => [qw{ varg_compare }],

    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Varg qw{ varg_compare };

diag(   q{Test varg_compare from Varg.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ varg compare };

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
    infile_path_truth => {
        input           => q{test_infile_truth.vcf},
        expected_output => q{--truth-set} . $SPACE . q{test_infile_truth.vcf}
    },
    infile_path_vcf => {
        input           => q{test_infile_vcf.vcf},
        expected_output => q{--variants} . $SPACE . q{test_infile_vcf.vcf}
    }
);

my %specific_argument = (
    infile_path_truth => {
        input           => q{test_infile_truth.vcf},
        expected_output => q{--truth-set} . $SPACE . q{test_infile_truth.vcf}
    },
    infile_path_vcf => {
        input           => q{test_infile_vcf.vcf},
        expected_output => q{--variants} . $SPACE . q{test_infile_vcf.vcf}
    }
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&varg_compare;

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
