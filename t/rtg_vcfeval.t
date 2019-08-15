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
use Modern::Perl qw{ 2014 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Commands qw{ test_function };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.03;

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
        q{MIP::Program::Rtg}   => [qw{ rtg_vcfeval }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Rtg qw{ rtg_vcfeval };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test rtg_vcfeval from Rtg.pm v}
      . $MIP::Program::Rtg::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ rtg vcfeval };

my %base_argument = (
    FILEHANDLE => {
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
    baselinefile_path => {
        input           => catfile(qw{path to baselinefile}),
        expected_output => q{--baseline=} . catfile(qw{path to baselinefile}),
    },
    callfile_path => {
        input           => catfile(qw{path to callfile}),
        expected_output => q{--calls=} . catfile(qw{path to callfile}),
    },
    eval_region_file_path => {
        input           => catfile(qw{path to eval_regionsfile}),
        expected_output => q{--evaluation-regions=}
          . catfile(qw{path to eval_regionsfile}),
    },
    outputdirectory_path => {
        input           => catfile(qw{path to outputdirectory_path}),
        expected_output => q{--output=} . catfile(qw{path to outputdirectory_path}),
    },
    sdf_template_file_path => {
        input           => catfile(qw{path to sdf_template_file_path}),
        expected_output => q{--template=} . catfile(qw{path to sdf_template_file_path}),
    },
);

my %specific_argument = (
    all_record => {
        input           => 1,
        expected_output => q{--all-records},
    },
    baselinefile_path => {
        input           => catfile(qw{path to baselinefile}),
        expected_output => q{--baseline=} . catfile(qw{path to baselinefile}),
    },
    bed_regionsfile_path => {
        input           => catfile(qw{path to bed_regionsfile_path}),
        expected_output => q{--bed-regions=} . catfile(qw{path to bed_regionsfile_path}),
    },
    callfile_path => {
        input           => catfile(qw{path to callfile}),
        expected_output => q{--calls=} . catfile(qw{path to callfile}),
    },
    eval_region_file_path => {
        input           => catfile(qw{path to eval_regionsfile}),
        expected_output => q{--evaluation-regions=}
          . catfile(qw{path to eval_regionsfile}),
    },
    outputdirectory_path => {
        input           => catfile(qw{path to outputdirectory_path}),
        expected_output => q{--output=} . catfile(qw{path to outputdirectory_path}),
    },
    output_mode => {
        input           => q{annotate},
        expected_output => q{--output-mode=annotate},
    },
    sample_id => {
        input           => q{sample_1},
        expected_output => q{--sample=sample_1},
    },
    sdf_template_file_path => {
        input           => catfile(qw{path to sdf_template_file_path}),
        expected_output => q{--template=} . catfile(qw{path to sdf_template_file_path}),
    },
    thread_number => {
        input           => 2,
        expected_output => q{--threads=2},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&rtg_vcfeval;

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
