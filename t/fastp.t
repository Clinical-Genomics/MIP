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
    my %perl_module = ( q{MIP::Program::Fastp} => [qw{ fastp }], );

    test_import( { perl_module_href => \%perl_module, } );
}

## Constants
Readonly my $REQUIRED_LENGTH => 50;
Readonly my $THREADS         => 4;

use MIP::Program::Fastp qw{ fastp };

diag(   q{Test fastp from Fastp.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{fastp};

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
    first_infile_path => {
        input           => q{read_1.fq.gz},
        expected_output => q{--in1 read_1.fq.gz},
    },
    first_outfile_path => {
        input           => q{read_trim_1.fq.gz},
        expected_output => q{--out1 read_trim_1.fq.gz},
    },
    report_html => {
        input           => catfile(qw{ path to report.html }),
        expected_output => q{--html} . $SPACE . catfile(qw{ path to report.html }),
    },
    report_json => {
        input           => catfile(qw{ path to report.json }),
        expected_output => q{--json} . $SPACE . catfile(qw{ path to report.json }),
    },
);

my %specific_argument = (
    detect_pe_adapter => {
        input           => 1,
        expected_output => q{--detect_adapter_for_pe}
    },
    first_infile_path => {
        input           => q{read_1.fq.gz},
        expected_output => q{--in1 read_1.fq.gz},
    },
    first_outfile_path => {
        input           => q{read_trim_1.fq.gz},
        expected_output => q{--out1 read_trim_1.fq.gz},
    },
    interleaved_in => {
        input           => 1,
        expected_output => q{--interleaved_in},
    },
    length_required => {
        input           => $REQUIRED_LENGTH,
        expected_output => q{--length_required} . $SPACE . $REQUIRED_LENGTH,
    },
    low_complexity_filter => {
        input           => 1,
        expected_output => q{--low_complexity_filter},
    },
    overrepresentation_analysis => {
        input           => 1,
        expected_output => q{--overrepresentation_analysis},
    },
    second_infile_path => {
        input           => q{read_2.fq.gz},
        expected_output => q{--in2 read_2.fq.gz},
    },
    second_outfile_path => {
        input           => q{read_trim_2.fq.gz},
        expected_output => q{--out2 read_trim_2.fq.gz},
    },
    threads => {
        input           => $THREADS,
        expected_output => q{--thread 4}
    },
    trim_poly_g => {
        input           => 1,
        expected_output => q{--trim_poly_g},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&fastp;

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
