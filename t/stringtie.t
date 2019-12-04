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
        q{MIP::Program::Stringtie} => [qw{ stringtie }],
        q{MIP::Test::Fixtures}     => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Stringtie qw{ stringtie };

diag(   q{Test stringtie from Stringtie.pm v}
      . $MIP::Program::Stringtie::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $MIN_COVERAGE => 5;
Readonly my $THREADS      => 16;

## Base arguments
my @function_base_commands = qw{ stringtie };

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
    gtf_reference_path => {
        input           => catfile(qw{ path to gtf }),
        expected_output => q{-G} . $SPACE . catfile(qw{ path to gtf }),
    },
    infile_path => {
        input           => catfile(qw{ path to infile.bam }),
        expected_output => catfile(qw{ path to infile.bam }),
    },
    outfile_path => {
        input           => catfile(qw{ path to gtf }),
        expected_output => q{-G} . $SPACE . catfile(qw{ path to gtf }),
    },
);

my %specific_argument = (
    cov_ref_transcripts_outfile_path => {
        input           => catfile(qw{ path to cov.gtf }),
        expected_output => q{-C} . $SPACE . catfile(qw{ path to cov.gtf }),
    },
    gene_abundance_outfile_path => {
        input           => catfile(qw{ path to abound.tab }),
        expected_output => q{-A} . $SPACE . catfile(qw{ path to abound.tab }),
    },
    gtf_reference_path => {
        input           => catfile(qw{ path to gtf }),
        expected_output => q{-G} . $SPACE . catfile(qw{ path to gtf }),
    },
    infile_path => {
        input           => catfile(qw{ path to infile.bam }),
        expected_output => catfile(qw{ path to infile.bam }),
    },
    junction_reads => {
        input           => q{2.5},
        expected_output => q{-j 2.5},
    },
    library_type => {
        input           => q{forward_stranded},
        expected_output => q{--fr},
    },
    minimum_coverage => {
        input           => $MIN_COVERAGE,
        expected_output => q{-c} . $SPACE . $MIN_COVERAGE,
    },
    outfile_path => {
        input           => catfile(qw{ path to gtf }),
        expected_output => q{-G} . $SPACE . catfile(qw{ path to gtf }),
    },
    threads => {
        input           => $THREADS,
        expected_output => q{-p} . $SPACE . $THREADS,
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&stringtie;

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
