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
        q{MIP::Program::Variantcalling::Gffcompare} => [qw{ gffcompare }],
        q{MIP::Test::Fixtures}                      => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Variantcalling::Gffcompare qw{ gffcompare };

diag(   q{Test gffcompare from Gffcompare.pm v}
      . $MIP::Program::Variantcalling::Gffcompare::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ gffcompare };

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
        expected_output => q{-r} . $SPACE . catfile(qw{ path to gtf }),
    },
    infile_paths_ref => {
        inputs_ref => [ catfile(qw{ path infile.bam }), catfile(qw{ path infile2.bam }) ],
        expected_output =>
          [ catfile(qw{ path infile.bam }), catfile(qw{ path infile2.bam }) ],
    },
    outfile_path_prefix => {
        input           => catfile(qw{ path to prefix }),
        expected_output => q{-o} . $SPACE . catfile(qw{ path to prefix }),
    },
);

my %specific_argument = (
    genome_sequence_path => {
        input           => catfile(qw{ path to fasta }),
        expected_output => q{-s} . $SPACE . catfile(qw{ path to fasta }),
    },
    gtf_reference_path => {
        input           => catfile(qw{ path to gtf }),
        expected_output => q{-r} . $SPACE . catfile(qw{ path to gtf }),
    },
    ignore_non_overlapping_ref => {
        input           => 1,
        expected_output => q{-R},
    },
    infile_paths_ref => {
        inputs_ref => [ catfile(qw{ path infile.bam }), catfile(qw{ path infile2.bam }) ],
        expected_output => catfile(qw{ path infile.bam })
          . $SPACE
          . catfile(qw{ path infile2.bam }),
    },
    outfile_path_prefix => {
        input           => catfile(qw{ path to prefix }),
        expected_output => q{-o} . $SPACE . catfile(qw{ path to prefix }),
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&gffcompare;

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
