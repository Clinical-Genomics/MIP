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
our $VERSION = 1.01;

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
        q{MIP::Program::Variantcalling::Picardtools} =>
          [qw{ picardtools_genotypeconcordance }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Variantcalling::Picardtools qw{ picardtools_genotypeconcordance };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test picardtools_genotypeconcordance from Variantcalling::Picardtools.pm v}
      . $MIP::Program::Variantcalling::Picardtools::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $MIN_GQ    => 20;
Readonly my $MIN_DEPTH => 10;

## Base arguments
my @function_base_commands = qw{ GenotypeConcordance };

my %base_argument = (
    stderrfile_path => {
        input           => q{stderrfile.test},
        expected_output => q{2> stderrfile.test},
    },
    FILEHANDLE => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    infile_path => {
        input           => catfile(qw{ dir infile.vcf }),
        expected_output => q{CALL_VCF=} . catfile(qw{ dir infile.vcf }),
    },
    truth_file_path => {
        input           => catfile(qw{ dir truth.vcf }),
        expected_output => q{TRUTH_VCF=} . catfile(qw{ dir truth.vcf }),
    },
    outfile_prefix_path => {
        input           => catfile(qw{ dir outfile_prefix.vcf }),
        expected_output => q{OUTPUT=} . catfile(qw{ dir outfile_prefix.vcf }),
    },
    truth_sample => {
        input           => catfile(qw{ dir truth_sample.vcf }),
        expected_output => q{TRUTH_SAMPLE=} . catfile(qw{ dir truth_sample.vcf }),
    },
    call_sample => {
        input           => catfile(qw{ dir call_sample.vcf }),
        expected_output => q{CALL_SAMPLE=} . catfile(qw{ dir call_sample.vcf }),
    },
    referencefile_path => {
        input           => catfile(qw{ references grch37_homo_sapiens_-d5-.fasta }),
        expected_output => q{R=}
          . catfile(qw{ references grch37_homo_sapiens_-d5-.fasta }),
    },
);

my %specific_argument = (
    intervals_ref => {
        inputs_ref      => [qw{ chr1 chr2 }],
        expected_output => q{INTERVALS=chr1 INTERVALS=chr2},
    },
    min_genotype_quality => {
        input           => $MIN_GQ,
        expected_output => q{MIN_GQ=} . $MIN_GQ,
    },
    min_depth => {
        input           => $MIN_DEPTH,
        expected_output => q{MIN_DP=} . $MIN_DEPTH,
    },
    infile_path => {
        input           => catfile(qw{ dir infile.vcf }),
        expected_output => q{CALL_VCF=} . catfile(qw{ dir infile.vcf }),
    },
    truth_file_path => {
        input           => catfile(qw{ dir truth.vcf }),
        expected_output => q{TRUTH_VCF=} . catfile(qw{ dir truth.vcf }),
    },
    outfile_prefix_path => {
        input           => catfile(qw{ dir outfile_prefix.vcf }),
        expected_output => q{OUTPUT=} . catfile(qw{ dir outfile_prefix.vcf }),
    },
    truth_sample => {
        input           => catfile(qw{ dir truth_sample.vcf }),
        expected_output => q{TRUTH_SAMPLE=} . catfile(qw{ dir truth_sample.vcf }),
    },
    call_sample => {
        input           => catfile(qw{ dir call_sample.vcf }),
        expected_output => q{CALL_SAMPLE=} . catfile(qw{ dir call_sample.vcf }),
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&picardtools_genotypeconcordance;

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

## Test specfic java core options
## Base arguments
@function_base_commands = qw{ java };

my %specific_java_argument = (
    java_jar => {
        input           => q{picard.jar},
        expected_output => q{-jar picard.jar},
    },
);

## Test both base and function specific arguments
@arguments = ( \%specific_java_argument );

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
