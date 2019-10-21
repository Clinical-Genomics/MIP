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
our $VERSION = 1.00;

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
        q{MIP::Program::Variantcalling::Gatk} => [qw{ gatk_varianteval }],
        q{MIP::Test::Fixtures}                => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Variantcalling::Gatk qw{ gatk_varianteval };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test gatk_varianteval from Gatk v}
      . $MIP::Program::Variantcalling::Gatk::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ gatk VariantEval };

my %base_argument = (
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    infile_paths_ref => {
        inputs_ref      => [qw{ var_1.vcf var_2.vcf var_3.vcf }],
        expected_output => q{--eval var_1.vcf --eval var_2.vcf --eval var_3.vcf},
    },
    outfile_path => {
        input           => catfile(qw{ path_to_analysis_dir outfile.vcf }),
        expected_output => q{--output}
          . $SPACE
          . catfile(qw{ path_to_analysis_dir outfile.vcf }),
    },
    referencefile_path => {
        input           => catfile(qw{reference_dir human_genome_build.fasta }),
        expected_output => q{--reference}
          . $SPACE
          . catfile(qw{reference_dir human_genome_build.fasta }),
    },
);

my %specific_argument = (
    dbsnp_file_path => {
        input           => catfile(qw{ dir grch37_dbsnp_-138_esa_129-.vcf}),
        expected_output => q{--dbsnp}
          . $SPACE
          . catfile(qw{ dir grch37_dbsnp_-138_esa_129-.vcf}),
    },
    indel_gold_standard_file_path => {
        input => catfile(qw{ dir grch37_mills_and_1000g_indels_-gold_standard-.vcf}),
        expected_output => q{--gold-standard}
          . $SPACE
          . catfile(qw{ dir grch37_mills_and_1000g_indels_-gold_standard-.vcf}),
    },
    infile_paths_ref => {
        inputs_ref      => [qw{ var_1.vcf var_2.vcf var_3.vcf }],
        expected_output => q{--eval var_1.vcf --eval var_2.vcf --eval var_3.vcf},
    },
    intervals_ref => {
        inputs_ref      => [qw{ chr1 chr2 chr3 }],
        expected_output => q{--intervals chr1 --intervals chr2 --intervals chr3},
    },
    outfile_path => {
        input           => catfile(qw{ path_to_analysis_dir outfile.vcf }),
        expected_output => q{--output}
          . $SPACE
          . catfile(qw{ path_to_analysis_dir outfile.vcf }),
    },
    pedigree => {
        input           => catfile(qw{ dir pedigree.fam }),
        expected_output => q{--pedigree} . $SPACE . catfile(qw{ dir pedigree.fam }),
    },
    verbosity => {
        input           => q{INFO},
        expected_output => q{--verbosity INFO},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&gatk_varianteval;

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
