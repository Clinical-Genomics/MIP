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


## Constants
Readonly my $COVERAGE => 90;

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Program::Gatk}  => [qw{ gatk_combinevariants }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Gatk qw{ gatk_combinevariants };

diag(   q{Test gatk_combinevariants from Gatk}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ gatk3 };

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
        inputs_ref => [qw{ var_1.vcf var_2.vcf var_3.vcf }],
        expected_output =>
          q{--variant: var_1.vcf --variant: var_2.vcf --variant: var_3.vcf},
    },
    outfile_path => {
        input           => catfile(qw{ dir outfile.vcf }),
        expected_output => q{--outputFile} . $SPACE . catfile(qw{ dir outfile.vcf }),
    },
    referencefile_path => {
        input           => catfile(qw{reference_dir human_genome_build.fasta }),
        expected_output => q{--reference_sequence}
          . $SPACE
          . catfile(qw{reference_dir human_genome_build.fasta }),
    },
);

my %specific_argument = (
    downsample_to_coverage => {
        input           => $COVERAGE,
        expected_output => q{--downsample_to_coverage } . $COVERAGE,
    },
    exclude_nonvariants => {
        input           => 1,
        expected_output => q{--excludeNonVariants},
    },
    gatk_disable_auto_index_and_file_lock => {
        input           => 1,
        expected_output => q{--disable_auto_index_creation_and_locking_when_reading_rods},
    },
    genotype_merge_option => {
        input           => q{PRIORITIZE},
        expected_output => q{--genotypemergeoption PRIORITIZE},
    },
    intervals_ref => {
        inputs_ref      => [qw{ chr1 chr2 chr3 }],
        expected_output => q{--intervals chr1 --intervals chr2 --intervals chr3},
    },
    logging_level => {
        input           => q{INFO},
        expected_output => q{--logging_level INFO},
    },
    pedigree => {
        input           => catfile(qw{ dir pedigree.fam }),
        expected_output => q{--pedigree} . $SPACE . catfile(qw{ dir pedigree.fam }),
    },
    pedigree_validation_type => {
        input           => q{SILENT},
        expected_output => q{--pedigreeValidationType SILENT},
    },
    prioritize_caller => {
        input           => q{priority_list},
        expected_output => q{--rod_priority_list priority_list},
    },

);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&gatk_combinevariants;

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

## Base arguments
@function_base_commands = qw{ gatk3 java };

my %specific_java_argument = (
    java_jar => {
        input           => q{gatk.jar},
        expected_output => q{-jar gatk.jar},
    },
);

## Test both base and function specific arguments
@arguments = ( \%specific_java_argument );

ARGUMENT_HASH_REF:
foreach my $argument_href (@arguments) {
    my @commands = test_function(
        {
            argument_href              => $argument_href,
            base_commands_index        => 1,
            do_test_base_command       => 1,
            function_base_commands_ref => \@function_base_commands,
            module_function_cref       => $module_function_cref,
            required_argument_href     => \%required_argument,
        }
    );
}
done_testing();
