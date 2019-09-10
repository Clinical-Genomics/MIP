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
our $VERSION = 1.01;

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
        q{MIP::Program::Base::Gatk} => [qw{ gatk_base }],
        q{MIP::Test::Fixtures}      => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Base::Gatk qw{ gatk_base };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test gatk_base from Base::Gatk.pm v}
      . $MIP::Program::Base::Gatk::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $DOWNSAMPLE_TO_COVERAGE               => 1000;
Readonly my $STATIC_QUANTIZED_QUALS_MINIMUM_LEVEL => 10;
Readonly my $STATIC_QUANTIZED_QUALS_MAX_LEVEL     => 20;

## Base arguments
my @function_base_commands = qw{ --analysis_type HaplotypeCaller };

my %base_argument = (
    FILEHANDLE => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    analysis_type => {
        input           => q{HaplotypeCaller},
        expected_output => q{--analysis_type HaplotypeCaller},
    },
    referencefile_path => {
        input           => catfile(qw{ reference_dir human_genome_build.fasta }),
        expected_output => q{--reference_sequence}
          . catfile(qw{reference_dir human_genome_build.fasta }),
    },
);

my %specific_argument = (
    analysis_type => {
        input           => q{HaplotypeCaller},
        expected_output => qw{ HaplotypeCaller },
    },
    base_quality_score_recalibration_file => {
        input           => catfile(qw{ dir infile.bsqr }),
        expected_output => q{--BQSR } . catfile(qw{ dir infile.bsqr }),
    },
    disable_indel_qual => {
        input           => 1,
        expected_output => q{--disable_indel_quals},
    },
    downsample_to_coverage => {
        input           => $DOWNSAMPLE_TO_COVERAGE,
        expected_output => q{--downsample_to_coverage } . $DOWNSAMPLE_TO_COVERAGE,
    },
    gatk_disable_auto_index_and_file_lock => {
        input           => 1,
        expected_output => q{--disable_auto_index_creation_and_locking_when_reading_rods},
    },
    intervals_ref => {
        inputs_ref      => [qw{ chr1 chr2}],
        expected_output => q{--intervals chr1 --intervals chr2},
    },
    num_cpu_threads_per_data_thread => {
        input           => 2,
        expected_output => q{--num_cpu_threads_per_data_thread 2},
    },
    logging_level => {
        input           => q{INFO},
        expected_output => q{--logging_level INFO},
    },
    pedigree => {
        input           => catfile(qw{ dir ped.fam }),
        expected_output => q{--pedigree } . catfile(qw{ dir ped.fam }),
    },
    pedigree_validation_type => {
        input           => q{SILENT},
        expected_output => q{--pedigreeValidationType SILENT},
    },
    read_filters_ref => {
        inputs_ref      => [qw{ MalformedRead BadCigar}],
        expected_output => q{--read_filter MalformedRead --read_filter BadCigar},
    },
    referencefile_path => {
        input           => catfile(qw{reference_dir human_genome_build.fasta }),
        expected_output => q{--reference_sequence }
          . catfile(qw{reference_dir human_genome_build.fasta }),
    },
    static_quantized_quals_ref => {
        inputs_ref =>
          [ $STATIC_QUANTIZED_QUALS_MINIMUM_LEVEL, $STATIC_QUANTIZED_QUALS_MAX_LEVEL ],
        expected_output => q{--static_quantized_quals }
          . $STATIC_QUANTIZED_QUALS_MINIMUM_LEVEL
          . q{ --static_quantized_quals }
          . $STATIC_QUANTIZED_QUALS_MAX_LEVEL,
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&gatk_base;

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
