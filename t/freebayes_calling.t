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
Readonly my $COMMA             => q{,};
Readonly my $SPACE             => q{ };
Readonly my $USE_BEST_N_ALLELS => q{4};

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Program::Variantcalling::Freebayes} => [qw{ freebayes_calling }],
        q{MIP::Test::Fixtures}                     => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Variantcalling::Freebayes qw{ freebayes_calling };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test freebayes_calling from Freebayes v}
      . $MIP::Program::Variantcalling::Freebayes::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ freebayes };

my %base_argument = (
    stderrfile_path => {
        input           => q{stderrfile.test},
        expected_output => q{2> stderrfile.test},
    },
    FILEHANDLE => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
    stdoutfile_path => {
        input           => q{stdoutfile_path.test},
        expected_output => q{1> stdoutfile_path.test},
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    infile_paths_ref => {
        inputs_ref      => [qw{ var_1.vcf var_2.vcf var_3.vcf }],
        expected_output => q{var_1.vcf var_2.vcf var_3.vcf},
    },
    referencefile_path => {
        input           => catfile(qw{reference_dir human_genome_build.fasta }),
        expected_output => q{--fasta-reference}
          . $SPACE
          . catfile(qw{reference_dir human_genome_build.fasta }),
    },
);

my %specific_argument = (
    apply_standard_filter => {
        input           => 1,
        expected_output => q{--standard-filters},
    },
    calculate_genotype_quality => {
        input           => 1,
        expected_output => q{--genotype-qualities},
    },
    use_best_n_alleles => {
        input           => $USE_BEST_N_ALLELS,
        expected_output => q{--use-best-n-alleles}
          . $SPACE
          . $USE_BEST_N_ALLELS,
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&freebayes_calling;

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

done_testing();
