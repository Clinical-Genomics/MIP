#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname  };
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
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.03;

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
    my %perl_module = ( q{MIP::Test::Fixtures} => [qw{ test_standard_cli }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Variantcalling::Gatk qw{ gatk_selectvariants };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test gatk_selectvariants from Variantcalling::Gatk.pm v}
      . $MIP::Program::Variantcalling::Gatk::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ gatk SelectVariants };

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
        input           => catfile(qw{ my variants.vcf }),
        expected_output => q{--variant } . catfile(qw{ my variants.vcf }),
    },
    outfile_path => {
        input           => catfile(qw{ my outfile }),
        expected_output => q{--output } . catfile(qw{ my outfile }),
    },
);

my %specific_argument = (
    exclude_non_variants => {
        input           => 1,
        expected_output => q{--exclude-non-variants},
    },
    infile_path => {
        input           => catfile(qw{ my variants.vcf }),
        expected_output => q{--variant } . catfile(qw{ my variants.vcf }),
    },
    outfile_path => {
        input           => catfile(qw{ my outfile }),
        expected_output => q{--output } . catfile(qw{ my outfile }),
    },
    restrict_alleles_to => {
        input           => q{BIALLELIC},
        expected_output => q{--restrict-alleles-to BIALLELIC},
    },
    sample_names_ref => {
        inputs_ref      => [qw{ mother child }],
        expected_output => q{--sample-name mother --sample-name child},
    },
    select_type_to_include_ref => {
        inputs_ref      => [qw{ INDEL SNP }],
        expected_output => q{--select-type-to-include INDEL --select-type-to-include SNP},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&gatk_selectvariants;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

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
