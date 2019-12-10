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

## Constants
Readonly my $HOM_REF_GENOTYPES_IN_CALL_SET => 7854;

BEGIN {
    use MIP::Test::Fixtures qw{ test_import };
    ### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Program::Gatk}  => [qw{ gatk_calculategenotypeposteriors }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Gatk qw{ gatk_calculategenotypeposteriors };

diag(   q{Test gatk_calculategenotypeposteriors from Variantcalling::Gatk.pm v}
      . $MIP::Program::Gatk::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ gatk CalculateGenotypePosteriors };

my %base_argument = (
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
    stderrfile_path => {
        input           => q{stderrfile.test},
        expected_output => q{2> stderrfile.test},
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    infile_path => {
        input           => catfile(qw{ my case.vcf  }),
        expected_output => q{--variant } . catfile(qw{ my case.vcf }),
    },
    outfile_path => {
        input           => catfile(qw{ my case_refined.vcf }),
        expected_output => q{--output } . catfile(qw{ my case_refined.vcf }),
    },
);

my %specific_argument = (
    infile_path => {
        input           => catfile(qw{ my case.vcf  }),
        expected_output => q{--variant } . catfile(qw{ my case.vcf }),
    },
    num_ref_samples_if_no_call => {
        input           => $HOM_REF_GENOTYPES_IN_CALL_SET,
        expected_output => q{--num-reference-samples-if-no-call }
          . $HOM_REF_GENOTYPES_IN_CALL_SET,
    },
    outfile_path => {
        input           => catfile(qw{ my case_refined.vcf }),
        expected_output => q{--output } . catfile(qw{ my case_refined.vcf }),
    },
    supporting_callset_file_path => {
        input           => catfile(qw{ my support.vcf }),
        expected_output => q{--supporting-callsets } . catfile(qw{ my support.vcf }),
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&gatk_calculategenotypeposteriors;

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
