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
our $VERSION = 1.04;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $GAUSSIANS    => 8;
Readonly my $MAX_ATTEMPTS => 2;

BEGIN {
    use MIP::Test::Fixtures qw{ test_import };
    ### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Program::Gatk}  => [qw{ gatk_variantrecalibrator }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Gatk qw{ gatk_variantrecalibrator };

diag(   q{Test gatk_variantrecalibrator from Gatk.pm v}
      . $MIP::Program::Gatk::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ gatk VariantRecalibrator };

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
    annotations_ref => {
        inputs_ref      => [qw{ MQRankSum QD }],
        expected_output => q{--use-annotation MQRankSum --use-annotation: QD},
    },
    infile_path => {
        input           => catfile(qw{ my case.vcf  }),
        expected_output => q{--variant } . catfile(qw{ my case.vcf }),
    },
    outfile_path => {
        input           => catfile(qw{ my output.recal }),
        expected_output => q{--output } . catfile(qw{ my output.recal }),
    },
    resources_ref => {
        inputs_ref      => [qw{ resource_1 resource_2 }],
        expected_output => q{--resource:resource_1 --resource:resource_2},
    },
    tranches_file_path => {
        input           => catfile(qw{ my output.tranches }),
        expected_output => q{--tranches-file } . catdir(qw{ my output.tranches }),
    },
);

my %specific_argument = (
    infile_path => {
        input           => catfile(qw{ my case.vcf }),
        expected_output => q{--variant } . catfile(qw{ my case.vcf }),
    },
    max_attempts => {
        input           => $MAX_ATTEMPTS,
        expected_output => q{--max-attempts } . $MAX_ATTEMPTS,
    },
    max_gaussian_level => {
        input           => $GAUSSIANS,
        expected_output => q{--max-gaussians } . $GAUSSIANS,
    },
    mode => {
        input           => q{SNP},
        expected_output => q{--mode SNP},
    },
    outfile_path => {
        input           => catfile(qw{ my output.recal }),
        expected_output => q{--output } . catfile(qw{ my output.recal }),
    },
    rscript_file_path => {
        input           => catfile(qw{ my plots.R }),
        expected_output => q{--rscript-file } . catfile(qw{ my plots.R }),
    },
    resources_ref => {
        inputs_ref      => [qw{ resource_1 resource_2 }],
        expected_output => q{--resource:resource_1 --resource:resource_2},
    },
    tranches_file_path => {
        input           => catfile(qw{ my output.tranches }),
        expected_output => q{--tranches-file } . catdir(qw{ my output.tranches }),
    },
    ts_tranches_ref => {
        inputs_ref      => [qw{ 100 99.9 }],
        expected_output => q{-tranche 100 -tranche 99.9},
    },
    trust_all_polymorphic => {
        input           => 1,
        expected_output => q{--trust-all-polymorphic},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&gatk_variantrecalibrator;

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
