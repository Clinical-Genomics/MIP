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
        q{MIP::Program::Plink} => [qw{ plink_variant_pruning }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

## Constants
Readonly my $INDEP_WINDOW_SIZE   => 50;
Readonly my $INDEP_STEP_SIZE     => 5;
Readonly my $INDEP_VIF_THRESHOLD => 2;

use MIP::Program::Plink qw{ plink_variant_pruning };

diag(   q{Test plink_variant_pruning from MIP::Program::Plink v}
      . $MIP::Program::Plink::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ plink2 };

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
    const_fid => {
        input           => q{case_id},
        expected_output => q{--const-fid} . $SPACE . q{case_id},
    },
    indep => {
        input           => 1,
        expected_output => q{--indep},
    },
    indep_step_size => {
        input           => $INDEP_STEP_SIZE,
        expected_output => $INDEP_STEP_SIZE,
    },
    indep_vif_threshold => {
        input           => $INDEP_VIF_THRESHOLD,
        expected_output => $INDEP_VIF_THRESHOLD,
    },
    indep_window_size => {
        input           => $INDEP_WINDOW_SIZE,
        expected_output => $INDEP_WINDOW_SIZE,
    },

    outfile_prefix => {
        input           => catfile(qw{ temp_directory $case_id _data }),
        expected_output => q{--out}
          . $SPACE
          . catfile(qw{ temp_directory $case_id _data }),
    },
    set_missing_var_ids => {
        input           => q?@:#[hg19]\$1,\$2?,
        expected_output => q{--set-missing-var-ids} . $SPACE . q?@:#[hg19]\$1,\$2?,
    },
    vcffile_path => {
        input           => catfile(qw{ dir infile.vcf }),
        expected_output => q{--vcf} . $SPACE . catfile(qw{ dir infile.vcf }),
    },
);

my %specific_argument = (
    make_bed => {
        input           => 1,
        expected_output => q{--make-bed},
    },
    vcf_half_call => {
        input           => q{haploid},
        expected_output => q{--vcf-half-call haploid},
    },
    vcf_require_gt => {
        input           => 1,
        expected_output => q{--vcf-require-gt},
    },
);

# Coderef - enables generalized use of generate call
my $module_function_cref = \&plink_variant_pruning;

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
