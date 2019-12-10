#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };
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
        q{MIP::Program::Bcftools} => [qw{ bcftools_filter }],
        q{MIP::Test::Fixtures}    => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

## Constants
Readonly my $SNP_GAP_FILTER_DISTANCE   => 50;
Readonly my $INDEL_GAP_FILTER_DISTANCE => 100;

use MIP::Program::Bcftools qw{ bcftools_filter };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test bcftools_filter from Bcftools.pm v}
      . $MIP::Program::Bcftools::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ bcftools };

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
my %required_argument = ();

my %specific_argument = (
    exclude => {
        input           => q{%QUAL<10 || (RPB<0.1 && %QUAL<15)},
        expected_output => q{--exclude} . $SPACE . q{%QUAL<10 || (RPB<0.1 && %QUAL<15)},
    },
    filter_mode => {
        input           => q{+},
        expected_output => q{--mode} . $SPACE . q{+},
    },
    include => {
        input           => q{INFO/CSQ[*]~":p[.]"},
        expected_output => q{--include} . $SPACE . q{INFO/CSQ[*]~":p[.]"},
    },
    indel_gap => {
        input           => $INDEL_GAP_FILTER_DISTANCE,
        expected_output => q{--IndelGap} . $SPACE . $INDEL_GAP_FILTER_DISTANCE,
    },
    infile_path => {
        input           => q{infile.test},
        expected_output => q{infile.test},
    },
    snp_gap => {
        input           => $SNP_GAP_FILTER_DISTANCE,
        expected_output => q{--SnpGap} . $SPACE . $SNP_GAP_FILTER_DISTANCE,
    },
    soft_filter => {
        input           => q{LowQual},
        expected_output => q{--soft-filter} . $SPACE . q{LowQual},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&bcftools_filter;

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
