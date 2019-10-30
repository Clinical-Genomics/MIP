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
our $VERSION = 1.03;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $VEP_VERSION => 91;

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Program::Variantcalling::Vep} => [qw{ variant_effect_predictor_install }],
        q{MIP::Test::Fixtures}               => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Variantcalling::Vep qw{ variant_effect_predictor_install };

diag(   q{Test variant_effect_predictor_install from Vep.pm v}
      . $MIP::Program::Variantcalling::Vep::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ INSTALL.pl };

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
my %required_argument = (
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

my %specific_argument = (
    assembly => {
        input           => q{GRCh37},
        expected_output => q{--ASSEMBLY GRCh37},
    },
    auto => {
        input           => q{alcf},
        expected_output => q{--AUTO alcf},
    },
    cache_directory => {
        input           => catdir(qw{ ensembl cache }),
        expected_output => q{--CACHEDIR } . catdir(qw{ ensembl cache }),
    },
    cache_version => {
        input           => $VEP_VERSION,
        expected_output => q{--CACHE_VERSION} . $SPACE . $VEP_VERSION,
    },
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
    no_update => {
        input           => 1,
        expected_output => q{--NO_UPDATE},
    },
    no_htslib => {
        input           => 1,
        expected_output => q{--NO_HTSLIB},
    },
    plugins_ref => {
        inputs_ref      => [qw{ LofTool maxEntScan}],
        expected_output => q{--PLUGINS LofTool,maxEntScan},
    },
    species_ref => {
        inputs_ref      => [qw{ homo_spaiens }],
        expected_output => q{--SPECIES homo_spaiens},
    },
    version => {
        input           => $VEP_VERSION,
        expected_output => q{--VERSION} . $SPACE . $VEP_VERSION,
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&variant_effect_predictor_install;

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
