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
        q{MIP::Program::Conda} => [qw{ conda_create }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Conda qw{conda_create};
use MIP::Test::Commands qw{test_function};

diag(   q{Test conda_create }
      . $MIP::Program::Conda::VERSION
      . q{, Perl }
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ conda create };

my %base_argument = (
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

## Can be duplicated with %base and/or %specific to enable testing of each individual argument
my %required_argument = (
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

my %specific_argument = (
    conda_channels_ref => {
        inputs_ref      => [qw{ bioconda conda-forge }],
        expected_output => q{--channel bioconda --channel conda-forge},
    },
    env_name => {
        input           => q{test_env},
        expected_output => q{--name test_env},
    },
    no_confirmation => {
        input           => 1,
        expected_output => q{--yes},
    },
    packages_ref => {
        inputs_ref      => [qw{ test_package_1 test_package_2}],
        expected_output => q{test_package_1 test_package_2},
    },
    quiet => {
        input           => 1,
        expected_output => q{--quiet},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&conda_create;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

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
