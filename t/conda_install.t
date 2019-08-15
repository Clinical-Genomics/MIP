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
use MIP::Test::Commands qw{ test_function };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.05;

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
    my %perl_module = (
        q{MIP::Package_manager::Conda} => [qw{ conda_install }],
        q{MIP::Test::Fixtures}         => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Package_manager::Conda qw{conda_install};

diag(   q{Test conda_install from Conda.pm v}
      . $MIP::Package_manager::Conda::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ conda install };

my %base_argument = (
    FILEHANDLE => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

## Can be duplicated with %base and/or %specific to enable testing of each individual argument
my %required_argument = (
    packages_ref => {
        inputs_ref      => [qw{ test_package_1=1.2.3 test_package_2=1.2 }],
        expected_output => q{test_package_1=1.2.3 test_package_2=1.2},
    },
    FILEHANDLE => {
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
    no_update_dep => {
        input           => 1,
        expected_output => q{--no-update-deps},
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
my $module_function_cref = \&conda_install;

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
