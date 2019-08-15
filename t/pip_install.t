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
our $VERSION = 1.02;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $SPACE   => q{ };
Readonly my $NEWLINE => qq{\n};
Readonly my $COMMA   => q{,};

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Package_manager::Pip} => [qw{ pip_install }],
        q{MIP::Test::Fixtures}       => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Package_manager::Pip qw{ pip_install };

diag(   q{Test pip_install from Pip.pm v}
      . $MIP::Package_manager::Pip::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ pip install };

my %base_argument = (
    FILEHANDLE => {
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

my %specific_argument = (
    editable => {
        input           => catdir(qw{ test path }),
        expected_output => q{--editable test/path},
    },
    FILEHANDLE => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
    packages_ref => {
        inputs_ref      => [qw{ test_package_1 test_package_2 }],
        expected_output => q{test_package_1 test_package_2},
    },
    no_cache_dir => {
        input           => 1,
        expected_output => q{--no-cache-dir},
    },
    quiet => {
        input           => 1,
        expected_output => q{--quiet},
    },
    requirement => {
        input           => q{test_file.txt},
        expected_output => q{--requirement test_file.txt},
    },
    verbose => {
        input           => 1,
        expected_output => q{--verbose},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&pip_install;

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
        }
    );
}

## When given the python_module option
my @commands = pip_install( { python_module => 1, } );
## Then prepend python -m to base command
is( $commands[0], q{python -m}, q{Argument: python_module} );

done_testing();
