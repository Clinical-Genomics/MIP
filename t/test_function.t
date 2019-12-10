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
use autodie qw { :all };
use Modern::Perl qw{ 2018 };
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
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
        q{MIP::Program::Bcftools} => [qw{ bcftools_view }],
        q{MIP::Test::Commands}    => [qw{ test_command test_function }],
        q{MIP::Test::Fixtures}    => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Bcftools qw{ bcftools_view };
use MIP::Test::Commands qw{ test_command test_function };

diag(   q{Test test_function from Commands.pm v}
      . $MIP::Test::Commands::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a test_command, when no expected output and a faulty base command

## Base arguments
my @function_base_commands = qw{ bcftools does_not_exists };

Readonly my $MAX_ALLELES => 2;
Readonly my $MIN_ALLELES => 2;

my %base_argument = (
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = ();

my %specific_argument = (
    apply_filters_ref => {
        inputs_ref      => [qw{ filter1 filter2 }],
        expected_output => q{not_the_expected_return},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&bcftools_view;

trap {
    test_function(
        {
            argument_href              => \%specific_argument,
            do_test_base_command       => 1,
            function_base_commands_ref => \@function_base_commands,
            module_function_cref       => $module_function_cref,
            required_argument_href     => \%required_argument,
            is_self_testing            => 1,
        }
    )
};
## Then throw error
like( $trap->stderr, qr/Command\sline\sdoes\snot/xms, q{Throw error message} );

## Given a scalar input when base command exists
@function_base_commands = qw{ test command };
%required_argument      = (
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
    hash_arg_href => {
        input_href => {
            water => q{wet},
            fire  => q{hot},
        },

        # Always sorted to an alphabetical order according to ASCII table
        expected_output => q{--hash_arg fire=hot --hash_arg water=wet},
    },
);
%specific_argument = (
    array_args_ref => {
        inputs_ref      => [qw{ test_value_1 test_value_2 }],
        expected_output => q{--array_args test_value_1 --array_args test_value_2},
    },
    hash_arg_href => {
        input_href => {
            water => q{wet},
            fire  => q{hot},
        },

        # Always sorted to an alphabetical order according to ASCII table
        expected_output => q{--hash_arg fire=hot --hash_arg water=wet},
    },
    scalar_arg => {
        input           => q{test_scalar},
        expected_output => q{--scalar_arg} . $SPACE . q{test_scalar},
    },
);

# Coderef - enables generalized use of generate call
$module_function_cref = \&test_command;

## Test both base and function specific arguments
my @arguments = ( \%required_argument, \%specific_argument );

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

## Given no filehandle and hash argument
@arguments = ( \%specific_argument );

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

done_testing();
