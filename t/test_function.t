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
use Modern::Perl qw{ 2014 };
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

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
        q{MIP::Program::Variantcalling::Bcftools} => [qw{ bcftools_view }],
        q{MIP::Test::Commands}                    => [qw{ test_function }],
        q{MIP::Test::Fixtures}                    => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Variantcalling::Bcftools qw{ bcftools_view };
use MIP::Test::Commands qw{ test_function };

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
    FILEHANDLE => {
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
say STDERR $trap->stderr;
## Then throw error
like( $trap->stderr, qr/Command \s+ line \s+ does \s+ not/xms, q{Throw error message} );

done_testing();
