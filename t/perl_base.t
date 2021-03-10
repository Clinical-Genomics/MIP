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

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = ( q{MIP::Language::Perl} => [qw{ perl_base }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Language::Perl qw{ perl_base };

diag(   q{Test perl_base from Perl.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ perl };

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument

my %specific_argument = (
    autosplit => {
        input           => 1,
        expected_output => q{-a},
    },
    command_line => {
        input           => 1,
        expected_output => q{-e},
    },
    inplace => {
        input           => 1,
        expected_output => q{-i},
    },
    n => {
        input           => 1,
        expected_output => q{-n},
    },
    print => {
        input           => 1,
        expected_output => q{-p},
    },
    print_newline => {
        input           => 1,
        expected_output => q{-l},
    },
    use_container => => {
        input           => 1,
        expected_output => q{perl},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&perl_base;

## Test both base and function specific arguments
my @arguments = ( \%specific_argument );

ARGUMENT_HASH_REF:
foreach my $argument_href (@arguments) {

    test_function(
        {
            argument_href              => $argument_href,
            do_test_base_command       => 1,
            function_base_commands_ref => \@function_base_commands,
            module_function_cref       => $module_function_cref,
        }
    );
}

done_testing();
