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
our $VERSION = 1.00;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COMMA        => q{,};
Readonly my $DOUBLE_QUOTE => q{"};
Readonly my $SPACE        => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Gnu::Findutils} => [qw{ xargs }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Gnu::Findutils qw{xargs};
use MIP::Test::Commands qw{test_function};

diag(   q{Test xargs from Findutils v}
      . $MIP::Gnu::Findutils::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $MAX_PROCESSES => 16;
Readonly my $MAX_ARGUMENTS => 1;

## Base arguments
my @function_base_commands = qw{ xargs };

my %base_argument = (
    FILEHANDLE => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    shell_commands_ref => {
        inputs_ref      => [qw{java -jar}],
        expected_output => q{java -jar},
    },
    FILEHANDLE => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

my %specific_argument = (
    max_args => {
        input           => $MAX_ARGUMENTS,
        expected_output => q{-n} . $SPACE . $MAX_ARGUMENTS,
    },
    max_procs => {
        input           => $MAX_PROCESSES,
        expected_output => q{-P} . $SPACE . $MAX_PROCESSES,
    },
    null_character => {
        input           => 1,
        expected_output => q{-0},
    },
    placeholder_symbol => {
        input           => q?{}?,
        expected_output => q?{}? . $DOUBLE_QUOTE . $SPACE,
    },
    replace_str => {
        input           => 1,
        expected_output => q{-i},
    },
    shell_commands_ref => {
        inputs_ref      => [qw{java -jar}],
        expected_output => q{java -jar},
    },
    verbose => {
        input           => 1,
        expected_output => q{--verbose},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&xargs;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

HASHES_OF_ARGUMENTS:
foreach my $argument_href (@arguments) {
    my @commands = test_function(
        {
            argument_href              => $argument_href,
            required_argument_href     => \%required_argument,
            module_function_cref       => $module_function_cref,
            function_base_commands_ref => \@function_base_commands,
            do_test_base_command       => 1,
        }
    );
}

done_testing();
