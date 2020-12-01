#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile devnull };
use FindBin qw{ $Bin };
use List::Util qw{ any };
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
use MIP::Constants qw{ $COMMA $NEWLINE $SPACE };
use MIP::Test::Commands qw{ test_function };
use MIP::Test::Fixtures qw{ test_standard_cli };
use MIP::Test::Writefile qw{ test_write_to_file };

my $VERBOSE = 1;
our $VERSION = 1.02;

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
        q{MIP::Program::Gnu::Bash} => [qw{ gnu_set }],
        q{MIP::Test::Fixtures}     => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Gnu::Bash qw{ gnu_set };

diag(   q{Test gnu_set from Bash.pm v}
      . $MIP::Program::Gnu::Bash::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Specific arguments
my %argument = (
    set_errexit => {
        input           => 1,
        expected_output => q{set -e},
    },
    set_nounset => {
        input           => 1,
        expected_output => q{set -u},
    },
    set_pipefail => {
        input           => 1,
        expected_output => q{set -o pipefail},
    },
    unset_errexit => {
        input           => 1,
        expected_output => q{set +e},
    },
);

my @commands = gnu_set(
    {
        set_errexit   => $argument{set_errexit}{input},
        set_nounset   => $argument{set_nounset}{input},
        set_pipefail  => $argument{set_pipefail}{input},
        unset_errexit => $argument{unset_errexit}{input},
    }
);

## Testing return of commands
foreach my $key ( keys %argument ) {

    # Alias expected output
    my $expected_output = $argument{$key}{expected_output};

    ok( ( any { $_ eq $expected_output } @commands ), 'Argument: ' . $key );
}

## Testing write to file

## Base arguments
my @function_base_commands = qw{ set };

# Fake arguments
my @args = (
    filehandle  => undef,
    set_errexit => $argument{set_errexit}{input},
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&gnu_set;

test_write_to_file(
    {
        args_ref             => \@args,
        base_commands_ref    => \@function_base_commands,
        module_function_cref => $module_function_cref,
        separator            => $NEWLINE,
    }
);

done_testing();
