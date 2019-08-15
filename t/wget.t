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
use MIP::Constants qw{ $COMMA $EQUALS $SPACE };
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
        q{MIP::Program::Download::Wget} => [qw{ wget }],
        q{MIP::Test::Fixtures}          => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Download::Wget qw{ wget };

diag(   q{Test wget from Wget.pm v}
      . $MIP::Program::Download::Wget::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ wget };

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

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    FILEHANDLE => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
    url => {
        input           => q{https://www.gnu.org/software/wget},
        expected_output => q{https://www.gnu.org/software/wget},
    },
);

my %specific_argument = (
    continue => {
        input           => 1,
        expected_output => q{--continue},
    },
    outfile_path => {
        input           => catdir(qw{ outdir test }),
        expected_output => q{-O} . $SPACE . catdir(qw{ outdir test }),
    },
    quiet => {
        input           => 1,
        expected_output => q{--quiet},
    },
    read_timeout => {
        input           => 1,
        expected_output => q{--read-timeout} . $EQUALS . q{1},
    },
    retry_connrefused => {
        input           => 1,
        expected_output => q{--retry-connrefused},
    },
    timeout => {
        input           => 1,
        expected_output => q{--timeout} . $EQUALS . q{1},
    },
    tries => {
        input           => 1,
        expected_output => q{--tries} . $EQUALS . q{1},
    },
    url => {
        input           => q{https://www.gnu.org/software/wget},
        expected_output => q{https://www.gnu.org/software/wget},
    },
    user => {
        input           => q{superman},
        expected_output => q{--user} . $EQUALS . q{superman},
    },
    verbose => {
        input           => 1,
        expected_output => q{--verbose},
    },
    wait_retry => {
        input           => 1,
        expected_output => q{--waitretry} . $EQUALS . q{1},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&wget;

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
