#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
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
        q{MIP::Program::Bedtools} => [qw{ bedtools_makewindows }],
        q{MIP::Test::Fixtures}    => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Bedtools qw{ bedtools_makewindows };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test bedtools_makewindows from Bedtools.pm v}
      . $MIP::Program::Bedtools::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

Readonly my $WINDOW_SIZE => 200_000;
Readonly my $STEP_SIZE   => 199_750;

## Base arguments
my @function_base_commands = qw{ bedtools makewindows };

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

## Can be duplicated with %base and/or %specific to enable testing of each individual argument
my %required_argument = (
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
    infile_bed_path => {
        input           => catfile(qw{ path to infile_bed }),
        expected_output => q{-b} . $SPACE . catfile(qw{ path to infile_bed }),
    },
    window_size => {
        input           => $WINDOW_SIZE,
        expected_output => q{-w} . $SPACE . $WINDOW_SIZE,
    },
);

## Specific arguments
my %specific_argument = (
    infile_bed_path => {
        input           => catfile(qw{ path to infile_bed }),
        expected_output => q{-b} . $SPACE . catfile(qw{ path to infile_bed }),
    },
    step_size => {
        input           => $STEP_SIZE,
        expected_output => q{-s} . $SPACE . $STEP_SIZE,
    },
    window_size => {
        input           => $WINDOW_SIZE,
        expected_output => q{-w} . $SPACE . $WINDOW_SIZE,
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&bedtools_makewindows;

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
