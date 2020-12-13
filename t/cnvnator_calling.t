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
    my %perl_module = (
        q{MIP::Program::Cnvnator} => [qw{ cnvnator_calling }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

## Constants
Readonly my $BIN_SIZE => 1000;

use MIP::Program::Cnvnator qw{ cnvnator_calling };

diag(   q{Test cnvnator_calling from Cnvnator}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ cnvnator };

my %base_argument = (
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
    stderrfile_path => {
        input           => q{stderrfile.test},
        expected_output => q{2> stderrfile.test},
    },
    stdoutfile_path => {
        input           => q{stdoutfile_path.test},
        expected_output => q{1> stdoutfile_path.test},
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    infile_path => {
        input           => q{infile_path},
        expected_output => q{infile_path},
    },
);

my %specific_argument = (
    cnv_bin_size => {
        input           => $BIN_SIZE,
        expected_output => q{-call} . $SPACE . $BIN_SIZE,
    },
    regions_ref => {
        inputs_ref      => [qw{ 1 2 3 }],
        expected_output => q{-chrom 1 2 3},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&cnvnator_calling;

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
