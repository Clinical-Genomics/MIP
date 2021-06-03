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
    my %perl_module = ( q{MIP::Program::Pigz} => [qw{ pigz }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Pigz qw{ pigz };

diag(   q{Test pigz from Pigz}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ pigz };

my %base_argument = (
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

my %specific_argument = (
    decompress => {
        input           => q{decompress},
        expected_output => q{--decompress},
    },
    infile_path => {
        input           => q{infile_path},
        expected_output => q{infile_path},
    },
    outfile_path => {
        input           => q{outfile_path},
        expected_output => q{> outfile_path},
    },
    processes => {
        input           => 2,
        expected_output => q{--processes 2},
    },
    quiet => {
        input           => 1,
        expected_output => q{--quiet},
    },
    stderrfile_path => {
        input           => q{stderrfile.test},
        expected_output => q{2> stderrfile.test},
    },
    stderrfile_path_append => {
        input           => q{stderrfile_path_append},
        expected_output => q{2>> stderrfile_path_append},
    },
    stdout => {
        input           => q{stdout},
        expected_output => q{--stdout},
    },
    verbose => {
        input           => 1,
        expected_output => q{--verbose},
    },

);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&pigz;

my @commands = test_function(
    {
        argument_href              => \%specific_argument,
        do_test_base_command       => 1,
        function_base_commands_ref => \@function_base_commands,
        module_function_cref       => $module_function_cref,
    }
);

done_testing();
