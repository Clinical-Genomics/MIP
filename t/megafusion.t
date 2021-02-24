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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Commands qw{ test_function };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = ( q{MIP::Program::Megafusion} => [qw{ megafusion }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Megafusion qw{ megafusion };

diag(   q{Test megafusion from Megafusion.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ MegaFusion.py };

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

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    config_file_path => {
        input           => catfile(qw{ path to arriba.json }),
        expected_output => q{--json} . $SPACE . catfile(qw{ path to arriba.json }),
    },
    infile_path => {
        input           => catfile(qw{ path to arriba_fusions.tsv }),
        expected_output => q{--fusion} . $SPACE . catfile(qw{ path to arriba_fusions.tsv }),
    },
);

my %specific_argument = (
    config_file_path => {
        input           => catfile(qw{ path to arriba.json }),
        expected_output => q{--json} . $SPACE . catfile(qw{ path to arriba.json }),
    },
    infile_path => {
        input           => catfile(qw{ path to arriba_fusions.tsv }),
        expected_output => q{--fusion} . $SPACE . catfile(qw{ path to arriba_fusions.tsv }),
    },
    sample_id => {
        input           => q{sample_id},
        expected_output => q{--sample} . $SPACE . q{sample_id},
    },
);

my $module_function_cref = \&megafusion;

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
