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

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = ( q{MIP::Program::Svdb} => [qw{ svdb_merge }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Svdb qw{ svdb_merge };

diag(   q{Test svdb_merge from Svdb.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $BND_DISTANCE  => 2_000;
Readonly my $EVENT_OVERLAP => -1;

## Base arguments
my @function_base_commands = qw{ svdb --merge };

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
    infile_paths_ref => {
        inputs_ref      => [ catfile(qw{ a test infile_1 }), catfile(qw{ a test infile_2 }) ],
        expected_output => q{--vcf}
          . $SPACE
          . catfile(qw{ a test infile_1 })
          . $SPACE
          . catfile(qw{ a test infile_2 }),
    },
);

my %specific_argument = (
    bnd_distance => {
        input           => $BND_DISTANCE,
        expected_output => q{--bnd_distance} . $SPACE . $BND_DISTANCE,
    },
    infile_paths_ref => {
        inputs_ref      => [ catfile(qw{ a test infile_1 }), catfile(qw{ a test infile_2 }) ],
        expected_output => q{--vcf}
          . $SPACE
          . catfile(qw{ a test infile_1 })
          . $SPACE
          . catfile(qw{ a test infile_2 }),
    },
    notag => {
        input           => q{1},
        expected_output => q{--notag},
    },
    overlap => {
        input           => $EVENT_OVERLAP,
        expected_output => q{--overlap} . $SPACE . $EVENT_OVERLAP,
    },
    pass_only => {
        input           => 1,
        expected_output => q{--pass_only},
    },
    priority => {
        input           => q{manta,delly,cnvnator,tiddit},
        expected_output => q{--priority} . $SPACE . q{manta,delly,cnvnator,tiddit},
    },
    same_order => {
        input           => q{1},
        expected_output => q{--same_order},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&svdb_merge;

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
