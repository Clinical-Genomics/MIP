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
        q{MIP::Program::Dragen} => [qw{ dragen_build_hash_table }],
        q{MIP::Test::Fixtures}  => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Dragen qw{ dragen_build_hash_table };

diag(   q{Test dragen_build_hash_table from Dragen.pm v}
      . $MIP::Program::Dragen::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $THREAD_NUMBER => 2;

## Base arguments
my @function_base_commands = qw{ dragen };

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
    outdirectory_path => {
        input           => catdir(qw{ outdirectory path }),
        expected_output => q{--output-directory}
          . $SPACE
          . catdir(qw{ outdirectory path }),
    },
    reference_genome_file_path => {
        input           => catfile(qw{a reference_file.fasta}),
        expected_output => q{--ht-reference}
          . $SPACE
          . catfile(qw{a reference_file.fasta}),
    },
);

my %specific_argument = (
    build_hash_table => {
        input           => 1,
        expected_output => q{--build-hash-table} . $SPACE . q{true},
    },
    enable_cnv => {
        input           => 1,
        expected_output => q{--enable-cnv} . $SPACE . q{true},
    },
    ht_alt_liftover_file_path => {
        input           => catfile(qw{a reference_liftover_file.fasta}),
        expected_output => q{--ht-alt-liftover}
          . $SPACE
          . catfile(qw{a reference_liftover_file.fasta}),
    },
    ht_decoys_file_path => {
        input           => catfile(qw{a decoy_file.fasta}),
        expected_output => q{--ht-decoys} . $SPACE . catfile(qw{a decoy_file.fasta}),
    },
    outdirectory_path => {
        input           => catdir(qw{ outdirectory path }),
        expected_output => q{--output-directory}
          . $SPACE
          . catdir(qw{ outdirectory path }),
    },
    reference_genome_file_path => {
        input           => catfile(qw{a reference_file.fasta}),
        expected_output => q{--ht-reference}
          . $SPACE
          . catfile(qw{a reference_file.fasta}),
    },
    thread_number => {
        input           => $THREAD_NUMBER,
        expected_output => q{--ht-num-threads} . $SPACE . $THREAD_NUMBER,
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&dragen_build_hash_table;

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
