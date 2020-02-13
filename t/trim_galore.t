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
use MIP::Test::Fixtures qw{ test_standard_cli };

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
        q{MIP::Program::Trim_galore} => [qw{ trim_galore }],
        q{MIP::Test::Fixtures}       => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Trim_galore qw{ trim_galore };

diag(   q{Test trim_galore from Trim_galore.pm v}
      . $MIP::Program::Trim_galore::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

Readonly my $CORES => 12;

## Base arguments
my @function_base_commands = qw{ trim_galore };

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
        inputs_ref      => [qw{ file_1.fastq file_2.fastq }],
        expected_output => q{file_1.fastq} . $SPACE . q{file_2.fastq},
    },
);

my %specific_argument = (
    cores => {
        input           => $CORES,
        expected_output => q{--cores} . $SPACE . $CORES,
    },
    fastqc => {
        input           => 1,
        expected_output => q{--fastqc},
    },
    gzip_output => {
        input           => 1,
        expected_output => q{--gzip},
    },
    infile_paths_ref => {
        inputs_ref      => [qw{ file_1.fastq file_2.fastq }],
        expected_output => q{file_1.fastq} . $SPACE . q{file_2.fastq},
    },
    outdir_path => {
        input           => catdir(qw{ output dir }),
        expected_output => q{--output_dir} . $SPACE . catdir(qw{ output dir }),
    },
    paired_reads => {
        input           => 1,
        expected_output => q{--paired},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&trim_galore;

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
