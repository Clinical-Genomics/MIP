#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname  };
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
Readonly my $COMMA             => q{,};
Readonly my $EXPECT_VALUE      => 1e-2;
Readonly my $MAX_TARGET_SEQS   => 1001;
Readonly my $MAX_THREAD_NUMBER => 16;
Readonly my $SPACE             => q{ };
Readonly my $TABULAR           => 6;
Readonly my $WORD_SIZE         => 11;

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Program::Alignment::Blast} => [qw{ blast_blastn }],
        q{MIP::Test::Fixtures}            => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Alignment::Blast qw{ blast_blastn };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test blast_blastn from Blast.pm v}
      . $MIP::Program::Alignment::Blast::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ blastn };

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
    database_name => {
        input           => catfile(qw{ a test database.fa }),
        expected_output => q{-db} . $SPACE . catfile(qw{ a test database.fa }),
    },
    query_file_path => {
        input           => catfile(qw{ a test query_file.fa }),
        expected_output => q{-query} . $SPACE . catfile(qw{ a test query_file.fa }),
    },
);

my %specific_argument = (
    database_name => {
        input           => catfile(qw{ a test database.fa }),
        expected_output => q{-db} . $SPACE . catfile(qw{ a test database.fa }),
    },
    evalue => {
        input           => $EXPECT_VALUE,
        expected_output => q{-evalue} . $SPACE . $EXPECT_VALUE,
    },
    lcase_masking => {
        input           => 1,
        expected_output => q{-lcase_masking},
    },
    max_target_seqs => {
        input           => $MAX_TARGET_SEQS,
        expected_output => q{-max_target_seqs} . $SPACE . $MAX_TARGET_SEQS,
    },
    output_format => {
        input           => $TABULAR,
        expected_output => q{-outfmt} . $SPACE . $TABULAR,
    },
    query_file_path => {
        input           => catfile(qw{ a test query_file.fa }),
        expected_output => q{-query} . $SPACE . catfile(qw{ a test query_file.fa }),
    },
    thread_number => {
        input           => $MAX_THREAD_NUMBER,
        expected_output => q{-num_threads} . $SPACE . $MAX_THREAD_NUMBER,
    },
    word_size => {
        input           => $WORD_SIZE,
        expected_output => q{-word_size} . $SPACE . $WORD_SIZE,
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&blast_blastn;

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
