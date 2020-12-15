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
use MIP::Constants qw{ $COMMA $DOUBLE_QUOTE $SPACE $TAB };
use MIP::Test::Commands qw{ test_function };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Program::Bwa}   => [qw{ run_bwamem }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Bwa qw{ run_bwamem };

diag(   q{Test run_bwamem from Bwa.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ run-bwamem };

## Read group header line
my @read_group_headers = (
    $DOUBLE_QUOTE . q{@RG},
    q{ID:} . q{1_140128_H8AHNADXX_1-1-2A_GATCAG_1} . $TAB,
    q{SM:} . q{1-1-2A},
    q{PL:} . q{ILLUMINA} . $DOUBLE_QUOTE,
);

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
);

## Can be duplicated with %base and/or %specific to enable testing of each individual argument
my %required_argument = (
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
    idxbase => {
        input           => q{grch37_homo_sapiens_-d5-.fasta},
        expected_output => q{grch37_homo_sapiens_-d5-.fasta},
    },
    infile_path => {
        input           => q{test_infile.fastq},
        expected_output => q{test_infile.fastq},
    },
    outfiles_prefix_path => {
        input           => q{test_outfile.bam},
        expected_output => q{-o test_outfile.bam},
    },
);

my %specific_argument = (
    hla_typing => {
        input           => 1,
        expected_output => q{-H},
    },
    infile_path => {
        input           => q{test_infile_1.fastq},
        expected_output => q{test_infile_1.fastq},
    },
    outfiles_prefix_path => {
        input           => q{test_outfile.bam},
        expected_output => q{-o test_outfile.bam},
    },
    read_group_header => {
        input           => ( join $SPACE,                  @read_group_headers ),
        expected_output => ( q{-R} . $SPACE . join $SPACE, @read_group_headers ),
    },
    second_infile_path => {
        input           => q{test_infile_2.fastq},
        expected_output => q{test_infile_2.fastq},
    },
    stderrfile_path => {
        input           => q{stderrfile.test},
        expected_output => q{2> stderrfile.test},
    },
    stderrfile_path_append => {
        input           => q{stderrfile.test},
        expected_output => q{2>> stderrfile.test},
    },
    thread_number => {
        input           => 2,
        expected_output => q{-t 2},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&run_bwamem;

## Test both base and function specific arguments
my @arguments = ( \%required_argument, \%specific_argument );

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
