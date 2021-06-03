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
    my %perl_module = (
        q{MIP::Program::Picardtools} => [qw{ picardtools_collectrnaseqmetrics }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Picardtools qw{ picardtools_collectrnaseqmetrics };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test picardtools_collectrnaseqmetrics from Picardtools.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ picard CollectRnaSeqMetrics };

my %base_argument = (
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    gene_annotation_file_path => {
        input =>
          catfile(qw{ references grch37_gencode_annotation_reformated-v31-.refFlat }),
        expected_output => q{-REF_FLAT}
          . $SPACE
          . catfile(qw{ references grch37_gencode_annotation_reformated-v31-.refFlat }),
    },
    infile_path => {
        input           => catfile(qw{ indirectory infile_1 }),
        expected_output => q{-INPUT} . $SPACE . catfile(qw{ indirectory infile_1 }),
    },
    outfile_path => {
        input           => catfile(qw{ out_directory outfile }),
        expected_output => q{-OUTPUT} . $SPACE . catdir(qw{ out_directory outfile }),
    },
    strand_specificity => {
        input           => q{NONE},
        expected_output => q{-STRAND_SPECIFICITY NONE},
    },
);

my %specific_argument = (
    chart_outfile_path => {
        input           => catfile(qw{ out chart.pdf}),
        expected_output => q{-CHART_OUTPUT} . $SPACE . catfile(qw{ out chart.pdf}),
    },
    gene_annotation_file_path => {
        input =>
          catfile(qw{ references grch37_gencode_annotation_reformated-v31-.refFlat }),
        expected_output => q{-REF_FLAT}
          . $SPACE
          . catfile(qw{ references grch37_gencode_annotation_reformated-v31-.refFlat }),
    },
    infile_path => {
        input           => catfile(qw{ indirectory infile_1 }),
        expected_output => q{-INPUT} . $SPACE . catfile(qw{ indirectory infile_1 }),
    },
    outfile_path => {
        input           => catfile(qw{ out_directory outfile }),
        expected_output => q{-OUTPUT} . $SPACE . catfile(qw{ out_directory outfile }),
    },
    rrna_intervals_file_path => {
        input           => catfile(qw{ references rrna_intervals.interval_list }),
        expected_output => q{-RIBOSOMAL_INTERVALS}
          . $SPACE
          . catfile(qw{  references rrna_intervals.interval_list }),
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&picardtools_collectrnaseqmetrics;

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
