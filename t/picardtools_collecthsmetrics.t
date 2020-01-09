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
use autodie qw { :all };
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.03;

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
        q{MIP::Program::Picardtools} => [qw{ picardtools_collecthsmetrics }],
        q{MIP::Test::Fixtures}       => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Picardtools qw{ picardtools_collecthsmetrics };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test picardtools_collecthsmetrics from Picardtools.pm v}
      . $MIP::Program::Picardtools::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ picard CollectHsMetrics };

my %base_argument = (
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    bait_interval_file_paths_ref => {
        inputs_ref      => [ catfile(qw{ indirectory exome_padded_interval_list_1 }) ],
        expected_output => q{-BAIT_INTERVALS}
          . $SPACE
          . catfile(qw{ indirectory exome_padded_interval_list_1 }),
    },
    target_interval_file_paths_ref => {
        inputs_ref      => [ catfile(qw{ indirectory exome_interval_list_1 }) ],
        expected_output => q{-TARGET_INTERVALS}
          . $SPACE
          . catfile(qw{ indirectory exome_interval_list_1 }),
    },
    infile_path => {
        input           => catfile(qw{ indirectory infile_1 }),
        expected_output => q{-INPUT} . $SPACE . catfile(qw{ indirectory infile_1 }),
    },
    outfile_path => {
        input           => catfile(qw{ out_directory outfile }),
        expected_output => q{-OUTPUT} . $SPACE . catfile(qw{ out_directory outfile }),
    },
    referencefile_path => {
        input           => catfile(qw{ references grch37_homo_sapiens_-d5-.fasta }),
        expected_output => q{-R}
          . $SPACE
          . catfile(qw{ references grch37_homo_sapiens_-d5-.fasta }),
    },
);

my %specific_argument = (
    create_index => {
        input           => q{true},
        expected_output => q{-CREATE_INDEX true},
    },
    bait_interval_file_paths_ref => {
        inputs_ref      => [ catfile(qw{ indirectory exome_padded_interval_list_1 }) ],
        expected_output => q{-BAIT_INTERVALS}
          . $SPACE
          . catfile(qw{ indirectory exome_padded_interval_list_1 }),
    },
    target_interval_file_paths_ref => {
        inputs_ref      => [ catfile(qw{ indirectory exome_interval_list_1 }) ],
        expected_output => q{-TARGET_INTERVALS}
          . $SPACE
          . catfile(qw{ indirectory exome_interval_list_1 }),
    },
    infile_path => {
        input           => q{infile_1},
        expected_output => q{-INPUT infile_1},
    },
    outfile_path => {
        input           => catfile(qw{ out_directory outfile }),
        expected_output => q{-OUTPUT} . $SPACE . catfile(qw{ out_directory outfile }),
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&picardtools_collecthsmetrics;

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
