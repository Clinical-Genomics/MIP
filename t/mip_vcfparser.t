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
our $VERSION = 1.04;

## Constants
Readonly my $PADDING => 50;

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
        q{MIP::Program::Mip}   => [qw{ mip_vcfparser }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Mip qw{ mip_vcfparser };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test mip_vcfparser from Mip.pm v}
      . $MIP::Program::Mip::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ mip vcfparser };

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
    infile_path => {
        input => catfile(
            qw{ file_path_prefix_contig_analysis-type_suffix.annotation_infile_number}),
        expected_output => catfile(
            qw{ file_path_prefix_contig_analysis-type_suffix.annotation_infile_number}),
    },
);

my %specific_argument = (
    log_file_path => {
        input           => catfile(qw{ a dir vcfparser_contig.log}),
        expected_output => q{--log_file}
          . $SPACE
          . catfile(qw{ a dir vcfparser_contig.log}),
    },
    padding => {
        input           => $PADDING,
        expected_output => q{--padding} . $SPACE . $PADDING,
    },
    parse_vep => {
        input           => 1,
        expected_output => q{--parse_vep},
    },
    per_gene => {
        input           => 1,
        expected_output => q{--per_gene},
    },
    pli_values_file_path => {
        input           => catfile(qw{a dir plifile_path}),
        expected_output => q{--pli_values_file } . catfile(qw{a dir plifile_path}),
    },
    range_feature_annotation_columns_ref => {
        inputs_ref => [qw{ feature_anno1 feature_anno2 }],
        expected_output =>
          q{--range_feature_annotation_columns feature_anno1,feature_anno2},
    },
    range_feature_file_path => {
        input           => q{sv_vcfparser_range_feature_file},
        expected_output => q{--range_feature_file sv_vcfparser_range_feature_file},
    },
    select_feature_annotation_columns_ref => {
        inputs_ref => [qw{ feature_anno1 feature_anno2 }],
        expected_output =>
          q{--select_feature_annotation_columns feature_anno1,feature_anno2},
    },
    select_feature_file_path => {
        input           => catfile(qw{ active_parameter_href vcfparser_select_file }),
        expected_output => q{--select_feature_file}
          . $SPACE
          . catfile(qw{ active_parameter_href vcfparser_select_file }),
    },
    select_feature_matching_column => {
        input           => 2,
        expected_output => q{--select_feature_matching_column 2},
    },
    select_outfile => {
        input           => catdir(qw{ outfile_path prefix_contig.selectedsuffix }),
        expected_output => q{--select_outfile}
          . $SPACE
          . catdir(qw{ outfile_path prefix_contig.selectedsuffix }),
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&mip_vcfparser;

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
