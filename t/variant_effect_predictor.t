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
our $VERSION = 1.01;

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
        q{MIP::Program::Variantcalling::Vep} => [qw{ variant_effect_predictor }],
        q{MIP::Test::Fixtures}               => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Variantcalling::Vep qw{ variant_effect_predictor };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test variant_effect_predictor from Vep.pm v}
      . $MIP::Program::Variantcalling::Vep::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $VARIANT_BUFFERT_SIZE => 20_000;

## Base arguments
my @function_base_commands = qw{ vep };

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
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

my %specific_argument = (
    assembly => {
        input           => q{GRCh37},
        expected_output => q{--assembly} . $SPACE . q{GRCh37},
    },
    buffer_size => {
        input           => $VARIANT_BUFFERT_SIZE,
        expected_output => q{--buffer_size} . $SPACE . $VARIANT_BUFFERT_SIZE,
    },
    cache_directory => {
        input           => catdir( q{test_dir}, q{test_cache_dir} ),
        expected_output => q{--dir_cache}
          . $SPACE
          . catdir( q{test_dir}, q{test_cache_dir} ),
    },
    custom_annotations_ref => {
        inputs_ref => [
            (
                q{path,key,file_type,annotation_type,force_report_coordinates},
                q{path_1key_1,file_type_1,annotation_type_1,force_report_coordinates_1}
            )
        ],
        expected_output =>
          q{--custom path,key,file_type,annotation_type,force_report_coordinates}
          . $SPACE
          . q{--custom path_1key_1,file_type_1,annotation_type_1,force_report_coordinates_1},
    },
    distance => {
        input           => 10,
        expected_output => q{--distance} . $SPACE . q{10},
    },
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
    fork => {
        input           => 1,
        expected_output => q{--fork} . $SPACE . q{1},
    },
    infile_format => {
        input           => q{vcf},
        expected_output => q{--format} . $SPACE . q{vcf},
    },
    infile_path => {
        input           => catfile( q{test_dir}, q{infile.vcf} ),
        expected_output => q{--input_file}
          . $SPACE
          . catfile( q{test_dir}, q{infile.vcf} ),
    },
    max_sv_size => {
        input           => 1,
        expected_output => q{--max_sv_size} . $SPACE . 1,
    },
    outfile_format => {
        input           => q{vcf},
        expected_output => q{--} . q{vcf},
    },
    outfile_path => {
        input           => catfile( q{test_dir}, q{infile.vcf} ),
        expected_output => q{--output_file}
          . $SPACE
          . catfile( q{test_dir}, q{infile.vcf} ),
    },
    plugins_dir_path => {
        input           => catdir(qw{ test_dir plugins }),
        expected_output => q{--dir_plugins} . $SPACE . catdir(qw{ test_dir plugins }),
    },
    plugins_ref => {
        inputs_ref      => [qw{ LoFtool LoF }],
        expected_output => q{--plugin LoFtool} . $SPACE . q{--plugin LoF},
    },
    reference_path => {
        input           => catfile( q{test_dir}, q{hum_ref.pl} ),
        expected_output => q{--fasta} . $SPACE . catfile( q{test_dir}, q{hum_ref.pl} ),
    },
    regions_ref => {
        inputs_ref      => [qw{ 1 2 }],
        expected_output => q{--chr} . $SPACE . q{1,2},
    },
    synonyms_file_path => {
        input           => catfile(qw{a synonym_file.tsv}),
        expected_output => q{--synonyms} . $SPACE . catfile(qw{a synonym_file.tsv}),
    },
    vep_features_ref => {
        inputs_ref      => [qw{ tsl hgvs}],
        expected_output => q{--tsl} . $SPACE . q{--hgvs},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&variant_effect_predictor;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

HASHES_OF_ARGUMENTS:
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
