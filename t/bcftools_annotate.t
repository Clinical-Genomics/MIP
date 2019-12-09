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
        q{MIP::Program::Bcftools} => [qw{ bcftools_annotate }],
        q{MIP::Test::Fixtures}    => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Bcftools qw{ bcftools_annotate };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test bcftools_annotate from Bcftools.pm v}
      . $MIP::Program::Bcftools::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ bcftools };

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
my %required_argument = ();

my %specific_argument = (
    annotations_file_path => {
        input           => q{annotation_file.tsv},
        expected_output => q{--annotations annotation_file.tsv},
    },
    columns_name => {
        input           => q{Chrom,Pos,Ref,Alt,-,CADD},
        expected_output => q{--columns Chrom,Pos,Ref,Alt,-,CADD},
    },
    headerfile_path => {
        input           => q{headerlines},
        expected_output => q{--header-lines headerlines},
    },
    infile_path => {
        input           => q{infile.test},
        expected_output => q{infile.test},
    },
    outfile_path => {
        input           => q{outfile.txt},
        expected_output => q{--output outfile.txt},
    },
    output_type => {
        input           => q{v},
        expected_output => q{--output-type v},
    },
    remove_ids_ref => {
        inputs_ref      => [qw{ idref1 idref2 idref3 }],
        expected_output => q{--remove idref1,idref2,idref3},
    },
    samples_file_path => {
        input           => q{samplesfile},
        expected_output => q{--samples-file samplesfile},
    },
    set_id => {
        input           => q{id},
        expected_output => q{--set-id id},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&bcftools_annotate;

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
