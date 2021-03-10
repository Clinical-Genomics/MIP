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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Commands qw{ test_function };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = ( q{MIP::Program::Pdfmerger} => [qw{ pdfmerger }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Pdfmerger qw{ pdfmerger };

diag(   q{Test pdfmerger from Pdfmerger.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ pdfmerger };

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
    stdinfile_path => {
        input           => q{stdinfile.test},
        expected_output => q{< stdinfile.test},
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
        inputs_ref      => [qw{ fusion_file1.pdf fusion_file2.pdf }],
        expected_output => q{--infile fusion_file1.pdf --infile fusion_file2.pdf},
    },
    orientation => {
        input           => q{landscape},
        expected_output => q{--orientation landscape},
    },
    outdir_path => {
        input           => catdir(qw{ outfolder path }),
        expected_output => q{--outfolder} . $SPACE . catdir(qw{ outfolder path }),
    },
    outfile_name => {
        input           => q{outfile.pdf},
        expected_output => q{--outfile outfile.pdf},
    },
);

my %specific_argument = (
    infile_paths_ref => {
        inputs_ref      => [qw{fusion_file1.pdf fusion_file2.pdf}],
        expected_output => q{--infile fusion_file1.pdf --infile fusion_file2.pdf},
    },
    orientation => {
        input           => q{landscape},
        expected_output => q{--orientation landscape},
    },
    outdir_path => {
        input           => catdir(qw{outfolder path}),
        expected_output => q{--outfolder} . $SPACE . catdir(qw{outfolder path}),
    },
    outfile_name => {
        input           => q{outfile.pdf},
        expected_output => q{--outfile outfile.pdf},
    },
    write_filenames => {
        input           => 1,
        expected_output => q{--write-filenames},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&pdfmerger;

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
