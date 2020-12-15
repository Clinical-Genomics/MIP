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
        q{MIP::Program::Manta} => [qw{ manta_config }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Manta qw{ manta_config };

diag(   q{Test manta_config from Manta}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ configManta.py };

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

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    infile_paths_ref => {
        inputs_ref      => [qw{ infile_path1 infile_path2 infile_path3 }],
        expected_output => q{--bam infile_path1 --bam infile_path2 --bam infile_path3},
    },
    referencefile_path => {
        input           => q{referencefile_path},
        expected_output => q{--referenceFasta referencefile_path},
    },
);

my %specific_argument = (
    call_regions_file_path => {
        input           => catfile(qw{a region_file.bed.gz}),
        expected_output => q{--callRegions} . $SPACE . catfile(qw{a region_file.bed.gz}),
    },
    exome_analysis => {
        input           => 1,
        expected_output => q{--exome},
    },
    outdirectory_path => {
        input           => q{outdirectory_path},
        expected_output => q{--runDir outdirectory_path},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&manta_config;

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
