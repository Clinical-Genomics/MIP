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
    my %perl_module = ( q{MIP::Program::Expansionhunter} => [qw{ expansionhunter }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Expansionhunter qw{ expansionhunter };

diag(   q{Test expansionhunter from Expansionhunter.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Module test constants
Readonly my $REGION_EXTENSION_LENGTH => 1500;

## Base arguments
my @function_base_commands = qw{ ExpansionHunter };

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
        input           => catfile(qw{ a test path }),
        expected_output => q{--reads} . $SPACE . catfile(qw{ a test path }),
    },
    outfile_path_prefix => {
        input           => catfile(qw{ a test outfile_path_prefix }),
        expected_output => q{--output-prefix} . $SPACE . catfile(qw{ a test outfile_path_prefix }),
    },
    reference_genome_path => {
        input           => catfile(qw{ a test fasta }),
        expected_output => q{--reference} . $SPACE . catfile(qw{ a test fasta }),
    },
    variant_catalog_file_path => {
        input           => catdir(qw{ a test variant_catalog_file.json }),
        expected_output => q{--variant-catalog}
          . $SPACE
          . catdir(qw{ a test variant_catalog_file.json }),
    },
);

my %specific_argument = (
    infile_path => {
        input           => catfile(qw{ a test path }),
        expected_output => q{--reads} . $SPACE . catfile(qw{ a test path }),
    },
    log_level => {
        input           => q{debug},
        expected_output => q{--log-level} . $SPACE . q{debug},
    },
    outfile_path_prefix => {
        input           => catfile(qw{ a test outfile_path_prefix }),
        expected_output => q{--output-prefix} . $SPACE . catfile(qw{ a test outfile_path_prefix }),
    },
    reference_genome_path => {
        input           => catfile(qw{ a test fasta }),
        expected_output => q{--reference} . $SPACE . catfile(qw{ a test fasta }),
    },
    region_extension_length => {
        input           => $REGION_EXTENSION_LENGTH,
        expected_output => q{--region-extension-length} . $SPACE . $REGION_EXTENSION_LENGTH,
    },
    sex => {
        input           => q{female},
        expected_output => q{--sex female},
    },
    threads => {
        input           => 2,
        expected_output => q{--threads 2},
    },
    variant_catalog_file_path => {
        input           => catdir(qw{ a test variant_catalog_file.json }),
        expected_output => q{--variant-catalog}
          . $SPACE
          . catdir(qw{ a test variant_catalog_file.json }),
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&expansionhunter;

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
