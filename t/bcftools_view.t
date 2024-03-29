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
use MIP::Constants qw{ $SPACE $COMMA };
use MIP::Test::Commands qw{ test_function };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = ( q{MIP::Program::Bcftools} => [qw{ bcftools_view }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Bcftools qw{ bcftools_view };

diag(   q{Test bcftools_view from Bcftools.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ bcftools view };

Readonly my $MAX_ALLELES => 2;
Readonly my $MIN_AC      => 1;
Readonly my $MIN_ALLELES => 2;

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
my %specific_argument = (
    apply_filters_ref => {
        inputs_ref      => [qw{ filter1 filter2 }],
        expected_output => q{--apply-filters filter1,filter2},
    },
    exclude => {
        input           => q{%QUAL<10 || (RPB<0.1 && %QUAL<15)},
        expected_output => q{--exclude %QUAL<10 || (RPB<0.1 && %QUAL<15)},
    },
    exclude_types_ref => {
        inputs_ref      => [qw{ type1 type2 }],
        expected_output => q{--exclude-types type1,type2},
    },
    genotype => {
        input           => q{het},
        expected_output => q{--genotype het},
    },
    header_only => {
        input           => 1,
        expected_output => q{--header-only},
    },
    include => {
        input           => q{INFO/CSQ[*]~":p[.]"},
        expected_output => q{--include INFO/CSQ[*]~":p[.]"},
    },
    infile_path => {
        input           => q{infile.test},
        expected_output => q{infile.test},
    },
    max_alleles => {
        input           => $MAX_ALLELES,
        expected_output => q{--max-alleles } . $MAX_ALLELES,
    },
    min_ac => {
        input           => $MIN_AC,
        expected_output => q{--min-ac } . $MIN_AC,
    },
    min_alleles => {
        input           => $MIN_ALLELES,
        expected_output => q{--min-alleles } . $MIN_ALLELES,
    },
    no_header => {
        input           => 1,
        expected_output => q{--no-header},
    },
    outfile_path => {
        input           => q{outfile.txt},
        expected_output => q{--output-file outfile.txt},
    },
    output_type => {
        input           => q{v},
        expected_output => q{--output-type v},
    },
    types => {
        input           => q{snps,indel},
        expected_output => q{--types snps,indel},
    },
    targets => {
        input           => q{MT},
        expected_output => q{--targets MT},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&bcftools_view;

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
        }
    );
}

done_testing();
