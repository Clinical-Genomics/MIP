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
        q{MIP::Program::Picardtools} => [qw{ picardtools_sortvcf }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Picardtools qw{ picardtools_sortvcf };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test picardtools_sortvcf from Picardtools.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $MIN_GQ    => 20;
Readonly my $MIN_DEPTH => 10;

## Base arguments
my @function_base_commands = qw{ picard SortVcf };

my %base_argument = (
    stderrfile_path => {
        input           => q{stderrfile.test},
        expected_output => q{2> stderrfile.test},
    },
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    infile_paths_ref => {
        inputs_ref =>
          [ catfile(qw{ dir infile_1.vcf }), catfile(qw{ dir infile_2.vcf }) ],
        expected_output => q{-INPUT}
          . $SPACE
          . catfile(qw{ dir infile_1.vcf })
          . $SPACE
          . q{-INPUT}
          . $SPACE
          . catfile(qw{ dir infile_2.vcf }),
    },
    sequence_dictionary => {
        input           => catfile(qw{ dir ref.dict }),
        expected_output => q{-SEQUENCE_DICTIONARY}
          . $SPACE
          . catfile(qw{ dir ref.dict  }),
    },
    outfile_path => {
        input           => catfile(qw{ dir outfile.vcf }),
        expected_output => q{-OUTPUT} . $SPACE . catfile(qw{ dir outfile.vcf }),
    },
    referencefile_path => {
        input           => catfile(qw{ references grch37_homo_sapiens_-d5-.fasta }),
        expected_output => q{-REFERENCE_SEQUENCE}
          . $SPACE
          . catfile(qw{ references grch37_homo_sapiens_-d5-.fasta }),
    },
);

my %specific_argument = (
    infile_paths_ref => {
        inputs_ref =>
          [ catfile(qw{ dir infile_1.vcf }), catfile(qw{ dir infile_2.vcf }) ],
        expected_output => q{-INPUT}
          . $SPACE
          . catfile(qw{ dir infile_1.vcf })
          . $SPACE
          . q{-INPUT}
          . $SPACE
          . catfile(qw{ dir infile_2.vcf }),
    },
    sequence_dictionary => {
        input           => catfile(qw{ dir ref.dict }),
        expected_output => q{-SEQUENCE_DICTIONARY}
          . $SPACE
          . catfile(qw{ dir ref.dict  }),
    },
    outfile_path => {
        input           => catfile(qw{ dir outfile.vcf }),
        expected_output => q{-OUTPUT} . $SPACE . catfile(qw{ dir outfile.vcf }),
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&picardtools_sortvcf;

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
