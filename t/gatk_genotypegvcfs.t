#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname  };
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
our $VERSION = 1.02;

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
        q{MIP::Program::Gatk}  => [qw{ gatk_genotypegvcfs }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Gatk qw{ gatk_genotypegvcfs };

diag(   q{Test gatk_genotypegvcfs from Gatk.pm v}
      . $MIP::Program::Gatk::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ gatk GenotypeGVCFs };

my %base_argument = (
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
    stderrfile_path => {
        input           => q{stderrfile.test},
        expected_output => q{2> stderrfile.test},
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    infile_path => {
        input           => q{gendb://} . catdir(qw{ path to my_db }),
        expected_output => q{--variant gendb://} . catdir(qw{ path to my_db }),
    },
    outfile_path => {
        input           => catfile(qw{ my outfile }),
        expected_output => q{--output } . catfile(qw{ my outfile }),
    },
    referencefile_path => {
        input           => catfile(qw{ my genome }),
        expected_output => q{--reference } . catdir(qw{ my genome }),
    },
);

my %specific_argument = (
    dbsnp_path => {
        input           => catfile(qw{ dir grch37_dbsnp_-138-.vcf }),
        expected_output => q{--dbsnp } . catfile(qw{ dir grch37_dbsnp_-138-.vcf }),
    },
    include_nonvariant_sites => {
        input           => 1,
        expected_output => q{--include-non-variant-sites},
    },
    infile_path => {
        input           => q{gendb://} . catdir(qw{ path to my_db }),
        expected_output => q{--variant gendb://} . catdir(qw{ path to my_db }),
    },
    outfile_path => {
        input           => catfile(qw{ my outfile }),
        expected_output => q{--output } . catfile(qw{ my outfile }),
    },
    referencefile_path => {
        input           => catfile(qw{ my genome }),
        expected_output => q{--reference } . catdir(qw{ my genome }),
    },
    use_new_qual_calculator => {
        input           => 1,
        expected_output => q{--use-new-qual-calculator},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&gatk_genotypegvcfs;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

ARGUMENT_HASH_REF:
foreach my $argument_href (@arguments) {
    my @commands = test_function(
        {
            argument_href              => $argument_href,
            base_commands_index        => 1,
            do_test_base_command       => 1,
            function_base_commands_ref => \@function_base_commands,
            module_function_cref       => $module_function_cref,
            required_argument_href     => \%required_argument,
        }
    );
}

done_testing();
