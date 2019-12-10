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
        q{MIP::Program::Gatk}  => [qw{ gatk_common_options }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Gatk qw{ gatk_common_options };

diag(   q{Test gatk_common_options from Gatk.pm v}
      . $MIP::Program::Gatk::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my @function_base_commands = qw{ program };

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    commands_ref => {
        inputs_ref      => [qw{ program }],
        expected_output => [qw{ program }],
    },
);

my %specific_argument = (
    intervals_ref => {
        inputs_ref      => [qw{ chr1 chr2 }],
        expected_output => q{--intervals chr1 --intervals chr2},
    },
    pedigree => {
        input           => catfile(qw{ a pedigree }),
        expected_output => q{--pedigree } . catfile(qw{ a pedigree }),
    },
    read_filters_ref => {
        inputs_ref      => [qw{ MalformedRead BadCigar}],
        expected_output => q{--read-filter MalformedRead --read-filter BadCigar},
    },
    referencefile_path => {
        input           => catfile(qw{reference_dir human_genome_build.fasta }),
        expected_output => q{--reference }
          . catfile(qw{reference_dir human_genome_build.fasta }),
    },
    temp_directory => {
        input           => catdir(qw{ a dir }),
        expected_output => q{--tmp-dir } . catdir(qw{ a dir }),
    },
    verbosity => {
        input           => q{INFO},
        expected_output => q{--verbosity INFO},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&gatk_common_options;

## Test both base and function specific arguments
my @arguments = ( \%specific_argument );

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
