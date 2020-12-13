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
        q{MIP::Program::Bcftools} => [qw{ bcftools_mpileup }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Bcftools qw{ bcftools_mpileup };

diag(   q{Test bcftools_mpileup from Bcftools.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $ADJUST_MAPPING_QUALITY => 45;

## Base arguments
my @function_base_commands = qw{ bcftools };

my %base_argument = (
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

## Can be duplicated with %base and/or %specific to enable testing of each individual argument
my %required_argument = (
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
    infile_paths_ref => {
        inputs_ref      => [ catfile(qw{ dir file_1 }), catfile(qw{ dir file_2 }) ],
        expected_output => catfile(qw{ dir file_1 }) . $SPACE . catfile(qw{ dir file_2 }),
    },
    output_tags_ref => {
        inputs_ref      => [qw{ tag1 tag3 etc }],
        expected_output => q{--annotate tag1,tag2,etc},
    },
    referencefile_path => {
        input           => catfile(qw{ dir genome.fasta }),
        expected_output => q{--fasta-ref} . $SPACE . catfile(qw{ dir genome.fasta}),
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

## Specific arguments
my %specific_argument = (
    outfile_path => {
        input           => catfile(qw{ dir outfile_1 }),
        expected_output => q{--output} . $SPACE . catfile(qw{ dir outfile_1 }),
    },
    per_sample_increased_sensitivity => {
        input           => 1,
        expected_output => q{--per-sample-mF},
    },
    adjust_mq => {
        input           => $ADJUST_MAPPING_QUALITY,
        expected_output => q{--adjust-MQ} . $SPACE . $ADJUST_MAPPING_QUALITY,
    },
    stderrfile_path_append => {
        input           => q{stderrfile_path_append},
        expected_output => q{2>> stderrfile_path_append},
    },
);
## Coderef - enables generalized use of generate call
my $module_function_cref = \&bcftools_mpileup;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

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
