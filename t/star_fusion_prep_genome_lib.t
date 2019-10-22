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
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.01;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $THREAD_NUMBER => 16;

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Program::Star_fusion} => [qw{ star_fusion_prep_genome_lib }],
        q{MIP::Test::Fixtures}       => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Star_fusion qw{ star_fusion_prep_genome_lib };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test star_fusion_prep_genome_lib from Star_fusion.pm v}
      . $MIP::Program::Star_fusion::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ prep_genome_lib.pl };

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
    blast_pairs_file_path => {
        input           => catfile(qw{ a test blast_pairs_file.tab }),
        expected_output => q{--blast_pairs}
          . $SPACE
          . catfile(qw{ a test blast_pairs_file.tab }),
    },
    gtf_path => {
        input           => catfile(qw{ a test transcripts_file.gtf }),
        expected_output => q{--gtf} . $SPACE . catfile(qw{ a test transcripts_file.gtf }),
    },
    output_dir_path => {
        input           => catdir(qw{ a test outdir }),
        expected_output => q{--output_dir} . $SPACE . catdir(qw{ a test outdir }),
    },
    referencefile_path => {
        input           => catfile(qw{ a test human_reference.fasta }),
        expected_output => q{--genome_fa}
          . $SPACE
          . catfile(qw{ a test  human_reference.fasta }),
    },
);

my %specific_argument = (
    blast_pairs_file_path => {
        input           => catfile(qw{ a test blast_pairs_file.tab }),
        expected_output => q{--blast_pairs}
          . $SPACE
          . catfile(qw{ a test blast_pairs_file.tab }),
    },
    gtf_path => {
        input           => catfile(qw{ a test transcripts_file.gtf }),
        expected_output => q{--gtf} . $SPACE . catfile(qw{ a test transcripts_file.gtf }),
    },
    output_dir_path => {
        input           => catdir(qw{ a test outdir }),
        expected_output => q{--output_dir} . $SPACE . catdir(qw{ a test outdir }),
    },
    referencefile_path => {
        input           => catfile(qw{ a test human_reference.fasta }),
        expected_output => q{--genome_fa}
          . $SPACE
          . catfile(qw{ a test  human_reference.fasta }),
    },
    thread_number => {
        input           => $THREAD_NUMBER,
        expected_output => q{--cpu} . $SPACE . $THREAD_NUMBER,
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&star_fusion_prep_genome_lib;

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
