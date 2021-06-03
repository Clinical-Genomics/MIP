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


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Program::Star_fusion} => [qw{ star_fusion_gtf_file_to_feature_seqs }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Star_fusion qw{ star_fusion_gtf_file_to_feature_seqs };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test star_fusion_gtf_file_to_feature_seqs from Star_fusion.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ gtf_file_to_feature_seqs.pl };

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
    gtf_path => {
        input           => catfile(qw{ a test transcripts_file.gtf }),
        expected_output => q{--gtf_file}
          . $SPACE
          . catfile(qw{ a test transcripts_file.gtf }),
    },
    referencefile_path => {
        input           => catfile(qw{ a test human_reference.fasta }),
        expected_output => q{--genome_fa}
          . $SPACE
          . catfile(qw{ a test  human_reference.fasta }),
    },
    seq_type => {
        input           => q{cDNA},
        expected_output => q{--seqType} . $SPACE . q{cDNA},
    },
);

my %specific_argument = (
    gtf_path => {
        input           => catfile(qw{ a test transcripts_file.gtf }),
        expected_output => q{--gtf_file}
          . $SPACE
          . catfile(qw{ a test transcripts_file.gtf }),
    },
    referencefile_path => {
        input           => catfile(qw{ a test human_reference.fasta }),
        expected_output => q{--genome_fa}
          . $SPACE
          . catfile(qw{ a test  human_reference.fasta }),
    },
    seq_type => {
        input           => q{cDNA},
        expected_output => q{--seqType} . $SPACE . q{cDNA},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&star_fusion_gtf_file_to_feature_seqs;

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
