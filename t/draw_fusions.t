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
use MIP::Constants qw{ $COMMA $EQUALS $SPACE };
use MIP::Test::Commands qw{ test_function };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Program::Arriba} => [qw{ draw_fusions }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Arriba qw{ draw_fusions };

diag(   q{Test draw_fusions from Arriba.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ draw_fusions.R };

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

my %required_argument = (
    alignment_file_path => {
        input           => catfile(qw{ my aligned.bam }),
        expected_output => q{--alignments} . $EQUALS . catfile(qw{ my aligned.bam }),
    },
    annotation_file_path => {
        input           => catfile(qw{ my annotations.gtf }),
        expected_output => q{--annotation} . $EQUALS . catfile(qw{ my annotations.gtf }),
    },
    fusion_file_path => {
        input           => catfile(qw{ my fusions.tsv }),
        expected_output => q{--fusions} . $EQUALS . catfile(qw{ my fusions.tsv }),
    },
    outfile_path => {
        input           => catfile(qw{ my fusions.pdf }),
        expected_output => q{--output} . $EQUALS . catfile(qw{ my fusions.pdf }),
    },
);

my %specific_argument = (
    alignment_file_path => {
        input           => catfile(qw{ my aligned.bam }),
        expected_output => q{--alignments} . $EQUALS . catfile(qw{ my aligned.bam }),
    },
    annotation_file_path => {
        input           => catfile(qw{ my annotations.gtf }),
        expected_output => q{--annotation} . $EQUALS . catfile(qw{ my annotations.gtf }),
    },
    cytoband_file_path => {
        input           => catfile(qw{ my cytobands.tsv }),
        expected_output => q{--cytobands} . $EQUALS . catfile(qw{ my cytobands.tsv }),
    },
    fusion_file_path => {
        input           => catfile(qw{ my fusions.tsv }),
        expected_output => q{--fusions} . $EQUALS . catfile(qw{ my fusions.tsv }),
    },
    outfile_path => {
        input           => catfile(qw{ my fusions.pdf }),
        expected_output => q{--output} . $EQUALS . catfile(qw{ my fusions.pdf }),
    },
    protein_domain_file_path => {
        input           => catfile(qw{ my protein_domains.gff }),
        expected_output => q{--proteinDomains}
          . $EQUALS
          . catfile(qw{ my protein_domains.gff }),
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&draw_fusions;

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
