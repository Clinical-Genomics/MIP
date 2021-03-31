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
    my %perl_module = ( q{MIP::Program::Arriba} => [qw{ arriba }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Arriba qw{ arriba };

diag(   q{Test arriba from Arriba.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ arriba };

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
    annotation_file_path => {
        input           => catfile(qw{ my annotations.gtf }),
        expected_output => q{-g} . $SPACE . catfile(qw{ my annotations.gtf }),
    },
    genome_file_path => {
        input           => catfile(qw{ my genome.fa }),
        expected_output => q{-a} . $SPACE . catfile(qw{ my genome.fa }),
    },
    infile_path => {
        input           => catfile(qw{ my aligned.bam }),
        expected_output => q{-x} . $SPACE . catfile(qw{ my aligned.bam }),
    },
    outfile_path => {
        input           => catfile(qw{ my fusions.tsv }),
        expected_output => q{-o} . $SPACE . catfile(qw{ my fusions.tsv }),
    },
);

my %specific_argument = (
    annotation_file_path => {
        input           => catfile(qw{ my annotations.gtf }),
        expected_output => q{-g} . $SPACE . catfile(qw{ my annotations.gtf }),
    },
    blacklist_file_path => {
        input           => catfile(qw{ my fusion blacklist.tsv }),
        expected_output => q{-b} . $SPACE . catfile(qw{ my fusion blacklist.tsv }),
    },
    dna_sv_file_path => {
        input           => catfile(qw{ my wgs_sv_calls.tsv }),
        expected_output => q{-d} . $SPACE . catfile(qw{ my wgs_sv_calls.tsv }),
    },
    discarded_fusion_file_path => {
        input           => catfile(qw{ my fusion.discarded.tsv }),
        expected_output => q{-O} . $SPACE . catfile(qw{ my fusion.discarded.tsv }),
    },
    genome_file_path => {
        input           => catfile(qw{ my genome.fa }),
        expected_output => q{-a} . $SPACE . catfile(qw{ my genome.fa }),
    },
    infile_path => {
        input           => catfile(qw{ my aligned.bam }),
        expected_output => q{-x} . $SPACE . catfile(qw{ my aligned.bam }),
    },
    known_fusion_file_path => {
        input           => catfile(qw{ my known_fusions.tsv}),
        expected_output => q{-k} . $SPACE . catfile(qw{ my known_fusions.tsv}),
    },
    outfile_path => {
        input           => catfile(qw{ my fusions.tsv }),
        expected_output => q{-o} . $SPACE . catfile(qw{ my fusions.tsv }),
    },
    protein_domain_file_path => {
        input           => catfile(qw{ my proteindomains.gff }),
        expected_output => q{-p} . $SPACE . catfile(qw{ my proteindomains.gff }),
    },
    tag_file_path => {
        input           => catfile(qw{ my tags.tsv }),
        expected_output => q{-t} . $SPACE . catfile(qw{ my tags.tsv }),
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&arriba;

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
