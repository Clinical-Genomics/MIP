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
        q{MIP::Program::Star}  => [qw{ star_genome_generate }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Star qw{ star_genome_generate };

diag(   q{Test star_genome_generate from Star.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $READ_LENGTH   => 150;
Readonly my $THREAD_NUMBER => 16;

my @function_base_commands = qw{ STAR --runMode genomeGenerate };

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
    fasta_path => {
        input           => catfile(qw{ dir test_file.fasta }),
        expected_output => q{--genomeFastaFiles}
          . $SPACE
          . catfile(qw{ dir test_file.fasta }),
    },
    genome_dir_path => {
        input           => catfile(qw{ dir genome_dir_path }),
        expected_output => q{--genomeDir} . $SPACE . catfile(qw{ dir genome_dir_path }),
    },
    gtf_path => {
        input           => catfile(qw{ dir test_gtf.gtf }),
        expected_output => q{--sjdbGTFfile} . $SPACE . catfile(qw{ dir test_gtf.gtf }),
    },
);

my %specific_argument = (
    fasta_path => {
        input           => catfile(qw{ dir test_file.fasta }),
        expected_output => q{--genomeFastaFiles}
          . $SPACE
          . catfile(qw{ dir test_file.fasta }),
    },
    genome_dir_path => {
        input           => catfile(qw{ dir genome_dir_path }),
        expected_output => q{--genomeDir} . $SPACE . catfile(qw{ dir genome_dir_path }),
    },
    gtf_path => {
        input           => catfile(qw{ dir test_gtf.gtf }),
        expected_output => q{--sjdbGTFfile} . $SPACE . catfile(qw{ dir test_gtf.gtf }),
    },
    read_length => {
        input           => $READ_LENGTH,
        expected_output => q{--sjdbOverhang} . $SPACE . $READ_LENGTH,
    },
    thread_number => {
        input           => $THREAD_NUMBER,
        expected_output => q{--runThreadN} . $SPACE . $THREAD_NUMBER,
    },
);

my $module_function_cref = \&star_genome_generate;

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
