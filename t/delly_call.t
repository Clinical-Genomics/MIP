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
Readonly my $MIN_MAP_QUAL       => q{20};
Readonly my $INSERT_SIZE_CUTOFF => q{15};

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Program::Variantcalling::Delly} => [qw{ delly_call }],
        q{MIP::Test::Fixtures}                 => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Variantcalling::Delly qw{ delly_call };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test delly_call from Delly.pm v}
      . $MIP::Program::Variantcalling::Delly::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ delly call };

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
    infile_path => {
        input           => catfile(qw{ file_path prefix_contig_suffix }),
        expected_output => catfile(qw{ file_path prefix_contig_suffix }),
    },
    referencefile_path => {
        input           => catfile(qw{reference_dir human_genome_build.fasta }),
        expected_output => q{--genome}
          . $SPACE
          . catfile(qw{reference_dir human_genome_build.fasta }),
    },
    sv_type => {
        input           => q{DEL},
        expected_output => q{--type DEL},
    },
);

my %specific_argument = (
    exclude_file_path => {
        input           => q{delly_exclude_file},
        expected_output => q{--exclude delly_exclude_file},
    },
    genotypefile_path => {
        input => catfile(qw{ outfile_path prefix_contig_infile_suffix_sv-type_suffix }),
        expected_output => q{--vcffile}
          . $SPACE
          . catfile(qw{ outfile_path prefix_contig_infile_suffix_sv-type_suffix }),
    },
    mad_cutoff => {
        input           => $INSERT_SIZE_CUTOFF,
        expected_output => q{--mad-cutoff} . $SPACE . $INSERT_SIZE_CUTOFF,
    },
    mapping_qual => {
        input           => $MIN_MAP_QUAL,
        expected_output => q{--map-qual} . $SPACE . $MIN_MAP_QUAL,
    },
    small_indel => {
        input           => 1,
        expected_output => q{--i},
    },
    outfile_path => {
        input => catfile(qw{ outfile_path prefix_contig_infile_suffix_sv-type }),
        expected_output => q{--outfile}
          . $SPACE
          . catfile(qw{ outfile_path prefix_contig_infile_suffix_sv-type }),
    },

);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&delly_call;

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
