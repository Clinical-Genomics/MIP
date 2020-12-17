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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Commands qw{ test_function };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Program::Deeptrio} => [qw{ deeptrio }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Deeptrio qw{ deeptrio };

diag(   q{Test deeptrio from Deeptrio.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ run_deeptrio };

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
    model_type          => {
        input           => q{WES},
        expected_output => q{--model_type WES},
    },
    num_shards          => {
        input           => q{36},
        expected_output => q{--num_shards 36},
    },
    output_gvcf_child   => {
        input           => catfile(q{ dir outfile.g.vcf }),
        expected_output => q{--output_gvcf_child } . catfile(q{ dir outfile.g.vcf }),
    },
    output_gvcf_parent1 => {
        input           => catfile(q{ dir outfile.g.vcf }),
        expected_output => q{--output_gvcf_parent1 } . catfile(q{ dir outfile.g.vcf }),
    },
    output_vcf_child    => {
        input           => catfile(q{ dir outfile.vcf }),
        expected_output => q{--output_vcf_child } . catfile(q{ dir outfile.vcf }),
    },
    output_vcf_parent1  => {
        input           => catfile(q{ dir outfile.vcf }),
        expected_output => q{--output_vcf_parent1 }  . catfile(q{ dir outfile.vcf }),
    },
    reads_child         => {
        input           => catfile(qw{ dir infile.bam }),
        expected_output => q{--reads_child } . catfile(qw{ dir infile.bam }) ,
    },
    reads_parent1       => {
        input           => catfile(qw{ dir infile.bam }),
        expected_output => q{--reads_parent1 } . catfile(qw{ dir infile.bam }),
    },
    referencefile_path  => {
        input           => q{test},
        expected_output => q{--ref test},
    },
    sample_name_child   => {
        input           => q{child},
        expected_output => q{--sample_name_child } . q{child},
    },
    sample_name_parent1 => {
        input           => q{father},
        expected_output => q{--sample_name_parent1 } . q{father},
    },
);

my %specific_argument = (
    bedfile => {
        input           => catfile(qw{ dir infile.bed }),
        expected_output => q{--regions } . catfile(qw{ dir infile.bed }),
    },
    output_gvcf_parent2 => {
        input           => catfile(q{ dir outfile.g.vcf }),
        expected_output => q{--output_gvcf_parent2 } . catfile(q{ dir outfile.g.vcf }),
    },
    output_vcf_parent2  => {
        input => catfile(q{ dir outfile.vcf }),
        expected_output => q{--output_vcf_parent2 }  . catfile(q{ dir outfile.vcf }),
    },
    reads_parent2 => {
        input           => catfile(qw{ dir infile.bam }),
        expected_output => q{--reads_parent2 } . catfile(qw{ dir infile.bam }),
    },
    sample_name_parent2 => {
        input => q{mother},
        expected_output => q{--sample_name_parent2 } . q{mother},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&deeptrio;

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
