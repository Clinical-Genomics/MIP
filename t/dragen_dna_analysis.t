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
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.01;

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
        q{MIP::Program::Dragen} => [qw{ dragen_dna_analysis }],
        q{MIP::Test::Fixtures}  => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Dragen qw{ dragen_dna_analysis };

diag(   q{Test dragen_dna_analysis from Dragen.pm v}
      . $MIP::Program::Dragen::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $THREAD_NUMBER => 2;

## Base arguments
my @function_base_commands = qw{ dragen };

my %base_argument = (
    FILEHANDLE => {
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
    dragen_hash_ref_dir_path => {
        input           => catdir(qw{a dragen_ref_hash_table_dir}),
        expected_output => q{--ref-dir}
          . $SPACE
          . catdir(qw{a dragen_ref_hash_table_dir}),
    },
    outdirectory_path => {
        input           => catdir(qw{ outdirectory path }),
        expected_output => q{--output-directory}
          . $SPACE
          . catdir(qw{ outdirectory path }),
    },
    outfile_prefix => {
        input           => q{sample_file_prefix},
        expected_output => q{--output-file-prefix} . $SPACE . q{sample_file_prefix},
    },
);

my %specific_argument = (
    alignment_output_format => {
        input           => q{BAM},
        expected_output => q{--output-format} . $SPACE . q{BAM},
    },
    alt_aware => {
        input           => 1,
        expected_output => q{--alt-aware} . $SPACE . q{true},
    },
    bam_file_path => {
        input           => catfile(qw{ a sample.bam }),
        expected_output => q{-b} . $SPACE . catfile(qw{ a sample.bam }),
    },
    cnv_enable_self_normalization => {
        input           => 1,
        expected_output => q{--cnv-enable-self-normalization} . $SPACE . q{true},
    },
    combine_samples_by_name => {
        input           => 1,
        expected_output => q{--combine-samples-by-name} . $SPACE . q{true},
    },
    dbsnp_file_path => {
        input           => catfile(qw{ a dbsnp_file_path.vcf.gz }),
        expected_output => q{--dbsnp} . $SPACE . catfile(qw{ a dbsnp_file_path.vcf.gz }),
    },
    disable_vcf_compression => {
        input           => 1,
        expected_output => q{--enable-vcf-compression} . $SPACE . q{false},
    },
    dragen_hash_ref_dir_path => {
        input           => catdir(qw{a dragen_ref_hash_table_dir}),
        expected_output => q{--ref-dir}
          . $SPACE
          . catdir(qw{a dragen_ref_hash_table_dir}),
    },
    enable_duplicate_marking => {
        input           => 1,
        expected_output => q{--enable-bam-indexing} . $SPACE . q{true},
    },
    enable_cnv => {
        input           => 1,
        expected_output => q{--enable-cnv} . $SPACE . q{true},
    },
    enable_duplicate_marking => {
        input           => 1,
        expected_output => q{--enable-duplicate-marking} . $SPACE . q{true},
    },
    enable_combinegvcfs => {
        input           => 1,
        expected_output => q{--enable-combinegvcfs} . $SPACE . q{true},
    },
    enable_joint_genotyping => {
        input           => 1,
        expected_output => q{--enable-joint-genotyping} . $SPACE . q{true},
    },
    enable_map_align => {
        input           => 1,
        expected_output => q{--enable-map-align} . $SPACE . q{true},
    },
    enable_map_align_output => {
        input           => 1,
        expected_output => q{--enable-map-align-output} . $SPACE . q{true},
    },
    enable_multi_sample_gvcf => {
        input           => 1,
        expected_output => q{--enable-multi-sample-gvcf} . $SPACE . q{true},
    },
    enable_sampling => {
        input           => 1,
        expected_output => q{--enable-sampling} . $SPACE . q{true},
    },
    enable_sort => {
        input           => 1,
        expected_output => q{--enable-sort} . $SPACE . q{true},
    },
    enable_variant_caller => {
        input           => 1,
        expected_output => q{--enable-variant-caller} . $SPACE . q{true},
    },
    fastq_infile_path => {
        input           => catfile(qw{ a read_1.fastq.gz }),
        expected_output => q{-1} . $SPACE . catfile(qw{ a read_1.fastq.gz }),
    },
    fastq_list_all_samples => {
        input           => 1,
        expected_output => q{--fastq-list-all-samples} . $SPACE . q{true},
    },
    fastq_list_file_path => {
        input           => catfile(qw{ a fastq_list.csv }),
        expected_output => q{--fastq-list} . $SPACE . catfile(qw{ a fastq_list.csv }),
    },
    fastq_list_sample_id => {
        input           => q{sample-1},
        expected_output => q{--fastq-list-sample-id} . $SPACE . q{sample-1},
    },
    force => {
        input           => 1,
        expected_output => q{--force},
    },
    is_fastq_interleaved => {
        input           => 1,
        expected_output => q{--interleaved},
    },
    outdirectory_path => {
        input           => catdir(qw{ outdirectory path }),
        expected_output => q{--output-directory}
          . $SPACE
          . catdir(qw{ outdirectory path }),
    },
    outfile_prefix => {
        input           => q{sample_file_prefix},
        expected_output => q{--output-file-prefix} . $SPACE . q{sample_file_prefix},
    },
    pedigree_file_path => {
        input           => catfile(qw{ a family.ped }),
        expected_output => q{--vc-pedigree} . $SPACE . catfile(qw{ a family.ped }),
    },
    sample_gvcf_file_paths_ref => {
        inputs_ref      => [qw{ sample-1.gvcf sample-2.gvcf }],
        expected_output => q{--variant sample-1.gvcf --variant sample-2.gvcf},
    },
    sample_id => {
        input           => q{sample-1},
        expected_output => q{--vc-sample-name} . $SPACE . q{sample-1},
    },
    second_fastq_infile_path => {
        input           => catfile(qw{ a read_2.fastq.gz }),
        expected_output => q{-2} . $SPACE . catfile(qw{ a read_2.fastq.gz }),
    },
    vc_emit_ref_confidence => {
        input           => q{GVCF},
        expected_output => q{--vc-emit-ref-confidence} . $SPACE . q{GVCF},
    },
    vc_enable_gatk_acceleration => {
        input           => 1,
        expected_output => q{--vc-enable-gatk-acceleration} . $SPACE . q{true},
    },
    vc_target_bed_file_path => {
        input           => catfile(qw{ a targets.bed }),
        expected_output => q{--vc-target-bed} . $SPACE . catfile(qw{ a targets.bed }),
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&dragen_dna_analysis;

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
