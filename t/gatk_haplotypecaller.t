#!/usr/bin/env perl

use Modern::Perl qw{ 2014 };
use warnings qw{ FATAL utf8 };
use autodie;
use 5.018;
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };

use FindBin qw{ $Bin };
use File::Basename qw{ dirname basename };
use File::Spec::Functions qw{ catdir catfile };
use Getopt::Long;
use Test::More;
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = 1.0.0;

## Constants
Readonly my $SPACE   => q{ };
Readonly my $NEWLINE => qq{\n};
Readonly my $COMMA   => q{,};

### User Options
GetOptions(

    # Display help text
    q{h|help} => sub {
        done_testing();
        say {*STDOUT} $USAGE;
        exit;
    },

    # Display version number
    q{v|version} => sub {
        done_testing();
        say {*STDOUT} $NEWLINE
          . basename($PROGRAM_NAME)
          . $SPACE
          . $VERSION
          . $NEWLINE;
        exit;
    },
    q{vb|verbose} => $VERBOSE,
  )
  or (
    done_testing(),
    help(
        {
            USAGE     => $USAGE,
            exit_code => 1,
        }
    )
  );

BEGIN {

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module;

    $perl_module{q{MIP::Script::Utils}} = [qw{ help }];

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

## Modules
    my @modules = (q{MIP::Program::Alignment::Gatk});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Program::Alignment::Gatk qw{ gatk_haplotypecaller };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test gatk_haplotypecaller from Alignment::Gatk.pm v}
      . $MIP::Program::Alignment::Gatk::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $MITOCHONDRIA_PLOIDY                           => 3;
Readonly my $SAMPLE_PLOIDY                                 => 3;
Readonly my $STANDARD_MIN_CONFIDENCE_THRESHOLD_FOR_CALLING => 10;
Readonly my $VARIANT_INDEX_PARAMETER                       => 128_000;

## Base arguments
my $function_base_command = q{--analysis_type HaplotypeCaller};

my %base_argument = (
    stderrfile_path => {
        input           => q{stderrfile.test},
        expected_output => q{2> stderrfile.test},
    },
    FILEHANDLE => {
        input           => undef,
        expected_output => $function_base_command,
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    intervals_ref => {
        inputs_ref      => [qw{ chr1 chr2}],
        expected_output => q{--intervals chr1 --intervals chr2},
    },
    referencefile_path => {
        input           => catfile(qw{reference_dir human_genome_build.fasta }),
        expected_output => q{--reference_sequence}
          . catfile(qw{reference_dir human_genome_build.fasta }),
    },
    annotations_ref => {
        inputs_ref => [qw{ BaseQualityRankSumTest ChromosomeCounts }],
        expected_output =>
          q{--annotation BaseQualityRankSumTest --annotation ChromosomeCounts},
    },
    infile_path => {
        input           => catfile(qw{ dir infile.bam }),
        expected_output => q{--input_file } . catfile(qw{ dir infile.bam }),
    },
    outfile_path => {
        input           => catfile(qw{ dir outfile.bam }),
        expected_output => q{--out } . catfile(qw{ dir outfile.bam }),
    },
);

my %specific_argument = (
    intervals_ref => {
        inputs_ref      => [qw{ chr1 chr2}],
        expected_output => q{--intervals chr1 --intervals chr2},
    },
    referencefile_path => {
        input           => catfile(qw{reference_dir human_genome_build.fasta }),
        expected_output => q{--reference_sequence }
          . catfile(qw{ reference_dir human_genome_build.fasta }),
    },
    annotations_ref => {
        inputs_ref => [qw{ BaseQualityRankSumTest ChromosomeCounts }],
        expected_output =>
          q{--annotation BaseQualityRankSumTest --annotation ChromosomeCounts},
    },
    infile_path => {
        input           => catfile(qw{ dir infile.bam }),
        expected_output => q{--input_file } . catfile(qw{ dir infile.bam }),
    },
    outfile_path => {
        input           => catfile(qw{ dir outfile.bam }),
        expected_output => q{--out } . catfile(qw{ dir outfile.bam }),
    },
    dbsnp => {
        input           => catfile(qw{ dir GRCh37_dbsnp_-138-.vcf }),
        expected_output => q{--dbsnp }
          . catfile(qw{ dir GRCh37_dbsnp_-138-.vcf }),
    },
    sample_ploidy => {
        input           => $MITOCHONDRIA_PLOIDY,
        expected_output => q{--sample_ploidy} . $SPACE . $MITOCHONDRIA_PLOIDY,
    },
    standard_min_confidence_threshold_for_calling => {
        input           => $STANDARD_MIN_CONFIDENCE_THRESHOLD_FOR_CALLING,
        expected_output => q{--standard_min_confidence_threshold_for_calling }
          . $STANDARD_MIN_CONFIDENCE_THRESHOLD_FOR_CALLING,
    },
    dont_use_soft_clipped_bases => {
        input           => 1,
        expected_output => q{--dontUseSoftClippedBases},
    },
    pcr_indel_model => {
        input           => q{NONE},
        expected_output => q{--pcr_indel_model NONE},
    },
    variant_index_parameter => {
        input           => $VARIANT_INDEX_PARAMETER,
        expected_output => q{--variant_index_parameter }
          . $VARIANT_INDEX_PARAMETER,
    },
    emit_ref_confidence => {
        input           => q{GVCF},
        expected_output => q{--emitRefConfidence GVCF},
    },
    variant_index_type => {
        input           => q{LINEAR},
        expected_output => q{--variant_index_type LINEAR},
    },
    sample_ploidy => {
        input           => $SAMPLE_PLOIDY,
        expected_output => q{--sample_ploidy} . $SPACE . $SAMPLE_PLOIDY,
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&gatk_haplotypecaller;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

ARGUMENT_HASH_REF:
foreach my $argument_href (@arguments) {
    my @commands = test_function(
        {
            argument_href          => $argument_href,
            required_argument_href => \%required_argument,
            module_function_cref   => $module_function_cref,
            function_base_command  => $function_base_command,
            do_test_base_command   => 1,
        }
    );
}

done_testing();

######################
####SubRoutines#######
######################

sub build_usage {

## Function  : Build the USAGE instructions
## Returns   :
## Arguments : $program_name => Name of the script

    my ($arg_href) = @_;

    ## Default(s)
    my $program_name;

    my $tmpl = {
        program_name => {
            default     => basename($PROGRAM_NAME),
            strict_type => 1,
            store       => \$program_name,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return <<"END_USAGE";
 $program_name [options]
    -vb/--verbose Verbose
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}
