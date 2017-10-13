#!/usr/bin/env perl

use Modern::Perl qw{ 2014 };
use warnings qw{ FATAL utf8 };
use autodie;
use 5.018;    #Require at least perl 5.18
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };

use FindBin qw{ $Bin };    #Find directory of script
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
Readonly my $SPACE       => q{ };
Readonly my $NEWLINE     => qq{\n};
Readonly my $COMMA       => q{,};
Readonly my $COVERAGE    => 90;
Readonly my $TS_FILTER   => 99;
Readonly my $CPU_THREADS => 8;

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
    my @modules = (q{MIP::Program::Variantcalling::Gatk});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Program::Variantcalling::Gatk qw{ gatk_applyrecalibration };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test gatk_applyrecalibration from Gatk v}
      . $MIP::Program::Variantcalling::Gatk::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my $function_base_command = q{--analysis_type ApplyRecalibration};

my %base_argument = (
    FILEHANDLE => {
        input           => undef,
        expected_output => $function_base_command,
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    infile_path => {
        input           => catfile(qw{ path_to_analysis_dir infile.vcf }),
        expected_output => q{--input} . $SPACE
          . catfile(qw{ path_to_analysis_dir infile.vcf }),
    },
    outfile_path => {
        input           => catfile(qw{ path_to_analysis_dir outfile.vcf }),
        expected_output => q{--out} . $SPACE
          . catfile(qw{ path_to_analysis_dir outfile.vcf }),
    },
    referencefile_path => {
        input           => catfile(qw{reference_dir human_genome_build.fasta }),
        expected_output => q{--reference_sequence} . $SPACE
          . catfile(qw{reference_dir human_genome_build.fasta }),
    },
    mode => {
        inputs_ref      => q{SNP},
        expected_output => q{--mode SNP},
    },
);

my %specific_argument = (
    downsample_to_coverage => {
        input           => $COVERAGE,
        expected_output => q{--downsample_to_coverage } . $COVERAGE,
    },
    gatk_disable_auto_index_and_file_lock => {
        input => 1,
        expected_output =>
          q{--disable_auto_index_creation_and_locking_when_reading_rods},
    },
    logging_level => {
        input           => q{INFO},
        expected_output => q{--logging_level INFO},
    },
    intervals_ref => {
        inputs_ref => [qw{ chr1 chr2 chr3 }],
        expected_output =>
          q{--intervals chr1 --intervals chr2 --intervals chr3},
    },
    pedigree => {
        input           => catfile(qw{ dir pedigree.fam }),
        expected_output => q{--pedigree} . $SPACE . catfile(qw{ dir pedigree.fam }),
    },
    pedigree_validation_type => {
        input           => q{SILENT},
        expected_output => q{--pedigreeValidationType SILENT},
    },
    read_filters_ref => {
        inputs_ref => [qw{ MaxReadLength MappingQuality }],
        expected_output =>
          q{--read_filter MaxReadLength --read_filter MappingQuality},
    },
    static_quantized_quals_ref => {
        inputs_ref => [qw{ 10 20 }],
        expected_output =>
          q{--static_quantized_quals 10 --static_quantized_quals 20},
    },
    base_quality_score_recalibration_file => {
        input           => catfile(qw{ dir recalibration_report.grp }),
        expected_output => q{--BQSR} . $SPACE
          . catfile(qw{ dir recalibration_report.grp }),
    },
    recal_file_path => {
        input           => catfile(qw{ dir recal_file.intervals }),
        expected_output => q{--recal_file} . $SPACE
          . catfile(qw{ dir recal_file.intervals }),
    },
    tranches_file_path => {
        input           => catfile(qw{ dir tranchesfile.tranches }),
        expected_output => q{--tranches_file} . $SPACE
          . catfile(qw{ dir tranchesfile.tranches }),
    },
    num_cpu_threads_per_data_thread => {
        input           => $CPU_THREADS,
        expected_output => q{--num_cpu_threads_per_data_thread } . $CPU_THREADS,
    },
    disable_indel_qual => {
        input           => 1,
        expected_output => q{--disable_indel_quals},
    },
    ts_filter_level => {
        input           => 99,
        expected_output => q{--ts_filter_level } . $TS_FILTER,
    },

);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&gatk_applyrecalibration;

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

## build_usage

## Function  : Build the USAGE instructions
## Returns   : ""
## Arguments : $program_name
##           : $program_name => Name of the script

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
