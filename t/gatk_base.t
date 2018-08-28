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
    my @modules = (q{MIP::Program::Base::Gatk});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Program::Base::Gatk qw{ gatk_base };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test gatk_base from Base::Gatk.pm v}
      . $MIP::Program::Base::Gatk::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $DOWNSAMPLE_TO_COVERAGE               => 1000;
Readonly my $STATIC_QUANTIZED_QUALS_MINIMUM_LEVEL => 10;
Readonly my $STATIC_QUANTIZED_QUALS_MAX_LEVEL     => 20;

## Base arguments
my @function_base_commands = qw{ --analysis_type HaplotypeCaller };

my %base_argument = (
    FILEHANDLE => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    analysis_type => {
        input           => q{HaplotypeCaller},
        expected_output => q{--analysis_type HaplotypeCaller},
    },
    referencefile_path => {
        input           => catfile(qw{reference_dir human_genome_build.fasta }),
        expected_output => q{--reference_sequence}
          . catfile(qw{reference_dir human_genome_build.fasta }),
    },
);

my %specific_argument = (
    intervals_ref => {
        inputs_ref      => [qw{ chr1 chr2}],
        expected_output => q{--intervals chr1 --intervals chr2},
    },
    read_filters_ref => {
        inputs_ref => [qw{ MalformedRead BadCigar}],
        expected_output =>
          q{--read_filter MalformedRead --read_filter BadCigar},
    },
    analysis_type => {
        input           => q{HaplotypeCaller},
        expected_output => q{--analysis_type HaplotypeCaller},
    },
    logging_level => {
        input           => q{INFO},
        expected_output => q{--logging_level INFO},
    },
    pedigree_validation_type => {
        input           => q{SILENT},
        expected_output => q{--pedigreeValidationType SILENT},
    },
    pedigree => {
        input           => catfile(qw{ dir ped.fam }),
        expected_output => q{--pedigree } . catfile(qw{ dir ped.fam }),
    },
    num_cpu_threads_per_data_thread => {
        input           => 2,
        expected_output => q{--num_cpu_threads_per_data_thread 2},
    },
    downsample_to_coverage => {
        input           => $DOWNSAMPLE_TO_COVERAGE,
        expected_output => q{--downsample_to_coverage }
          . $DOWNSAMPLE_TO_COVERAGE,
    },
    gatk_disable_auto_index_and_file_lock => {
        input => 1,
        expected_output =>
          q{--disable_auto_index_creation_and_locking_when_reading_rods},
    },
    referencefile_path => {
        input           => catfile(qw{reference_dir human_genome_build.fasta }),
        expected_output => q{--reference_sequence }
          . catfile(qw{reference_dir human_genome_build.fasta }),
    },
    base_quality_score_recalibration_file => {
        input           => catfile(qw{ dir infile.bsqr }),
        expected_output => q{--BQSR } . catfile(qw{ dir infile.bsqr }),
    },
    disable_indel_qual => {
        input           => 1,
        expected_output => q{--disable_indel_quals},
    },
    static_quantized_quals_ref => {
        inputs_ref => [
            $STATIC_QUANTIZED_QUALS_MINIMUM_LEVEL,
            $STATIC_QUANTIZED_QUALS_MAX_LEVEL
        ],
        expected_output => q{--static_quantized_quals }
          . $STATIC_QUANTIZED_QUALS_MINIMUM_LEVEL
          . q{ --static_quantized_quals }
          . $STATIC_QUANTIZED_QUALS_MAX_LEVEL,
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&gatk_base;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

ARGUMENT_HASH_REF:
foreach my $argument_href (@arguments) {
    my @commands = test_function(
        {
            argument_href              => $argument_href,
            required_argument_href     => \%required_argument,
            module_function_cref       => $module_function_cref,
            function_base_commands_ref => \@function_base_commands,
            do_test_base_command       => 1,
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
