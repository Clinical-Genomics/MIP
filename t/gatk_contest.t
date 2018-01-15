#!/usr/bin/env perl

use 5.018;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname  };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use Getopt::Long;
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Modern::Perl qw{ 2014 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = 1.0.0;

## Constants
Readonly my $COMMA   => q{,};
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };

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
    my %perl_module = ( q{MIP::Script::Utils} => [qw{ help }], );

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

## Modules
    my @modules = (q{MIP::Program::Qc::Contest});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Program::Qc::Contest qw{ gatk_contest };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test gatk_contest from Contest.pm v}
      . $MIP::Program::Qc::Contest::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my $function_base_command = q{--analysis_type ContEst};

my %base_argument = (
    FILEHANDLE => {
        input           => undef,
        expected_output => $function_base_command,
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    infile_eval => {
        input           => catfile(qw{path to input_bam_eval}),
        expected_output => q{--input_file:eval}
          . $SPACE
          . catfile(qw{path to input_bam_eval}),
    },
    infile_genotype => {
        input           => catfile(qw{path to input_bam_genotype}),
        expected_output => q{--input_file:genotype}
          . $SPACE
          . catfile(qw{path to input_bam_genotype}),
    },
    outfile_path => {
        input           => catfile(qw{path to output_file}),
        expected_output => q{--out} . $SPACE . catfile(qw{path to output_file}),
    },
    pop_vcffile_path => {
        input           => catfile(qw{path to population_vcf}),
        expected_output => q{--popfile}
          . $SPACE
          . catfile(qw{path to population_vcf}),
    },
    referencefile_path => {
        input           => catfile(qw{path to human_fasta_file}),
        expected_output => q{--reference_sequence}
          . $SPACE
          . catfile(qw{path to human_fasta_file}),
    },
);

my %specific_argument = (
    infile_eval => {
        input           => catfile(qw{path to input_bam_eval}),
        expected_output => q{--input_file:eval}
          . $SPACE
          . catfile(qw{path to input_bam_eval}),
    },
    infile_genotype => {
        input           => catfile(qw{path to input_bam_genotype}),
        expected_output => q{--input_file:genotype}
          . $SPACE
          . catfile(qw{path to input_bam_genotype}),
    },
    min_genotype_ratio => {
        input           => q{0.95},
        expected_output => q{--min_genotype_ratio} . $SPACE . q{0.95},
    },
    pop_vcffile_path => {
        input           => catfile(qw{path to population_vcf}),
        expected_output => q{--popfile}
          . $SPACE
          . catfile(qw{path to population_vcf}),
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&gatk_contest;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

ARGUMENT_HASH_REF:
foreach my $argument_href (@arguments) {
    my @commands = test_function(
        {
            argument_href          => $argument_href,
            do_test_base_command   => 1,
            function_base_command  => $function_base_command,
            module_function_cref   => $module_function_cref,
            required_argument_href => \%required_argument,
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
            store       => \$program_name,
            strict_type => 1,
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
