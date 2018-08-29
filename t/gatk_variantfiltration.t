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
use List::Util qw{ any };
use Test::More;
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = 1.0.0;

## Constants
Readonly my $COMMA        => q{,};
Readonly my $DOUBLE_QOUTE => q{"};
Readonly my $NEWLINE      => qq{\n};
Readonly my $SPACE        => q{ };

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

use MIP::Program::Variantcalling::Gatk qw{ gatk_variantfiltration };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test gatk_variantfiltration from Varinatcalling::Gatk.pm v}
      . $MIP::Program::Variantcalling::Gatk::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $CLUSTER_SIZE        => 3;
Readonly my $CLUSTER_WINDOW_SIZE => 35;

## Base arguments
my @function_base_commands = qw{ --analysis_type VariantFiltration };

my %base_argument = (
    stderrfile_path => {
        input           => q{stderrfile.test},
        expected_output => q{2> stderrfile.test},
    },
    FILEHANDLE => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    infile_path => {
        input           => catfile(qw{ dir infile.vcf }),
        expected_output => q{--variant } . catfile(qw{ dir infile.vcf }),
    },
    outfile_path => {
        input           => catfile(qw{ dir outfile.vcf }),
        expected_output => q{--out } . catfile(qw{ dir outfile.vcf }),
    },
    referencefile_path => {
        input           => catfile(qw{reference_dir human_genome_build.fasta }),
        expected_output => q{--reference_sequence}
          . catfile(qw{reference_dir human_genome_build.fasta }),
    },
);

my %specific_argument = (
    cluster_size => {
        input           => $CLUSTER_SIZE,
        expected_output => q{--clusterSize} . $SPACE . $CLUSTER_SIZE,
    },
    infile_path => {
        input           => catfile(qw{ dir infile.vcf }),
        expected_output => q{--variant } . catfile(qw{ dir infile.vcf }),
    },
    outfile_path => {
        input           => catfile(qw{ dir outfile.vcf }),
        expected_output => q{--out } . catfile(qw{ dir outfile.vcf }),
    },
    referencefile_path => {
        input           => catfile(qw{reference_dir human_genome_build.fasta }),
        expected_output => q{--reference_sequence }
          . catfile(qw{ reference_dir human_genome_build.fasta }),
    },
    cluster_window_size => {
        input           => $CLUSTER_WINDOW_SIZE,
        expected_output => q{--clusterWindowSize}
          . $SPACE
          . $CLUSTER_WINDOW_SIZE,
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&gatk_variantfiltration;

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

## Special case for hash input.
## Given a hash with filter names and filter expressions
my %filter = ( filter_1 => q{test_1 > 1}, );

## When the subroutine is executed
my @commands = gatk_variantfiltration(
    {
        infile_path        => catfile(qw{ a infile_path }),
        outfile_path       => catfile(qw{ a outfile_path }),
        referencefile_path => catfile(qw{ a genome }),
        filter_href        => \%filter,
    }
);

## Then the string below should exactly match one of the strings in the command array.
my $expected_output =
    q{--filterName filter_1}
  . $SPACE
  . q{--filterExpression}
  . $SPACE
  . $DOUBLE_QOUTE
  . $filter{filter_1}
  . $DOUBLE_QOUTE;
my $filter_match = grep { $_ eq $expected_output } @commands;
ok( $filter_match, q{Argument: filter_href} );

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
