#!/usr/bin/env perl

use Modern::Perl qw{ 2018 };
use warnings qw{ FATAL utf8 };
use autodie;
use 5.026;
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };

use FindBin qw{ $Bin };
use File::Basename qw{ dirname basename };
use File::Spec::Functions qw{ catdir };
use Getopt::Long;
use Test::More;

## CPANM
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = 1.0.1;

## Constants
Readonly my $SPACE                     => q{ };
Readonly my $NEWLINE                   => qq{\n};
Readonly my $COMMA                     => q{,};
Readonly my $SNP_GAP_FILTER_DISTANCE   => 50;
Readonly my $INDEL_GAP_FILTER_DISTANCE => 100;

###User Options
GetOptions(
    q{h|help} => sub {
        done_testing();
        say {*STDOUT} $USAGE;
        exit;
    },    #Display help text
    q{v|version} => sub {
        done_testing();
        say {*STDOUT} $NEWLINE . basename($PROGRAM_NAME) . $SPACE . $VERSION . $NEWLINE;
        exit;
    },    #Display version number
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
    my @modules = (q{MIP::Program::Variantcalling::Bcftools});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Program::Variantcalling::Bcftools qw{ bcftools_filter };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test bcftools_filter from Bcftools.pm v}
      . $MIP::Program::Variantcalling::Bcftools::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ bcftools };

my %base_argument = (
    stdoutfile_path => {
        input           => q{stdoutfile.test},
        expected_output => q{1> stdoutfile.test},
    },
    stderrfile_path => {
        input           => q{stderrfile.test},
        expected_output => q{2> stderrfile.test},
    },
    stderrfile_path_append => {
        input           => q{stderrfile.test},
        expected_output => q{2>> stderrfile.test},
    },
    FILEHANDLE => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = ();

my %specific_argument = (
    infile_path => {
        input           => q{infile.test},
        expected_output => q{infile.test},
    },
    exclude => {
        input           => q{%QUAL<10 || (RPB<0.1 && %QUAL<15)},
        expected_output => q{--exclude} . $SPACE . q{%QUAL<10 || (RPB<0.1 && %QUAL<15)},
    },
    filter_mode => {
        input           => q{+},
        expected_output => q{--mode} . $SPACE . q{+},
    },
    include => {
        input           => q{INFO/CSQ[*]~":p[.]"},
        expected_output => q{--include} . $SPACE . q{INFO/CSQ[*]~":p[.]"},
    },
    soft_filter => {
        input           => q{LowQual},
        expected_output => q{--soft-filter} . $SPACE . q{LowQual},
    },
    snp_gap => {
        input           => $SNP_GAP_FILTER_DISTANCE,
        expected_output => q{--SnpGap} . $SPACE . $SNP_GAP_FILTER_DISTANCE,
    },
    indel_gap => {
        input           => $INDEL_GAP_FILTER_DISTANCE,
        expected_output => q{--IndelGap} . $SPACE . $INDEL_GAP_FILTER_DISTANCE,
    },

);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&bcftools_filter;

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
