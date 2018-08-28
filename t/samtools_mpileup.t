#!/usr/bin/env perl

#### Copyright 2017 Henrik Stranneheim

use Modern::Perl qw{2014};
use warnings qw{FATAL utf8};
use autodie;
use 5.018;    #Require at least perl 5.18
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{-no_match_vars};
use Params::Check qw{check allow last_error};

use FindBin qw{$Bin};    #Find directory of script
use File::Basename qw{dirname basename};
use File::Spec::Functions qw{catdir};
use Getopt::Long;
use Test::More;
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use MIP::Script::Utils qw{help};

## Constants
Readonly my $SPACE   => q{ };
Readonly my $NEWLINE => qq{\n};

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '1.0.0';

###User Options
GetOptions(
    'h|help' => sub {
        done_testing();
        say {*STDOUT} $USAGE;
        exit;
    },    #Display help text
    'v|version' => sub {
        done_testing();
        say {*STDOUT} basename($PROGRAM_NAME) . $SPACE . $VERSION, $NEWLINE;

        exit;
    },    #Display version number
    'vb|verbose' => $VERBOSE,
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

    my %perl_module;

    $perl_module{'MIP::Script::Utils'} = [qw{help}];

  PERL_MODULES:
    while ( my ( $module, $module_import ) = each %perl_module ) {

        use_ok( $module, @{$module_import} )
          or BAIL_OUT 'Cannot load ' . $module;
    }

    my @modules = ('MIP::Program::Alignment::Samtools');

  MODULES:
    for my $module (@modules) {

        require_ok($module) or BAIL_OUT 'Cannot load ' . $module;
    }
}

use MIP::Program::Alignment::Samtools qw{samtools_mpileup};
use MIP::Test::Commands qw{test_function};

diag(
"Test samtools_mpileup $MIP::Program::Alignment::Samtools::VERSION, Perl $^V, $EXECUTABLE_NAME"
);

## Base arguments
my @function_base_commands = qw{ samtools };

my %base_argument = (
    FILEHANDLE => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

## Can be duplicated with %base and/or %specific to enable testing of each individual argument
my %required_argument = (
    FILEHANDLE => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
    infile_paths_ref => {
        inputs_ref      => [],
        expected_output => q{},
    },
    output_tags_ref => {
        inputs_ref      => [qw{tag1 tag3 etc}],
        expected_output => q{--output-tags tag1,tag2,etc},
    },
    referencefile_path => {
        input           => q{fastaref},
        expected_output => q{--fasta-ref fastaref},
    },
);

## Specific arguments
my %specific_argument = (

    outfile_path => {
        input           => q{outpath},
        expected_output => q{--output outpath},
    },
    stderrfile_path => {
        input           => q{stderrfile.test},
        expected_output => q{2> stderrfile.test},
    },
    region => {
        input           => q{3:45000-67000},
        expected_output => q{--region 3:45000-67000},
    },
    output_bcf => {
        input           => q{1},
        expected_output => q{--BCF},
    },
    per_sample_increased_sensitivity => {
        input           => q{1},
        expected_output => q{--per-sample-mF},
    },
    adjust_mq => {
        input           => q{45},
        expected_output => q{--adjust-MQ 45},
    },
    stderrfile_path_append => {
        input           => q{stderrfile_path_append},
        expected_output => q{2>> stderrfile_path_append},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&samtools_mpileup;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

foreach my $argument_href (@arguments) {

    my @commands = test_function(
        {
            argument_href              => $argument_href,
            required_argument_href     => \%required_argument,
            module_function_cref       => $module_function_cref,
            function_base_commands_ref => \@function_base_commands,
        }
    );
}

done_testing();

######################
####SubRoutines#######
######################

sub build_usage {

##build_usage

##Function : Build the USAGE instructions
##Returns  : ""
##Arguments: $program_name
##         : $program_name => Name of the script

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

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    return <<"END_USAGE";
 $program_name [options]
    -vb/--verbose Verbose
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}
