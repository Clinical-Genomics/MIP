#!/usr/bin/env perl

use Modern::Perl qw{ 2018 };
use warnings qw{FATAL utf8};
use autodie;
use 5.026;    #Require at least perl 5.18
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

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '1.0.0';

## Constants
Readonly my $SPACE        => q{ };
Readonly my $TAB          => q{\t};
Readonly my $DOUBLE_QUOTE => q{"};
Readonly my $NEWLINE      => qq{\n};

###User Options
GetOptions(
    q{h|help} => sub {
        done_testing();
        say {*STDOUT} $USAGE;
        exit;
    },    #Display help text
    q{v|version} => sub {
        done_testing();
        say {*STDOUT} $NEWLINE, basename($PROGRAM_NAME), $SPACE, $VERSION, $NEWLINE;
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
##Modules with import
    my %perl_module;

    $perl_module{'MIP::Script::Utils'} = [qw{help}];

  PERL_MODULES:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load } . $module;
    }

## Modules
    my @modules = (q{MIP::Program::Alignment::Bwa});

  MODULES:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load } . $module;
    }
}

use MIP::Program::Alignment::Bwa qw{run_bwamem};
use MIP::Test::Commands qw{test_function};

diag(
    q{Test run_bwamem }
      . $MIP::Program::Alignment::Bwa::VERSION
      . q{, Perl}
      . $PERL_VERSION,
    $EXECUTABLE_NAME
);

## Base arguments
my @function_base_commands = qw{ run-bwamem };

## Read group header line
my @read_group_headers = (
    $DOUBLE_QUOTE . q{@RG},
    q{ID:} . q{1_140128_H8AHNADXX_1-1-2A_GATCAG_1} . $TAB,
    q{SM:} . q{1-1-2A},
    q{PL:} . q{ILLUMINA} . $DOUBLE_QUOTE,
);

my %base_argument = (
    stderrfile_path => {
        input           => q{stderrfile.test},
        expected_output => q{2> stderrfile.test},
    },
    stderrfile_path_append => {
        input           => q{stderrfile.test},
        expected_output => q{2>> stderrfile.test},
    },
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

## Can be duplicated with %base and/or %specific to enable testing of each individual argument
my %required_argument = (
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
    infile_path => {
        input           => q{test_infile.fastq},
        expected_output => q{test_infile.fastq},
    },
    idxbase => {
        input           => q{grch37_homo_sapiens_-d5-.fasta},
        expected_output => q{grch37_homo_sapiens_-d5-.fasta},
    },
    outfiles_prefix_path => {
        input           => q{test_outfile.bam},
        expected_output => q{-o test_outfile.bam},
    },
);

my %specific_argument = (
    thread_number => {
        input           => 2,
        expected_output => q{-t 2},
    },
    hla_typing => {
        input           => 1,
        expected_output => q{-H},
    },
    read_group_header => {
        input           => ( join $SPACE,                  @read_group_headers ),
        expected_output => ( q{-R} . $SPACE . join $SPACE, @read_group_headers ),
    },
    infile_path => {
        input           => q{test_infile_1.fastq},
        expected_output => q{test_infile_1.fastq},
    },
    second_infile_path => {
        input           => q{test_infile_2.fastq},
        expected_output => q{test_infile_2.fastq},
    },
    outfiles_prefix_path => {
        input           => q{test_outfile.bam},
        expected_output => q{-o test_outfile.bam},
    },
    stderrfile_path => {
        input           => q{stderrfile.test},
        expected_output => q{2> stderrfile.test},
    },
    stderrfile_path_append => {
        input           => q{stderrfile.test},
        expected_output => q{2>> stderrfile.test},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&run_bwamem;

## Test both base and function specific arguments
my @arguments = ( \%required_argument, \%specific_argument );

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

    check( $tmpl, $arg_href, 1 ) or croak qw(Could not parse arguments!);

    return <<"END_USAGE";
 $program_name [options]
    -vb/--verbose Verbose
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}
