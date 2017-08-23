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
use FileHandle;
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use Script::Utils qw{help};

## Constants
Readonly my $SPACE => q{ };
Readonly my $NEWLINE => qq{\n};

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '0.0.0';

###User Options
GetOptions(
    'h|help' => sub {
        done_testing();
        say {*STDOUT} $USAGE;
        exit;
    },    #Display help text
    'v|version' => sub {
        done_testing();
        say {*STDOUT} basename($PROGRAM_NAME) . $SPACE . $VERSION,
          $NEWLINE;

        exit;
    },    #Display version number
    'vb|verbose' => $VERBOSE,
  )
  or (
    done_testing(),
    Script::Utils::help(
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

    $perl_module{'Script::Utils'} = [qw{help}];
    
    #PERL MODULES
    while ( my ( $module, $module_import ) = each %perl_module ) {

        use_ok( $module, @{$module_import} )
          or BAIL_OUT 'Cannot load ' . $module;
    }

## Modules
    my @modules = ('MIP::Program::Alignment::Sambamba_view');

    #MODULES
    for my $module (@modules) {

        require_ok($module) or BAIL_OUT 'Cannot load ' . $module;
    }
}

use MIP::Program::Alignment::Sambamba_view qw{sambamba_view};
use MIP::Test::Commands qw{test_function};

diag(
"Test sambamba_view MIP::Program::Alignment::Sambamba_view::VERSION, Perl $^V, $EXECUTABLE_NAME"
);

## Default(s)
    my $with_header;
    my $show_progress;
    my $output_format;

    ## Flatten argument(s)
    my $regions_ref;
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $referencefile_path;
    my $stderrfile_path_append;


## Base arguments
my $function_base_command = q{sambamba};

my %base_argument = (
    FILEHANDLE => {
        input           => undef,
        expected_output => $function_base_command,
    },
);

## Can be duplicated with %base and/or %specific to enable testing of each individual argument
my %required_argument = (
    FILEHANDLE => {
        input           => undef,
        expected_output => $function_base_command,
    },
    infile_path => {
        input           => q{infile.test},
        expected_output => q{infile.test},
    },
);

## Specific arguments
my %specific_argument = (
    stderrfile_path => {
        input           => q{stderrfile.test},
        expected_output => q{2> stderrfile.test},
    },
    with_header => {
		input           => 1,
		expected_output => q{--with-header},
    },
    show_progress => {
		input           => 1,
		expected_output => q{--show-progress},
    },
    output_format => {
		input   => q{bam},
		expected_output => q{--format bam},
    },
    referencefile_path => {
        input           => q{pathToRef.test},
        expected_output => q{--ref-filename=pathToRef.test},
    
    },
    regions_ref => {
		inputs_ref           => [qw(1:1000000-2000000 2:1000-5000)],
		expected_output => q{1:1000000-2000000 2:1000-5000},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&sambamba_view;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

foreach my $argument_href (@arguments) {
	
    my @commands = test_function(
        {
            argument_href          => $argument_href,
            required_argument_href => \%required_argument,
            module_function_cref   => $module_function_cref,
            function_base_command  => $function_base_command,
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
