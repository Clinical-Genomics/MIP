#!/usr/bin/env perl

#### Copyright 2017 Henrik Stranneheim

use Modern::Perl qw{ 2018 };
use warnings qw(FATAL utf8);
use autodie;
use 5.026;    #Require at least perl 5.18
use utf8;
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use English qw(-no_match_vars);
use Params::Check qw(check allow last_error);

use FindBin qw($Bin);    #Find directory of script
use File::Basename qw(dirname basename);
use File::Spec::Functions qw(catdir);
use Getopt::Long;
use Test::More;

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use MIP::Script::Utils qw(help);

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '0.0.0';

###User Options
GetOptions(
    'h|help' => sub {
        done_testing();
        print {*STDOUT} $USAGE, "\n";
        exit;
    },    #Display help text
    'v|version' => sub {
        done_testing();
        print {*STDOUT} "\n" . basename($PROGRAM_NAME) . q{  } . $VERSION, "\n\n";
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
##Modules with import
    my %perl_module;

    $perl_module{'MIP::Script::Utils'} = [qw(help)];
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT 'Cannot load ' . $module;
    }

##Modules
    my @modules = ('MIP::Gnu::Coreutils');
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT 'Cannot load ' . $module;
    }
}

use MIP::Gnu::Coreutils qw(gnu_mkdir);
use MIP::Test::Commands qw(test_function);

diag("Test gnu_mkdir $MIP::Gnu::Coreutils::VERSION, Perl $^V, $EXECUTABLE_NAME");

## Base arguments
my @function_base_commands = 'mkdir';

my %base_argument = (
    stderrfile_path => {
        input           => 'stderrfile.test',
        expected_output => '2> stderrfile.test',
    },
    stderrfile_path_append => {
        input           => 'stderrfile.test',
        expected_output => '2>> stderrfile.test',
    },
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
);

## Can be duplicated with %base and/or %specific to enable testing of each individual argument
my %required_argument = (
    indirectory_path => {
        input           => '/test/path/',
        expected_output => '/test/path/',
    },
);

## Specific arguments
my %specific_argument = (
    indirectory_path => {
        input           => '/test/path/',
        expected_output => '/test/path/',
    },
    parents => {
        input           => 1,
        expected_output => '--parents',
    },
    verbose => {
        input           => 1,
        expected_output => '--verbose',
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&gnu_mkdir;

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

    check( $tmpl, $arg_href, 1 ) or croak qw(Could not parse arguments!);

    return <<"END_USAGE";
 $program_name [options]
    -vb/--verbose Verbose
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}
