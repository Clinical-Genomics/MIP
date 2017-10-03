#!/usr/bin/env perl

#### Copyright 2017 Henrik Stranneheim

use Modern::Perl qw(2014);
use warnings qw(FATAL utf8);
use autodie;
use 5.018;    #Require at least perl 5.18
use utf8;
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use English qw(-no_match_vars);
use Params::Check qw(check allow last_error);

use FindBin qw($Bin);    #Find directory of script
use File::Basename qw(dirname basename);
use File::Spec::Functions qw(catfile catdir devnull);
use Getopt::Long;
use Test::More;

## Third party module(s)
use List::Util qw(any);

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use MIP::Script::Utils qw(help);
use MIP::Test::Writefile qw(test_write_to_file);

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
        print {*STDOUT} "\n" . basename($PROGRAM_NAME) . q{  } . $VERSION,
          "\n\n";
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
    my @modules = ('MIP::Language::Shell');

    for my $module (@modules) {

        require_ok($module) or BAIL_OUT 'Cannot load ' . $module;
    }
}

use MIP::Language::Shell qw(build_shebang);

my $NEWLINE = q{\n};

diag(
"Test build_shebang $MIP::Language::Shell::VERSION, Perl $^V, $EXECUTABLE_NAME"
);

## Base arguments
my $batch_shebang = q{#!};

my %base_argument = (
    FILEHANDLE => {
        input           => undef,
        expected_output => $batch_shebang,
    },
);

my $bash_bin_path =
  catfile( dirname( dirname( devnull() ) ), qw(usr bin env bash) );

## Specific arguments
my %argument = (
    bash_bin_path => {
        input           => $bash_bin_path,
        expected_output => $batch_shebang . $bash_bin_path . q{ --login},
    },
    invoke_login_shell => {
        input           => 1,
        expected_output => $batch_shebang . $bash_bin_path . q{ --login},
    },
);

my @commands = build_shebang(
    {
        bash_bin_path   => $argument{bash_bin_path}{input},
        invoke_login_shell => $argument{invoke_login_shell}{input},
    }
);

## Testing return of commands
foreach my $key ( keys %argument ) {

    # Alias expeceted output
    my $expected_output = $argument{$key}{expected_output};

    ok( ( any { $_ eq $expected_output } @commands ), 'Argument: ' . $key );
}

## Testing write to file

# Fake arguments
my @args = (
    bash_bin_path => $bash_bin_path,
    FILEHANDLE    => undef,
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&build_shebang;

my $function_base_command = $batch_shebang . $bash_bin_path,;

test_write_to_file(
    {
        args_ref             => \@args,
        module_function_cref => $module_function_cref,
        base_command         => $function_base_command,
        separator            => $NEWLINE,
    }
);

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
