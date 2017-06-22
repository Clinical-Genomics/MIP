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

use Cwd;
use FindBin qw($Bin);    #Find directory of script
use File::Basename qw(dirname basename);
use File::Spec::Functions qw(catfile catdir devnull);
use Getopt::Long;
use Test::More;
use IPC::Cmd qw(can_run run);

## Third party module(s)
use List::Util qw(any);

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use Script::Utils qw(help);
use MIP::Gnu::Coreutils qw(gnu_mkdir gnu_rm);

our $USAGE = build_usage( {} );

my $VERBOSE = 0;
our $VERSION = '1.0.0';

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

    $perl_module{'Script::Utils'}           = [qw(help)];
    $perl_module{'MIP::Gnu::Coreutils'} = [qw(gnu_mkdir gnu_rm)];

    while ( my ( $module, $module_import ) = each %perl_module ) {

        use_ok( $module, @{$module_import} )
          or BAIL_OUT 'Cannot load ' . $module;
    }

    ## Modules
    my @modules = ('MIP::Language::Shell');

    for my $module (@modules) {

        require_ok($module) or BAIL_OUT 'Cannot load ' . $module;
    }
}

use MIP::Language::Shell qw(build_shebang enable_trap create_error_trap_function);
use MIP::Gnu::Bash qw(gnu_set);

my $NEWLINE = q{\n};

diag(
"Test create_error_trap_function $MIP::Language::Shell::VERSION, Perl $^V, $EXECUTABLE_NAME"
);

# Create anonymous filehandle
my $FILEHANDLE = IO::Handle->new();

# Shell create error trap file
my $bash_file_path = catfile( cwd(), 'test_create_error_trap_function.sh' );

# Temporary directory
my $temp_dir = catdir( cwd(), qw(test_dir .test_create_error_trap_function));

# Open filehandle
open $FILEHANDLE, '>', $bash_file_path
  or
  croak( q{Cannot write to '} . $bash_file_path . q{' :} . $OS_ERROR . "\n" );

## Write to bash file
_build_test_file_recipe(
    {
        FILEHANDLE     => $FILEHANDLE,
        temp_dir       => $temp_dir,
        bash_file_path => $bash_file_path,
    }
);
close $FILEHANDLE;

## Testing write to file
ok( -e $bash_file_path, 'Create bash' );

ok( can_run('bash'), 'Checking can run bash binary' );

my $cmds_ref = [ 'bash', $bash_file_path ];
my ( $success, $error_message, $full_buf_ref, $stdout_buf_ref, $stderr_buf_ref ) =
  run( command => $cmds_ref, verbose => $VERBOSE );

## Testing error trap function
ok($stderr_buf_ref->[1] =~/[:] Unknown Error - ExitCode[=]/, q{Performed error trap} );

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

sub _build_test_file_recipe {

##_build_test_file_recipe

##Function : Builds the test file for testing the housekeeping function
##Returns  : ""
##Arguments: $FILEHANDLE, $temp_dir, $bash_file_path
##         : $FILEHANDLE     => FILEHANDLE to write to
##         : $temp_dir       => Temporary directory to use for test
##         : $bash_file_path => Test file to write recipe to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $temp_dir;
    my $bash_file_path;

    my $tmpl = {
        FILEHANDLE     => { required => 1, store => \$FILEHANDLE },
        temp_dir       => { required => 1, store => \$temp_dir },
        bash_file_path => { required => 1, store => \$bash_file_path },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw(Could not parse arguments!);

    # Add bash shebang
    build_shebang(
        {
            FILEHANDLE  => $FILEHANDLE,
        }
    );

    ## Set shell attributes
    gnu_set({FILEHANDLE         => $FILEHANDLE,
	     set_errexit        => 1,
	     set_nounset        => 1,
	    });

    enable_trap(
        {
            FILEHANDLE         => $FILEHANDLE,
            trap_signals_ref   => ['DEBUG'],
            trap_function_call => q{previous_command="$BASH_COMMAND"},
        }
    );

    # Create housekeeping fucntion to remove temp_dir
    create_error_trap_function(
        {
            trap_function_name => 'error',
            FILEHANDLE         => $FILEHANDLE,
        }
    );

    # Remove batch file to make clean exit
    gnu_rm(
        {
            infile_path => $bash_file_path,
            FILEHANDLE  => $FILEHANDLE,
        }
    );

    # Create dir to test removal later
    gnu_mkdir(
        {
            indirectory_path => $temp_dir,
            parents          => 0,
            FILEHANDLE       => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} "\n";

    return;
}
