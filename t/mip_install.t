#!/usr/bin/env perl

### Copyright 2016 Henrik Stranneheim

use strict;
use warnings;
use warnings qw( FATAL utf8 );
use utf8;
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use English qw(-no_match_vars);

BEGIN {

    ## Special case to initiate testing
    my @modules = qw(Test::More);

    ## Evaluate that all modules required are installed
    check_modules( \@modules );

    sub check_modules {

        ##check_modules

        ##Function : Evaluate that all modules required are installed
        ##Returns  : ""
        ##Arguments: $modules_ref
        ##         : $modules_ref => Array of module names

        my ($modules_ref) = @_;

        foreach my $module ( @{$modules_ref} ) {

# Replace "::" with "/" since the automatic replacement magic only occurs for barewords.
            $module =~ s/::/\//gxms;

            # Add perl module ending for the same reason
            $module .= qw(.pm);

            if ( eval { require $module; 1 } ) {
            }
            if ($EVAL_ERROR) {

                warn qw(NOTE: ) . $module
                  . " not installed - Please install to run install.t.\n";
                warn "NOTE: Aborting!\n";
                exit 1;
            }
        }
        return;
    }

    ## Initate tests
    print {*STDOUT} "Initiate tests:\n";
    print {*STDOUT} "\nTesting default perl modules and selected functions\n\n";

    ## More proper testing
    use Test::More;

    ## Modules with import
    my %perl_module;

    $perl_module{charnames}               = [qw(:full :short)];
    $perl_module{'File::Basename'}        = [qw(dirname basename)];
    $perl_module{'File::Spec::Functions'} = [qw(catfile catdir devnull)];
    $perl_module{FindBin}                 = [qw($Bin)];
    $perl_module{'IPC::Cmd'}              = [qw(can_run run)];
    $perl_module{open}                    = [qw(:encoding(UTF-8) :std )];
    $perl_module{'Params::Check'}         = [qw(check allow last_error)];
    $perl_module{warnings}                = [qw(FATAL utf8 )];

    while ( my ( $module, $module_imports_ref ) = each %perl_module ) {

        use_ok( $module, @{$module_imports_ref} )
          or BAIL_OUT "Can't load $module";
    }

    ## Modules
    @modules = qw(Cwd Getopt::Long utf8 strict warnings);

    for my $module (@modules) {

        require_ok($module) or BAIL_OUT "Can't load $module";
    }
}

use Test::More;
use FindBin qw( $Bin );
use File::Basename qw( dirname basename );
use File::Spec::Functions qw( catfile catdir devnull );
use Params::Check qw[ check allow last_error ];
use IPC::Cmd qw[ can_run run ];
use Getopt::Long;
use Cwd;

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use Script::Utils qw( help );

our $USAGE = build_usage( {} );

my $conda_dir_path;
my $verbose = 1;
our $VERSION = '0.0.1';

###User Options
GetOptions(
    'cdp|conda_dir_path:s' => \$conda_dir_path,
    'vb|verbose'           => $verbose,
    'h|help'               => sub {
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

ok( can_run('conda'), 'Checking can run of conda binary' );

ok(
    catfile( dirname($Bin), 'mip_install.pl' ),
    'Locating install script in MIP dir'
);

my $install_script = catfile( dirname($Bin), 'mip_install.pl' );

## Test execution of install.pl
# Create array ref for cmd
my $cmds_ref = [ qw(perl), $install_script, qw(-sp mip_scripts) ];
if ($conda_dir_path) {

    push @{$cmds_ref}, '--conda_dir_path', $conda_dir_path;
}

my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run(
    command => $cmds_ref,
    verbose => $verbose
);

ok( $success, 'Executed install.pl' );

is( -e catfile( getcwd(), 'mip.sh' ), 1, 'Locating created mip.sh in MIP dir' );

##Clean-up
$cmds_ref = [ 'rm', catfile( getcwd(), 'mip.sh' ) ];
( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run(
    command => $cmds_ref,
    verbose => $verbose
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

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    return <<"END_USAGE";
 $program_name [options]
    -cdp/--conda_dir_path The conda directory path (Default: "HOME/miniconda")
    -vb/--verbose Verbose
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}

