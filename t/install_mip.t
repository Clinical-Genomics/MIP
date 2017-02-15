#!/usr/bin/env perl

###Copyright 2016 Henrik Stranneheim'

BEGIN {

    ## Special case to initiate testing
    my @modules = ("Test::More",
	);

    ## Evaluate that all modules required are installed
    &check_modules(\@modules);
    
    sub check_modules {
	
	##check_modules
	
	##Function : Evaluate that all modules required are installed 
	##Returns  : ""
	##Arguments: $modules_ref
	##         : $modules_ref => Array of module names
	
	my $modules_ref = $_[0];
	
	foreach my $module (@$modules_ref) {
	    
	    $module =~s/::/\//g;  #Replace "::" with "/" since the automatic replacement magic only occurs for barewords.
	    $module .= ".pm";  #Add perl module ending for the same reason
	    
	    eval { 
		
		require $module; 
	    };
	    if($@) {
		
		warn("NOTE: ".$module." not installed - Please install to run install.t.\n");
		warn("NOTE: Aborting!\n");
		exit 1;
	    }
	}
    }

    ##Initate tests
    print STDOUT "Initiate tests:\n";
    print STDOUT "\nTesting default perl modules and selected functions\n\n";


    ##More proper testing 
    use Test::More;

    ##Modules with import
    my %perl_module;

    $perl_module{charnames} = [qw(:full :short)];
    $perl_module{"File::Basename"} = [qw(dirname basename)];
    $perl_module{"File::Spec::Functions"} =  [qw(catfile catdir devnull)];
    $perl_module{FindBin} = [qw($Bin)];
    $perl_module{"IPC::Cmd"} = [qw(can_run run)];
    $perl_module{open} = [qw(:encoding(UTF-8) :std )];
    $perl_module{"Params::Check"} = [qw(check allow last_error)];
    $perl_module{warnings} = [qw(FATAL utf8 )];
    
    while (my ($module, $module_import) = each %perl_module) { 

	use_ok($module, @$module_import) or BAIL_OUT "Can't load $module";
    }
    
    ##Modules
    @modules = ("Cwd",
		"Getopt::Long",
		"utf8",
		"strict",
		"warnings",
	);

    for my $module (@modules) {

	require_ok($module) or BAIL_OUT "Can't load $module";
    }
}


use strict;
use warnings;
use warnings qw( FATAL utf8 );
use utf8;
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Test::More;
use FindBin qw($Bin);
use File::Basename qw(dirname basename);
use File::Spec::Functions qw(catfile catdir devnull);
use Params::Check qw[check allow last_error];
use IPC::Cmd qw[can_run run];
use Getopt::Long;
use Cwd;

##MIPs lib/
use lib catdir(dirname($Bin), "lib");
use Script::Utils qw(help);

our $USAGE;

BEGIN {

    $USAGE =
	basename($0).qq{
           -vb/--verbose Verbose
           -h/--help Display this help message   
           -v/--version Display version
        };    
}

my $verbose = 1;
my $install_version = "0.0.0";

###User Options
GetOptions('vb|verbose' => $verbose,
	   'h|help' => sub { done_testing(); print STDOUT $USAGE, "\n"; exit;},  #Display help text
	   'v|version' => sub { done_testing(); print STDOUT "\n".basename($0)." ".$install_version, "\n\n"; exit;},  #Display version number
    ) or done_testing(), Script::Utils::help({USAGE => $USAGE,
					      exit_code => 1,
					     });;

ok(can_run("conda")), "Checking can run of conda binary");

ok(catfile(dirname($Bin), "install_mip.pl"), "Locating install script in MIP dir");

my $install_script = catfile(dirname($Bin), "install_mip.pl");

## Test execution of install.pl
my $cmds_ref = ["perl", $install_script, "-sp", "mip_scripts"];  # Create array ref for cmd
my( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
    run( command => $cmds_ref, verbose => $verbose );
ok($success, "Executed install.pl");

is(-e catfile(getcwd(), "mip.sh"), 1, "Locating created mip.sh in MIP dir");

##Clean-up
$cmds_ref = ["rm", catfile(getcwd(), "mip.sh")];
( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
    run( command => $cmds_ref, verbose => $verbose );

done_testing();


######################
####SubRoutines#######
######################



