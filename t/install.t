#!/usr/bin/env perl

###Copyright 2016 Henrik Stranneheim'

BEGIN {

    ## Special case to initiate testing
    my @modules = ("Test::More",
	);

    ## Evaluate that all modules required are installed
    &EvalModules(\@modules);
    
    sub EvalModules {
	
	##EvalModules
	
	##Function : Evaluate that all modules required are installed 
	##Returns  : ""
	##Arguments: $modules_ref
	##         : $modules_ref => Array of module names
	
	my $modules_ref = $_[0];
	
	foreach my $module (@{$modules_ref}) {
	    
	    $module =~s/::/\//g;  #Replace "::" with "/" since the automatic replacement magic only occurs for barewords.
	    $module .= ".pm";  #Add perl module ending for the same reason
	    
	    eval { 
		
		require $module; 
	    };
	    if($@) {
		
		warn("NOTE: ".$module." not installed - Please install to run mip tests.\n");
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

    $perl_module{charnames} = qw( :full :short );
    $perl_module{File::Basename} = qw(dirname);
    $perl_module{File::Spec::Functions} =  qw(catfile catdir devnull);
    $perl_module{FindBin} = qw($Bin);
    $perl_module{IPC::Cmd} = qw[can_run run];
    $perl_module{open} = qw( :encoding(UTF-8) :std );
    $perl_module{Params::Check} = qw[check allow last_error];
    $perl_module{warnings} = qw( FATAL utf8 );
    
    while (my ($module, $module_import) = each %perl_module) { 

	use_ok($module, $module_import) or BAIL_OUT "Can't load $module";
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
use File::Basename qw(dirname);
use File::Spec::Functions qw(catfile catdir devnull);
use Params::Check qw[check allow last_error];
use IPC::Cmd qw[can_run run];
use Getopt::Long;
use Cwd;

our $USAGE;

BEGIN {
    $USAGE =
	qq{install.t
           -h/--help Display this help message   
           -v/--version Display version
        };    
}

my $install_version = "0.0.0";

###User Options
GetOptions('h|help' => sub { print STDOUT $USAGE, "\n"; exit;},  #Display help text
	   'v|version' => sub { print STDOUT "\ntest.t ".$install_version, "\n\n"; exit;},  #Display version number
);

ok(check_command_in_path({program => "conda"}), "Checking execution of conda binary");

ok(catfile(dirname($Bin), "install.pl"), "Locating install script in MIP dir");

my $install_script = catfile(dirname($Bin), "install.pl");

## Test execution of install.pl
my $verbose = 0;

my $cmds_ref = ["perl", $install_script, "-sp", "mip_scripts"];  # Create array ref for cmd
my( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
    run( command => $cmds_ref, verbose => $verbose );
ok($success, "Execute install.pl");

is(-e catfile(getcwd(), "mip.sh"), 1, "Locating created mip.sh in MIP dir");

##Clean-up
$cmds_ref = ["rm", catfile(getcwd(), "mip.sh")];
( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
    run( command => $cmds_ref, verbose => $verbose );

done_testing();


######################
####SubRoutines#######
######################


sub check_command_in_path {

##check_command_in_path

##Function : Checking command in your path and executable
##Returns  : ""
##Arguments: $program
##         : $program => Program to check

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $program;

    my $tmpl = { 
	program => { required => 1, defined => 1, strict_type => 1, store => \$program},
    };
        
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    if(can_run($program)) {  #IPC::Cmd

	return 1;
    }
    else {
	
	return 0;
    }
}



