#!/usr/bin/env perl

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

ok(check_command_in_path({program => "conda"}), "Checking execution of conda binary");

ok(catfile(dirname($Bin), "install.pl"), "Locating install script in MIP dir");

my $install_script = catfile(dirname($Bin), "install.pl");

## Test execution of install.pl
my $verbose = 0;

my $cmds_ref = ["perl", $install_script, "-sp", "mip_scripts"];  # Create array ref for cmd
my( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
    run( command => $cmds_ref, verbose => $verbose );
ok($success, "Execute install.pl");

is(-e catfile(dirname($Bin), "t", "mip.sh"), 1, "Locating created mip.sh in MIP dir");

##Clean-up
$cmds_ref = ["rm", catfile(dirname($Bin), "t", "mip.sh")];
( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
    run( command => $cmds_ref, verbose => $verbose );
#ok(check_command_in_path({program => "pip"}), "Checking execution of pip binary");
#ok(check_command_in_path({program => "cpanm"}), "Checking execution of cpanm binary");

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



