#!/usr/bin/env perl

###Copyright 2016 Henrik Stranneheim

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
	qq{mip.t
           -h/--help Display this help message   
           -v/--version Display version
        };    
}

my $mip_version = "0.0.0";

###User Options
GetOptions('h|help' => sub { print STDOUT $USAGE, "\n"; exit;},  #Display help text
	   'v|version' => sub { print STDOUT "\ntest.t ".$mip_version, "\n\n"; exit;},  #Display version number
);

ok(check_command_in_path({program => "mip.pl"}), "Checking can run mip.pl");

my $mip_script = catfile(dirname($Bin), "mip.pl");

## Test execution of mip.pl
my $verbose = 1;

my $cmds_ref = [$mip_script,
		"-f", "643594-miptest",
		"-c", catfile(dirname($Bin), "templates", "mip_travis_config.yaml"),
		"-ifd", catfile("data", "643594-miptest", "test_data", "ADM1059A1", "fastq=ADM1059A1"),
		"-ifd", catfile("data", "643594-miptest", "test_data", "ADM1059A2", "fastq=ADM1059A2"),
		"-ifd", catfile("data", "643594-miptest", "test_data", "ADM1059A3", "fastq=ADM1059A3"),
		"-rio", "1",
		"-dra", "2",
		"-pvep", "0",
		"-psvv", "0",
    ];  # Create array ref for cmd

my( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
    run( command => $cmds_ref, verbose => $verbose );
ok($success, "Executed mip.pl");

##Clean-up
$cmds_ref = ["rm", catfile(getcwd(), "data", "643594-miptest", "analysis")];
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
