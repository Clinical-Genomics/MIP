#!/usr/bin/env perl

###Copyright 2016 Henrik Stranneheim

use Modern::Perl '2014';
use warnings qw( FATAL utf8 );
use autodie;
use v5.18;  #Require at least perl 5.18
use utf8;
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );

use FindBin qw($Bin);  #Find directory of script
use File::Basename qw(dirname basename);
use File::Spec::Functions qw(catdir catfile devnull);
use Getopt::Long;

##MIPs lib/
use lib catdir(dirname($Bin), "lib");
use Check::Check_modules qw(check_modules);
use Script::Utils qw(help);

our $USAGE;

BEGIN {

    ## Special case to initiate testing
    my @modules = ("Test::More",
	);

    # Evaluate that all modules required are installed
    Check::Check_modules::check_modules({modules_ref => \@modules,
					 program_name => $0,
					});

    ##Initate tests
    print STDOUT "Initiate tests:\n";
    print STDOUT "\nTesting perl modules and selected functions\n\n";


    ##More proper testing 
    use Test::More;

    ##Modules with import
    my %perl_module;

    $perl_module{autodie} = [qw(open close :all)];
    $perl_module{charnames} = [qw(:full :short)];
    $perl_module{Cwd} = [qw(abs_path)];
    $perl_module{"File::Basename"} = [qw(dirname basename)];
    $perl_module{"File::Path"} = [qw(make_path remove_tree)];
    $perl_module{"File::Spec::Functions"} =  [qw(catfile catdir devnull)];
    $perl_module{FindBin} = [qw($Bin)];
    $perl_module{"List::Util"} = [qw(any all uniq)];
    $perl_module{"IPC::Cmd"} = [qw(can_run run)];
    $perl_module{"Modern::Perl"} = [qw(2014)];
    $perl_module{open} = [qw(:encoding(UTF-8) :std)];
    $perl_module{"Params::Check"} = [qw(check allow last_error)];
    $perl_module{warnings} = [qw(FATAL utf8)];
    
    while (my ($module, $module_import) = each %perl_module) { 

	use_ok($module, @$module_import) or BAIL_OUT "Can't load $module";
    }
    
    ##Modules
    @modules = ("Cwd",
		"Getopt::Long",
		"IO::Handle",
		"IPC::System::Simple",
		"Log::Log4perl",
		"Path::Iterator::Rule",
		"POSIX",
		"strict",
		"TAP::Harness",
		"Time::Piece",
		"utf8",
		"YAML",
		"warnings",
	);

    for my $module (@modules) {

	require_ok($module) or BAIL_OUT "Can't load $module";
    }

    $USAGE =
	basename($0).qq{
           -h/--help Display this help message   
           -v/--version Display version
        };    
}



my $run_tests_version = "0.0.0";

###User Options
GetOptions('h|help' => sub { done_testing(); print STDOUT $USAGE, "\n"; exit;},  #Display help text
	   'v|version' => sub { done_testing(); print STDOUT "\n".basename($0)." ".$run_tests_version, "\n\n"; exit;},  #Display version number
    ) or done_testing(), Script::Utils::help({USAGE => $USAGE,
					      exit_code => 1,
					     });

use TAP::Harness;
use Cwd;


## Test central perl modules and import functions
test_modules();

mip_scripts();

done_testing(); # Reached the end safely

my %args = (verbosity => 1);

my $harness = TAP::Harness->new( \%args );
my @tests = (catfile(getcwd(), "mip.t"),
    );
$harness->runtests(@tests);

######################
####SubRoutines#######
######################

sub test_modules {

##test_modules
    
##Function : Test perl modules and functions 
##Returns  : ""
##Arguments: 
##         : 

    use Cwd;
    ok(getcwd(), "Cwd: Locate current working directory");

    use FindBin qw($Bin);
    ok(defined($Bin),"FindBin: Locate directory of script");

    use File::Basename qw(dirname);
    ok(dirname($Bin), "File::Basename qw(dirname): Strip the last part of directory");
 
    use File::Spec::Functions qw(catfile catdir devnull);
    ok(catdir(dirname($Bin), "t"),"File::Spec::Functions qw(catdir): Concatenate directories");
    ok(catfile($Bin, "run_tests.t"),"File::Spec::Functions qw(catfile): Concatenate files");
    ok(catfile(dirname(devnull()), "stdout"), "File::Spec::Functions qw(devnull): Use devnull");
 
    use File::Basename qw(basename);
    my $file_path = catfile($Bin, "run_tests.t");
    ok(basename($file_path), "File::Basename qw(basename): Strip directories");

    use Cwd qw(abs_path);
    ok(abs_path( catfile($Bin, "run_tests.pl") ), "Cwd_abs_path: Add absolute path");

    use File::Path qw(make_path remove_tree);
    ok(make_path("TEST"), "File::Path_make_path: Create path");
    ##Clean-up
    ok(remove_tree("TEST"), "File::Path_remove_tree: Remove path");

    ##MIPs lib/
    use lib catdir(dirname($Bin), "lib");
    use File::Format::Yaml qw(load_yaml);
    use YAML;    
    my $yaml_file = catdir(dirname($Bin), "templates", "643594-miptest_pedigree.yaml");
    ok( -f $yaml_file, "YAML: File= $yaml_file in MIP/templates directory");
    
    my $yaml = File::Format::Yaml::load_yaml({yaml_file => $yaml_file,
					     });
    ok( defined $yaml, "YAML: Load File" );  #Check that we got something
    ok(Dump( $yaml ), "YAML: Dump file");
    
    use Params::Check qw[check allow last_error];
    use Log::Log4perl;
    ## Creates log
    my $log_file = catdir(dirname($Bin), "templates", "mip_log.yaml");
    ok( -f $log_file, "Log::Log4perl: File= $log_file in MIP directory");

    use MIP_log::Log4perl qw(initiate_logger);
    ## Creates log object
    my $log = MIP_log::Log4perl::initiate_logger({categories_ref => ["TRACE", "ScreenApp"],
						  file_path_ref => \$log_file,
						  log_name => "Run_tests",
						 });
    
    ok($log->info("1"), "Log::Log4perl: info");
    ok($log->warn("1"), "Log::Log4perl: warn");
    ok($log->error("1"), "Log::Log4perl: error");
    ok($log->fatal("1"), "Log::Log4perl: fatal");

    use Getopt::Long;
    push(@ARGV, ("-verbose", "2"));
    my $verbose = 1;
    ok(GetOptions("verbose:n"  => \$verbose), "Getopt::Long: Get options call");
    ok ($verbose == 2, "Getopt::Long: Get options modified");

    ## Check time
    use Time::Piece;
    my $dateTime = localtime;
    ok($dateTime, "localtime = $dateTime");
    my $dateTimeStamp = $dateTime->datetime;
    ok($dateTimeStamp, "datetime = $dateTime");
    my $date = $dateTime->ymd;
    ok($date, "ymd = $date");

    ## Locate name of script
    my $script = (`basename $0`);
    ok((`basename $0`), "Detect script name = $script");

    ## Execution of programs
    use IPC::Cmd qw[can_run run];
    ok(can_run("perl"), "Can run IPC::Cmd");
    ok(my $bool = IPC::Cmd->can_capture_buffer, "IPC::Cmd can capture buffer");
}


sub mip_scripts{

##mip_scripts
    
##Function : Test MIP kit completion
##Returns  : ""
##Arguments: 
##         :

    my @mip_scripts = ("calculate_af.pl",
		       "download_reference.pl",
		       "mip_install.pl",
		       "max_af.pl",
		       "mip.pl",
		       "qccollect.pl",
		       "vcfparser.pl",
	);

    foreach my $script (@mip_scripts) {
	
	is(-e catfile(dirname($Bin), $script), 1, "Found MIP file: ".$script);
    }

    my %mip_sub_scripts;
    $mip_sub_scripts{"definitions"} = ["define_download_references.yaml",
				       "define_parameters.yaml",
	];
    $mip_sub_scripts{"t"} = ["mip_install.t",
			     "mip.t",
			     "run_tests.t",
			     "mip_analysis.t",
			     
	];
    $mip_sub_scripts{"templates"} = ["mip_config.yaml",
				     "mip_travis_config.yaml",
				     "643594-miptest_pedigree.yaml",
				     "mip_log.yaml",
	];

    foreach my $directory (keys %mip_sub_scripts) {
	
	foreach my $script (@{$mip_sub_scripts{$directory}}) {
	    
	    is(-e catfile(dirname($Bin), $directory, $script), 1, "Found MIP file: ".$script);
	}
    }
    my @mip_directories = ("lib",
			   catdir("t", "data"),
	);
    foreach my $directory (@mip_directories) {
	
	is(-e catfile(dirname($Bin), $directory), 1, "Found MIP dir: ".$directory);
    }
}
