#!/usr/bin/env perl

###Copyright 2016 Henrik Stranneheim

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
    print STDOUT "\nTesting perl modules and selected functions\n\n";


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
		"IO::Handle",
		"utf8",
		"strict",
		"TAP::Harness",
	);

    for my $module (@modules) {
	require_ok($module) or BAIL_OUT "Can't load $module";
    }
}

use Getopt::Long;

our $USAGE;

BEGIN {
    $USAGE =
	qq{run_tests.t
           -h/--help Display this help message   
           -v/--version Display version
        };    
}

my $run_tests_version = "0.0.0";

###User Options
GetOptions('h|help' => sub { print STDOUT $USAGE, "\n"; exit;},  #Display help text
	   'v|version' => sub { print STDOUT "\ntest.t ".$run_tests_version, "\n\n"; exit;},  #Display version number
);

use TAP::Harness;

## Test central perl modules and import functions
test_modules();

mip_scripts();

done_testing(); # Reached the end safely

my %args = (verbosity => 1);

my $harness = TAP::Harness->new( \%args );
my @tests = [ "install.t" ];
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

    use FindBin qw($Bin); #Find directory of script

    ok(defined($Bin),"FindBin: Locate directory of script");

    use File::Basename qw(dirname);  #Strip the last part of directory

    ok(dirname($Bin), "File::Basename qw(dirname): Strip the last part of directory");
 
    use File::Spec::Functions qw(catfile catdir);

    ok(catdir(dirname($Bin), "t"),"File::Spec::Functions qw(catdir): Concatenate directories");
    ok(catfile($Bin, "run_tests.pl"),"File::Spec::Functions qw(catfile): Concatenate files");
    
    use YAML;
    
    my $yaml_file = catdir(dirname($Bin), "templates", "118_pedigree.yaml");
    ok( -f $yaml_file,"YAML: File= $yaml_file in MIP/templates directory");
    
    my $yaml = YAML::LoadFile($yaml_file);  #Create an object
    ok( defined $yaml,"YAML: Load File" );  #Check that we got something
    ok(Dump( $yaml ),"YAML: Dump file");

    use Log::Log4perl;
    use Params::Check qw[check allow last_error];
    ## Creates log
    my $log_file = catdir(dirname($Bin), "templates", "mip_log.yaml");
    ok( -f $log_file,"Log::Log4perl: File= $log_file in MIP directory");

    ## Create log4perl config file
    my $config = &create_log4perl_congfig({file_path_ref => \$log_file});
    
    ok(Log::Log4perl->init(\$config), "Log::Log4perl: Initate");
    ok(Log::Log4perl->get_logger("MIP_logger"), "Log::Log4perl: Get logger");
    
    my $logger = Log::Log4perl->get_logger("MIP_logger");
    ok($logger->info("1"), "Log::Log4perl: info");
    ok($logger->warn("1"), "Log::Log4perl: warn");
    ok($logger->error("1"), "Log::Log4perl: error");
    ok($logger->fatal("1"), "Log::Log4perl: fatal");

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
    ok($bool = IPC::Cmd->can_capture_buffer, "IPC::Cmd can capture buffer");
}

sub create_log4perl_congfig {

##create_log4perl_congfig

##Function : Create log4perl config file.
##Returns  : "$config"
##Arguments: $file_path_ref
##         : $file_path_ref => log4perl config file path {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path_ref;

    my $tmpl = {
	file_path_ref => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$file_path_ref},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $config = q?
        log4perl.category.MIP_logger = TRACE, LogFile, ScreenApp
        log4perl.appender.LogFile = Log::Log4perl::Appender::File
        log4perl.appender.LogFile.filename = ?.$$file_path_ref.q?
        log4perl.appender.LogFile.layout=PatternLayout
        log4perl.appender.LogFile.layout.ConversionPattern = [%p] %d %c - %m%n

        log4perl.appender.ScreenApp = Log::Log4perl::Appender::Screen
        log4perl.appender.ScreenApp.layout = PatternLayout
        log4perl.appender.ScreenApp.layout.ConversionPattern = [%p] %d %c - %m%n
        ?;
    return $config;
}


sub LoadYAML {
 
##LoadYAML
    
##Function : Loads a YAML file into an arbitrary hash and returns it. Note: Currently only supports hashreferences and hashes and no mixed entries.
##Returns  : %yaml_hash
##Arguments: $yaml_file
##         : $yaml_file => The yaml file to load

    my ($arg_hef) = @_;

    ##Flatten argument(s)
    my $yaml_file;

    my $tmpl = { 
	yaml_file => { required => 1, defined => 1, strict_type => 1, store => \$yaml_file},
    };

    check($tmpl, $arg_hef, 1) or die qw[Could not parse arguments!];

    my %yaml_hash;

    open (my $YAML, "<", $yaml_file) or die "can't open ".$yaml_file.":".$!, "\n";  #Log4perl not initialised yet, hence no logdie
    local $YAML::QuoteNumericStrings = 1;  #Force numeric values to strings in YAML representation
    %yaml_hash = %{ YAML::LoadFile($yaml_file) };  #Load hashreference as hash
        
    close($YAML);

    return %yaml_hash;
}


sub mip_scripts{

##mip_scripts
    
##Function : Test MIP kit completion
##Returns  : ""
##Arguments: 
##         :

    my @mip_scripts = ("calculate_af.pl",
		       "max_af.pl",
		       "mip.pl",
		       "qccollect.pl",
		       "vcfparser.pl",
		       "install.pl",
	);

    foreach my $script (@mip_scripts) {
	
	is(-e catfile(dirname($Bin), $script), 1, "Found MIP file: ".$script);
    }

    my %mip_sub_scripts;
    $mip_sub_scripts{"definitions"} = ["define_parameters.yaml"];
    $mip_sub_scripts{"t"} = ["run_tests.t",
			     "test.t",
			     "install.t",
	];
    $mip_sub_scripts{"templates"} = ["mip_config.yaml",
				     "118_pedigree.yaml",
				     "mip_log.yaml",
	];

    foreach my $directory (keys %mip_sub_scripts) {
	
	foreach my $script (@{$mip_sub_scripts{$directory}}) {
	    
	    is(-e catfile(dirname($Bin), $directory, $script), 1, "Found MIP file: ".$script);
	}
    }
}
