#!/usr/bin/env perl

use Modern::Perl '2014';
use warnings qw( FATAL utf8 );
use autodie qw(open close :all);
use v5.18;  #Require at least perl 5.18
use utf8;
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );

use Cwd;
use Cwd qw(abs_path);
use File::Basename qw(dirname basename);
use File::Spec::Functions qw(catfile catdir devnull);
use FindBin qw($Bin); #Find directory of script
use Getopt::Long;
use IO::Handle;
use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case

## Third party module(s)
use List::Util qw(any);

##MIPs lib/
use lib catdir($Bin, "lib");
use Check::Check_modules qw(check_modules);
use File::Format::Shell qw(create_bash_file);
use File::Format::Yaml qw(load_yaml);
use MIP_log::Log4perl qw(initiate_logger);
use Script::Utils qw(help set_default_array_parameters);

our $USAGE;

BEGIN {

    my @modules = ("Modern::Perl",
		   "autodie",
		   "YAML",
		   "File::Format::Yaml",
		   "Log::Log4perl",
		   "MIP_log::Log4perl",
		   "Script::Utils",
	);

    ## Evaluate that all modules required are installed
    Check::Check_modules::check_modules({modules_ref => \@modules,
					 program_name => "download_reference",
					});

    $USAGE =
	basename($0).qq{ [options]
           -rd/--reference_dir Reference(s) directory (Default: "")
           -r/--reference Reference to download (e.g. 'clinvar=20170104')
           -rg/--reference_genome_versions Reference versions to download ((Default: ["GRCh37", "hg38"]))
           -l/--log_file Log file (Default: "download_reference.log")
           -h/--help Display this help message
           -ver/--version Display version
           -v/--verbose Set verbosity
        };
}

### Set parameter default

## Loads a YAML file into an arbitrary hash and returns it.
my %parameter = load_yaml({yaml_file => catfile($Bin, "definitions", "define_download_references.yaml"),
			  });

## Set parameter default
$parameter{reference_dir} = cwd();

## Define default parameters
my %array_parameter;
$array_parameter{reference_genome_versions}{default} = ["GRCh37", "hg38"];

my $download_reference_version = "0.0.1";


###User Options
GetOptions('rd|reference_dir:s' => \$parameter{reference_dir},  #MIPs reference directory
	   'r|reference:s' => \%{ $parameter{cmd_reference} },
	   'rg|reference_genome_versions:s' => \@{ $parameter{reference_genome_versions} },
	   'l|log_file:s' => \$parameter{log_file},
	   'h|help' => sub { print STDOUT $USAGE, "\n"; exit;},  #Display help text
	   'ver|version' => sub { print STDOUT "\n".basename($0)." ".$download_reference_version, "\n\n"; exit;},  #Display version number
	   'v|verbose' => \$parameter{verbose},
    ) or Script::Utils::help({USAGE => $USAGE,
			      exit_code => 1,
			     });

## Creates log object
my $log = MIP_log::Log4perl::initiate_logger({file_path_ref => \$parameter{log_file},
					      log_name => "Download_reference",
					     });
check_user_reference({cmd_reference_ref => \%{ $parameter{cmd_reference} },
		      reference_ref => \% { $parameter{reference} },
		     });


## Set default for array parameters
Script::Utils::set_default_array_parameters({parameter_href => \%parameter,
					     array_parameter_href => \%array_parameter,
					    });

## Change relative path to absolute path for certain parameters
update_to_absolute_path({parameter_href => \%parameter,
			});

##########
###MAIN###
##########


## Create bash file for writing install instructions
my $BASHFILEHANDLE = File::Format::Shell::create_bash_file({file_name => "download_reference.sh",
							    directory_remove => ".download_reference",
							    log => $log,
							   });
$log->info("Will write install instructions to '".catfile(cwd(), "download_reference.sh"), "'\n");

references({parameter_href => \%parameter,
	    FILEHANDLE => $BASHFILEHANDLE,
	   });

###SubRoutines###


sub references {

##references

##Function : Install references
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE, $verbose
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to
##         : $verbose        => Verbosity

    my ($arg_href) = @_;

    ## Default(s)
    my $verbose;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
	verbose => { default => $arg_href->{parameter_href}{verbose},
		     strict_type => 1, store => \$verbose},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    use MIP::Gnu::Bash qw(gnu_cd);
    use MIP::Gnu:Coreutils qw (gnu_mkdir);

    ## Retrieve logger object now that log_file has been set
    my $log = Log::Log4perl->get_logger("Download_reference");

    my $pwd = cwd();

    print $FILEHANDLE "## Create reference directory\n";
    gnu_mkdir({indirectory_path => $parameter_href->{reference_dir},
				    parents => 1,
				    FILEHANDLE => $FILEHANDLE,
				   });
    print $FILEHANDLE "\n\n";

    ## Since all commands should assume working directory to be the reference directory
    gnu_cd({directory_path => $parameter_href->{reference_dir},
	FILEHANDLE => $FILEHANDLE,
       });
    print $FILEHANDLE "\n\n";

  REFERENCE:
    while (my ($reference_id, $versions_ref) = each(%{ $parameter_href->{reference} }) ) {

	## Remodel depending on if "--reference" was used or not as the user info is stored as a scalar per reference_id while yaml is stored as arrays per reference_id

	my @reference_versions;
	if (ref($versions_ref) eq "ARRAY") {

	    @reference_versions = @$versions_ref;
	}
	else {

	    push (@reference_versions, $versions_ref);
	}
	
      REFERENCE_VERSION:
	foreach my $reference_version (@reference_versions) {

	  GENOME_VERSION:
	    foreach my $genome_version (@{ $parameter_href->{reference_genome_versions} })    {

		$genome_version = lc($genome_version);
		my $reference_href = $parameter_href->{$reference_id}{$genome_version}{$reference_version};

		if( (exists($parameter_href->{$reference_id}{$genome_version}))
		    && (exists($parameter_href->{$reference_id}{$genome_version}{$reference_version})) ) {

		    ## Build file name and path
		    my $outfile_name = $reference_href->{outfile};
		    my $outfile_path = catfile($parameter_href->{reference_dir}, $outfile_name);
	    
		    ## Check if reference already exists in reference directory
		    if (! -f $outfile_path) {
		    
			$log->warn("Cannot find reference file:".$outfile_path, "\n");
			$log->warn("Will try to download: ".$reference_id." version: ".$reference_version, "\n");

			## Potential download files
			my @file_keys = ("file",
					 "file_check",
					 "file_index",
					 "file_index_check");

		      REFERENCE_FILES:
			foreach my $key (@file_keys) {

			    ## Install reference
			    if (exists($reference_href->{$key})) {
				
				my $file = $reference_href->{$key};
				my $outfile = $reference_href->{"out".$key};
				my $outfile_path = catfile($parameter_href->{reference_dir}, $outfile);
				
				download({parameter_href => $parameter_href,
					  FILEHANDLE => $FILEHANDLE,
					  url => $reference_href->{url_prefix}.$file,
					  outfile_path => $outfile_path,
					  file_id => $reference_id,
					 });
				
				## Check if file needs to be decompress and write decompression if so
				decompress_file({parameter_href => $parameter_href,
						 FILEHANDLE => $FILEHANDLE,
						 outfile_path => $outfile_path,
						 file_decompress => $reference_href->{"out".$key."_decompress"},
						});

				## Check file integrity of file
				check_file({FILEHANDLE => $FILEHANDLE,
					    outfile_path => $outfile_path,
					    outfile_path_check => $outfile_path,
					    check_method => $reference_href->{"out".$key."_method"},
					   });
			    }
			}

			## Process reference with commands
			if ( ( map { $_ =~/_command$/ } keys %$reference_href ) ) { #If key contains command

			    ## Reformat command
			    if (exists($reference_href->{outfile_reformat_command})) {
				
				print $FILEHANDLE $reference_href->{outfile_reformat_command}, "\n\n";
			    }

			    ## Bgzip command
			    if (exists($reference_href->{outfile_bgzip_command})) {
				
				print $FILEHANDLE $reference_href->{outfile_bgzip_command}, "\n\n";
			    }

			    ## Tabix command
			    if (exists($reference_href->{outfile_tabix_command})) {
				
				print $FILEHANDLE $reference_href->{outfile_tabix_command}, "\n\n";
			    }
			}
		    }
		}
	    }
	}
    }

    ## Move back to original
    cd({directory_path => $pwd,
	FILEHANDLE => $FILEHANDLE,
       });
    print $FILEHANDLE "\n\n";
}


sub download {

##download

##Function : Downloads files
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE, $url, $outfile_path, $file_id, $program, $quiet, $verbose
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to
##         : $url            => Url to use for download
##         : $outfile_path   => Outfile path 
##         : $program        => Program to use for download
##         : $file_id        => File id
##         : $quiet          => Quiet (no output)
##         : $verbose        => Verbosity

    my ($arg_href) = @_;

    ## Default(s)
    my $program;
    my $quiet;
    my $verbose;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;
    my $url;
    my $outfile_path;
    my $file_id;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
	url => { required => 1, defined => 1, strict_type => 1, store => \$url},
	outfile_path => { required => 1, defined => 1, strict_type => 1, store => \$outfile_path},
	file_id => { required => 1, defined => 1, strict_type => 1, store => \$file_id},
	program => { default => "wget",
		     allow => ["wget"],
		     strict_type => 1, store => \$program},
	quiet => { default => 1,
		   allow => [0, 1],
		   strict_type => 1, store => \$quiet},
	verbose => { default => $arg_href->{parameter_href}{verbose},
		     strict_type => 1, store => \$verbose},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Download
    print $FILEHANDLE "## Download ".$file_id, "\n";

    if($program eq "wget") {
	
	use Program::Download::Wget qw(wget);

	wget({url => $url,
	      FILEHANDLE => $FILEHANDLE,
	      quiet => $quiet,
	      verbose => $verbose,
	      outfile_path => $outfile_path,
	     });
	say $FILEHANDLE "\n";
    }
}


sub remove_file_ending {

##remove_file_ending

##Function : Removes ".file_ending" in filename.file_ending(.gz)
##Returns  : File name with supplied $file_ending or $file_ending(.gz) removed
##Arguments: $file_name_ref, $file_ending
##         : $file_name_ref => File name {REF}
##         : $file_ending   => File ending to be removed

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_name_ref;
    my $file_ending;

    my $tmpl = {
	file_name_ref => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$file_name_ref},
	file_ending => { required => 1, defined => 1, strict_type => 1, store => \$file_ending},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $file_name_noending;

    if ( (defined($$file_name_ref)) && ($$file_name_ref =~/(\S+)($file_ending$|$file_ending.gz$)/) ) {

	$file_name_noending = $1;
    }
    return $file_name_noending;
}


sub decompress_file {

##decompress_file

##Function : Check if file needs to be decompress and write decompression if so
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE, $outfile_path, $file_decompress
##         : $parameter_href   => Holds all parameters
##         : $FILEHANDLE       => Filehandle to write to
##         : $outfile_path     => Outfile path
##         : $file_decompress  => Decompress the downloaded file
 
    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;
    my $outfile_path;
    my $file_decompress;

    my $tmpl = { 
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
	outfile_path => { required => 1, defined => 1, strict_type => 1, store => \$outfile_path},
	file_decompress => { strict_type => 1, store => \$file_decompress},
    };
     
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    if (defined($outfile_path)) {

	if ( (defined($file_decompress)) && ($file_decompress eq "gzip") ) {

	    ## Removes ".file_ending" in filename.FILENDING(.gz)
	    my $outfile_path_no_ending = remove_file_ending({file_name_ref => \$outfile_path,
							 file_ending => ".gz",
							});
	
	    print $FILEHANDLE "gzip ";
	    print $FILEHANDLE "-f ";
	    print $FILEHANDLE "--quiet ";
	    print $FILEHANDLE "-d ";
	    print $FILEHANDLE "-c ";
	    print $FILEHANDLE $outfile_path." ";
	    print $FILEHANDLE "> ".$outfile_path_no_ending, "\n\n";
	}

	if ( (defined($file_decompress)) && ($file_decompress eq "unzip") ) {
	
	    print $FILEHANDLE "unzip ";
	    print $FILEHANDLE "-d ".$parameter_href->{reference_dir}." ";
	    print $FILEHANDLE $outfile_path, "\n\n";;
	}
    }
}


sub check_file {

##check_file

##Function : Check file integrity of file
##Returns  : ""
##Arguments: $FILEHANDLE, $outfile_path, $outfile_path_check, $check_method
##         : $FILEHANDLE          => Filehandle to write to
##         : $outfile_path        => Outfile path
##         : $outfile_path_check  => File to check
##         : check_method         => Method to perform file check
 
    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $outfile_path;
    my $outfile_path_check;
    my $check_method;

    my $tmpl = { 
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
	outfile_path => { required => 1, defined => 1, strict_type => 1, store => \$outfile_path},
	outfile_path_check => { strict_type => 1, store => \$outfile_path_check},
	check_method => { strict_type => 1, store => \$check_method},
    };
     
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    use Program::Gnu::Coreutils qw(rm);

    if( (defined($check_method)) && ($check_method eq "md5sum") ) {

	## Removes ".file_ending" in filename.FILENDING(.gz)
	my $outfile_path_no_ending = remove_file_ending({file_name_ref => \$outfile_path,
							 file_ending => ".md5",
							});
	if (defined($outfile_path_no_ending)) {
	    
	    my $perl_regexp = q?perl -nae 'print $F[0]."  ?.$outfile_path_no_ending.q?" ' ?.$outfile_path_check;
	    print $FILEHANDLE $perl_regexp." > md5sum_check.txt", "\n\n";
	    print $FILEHANDLE "md5sum ";
	    print $FILEHANDLE "-c md5sum_check.txt", "\n\n";

	    ## Clean-up
	    rm({infile_path => "md5sum_check.txt",
		force => 1,
		FILEHANDLE => $FILEHANDLE,
	       });
	    print $FILEHANDLE "\n\n";
	}
    }
}


sub update_to_absolute_path {

##update_to_absolute_path

##Function : Change relative path to absolute path for certain parameter_names
##Returns  : ""
##Arguments: $parameter_href
##         : $parameter_href => The parameter hash {REF}

    my ($arg_href) = @_;

    ##Flatten argument(s)
    my $parameter_href;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    use File::Parse::Parse qw(find_absolute_path);

    ## Retrieve logger object now that log_file has been set
    my $log = Log::Log4perl->get_logger("Download_reference");
    
    foreach my $parameter_name (@{ $parameter_href->{absolute_paths} }) {

	if(defined($parameter{$parameter_name})) {

	    if (ref($parameter_href->{$parameter_name}) eq "ARRAY") {  #Array reference
	    
		foreach my $parameter_value (@{ $parameter_href->{$parameter_name} }) {
		
		    ## Replace original input with abolute path for supplied path or croaks and exists if path does not exists
		    $parameter_value = find_absolute_path({path => $parameter_value,
							   parameter_name => $parameter_name,
							   log => $log,
							  });
		}
	    }
	    elsif (ref($parameter_href->{$parameter_name}) eq "HASH") {  #Hash reference

		foreach my $key (keys %{ $parameter_href->{$parameter_name} }) {  #Cannot use each since we are updating key

		    ## Find aboslute path for supplied path or croaks and exists if path does not exists
		    my $updated_key = find_absolute_path({path => $key,
							  parameter_name => $parameter_name,
							  log => $log,
							 });
		    $parameter_href->{$parameter_name}{$updated_key} = delete($parameter_href->{$parameter_name}{$key});
		}
	    }
	    else {  #Scalar - not a reference
		
		## Find aboslute path for supplied path or croaks and exists if path does not exists
		$parameter_href->{$parameter_name} = find_absolute_path({path => $parameter_href->{$parameter_name},
									 parameter_name => $parameter_name,
									 log => $log,
									});
	    }
	}
    }
}


sub check_user_reference {

##check_user_reference

##Function : Check that the user supplied reference id and version
##Returns  : ""
##Arguments: $cmd_reference_ref, $reference_ref,
##         : $cmd_reference_ref => User supplied reference id and version
##         : $reference_ref     => Defined reference id and version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $cmd_reference_ref;
    my $reference_ref;

    my $tmpl = { 
	cmd_reference_ref => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$cmd_reference_ref},
	reference_ref => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$reference_ref},
    };
     
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    while (my ($reference_id, $version) = each(%$cmd_reference_ref) ) {

	if(! exists($reference_ref->{$reference_id})) {

	    $log->fatal("Cannot find reference key:".$reference_id, "\n");
	    exit;
	}
	elsif (! any {$_ eq $version} @{ $reference_ref->{$reference_id} }) {  #If element is part of array

	    $log->fatal("Cannot find version key: ".$version." reference key:".$reference_id, "\n");
	    exit;
	}
    }
    
}
