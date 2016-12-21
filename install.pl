#!/usr/bin/env perl

##Assumes you have a working conda installation

use strict;
use warnings;
use warnings qw( FATAL utf8 );
use utf8;
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );

use Getopt::Long;
use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case
use Cwd;
use FindBin qw($Bin); #Find directory of script
use IO::Handle;
use File::Basename qw(dirname);
use File::Spec::Functions qw(catfile catdir devnull);

our $USAGE;

BEGIN {
    $USAGE =
	qq{install.pl [options]
           -env/--conda_environment Conda environment (Default: "mip")
           -cdp/--conda_path The conda path (Default: "HOME/miniconda")
           -cdu/--conda_update Update conda before installing (Supply flag to enable)
           -bvc/--bioconda Set the module version of the programs that can be installed with bioconda (e.g. 'bwa=0.7.12')
           -pip/--pip Set the module version of the programs that can be installed with pip (e.g. 'genmod=3.5.9')
           -pyv/--python_version Set the env python version (Default: "2.7")

           ## SHELL
           -pei/--perl_install Install perl (Supply flag to enable)
           -pev/--perl_version Set the perl version (defaults: "5.18.2")
           -pevs/--perl_skip_test Skip "tests" in perl installation
           -pm/--perl_modules Set the perl modules to be installed via cpanm (Default: ["Modern::Perl", "List::Util", "IPC::System::Simple", "Path::Iterator::Rule", "YAML", "Log::Log4perl", "Set::IntervalTree", "Net::SSLeay",P, "LWP::Simple", "LWP::Protocol::https", "Archive::Zip", "Archive::Extract", "DBI","JSON", "DBD::mysql", "CGI", "Sereal::Encoder", "Sereal::Decoder", "Bio::Root::Version", "Module::Build"])
           -pmf/--perl_modules_force Force installation of perl modules
           -pic/--picardtools Set the picardtools version (Default: "2.3.0"),
           -sbb/--sambamba Set the sambamba version (Default: "0.6.1")
           -vct/--vcftools Set the vcftools version (Default: "0.1.14")
           -bet/--bedtools Set the bedtools version (Default: "2.25.0")
           -vt/--vt Set the vt version (Default: "0.57")
           -plk/--plink  Set the plink version (Default: "160224")
           -snpg/--snpeff_genome_versions Set the snpEff genome version (Default: ["GRCh37.75"])
           -vep/--varianteffectpredictor Set the VEP version (Default: "87")
           -vepa/--vep_auto_flag Set the VEP auto installer flags
	   -vepc/--vep_cache_dir Specify the cache directory to use (whole path; defaults to "~/miniconda/envs/conda_environment/ensembl-tools-release-varianteffectpredictorVersion/cache")
           -vepa/--vep_assemblies Select the assembly version (Default: ["GRCh37"])
           -vepp/--vep_plugin Supply a comma separated list of VEP plugins (Default: "UpDownDistance,LoFtool,Lof")
           -rhc/--rhocall Set the rhocall version (Default: "0.1")
           -rhcp/--rhocall_path Set the path to where to install rhocall (Defaults: "HOME/rhocall")

           ## Utility
           -psh/--prefer_shell Shell will be used for overlapping shell and biconda installations (Supply flag to enable)
           -ppd/--print_parameters_default Print the parameter defaults
           -nup/--noupdate Do not update already installed programs (Supply flag to enable)
           -sp/--select_programs Install supplied programs e.g. -sp perl -sp bedtools (Default: "";)
           -h/--help Display this help message
           -v/--version Display version
        };
}


my %parameter;

### Set parameter default

##Conda
$parameter{conda_environment} = "mip";
$parameter{conda_path} = catdir($ENV{HOME}, "miniconda");
$parameter{python_version} = "2.7";

$parameter{bioconda}{bwa} = "0.7.15";
$parameter{bioconda}{bwakit} = "0.7.12";
$parameter{bioconda_bwakit_patch} = "-0";  #For correct softlinking in share and bin in conda env
$parameter{bioconda}{fastqc} = "0.11.5";
$parameter{bioconda}{cramtools} = "3.0.b47";
$parameter{bioconda}{samtools} = "1.3.1";
$parameter{bioconda}{bcftools} = "1.3.1";
$parameter{bioconda}{snpeff} = "4.2";
$parameter{bioconda_snpeff_patch} = "-0";  #For correct softlinking in share and bin in conda env
$parameter{bioconda}{picard} = "2.5.0";
$parameter{bioconda_picard_patch} = "-1";  #For correct softlinking in share and bin in conda env
$parameter{bioconda}{mosaik} = "2.2.26";
$parameter{bioconda}{htslib} = "1.3.1";
$parameter{bioconda}{bedtools} = "2.26.0";
$parameter{bioconda}{vt} = "2015.11.10";
$parameter{bioconda}{sambamba} = "0.6.3";
$parameter{bioconda}{freebayes} = "1.0.2.0";
$parameter{bioconda}{delly} = "0.7.2";
$parameter{bioconda}{manta} = "1.0.0";
$parameter{bioconda_manta_patch} = "-0";
$parameter{bioconda}{multiqc} = "0.8dev0";
$parameter{bioconda}{plink2} = "1.90b3.35";
$parameter{bioconda}{vcfanno} = "0.1.0";
$parameter{bioconda}{gcc} = "4.8.5";  #Required for CNVnator
#$parameter{bioconda}{cmake} = "3.3.1";
#$parameter{bioconda}{boost} = "1.57.0";
#$parameter{bioconda_boost_patch} = "-4";


##Perl Modules
$parameter{perl_version} = "5.18.2";

## PIP
$parameter{pip}{genmod} = "3.5.9";
$parameter{pip}{variant_integrity} = "0.0.4";
$parameter{pip}{chanjo} = "4.0.0";
$parameter{pip}{cosmid} = "0.4.9.1";
$parameter{pip}{'python-Levenshtein'} = "0.12.0";

## Programs currently installable by SHELL
$parameter{mip_scripts} = "Your current MIP version";
$parameter{picardtools} = "2.3.0";
$parameter{sambamba} = "0.6.1";
$parameter{vcftools} = "0.1.14";
$parameter{bedtools} = "2.25.0";
$parameter{vt} = "gitRepo";
$parameter{plink2} = "160316";
$parameter{snpeff} = "v4_2";
$parameter{varianteffectpredictor} = "87";
$parameter{vep_auto_flag} = "alcf";
$parameter{vep_plugin} = "UpDownDistance,LoFtool,Lof";
$parameter{rhocall} = "0.3";
$parameter{rhocall_path} = catdir($ENV{HOME}, "rhocall");

#$parameter{cnvnator} = "0.3.2";
#$parameter{findtranslocations} = "0";

my $install_version = "1.0.0";

###User Options
GetOptions('env|conda_environment:s'  => \$parameter{conda_environment},
	   'cdp|conda_path:s' => \$parameter{conda_path},
	   'cdu|conda_update' => \$parameter{conda_update},
	   'bcv|bioconda=s' => \%{ $parameter{bioconda} },
	   'pip|pip=s' => \%{ $parameter{pip} },
	   'pyv|python_version=s' => \$parameter{python_version},
	   'pev|perl_version=s' => \$parameter{perl_version},
	   'pei|perl_install' => \$parameter{perl_install},
	   'pevs|perl_skip_test' => \$parameter{perl_skip_test},
	   'pm|perl_modules:s' => \@{ $parameter{perl_modules} },  #Comma separated list
           'pmf|perl_modules_force' =>  \$parameter{perl_modules_force},
	   'pic|picardtools:s' => \$parameter{picardtools},
	   'sbb|sambamba:s' => \$parameter{sambamba},
	   'vct|vcftools:s' => \$parameter{vcftools},
	   'bet|bedtools:s' =>\$parameter{bedtools},
	   'vt|vt:s' => \$parameter{vt},
	   'plk|plink2:s' => \$parameter{plink2},
	   'snpg|snpeff_genome_versions:s' => \@{ $parameter{snpeff_genome_versions} },
	   'vep|varianteffectpredictor:s' => \$parameter{varianteffectpredictor},
	   'vepai|vep_auto_flag:s' => \$parameter{vep_auto_flag},
	   'vepc|vep_cache_dir:s' => \$parameter{vep_cache_dir},  #path to vep cache dir
	   'vepa|vep_assemblies:s' => \@{ $parameter{vep_assemblies} },  #Select assembly version to use
	   'vepp|vep_plugin:s' => \$parameter{vep_plugin},  #Comma sep string
	   'rhc|rhocall:s' => \$parameter{rhocall},
	   'rhcp|rhocall_path:s' => \$parameter{rhocall_path},
#	   'cnv|cnvnator:s' => \$parameter{cnvnator},
#	   'ftr|findtranslocations:s' => \$parameter{findtranslocations},
	   'psh|prefer_shell' => \$parameter{prefer_shell},  # Shell will be used for overlapping shell and biconda installations
	   'ppd|print_parameters_default' => sub { print_parameters({parameter_href => \%parameter}); exit;},  #Display parameter defaults
	   'nup|noupdate' => \$parameter{noupdate},
	   'sp|select_programs:s' => \@{ $parameter{select_programs} },  #Comma sep string
	   'h|help' => sub { print STDOUT $USAGE, "\n"; exit;},  #Display help text
	   'v|version' => sub { print STDOUT "\ninstall.pl ".$install_version, "\n\n"; exit;},  #Display version number
    ) or help({USAGE => $USAGE,
	       exit_code => 1,
	      });

## Update default parameter dependent on other parameters
if (! $parameter{vep_cache_dir}) {

    $parameter{vep_cache_dir} = catdir($parameter{conda_path}, "envs", $parameter{conda_environment}, "ensembl-tools-release-".$parameter{varianteffectpredictor}, "cache");  #Cache directory;)
}

## Set default for array parameters
set_default_array_parameters({parameter_href => \%parameter,
			     });

##########
###MAIN###
##########


#my $LOGFILEHANDLE = &OpenLogFile({file_name => "MIP_installation.log",
#				 });

## Create bash file for writing install instructions
my $BASHFILEHANDLE = create_bash_file({file_name => "mip.sh",
				      });

## Check existance of conda environment
check_conda({parameter_href => \%parameter,
	     FILEHANDLE => $BASHFILEHANDLE,
	    });


## Create Conda environment if required
create_conda_environment({parameter_href => \%parameter,
			  FILEHANDLE => $BASHFILEHANDLE,
			 });


## Install modules into Conda environment using channel Bioconda
install_bioconda_modules({parameter_href => \%parameter,
			  FILEHANDLE => $BASHFILEHANDLE,
			 });

if (@{ $parameter{select_programs} }) {

    if ( ( grep {$_ eq "perl"} @{ $parameter{select_programs} } ) ) { #If element is part of array

	perl({parameter_href => \%parameter,
	      FILEHANDLE => $BASHFILEHANDLE,
	     });
    }
}
else {

    perl({parameter_href => \%parameter,
	  FILEHANDLE => $BASHFILEHANDLE,
	 });
}


pip_install({parameter_href => \%parameter,
	     FILEHANDLE => $BASHFILEHANDLE,
	    });

if ($parameter{prefer_shell}) {

    if (@{ $parameter{select_programs} }) {

	if ( ( grep {$_ eq "picardtools"} @{ $parameter{select_programs} } ) ) { #If element is part of array

	    picardtools({parameter_href => \%parameter,
			 FILEHANDLE => $BASHFILEHANDLE,
			});
	}
	if ( ( grep {$_ eq "sambamba"} @{ $parameter{select_programs} } ) ) { #If element is part of array

	    sambamba({parameter_href => \%parameter,
		      FILEHANDLE => $BASHFILEHANDLE,
		     });
	}
	if ( ( grep {$_ eq "bedtools"} @{ $parameter{select_programs} } ) ) { #If element is part of array

	    bedtools({parameter_href => \%parameter,
		      FILEHANDLE => $BASHFILEHANDLE,
		     });
	}
	if ( ( grep {$_ eq "vt"} @{ $parameter{select_programs} } ) ) { #If element is part of array

	    vt({parameter_href => \%parameter,
		FILEHANDLE => $BASHFILEHANDLE,
	       });
	}
	if ( ( grep {$_ eq "snpeff"} @{ $parameter{select_programs} } ) ) { #If element is part of array

	    snpeff({parameter_href => \%parameter,
		    FILEHANDLE => $BASHFILEHANDLE,
		   });
	}
	if ( ( grep {$_ eq "plink2"} @{ $parameter{select_programs} } ) ) { #If element is part of array

	    plink2({parameter_href => \%parameter,
		    FILEHANDLE => $BASHFILEHANDLE,
		   });
	}
	if ( ( grep {$_ eq "rhocall"} @{ $parameter{select_programs} } ) ) { #If element is part of array

	    rhocall({parameter_href => \%parameter,
		     FILEHANDLE => $BASHFILEHANDLE,
		    });
	}
    }
    else {

	picardtools({parameter_href => \%parameter,
		     FILEHANDLE => $BASHFILEHANDLE,
		    });

	sambamba({parameter_href => \%parameter,
		  FILEHANDLE => $BASHFILEHANDLE,
		 });

	bedtools({parameter_href => \%parameter,
		  FILEHANDLE => $BASHFILEHANDLE,
		 });

	vt({parameter_href => \%parameter,
	    FILEHANDLE => $BASHFILEHANDLE,
	   });

	snpeff({parameter_href => \%parameter,
		FILEHANDLE => $BASHFILEHANDLE,
	       });

	plink2({parameter_href => \%parameter,
		FILEHANDLE => $BASHFILEHANDLE,
	       });

	rhocall({parameter_href => \%parameter,
		 FILEHANDLE => $BASHFILEHANDLE,
		});
    }
}

if (@{ $parameter{select_programs} }) {

    if ( ( grep {$_ eq "vcftools"} @{ $parameter{select_programs} } ) ) { #If element is part of array

	vcftools({parameter_href => \%parameter,
		  FILEHANDLE => $BASHFILEHANDLE,
		 });
    }
    if ( ( grep {$_ eq "mip_scripts"} @{ $parameter{select_programs} } ) ) { #If element is part of array

	mip_scripts({parameter_href => \%parameter,
		     FILEHANDLE => $BASHFILEHANDLE,
		    });
    }
    if ( ( grep {$_ eq "varianteffectpredictor"} @{ $parameter{select_programs} } ) ) { #If element is part of array

	varianteffectpredictor({parameter_href => \%parameter,
				FILEHANDLE => $BASHFILEHANDLE,
			       });
    }
#    if ( ( grep {$_ eq "cnvnator"} @{ $parameter{select_programs} } ) ) { #If element is part of array

#	cnvnator({parameter_href => \%parameter,
#		   FILEHANDLE => $BASHFILEHANDLE,
#		  });
#   }
#    if ( ( grep {$_ eq "findtranslocations"} @{ $parameter{select_programs} } ) ) { #If element is part of array

#	findtranslocations({parameter_href => \%parameter,
#			     FILEHANDLE => $BASHFILEHANDLE,
#			    });
#    }
}
else {

    vcftools({parameter_href => \%parameter,
	      FILEHANDLE => $BASHFILEHANDLE,
	     });

    mip_scripts({parameter_href => \%parameter,
		 FILEHANDLE => $BASHFILEHANDLE,
		});

    varianteffectpredictor({parameter_href => \%parameter,
			    FILEHANDLE => $BASHFILEHANDLE,
			   });

#    cnvnator({parameter_href => \%parameter,
#	       FILEHANDLE => $BASHFILEHANDLE,
#	      });

#    findtranslocations({parameter_href => \%parameter,
#			 FILEHANDLE => $BASHFILEHANDLE,
#			});
}

close($BASHFILEHANDLE);
#close($LOGFILEHANDLE);

###SubRoutines###

sub set_default_array_parameters {

##set_default_array_parameters

##Function : Set default for array parameters
##Returns  : ""
##Arguments: $parameter_href
##         : $parameter_href => Holds all parameters

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my %array_parameter;
    $array_parameter{vep_assemblies}{default} = ["GRCh37"];
    $array_parameter{snpeff_genome_versions}{default} = ["GRCh37.75"];  #GRCh38.82 but check current on the snpEff sourceForge
    $array_parameter{perl_modules}{default} = ["Modern::Perl",  #MIP
					       "IPC::System::Simple",  #MIP
					       "Path::Iterator::Rule",  #MIP
					       "YAML",  #MIP
					       "Log::Log4perl",  #MIP
					       "List::Util",  #MIP
					       "Set::IntervalTree",  # MIP/vcfParser.pl
					       "Net::SSLeay",  # VEP
					       "LWP::Simple",  # VEP
					       "LWP::Protocol::https",  # VEP
					       "PerlIO::gzip",  #VEP
                                               "IO::Uncompress::Gunzip",  #VEP
                                               "HTML::Lint",  #VEP
                                               "Archive::Zip",  # VEP
					       "Archive::Extract",  #VEP
					       "DBI",  # VEP
					       "JSON",  # VEP
					       "DBD::mysql",  # VEP
					       "CGI",  # VEP
					       "Sereal::Encoder",  # VEP
					       "Sereal::Decoder",  # VEP
					       "Bio::Root::Version",  #VEP
					       "Module::Build", #VEP
                                               "File::Copy::Recursive", #VEP
	];

    foreach my $parameter_name (keys %array_parameter) {

	if (! @{ $parameter_href->{$parameter_name} }) {  #Unless parameter was supplied on cmd

	    $parameter_href->{$parameter_name} = $array_parameter{$parameter_name}{default};
	}
    }
}

sub create_bash_file {

##create_bash_file

##Function : Create bash file for writing install instructions
##Returns  : ""
##Arguments: $file_name, install_directory
##         : $file_name        => File name
##         : install_directory => The temporary installation directory 


    my ($arg_href) = @_;

    ## Default(s)
    my $install_directory;

    ## Flatten argument(s)
    my $file_name;

    my $tmpl = {
	file_name => { required => 1, defined => 1, strict_type => 1, store => \$file_name},
	install_directory => { default => ".MIP",
			       allow => qr/^\.\S+$/,
			       strict_type => 1, store => \$install_directory},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $FILEHANDLE = IO::Handle->new();  #Create anonymous filehandle
    my $pwd = cwd();

    ## Open batch file
    open ($FILEHANDLE, ">", catfile($pwd, $file_name)) or die("Cannot write to '".catfile($pwd, $file_name)."' :".$!."\n");

    print $FILEHANDLE "#!".catfile( dirname( dirname( devnull() ) ) ).catfile("usr", "bin", "env", "bash"), "\n\n";

    ## Create housekeeping function and trap
    say $FILEHANDLE q?finish() {?, "\n";
    say $FILEHANDLE "\t".q?## Perform exit housekeeping?;
    say $FILEHANDLE "\t".q?rm -rf ?.$install_directory;

    say $FILEHANDLE q?}?;
    say $FILEHANDLE q?trap finish EXIT TERM INT?, "\n";

    ## Create error handling function and trap
    say $FILEHANDLE q?error() {?, "\n";
    say $FILEHANDLE "\t".q?## Display error message and exit?;
    say $FILEHANDLE "\t".q{ret="$?"};
    say $FILEHANDLE "\t".q?echo "${PROGNAME}: ${1:-"Unknown Error - ExitCode="$ret}" 1>&2?, "\n";
    say $FILEHANDLE "\t".q?exit 1?;

    say $FILEHANDLE q?}?;
    say $FILEHANDLE q?trap error ERR?, "\n";

    print STDOUT "Will write install instructions to '".catfile($pwd, $file_name), "'\n";

    return $FILEHANDLE;
}

sub OpenLogFile {

##OpenLogFile

##Function : Open log file
##Returns  : ""
##Arguments: $file_name
##         : $file_name => File name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_name;

    my $tmpl = {
	file_name => { required => 1, defined => 1, strict_type => 1, store => \$file_name},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $FILEHANDLE = IO::Handle->new();  #Create anonymous filehandle

    ## Open batch file
    open ($FILEHANDLE, ">", catfile($file_name)) or die("Cannot write to '".catfile($file_name)."' :".$!."\n");

    return $FILEHANDLE;
}

sub print_parameters {

##print_parameters

##Function : Print all parameters and the default values
##Returns  : ""
##Arguments: $parameter_href
##         : $parameter_href => Holds all parameters

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Set default for array parameters
    set_default_array_parameters({parameter_href => $parameter_href,
				 });

    foreach my $key (keys %{$parameter_href}) {

	if (ref($parameter_href->{$key})!~/ARRAY|HASH/) {

	    print STDOUT $key." ";
	    if ($parameter_href->{$key}) {

		print $parameter_href->{$key}, "\n";
	    }
	    else {  ##Boolean value

		print "0", "\n";
	    }
	}
	elsif (ref($parameter_href->{$key})=~/HASH/) {

	    foreach my $program (keys %{$parameter_href->{$key}}) {

		print STDOUT $key." ".$program.": ".$parameter_href->{$key}{$program}, "\n";
	    }
	}
	elsif (ref($parameter_href->{$key})=~/ARRAY/)  {

	    print STDOUT $key.": ".join(" ", @{$parameter_href->{$key}}), "\n";
	}
    }
}


sub check_conda {

##check_conda

##Function : Check existance of conda environment
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $program = "conda";

    if ($ENV{PATH}=~/conda/) {

	print STDERR "Program check: ".$program." installed\n";
    }
    else {

	print STDERR "Could not detect ".$program." in your PATH\n";
	exit 1;
    }

    ##Deactivate any activate env prior to installation
    my $detect_active_conda_env = q?perl -nae 'if( ($_!~/^root/) && ($_=~/\*/) ) {print $F[0]}'?;
    my $ret = `conda info --envs | $detect_active_conda_env`;

    if ($ret) {

	print STDOUT "Found activated conda env: ".$ret."\n";
	print STDOUT "Please exit conda env: ".$ret." before executing install script\n";
	exit 1;
    }

    ## Check Conda path
    if (! -d $parameter_href->{conda_path}) {

	print STDERR "Could not find miniconda directory in: ".catdir($parameter_href->{conda_path}), "\n";
	exit 1;
    }

    print STDERR "Writting install instructions for Conda packages\n";

    ## Update Conda
    if ($parameter_href->{conda_update}) {

	print $FILEHANDLE "### Update Conda\n";
	print $FILEHANDLE "conda update -y conda ";
	print $FILEHANDLE "\n\n";
    }
}

sub create_conda_environment {

##create_conda_environment

##Function : Create Conda environment if required
##Returns  : ""
##Arguments: $parameter_href
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Check Conda environment
    if (! -d catdir($parameter{conda_path}, "envs", $parameter{conda_environment})) {

	## Create conda environment
	print $FILEHANDLE "### Creating Conda Environment and install: ".$parameter_href->{conda_environment}, "\n";
	print $FILEHANDLE "conda create -n ".$parameter_href->{conda_environment}." ";
	print $FILEHANDLE "-y ";
	print $FILEHANDLE "pip ";
	print $FILEHANDLE "python=".$parameter_href->{python_version}." ";
	print $FILEHANDLE "\n\n";
    }
}

sub install_bioconda_modules {

##install_bioconda_modules

##Function : Install modules into Conda environment using channel Bioconda
##Returns  : ""
##Arguments: $parameter_href
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Install into conda environment using bioconda channel
    print $FILEHANDLE "### Installing into Conda Environment: ".$parameter_href->{conda_environment}, "\n";
    print $FILEHANDLE "conda install ";
    print $FILEHANDLE "-n ".$parameter_href->{conda_environment}." ";
    print $FILEHANDLE "-y ";
    print $FILEHANDLE "-c bioconda ";

    ## Install all bioconda packages
    foreach my $program (keys %{$parameter_href->{bioconda}}) {

	print $FILEHANDLE $program."=".$parameter_href->{bioconda}{$program}." ";
    }

    print $FILEHANDLE "\n\n";

    ## Custom
    foreach my $program (keys %{$parameter_href->{bioconda}}) {

	if ($program eq "bwakit") {

	    ## Define binaries
	    my @bwakit_binaries = ("k8",
				   "seqtk",
				   "bwa-postalt.js",
				   "run-HLA",
				   "typeHLA.sh",
				   "fermi2",
				   "fermi2.pl",
				   "ropebwt2",
				   "typeHLA-selctg.js",
				   "typeHLA.js"
		);

	    foreach my $binary (@bwakit_binaries) {

		create_softlink({parameter_href => $parameter_href,
				 FILEHANDLE => $BASHFILEHANDLE,
				 binary => catfile($parameter{conda_path}, "envs", $parameter{conda_environment}, "share", "bwakit-".$parameter_href->{bioconda}{bwakit}.$parameter_href->{bioconda_bwakit_patch}, $binary),
				 softlink => $binary,
				});
	    }

	    print $BASHFILEHANDLE "cp -rf ".catdir($parameter{conda_path}, "envs", $parameter{conda_environment}, "share", "bwakit-".$parameter_href->{bioconda}{bwakit}.$parameter_href->{bioconda_bwakit_patch}, "resource-human-HLA")." ";
	    print $BASHFILEHANDLE catdir($parameter{conda_path}, "envs", $parameter{conda_environment}, "bin"), "\n\n";
	}
	if ($program eq "picard") {

	    create_softlink({parameter_href => $parameter_href,
			     FILEHANDLE => $BASHFILEHANDLE,
			     binary => catfile($parameter{conda_path}, "envs", $parameter{conda_environment}, "share", "picard-".$parameter_href->{bioconda}{picard}.$parameter_href->{bioconda_picard_patch}, "picard.jar"),
			     softlink => "picard.jar",
			    });
	}
	if ($program eq "snpeff") {

	    ## Define binaries
	    my @snpeff_binaries = ("snpEff.jar",
				   "SnpSift.jar",
				   "snpEff.config",
		);

	    foreach my $binary (@snpeff_binaries) {

		create_softlink({parameter_href => $parameter_href,
				 FILEHANDLE => $BASHFILEHANDLE,
				 binary => catfile($parameter{conda_path}, "envs", $parameter{conda_environment}, "share", "snpeff-".$parameter_href->{bioconda}{snpeff}.$parameter_href->{bioconda_snpeff_patch}, $binary),
				 softlink => $binary,
				});
	    }

	    foreach my $genome_version (@{ $parameter_href->{snpeff_genome_versions} }) {

		## Check and if required add the vertebrate mitochondrial codon table to snpeff config
		check_mt_codon_table({parameter_href => $parameter_href,
				      FILEHANDLE => $BASHFILEHANDLE,
				      share_dir => catdir($parameter{conda_path}, "envs", $parameter_href->{conda_environment}, "share", "snpeff-".$parameter_href->{bioconda}{snpeff}.$parameter_href->{bioconda_snpeff_patch}),
				      config_file => "snpEff.config",
				      genome_version_ref => \$genome_version,
				     });

		unless (-d catdir($parameter{conda_path}, "envs", $parameter_href->{conda_environment}, "share", "snpeff-".$parameter_href->{bioconda}{snpeff}.$parameter_href->{bioconda_snpeff_patch}, "data", $genome_version)) {

		    ## Write instructions to download snpeff database. This is done by install script to avoid race conditin when doing first analysis run in MIP
		    snpeff_download({parameter_href => $parameter_href,
				     FILEHANDLE => $BASHFILEHANDLE,
				     genome_version_ref => \$genome_version,
				    });
		}
	    }
	}
	if ($program eq "manta") {

	    my @manta_binaries = ("configManta.py",
				  "configManta.py.ini",
		);

	    foreach my $binary (@manta_binaries) {

		create_softlink({parameter_href => $parameter_href,
				 FILEHANDLE => $BASHFILEHANDLE,
				 binary => catfile($parameter{conda_path}, "envs", $parameter{conda_environment}, "share", "manta-".$parameter_href->{bioconda}{manta}.$parameter_href->{bioconda_manta_patch}, "bin", $binary),
				 softlink => $binary,
				});
	    }

	    ## Make file executable
	    enable_executable({parameter_href => $parameter_href,
			       FILEHANDLE => $BASHFILEHANDLE,
			       binary => q?configManta.py?,
			      });
	}
    }
}


sub perl {

##perl

##Function : Installs perl
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();

    if ($ENV{PATH}=~/perl-$parameter_href->{perl_version}/) {

	if ($parameter_href->{noupdate}) {

	    print STDERR "Found perl-".$parameter_href->{perl_version}.". in your path", "\n";
	    print STDERR q?Skipping writting installation for perl-?.$parameter_href->{perl_version}, "\n";
	}
	else {

	    if ($parameter_href->{perl_install}) {

		## Removing specific Perl version
		print $FILEHANDLE "### Removing specific perl version\n";
		print $FILEHANDLE q?rm -rf $HOME/perl-?.$parameter_href->{perl_version};
		print $FILEHANDLE "\n\n";

		install_perl_cpnam({parameter_href => $parameter_href,
				    FILEHANDLE => $BASHFILEHANDLE,
				   });
	    }

	    perl_modules({parameter_href => $parameter_href,
			  FILEHANDLE => $BASHFILEHANDLE,
			 });
	}
    }
    else {

	if ($parameter_href->{perl_install}) {

	    install_perl_cpnam({parameter_href => $parameter_href,
				FILEHANDLE => $BASHFILEHANDLE,
				path => 1,
			       });
	}

	perl_modules({parameter_href => $parameter_href,
		      FILEHANDLE => $BASHFILEHANDLE,
		     });
    }
}


sub install_perl_cpnam {

##install_perl_cpnam

##Function : Install perl CPANM
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to
##         : $path           => Export path if provided {Optional}

    my ($arg_href) = @_;

    ## Default(s)
    my $path;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
	path => { default => 0,
		  allow => [0, 1],
		  strict_type => 1, store => \$path},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();

    print STDERR "Writting install instructions for perl and Cpanm\n";

    ## Install specific perl version
    print $FILEHANDLE "### Install specific perl version\n";

    ## Move to Home
    print $FILEHANDLE "## Move HOME\n";
    print $FILEHANDLE q?cd $HOME?;
    print $FILEHANDLE "\n\n";

    ## Download
    print $FILEHANDLE "## Download perl\n";
    print $FILEHANDLE "wget --quiet http://www.cpan.org/src/5.0/perl-".$parameter_href->{perl_version}.".tar.gz ";
    print $FILEHANDLE "-O perl-".$parameter_href->{perl_version}.".tar.gz";  #Dowload outfile
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "tar xzf perl-".$parameter_href->{perl_version}.".tar.gz";
    print $FILEHANDLE "\n\n";

    ## Move to perl directory
    print $FILEHANDLE "## Move to perl directory\n";
    print $FILEHANDLE "cd perl-".$parameter_href->{perl_version};
    print $FILEHANDLE "\n\n";

    ## Configure
    print $FILEHANDLE "## Configure\n";
    print $FILEHANDLE q?./Configure -des -Dprefix=$HOME/perl-?.$parameter_href->{perl_version}, "\n";
    print $FILEHANDLE "make", "\n";

    if (! $parameter{perl_skip_test}) {

	print $FILEHANDLE "make test", "\n";
    }
    print $FILEHANDLE "make install", "\n\n";

    if ($path) {

	## Export path
	print $FILEHANDLE "## Export path\n";
	print $FILEHANDLE q?echo 'export PATH=$HOME/perl-?.$parameter_href->{perl_version}.q?/:$PATH' >> ~/.bashrc?;
	print $FILEHANDLE "\n\n";
	print $FILEHANDLE q?export PATH=$HOME/perl-?.$parameter_href->{perl_version}.q?/:$PATH?;  #Use newly installed perl
	print $FILEHANDLE "\n\n";
    }

    ## Remove tar file
    print $FILEHANDLE "## Remove tar file\n";
    print $FILEHANDLE "cd && rm perl-".$parameter_href->{perl_version}.".tar.gz";
    print $FILEHANDLE "\n\n";

    ## Move to back
    print $FILEHANDLE "## Move to original working directory\n";
    print $FILEHANDLE "cd ".$pwd;
    print $FILEHANDLE "\n\n";

    print $FILEHANDLE q?echo 'eval `perl -I ~/perl-?.$parameter_href->{perl_version}.q?/lib/perl5/ -Mlocal::lib=~/perl-?.$parameter_href->{perl_version}.q?/`' >> ~/.bash_profile ?;  #Add at start-up
    print $FILEHANDLE "\n\n";
    print $FILEHANDLE q?echo 'export PERL_UNICODE=SAD' >> ~/.bash_profile ?;  #Add at start-up
    print $FILEHANDLE "\n\n";

    ## Install perl modules via cpanm
    print $FILEHANDLE "## Install cpanm\n";
    print $FILEHANDLE q?wget -O- http://cpanmin.us | perl - -l $HOME/perl-?.$parameter_href->{perl_version}.q?/bin App::cpanminus --local-lib=~/perl-?.$parameter_href->{perl_version}.q?/ local::lib ?;
    print $FILEHANDLE "\n\n";

    ## Use newly installed perl
    print $FILEHANDLE q?eval `perl -I ~/perl-?.$parameter_href->{perl_version}.q?/lib/perl5/ -Mlocal::lib=~/perl-?.$parameter_href->{perl_version}.q?/` ?;
    print $FILEHANDLE "\n\n";

    ## Use newly installed perl
    print $FILEHANDLE q?PERL5LIB=~/perl-?.$parameter_href->{perl_version}.q?/lib/perl5?;
    print $FILEHANDLE "\n\n";
}


sub perl_modules {

##perl_modules

##Function : Install perl modules via cpanm
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Install perl modules via cpanm
    print $FILEHANDLE "## Install perl modules via cpanm\n";
    print $FILEHANDLE "cpanm ";

    if ($parameter{perl_modules_force}) {

	print $FILEHANDLE "--force ";
    }
    print $FILEHANDLE join(" ", @{ $parameter_href->{perl_modules} })." ";
    print $FILEHANDLE "\n\n";
}


sub pip_install {

##pip_install

##Function : Writes install instructions for pip packages
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    print STDERR "Writting install instructions for pip packages\n";

    ## Install PIP packages in conda environment
    print $FILEHANDLE "### Install PIP packages in conda environment: ".$parameter_href->{conda_environment}, "\n";

    ## Activate conda environment
    activate_conda_environment({parameter_href => $parameter_href,
				FILEHANDLE => $FILEHANDLE,
			       });

    ## Install PIP packages
    print $FILEHANDLE "## Install PIP packages\n";
    print $FILEHANDLE "pip install ";

    ## Install all PIP packages
    foreach my $program (keys %{ $parameter_href->{pip} }) {

	print $FILEHANDLE $program."==".$parameter_href->{pip}{$program}." ";
    }
    print $FILEHANDLE "\n\n";

    ## Deactivate conda environment
    deactivate_conda_environment({FILEHANDLE => $FILEHANDLE,
				 });
}


sub picardtools {

##picardtools

##Function : Install picardtools
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Check if the binary of the program being installed already exists
    if (check_conda_bin_file_exists({parameter_href => $parameter_href,
				     program_name => "picard.jar",
				    })) {  # Assumes that picard.jar is there as well then

	return
    }

    ## Install picard
    print $FILEHANDLE "### Install Picard\n";

    ## Create the temporary install directory
    create_install_dir({FILEHANDLE => $FILEHANDLE,
		       });

    ## Download
    print $FILEHANDLE "## Download Picard\n";
    print $FILEHANDLE "wget --quiet https://github.com/broadinstitute/picard/releases/download/".$parameter_href->{picardtools}."/picard-tools-".$parameter_href->{picardtools}.".zip ";
    print $FILEHANDLE "-O picard-tools-".$parameter_href->{picardtools}.".zip";  #Download outfile
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "unzip picard-tools-".$parameter_href->{picardtools}.".zip";
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    if (-d catdir($parameter{conda_path}, "envs", $parameter_href->{conda_environment}, "share", "picard-tools-".$parameter_href->{picardtools})) {

	print $FILEHANDLE "rm -rf ".catdir($parameter{conda_path}, "envs", $parameter_href->{conda_environment}, "share", "picard-tools-".$parameter_href->{picardtools});
	print $FILEHANDLE "\n\n";
    }

    print $FILEHANDLE "## Make available from conda environment\n";
    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?picard-tools-?.$parameter_href->{picardtools}.q? ?.catdir($parameter{conda_path}, "envs", $parameter_href->{conda_environment}, "share");
    print $FILEHANDLE "\n\n";

    create_softlink({parameter_href => $parameter_href,
		     FILEHANDLE => $FILEHANDLE,
		     binary => catfile($parameter{conda_path}, "envs", $parameter{conda_environment}, "share", "picard-tools-".$parameter_href->{picardtools}, "picard.jar"),
		     softlink => "picard.jar",
		    });

    ## Remove the temporary install directory
    remove_install_dir({FILEHANDLE => $FILEHANDLE,
			pwd => $pwd,
		       });
}


sub sambamba {

##sambamba

##Function : Install sambamba
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Check if the binary of the program being installed already exists
    if (check_conda_bin_file_exists({parameter_href => $parameter_href,
				     program_name => "sambamba",
				     program_version => $parameter_href->{sambamba},
				    }) ) {
	return
    }

    ## Install sambamba
    print $FILEHANDLE "### Install sambamba\n";

    ## Create the temporary install directory
    create_install_dir({FILEHANDLE => $FILEHANDLE,
		       });

    ## Download
    print $FILEHANDLE "## Download sambamba release\n";
    print $FILEHANDLE q?wget --quiet https://github.com/lomereiter/sambamba/releases/download/v?.$parameter_href->{sambamba}.q?/sambamba_v?.$parameter_href->{sambamba}.q?_linux.tar.bz2 ?;
    print $FILEHANDLE "-O sambamba_v".$parameter_href->{sambamba}."_linux.tar.bz2";  #Download outfile
    print $FILEHANDLE "\n\n";

    ## Decompress
    print $FILEHANDLE "## Decompress sambamba file\n";
    print $FILEHANDLE "bzip2 ";
    print $FILEHANDLE "-f ";  #Force
    print $FILEHANDLE "-d ";  #Decompress
    print $FILEHANDLE "sambamba_v".$parameter_href->{sambamba}."_linux.tar.bz2";
    print $FILEHANDLE "\n\n";

    ## Extract files
    print $FILEHANDLE "## Extract files\n";
    print $FILEHANDLE "tar xvf sambamba_v".$parameter_href->{sambamba}."_linux.tar";
    print $FILEHANDLE "\n\n";

    ## Make executable
    print $FILEHANDLE "## Make executable\n";
    print $FILEHANDLE "chmod 755 ";
    print $FILEHANDLE "sambamba_v".$parameter_href->{sambamba};
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    print $FILEHANDLE "## Make available from conda environment\n";
    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?sambamba_v?.$parameter_href->{sambamba}.q? ?.$parameter{conda_path}.q?/envs/?.$parameter_href->{conda_environment}.q?/bin/?;
    print $FILEHANDLE "\n\n";

    create_softlink({parameter_href => $parameter_href,
		     FILEHANDLE => $BASHFILEHANDLE,
		     binary => "sambamba_v".$parameter_href->{bioconda}{sambamba},
		     softlink => "sambamba",
		    });

    ## Remove the temporary install directory
    remove_install_dir({FILEHANDLE => $FILEHANDLE,
			pwd => $pwd,
		       });

}


sub vcftools {

##vcftools

##Function : Install vcftools
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Check if the binary of the program being installed already exists
    if(check_conda_bin_file_exists({parameter_href => $parameter_href,
				    program_name => "vcftools",
				   })) {

	return
    }

    ## Install vcftools
    print $FILEHANDLE "### Install vcftools\n";

    ## Create the temporary install directory
    create_install_dir({FILEHANDLE => $FILEHANDLE,
		       });

    ## Download
    print $FILEHANDLE "## Download vcftools\n";
    print $FILEHANDLE "wget --quiet https://github.com/vcftools/vcftools/releases/download/v".$parameter_href->{vcftools}."/vcftools-".$parameter_href->{vcftools}.".tar.gz ";
    print $FILEHANDLE "-O vcftools-".$parameter_href->{vcftools}.".tar.gz";  #Download outfile
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "tar xvf vcftools-".$parameter_href->{vcftools}.".tar.gz";
    print $FILEHANDLE "\n\n";

    ## Export PERL5LIB environment variable
    print $FILEHANDLE "## Export PERL5LIB environment variable\n";
    print $FILEHANDLE q?export PERL5LIB=?.$Bin.q?/vcftools-?.$parameter_href->{vcftools}.q?/src/perl/?;
    print $FILEHANDLE "\n\n";

    ## Move to vcftools directory
    print $FILEHANDLE "## Move to vcftools directory\n";
    print $FILEHANDLE "cd vcftools-".$parameter_href->{vcftools};
    print $FILEHANDLE "\n\n";

    ## Configure
    my $filePath = $parameter{conda_path}.q?/envs/?.$parameter_href->{conda_environment};

    print $FILEHANDLE "## Configure", "\n";
    print $FILEHANDLE q?./configure --prefix=?.$filePath, "\n";
    print $FILEHANDLE "make", "\n";
    print $FILEHANDLE "make install", "\n";
    print $FILEHANDLE "\n\n";

    ## Move perl Module
    print $FILEHANDLE "## Move perl Module\n";
    print $FILEHANDLE q?cp src/perl/Vcf.pm $HOME/perl-?.$parameter_href->{perl_version}.q?/lib/perl5/?;
    print $FILEHANDLE "\n\n";

    ## Remove the temporary install directory
    remove_install_dir({FILEHANDLE => $FILEHANDLE,
			pwd => $pwd,
		       });

    ## Reset perl envionment
    print $FILEHANDLE q?PERL5LIB=~/perl-?.$parameter_href->{perl_version}.q?/lib/perl5?;
    print $FILEHANDLE "\n\n";
}


sub bedtools {

##bedtools

##Function : Install bedtools
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Check if the binary of the program being installed already exists
    if(check_conda_bin_file_exists({parameter_href => $parameter_href,
				    program_name => "bedtools",
				   })) {

	return
    }

    my $bedtools_main_version = substr($parameter_href->{bedtools}, 0, 1);

    ## Install bedtools
    print $FILEHANDLE "### Install bedtools\n";

    ## Create the temporary install directory
    create_install_dir({FILEHANDLE => $FILEHANDLE,
		       });

    ## Download
    print $FILEHANDLE "## Download bedtools\n";
    print $FILEHANDLE "wget --quiet https://github.com/arq5x/bedtools".$bedtools_main_version."/releases/download/v".$parameter_href->{bedtools}."/bedtools-".$parameter_href->{bedtools}.".tar.gz ";
    print $FILEHANDLE "-O bedtools-".$parameter_href->{bedtools}.".tar.gz";  #Download outfile
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "tar xvf bedtools-".$parameter_href->{bedtools}.".tar.gz";
    print $FILEHANDLE "\n\n";

    ## Move to bedtools directory
    print $FILEHANDLE "## Move to bedtools directory\n";
    print $FILEHANDLE "cd bedtools".$bedtools_main_version;
    print $FILEHANDLE "\n\n";

    print $FILEHANDLE "make";
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    print $FILEHANDLE "## Make available from conda environment\n";
    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?./bin/* ?.$parameter{conda_path}.q?/envs/?.$parameter_href->{conda_environment}.q?/bin/?;
    print $FILEHANDLE "\n\n";

    ## Remove the temporary install directory
    remove_install_dir({FILEHANDLE => $FILEHANDLE,
			pwd => $pwd,
		       });
}


sub vt {

##vt

##Function : Install vt
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Check if the binary of the program being installed already exists
    if (check_conda_bin_file_exists({parameter_href => $parameter_href,
				     program_name => "vt",
				    })) {

	return
    }

    ## Install vt
    print $FILEHANDLE "### Install VT\n";

    ## Create the temporary install directory
    create_install_dir({FILEHANDLE => $FILEHANDLE,
		       });

    ## Download
    print $FILEHANDLE "## Download VT\n";

    print $FILEHANDLE "git clone https://github.com/atks/vt.git ";
    print $FILEHANDLE "\n\n";

    ## Move to vt directory
    print $FILEHANDLE "## Move to vt directory\n";
    print $FILEHANDLE "cd vt ";
    print $FILEHANDLE "\n\n";

    ## Configure
    print $FILEHANDLE "## Configure\n";
    print $FILEHANDLE "make", "\n";
    print $FILEHANDLE "make test";
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    print $FILEHANDLE "## Make available from conda environment\n";
    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?vt ?.$parameter{conda_path}.q?/envs/?.$parameter_href->{conda_environment}.q?/bin/?;
    print $FILEHANDLE "\n\n";

    ## Remove the temporary install directory
    remove_install_dir({FILEHANDLE => $FILEHANDLE,
			pwd => $pwd,
		       });
}


sub plink2 {

##plink2

##Function : Install plink2
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Check if the binary of the program being installed already exists
    if (check_conda_bin_file_exists({parameter_href => $parameter_href,
				     program_name => "plink",
				    })) {

	return
    }

    ## Install Plink
    print $FILEHANDLE "### Install Plink\n";

    ## Create the temporary install directory
    create_install_dir({FILEHANDLE => $FILEHANDLE,
		       });

    ## Download
    print $FILEHANDLE "## Download Plink\n";
    print $FILEHANDLE "wget --quiet https://www.cog-genomics.org/static/bin/plink".$parameter_href->{plink2}."/plink_linux_x86_64.zip ";
    print $FILEHANDLE "-O plink-".$parameter_href->{plink2}."-x86_64.zip";  #Download outfile
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "unzip plink-".$parameter_href->{plink2}."-x86_64.zip";
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    print $FILEHANDLE "## Make available from conda environment\n";
    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?plink ?.$parameter{conda_path}.q?/envs/?.$parameter_href->{conda_environment}.q?/bin/plink2?;
    print $FILEHANDLE "\n\n";

    ## Remove the temporary install directory
    remove_install_dir({FILEHANDLE => $FILEHANDLE,
			pwd => $pwd,
		       });
}


sub snpeff {

##snpeff

##Function : Install snpeff
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Check if the binary of the program being installed already exists
    if (check_conda_bin_file_exists({parameter_href => $parameter_href,
				     program_name => "snpEff.jar",
				    })) {  # Assumes that SnpSift.jar is there as well then

	return
    }

    ## Install snpeff
    print $FILEHANDLE "### Install snpeff\n";

    ## Create the temporary install directory
    create_install_dir({FILEHANDLE => $FILEHANDLE,
		       });

    ## Download
    print $FILEHANDLE "## Download snpeff\n";
    print $FILEHANDLE "wget --quiet http://sourceforge.net/projects/snpeff/files/snpEff_".$parameter_href->{snpeff}."_core.zip/download ";
    print $FILEHANDLE "-O snpEff_".$parameter_href->{snpeff}."_core.zip";  #Download outfile
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "unzip snpEff_".$parameter_href->{snpeff}."_core.zip";
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    if (-d $parameter{conda_path}.q?/envs/?.$parameter_href->{conda_environment}.q?/share/snpEff.?.$parameter_href->{snpeff}) {

	print $FILEHANDLE "rm -rf ".$parameter{conda_path}.q?/envs/?.$parameter_href->{conda_environment}.q?/share/snpEff.?.$parameter_href->{snpeff};
	print $FILEHANDLE "\n\n";
    }

    print $FILEHANDLE "## Make available from conda environment\n";
    print $FILEHANDLE "mkdir -p ".$parameter{conda_path}.q?/envs/?.$parameter_href->{conda_environment}.q?/share/snpEff.?.$parameter_href->{snpeff};
    print $FILEHANDLE "\n\n";

    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?snpEff/*.jar ?.$parameter{conda_path}.q?/envs/?.$parameter_href->{conda_environment}.q?/share/snpEff.?.$parameter_href->{snpeff}.q?/?;
    print $FILEHANDLE "\n\n";

    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?snpEff/snpEff.config ?.$parameter{conda_path}.q?/envs/?.$parameter_href->{conda_environment}.q?/share/snpEff.?.$parameter_href->{snpeff}.q?/?;
    print $FILEHANDLE "\n\n";

    ## Define binaries
    my @snpeff_binaries = ("snpEff.jar",
			   "SnpSift.jar",
			   "snpEff.config",
	);

    foreach my $binary (@snpeff_binaries) {

	create_softlink({parameter_href => $parameter_href,
			 FILEHANDLE => $BASHFILEHANDLE,
			 binary => catfile($parameter{conda_path}, "envs", $parameter{conda_environment}, "share", "snpEff.".$parameter_href->{snpeff}, $binary),
			 softlink => $binary,
			});
    }

    foreach my $genome_version (@{$parameter_href->{snpeff_genome_versions}}) {

	## Check and if required add the vertebrate mitochondrial codon table to snpeff config
	check_mt_codon_table({parameter_href => $parameter_href,
			      FILEHANDLE => $BASHFILEHANDLE,
			      share_dir => catdir($parameter{conda_path}, "envs", $parameter_href->{conda_environment}, "share", "snpEff.".$parameter_href->{snpeff}),
			      config_file => "snpEff.config",
			      genome_version_ref => \$genome_version,
			     });

	unless (-d catdir($parameter{conda_path}, "envs", $parameter_href->{conda_environment}, "share", "snpEff.".$parameter_href->{snpeff}, "data", $genome_version)) {

	    ## Write instructions to download snpeff database. This is done by install script to avoid race conditin when doing first analysis run in MIP
	    snpeff_download({parameter_href => $parameter_href,
			     FILEHANDLE => $BASHFILEHANDLE,
			     genome_version_ref => \$genome_version,
			    });
	}
    }

    ## Remove the temporary install directory
    remove_install_dir({FILEHANDLE => $FILEHANDLE,
			pwd => $pwd,
		       });
}


sub varianteffectpredictor {

##varianteffectpredictor

##Function : Install varianteffectpredictor
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();

    my $miniconda_bin_dir = catdir($parameter{conda_path}, "envs", $parameter_href->{conda_environment}, "ensembl-tools-release-".$parameter_href->{varianteffectpredictor});

    if (-d $miniconda_bin_dir) {

	print STDERR q?Found varianteffectpredictor in miniconda directory: ?.$miniconda_bin_dir, "\n";

	if ($parameter_href->{noupdate}) {

	    print STDERR "Skipping writting installation process for varianteffectpredictor\n";
	    return
	}
	else {

	    ## Removing varianteffectpredictor
	    print $FILEHANDLE "### Removing varianteffectpredictor\n";
	    print $FILEHANDLE q?rm -rf ?.$miniconda_bin_dir;
	    print $FILEHANDLE "\n\n";
	}
    }
    else {

	print STDERR "Writting install instructions for varianteffectpredictor\n";
    }

    ## Install VEP
    print $FILEHANDLE "### Install varianteffectpredictor\n";

    ## Activate conda environment
    activate_conda_environment({parameter_href => $parameter_href,
				FILEHANDLE => $FILEHANDLE,
			       });

    ##Make sure that the cache directory exists
    print $FILEHANDLE "mkdir -p ".$parameter_href->{vep_cache_dir}." ";  #Cache directory
    print $FILEHANDLE "\n\n";

    ## Move to miniconda environment
    print $FILEHANDLE q?cd ?.catdir($parameter{conda_path}, "envs", $parameter_href->{conda_environment});
    print $FILEHANDLE "\n\n";

    ## Download
    print $FILEHANDLE "## Download VEP\n";
    print $FILEHANDLE "wget --quiet https://github.com/Ensembl/ensembl-tools/archive/release/".$parameter_href->{varianteffectpredictor}.".zip ";
    print $FILEHANDLE "-O VariantEffectPredictor-".$parameter_href->{varianteffectpredictor}.".zip";  #Download outfile
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "unzip VariantEffectPredictor-".$parameter_href->{varianteffectpredictor}.".zip";
    print $FILEHANDLE "\n\n";

    ## Move to VariantEffectPredictor directory
    print $FILEHANDLE "## Move to VariantEffectPredictor directory\n";
    print $FILEHANDLE "cd ".catdir("ensembl-tools-release-".$parameter_href->{varianteffectpredictor}, "scripts", "variant_effect_predictor");
    print $FILEHANDLE "\n\n";

    ## Install VEP
    print $FILEHANDLE "## Install VEP\n";
    print $FILEHANDLE "perl INSTALL.pl ";
    print $FILEHANDLE "--AUTO ".$parameter{vep_auto_flag};  #a (API), l (FAIDX/htslib), c (cache), f (FASTA)

    if ( (defined($parameter{vep_plugin})) && ($parameter{vep_plugin} ne 0)) {

	print $FILEHANDLE "p ";  #p (plugins)
	print $FILEHANDLE "-g ".$parameter{vep_plugin}." ";  #Plugins in comma sep string
    }
    else {

	print $FILEHANDLE " ";  #Add whitespace if no plugins
    }

    print $FILEHANDLE "-c ".$parameter_href->{vep_cache_dir}." ";  #Cache directory
    print $FILEHANDLE "-s homo_sapiens ";

    ## Only install first assembly version since VEP install cannot handle multiple versions at the same time
    print $FILEHANDLE "--ASSEMBLY ".$parameter_href->{vep_assemblies}[0]." ";
    print $FILEHANDLE "\n\n";

    if (scalar( @{$parameter_href->{vep_assemblies}} ) > 1 ) {

	for (my $assembly_version=1;$assembly_version<scalar(@{$parameter_href->{vep_assemblies}});$assembly_version++) {

	    print $FILEHANDLE "## Install additional VEP cache assembly version\n";
	    print $FILEHANDLE "perl INSTALL.pl ";
	    print $FILEHANDLE "--AUTO cf ";  #a (API), l (FAIDX/htslib), c (cache), f (FASTA), p (plugins)
	    print $FILEHANDLE "-c ".$parameter_href->{vep_cache_dir}." ";  #Cache directory
	    print $FILEHANDLE "-s homo_sapiens ";

	    ## Only install first assembly version since VEP install cannot handle multiple versions at the same time
	    print $FILEHANDLE "--ASSEMBLY ".$parameter_href->{vep_assemblies}[$assembly_version]." ";
	    print $FILEHANDLE "\n\n";
	}
    }

    if ( defined($parameter{vep_plugin}) && ($parameter{vep_plugin}=~/Loftool/) ) {

	##Add LofTool required text file
	print $FILEHANDLE "##Add LofTool required text file\n";
	print $FILEHANDLE "wget https://raw.githubusercontent.com/Ensembl/VEP_plugins/master/LoFtool_scores.txt ";
	print $FILEHANDLE q?-O $HOME/.vep/Plugins/LoFtool_scores.txt ?;
	print $FILEHANDLE "\n\n";
    }

    if ( defined($parameter{vep_plugin}) && ($parameter{vep_plugin}=~/Lof/) ) {

	##Add Lof required perl splice script
	print $FILEHANDLE "##Add Lof required perl splice script\n";
	print $FILEHANDLE "wget https://raw.githubusercontent.com/konradjk/loftee/master/splice_module.pl ";
	print $FILEHANDLE q?-O $HOME/.vep/Plugins/splice_module.pl ?;
	print $FILEHANDLE "\n\n";

	##Add Lof optional human_ancestor_fa
	print $FILEHANDLE "##Add Lof optional human_ancestor_fa\n";
	print $FILEHANDLE "wget https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz ";
	print $FILEHANDLE q?-O ?.catfile($parameter_href->{vep_cache_dir}, "human_ancestor.fa.gz")." ";
	print $FILEHANDLE "\n\n";

	##Uncompress
	print $FILEHANDLE "bgzip -d ".catfile($parameter_href->{vep_cache_dir}, "human_ancestor.fa.gz")." ";
	print $FILEHANDLE "\n\n";

	print $FILEHANDLE "wget https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz.fai ";
	print $FILEHANDLE q?-O ?.catfile($parameter_href->{vep_cache_dir}, "human_ancestor.fa.fai")." ";
	print $FILEHANDLE "\n\n";
    }

    ## Clean up
    print $FILEHANDLE "## Clean up\n";
    print $FILEHANDLE q?cd ?.catdir($parameter{conda_path}, "envs", $parameter_href->{conda_environment});
    print $FILEHANDLE "\n\n";

    print $FILEHANDLE "rm -rf VariantEffectPredictor-".$parameter_href->{varianteffectpredictor}.".zip";;
    print $FILEHANDLE "\n\n";

    ## Moving up
    print $FILEHANDLE "## Moving back to original working directory\n";
    print $FILEHANDLE "cd ".$pwd;  #Go back to subroutine origin
    print $FILEHANDLE "\n\n";

    ## Deactivate conda environment
    deactivate_conda_environment({FILEHANDLE => $FILEHANDLE,
				 });
}


sub cnvnator {

##cnvnator

##Function : Install cnvnator
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Check if the binary of the program being installed already exists
    if (check_conda_bin_file_exists({parameter_href => $parameter_href,
				     program_name => "cnvnator",
				    })) {

	return
    }

    my $miniconda_bin_dir = catdir($parameter{conda_path}, "envs", $parameter_href->{conda_environment}, "root");

    if (-d $miniconda_bin_dir) {

	print STDERR q?Found Root in miniconda directory: ?.$miniconda_bin_dir, "\n";

	if ($parameter_href->{noupdate}) {

	    print STDERR "Skipping writting installation process for Root\n";
	    return
	}
	else {

	    ## Removing Root
	    print $FILEHANDLE "### Removing Root\n";
	    print $FILEHANDLE q?rm -rf ?.$miniconda_bin_dir;
	    print $FILEHANDLE "\n\n";
	}
    }
    else {

	print STDERR "Writting install instructions for Root\n";
    }

    ## Install Root
    print $FILEHANDLE "### Install cnvnator/Root\n";

    ## Move to miniconda environment
    print $FILEHANDLE q?cd ?.catdir($parameter{conda_path}, "envs", $parameter_href->{conda_environment});
    print $FILEHANDLE "\n\n";

    ## Download
    print $FILEHANDLE "## Download Root\n";

    print $FILEHANDLE "wget --quiet https://root.cern.ch/download/root_v5.34.34.Linux-slc6-x86_64-gcc4.4.tar.gz ";  #Currently hardcoded
    print $FILEHANDLE "-O root_v5.34.34.Linux-slc6-x86_64-gcc4.4.tar.gz ";  #Download outfile
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "tar xvf root_v5.34.34.Linux-slc6-x86_64-gcc4.4.tar.gz ";
    print $FILEHANDLE "\n\n";

    unless ($ENV{PATH}=~/$parameter{conda_path}\/envs\/$parameter_href->{conda_environment}\/root\/bin/) {

	## Export path
	print $FILEHANDLE "## Export path\n";
	print $FILEHANDLE q?echo 'source ?.$parameter{conda_path}.q?/envs/?.$parameter_href->{conda_environment}.q?/root/bin/thisroot.sh' >> ~/.bashrc?;
	print $FILEHANDLE "\n\n";

	## Use newly installed root
	print $FILEHANDLE q?source ?.$parameter{conda_path}.q?/envs/?.$parameter_href->{conda_environment}.q?/root/bin/thisroot.sh ?;
	print $FILEHANDLE "\n\n";
    }

    ## Moving up
    print $FILEHANDLE "## Moving back to original working directory\n";
    print $FILEHANDLE "cd ".$pwd;  #Go back to subroutine origin
    print $FILEHANDLE "\n\n";

    ## Install CNVNator
    print $FILEHANDLE "### Install cnvnator\n";

    ## Activate conda environment
    activate_conda_environment({parameter_href => $parameter_href,
				FILEHANDLE => $FILEHANDLE,
			       });

    ## Create the temporary install directory
    create_install_dir({FILEHANDLE => $FILEHANDLE,
		       });

    ## Download
    print $FILEHANDLE "## Download CNVNator\n";
    print $FILEHANDLE "wget --quiet https://github.com/abyzovlab/CNVnator/releases/download/v".$parameter_href->{cnvnator}."/CNVnator_v".$parameter_href->{cnvnator}.".zip ";
    print $FILEHANDLE "-O CNVnator_v".$parameter_href->{cnvnator}.".zip";  #Download outfile
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "unzip CNVnator_v".$parameter_href->{cnvnator}.".zip";
    print $FILEHANDLE "\n\n";

    ## Move to CNVnator directory
    print $FILEHANDLE "## Move to CNVnator directory\n";
    print $FILEHANDLE "cd ".catdir("CNVnator_v".$parameter_href->{cnvnator}, "src", "samtools");
    print $FILEHANDLE "\n\n";

    ## Configure
    print $FILEHANDLE "## Configure CNVnator samTools specific version\n";
    print $FILEHANDLE "make";
    print $FILEHANDLE "\n\n";

    print $FILEHANDLE "## Move to CNVnator directory\n";
    print $FILEHANDLE "cd ..";
    print $FILEHANDLE "\n";

    print $FILEHANDLE "make";
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    print $FILEHANDLE "## Make available from conda environment\n";
    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?cnvnator ?.catdir($parameter{conda_path}, "envs", $parameter_href->{conda_environment}, "bin");
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    print $FILEHANDLE "## Make available from conda environment\n";
    print $FILEHANDLE "cd ..";
    print $FILEHANDLE "\n";
    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?cnvnator2VCF.pl ?.catdir($parameter{conda_path}, "envs", $parameter_href->{conda_environment}, "bin");
    print $FILEHANDLE "\n\n";

    print $FILEHANDLE "## Make executable from conda environment\n";
    print $FILEHANDLE "chmod +x ";
    print $FILEHANDLE catfile($parameter{conda_path}, "envs", $parameter_href->{conda_environment}, "bin", "cnvnator2VCF.pl");
    print $FILEHANDLE "\n\n";

    ## Remove the temporary install directory
    remove_install_dir({FILEHANDLE => $FILEHANDLE,
			pwd => $pwd,
		       });
    ## Deactivate conda environment
    deactivate_conda_environment({FILEHANDLE => $FILEHANDLE,
				 });
}


sub findtranslocations {

##findtranslocations

##Function : Install findtranslocations
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Check if the binary of the program being installed already exists
    if (check_conda_bin_file_exists({parameter_href => $parameter_href,
				     program_name => "findtranslocations",
				    })) {

	return
    }

    ## Install findtranslocations
    print $FILEHANDLE "### Install findtranslocations\n";

    ## Activate conda environment
    activate_conda_environment({parameter_href => $parameter_href,
				FILEHANDLE => $FILEHANDLE,
			       });

    ## Add to bashrc
    unless (-d catdir($parameter{conda_path}, "envs", $parameter_href->{conda_environment}, "FindTranslocations", "bin")) {

	## Export path
	print $FILEHANDLE "## Export to bashrc\n";
	print $FILEHANDLE q?printf '\nif [ -f ?.$parameter{conda_path}.q?/envs/?.$parameter_href->{conda_environment}.q?/FindTranslocations/bin/FindTranslocations ]; then\n?;
	print $FILEHANDLE q?\t\texport LD_LIBRARY_PATH=$LD_LIBRARY_PATH:?.$parameter{conda_path}.q?/pkgs/boost-?.$parameter{bioconda}{boost}.$parameter{bioconda_boost_patch}.q?/lib\n?;
	print $FILEHANDLE q?fi\n\n' >> ~/.bashrc?;
	print $FILEHANDLE "\n\n";
    }

    ## Move to miniconda environment
    print $FILEHANDLE "cd ".catdir($parameter{conda_path}, "envs", $parameter_href->{conda_environment});
    print $FILEHANDLE "\n\n";

    ## Download
    print $FILEHANDLE "## Download FindTranslocations\n";
    print $FILEHANDLE "wget --quiet https://github.com/J35P312/FindTranslocations/archive/version_".$parameter_href->{findtranslocations}.".zip ";
    print $FILEHANDLE "-O FindTranslocations-".$parameter_href->{findtranslocations}.".zip";  #Download outfile
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "rm -rf FindTranslocations";
    print $FILEHANDLE "\n\n";
    print $FILEHANDLE "unzip FindTranslocations-".$parameter_href->{findtranslocations}.".zip ";
    print $FILEHANDLE "\n\n";
    print $FILEHANDLE "mv FindTranslocations-version_".$parameter_href->{findtranslocations}." ";
    print $FILEHANDLE "FindTranslocations ";
    print $FILEHANDLE "\n\n";

    ## Move to FindTranslocations directory
    print $FILEHANDLE "## Move to FindTranslocations directory\n";
    print $FILEHANDLE "cd FindTranslocations";
    print $FILEHANDLE "\n\n";
    print $FILEHANDLE "mkdir -p build";
    print $FILEHANDLE "\n\n";
    print $FILEHANDLE "cd build";
    print $FILEHANDLE "\n\n";

    ## Configure
    print $FILEHANDLE "## Configure\n";
    print $FILEHANDLE "cmake .. -DBoost_NO_BOOST_CMAKE=ON";
    print $FILEHANDLE "\n\n";

    print $FILEHANDLE "make";
    print $FILEHANDLE "\n\n";

    print $FILEHANDLE "cd ../bin";
    print $FILEHANDLE "\n\n";
    print $FILEHANDLE "chmod a+x FindTranslocations";
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    my $cwd = cwd();
    print $FILEHANDLE "## Make available from conda environment\n";
    print $FILEHANDLE "ln -f -s  ";
    print $FILEHANDLE catfile($parameter{conda_path}, "envs", $parameter_href->{conda_environment}, "FindTranslocations", "bin", "FindTranslocations")." ".catdir($parameter{conda_path}, "envs", $parameter_href->{conda_environment}, "bin");
    print $FILEHANDLE "\n\n";

    ## Clean-up
    print $FILEHANDLE "cd ".catdir($parameter{conda_path}, "envs", $parameter_href->{conda_environment});
    print $FILEHANDLE "\n\n";
    print $FILEHANDLE "rm -rf FindTranslocations-".$parameter_href->{findtranslocations}.".zip";
    print $FILEHANDLE "\n\n";

    ## Deactivate conda environment
    deactivate_conda_environment({FILEHANDLE => $FILEHANDLE,
				 });
}


sub mip_scripts {

##mip_scripts

##Function : Install mip_scripts
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Define MIP scripts and yaml files
    my @mip_scripts = ("calculate_af.pl",
		       "max_af.pl",
		       "mip.pl",
		       "qccollect.pl",
		       "vcfparser.pl",
		       "install.pl",
	);
    my %mip_sub_scripts;
    $mip_sub_scripts{"definitions"} = ["define_parameters.yaml"];
    $mip_sub_scripts{"t"} = ["test.t"];
    $mip_sub_scripts{"templates"} = ["mip_config.yaml"];

    ## Check if the binary of the program being installed already exists
    if (check_conda_bin_file_exists({parameter_href => $parameter_href,
				     program_name => "mip.pl",
				    })) {  #Proxy for all

	return
    }

    ## Install mip_scripts
    print $FILEHANDLE "### Install mip_scripts\n";

    ## Create directories
    print $FILEHANDLE "## Create directories\n";
    foreach my $directory (keys %mip_sub_scripts) {

	print $FILEHANDLE "mkdir -p ";
	print $FILEHANDLE catdir($parameter{conda_path}, "envs", $parameter_href->{conda_environment}, "bin", $directory);
	print $FILEHANDLE "\n\n";
    }

    ## Copy mip scripts and sub scripts to conda env and make executable
    print $FILEHANDLE "## Copy mip scripts and subdirectory scripts to conda env and make executable\n\n";
    foreach my $script (@mip_scripts) {

	print $FILEHANDLE "cp ";
	print $FILEHANDLE catfile($Bin, $script)." ";
	print $FILEHANDLE catdir($parameter{conda_path}, "envs", $parameter_href->{conda_environment}, "bin"), "\n";
	print $FILEHANDLE "chmod a+x ".catfile($parameter{conda_path}, "envs", $parameter_href->{conda_environment}, "bin", $script);
	print $FILEHANDLE "\n\n";
    }

    foreach my $directory (keys %mip_sub_scripts) {

	foreach my $script (@{$mip_sub_scripts{$directory}}) {

	    print $FILEHANDLE "cp ";
	    print $FILEHANDLE catfile($Bin, $directory, $script)." ";
	    print $FILEHANDLE catdir($parameter{conda_path}, "envs", $parameter_href->{conda_environment}, "bin", $directory), "\n";
	    print $FILEHANDLE "chmod a+x ".catfile($parameter{conda_path}, "envs", $parameter_href->{conda_environment}, "bin", $directory, $script);
	    print $FILEHANDLE "\n\n";
	}
    }
}


sub rhocall {

##rhocall

##Function : Install rhocall
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Check if the binary of the program being installed already exists
    if (check_conda_bin_file_exists({parameter_href => $parameter_href,
				     program_name => "rhocall",
				    })) {

	return
    }

    ## Activate conda environment
    activate_conda_environment({parameter_href => $parameter_href,
				FILEHANDLE => $FILEHANDLE,
			       });

    ## Install rhocall
    print $FILEHANDLE "### Install rhocall\n";

    ## Create the temporary install directory
    create_install_dir({FILEHANDLE => $FILEHANDLE,
			install_directory => $parameter{rhocall_path},
		       });

    ## Download
    print $FILEHANDLE "## Download rhocall\n";
    print $FILEHANDLE "wget --quiet https://github.com/dnil/rhocall/archive/".$parameter_href->{rhocall}.".zip ";
    print $FILEHANDLE "-O rhocall-".$parameter_href->{rhocall}.".zip";  #Download outfile
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "unzip rhocall-".$parameter_href->{rhocall}.".zip";
    print $FILEHANDLE "\n\n";

    ## Move to rhocall directory
    print $FILEHANDLE "## Move to rhocall directory\n";
    print $FILEHANDLE "cd rhocall-".$parameter_href->{rhocall};
    print $FILEHANDLE "\n\n";

    ## Configure
    print $FILEHANDLE "## Configure\n";
    print $FILEHANDLE "pip install numpy Cython", "\n";
    print $FILEHANDLE "pip install -r requirements.txt", "\n";
    print $FILEHANDLE "pip install -e .";
    print $FILEHANDLE "\n\n";

    ## Moving up
    print $FILEHANDLE "## Moving back to original working directory\n";
    print $FILEHANDLE "cd ".$pwd;  #Go back to subroutine origin
    print $FILEHANDLE "\n\n";

    ## Deactivate conda environment
    deactivate_conda_environment({FILEHANDLE => $FILEHANDLE,
				 });
}


sub activate_conda_environment {

##activate_conda_environment

##Function : Activate conda environment
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE       => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Activate conda environment
    print $FILEHANDLE "## Activate conda environment\n";
    print $FILEHANDLE "source activate ".$parameter_href->{conda_environment}." ";
    print $FILEHANDLE "\n\n";
}


sub deactivate_conda_environment {

##deactivate_conda_environment

##Function : Deactivate conda environment
##Returns  : ""
##Arguments: $FILEHANDLE
##         : $FILEHANDLE => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;

    my $tmpl = {
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Deactivate conda environment
    print $FILEHANDLE "## Deactivate conda environment\n";
    print $FILEHANDLE "source deactivate ";
    print $FILEHANDLE "\n\n";
}


sub remove_install_dir {

##remove_install_dir

##Function : Remove the temporary install directory
##Returns  : ""
##Arguments: $FILEHANDLE, $pwd, $install_directory
##         : $FILEHANDLE        => FILEHANDLE to write to
##         : $pwd               => The original working directory
##         : $install_directory => Temporary installation directory

    my ($arg_href) = @_;

    ## Default(s)
    my $install_directory;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $pwd;

    my $tmpl = {
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
	pwd => { required => 1, defined => 1, strict_type => 1, store => \$pwd},
	install_directory => { default => ".MIP",
			       allow => qr/^\.\S+$/,
			       strict_type => 1, store => \$install_directory},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Moving up
    print $FILEHANDLE "## Moving back to original working directory\n";
    print $FILEHANDLE "cd ".$pwd;  #Go back to subroutine origin
    print $FILEHANDLE "\n\n";

    ## Clean up
    print $FILEHANDLE "## Clean up\n";
    print $FILEHANDLE "rm -rf ".$install_directory;
    print $FILEHANDLE "\n\n";
}

sub create_install_dir {

##create_install_dir

##Function : Create the temporary install directory
##Returns  : ""
##Arguments: $FILEHANDLE, $install_directory
##         : $FILEHANDLE        => FILEHANDLE to write to
##         : $install_directory => Temporary installation directory

    my ($arg_href) = @_;

    ## Default(s)
    my $install_directory;

    ## Flatten argument(s)
    my $FILEHANDLE;

    my $tmpl = {
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
	install_directory => { default => ".MIP",
			       strict_type => 1, store => \$install_directory},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Create temp install directory
    print $FILEHANDLE "## Create temp install directory\n";
    print $FILEHANDLE "mkdir -p ".$install_directory, "\n";
    print $FILEHANDLE "cd ".$install_directory;
    print $FILEHANDLE "\n\n";
}


sub check_conda_bin_file_exists {

##check_conda_bin_file_exists

##Function : Check if the binary of the program being installed already exists
##Returns  : ""
##Arguments: $parameter_href, $program_name, $program_version
##         : $parameter_href  => Holds all parameters
##         : $program_name    => Program name
##         : $program_version => Program version

    my ($arg_href) = @_;

    ## Default(s)
    my $program_version;

    ## Flatten argument(s)
    my $parameter_href;
    my $program_name;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
	program_name => { required => 1, defined => 1, strict_type => 1, store => \$program_name},
	program_version => { default => 0,
			     strict_type => 1, store => \$program_version},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $miniconda_bin_file;

    if ($program_version) {

	$miniconda_bin_file = catfile($parameter{conda_path}, "envs", $parameter_href->{conda_environment}, "bin", $program_name.$program_version);
    }
    else {

	$miniconda_bin_file = catfile($parameter{conda_path}, "envs", $parameter_href->{conda_environment}, "bin", $program_name);
    }

    if (-f $miniconda_bin_file) {

	if ($program_version) {

	    print STDERR q?Found ?.$program_name.q? version ?.$program_version.q? in miniconda directory: ?.catdir($parameter{conda_path}, "envs", $parameter_href->{conda_environment}, "bin"), "\n";

	    if ($parameter_href->{noupdate}) {

		print STDERR q?Skipping writting installation process for ?.$program_name.q? ?.$program_version, "\n";
		return 1;
	    }
	    print STDERR "Writting install instructions for ".$program_name, "\n";
	}
	else {

	    print STDERR q?Found ?.$program_name.q? in miniconda directory: ?.catdir($parameter{conda_path}, "envs", $parameter_href->{conda_environment}, "bin"), "\n";

	    if ($parameter_href->{noupdate}) {

		print STDERR q?Skipping writting installation process for ?.$program_name, "\n";
		return 1;
	    }
	    print STDERR "Writting install instructions for ".$program_name, "\n";
	}
	return 0;
    }
    else {

	print STDERR "Writting install instructions for ".$program_name, "\n";
	return 0;
    }
}


sub create_softlink {

##create_softlink

##Function : Create softlink
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE, $binary, $softlink
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => FILEHANDLE to write to
##         : $binary         => The binary file
##         : $softlink       => The name of the softlink

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;
    my $binary;
    my $softlink;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
	binary => { required => 1, defined => 1, strict_type => 1, store => \$binary},
	softlink => { required => 1, defined => 1, strict_type => 1, store => \$softlink},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Add softlink
    print $FILEHANDLE "## Move to directory and create softlink\n";
    print $FILEHANDLE "cd ".catdir($parameter{conda_path}, "envs", $parameter_href->{conda_environment}, "bin"), "\n";

    print $FILEHANDLE "ln -f -s ";
    print $FILEHANDLE $binary.q? ?.$softlink;
    print $FILEHANDLE "\n\n";

    ## Move to back
    print $FILEHANDLE "## Move to original working directory\n";
    print $FILEHANDLE "cd ".$pwd;
    print $FILEHANDLE "\n\n";
}


sub enable_executable {

##enable_executable

##Function : Make file executable
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE, $binary
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => FILEHANDLE to write to
##         : $binary         => The binary file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;
    my $binary;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
	binary => { required => 1, defined => 1, strict_type => 1, store => \$binary},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Enable executable
    print $FILEHANDLE "## Enable executable\n";
    print $FILEHANDLE "cd ".catdir($parameter{conda_path}, "envs", $parameter_href->{conda_environment}, "bin"), "\n";

    print $FILEHANDLE "chmod a+x ";
    print $FILEHANDLE $binary.q? ?;
    print $FILEHANDLE "\n\n";

    ## Move to back
    print $FILEHANDLE "## Move to original working directory\n";
    print $FILEHANDLE "cd ".$pwd;
    print $FILEHANDLE "\n\n";
}


sub check_mt_codon_table {

##check_mt_codon_table

##Function : Check and if required add the vertebrate mitochondrial codon table to snpeff config
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE, $share_dir, $config_file, $genome_version_ref
##         : $parameter_href     => Holds all parameters
##         : $FILEHANDLE         => FILEHANDLE to write to
##         : $share_dir          => The conda env shared directory
##         : $config_file        => The config config_file
##         : $genome_version_ref => snpeff genome version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;
    my $share_dir;
    my $config_file;
    my $genome_version_ref;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
	share_dir => { required => 1, defined => 1, strict_type => 1, store => \$share_dir},
	config_file => { required => 1, defined => 1, strict_type => 1, store => \$config_file},
	genome_version_ref => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$genome_version_ref},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();
    my $detect_regexp = q?perl -nae 'if($_=~/?.$$genome_version_ref.q?.MT.codonTable/) {print 1}' ?;
    my $add_regexp = q?perl -nae 'if($_=~/?.$$genome_version_ref.q?.reference/) {print $_; print "?.$$genome_version_ref.q?.MT.codonTable : Vertebrate_Mitochondrial\n"} else {print $_;}' ?;
    my $ret;

    if (-f catfile($share_dir, $config_file)) {

	$ret = `$detect_regexp $share_dir/$config_file`;
    }
    if (!$ret) {  #No MT.codonTable in config

	print $FILEHANDLE q?## Adding ?.$$genome_version_ref.q?.MT.codonTable : Vertebrate_Mitochondrial to ?.$share_dir.$config_file, "\n";

	## Add MT.codon Table to config
	print $FILEHANDLE $add_regexp." ".catfile($share_dir, $config_file)." > ".catfile($share_dir, $config_file.".tmp"), "\n";
	print $FILEHANDLE "mv ".catfile($share_dir, $config_file.".tmp")." ".catfile($share_dir, $config_file);
	print $FILEHANDLE "\n\n";

    }
    else {

	print STDERR  "Found MT.codonTable in ".catfile($share_dir, "snpEff.config").". Skipping addition to snpEff config", "\n";
    }
}


sub snpeff_download {

##snpeff_download

##Function : Write instructions to download snpeff database. This is done by install script to avoid race conditin when doing first analysis run in MIP
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE, $genome_version_ref
##         : $parameter_href     => Holds all parameters
##         : $FILEHANDLE         => FILEHANDLE to write to
##         : $genome_version_ref => snpeff genome version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;
    my $genome_version_ref;

    my $tmpl = {
	parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameter_href},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
	genome_version_ref => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$genome_version_ref},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Activate conda environment
    activate_conda_environment({parameter_href => $parameter_href,
				FILEHANDLE => $FILEHANDLE,
			       });

    print $FILEHANDLE "java -Xmx2g ";
    print $FILEHANDLE q?-jar ?.catfile($parameter_href->{conda_path}, "envs", $parameter_href->{conda_environment}, "bin", "snpEff.jar")." ";
    print $FILEHANDLE "download ";
    print $FILEHANDLE " -v ";
    print $FILEHANDLE $$genome_version_ref." ";
    print $FILEHANDLE q?-c ?.catfile($parameter_href->{conda_path}, "envs", $parameter_href->{conda_environment}, "bin", "snpEff.config")." ";
    print $FILEHANDLE "\n\n";

    ## Deactivate conda environment
    deactivate_conda_environment({FILEHANDLE => $FILEHANDLE,
				 });

}


sub help {

##help

##Function : Print help text and exit with supplied exit code
##Returns  : ""
##Arguments: $USAGE, $exit_code
##         : $USAGE     => Help text
##         : $exit_code => Exit code

    my ($arg_href) = @_;

    ## Default(s)
    my $exit_code;

    ## Flatten argument(s)
    my $USAGE;

    my $tmpl = {
	USAGE => { required => 1, defined => 1, strict_type => 1, store => \$USAGE},
	exit_code => { default => 0,
		       allow => qr/^\d+$/,
		       strict_type => 1, store => \$exit_code},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    print STDOUT $USAGE, "\n";
    exit $exit_code;
}
