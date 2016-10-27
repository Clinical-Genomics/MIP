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
           -env/--condaEnvironment Conda environment (Default: "mip")
           -cdp/--condaPath The conda path (Default: "HOME/miniconda")
           -cdu/--condaUpDate Update conda before installing (Supply flag to enable)
           -bvc/--bioConda Set the module version of the programs that can be installed with bioConda (e.g. 'bwa=0.7.12')
           -pip/--pip Set the module version of the programs that can be installed with pip (e.g. 'genmod=3.5.9')

           ## SHELL
           -pei/--perlInstall Install perl (Supply flag to enable)
           -per/--perl Set the perl version (defaults: "5.18.2")
           -pm/--perlModules Set the perl modules to be installed via cpanm (Default: ["Modern::Perl", "IPC::System::Simple", "Path::Iterator::Rule", "YAML", "Log::Log4perl", "Set::IntervalTree", "Net::SSLeay",P, "LWP::Simple", "LWP::Protocol::https", "Archive::Zip", "Archive::Extract", "DBI","JSON", "DBD::mysql", "CGI", "Sereal::Encoder", "Sereal::Decoder", "Bio::Root::Version", "Module::Build"])
           -pic/--picardTools Set the picardTools version (Default: "2.3.0"),
           -sbb/--sambamba Set the sambamba version (Default: "0.6.1")
           -vct/--vcfTools Set the vcftools version (Default: "0.1.14")
           -bet/--bedTools Set the bedtools version (Default: "2.25.0")
           -vt/--vt Set the vt version (Default: "0.57")
           -plk/--plink  Set the plink version (Default: "160224")
           -snpg/--snpEffGenomeVersions Set the snpEff genome version (Default: ["GRCh37.75"])
           -vep/--variantEffectPredictor Set the VEP version (Default: "85")
	   -vepc/--vepDirectoryCache Specify the cache directory to use (whole path; defaults to "~/miniconda/envs/condaEnvironment/ensembl-tools-release-variantEffectPredictorVersion/cache")
           -vepa/--vepAssemblies Select the assembly version (Default: ["GRCh37"])
           -vepp/--variantEffectPredictorPlugin Supply a comma separated list of VEP plugins (Default: "UpDownDistance,LoFtool,LoF")

           ## Utility
           -psh/--preferShell Shell will be used for overlapping shell and biconda installations (Supply flag to enable)
           -ppd/--printParameterDefaults Print the parameter defaults
           -nup/--noUpdate Do not update already installed programs (Supply flag to enable)
           -sp/--selectPrograms Install supplied programs e.g. -sp perl -sp bedTools (Default: "";)
           -h/--help Display this help message   
           -v/--version Display version
        };    
}
my $installDirectory = ".MIP";

my %parameter; 

### Set parameter default

##Conda
$parameter{condaEnvironment} = "mip";
$parameter{condaPath} = catdir($ENV{HOME}, "miniconda");

$parameter{bioConda}{bwa} = "0.7.15";
$parameter{bioConda}{bwakit} = "0.7.12";
$parameter{bioCondaBwakitPatch} = "-0";  #For correct softlinking in share and bin in conda env
$parameter{bioConda}{fastqc} = "0.11.5";
$parameter{bioConda}{cramtools} = "3.0.b47";
$parameter{bioConda}{samtools} = "1.3.1";
$parameter{bioConda}{bcftools} = "1.3.1";
$parameter{bioConda}{snpeff} = "4.2";
$parameter{bioCondaSnpeffPatch} = "-0";  #For correct softlinking in share and bin in conda env
$parameter{bioConda}{picard} = "2.5.0";
$parameter{bioCondaPicardPatch} = "-1";  #For correct softlinking in share and bin in conda env
$parameter{bioConda}{mosaik} = "2.2.26";
$parameter{bioConda}{htslib} = "1.3.1";
$parameter{bioConda}{bedtools} = "2.26.0";
$parameter{bioConda}{vt} = "2015.11.10";
$parameter{bioConda}{sambamba} = "0.6.3";
$parameter{bioConda}{freebayes} = "1.0.2.0";
$parameter{bioConda}{delly} = "0.7.2";
$parameter{bioConda}{manta} = "1.0.0";
$parameter{bioCondaMantaPatch} = "-0";
$parameter{bioConda}{multiqc} = "0.8dev0";
$parameter{bioConda}{plink2} = "1.90b3.35";
$parameter{bioConda}{vcfanno} = "0.1.0";
$parameter{bioConda}{gcc} = "4.8.5";  #Required for CNVnator
#$parameter{bioConda}{cmake} = "3.3.1";
#$parameter{bioConda}{boost} = "1.57.0";
#$parameter{bioCondaBoostPatch} = "-4";


##Perl Modules
$parameter{perl} = "5.18.2";

## PIP
$parameter{pip}{genmod} = "3.5.9";
$parameter{pip}{variant_integrity} = "0.0.4";
$parameter{pip}{chanjo} = "3.4.1";
$parameter{pip}{cosmid} = "0.4.9.1";
$parameter{pip}{'python-Levenshtein'} = "0.12.0";

## Programs currently installable by SHELL
$parameter{MIPScripts} = "Your current MIP version";
$parameter{picardTools} = "2.3.0";
$parameter{sambamba} = "0.6.1";
$parameter{vcfTools} = "0.1.14";
$parameter{bedTools} = "2.25.0";
$parameter{vt} = "gitRepo";
$parameter{plink2} = "160316";
$parameter{snpEff} = "v4_2";
$parameter{variantEffectPredictor} = "85";
$parameter{variantEffectPredictorPlugin} = "UpDownDistance,LoFtool,LoF";
#$parameter{CNVnator} = "0.3.2";
#$parameter{FindTranslocations} = "0";

my $installVersion = "0.0.6";

###User Options
GetOptions('env|condaEnvironment:s'  => \$parameter{condaEnvironment},
	   'cdp|condaPath:s' => \$parameter{condaPath},
	   'cdu|condaUpDate' => \$parameter{condaUpDate},
	   'bcv|bioConda=s' => \%{$parameter{bioConda}},
	   'pip|pip=s' => \%{$parameter{pip}},
	   'per|perl=s' => \$parameter{perl},
	   'pei|perlInstall' => \$parameter{perlInstall},
	   'pm|perlModules:s' => \@{$parameter{perlModules}},  #Comma separated list
	   'pic|picardTools:s' => \$parameter{picardTools},
	   'sbb|sambamba:s' => \$parameter{sambamba},
	   'vct|vcfTools:s' => \$parameter{vcfTools},
	   'bet|bedTools:s' =>\$parameter{bedTools}, 
	   'vt|vt:s' => \$parameter{vt},
	   'plk|plink2:s' => \$parameter{plink2},
	   'snpg|snpEffGenomeVersions:s' => \@{$parameter{snpEffGenomeVersions}},
	   'vep|variantEffectPredictor:s' => \$parameter{variantEffectPredictor},
	   'vepc|vepDirectoryCache:s' => \$parameter{vepDirectoryCache},  #path to vep cache dir
	   'vepa|vepAssemblies:s' => \@{$parameter{vepAssemblies}},  #Select assembly version to use
	   'vepp|variantEffectPredictorPlugin:s' => \$parameter{variantEffectPredictorPlugin},  #Comma sep string
#	   'cnv|CNVnator:s' => \$parameter{CNVnator},
#	   'ftr|FindTranslocations:s' => \$parameter{FindTranslocations},
	   'psh|preferShell' => \$parameter{preferShell},  # Shell will be used for overlapping shell and biconda installations
	   'ppd|printParameterDefaults' => sub { &PrintParameters({parameterHashRef => \%parameter}); exit;},  #Display parameter defaults
	   'nup|noUpdate' => \$parameter{noUpdate},
	   'sp|selectPrograms:s' => \@{$parameter{selectPrograms}},  #Comma sep string
	   'h|help' => sub { print STDOUT $USAGE, "\n"; exit;},  #Display help text
	   'v|version' => sub { print STDOUT "\ninstall.pl ".$installVersion, "\n\n"; exit;},  #Display version number
    ) or &Help({USAGE => $USAGE,
		exitCode => 1,
	       });

## Update default parameter dependent on other parameters
if (! $parameter{vepDirectoryCache}) {

    $parameter{vepDirectoryCache} = catdir($parameter{condaPath}, "envs", $parameter{condaEnvironment}, "ensembl-tools-release-".$parameter{variantEffectPredictor}, "cache");  #Cache directory;)
}

## Set default for array parameters
&SetDefaultArrayParameters({parameterHashRef => \%parameter,
			   });

###MAIN###

#my $LOGFILEHANDLE = &OpenLogFile({fileName => "MIP_installation.log",
#				 });

## Create bash file for writing install instructions
my $BASHFILEHANDLE = &CreateBashFile({fileName => "mip.sh",
				     });

## Check existance of conda environment
&CheckConda({parameterHashRef => \%parameter,
	     FILEHANDLE => $BASHFILEHANDLE,
	    });

## Create Conda environment if required
&CreateCondaEnvironment({parameterHashRef => \%parameter,
			 FILEHANDLE => $BASHFILEHANDLE,
			});

## Install modules into Conda environment using channel Bioconda
&InstallBiocondaModules({parameterHashRef => \%parameter,
			 FILEHANDLE => $BASHFILEHANDLE,
			});

if (@{$parameter{selectPrograms}}) {

    if ( ( grep {$_ eq "perl"} @{$parameter{selectPrograms}} ) ) { #If element is part of array

	&Perl({parameterHashRef => \%parameter,
	       FILEHANDLE => $BASHFILEHANDLE,
	      });
    }
}
else {
    
    &Perl({parameterHashRef => \%parameter,
	   FILEHANDLE => $BASHFILEHANDLE,
	  });
}


&PipInstall({parameterHashRef => \%parameter,
	     FILEHANDLE => $BASHFILEHANDLE,
	    });

if ($parameter{preferShell}) {

    if (@{$parameter{selectPrograms}}) {
	
	if ( ( grep {$_ eq "picardTools"} @{$parameter{selectPrograms}} ) ) { #If element is part of array
	    
	    &PicardTools({parameterHashRef => \%parameter,
			  FILEHANDLE => $BASHFILEHANDLE,
			 });
	}
	if ( ( grep {$_ eq "sambamba"} @{$parameter{selectPrograms}} ) ) { #If element is part of array

	    &Sambamba({parameterHashRef => \%parameter,
		       FILEHANDLE => $BASHFILEHANDLE,
		      });
	}
	if ( ( grep {$_ eq "bedTools"} @{$parameter{selectPrograms}} ) ) { #If element is part of array

	    &BedTools({parameterHashRef => \%parameter,
		       FILEHANDLE => $BASHFILEHANDLE,
		      });
	}
	if ( ( grep {$_ eq "vt"} @{$parameter{selectPrograms}} ) ) { #If element is part of array

	    &VT({parameterHashRef => \%parameter,
		 FILEHANDLE => $BASHFILEHANDLE,
		});
	}
	if ( ( grep {$_ eq "snpEff"} @{$parameter{selectPrograms}} ) ) { #If element is part of array

	    &SnpEff({parameterHashRef => \%parameter,
		     FILEHANDLE => $BASHFILEHANDLE,
		    });
	}
	if ( ( grep {$_ eq "plink2"} @{$parameter{selectPrograms}} ) ) { #If element is part of array
	
	    &Plink2({parameterHashRef => \%parameter,
		     FILEHANDLE => $BASHFILEHANDLE,
		    });
	}
    }
    else {
	
	&PicardTools({parameterHashRef => \%parameter,
		      FILEHANDLE => $BASHFILEHANDLE,
		     });
	
	&Sambamba({parameterHashRef => \%parameter,
		   FILEHANDLE => $BASHFILEHANDLE,
		  });
	
	&BedTools({parameterHashRef => \%parameter,
		   FILEHANDLE => $BASHFILEHANDLE,
		  });
	
	&VT({parameterHashRef => \%parameter,
	     FILEHANDLE => $BASHFILEHANDLE,
	    });
	
	&SnpEff({parameterHashRef => \%parameter,
		 FILEHANDLE => $BASHFILEHANDLE,
		});

	&Plink2({parameterHashRef => \%parameter,
		 FILEHANDLE => $BASHFILEHANDLE,
		});
    }
}

if (@{$parameter{selectPrograms}}) {
    
    if ( ( grep {$_ eq "vcfTools"} @{$parameter{selectPrograms}} ) ) { #If element is part of array
	
	&VcfTools({parameterHashRef => \%parameter,
		   FILEHANDLE => $BASHFILEHANDLE,
		  });
    }
    if ( ( grep {$_ eq "MIPScripts"} @{$parameter{selectPrograms}} ) ) { #If element is part of array

	&MIPScripts({parameterHashRef => \%parameter,
		     FILEHANDLE => $BASHFILEHANDLE,
		    });
    }
    if ( ( grep {$_ eq "variantEffectPredictor"} @{$parameter{selectPrograms}} ) ) { #If element is part of array
	
	&VariantEffectPredictor({parameterHashRef => \%parameter,
				 FILEHANDLE => $BASHFILEHANDLE,
				});
    }
#    if ( ( grep {$_ eq "CNVnator"} @{$parameter{selectPrograms}} ) ) { #If element is part of array
	
#	&CNVnator({parameterHashRef => \%parameter,
#		   FILEHANDLE => $BASHFILEHANDLE,
#		  });
#   }
#    if ( ( grep {$_ eq "FindTranslocations"} @{$parameter{selectPrograms}} ) ) { #If element is part of array
	
#	&FindTranslocations({parameterHashRef => \%parameter,
#			     FILEHANDLE => $BASHFILEHANDLE,
#			    });
#    }
}
else {
    
    &VcfTools({parameterHashRef => \%parameter,
	       FILEHANDLE => $BASHFILEHANDLE,
	      });

    &MIPScripts({parameterHashRef => \%parameter,
		 FILEHANDLE => $BASHFILEHANDLE,
		});
    
    &VariantEffectPredictor({parameterHashRef => \%parameter,
			     FILEHANDLE => $BASHFILEHANDLE,
			    });

#    &CNVnator({parameterHashRef => \%parameter,
#	       FILEHANDLE => $BASHFILEHANDLE,
#	      });

#    &FindTranslocations({parameterHashRef => \%parameter,
#			 FILEHANDLE => $BASHFILEHANDLE,
#			});
}

close($BASHFILEHANDLE);
#close($LOGFILEHANDLE);

###SubRoutines###

sub SetDefaultArrayParameters {
    
##SetDefaultArrayParameters
    
##Function : Set default for array parameters
##Returns  : ""
##Arguments: $parameterHashRef
##         : $parameterHashRef => Holds all parameters
    
    my ($argHashRef) = @_;
    
    ## Flatten argument(s)
    my $parameterHashRef;

    my $tmpl = { 
	parameterHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameterHashRef},
    };
        
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my %arrayParameter;
    $arrayParameter{vepAssemblies}{default} = ["GRCh37"];
    $arrayParameter{snpEffGenomeVersions}{default} = ["GRCh37.75"];  #GRCh38.82 but check current on the snpEff sourceForge
    $arrayParameter{perlModules}{default} = ["Modern::Perl",  #MIP
					     "IPC::System::Simple",  #MIP
					     "Path::Iterator::Rule",  #MIP
					     "YAML",  #MIP
					     "Log::Log4perl",  #MIP
					     "Set::IntervalTree",  # MIP/vcfParser.pl
					     "Net::SSLeay",  # VEP
					     "LWP::Simple",  # VEP
					     "LWP::Protocol::https",  # VEP
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
	];
    
    foreach my $parameterName (keys %arrayParameter) {
	
	if (! @{${$parameterHashRef}{$parameterName}}) {  #Unless parameter was supplied on cmd
	    
	    ${$parameterHashRef}{$parameterName} = $arrayParameter{$parameterName}{default};
	}
    }
}

sub CreateBashFile {

##CreateBashFile
    
##Function : Create bash file for writing install instructions
##Returns  : ""
##Arguments: $fileName
##         : $fileName => File name 

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $fileName;

    my $tmpl = { 
	fileName => { required => 1, defined => 1, strict_type => 1, store => \$fileName},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my $FILEHANDLE = IO::Handle->new();  #Create anonymous filehandle
    my $pwd = cwd();

    ## Open batch file
    open ($FILEHANDLE, ">", catfile($pwd, $fileName)) or die("Cannot write to '".catfile($pwd, $fileName)."' :".$!."\n");

    print $FILEHANDLE "#!".catfile( dirname( dirname( devnull() ) ) ).catfile("usr", "bin", "env", "bash"), "\n\n";
 
    print STDOUT "Will write install instructions to '".catfile($pwd, $fileName), "'\n";

   return $FILEHANDLE;
}

sub OpenLogFile {

##OpenLogFile
    
##Function : Open log file
##Returns  : ""
##Arguments: $fileName
##         : $fileName => File name 

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $fileName;

    my $tmpl = { 
	fileName => { required => 1, defined => 1, strict_type => 1, store => \$fileName},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my $FILEHANDLE = IO::Handle->new();  #Create anonymous filehandle
    
    ## Open batch file
    open ($FILEHANDLE, ">", catfile($fileName)) or die("Cannot write to '".catfile($fileName)."' :".$!."\n");

    return $FILEHANDLE;
}

sub PrintParameters {
    
##PrintParameters
    
##Function : Print all parameters and the default values
##Returns  : ""
##Arguments: $parameterHashRef
##         : $parameterHashRef => Holds all parameters

    my ($argHashRef) = @_;
    
    ## Flatten argument(s)
    my $parameterHashRef;
    
    my $tmpl = { 
	parameterHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameterHashRef},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];
    
    ## Set default for array parameters
    &SetDefaultArrayParameters({parameterHashRef => $parameterHashRef,
			       });
    
    foreach my $key (keys %{$parameterHashRef}) {
	
	if (ref(${$parameterHashRef}{$key})!~/ARRAY|HASH/) {

	    print STDOUT $key." ";
	    if (${$parameterHashRef}{$key}) {

		print ${$parameterHashRef}{$key}, "\n";
	    }
	    else {  ##Boolean value

		print "0", "\n";
	    }
	}
	elsif (ref(${$parameterHashRef}{$key})=~/HASH/) {
	    
	    foreach my $program (keys %{${$parameterHashRef}{$key}}) {
		
		print STDOUT $key." ".$program.": ".${$parameterHashRef}{$key}{$program}, "\n";
	    }
	}
	elsif (ref(${$parameterHashRef}{$key})=~/ARRAY/)  {

	    print STDOUT $key.": ".join(" ", @{${$parameterHashRef}{$key}}), "\n";
	}
    }
}


sub CheckConda {

##CheckConda
    
##Function : Check existance of conda environment
##Returns  : ""
##Arguments: $parameterHashRef
##         : $parameterHashRef => Holds all parameters
##         : $FILEHANDLE       => Filehandle to write to

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $parameterHashRef;
    my $FILEHANDLE;

    my $tmpl = { 
	parameterHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameterHashRef},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my $program = "conda";

    if ($ENV{PATH}=~/conda/) {

	print STDERR "ProgramCheck: ".$program." installed\n";
    }
    else {
	
	print STDERR "Could not detect ".$program." in your PATH\n";
	exit 1;
    }

    ##Deactivate any activate env prior to installation
    my $detectActiveCondaEnv = q?perl -nae 'if( ($_!~/^root/) && ($_=~/\*/) ) {print $F[0]}'?;
    my $ret = `conda info --envs | $detectActiveCondaEnv`;

    if ($ret) {

	print STDOUT "Found activated conda env: ".$ret."\n";
	print STDOUT "Please exit conda env: ".$ret." before executing install script\n";
	exit 1;
    }

    ## Check Conda path
    if (! -d ${$parameterHashRef}{condaPath}) {

	print STDERR "Could not find miniconda directory in: ".catdir(${$parameterHashRef}{condaPath}), "\n";
	exit 1;
    }

    print STDERR "Writting install instructions for Conda packages\n";

    ## Update Conda
    if (${$parameterHashRef}{condaUpDate}) {

	print $FILEHANDLE "### Update Conda\n";
	print $FILEHANDLE "conda update -y conda ";
	print $FILEHANDLE "\n\n";
    }
}

sub CreateCondaEnvironment {

##CreateCondaEnvironment
    
##Function : Create Conda environment if required
##Returns  : ""
##Arguments: $parameterHashRef
##         : $parameterHashRef => Holds all parameters
##         : $FILEHANDLE       => Filehandle to write to

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $parameterHashRef;
    my $FILEHANDLE;

    my $tmpl = { 
	parameterHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameterHashRef},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    ## Check Conda environment
    if (! -d catdir($parameter{condaPath}, "envs", $parameter{condaEnvironment})) {

	## Create conda environment
	print $FILEHANDLE "### Creating Conda Environment and install: ".${$parameterHashRef}{condaEnvironment}, "\n";
	print $FILEHANDLE "conda create -n ".${$parameterHashRef}{condaEnvironment}." ";
	print $FILEHANDLE "-y ";
	print $FILEHANDLE "pip ";
	print $FILEHANDLE "\n\n";
    }
}

sub InstallBiocondaModules {
    
##InstallBiocondaModules
    
##Function : Install modules into Conda environment using channel Bioconda
##Returns  : ""
##Arguments: $parameterHashRef
##         : $parameterHashRef => Holds all parameters
##         : $FILEHANDLE       => Filehandle to write to

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $parameterHashRef;
    my $FILEHANDLE;

    my $tmpl = { 
	parameterHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameterHashRef},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    ## Install into conda environment using bioconda channel
    print $FILEHANDLE "### Installing into Conda Environment: ".${$parameterHashRef}{condaEnvironment}, "\n";
    print $FILEHANDLE "conda install ";	
    print $FILEHANDLE "-n ".${$parameterHashRef}{condaEnvironment}." ";
    print $FILEHANDLE "-y ";
    print $FILEHANDLE "-c bioconda ";
    
    ## Install all bioConda packages
    foreach my $program (keys %{${$parameterHashRef}{bioConda}}) {
	    
	print $FILEHANDLE $program."=".${$parameterHashRef}{bioConda}{$program}." ";
    }

    print $FILEHANDLE "\n\n";

    ## Custom 
    foreach my $program (keys %{${$parameterHashRef}{bioConda}}) {

	if ($program eq "bwakit") {

	    ## Define binaries
	    my @bwaKitBinaries = ("k8",
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

	    foreach my $binary (@bwaKitBinaries) {
		
		&CreateSoftLink({parameterHashRef => $parameterHashRef,
				 FILEHANDLE => $BASHFILEHANDLE,
				 binary => catfile($parameter{condaPath}, "envs", $parameter{condaEnvironment}, "share", "bwakit-".${$parameterHashRef}{bioConda}{bwakit}.${$parameterHashRef}{bioCondaBwakitPatch}, $binary),
				 softLink => $binary,
				});
	    }

	    print $BASHFILEHANDLE "cp -rf ".catdir($parameter{condaPath}, "envs", $parameter{condaEnvironment}, "share", "bwakit-".${$parameterHashRef}{bioConda}{bwakit}.${$parameterHashRef}{bioCondaBwakitPatch}, "resource-human-HLA")." ";
	    print $BASHFILEHANDLE catdir($parameter{condaPath}, "envs", $parameter{condaEnvironment}, "bin"), "\n\n";
	}
	if ($program eq "picard") {

	    &CreateSoftLink({parameterHashRef => $parameterHashRef,
			     FILEHANDLE => $BASHFILEHANDLE,
			     binary => catfile($parameter{condaPath}, "envs", $parameter{condaEnvironment}, "share", "picard-".${$parameterHashRef}{bioConda}{picard}.${$parameterHashRef}{bioCondaPicardPatch}, "picard.jar"),
			     softLink => "picard.jar",
			    });
	}
	if ($program eq "snpeff") {
	    
	    ## Define binaries
	    my @snpEffBinaries = ("snpEff.jar",
				  "SnpSift.jar",
				  "snpEff.config",
		);

	    foreach my $binary (@snpEffBinaries) {
		
		&CreateSoftLink({parameterHashRef => $parameterHashRef,
				 FILEHANDLE => $BASHFILEHANDLE,
				 binary => catfile($parameter{condaPath}, "envs", $parameter{condaEnvironment}, "share", "snpeff-".${$parameterHashRef}{bioConda}{snpeff}.${$parameterHashRef}{bioCondaSnpeffPatch}, $binary),
				 softLink => $binary,
				});
	    }

	    foreach my $genomeVersion (@{${$parameterHashRef}{snpEffGenomeVersions}}) {

		## Check and if required add the vertebrate mitochondrial codon table to SnpEff config
		&CheckMTCodonTable({parameterHashRef => $parameterHashRef,
				    FILEHANDLE => $BASHFILEHANDLE,
				    shareDirectory => catdir($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment}, "share", "snpeff-".${$parameterHashRef}{bioConda}{snpeff}.${$parameterHashRef}{bioCondaSnpeffPatch}),
				    configFile => "snpEff.config",
				    genomeVersionRef => \$genomeVersion,
				   });

		unless (-d catdir($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment}, "share", "snpeff-".${$parameterHashRef}{bioConda}{snpeff}.${$parameterHashRef}{bioCondaSnpeffPatch}, "data", $genomeVersion)) {
		    
		    ## Write instructions to download SnpEff database. This is done by install script to avoid race conditin when doing first analysis run in MIP
		    &SnpEffDownload({parameterHashRef => $parameterHashRef,
				     FILEHANDLE => $BASHFILEHANDLE,
				     genomeVersionRef => \$genomeVersion,
				    });
		}	
	    }
	}
	if ($program eq "manta") {
	    
	    my @mantaBinaries = ("configManta.py",
				 "configManta.py.ini",
		);

	    foreach my $binary (@mantaBinaries) {

		&CreateSoftLink({parameterHashRef => $parameterHashRef,
				 FILEHANDLE => $BASHFILEHANDLE,
				 binary => catfile($parameter{condaPath}, "envs", $parameter{condaEnvironment}, "share", "manta-".${$parameterHashRef}{bioConda}{manta}.${$parameterHashRef}{bioCondaMantaPatch}, "bin", $binary),
				 softLink => $binary,
				});
	    }
	    
	    ## Make file executable
	    &EnableExecutable({parameterHashRef => $parameterHashRef,
			       FILEHANDLE => $BASHFILEHANDLE,
			       binary => q?configManta.py?,
			      });
	}
    }
}


sub Perl {

##Perl
    
##Function : Installs perl
##Returns  : ""
##Arguments: $parameterHashRef, $FILEHANDLE
##         : $parameterHashRef => Holds all parameters
##         : $FILEHANDLE       => Filehandle to write to

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $parameterHashRef;
    my $FILEHANDLE;

    my $tmpl = { 
	parameterHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameterHashRef},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];
    
    my $pwd = cwd();
    
    if ($ENV{PATH}=~/perl-${$parameterHashRef}{perl}/) {
	
	if (${$parameterHashRef}{noUpdate}) {
	    
	    print STDERR "Found perl-".${$parameterHashRef}{perl}.". in your path", "\n";
	    print STDERR q?Skipping writting installation for perl-?.${$parameterHashRef}{perl}, "\n";  
	}
	else {

	    if (${$parameterHashRef}{perlInstall}) {
		
		## Removing specific Perl version
		print $FILEHANDLE "### Removing specific Perl version\n";
		print $FILEHANDLE q?rm -rf $HOME/perl-?.${$parameterHashRef}{perl};
		print $FILEHANDLE "\n\n";
		
		&InstallPerlCpnam({parameterHashRef => $parameterHashRef,
				   FILEHANDLE => $BASHFILEHANDLE,
				  }); 
	    }
		
	    &PerlModules({parameterHashRef => $parameterHashRef,
			  FILEHANDLE => $BASHFILEHANDLE,
			 });
	}
    }
    else {
    	
	if (${$parameterHashRef}{perlInstall}) {
	    
	    &InstallPerlCpnam({parameterHashRef => $parameterHashRef,
			       FILEHANDLE => $BASHFILEHANDLE,
			       path => 1,
			      });
	}

	&PerlModules({parameterHashRef => $parameterHashRef,
		      FILEHANDLE => $BASHFILEHANDLE,
		     });
    }
}


sub InstallPerlCpnam {

##InstallPerlCpnam
    
##Function : Install perl CPANM 
##Returns  : ""
##Arguments: $parameterHashRef, $FILEHANDLE
##         : $parameterHashRef => Holds all parameters
##         : $FILEHANDLE       => Filehandle to write to
##         : $path             => Export path if provided {Optional}

    my ($argHashRef) = @_;

    ## Default(s)
    my $path;

    ## Flatten argument(s)
    my $parameterHashRef;
    my $FILEHANDLE;

    my $tmpl = { 
	parameterHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameterHashRef},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
	path => { default => 0,
		  allow => [0, 1],
		  strict_type => 1, store => \$path},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();

    print STDERR "Writting install instructions for Perl and Cpanm\n";
    
    ## Install specific Perl version
    print $FILEHANDLE "### Install specific Perl version\n";
    
    ## Move to Home
    print $FILEHANDLE "## Move HOME\n";
    print $FILEHANDLE q?cd $HOME?;
    print $FILEHANDLE "\n\n";
    
    ## Download
    print $FILEHANDLE "## Download Perl\n";
    print $FILEHANDLE "wget --quiet http://www.cpan.org/src/5.0/perl-".${$parameterHashRef}{perl}.".tar.gz ";
    print $FILEHANDLE "-O perl-".${$parameterHashRef}{perl}.".tar.gz";  #Dowload outfile
    print $FILEHANDLE "\n\n";
    
    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "tar xzf perl-".${$parameterHashRef}{perl}.".tar.gz";
    print $FILEHANDLE "\n\n";
    
    ## Move to perl directory
    print $FILEHANDLE "## Move to perl directory\n";
    print $FILEHANDLE "cd perl-".${$parameterHashRef}{perl};
    print $FILEHANDLE "\n\n";
    
    ## Configure
    print $FILEHANDLE "## Configure\n";
    print $FILEHANDLE q?./Configure -des -Dprefix=$HOME/perl-?.${$parameterHashRef}{perl}, "\n";    
    print $FILEHANDLE "make", "\n";
    print $FILEHANDLE "make test", "\n";    
    print $FILEHANDLE "make install", "\n\n";
    
    if ($path) {
	
	## Export path
	print $FILEHANDLE "## Export path\n";
	print $FILEHANDLE q?echo 'export PATH=$HOME/perl-?.${$parameterHashRef}{perl}.q?/:$PATH' >> ~/.bashrc?;
	print $FILEHANDLE "\n\n";
	print $FILEHANDLE q?export PATH=$HOME/perl-?.${$parameterHashRef}{perl}.q?/:$PATH?;  #Use newly installed perl
	print $FILEHANDLE "\n\n";
    }

    ## Remove tar file
    print $FILEHANDLE "## Remove tar file\n";
    print $FILEHANDLE "cd && rm perl-".${$parameterHashRef}{perl}.".tar.gz";
    print $FILEHANDLE "\n\n";
    
    ## Move to back
    print $FILEHANDLE "## Move to original working directory\n";
    print $FILEHANDLE "cd ".$pwd;
    print $FILEHANDLE "\n\n";

    print $FILEHANDLE q?echo 'eval `perl -I ~/perl-?.${$parameterHashRef}{perl}.q?/lib/perl5/ -Mlocal::lib=~/perl-?.${$parameterHashRef}{perl}.q?/`' >> ~/.bash_profile ?;  #Add at start-up
    print $FILEHANDLE "\n\n";
    print $FILEHANDLE q?echo 'export PERL_UNICODE=SAD' >> ~/.bash_profile ?;  #Add at start-up
    print $FILEHANDLE "\n\n";

    ## Install Perl modules via cpanm
    print $FILEHANDLE "## Install cpanm\n";
    print $FILEHANDLE q?wget -O- http://cpanmin.us | perl - -l $HOME/perl-?.${$parameterHashRef}{perl}.q?/bin App::cpanminus --local-lib=~/perl-?.${$parameterHashRef}{perl}.q?/ local::lib ?;
    print $FILEHANDLE "\n\n";

    ## Use newly installed perl
    print $FILEHANDLE q?eval `perl -I ~/perl-?.${$parameterHashRef}{perl}.q?/lib/perl5/ -Mlocal::lib=~/perl-?.${$parameterHashRef}{perl}.q?/` ?;
    print $FILEHANDLE "\n\n";

    ## Use newly installed perl
    print $FILEHANDLE q?PERL5LIB=~/perl-?.${$parameterHashRef}{perl}.q?/lib/perl5?;
    print $FILEHANDLE "\n\n";
}
    

sub PerlModules {

##PerlModules
    
##Function : Install Perl modules via cpanm
##Returns  : ""
##Arguments: $parameterHashRef, $FILEHANDLE
##         : $parameterHashRef => Holds all parameters
##         : $FILEHANDLE       => Filehandle to write to

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $parameterHashRef;
    my $FILEHANDLE;

    my $tmpl = { 
	parameterHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameterHashRef},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    ## Install Perl modules via cpanm
    print $FILEHANDLE "## Install Perl modules via cpanm\n";
    print $FILEHANDLE "cpanm ";
    print $FILEHANDLE join(" ", @{${$parameterHashRef}{perlModules}})." ";
    print $FILEHANDLE "\n\n";
}


sub PipInstall {

##PipInstall
    
##Function : Writes install instructions for pip packages
##Returns  : ""
##Arguments: $parameterHashRef, $FILEHANDLE
##         : $parameterHashRef => Holds all parameters
##         : $FILEHANDLE       => Filehandle to write to

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $parameterHashRef;
    my $FILEHANDLE;

    my $tmpl = { 
	parameterHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameterHashRef},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    print STDERR "Writting install instructions for pip packages\n";

    ## Install PIP packages in conda environment
    print $FILEHANDLE "### Install PIP packages in conda environment: ".${$parameterHashRef}{condaEnvironment}, "\n";

    ## Activate conda environment
    &ActivateCondaEnvironment({parameterHashRef => $parameterHashRef,
			       FILEHANDLE => $FILEHANDLE,
			      });

    ## Install PIP packages
    print $FILEHANDLE "## Install PIP packages\n";
    print $FILEHANDLE "pip install ";

    ## Install all PIP packages
    foreach my $program (keys %{${$parameterHashRef}{pip}}) {

	print $FILEHANDLE $program."==".${$parameterHashRef}{pip}{$program}." ";
    }
    print $FILEHANDLE "\n\n";

    ## Deactivate conda environment
    &DeactivateCondaEnvironment({FILEHANDLE => $FILEHANDLE,
				});
}


sub PicardTools {

##PicardTools
    
##Function : Install PicardTools
##Returns  : ""
##Arguments: $parameterHashRef, $FILEHANDLE
##         : $parameterHashRef => Holds all parameters
##         : $FILEHANDLE       => Filehandle to write to

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $parameterHashRef;
    my $FILEHANDLE;

    my $tmpl = { 
	parameterHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameterHashRef},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Check if the binary of the program being installed already exists
    if (&CheckCondaBinFileExists({parameterHashRef => $parameterHashRef,
				  programName => "picard.jar",
				 })) {  # Assumes that picard.jar is there as well then
	
	return
    }

    ## Install picard
    print $FILEHANDLE "### Install Picard\n";

    ## Create the temporary install directory
    &CreateInstallDirectory({FILEHANDLE => $FILEHANDLE,
			    });
    
    ## Download
    print $FILEHANDLE "## Download Picard\n";
    print $FILEHANDLE "wget --quiet https://github.com/broadinstitute/picard/releases/download/".${$parameterHashRef}{picardTools}."/picard-tools-".${$parameterHashRef}{picardTools}.".zip ";
    print $FILEHANDLE "-O picard-tools-".${$parameterHashRef}{picardTools}.".zip";  #Download outfile
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "unzip picard-tools-".${$parameterHashRef}{picardTools}.".zip";
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    if (-d catdir($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment}, "share", "picard-tools-".${$parameterHashRef}{picardTools})) {

	print $FILEHANDLE "rm -rf ".catdir($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment}, "share", "picard-tools-".${$parameterHashRef}{picardTools});
	print $FILEHANDLE "\n\n";
    }

    print $FILEHANDLE "## Make available from conda environment\n";
    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?picard-tools-?.${$parameterHashRef}{picardTools}.q? ?.catdir($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment}, "share");
    print $FILEHANDLE "\n\n";

    &CreateSoftLink({parameterHashRef => $parameterHashRef,
		     FILEHANDLE => $FILEHANDLE,
		     binary => catfile($parameter{condaPath}, "envs", $parameter{condaEnvironment}, "share", "picard-tools-".${$parameterHashRef}{picardTools}, "picard.jar"),
		     softLink => "picard.jar",
		    });
    
    ## Remove the temporary install directory
    &RemoveInstallDirectory({FILEHANDLE => $FILEHANDLE,
			     pwd => $pwd,
			    });
}


sub Sambamba {

##Sambamba
    
##Function : Install sambamba
##Returns  : ""
##Arguments: $parameterHashRef, $FILEHANDLE
##         : $parameterHashRef => Holds all parameters
##         : $FILEHANDLE       => Filehandle to write to

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $parameterHashRef;
    my $FILEHANDLE;

    my $tmpl = { 
	parameterHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameterHashRef},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Check if the binary of the program being installed already exists
    if (&CheckCondaBinFileExists({parameterHashRef => $parameterHashRef,
				  programName => "sambamba",
				  programVersion => ${$parameterHashRef}{sambamba},
				 }) ) {
	return 
    }

    ## Install sambamba
    print $FILEHANDLE "### Install sambamba\n";

    ## Create the temporary install directory
    &CreateInstallDirectory({FILEHANDLE => $FILEHANDLE,
			    });

    ## Download
    print $FILEHANDLE "## Download sambamba release\n";
    print $FILEHANDLE q?wget --quiet https://github.com/lomereiter/sambamba/releases/download/v?.${$parameterHashRef}{sambamba}.q?/sambamba_v?.${$parameterHashRef}{sambamba}.q?_linux.tar.bz2 ?;
    print $FILEHANDLE "-O sambamba_v".${$parameterHashRef}{sambamba}."_linux.tar.bz2";  #Download outfile
    print $FILEHANDLE "\n\n";

    ## Decompress
    print $FILEHANDLE "## Decompress sambamba file\n";
    print $FILEHANDLE "bzip2 ";
    print $FILEHANDLE "-f ";  #Force
    print $FILEHANDLE "-d ";  #Decompress
    print $FILEHANDLE "sambamba_v".${$parameterHashRef}{sambamba}."_linux.tar.bz2";
    print $FILEHANDLE "\n\n";

    ## Extract files
    print $FILEHANDLE "## Extract files\n";
    print $FILEHANDLE "tar xvf sambamba_v".${$parameterHashRef}{sambamba}."_linux.tar";
    print $FILEHANDLE "\n\n";

    ## Make executable
    print $FILEHANDLE "## Make executable\n";
    print $FILEHANDLE "chmod 755 ";
    print $FILEHANDLE "sambamba_v".${$parameterHashRef}{sambamba};
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    print $FILEHANDLE "## Make available from conda environment\n";
    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?sambamba_v?.${$parameterHashRef}{sambamba}.q? ?.$parameter{condaPath}.q?/envs/?.${$parameterHashRef}{condaEnvironment}.q?/bin/?;
    print $FILEHANDLE "\n\n";

    &CreateSoftLink({parameterHashRef => $parameterHashRef,
		     FILEHANDLE => $BASHFILEHANDLE,
		     binary => "sambamba_v".${$parameterHashRef}{bioConda}{sambamba},
		     softLink => "sambamba",
		    });
    
    ## Remove the temporary install directory
    &RemoveInstallDirectory({FILEHANDLE => $FILEHANDLE,
			     pwd => $pwd,
			    });

}


sub VcfTools {

##VcfTools
    
##Function : Install VcfTools
##Returns  : ""
##Arguments: $parameterHashRef, $FILEHANDLE
##         : $parameterHashRef => Holds all parameters
##         : $FILEHANDLE       => Filehandle to write to

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $parameterHashRef;
    my $FILEHANDLE;

    my $tmpl = { 
	parameterHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameterHashRef},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Check if the binary of the program being installed already exists
    if(&CheckCondaBinFileExists({parameterHashRef => $parameterHashRef,
				 programName => "vcftools",
				})) {
    
	return
    }

    ## Install vcfTools
    print $FILEHANDLE "### Install vcfTools\n";

    ## Create the temporary install directory
    &CreateInstallDirectory({FILEHANDLE => $FILEHANDLE,
			    });
    
    ## Download
    print $FILEHANDLE "## Download vcfTools\n";
    print $FILEHANDLE "wget --quiet https://github.com/vcftools/vcftools/releases/download/v".${$parameterHashRef}{vcfTools}."/vcftools-".${$parameterHashRef}{vcfTools}.".tar.gz ";
    print $FILEHANDLE "-O vcftools-".${$parameterHashRef}{vcfTools}.".tar.gz";  #Download outfile
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "tar xvf vcftools-".${$parameterHashRef}{vcfTools}.".tar.gz";
    print $FILEHANDLE "\n\n";

    ## Export PERL5LIB environment variable
    print $FILEHANDLE "## Export PERL5LIB environment variable\n";
    print $FILEHANDLE q?export PERL5LIB=?.$Bin.q?/vcftools-?.${$parameterHashRef}{vcfTools}.q?/src/perl/?;
    print $FILEHANDLE "\n\n";

    ## Move to vcfTools directory
    print $FILEHANDLE "## Move to vcfTools directory\n";
    print $FILEHANDLE "cd vcftools-".${$parameterHashRef}{vcfTools};
    print $FILEHANDLE "\n\n";

    ## Configure
    my $filePath = $parameter{condaPath}.q?/envs/?.${$parameterHashRef}{condaEnvironment};

    print $FILEHANDLE "## Configure", "\n";
    print $FILEHANDLE q?./configure --prefix=?.$filePath, "\n";
    print $FILEHANDLE "make", "\n";
    print $FILEHANDLE "make install", "\n";
    print $FILEHANDLE "\n\n";

    ## Move Perl Module
    print $FILEHANDLE "## Move Perl Module\n";
    print $FILEHANDLE q?cp src/perl/Vcf.pm $HOME/perl-?.${$parameterHashRef}{perl}.q?/lib/perl5/?;
    print $FILEHANDLE "\n\n";

    ## Remove the temporary install directory
    &RemoveInstallDirectory({FILEHANDLE => $FILEHANDLE,
			     pwd => $pwd,
			    });

    ## Reset perl envionment
    print $FILEHANDLE q?PERL5LIB=~/perl-?.${$parameterHashRef}{perl}.q?/lib/perl5?;
    print $FILEHANDLE "\n\n";
}


sub BedTools {

##BedTools
    
##Function : Install BedTools
##Returns  : ""
##Arguments: $parameterHashRef, $FILEHANDLE
##         : $parameterHashRef => Holds all parameters
##         : $FILEHANDLE       => Filehandle to write to

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $parameterHashRef;
    my $FILEHANDLE;

    my $tmpl = { 
	parameterHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameterHashRef},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];
    
    my $pwd = cwd();
    
    ## Check if the binary of the program being installed already exists
    if(&CheckCondaBinFileExists({parameterHashRef => $parameterHashRef,
				 programName => "bedtools",
				})) {
	
	return
    }

    my $bedToolsMainVersion = substr(${$parameterHashRef}{bedTools}, 0, 1);

    ## Install bedTools
    print $FILEHANDLE "### Install bedTools\n";
    
    ## Create the temporary install directory
    &CreateInstallDirectory({FILEHANDLE => $FILEHANDLE,
			    });
    
    ## Download
    print $FILEHANDLE "## Download bedTools\n";
    print $FILEHANDLE "wget --quiet https://github.com/arq5x/bedtools".$bedToolsMainVersion."/releases/download/v".${$parameterHashRef}{bedTools}."/bedtools-".${$parameterHashRef}{bedTools}.".tar.gz ";
    print $FILEHANDLE "-O bedtools-".${$parameterHashRef}{bedTools}.".tar.gz";  #Download outfile
    print $FILEHANDLE "\n\n";
    
    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "tar xvf bedtools-".${$parameterHashRef}{bedTools}.".tar.gz";
    print $FILEHANDLE "\n\n";

    ## Move to bedtools directory
    print $FILEHANDLE "## Move to bedtools directory\n";
    print $FILEHANDLE "cd bedtools".$bedToolsMainVersion;
    print $FILEHANDLE "\n\n";

    print $FILEHANDLE "make";
    print $FILEHANDLE "\n\n";
       
    ## Make available from conda environment
    print $FILEHANDLE "## Make available from conda environment\n";
    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?./bin/* ?.$parameter{condaPath}.q?/envs/?.${$parameterHashRef}{condaEnvironment}.q?/bin/?;
    print $FILEHANDLE "\n\n";
    
    ## Remove the temporary install directory
    &RemoveInstallDirectory({FILEHANDLE => $FILEHANDLE,
			     pwd => $pwd,
			    });
}


sub VT {

##VT
    
##Function : Install VT
##Returns  : ""
##Arguments: $parameterHashRef, $FILEHANDLE
##         : $parameterHashRef => Holds all parameters
##         : $FILEHANDLE       => Filehandle to write to

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $parameterHashRef;
    my $FILEHANDLE;

    my $tmpl = { 
	parameterHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameterHashRef},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Check if the binary of the program being installed already exists
    if (&CheckCondaBinFileExists({parameterHashRef => $parameterHashRef,
				  programName => "vt",
				 })) {
    
	return
    }

    ## Install VT
    print $FILEHANDLE "### Install VT\n";

    ## Create the temporary install directory
    &CreateInstallDirectory({FILEHANDLE => $FILEHANDLE,
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
    print $FILEHANDLE q?vt ?.$parameter{condaPath}.q?/envs/?.${$parameterHashRef}{condaEnvironment}.q?/bin/?;
    print $FILEHANDLE "\n\n";

    ## Remove the temporary install directory
    &RemoveInstallDirectory({FILEHANDLE => $FILEHANDLE,
			     pwd => $pwd,
			    });
}


sub Plink2 {

##Plink2
    
##Function : Install Plink2
##Returns  : ""
##Arguments: $parameterHashRef, $FILEHANDLE
##         : $parameterHashRef => Holds all parameters
##         : $FILEHANDLE       => Filehandle to write to

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $parameterHashRef;
    my $FILEHANDLE;

    my $tmpl = { 
	parameterHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameterHashRef},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Check if the binary of the program being installed already exists
    if (&CheckCondaBinFileExists({parameterHashRef => $parameterHashRef,
				  programName => "plink",
				 })) {

	return
    }

    ## Install Plink
    print $FILEHANDLE "### Install Plink\n";

    ## Create the temporary install directory
    &CreateInstallDirectory({FILEHANDLE => $FILEHANDLE,
			    });
    
    ## Download
    print $FILEHANDLE "## Download Plink\n";
    print $FILEHANDLE "wget --quiet https://www.cog-genomics.org/static/bin/plink".${$parameterHashRef}{plink2}."/plink_linux_x86_64.zip ";
    print $FILEHANDLE "-O plink-".${$parameterHashRef}{plink2}."-x86_64.zip";  #Download outfile
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "unzip plink-".${$parameterHashRef}{plink2}."-x86_64.zip";
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    print $FILEHANDLE "## Make available from conda environment\n";
    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?plink ?.$parameter{condaPath}.q?/envs/?.${$parameterHashRef}{condaEnvironment}.q?/bin/plink2?;
    print $FILEHANDLE "\n\n";

    ## Remove the temporary install directory
    &RemoveInstallDirectory({FILEHANDLE => $FILEHANDLE,
			     pwd => $pwd,
			    });
}


sub SnpEff {

##SnpEff
    
##Function : Install SnpEff
##Returns  : ""
##Arguments: $parameterHashRef, $FILEHANDLE
##         : $parameterHashRef => Holds all parameters
##         : $FILEHANDLE       => Filehandle to write to

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $parameterHashRef;
    my $FILEHANDLE;

    my $tmpl = { 
	parameterHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameterHashRef},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Check if the binary of the program being installed already exists
    if (&CheckCondaBinFileExists({parameterHashRef => $parameterHashRef,
				  programName => "snpEff.jar",
				 })) {  # Assumes that SnpSift.jar is there as well then

	return
    }

    ## Install SnpEff
    print $FILEHANDLE "### Install SnpEff\n";

    ## Create the temporary install directory
    &CreateInstallDirectory({FILEHANDLE => $FILEHANDLE,
			    });
    
    ## Download
    print $FILEHANDLE "## Download SnpEff\n";
    print $FILEHANDLE "wget --quiet http://sourceforge.net/projects/snpeff/files/snpEff_".${$parameterHashRef}{snpEff}."_core.zip/download ";
    print $FILEHANDLE "-O snpEff_".${$parameterHashRef}{snpEff}."_core.zip";  #Download outfile
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "unzip snpEff_".${$parameterHashRef}{snpEff}."_core.zip";
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    if (-d $parameter{condaPath}.q?/envs/?.${$parameterHashRef}{condaEnvironment}.q?/share/snpEff.?.${$parameterHashRef}{snpEff}) {
	
	print $FILEHANDLE "rm -rf ".$parameter{condaPath}.q?/envs/?.${$parameterHashRef}{condaEnvironment}.q?/share/snpEff.?.${$parameterHashRef}{snpEff};
	print $FILEHANDLE "\n\n";
    }

    print $FILEHANDLE "## Make available from conda environment\n";
    print $FILEHANDLE "mkdir -p ".$parameter{condaPath}.q?/envs/?.${$parameterHashRef}{condaEnvironment}.q?/share/snpEff.?.${$parameterHashRef}{snpEff};
    print $FILEHANDLE "\n\n";

    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?snpEff/*.jar ?.$parameter{condaPath}.q?/envs/?.${$parameterHashRef}{condaEnvironment}.q?/share/snpEff.?.${$parameterHashRef}{snpEff}.q?/?;
    print $FILEHANDLE "\n\n";

    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?snpEff/snpEff.config ?.$parameter{condaPath}.q?/envs/?.${$parameterHashRef}{condaEnvironment}.q?/share/snpEff.?.${$parameterHashRef}{snpEff}.q?/?;
    print $FILEHANDLE "\n\n";

    ## Define binaries
    my @snpEffBinaries = ("snpEff.jar",
			  "SnpSift.jar",
			  "snpEff.config",
	);
    
    foreach my $binary (@snpEffBinaries) {
	
	&CreateSoftLink({parameterHashRef => $parameterHashRef,
			 FILEHANDLE => $BASHFILEHANDLE,
			 binary => catfile($parameter{condaPath}, "envs", $parameter{condaEnvironment}, "share", "snpEff.".${$parameterHashRef}{snpEff}, $binary),
			 softLink => $binary,
			});
    }

    foreach my $genomeVersion (@{${$parameterHashRef}{snpEffGenomeVersions}}) {

	## Check and if required add the vertebrate mitochondrial codon table to SnpEff config
	&CheckMTCodonTable({parameterHashRef => $parameterHashRef,
			    FILEHANDLE => $BASHFILEHANDLE,
			    shareDirectory => catdir($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment}, "share", "snpEff.".${$parameterHashRef}{snpEff}),
			    configFile => "snpEff.config",
			    genomeVersionRef => \$genomeVersion,
			   });
	
	unless (-d catdir($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment}, "share", "snpEff.".${$parameterHashRef}{snpEff}, "data", $genomeVersion)) {
	    
	    ## Write instructions to download SnpEff database. This is done by install script to avoid race conditin when doing first analysis run in MIP
	    &SnpEffDownload({parameterHashRef => $parameterHashRef,
			     FILEHANDLE => $BASHFILEHANDLE,
			     genomeVersionRef => \$genomeVersion,
			    });
	}
    }

    ## Remove the temporary install directory
    &RemoveInstallDirectory({FILEHANDLE => $FILEHANDLE,
			     pwd => $pwd,
			    });
}


sub VariantEffectPredictor {

##VariantEffectPredictor
    
##Function : Install VariantEffectPredictor
##Returns  : ""
##Arguments: $parameterHashRef, $FILEHANDLE
##         : $parameterHashRef => Holds all parameters
##         : $FILEHANDLE       => Filehandle to write to

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $parameterHashRef;
    my $FILEHANDLE;

    my $tmpl = { 
	parameterHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameterHashRef},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();
    
    my $minicondaBinDirectory = catdir($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment}, "ensembl-tools-release-".${$parameterHashRef}{variantEffectPredictor});

    if (-d $minicondaBinDirectory) {

	print STDERR q?Found VariantEffectPredictor in miniconda directory: ?.$minicondaBinDirectory, "\n";
	
	if (${$parameterHashRef}{noUpdate}) {

	    print STDERR "Skipping writting installation process for VariantEffectPredictor\n";  	    
	    return
	}
	else {

	    ## Removing VariantEffectPredictor
	    print $FILEHANDLE "### Removing VariantEffectPredictor\n";
	    print $FILEHANDLE q?rm -rf ?.$minicondaBinDirectory;
	    print $FILEHANDLE "\n\n";
	}
    }
    else {
	
	print STDERR "Writting install instructions for VariantEffectPredictor\n";
    }

    ## Install VEP
    print $FILEHANDLE "### Install VariantEffectPredictor\n";

    ## Activate conda environment
    &ActivateCondaEnvironment({parameterHashRef => $parameterHashRef,
			       FILEHANDLE => $FILEHANDLE,
			      });

    ##Make sure that the cache directory exists
    print $FILEHANDLE "mkdir -p ".${$parameterHashRef}{vepDirectoryCache}." ";  #Cache directory
    print $FILEHANDLE "\n\n";

    ## Move to miniconda environment
    print $FILEHANDLE q?cd ?.catdir($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment});
    print $FILEHANDLE "\n\n";

    ## Download
    print $FILEHANDLE "## Download VEP\n";
    print $FILEHANDLE "wget --quiet https://github.com/Ensembl/ensembl-tools/archive/release/".${$parameterHashRef}{variantEffectPredictor}.".zip ";
    print $FILEHANDLE "-O VariantEffectPredictor-".${$parameterHashRef}{variantEffectPredictor}.".zip";  #Download outfile
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "unzip VariantEffectPredictor-".${$parameterHashRef}{variantEffectPredictor}.".zip";
    print $FILEHANDLE "\n\n";    

    ## Move to VariantEffectPredictor directory
    print $FILEHANDLE "## Move to VariantEffectPredictor directory\n";
    print $FILEHANDLE "cd ".catdir("ensembl-tools-release-".${$parameterHashRef}{variantEffectPredictor}, "scripts", "variant_effect_predictor");
    print $FILEHANDLE "\n\n";

    ## Install VEP
    print $FILEHANDLE "## Install VEP\n";
    print $FILEHANDLE "perl INSTALL.pl ";
    print $FILEHANDLE "--AUTO alcfp ";  #a (API), l (FAIDX/htslib), c (cache), f (FASTA), p (plugins)
    print $FILEHANDLE "-g ".$parameter{variantEffectPredictorPlugin}." ";  #Plugins 
    print $FILEHANDLE "-c ".${$parameterHashRef}{vepDirectoryCache}." ";  #Cache directory
    print $FILEHANDLE "-s homo_sapiens ";

    ## Only install first assembly version since VEP install cannot handle multiple versions at the same time
    print $FILEHANDLE "--ASSEMBLY ".${$parameterHashRef}{vepAssemblies}[0]." ";
    print $FILEHANDLE "\n\n";
    
    if (scalar( @{${$parameterHashRef}{vepAssemblies}} ) > 1 ) {

	for (my $assemblyVersion=1;$assemblyVersion<scalar(@{${$parameterHashRef}{vepAssemblies}});$assemblyVersion++) {

	    print $FILEHANDLE "## Install additional VEP cache assembly version\n";
	    print $FILEHANDLE "perl INSTALL.pl ";
	    print $FILEHANDLE "--AUTO cf ";  #a (API), l (FAIDX/htslib), c (cache), f (FASTA), p (plugins)
	    print $FILEHANDLE "-c ".${$parameterHashRef}{vepDirectoryCache}." ";  #Cache directory
	    print $FILEHANDLE "-s homo_sapiens ";
	    
	    ## Only install first assembly version since VEP install cannot handle multiple versions at the same time
	    print $FILEHANDLE "--ASSEMBLY ".${$parameterHashRef}{vepAssemblies}[$assemblyVersion]." ";
	    print $FILEHANDLE "\n\n";
	}
    }

    ##Add LofTool required text file
    print $FILEHANDLE "##Add LofTool required text file\n";
    print $FILEHANDLE "wget https://raw.githubusercontent.com/Ensembl/VEP_plugins/master/LoFtool_scores.txt ";
    print $FILEHANDLE q?-O $HOME/.vep/Plugins/LoFtool_scores.txt ?;
    print $FILEHANDLE "\n\n";

    ##Add Lof required perl splice script
    print $FILEHANDLE "##Add Lof required perl splice script\n";
    print $FILEHANDLE "wget https://raw.githubusercontent.com/konradjk/loftee/master/splice_module.pl ";
    print $FILEHANDLE q?-O $HOME/.vep/Plugins/splice_module.pl ?;
    print $FILEHANDLE "\n\n";

    ##Add Lof optional human_ancestor_fa
    print $FILEHANDLE "##Add Lof optional human_ancestor_fa\n";
    print $FILEHANDLE "wget https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz ";
    print $FILEHANDLE q?-O ?.catfile(${$parameterHashRef}{vepDirectoryCache}, "human_ancestor.fa.gz")." ";
    print $FILEHANDLE "\n\n";

    ##Uncompress
    print $FILEHANDLE "bgzip -d ".catfile(${$parameterHashRef}{vepDirectoryCache}, "human_ancestor.fa.gz")." ";
    print $FILEHANDLE "\n\n";

    print $FILEHANDLE "wget https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz.fai ";
    print $FILEHANDLE q?-O ?.catfile(${$parameterHashRef}{vepDirectoryCache}, "human_ancestor.fa.fai")." ";
    print $FILEHANDLE "\n\n";

    ## Clean up
    print $FILEHANDLE "## Clean up\n";
    print $FILEHANDLE q?cd ?.catdir($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment});
    print $FILEHANDLE "\n\n";

    print $FILEHANDLE "rm -rf VariantEffectPredictor-".${$parameterHashRef}{variantEffectPredictor}.".zip";;
    print $FILEHANDLE "\n\n";

    ## Moving up
    print $FILEHANDLE "## Moving back to original working directory\n";
    print $FILEHANDLE "cd ".$pwd;  #Go back to subroutine origin
    print $FILEHANDLE "\n\n";

    ## Deactivate conda environment
    &DeactivateCondaEnvironment({FILEHANDLE => $FILEHANDLE,
				});
}


sub CNVnator {

##CNVnator
    
##Function : Install CNVnator
##Returns  : ""
##Arguments: $parameterHashRef, $FILEHANDLE
##         : $parameterHashRef => Holds all parameters
##         : $FILEHANDLE       => Filehandle to write to

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $parameterHashRef;
    my $FILEHANDLE;

    my $tmpl = { 
	parameterHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameterHashRef},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Check if the binary of the program being installed already exists
    if (&CheckCondaBinFileExists({parameterHashRef => $parameterHashRef,
				  programName => "cnvnator",
				 })) {

	return
    }

    my $minicondaBinDirectory = catdir($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment}, "root");

    if (-d $minicondaBinDirectory) {

	print STDERR q?Found Root in miniconda directory: ?.$minicondaBinDirectory, "\n";
	
	if (${$parameterHashRef}{noUpdate}) {

	    print STDERR "Skipping writting installation process for Root\n";  	    
	    return
	}
	else {

	    ## Removing Root
	    print $FILEHANDLE "### Removing Root\n";
	    print $FILEHANDLE q?rm -rf ?.$minicondaBinDirectory;
	    print $FILEHANDLE "\n\n";
	}
    }
    else {
	
	print STDERR "Writting install instructions for Root\n";
    }

    ## Install Root
    print $FILEHANDLE "### Install CNVnator/Root\n";

    ## Move to miniconda environment
    print $FILEHANDLE q?cd ?.catdir($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment});
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

    unless ($ENV{PATH}=~/$parameter{condaPath}\/envs\/${$parameterHashRef}{condaEnvironment}\/root\/bin/) {
	
	## Export path
	print $FILEHANDLE "## Export path\n";
	print $FILEHANDLE q?echo 'source ?.$parameter{condaPath}.q?/envs/?.${$parameterHashRef}{condaEnvironment}.q?/root/bin/thisroot.sh' >> ~/.bashrc?;
	print $FILEHANDLE "\n\n";

	## Use newly installed root
	print $FILEHANDLE q?source ?.$parameter{condaPath}.q?/envs/?.${$parameterHashRef}{condaEnvironment}.q?/root/bin/thisroot.sh ?;
	print $FILEHANDLE "\n\n";
    }
    
    ## Moving up
    print $FILEHANDLE "## Moving back to original working directory\n";
    print $FILEHANDLE "cd ".$pwd;  #Go back to subroutine origin
    print $FILEHANDLE "\n\n";

    ## Install CNVNator
    print $FILEHANDLE "### Install CNVnator\n";

    ## Activate conda environment
    &ActivateCondaEnvironment({parameterHashRef => $parameterHashRef,
			       FILEHANDLE => $FILEHANDLE,
			      });

    ## Create the temporary install directory
    &CreateInstallDirectory({FILEHANDLE => $FILEHANDLE,
			    });
    
    ## Download
    print $FILEHANDLE "## Download CNVNator\n";
    print $FILEHANDLE "wget --quiet https://github.com/abyzovlab/CNVnator/releases/download/v".${$parameterHashRef}{CNVnator}."/CNVnator_v".${$parameterHashRef}{CNVnator}.".zip ";
    print $FILEHANDLE "-O CNVnator_v".${$parameterHashRef}{CNVnator}.".zip";  #Download outfile
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "unzip CNVnator_v".${$parameterHashRef}{CNVnator}.".zip";
    print $FILEHANDLE "\n\n";

    ## Move to CNVnator directory
    print $FILEHANDLE "## Move to CNVnator directory\n";
    print $FILEHANDLE "cd ".catdir("CNVnator_v".${$parameterHashRef}{CNVnator}, "src", "samtools");
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
    print $FILEHANDLE q?cnvnator ?.catdir($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment}, "bin");
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    print $FILEHANDLE "## Make available from conda environment\n";
    print $FILEHANDLE "cd ..";
    print $FILEHANDLE "\n";
    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?cnvnator2VCF.pl ?.catdir($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment}, "bin");
    print $FILEHANDLE "\n\n";

    print $FILEHANDLE "## Make executable from conda environment\n";
    print $FILEHANDLE "chmod +x ";
    print $FILEHANDLE catfile($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment}, "bin", "cnvnator2VCF.pl");
    print $FILEHANDLE "\n\n";
    
    ## Remove the temporary install directory
    &RemoveInstallDirectory({FILEHANDLE => $FILEHANDLE,
			     pwd => $pwd,
			    });
    ## Deactivate conda environment
    &DeactivateCondaEnvironment({FILEHANDLE => $FILEHANDLE,
				});
}


sub FindTranslocations {

##FindTranslocations
    
##Function : Install FindTranslocations
##Returns  : ""
##Arguments: $parameterHashRef, $FILEHANDLE
##         : $parameterHashRef => Holds all parameters
##         : $FILEHANDLE       => Filehandle to write to

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $parameterHashRef;
    my $FILEHANDLE;

    my $tmpl = { 
	parameterHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameterHashRef},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Check if the binary of the program being installed already exists
    if (&CheckCondaBinFileExists({parameterHashRef => $parameterHashRef,
				  programName => "FindTranslocations",
				 })) {

	return
    }

    ## Install FindTranslocations
    print $FILEHANDLE "### Install FindTranslocations\n";

    ## Activate conda environment
    &ActivateCondaEnvironment({parameterHashRef => $parameterHashRef,
			       FILEHANDLE => $FILEHANDLE,
			      });

    ## Add to bashrc
    unless (-d catdir($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment}, "FindTranslocations", "bin")) {
	
	## Export path
	print $FILEHANDLE "## Export to bashrc\n";
	print $FILEHANDLE q?printf '\nif [ -f ?.$parameter{condaPath}.q?/envs/?.${$parameterHashRef}{condaEnvironment}.q?/FindTranslocations/bin/FindTranslocations ]; then\n?;
	print $FILEHANDLE q?\t\texport LD_LIBRARY_PATH=$LD_LIBRARY_PATH:?.$parameter{condaPath}.q?/pkgs/boost-?.$parameter{bioConda}{boost}.$parameter{bioCondaBoostPatch}.q?/lib\n?;
	print $FILEHANDLE q?fi\n\n' >> ~/.bashrc?;
	print $FILEHANDLE "\n\n";
    }

    ## Move to miniconda environment
    print $FILEHANDLE "cd ".catdir($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment});
    print $FILEHANDLE "\n\n";

    ## Download
    print $FILEHANDLE "## Download FindTranslocations\n";
    print $FILEHANDLE "wget --quiet https://github.com/J35P312/FindTranslocations/archive/version_".${$parameterHashRef}{FindTranslocations}.".zip ";
    print $FILEHANDLE "-O FindTranslocations-".${$parameterHashRef}{FindTranslocations}.".zip";  #Download outfile
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "rm -rf FindTranslocations";
    print $FILEHANDLE "\n\n";
    print $FILEHANDLE "unzip FindTranslocations-".${$parameterHashRef}{FindTranslocations}.".zip ";
    print $FILEHANDLE "\n\n";
    print $FILEHANDLE "mv FindTranslocations-version_".${$parameterHashRef}{FindTranslocations}." ";
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
    print $FILEHANDLE catfile($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment}, "FindTranslocations", "bin", "FindTranslocations")." ".catdir($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment}, "bin");
    print $FILEHANDLE "\n\n";    

    ## Clean-up
    print $FILEHANDLE "cd ".catdir($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment});
    print $FILEHANDLE "\n\n";
    print $FILEHANDLE "rm -rf FindTranslocations-".${$parameterHashRef}{FindTranslocations}.".zip";
    print $FILEHANDLE "\n\n";

    ## Deactivate conda environment
    &DeactivateCondaEnvironment({FILEHANDLE => $FILEHANDLE,
				});
}


sub MIPScripts {

##MIPScripts
    
##Function : Install MIPScripts 
##Returns  : ""
##Arguments: $parameterHashRef, $FILEHANDLE
##         : $parameterHashRef => Holds all parameters
##         : $FILEHANDLE       => Filehandle to write to

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $parameterHashRef;
    my $FILEHANDLE;

    my $tmpl = { 
	parameterHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameterHashRef},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();
    my $mipDirectory = $Bin;  #Alias

    ## Define MIP scripts and yaml files
    my @mipScripts = ("calculateAF.pl",
		      "maxAF.pl",
		      "mip.pl",
		      "qcCollect.pl",
		      "vcfParser.pl",
		      "install.pl",
	);
    my %mipSubScripts;
    $mipSubScripts{"definitions"} = ["defineParameters.yaml"];
    $mipSubScripts{"t"} = ["test.t"];
    $mipSubScripts{"templates"} = ["mip_config.yaml"];

    ## Check if the binary of the program being installed already exists
    if (&CheckCondaBinFileExists({parameterHashRef => $parameterHashRef,
				  programName => "mip.pl",
				 })) {  #Proxy for all 

	return
    }

    ## Install MIPScripts
    print $FILEHANDLE "### Install MIPScripts\n";
    
    ## Create directories
    print $FILEHANDLE "## Create directories\n";
    foreach my $directory (keys %mipSubScripts) {

	print $FILEHANDLE "mkdir -p ";
	print $FILEHANDLE catdir($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment}, "bin", $directory);
	print $FILEHANDLE "\n\n";
    }

    ## Copy mip scripts and sub scripts to conda env and make executable
    print $FILEHANDLE "## Copy mip scripts and subdirectory scripts to conda env and make executable\n\n";
    foreach my $script (@mipScripts) {
	    
	print $FILEHANDLE "cp ";
	print $FILEHANDLE catfile($Bin, $script)." ";
	print $FILEHANDLE catdir($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment}, "bin"), "\n";
	print $FILEHANDLE "chmod a+x ".catfile($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment}, "bin", $script);
	print $FILEHANDLE "\n\n";
    }

    foreach my $directory (keys %mipSubScripts) {

	foreach my $script (@{$mipSubScripts{$directory}}) {

	    print $FILEHANDLE "cp ";
	    print $FILEHANDLE catfile($Bin, $directory, $script)." ";
	    print $FILEHANDLE catdir($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment}, "bin", $directory), "\n";
	    print $FILEHANDLE "chmod a+x ".catfile($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment}, "bin", $directory, $script);
	    print $FILEHANDLE "\n\n";
	}
    }
}


sub ActivateCondaEnvironment {

##ActivateCondaEnvironment

##Function : Activate conda environment
##Returns  : ""
##Arguments: $parameterHashRef, $FILEHANDLE
##         : $parameterHashRef => Holds all parameters
##         : $FILEHANDLE       => Filehandle to write to

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $parameterHashRef;
    my $FILEHANDLE;

    my $tmpl = { 
	parameterHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameterHashRef},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    ## Activate conda environment
    print $FILEHANDLE "## Activate conda environment\n";
    print $FILEHANDLE "source activate ".${$parameterHashRef}{condaEnvironment}." ";
    print $FILEHANDLE "\n\n";
}


sub DeactivateCondaEnvironment {

##DeactivateCondaEnvironment
    
##Function : Deactivate conda environment
##Returns  : ""
##Arguments: $FILEHANDLE
##         : $FILEHANDLE => Filehandle to write to

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;

    my $tmpl = { 
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    ## Deactivate conda environment
    print $FILEHANDLE "## Deactivate conda environment\n";
    print $FILEHANDLE "source deactivate ";
    print $FILEHANDLE "\n\n";
}


sub RemoveInstallDirectory {

##RemoveInstallDirectory
    
##Function : Remove the temporary install directory
##Returns  : ""
##Arguments: $FILEHANDLE, $pwd
##         : $FILEHANDLE => FILEHANDLE to write to
##         : $pwd        => The original working directory

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $pwd;

    my $tmpl = { 
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
	pwd => { required => 1, defined => 1, strict_type => 1, store => \$pwd},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    ## Moving up
    print $FILEHANDLE "## Moving back to original working directory\n";
    print $FILEHANDLE "cd ".$pwd;  #Go back to subroutine origin
    print $FILEHANDLE "\n\n";

    ## Clean up
    print $FILEHANDLE "## Clean up\n";
    print $FILEHANDLE "rm -rf .MIP";
    print $FILEHANDLE "\n\n";
}

sub CreateInstallDirectory {

##CreateInstallDirectory
    
##Function : Create the temporary install directory
##Returns  : ""
##Arguments: $FILEHANDLE
##         : $FILEHANDLE => FILEHANDLE to write to

    my ($argHashRef) = @_;
    
    ## Flatten argument(s)
    my $FILEHANDLE;
    
    my $tmpl = { 
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];
    
    ## Create temp install directory
    print $FILEHANDLE "## Create temp install directory\n";
    print $FILEHANDLE "mkdir -p .MIP ", "\n";
    print $FILEHANDLE "cd .MIP";
    print $FILEHANDLE "\n\n";
}


sub CheckCondaBinFileExists {
    
##CheckCondaBinFileExists
    
##Function : Check if the binary of the program being installed already exists
##Returns  : ""
##Arguments: $parameterHashRef, $programName, $programVersion
##         : $parameterHashRef => Holds all parameters
##         : $programName      => Program name
##         : $programVersion   => Program version

    my ($argHashRef) = @_;

    ## Default(s)
    my $programVersion;

    ## Flatten argument(s)
    my $parameterHashRef;
    my $programName;
    
    my $tmpl = { 
	parameterHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameterHashRef},
	programName => { required => 1, defined => 1, strict_type => 1, store => \$programName},	
	programVersion => { default => 0,
			    strict_type => 1, store => \$programVersion},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my $minicondaBinFile;
    
    if ($programVersion) {

	$minicondaBinFile = catfile($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment}, "bin", $programName.$programVersion);
    }
    else {

	$minicondaBinFile = catfile($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment}, "bin", $programName);
    }

    if (-f $minicondaBinFile) {

	if ($programVersion) {
	    
	    print STDERR q?Found ?.$programName.q? version ?.$programVersion.q? in miniconda directory: ?.catdir($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment}, "bin"), "\n";
	    
	    if (${$parameterHashRef}{noUpdate}) {

		print STDERR q?Skipping writting installation process for ?.$programName.q? ?.$programVersion, "\n";  
		return 1;
	    }
	    print STDERR "Writting install instructions for ".$programName, "\n";
	}   
	else {

	    print STDERR q?Found ?.$programName.q? in miniconda directory: ?.catdir($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment}, "bin"), "\n";
	    
	    if (${$parameterHashRef}{noUpdate}) {
		
		print STDERR q?Skipping writting installation process for ?.$programName, "\n";  	    
		return 1;
	    }
	    print STDERR "Writting install instructions for ".$programName, "\n";
	}
	return 0;
    }
    else {
	
	print STDERR "Writting install instructions for ".$programName, "\n";
	return 0;
    }
}


sub CreateSoftLink {

##CreateSoftLink
    
##Function : Create softlink
##Returns  : ""
##Arguments: $parameterHashRef, $FILEHANDLE, $binary, $softLink
##         : $parameterHashRef => Holds all parameters
##         : $FILEHANDLE       => FILEHANDLE to write to
##         : $binary           => The binary file
##         : $softLink         => The name of the softlink

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $parameterHashRef;
    my $FILEHANDLE;
    my $binary;
    my $softLink;
    
    my $tmpl = { 
	parameterHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameterHashRef},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
	binary => { required => 1, defined => 1, strict_type => 1, store => \$binary},
	softLink => { required => 1, defined => 1, strict_type => 1, store => \$softLink},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Add softlink
    print $FILEHANDLE "## Move to directory and create softlink\n";
    print $FILEHANDLE "cd ".catdir($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment}, "bin"), "\n";

    print $FILEHANDLE "ln -f -s ";
    print $FILEHANDLE $binary.q? ?.$softLink;
    print $FILEHANDLE "\n\n";

    ## Move to back
    print $FILEHANDLE "## Move to original working directory\n";
    print $FILEHANDLE "cd ".$pwd;
    print $FILEHANDLE "\n\n";
}


sub EnableExecutable {

##EnableExecutable
    
##Function : Make file executable
##Returns  : ""
##Arguments: $parameterHashRef, $FILEHANDLE, $binary
##         : $parameterHashRef => Holds all parameters
##         : $FILEHANDLE       => FILEHANDLE to write to
##         : $binary           => The binary file

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $parameterHashRef;
    my $FILEHANDLE;
    my $binary;
    
     my $tmpl = { 
	parameterHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameterHashRef},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
	binary => { required => 1, defined => 1, strict_type => 1, store => \$binary},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Enable executable
    print $FILEHANDLE "## Enable executable\n";
    print $FILEHANDLE "cd ".catdir($parameter{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment}, "bin"), "\n";

    print $FILEHANDLE "chmod a+x ";
    print $FILEHANDLE $binary.q? ?;
    print $FILEHANDLE "\n\n";

    ## Move to back
    print $FILEHANDLE "## Move to original working directory\n";
    print $FILEHANDLE "cd ".$pwd;
    print $FILEHANDLE "\n\n";
}


sub CheckMTCodonTable {

##CheckMTCodonTable
    
##Function : Check and if required add the vertebrate mitochondrial codon table to SnpEff config
##Returns  : ""
##Arguments: $parameterHashRef, $FILEHANDLE, $shareDirectory, $configFile, $genomeVersionRef
##         : $parameterHashRef => Holds all parameters
##         : $FILEHANDLE       => FILEHANDLE to write to
##         : $shareDirectory   => The conda env shared directory
##         : $configFile       => The config configFile
##         : $genomeVersionRef => SnpEff genome version

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $parameterHashRef;
    my $FILEHANDLE;
    my $shareDirectory;
    my $configFile;
    my $genomeVersionRef;
    
    my $tmpl = { 
	parameterHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameterHashRef},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
	shareDirectory => { required => 1, defined => 1, strict_type => 1, store => \$shareDirectory},
	configFile => { required => 1, defined => 1, strict_type => 1, store => \$configFile},
	genomeVersionRef => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$genomeVersionRef},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my $pwd = cwd();
    my $detectRegExp = q?perl -nae 'if($_=~/?.$$genomeVersionRef.q?.MT.codonTable/) {print 1}' ?;
    my $addRegExp = q?perl -nae 'if($_=~/?.$$genomeVersionRef.q?.reference/) {print $_; print "?.$$genomeVersionRef.q?.MT.codonTable : Vertebrate_Mitochondrial\n"} else {print $_;}' ?;
    my $ret;

    if (-f catfile($shareDirectory, $configFile)) {

	$ret = `$detectRegExp $shareDirectory/$configFile`;
    }
    if (!$ret) {  #No MT.codonTable in config

	print $FILEHANDLE q?## Adding ?.$$genomeVersionRef.q?.MT.codonTable : Vertebrate_Mitochondrial to ?.$shareDirectory.$configFile, "\n";

	## Add MT.codon Table to config
	print $FILEHANDLE $addRegExp." ".catfile($shareDirectory, $configFile)." > ".catfile($shareDirectory, $configFile.".tmp"), "\n";
	print $FILEHANDLE "mv ".catfile($shareDirectory, $configFile.".tmp")." ".catfile($shareDirectory, $configFile);
	print $FILEHANDLE "\n\n";
	
    }
    else {

	print STDERR  "Found MT.codonTable in ".catfile($shareDirectory, "snpEff.config").". Skipping addition to snpEff config", "\n"; 
    }
}


sub SnpEffDownload {

##SnpEffDownload
    
##Function : Write instructions to download SnpEff database. This is done by install script to avoid race conditin when doing first analysis run in MIP
##Returns  : ""
##Arguments: $parameterHashRef, $FILEHANDLE, $genomeVersionRef
##         : $parameterHashRef => Holds all parameters
##         : $FILEHANDLE       => FILEHANDLE to write to
##         : $genomeVersionRef => SnpEff genome version

    my ($argHashRef) = @_;
    
    ## Flatten argument(s)
    my $parameterHashRef;
    my $FILEHANDLE;
    my $genomeVersionRef;

    my $tmpl = { 
	parameterHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$parameterHashRef},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
	genomeVersionRef => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$genomeVersionRef},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    ## Activate conda environment
    &ActivateCondaEnvironment({parameterHashRef => $parameterHashRef,
			       FILEHANDLE => $FILEHANDLE,
			      });
    
    print $FILEHANDLE "java -Xmx2g ";
    print $FILEHANDLE q?-jar ?.catfile(${$parameterHashRef}{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment}, "bin", "snpEff.jar")." ";
    print $FILEHANDLE "download ";
    print $FILEHANDLE " -v ";
    print $FILEHANDLE $$genomeVersionRef." ";
    print $FILEHANDLE q?-c ?.catfile(${$parameterHashRef}{condaPath}, "envs", ${$parameterHashRef}{condaEnvironment}, "bin", "snpEff.config")." ";
    print $FILEHANDLE "\n\n";

    ## Deactivate conda environment
    &DeactivateCondaEnvironment({FILEHANDLE => $FILEHANDLE,
				});
    
}


sub Help {

##Help
    
##Function : Print help text and exit with supplied exit code
##Returns  : ""
##Arguments: $USAGE, $exitCode
##         : $USAGE    => Help text
##         : $exitCode => Exit code

    my ($argHashRef) = @_;

    ## Default(s)
    my $exitCode;

    ## Flatten argument(s)
    my $USAGE;

    my $tmpl = { 
	USAGE => { required => 1, defined => 1, strict_type => 1, store => \$USAGE},
	exitCode => { default => 0,
		      allow => qr/^\d+$/,
		      strict_type => 1, store => \$exitCode},
    };
 
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    print STDOUT $USAGE, "\n";
    exit $exitCode;
}
