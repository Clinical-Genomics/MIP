#!/usr/bin/env perl

##Assumes you have a working conda installation

use strict;
use warnings;

use Getopt::Long;
use IO::File;
use Cwd;
use FindBin qw($Bin); #Find directory of script
use vars qw($USAGE);
use IPC::Cmd qw[can_run run];

## Third party module(s)
use List::MoreUtils qw(any);

BEGIN {
    $USAGE =
	qq{install.pl [options]
           -env/--condaEnvironment Conda environment (Default: "mip")
           -cdp/--condaPath The conda path (Default: "HOME/miniconda")
           -bvc/--bioConda Set the module version of the programs that can be installed with bioConda (e.g. 'bwa=0.7.12')
           -pip/--pip Set the module version of the programs that can be installed with pip (e.g. 'genmod=3.4.2')

           ## SHELL
           -pei/--perlInstall Install perl (defaults: "0" (=no))
           -per/--perl Set the perl version (defaults: "5.18.2")
           -pm/perlModules Set the perl modules to be installed via cpanm (comma sep)
           -pic/--picardTools Set the picardTools version (Default: "2.0.1"),
           -sbb/sambamba Set the sambamba version (Default: "0.5.9")
           -vct/--vcfTools Set the vcftools version (Default: "0.1.14")
           -bet/--bedTools Set the bedtools version (Default: "2.25.0")
           -vt/--vt Set the vt version (Default: "0.57")
           -plk/--plink  Set the plink version (Default: "160224")
           -vep/--variantEffectPredictor Set the VEP version (Default: "83")
	   -vepc/--vepDirectoryCache Specify the cache directory to use (whole path; defaults to "~/miniconda/envs/condaEnvironment/ensembl-tools-release-variantEffectPredictorVersion/cache")
           -vepp/--variantEffectPredictorPlugin Supply a comma separated list of VEP plugins (Default: "UpDownDistance")
           -cnv/--CNVnator Set the CNVnator version (Default: "0.3.2")
           -ftr/--FindTranslocations Set the FindTranslocations version (Default: "X.X.X")

           ## Utility
           -pbc/--preferBioConda Bioconda will used for overlapping shell and biconda installations (Default: "1" (=yes))
           -ppd/--printParameterDefaults Print the parameter defaults
           -u/--update Always install all programs (Default: "1" (=yes))
           -sp/--selectPrograms Install supplied programs e.g. -sp perl -sp bedTools (Default: "";)
           -h/--help Display this help message   
           -v/--version Display version
        };    
}
my ($installDirectory) = (".MIP");

my %parameter; 

### Set parameter default

$parameter{'update'} = 1;

##Conda
$parameter{'preferBioConda'} = 1;
$parameter{'condaEnvironment'} = "mip";
$parameter{'condaPath'} = $ENV{HOME}."/miniconda";

$parameter{'bioConda'}{'bwa'} = "0.7.12";
$parameter{'bioConda'}{'fastqc'} = "0.11.4";
$parameter{'bioConda'}{'cramtools'} = "3.0.b47";
$parameter{'bioConda'}{'samtools'} = "1.3";
$parameter{'bioConda'}{'bcftools'} = "1.3";
$parameter{'bioConda'}{'snpeff'} = "4.2";
$parameter{'bioCondaSnpeffPatch'} = "-0";  #For correct softlinking in share and bin in conda env
$parameter{'bioConda'}{'picard'} = "1.141";
$parameter{'bioCondaPicardPatch'} = "-1";  #For correct softlinking in share and bin in conda env
$parameter{'bioConda'}{'mosaik'} = "2.2.26";
$parameter{'bioConda'}{'htslib'} = "1.3";
$parameter{'bioConda'}{'bedtools'} = "2.25.0";
$parameter{'bioConda'}{'vt'} = "2015.11.10";
$parameter{'bioConda'}{'sambamba'} = "0.5.9";
$parameter{'bioConda'}{'freebayes'} = "1.0.2";
$parameter{'bioConda'}{'delly'} = "0.7.2";
$parameter{'bioConda'}{'manta'} = "0.29.3";
$parameter{'bioConda'}{'gcc'} = "4.8.5";


##Perl Modules
$parameter{'perlInstall'} = 0;
$parameter{'perl'} = "5.18.2";
$parameter{'perlModules'} = ["YAML",
			     "Log::Log4perl",
			     "List::MoreUtils",
			     "DateTime",
			     "DateTime::Format::ISO8601",
			     "DateTime::Format::HTTP",
			     "DateTime::Format::Mail",
			     "Set::IntervalTree",  # vcfParser
			     "Net::SSLeay",  # VEP
			     "LWP::Simple",  # VEP
			     "LWP::Protocol::https",  # VEP
			     "Archive::Zip",  # VEP
			     "DBI",  # VEP
			     "JSON",  # VEP
			     "DBD::mysql",  # VEP
			     "CGI",  # VEP
			     "Sereal::Encoder",  # VEP
			     "Sereal::Decoder",  # VEP
    ];

## PIP
$parameter{'pip'}{'genmod'} = "3.4.5";
$parameter{'pip'}{'chanjo'} = "3.3.1";
$parameter{'pip'}{'cosmid'} = "0.4.9.1";
$parameter{'pip'}{'python-Levenshtein'} = "0.12.0";

## Programs currently installable by SHELL
$parameter{'picardTools'} = "2.0.1";
$parameter{'sambamba'} = "0.5.9";
$parameter{'vcfTools'} = "0.1.14";
$parameter{'bedTools'} = "2.25.0";
$parameter{'vt'} = "gitRepo";
$parameter{'plink'} = "160224";
$parameter{'snpEff'} = "v4_2";
$parameter{'variantEffectPredictor'} = "83";
$parameter{'vepDirectoryCache'} = $parameter{'condaPath'}.q?/envs/?.$parameter{'condaEnvironment'}.q?/ensembl-tools-release-?.$parameter{'variantEffectPredictor'}.q?/cache?;  #Cache directory;
$parameter{'variantEffectPredictorPlugin'} = "UpDownDistance";
$parameter{'CNVnator'} = "0.3.2";
$parameter{'FindTranslocations'} = "X";

my $installVersion = "0.0.4";

###User Options
GetOptions('env|condaEnvironment:s'  => \$parameter{'condaEnvironment'},
	   'cdp|condaPath:s' => \$parameter{'condaPath'},
	   'bcv|bioConda=s' => \%{$parameter{'bioConda'}},
	   'pip|pip=s' => \%{$parameter{'pip'}},
	   'per|perl=s' => \$parameter{'perl'},
	   'pei|perlInstall:n' => \$parameter{'perlInstall'},
	   'pm|perlModules:s' => \@{$parameter{'perlModules'}},  #Comma separated list
	   'pic|picardTools:s' => \$parameter{'picardTools'},
	   'sbb|sambamba:s' => \$parameter{'sambamba'},
	   'vct|vcfTools:s' => \$parameter{'vcfTools'},
	   'bet|bedTools:s' =>\$parameter{'bedTools'}, 
	   'vt|vt:s' => \$parameter{'vt'},
	   'plk|plink:s' => \$parameter{'plink'},
	   'vep|variantEffectPredictor:s' => \$parameter{'variantEffectPredictor'},
	   'vepc|vepDirectoryCache:s' => \$parameter{'vepDirectoryCache'},  #path to vep cache dir
	   'vepp|variantEffectPredictorPlugin:s' => \$parameter{'variantEffectPredictorPlugin'},  #Comma sep string
	   'cnv|CNVnator:s' => \$parameter{'CNVnator'},
	   'ftr|FindTranslocations:s' => \$parameter{'FindTranslocations'},
	   'pbc|preferBioConda:s' => \$parameter{'preferBioConda'},  # Bioconda will used for overlapping shell and biconda installlations
	   'ppd|printParameterDefaults' => sub { &PrintParameters(\%parameter); exit;},  #Display parameter defaults
	   'u|update:n' => \$parameter{'update'},
	   'sp|selectPrograms:s' => \@{$parameter{'selectPrograms'}},  #Comma sep string
	   'h|help' => sub { print STDOUT $USAGE, "\n"; exit;},  #Display help text
	   'v|version' => sub { print STDOUT "\ninstall.pl ".$installVersion, "\n\n"; exit;},  #Display version number
    );

###MAIN###

#my $LOGFILEHANDLE = &OpenLogFile("MIP_installation.log");

my $BASHFILEHANDLE = &CreateBashFile("mip.sh");

&CreateConda(\%parameter, $BASHFILEHANDLE);

&CreateCondaEnvironment(\%parameter, $BASHFILEHANDLE);

if (scalar(@{$parameter{'selectPrograms'}}) > 0) {
    
    if ( ( any {$_ eq "perl"} @{$parameter{'selectPrograms'}} ) ) { #If element is part of array
	
	&Perl(\%parameter, $BASHFILEHANDLE);
    }
}
else {
    
    &Perl(\%parameter, $BASHFILEHANDLE);
}


&PipInstall(\%parameter, $BASHFILEHANDLE);

if ($parameter{'preferBioConda'} != 1) {

    if (scalar(@{$parameter{'selectPrograms'}}) > 0) {
	
	if ( ( any {$_ eq "picardTools"} @{$parameter{'selectPrograms'}} ) ) { #If element is part of array
	    
	    &PicardTools(\%parameter, $BASHFILEHANDLE);
	}
	if ( ( any {$_ eq "sambamba"} @{$parameter{'selectPrograms'}} ) ) { #If element is part of array

	    &Sambamba(\%parameter, $BASHFILEHANDLE);
	}
	if ( ( any {$_ eq "bedTools"} @{$parameter{'selectPrograms'}} ) ) { #If element is part of array

	    &BedTools(\%parameter, $BASHFILEHANDLE);
	}
	if ( ( any {$_ eq "vt"} @{$parameter{'selectPrograms'}} ) ) { #If element is part of array

	    &VT(\%parameter, $BASHFILEHANDLE);
	}
	if ( ( any {$_ eq "snpEff"} @{$parameter{'selectPrograms'}} ) ) { #If element is part of array

	    &SnpEff(\%parameter, $BASHFILEHANDLE);
	}
    }
    else {
	
	&PicardTools(\%parameter, $BASHFILEHANDLE);
	
	&Sambamba(\%parameter, $BASHFILEHANDLE);
	
	&BedTools(\%parameter, $BASHFILEHANDLE);
	
	&VT(\%parameter, $BASHFILEHANDLE);
	
	&SnpEff(\%parameter, $BASHFILEHANDLE);
    }
}

if (scalar(@{$parameter{'selectPrograms'}}) > 0) {
	
	if ( ( any {$_ eq "vcfTools"} @{$parameter{'selectPrograms'}} ) ) { #If element is part of array

	    &VcfTools(\%parameter, $BASHFILEHANDLE);
	}
	if ( ( any {$_ eq "plink"} @{$parameter{'selectPrograms'}} ) ) { #If element is part of array
	    
	    &Plink(\%parameter, $BASHFILEHANDLE);
	}
	if ( ( any {$_ eq "variantEffectPredictor"} @{$parameter{'selectPrograms'}} ) ) { #If element is part of array

	    &VariantEffectPredictor(\%parameter, $BASHFILEHANDLE);
	}
	if ( ( any {$_ eq "CNVnator"} @{$parameter{'selectPrograms'}} ) ) { #If element is part of array
	    
	    &CNVnator(\%parameter, $BASHFILEHANDLE);
	}
}
else {

    &VcfTools(\%parameter, $BASHFILEHANDLE);
    
    &Plink(\%parameter, $BASHFILEHANDLE);
    
    &VariantEffectPredictor(\%parameter, $BASHFILEHANDLE);

    &CNVnator(\%parameter, $BASHFILEHANDLE);
}

close($BASHFILEHANDLE);
#close($LOGFILEHANDLE);

###SubRoutines###

sub CreateBashFile {

    my $fileName = $_[0];

    my $FILEHANDLE = IO::Handle->new();  #Create anonymous filehandle
    my $pwd = cwd();

    ## Open batch file
    open ($FILEHANDLE, ">",$pwd."/".$fileName) or die("Can't write to '".$pwd."/".$fileName."' :".$!."\n");

    print $FILEHANDLE "#!/usr/bin/env bash", "\n\n";
 
    print STDOUT "Will write install instructions to '".$pwd."/".$fileName, "'\n";

   return $FILEHANDLE;
}

sub OpenLogFile {

    my $fileName = $_[0];

    my $FILEHANDLE = IO::Handle->new();  #Create anonymous filehandle
    
    ## Open batch file
    open ($FILEHANDLE, ">",$fileName) or die("Can't write to '".$fileName."' :".$!."\n");

    return $FILEHANDLE;
}

sub PrintParameters {
    
    my $parameterHashRef = $_[0];
    
    foreach my $key (keys %{$parameterHashRef}) {
	
	if (ref(${$parameterHashRef}{$key})!~/ARRAY|HASH/) {

	    print STDOUT $key." ".${$parameterHashRef}{$key}, "\n";
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


sub CreateConda {

    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];

    my $program = "conda";

    if(can_run($program)) {  #IPC::Cmd
	
	print STDERR "ProgramCheck: ".$program." installed", "\n";
    }
    else {
	
	print STDERR "Could not detect ".$program." in your PATH\n";
	exit 1;
    }

    ## Check Conda path
    if (! -d $parameter{'condaPath'}) {

	print STDERR "Could not find miniconda directory in: ".$parameter{'condaPath'}, "\n";
	exit 1;
    }

    print STDERR "Writting install instructions for Conda packages", "\n";

    ## Update Conda
    print $FILEHANDLE "### Update Conda\n";
    print $FILEHANDLE "conda update -y conda ";
    print $FILEHANDLE "\n\n";
}

sub CreateCondaEnvironment {

    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];

    ## Check Conda environment
    if (! -d $parameter{'condaPath'}."/envs/".$parameter{'condaEnvironment'}) {

	## Create conda environment
	print $FILEHANDLE "### Creating Conda Environment and install: ".${$parameterHashRef}{'condaEnvironment'}."\n";
	print $FILEHANDLE "conda create -n ".${$parameterHashRef}{'condaEnvironment'}." ";
	print $FILEHANDLE "-y ";
	print $FILEHANDLE "pip ";
	print $FILEHANDLE "\n\n";
    }
    
    ## Install into conda environment
    print $FILEHANDLE "### Installing into Conda Environment: ".${$parameterHashRef}{'condaEnvironment'}."\n";
    print $FILEHANDLE "conda install ";	
    print $FILEHANDLE "-n ".${$parameterHashRef}{'condaEnvironment'}." ";
    print $FILEHANDLE "-y ";
    print $FILEHANDLE "-c bioconda ";
    
    ## Install all bioConda packages
    foreach my $program (keys %{${$parameterHashRef}{'bioConda'}}) {
	    
	print $FILEHANDLE $program."=".${$parameterHashRef}{'bioConda'}{$program}." ";
    }

    print $FILEHANDLE "\n\n";

    ## Custom 
    foreach my $program (keys %{${$parameterHashRef}{'bioConda'}}) {


	if ($program eq "sambamba") {

	    &AddSoftLink({'parameterHashRef' => $parameterHashRef,
			  'FILEHANDLE' => $BASHFILEHANDLE,
			  'binary' => "sambamba",
			  'softLink' => "sambamba_v".${$parameterHashRef}{'bioConda'}{'sambamba'},
			 });
	}
	if ($program eq "picard") {

	    &AddSoftLink({'parameterHashRef' => $parameterHashRef,
			  'FILEHANDLE' => $BASHFILEHANDLE,
			  'binary' => q?../share/picard-?.${$parameterHashRef}{'bioConda'}{'picard'}.${$parameterHashRef}{'bioCondaPicardPatch'}.q?/picard.jar?,
			  'softLink' => "picard.jar",
			 });
	}
	if ($program eq "snpeff") {
	    
	    &AddSoftLink({'parameterHashRef' => $parameterHashRef,
			  'FILEHANDLE' => $BASHFILEHANDLE,
			  'binary' => q?../share/snpeff-?.${$parameterHashRef}{'bioConda'}{'snpeff'}.${$parameterHashRef}{'bioCondaSnpeffPatch'}.q?/snpEff.jar?,
			  'softLink' => "snpEff.jar",
			 });
	    
	    &AddSoftLink({'parameterHashRef' => $parameterHashRef,
			  'FILEHANDLE' => $BASHFILEHANDLE,
			  'binary' => q?../share/snpeff-?.${$parameterHashRef}{'bioConda'}{'snpeff'}.${$parameterHashRef}{'bioCondaSnpeffPatch'}.q?/SnpSift.jar?,
			  'softLink' => "SnpSift.jar",
			 });
	    &AddSoftLink({'parameterHashRef' => $parameterHashRef,
			  'FILEHANDLE' => $BASHFILEHANDLE,
			  'binary' => q?../share/snpeff-?.${$parameterHashRef}{'bioConda'}{'snpeff'}.${$parameterHashRef}{'bioCondaSnpeffPatch'}.q?/snpEff.config?,
			  'softLink' => "snpEff.config",
			 });
	    &CheckMTCodonTable({'parameterHashRef' => $parameterHashRef,
				'FILEHANDLE' => $BASHFILEHANDLE,
				'shareDirectory' => $parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/share/snpeff-?.${$parameterHashRef}{'bioConda'}{'snpeff'}.${$parameterHashRef}{'bioCondaSnpeffPatch'}.q?/?,
				'binary' => "snpEff.config",
			       });
	}
    }
}


sub Perl {
    
    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];
    
    my $pwd = cwd();
    
    if ($ENV{PATH}=~/perl-${$parameterHashRef}{'perl'}/) {
	
	if (${$parameterHashRef}{'update'} == 0) {
	    
	    print STDERR "Found perl-".${$parameterHashRef}{'perl'}.". in your path\n";
	    print STDERR q?Skipping writting installation for perl-?.${$parameterHashRef}{'perl'},"\n";  
	}
	else {

	    if (${$parameterHashRef}{'perlInstall'} == 1) {
	    
		## Removing specific Perl version
		print $FILEHANDLE "### Removing specific Perl version\n";
		print $FILEHANDLE q?rm -rf $HOME/perl-?.${$parameterHashRef}{'perl'};
		print $FILEHANDLE "\n\n";
		
		&InstallPerlCpnam($parameterHashRef, $BASHFILEHANDLE); 
	    }
	    
	    &PerlModules($parameterHashRef, $BASHFILEHANDLE,);
	}
    }
    else {
    	
	if (${$parameterHashRef}{'perlInstall'} == 1) {
	    
	    &InstallPerlCpnam($parameterHashRef, $BASHFILEHANDLE, "AddPath");
	}

	&PerlModules($parameterHashRef, $BASHFILEHANDLE,);
    }
}


sub InstallPerlCpnam {

    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];
    my $path = $_[2];
    
    my $pwd = cwd();

    print STDERR "Writting install instructions for Perl and Cpanm", "\n";
    
    ## Install specific Perl version
    print $FILEHANDLE "### Install specific Perl version\n";
    
    ## Move to Home
    print $FILEHANDLE "## Move HOME\n";
    print $FILEHANDLE q?cd $HOME?;
    print $FILEHANDLE "\n\n";
    
    ## Download
    print $FILEHANDLE "## Download Perl\n";
    print $FILEHANDLE "wget --quiet http://www.cpan.org/src/5.0/perl-".${$parameterHashRef}{'perl'}.".tar.gz ";
    print $FILEHANDLE "-O perl-".${$parameterHashRef}{'perl'}.".tar.gz";  #Dowload outfile
    print $FILEHANDLE "\n\n";
    
    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "tar xzf perl-".${$parameterHashRef}{'perl'}.".tar.gz";
    print $FILEHANDLE "\n\n";
    
    ## Move to perl directory
    print $FILEHANDLE "## Move to perl directory\n";
    print $FILEHANDLE "cd perl-".${$parameterHashRef}{'perl'};
    print $FILEHANDLE "\n\n";
    
    ## Configure
    print $FILEHANDLE "## Configure\n";
    print $FILEHANDLE q?./Configure -des -Dprefix=$HOME/perl-?.${$parameterHashRef}{'perl'};
    print $FILEHANDLE "\n";
    
    print $FILEHANDLE "make";
    print $FILEHANDLE "\n";
    
    print $FILEHANDLE "make test";
    print $FILEHANDLE "\n";
    
    print $FILEHANDLE "make install";
    print $FILEHANDLE "\n\n";
    
    if ($path) {
	
	## Export path
	print $FILEHANDLE "## Export path\n";
	print $FILEHANDLE q?echo 'export PATH=$HOME/perl-?.${$parameterHashRef}{'perl'}.q?/:$PATH' >> ~/.bashrc?;
	print $FILEHANDLE "\n\n";
	print $FILEHANDLE q?export PATH=$HOME/perl-?.${$parameterHashRef}{'perl'}.q?/:$PATH?;  #Use newly installed perl
	print $FILEHANDLE "\n\n";
    }

    ## Remove tar file
    print $FILEHANDLE "## Remove tar file\n";
    print $FILEHANDLE "cd && rm perl-".${$parameterHashRef}{'perl'}.".tar.gz";
    print $FILEHANDLE "\n\n";
    
    ## Move to back
    print $FILEHANDLE "## Move to original working directory\n";
    print $FILEHANDLE "cd ".$pwd;
    print $FILEHANDLE "\n\n";

    if ($path) {

	print $FILEHANDLE q?echo 'eval `perl -I ~/perl-?.${$parameterHashRef}{'perl'}.q?/lib/perl5/ -Mlocal::lib=~/perl-?.${$parameterHashRef}{'perl'}.q?/`' >> ~/.bash_profile ?;  #Add at start-up
	print $FILEHANDLE "\n\n";
    }

    ## Install Perl modules via cpanm
    print $FILEHANDLE "## Install cpanm\n";
    print $FILEHANDLE q?wget -O- http://cpanmin.us | perl - -l $HOME/perl-?.${$parameterHashRef}{'perl'}.q?/bin App::cpanminus --local-lib=~/perl-?.${$parameterHashRef}{'perl'}.q?/ local::lib ?;
    print $FILEHANDLE "\n\n";

    ## Use newly installed perl
    print $FILEHANDLE q?eval `perl -I ~/perl-?.${$parameterHashRef}{'perl'}.q?/lib/perl5/ -Mlocal::lib=~/perl-?.${$parameterHashRef}{'perl'}.q?/` ?;
    print $FILEHANDLE "\n\n";

    ## Use newly installed perl
    print $FILEHANDLE q?PERL5LIB=~/perl-?.${$parameterHashRef}{'perl'}.q?/lib/perl5?;
    print $FILEHANDLE "\n\n";
}
    

sub PerlModules {
    
    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];
    
    ## Install Perl modules via cpanm
    print $FILEHANDLE "## Install Perl modules via cpanm\n";
    print $FILEHANDLE "cpanm ";
    print $FILEHANDLE join(" ", @{${$parameterHashRef}{'perlModules'}})." ";
    print $FILEHANDLE "\n\n";
}


sub PipInstall {

    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];

    print STDERR "Writting install instructions for pip packages", "\n";

    ## Install PIP packages in conda environment
    print $FILEHANDLE "### Install PIP packages in conda environment: ".${$parameterHashRef}{'condaEnvironment'}."\n";
    &ActivateCondaEnvironment($parameterHashRef, $FILEHANDLE);

    ## Install PIP packages
    print $FILEHANDLE "## Install PIP packages\n";
    print $FILEHANDLE "pip install ";

    ## Install all PIP packages
    foreach my $program (keys %{${$parameterHashRef}{'pip'}}) {

	print $FILEHANDLE $program."==".${$parameterHashRef}{'pip'}{$program}." ";
    }
    print $FILEHANDLE "\n\n";

    &DeactivateCondaEnvironment($FILEHANDLE);
}


sub PicardTools {

    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];

    my $pwd = cwd();

    if (&CheckCondaBinFileExists($parameterHashRef, "picard.jar")) {  # Assumes that picard.jar is there as well then

	return
    }

    ## Install picard
    print $FILEHANDLE "### Install Picard\n\n";

    &CreateInstallDirectory($FILEHANDLE);
    
    ## Download
    print $FILEHANDLE "## Download Picard\n";
    print $FILEHANDLE "wget --quiet https://github.com/broadinstitute/picard/releases/download/".${$parameterHashRef}{'picardTools'}."/picard-tools-".${$parameterHashRef}{'picardTools'}.".zip ";
    print $FILEHANDLE "-O picard-tools-".${$parameterHashRef}{'picardTools'}.".zip";  #Download outfile
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "unzip picard-tools-".${$parameterHashRef}{'picardTools'}.".zip";
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    if (-d $parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/share/picard-tools-?.${$parameterHashRef}{'picardTools'}) {

	print $FILEHANDLE "rm -rf ".$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/share/picard-tools-?.${$parameterHashRef}{'picardTools'};
	print $FILEHANDLE "\n\n";
    }

    print $FILEHANDLE "## Make available from conda environment\n";
    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?picard-tools-?.${$parameterHashRef}{'picardTools'}.q? ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/share/?;
    print $FILEHANDLE "\n\n";

    &AddSoftLink({'parameterHashRef' => $parameterHashRef,
		  'FILEHANDLE' => $BASHFILEHANDLE,
		  'binary' => q?../share/picard-tools-?.${$parameterHashRef}{'picardTools'}.q?/picard.jar?,
		  'softLink' => "picard.jar",
		 });

    &CleanUpModuleInstall($FILEHANDLE, $pwd);
}


sub Sambamba {

    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];

    my $pwd = cwd();

    if (&CheckCondaBinFileExists($parameterHashRef, "sambamba", ${$parameterHashRef}{'sambamba'}) ) {

	return 
    }
    ## Install sambamba
    print $FILEHANDLE "### Install sambamba\n\n";

    &CreateInstallDirectory($FILEHANDLE);

    ## Download
    print $FILEHANDLE "## Download sambamba release\n";
    print $FILEHANDLE q?wget --quiet https://github.com/lomereiter/sambamba/releases/download/v?.${$parameterHashRef}{'sambamba'}.q?/sambamba_v?.${$parameterHashRef}{'sambamba'}.q?_linux.tar.bz2 ?;
    print $FILEHANDLE "-O sambamba_v".${$parameterHashRef}{'sambamba'}."_linux.tar.bz2";  #Download outfile
    print $FILEHANDLE "\n\n";

    ## Decompress
    print $FILEHANDLE "## Decompress sambamba file\n";
    print $FILEHANDLE "bzip2 ";
    print $FILEHANDLE "-f ";  #Force
    print $FILEHANDLE "-d ";  #Decompress
    print $FILEHANDLE "sambamba_v".${$parameterHashRef}{'sambamba'}."_linux.tar.bz2";
    print $FILEHANDLE "\n\n";

    ## Extract files
    print $FILEHANDLE "## Extract files\n";
    print $FILEHANDLE "tar xvf sambamba_v".${$parameterHashRef}{'sambamba'}."_linux.tar";
    print $FILEHANDLE "\n\n";

    ## Make executable
    print $FILEHANDLE "## Make executable\n";
    print $FILEHANDLE "chmod 755 ";
    print $FILEHANDLE "sambamba_v".${$parameterHashRef}{'sambamba'};
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    print $FILEHANDLE "## Make available from conda environment\n";
    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?sambamba_v?.${$parameterHashRef}{'sambamba'}.q? ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/?;
    print $FILEHANDLE "\n\n";

    &AddSoftLink({'parameterHashRef' => $parameterHashRef,
		  'FILEHANDLE' => $BASHFILEHANDLE,
		  'binary' => "sambamba_v".${$parameterHashRef}{'bioConda'}{'sambamba'},
		  'softLink' => "sambamba",
		 });

    &CleanUpModuleInstall($FILEHANDLE, $pwd);

}


sub VcfTools {

    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];

    my $pwd = cwd();

    if(&CheckCondaBinFileExists($parameterHashRef, "vcftools")) {
    
	return
    }

    ## Install vcfTools
    print $FILEHANDLE "### Install vcfTools\n\n";

    &CreateInstallDirectory($FILEHANDLE);
    
    ## Download
    print $FILEHANDLE "## Download vcfTools\n";
    print $FILEHANDLE "wget --quiet https://github.com/vcftools/vcftools/releases/download/v".${$parameterHashRef}{'vcfTools'}."/vcftools-".${$parameterHashRef}{'vcfTools'}.".tar.gz ";
    print $FILEHANDLE "-O vcftools-".${$parameterHashRef}{'vcfTools'}.".tar.gz";  #Download outfile
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "tar xvf vcftools-".${$parameterHashRef}{'vcfTools'}.".tar.gz";
    print $FILEHANDLE "\n\n";

    ## Export PERL5LIB environment variable
    print $FILEHANDLE "## Export PERL5LIB environment variable\n";
    print $FILEHANDLE q?export PERL5LIB=?.$Bin.q?/vcftools-?.${$parameterHashRef}{'vcfTools'}.q?/src/perl/?;
    print $FILEHANDLE "\n\n";

    ## Move to vcfTools directory
    print $FILEHANDLE "## Move to vcfTools directory\n";
    print $FILEHANDLE "cd vcftools-".${$parameterHashRef}{'vcfTools'};
    print $FILEHANDLE "\n\n";

    ## Configure
    my $filePath = $parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'};

    print $FILEHANDLE "## Configure\n";
    print $FILEHANDLE q?./configure --prefix=?.$filePath;
    print $FILEHANDLE "\n";

    print $FILEHANDLE "make";
    print $FILEHANDLE "\n";

    print $FILEHANDLE "make install";
    print $FILEHANDLE "\n\n";

    ## Move Perl Module
    print $FILEHANDLE "## Move Perl Module\n";
    print $FILEHANDLE q?cp src/perl/Vcf.pm $HOME/perl-?.${$parameterHashRef}{'perl'}.q?/lib/perl5/?;
    print $FILEHANDLE "\n\n";

    &CleanUpModuleInstall($FILEHANDLE, $pwd);

    ## Reset perl envionment
    print $FILEHANDLE q?PERL5LIB=~/perl-?.${$parameterHashRef}{'perl'}.q?/lib/perl5?;
    print $FILEHANDLE "\n\n";
}


sub BedTools {
    
    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];
    
    my $pwd = cwd();
    
    if(&CheckCondaBinFileExists($parameterHashRef, "bedtools")) {
	
	return
    }

    my $bedToolsMainVersion = substr(${$parameterHashRef}{'bedTools'}, 0, 1);

    ## Install bedTools
    print $FILEHANDLE "### Install bedTools\n\n";
    
    &CreateInstallDirectory($FILEHANDLE);
    
    ## Download
    print $FILEHANDLE "## Download bedTools\n";
    print $FILEHANDLE "wget --quiet https://github.com/arq5x/bedtools".$bedToolsMainVersion."/releases/download/v".${$parameterHashRef}{'bedTools'}."/bedtools-".${$parameterHashRef}{'bedTools'}.".tar.gz ";
    print $FILEHANDLE "-O bedtools-".${$parameterHashRef}{'bedTools'}.".tar.gz";  #Download outfile
    print $FILEHANDLE "\n\n";
    
    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "tar xvf bedtools-".${$parameterHashRef}{'bedTools'}.".tar.gz";
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
    print $FILEHANDLE q?./bin/* ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/?;
    print $FILEHANDLE "\n\n";
    
    &CleanUpModuleInstall($FILEHANDLE, $pwd);
}


sub VT {

    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];

    my $pwd = cwd();

    if (&CheckCondaBinFileExists($parameterHashRef, "vt")) {
    
	return
    }

    ## Install VT
    print $FILEHANDLE "### Install VT\n\n";

    &CreateInstallDirectory($FILEHANDLE);
    
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
    print $FILEHANDLE "make";
    print $FILEHANDLE "\n";

    print $FILEHANDLE "make test";
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    print $FILEHANDLE "## Make available from conda environment\n";
    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?vt ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/?;
    print $FILEHANDLE "\n\n";

    &CleanUpModuleInstall($FILEHANDLE, $pwd);
}


sub Plink {

    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];

    my $pwd = cwd();

    if (&CheckCondaBinFileExists($parameterHashRef, "plink")) {

	return
    }

    ## Install Plink
    print $FILEHANDLE "### Install Plink\n\n";

    &CreateInstallDirectory($FILEHANDLE);
    
    ## Download
    print $FILEHANDLE "## Download Plink\n";
    print $FILEHANDLE "wget --quiet https://www.cog-genomics.org/static/bin/plink".${$parameterHashRef}{'plink'}."/plink_linux_x86_64.zip ";
    print $FILEHANDLE "-O plink-".${$parameterHashRef}{'plink'}."-x86_64.zip";  #Download outfile
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "unzip plink-".${$parameterHashRef}{'plink'}."-x86_64.zip";
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    print $FILEHANDLE "## Make available from conda environment\n";
    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?plink ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/?;
    print $FILEHANDLE "\n\n";

    &CleanUpModuleInstall($FILEHANDLE, $pwd);
}


sub SnpEff {

    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];

    my $pwd = cwd();

    if (&CheckCondaBinFileExists($parameterHashRef, "snpEff.jar")) {  # Assumes that SnpSift.jar is there as well then

	return
    }

    ## Install SnpEff
    print $FILEHANDLE "### Install SnpEff\n\n";

    &CreateInstallDirectory($FILEHANDLE);
    
    ## Download
    print $FILEHANDLE "## Download SnpEff\n";
    print $FILEHANDLE "wget --quiet http://sourceforge.net/projects/snpeff/files/snpEff_".${$parameterHashRef}{'snpEff'}."_core.zip/download ";
    print $FILEHANDLE "-O snpEff_".${$parameterHashRef}{'snpEff'}."_core.zip";  #Download outfile
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "unzip snpEff_".${$parameterHashRef}{'snpEff'}."_core.zip";
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    if (-d $parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/share/snpEff.?.${$parameterHashRef}{'snpEff'}) {
	
	print $FILEHANDLE "rm -rf ".$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/share/snpEff.?.${$parameterHashRef}{'snpEff'};
	print $FILEHANDLE "\n\n";
    }

    print $FILEHANDLE "## Make available from conda environment\n";
    print $FILEHANDLE "mkdir -p ".$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/share/snpEff.?.${$parameterHashRef}{'snpEff'};
    print $FILEHANDLE "\n\n";

    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?snpEff/*.jar ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/share/snpEff.?.${$parameterHashRef}{'snpEff'}.q?/?;
    print $FILEHANDLE "\n\n";

    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?snpEff/snpEff.config ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/share/snpEff.?.${$parameterHashRef}{'snpEff'}.q?/?;
    print $FILEHANDLE "\n\n";

    &AddSoftLink({'parameterHashRef' => $parameterHashRef,
		  'FILEHANDLE' => $BASHFILEHANDLE,
		  'binary' => q?../share/snpEff.?.${$parameterHashRef}{'snpEff'}.q?/snpEff.jar?,
		  'softLink' => "snpEff.jar",
		 });
    
    &AddSoftLink({'parameterHashRef' => $parameterHashRef,
		  'FILEHANDLE' => $BASHFILEHANDLE,
		  'binary' => q?../share/snpEff.?.${$parameterHashRef}{'snpEff'}.q?/SnpSift.jar?,
		  'softLink' => "SnpSift.jar",
		 });

    &AddSoftLink({'parameterHashRef' => $parameterHashRef,
		  'FILEHANDLE' => $BASHFILEHANDLE,
		  'binary' => q?../share/snpEff.?.${$parameterHashRef}{'snpEff'}.q?/snpEff.config?,
		  'softLink' => "snpEff.config",
		 });
    &CheckMTCodonTable({'parameterHashRef' => $parameterHashRef,
			'FILEHANDLE' => $BASHFILEHANDLE,
			'shareDirectory' => $parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/share/snpEff.?.${$parameterHashRef}{'snpEff'}.q?/?,
			'binary' => "snpEff.config",
		       });
    &CleanUpModuleInstall($FILEHANDLE, $pwd);
}


sub VariantEffectPredictor {

    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];

    my $pwd = cwd();
    
    my $minicondaBinDirectory = $parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/ensembl-tools-release-?.${$parameterHashRef}{'variantEffectPredictor'};

    if (-d $minicondaBinDirectory) {

	print STDERR q?Found VariantEffectPredictor in miniconda directory: ?.$minicondaBinDirectory, "\n";
	
	if (${$parameterHashRef}{'update'} == 0) {

	    print STDERR "Skipping writting installation process for VariantEffectPredictor","\n";  	    
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
	
	print STDERR "Writting install instructions for VariantEffectPredictor", "\n";
    }

    ## Install VEP
    print $FILEHANDLE "### Install VariantEffectPredictor\n\n";

    &ActivateCondaEnvironment($parameterHashRef, $FILEHANDLE);

    ##Make sure that the cache directory exists
    print $FILEHANDLE "mkdir -p ".${$parameterHashRef}{'vepDirectoryCache'}." ";  #Cache directory
    print $FILEHANDLE "\n\n";

    ## Move to miniconda environment
    print $FILEHANDLE q?cd ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'};
    print $FILEHANDLE "\n\n";

    ## Download
    print $FILEHANDLE "## Download VEP\n";
    print $FILEHANDLE "wget --quiet https://github.com/Ensembl/ensembl-tools/archive/release/".${$parameterHashRef}{'variantEffectPredictor'}.".zip ";
    print $FILEHANDLE "-O VariantEffectPredictor-".${$parameterHashRef}{'variantEffectPredictor'}.".zip";  #Download outfile
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "unzip VariantEffectPredictor-".${$parameterHashRef}{'variantEffectPredictor'}.".zip";
    print $FILEHANDLE "\n\n";    

    ## Move to VariantEffectPredictor directory
    print $FILEHANDLE "## Move to VariantEffectPredictor directory\n";
    print $FILEHANDLE "cd ensembl-tools-release-".${$parameterHashRef}{'variantEffectPredictor'}."/scripts/variant_effect_predictor/";
    print $FILEHANDLE "\n\n";

    ## Install VEP
    print $FILEHANDLE "## Install VEP\n";
    print $FILEHANDLE "perl INSTALL.pl ";
    print $FILEHANDLE "--AUTO alcfp ";  #a (API), l (FAIDX/htslib), c (cache), f (FASTA), p (plugins)
    print $FILEHANDLE "-g ".$parameter{'variantEffectPredictorPlugin'}." ";  #Plugins 
    print $FILEHANDLE "-c ".${$parameterHashRef}{'vepDirectoryCache'}." ";  #Cache directory
    print $FILEHANDLE "-s homo_sapiens ";
    print $FILEHANDLE "--ASSEMBLY GRCh37 ";
    print $FILEHANDLE "\n\n";

    ## Clean up
    print $FILEHANDLE "## Clean up\n";
    print $FILEHANDLE q?cd ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'};
    print $FILEHANDLE "\n\n";

    print $FILEHANDLE "rm -rf VariantEffectPredictor-".${$parameterHashRef}{'variantEffectPredictor'}.".zip";;
    print $FILEHANDLE "\n\n";

    ## Moving up
    print $FILEHANDLE "## Moving back to original working directory\n";
    print $FILEHANDLE "cd ".$pwd;  #Go back to subroutine origin
    print $FILEHANDLE "\n\n";

    &DeactivateCondaEnvironment($FILEHANDLE);
}


sub CNVnator {

    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];

    my $pwd = cwd();

    if (&CheckCondaBinFileExists($parameterHashRef, "cnvnator")) {

	return
    }

    my $minicondaBinDirectory = $parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/root?;

    if (-d $minicondaBinDirectory) {

	print STDERR q?Found Root in miniconda directory: ?.$minicondaBinDirectory, "\n";
	
	if (${$parameterHashRef}{'update'} == 0) {

	    print STDERR "Skipping writting installation process for Root","\n";  	    
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
	
	print STDERR "Writting install instructions for Root", "\n";
    }

    ## Install Root
    print $FILEHANDLE "### Install CNVnator/Root\n\n";

    ## Move to miniconda environment
    print $FILEHANDLE q?cd ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'};
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

    unless ($ENV{PATH}=~/$parameter{'condaPath'}\/envs\/${$parameterHashRef}{'condaEnvironment'}\/root\/bin/) {
	
	## Export path
	print $FILEHANDLE "## Export path\n";
	print $FILEHANDLE q?echo 'source ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/root/bin/thisroot.sh' >> ~/.bashrc?;
	print $FILEHANDLE "\n\n";

	## Use newly installed root
	print $FILEHANDLE q?source ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/root/bin/thisroot.sh ?;
	print $FILEHANDLE "\n\n";
    }
    
    ## Moving up
    print $FILEHANDLE "## Moving back to original working directory\n";
    print $FILEHANDLE "cd ".$pwd;  #Go back to subroutine origin
    print $FILEHANDLE "\n\n";

    ## Install CNVNator
    print $FILEHANDLE "### Install CNVnator\n\n";

    &CreateInstallDirectory($FILEHANDLE);
    
    ## Download
    print $FILEHANDLE "## Download CNVNator\n";
    print $FILEHANDLE "wget --quiet https://github.com/abyzovlab/CNVnator/releases/download/v".${$parameterHashRef}{'CNVnator'}."/CNVnator_v".${$parameterHashRef}{'CNVnator'}.".zip ";
    print $FILEHANDLE "-O CNVnator_v".${$parameterHashRef}{'CNVnator'}.".zip";  #Download outfile
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "unzip CNVnator_v".${$parameterHashRef}{'CNVnator'}.".zip";
    print $FILEHANDLE "\n\n";

    ## Move to CNVnator directory
    print $FILEHANDLE "## Move to CNVnator directory\n";
    print $FILEHANDLE "cd CNVnator_v".${$parameterHashRef}{'CNVnator'}."/src/samtools/";
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
    print $FILEHANDLE q?cnvnator ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/?;
    print $FILEHANDLE "\n\n";

    &CleanUpModuleInstall($FILEHANDLE, $pwd);
}


sub FindTranslocations {

    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];

    my $pwd = cwd();

    if (&CheckCondaBinFileExists($parameterHashRef, "FindTranslocations")) {

	return
    }

    ## Install Plink
    print $FILEHANDLE "### Install FindTranslocations\n\n";

    &CreateInstallDirectory($FILEHANDLE);
    
    ## Download
    print $FILEHANDLE "## Download FindTranslocations\n";
    print $FILEHANDLE "wget --quiet https://www.cog-genomics.org/static/bin/plink".${$parameterHashRef}{'FindTranslocations'}."/plink_linux_x86_64.zip ";
    print $FILEHANDLE "-O FindTranslocations-".${$parameterHashRef}{'FindTranslocations'}."-x86_64.zip";  #Download outfile
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "unzip FindTranslocations-".${$parameterHashRef}{'plink'}."-x86_64.zip";
    print $FILEHANDLE "\n\n";

    ## Move to FindTranslocations directory
    print $FILEHANDLE "## Move to FindTranslocations directory\n";
    print $FILEHANDLE "cd FindTranslocations".${$parameterHashRef}{'CNVnator'};
    print $FILEHANDLE "\n";
    print $FILEHANDLE "mkdir -p build";
    print $FILEHANDLE "\n";
    print $FILEHANDLE "cd build";
    print $FILEHANDLE "\n";

    ## Configure
    print $FILEHANDLE "## Configure\n";
    print $FILEHANDLE "cmake .. -DBoost_NO_BOOST_CMAKE=ON";
    print $FILEHANDLE "\n\n";

    print $FILEHANDLE "make";
    print $FILEHANDLE "\n\n";

    print $FILEHANDLE "cd ../bin";
    print $FILEHANDLE "\n";
    print $FILEHANDLE "chmod a+x FindTranslocations";
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    print $FILEHANDLE "## Make available from conda environment\n";
    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?FindTranslocations ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/?;
    print $FILEHANDLE "\n\n";

    &CleanUpModuleInstall($FILEHANDLE, $pwd);
}



sub ActivateCondaEnvironment {

    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];

    ## Activate conda environment and install cpanm and MIP modules
    print $FILEHANDLE "## Activate conda environment\n";
    print $FILEHANDLE "source activate ".${$parameterHashRef}{'condaEnvironment'}." ";
    print $FILEHANDLE "\n\n";
}


sub DeactivateCondaEnvironment {

    my $FILEHANDLE = $_[0];

    ## Deactivate conda environment
    print $FILEHANDLE "## Deactivate conda environment\n";
    print $FILEHANDLE "source deactivate ";
    print $FILEHANDLE "\n\n";
}


sub CleanUpModuleInstall {

    my $FILEHANDLE = $_[0];
    my $pwd = $_[1];

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

    my $FILEHANDLE = $_[0];

    ## Create temp install directory
    print $FILEHANDLE "## Create temp install directory\n";
    print $FILEHANDLE "mkdir -p .MIP ";
    print $FILEHANDLE "\n";
    print $FILEHANDLE "cd .MIP";
    print $FILEHANDLE "\n\n";
}


sub CheckCondaBinFileExists {
    
    my $parameterHashRef = $_[0];
    my $programName = $_[1];
    my $programVersion = $_[2];
    
    my $minicondaBinFile = $parameter{'condaPath'};

    if ($programName eq "sambamba") {
	
	$programVersion = "_v".${$parameterHashRef}{'sambamba'};
    }
    
    if ($programVersion) {

	$minicondaBinFile .= q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/?.$programName.$programVersion;
    }
    else {

	$minicondaBinFile .= q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/?.$programName;
    }

    if (-f $minicondaBinFile) {

	if ($programVersion) {
	    
	    print STDERR q?Found ?.$programName.q? version ?.$programVersion.q? in miniconda directory: ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/?, "\n";
	    
	    if (${$parameterHashRef}{'update'} == 0) {

		print STDERR q?Skipping writting installation process for ?.$programName.q? ?.$programVersion,"\n";  
		return 1;
	    }
	    print STDERR "Writting install instructions for ".$programName, "\n";
	}   
	else {
	    
	    print STDERR q?Found ?.$programName.q? in miniconda directory: ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/?, "\n";
	    
	    if (${$parameterHashRef}{'update'} == 0) {
		
		print STDERR q?Skipping writting installation process for ?.$programName,"\n";  	    
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


sub AddSoftLink {

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $parameterHashRef = ${$argHashRef}{'parameterHashRef'};
    my $FILEHANDLE = ${$argHashRef}{'FILEHANDLE'};
    my $binary = ${$argHashRef}{'binary'};
    my $softLink = ${$argHashRef}{'softLink'} ;
    
    my $pwd = cwd();

    ## Add softlink
    print $FILEHANDLE "## Add softlink\n";
    print $FILEHANDLE "cd ".$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/?;
    print $FILEHANDLE "\n";

    print $FILEHANDLE "ln -f -s ";
    print $FILEHANDLE $binary.q? ?.$softLink;
    print $FILEHANDLE "\n\n";

    ## Move to back
    print $FILEHANDLE "## Move to original working directory\n";
    print $FILEHANDLE "cd ".$pwd;
    print $FILEHANDLE "\n\n";
}


sub CheckMTCodonTable {

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $parameterHashRef = ${$argHashRef}{'parameterHashRef'};
    my $FILEHANDLE = ${$argHashRef}{'FILEHANDLE'};
    my $shareDirectory = ${$argHashRef}{'shareDirectory'};
    my $binary = ${$argHashRef}{'binary'};
    
    my $pwd = cwd();

    my $detectRegExp = q?perl -nae 'if($_=~/GRCh37.75.MT.codonTable/) {print 1}' ?;
    my $addRegExp = q?perl -nae 'if($_=~/GRCh37.75.reference/) {print $_; print "GRCh37.75.MT.codonTable : Vertebrate_Mitochondrial\n"} else {print $_;}' ?;
    my $ret;

    if (-f $shareDirectory."/".$binary) {

	$ret = `$detectRegExp $shareDirectory/$binary`;
    }
    if (!$ret) {  #No MT.codonTable in config

	print $FILEHANDLE "## Adding GRCh37.75.MT.codonTable : Vertebrate_Mitochondrial to ".$shareDirectory.$binary;
	print $FILEHANDLE "\n";

	## Add MT.codon Table to config
	print $FILEHANDLE $addRegExp." ".$shareDirectory.$binary." > ".$shareDirectory.$binary.".tmp";
	print $FILEHANDLE "\n";
	print $FILEHANDLE "mv ".$shareDirectory.$binary.".tmp ".$shareDirectory.$binary;
	print $FILEHANDLE "\n\n";
	
    }
    else {

	print STDERR  "Found MT.codonTable in ".$shareDirectory."snpEff.config. Skipping addition to snpEff config\n"; 
    }
}
