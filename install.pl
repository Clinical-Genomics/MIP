#!/usr/bin/env perl

##Assumes you have a working conda installation

use strict;
use warnings;

use Getopt::Long;
use IO::File;
use Cwd;
use FindBin qw($Bin); #Find directory of script
use vars qw($USAGE);

BEGIN {
    $USAGE =
	qq{install.pl [options]
           -env/--condaEnvironment Conda environment (Default: "mip")
           -cdp/--condaPath The conda path (Default: "HOME/miniconda")
           -bvc/--bioConda Set the module version of the programs that can be installed with bioConda (e.g. 'bwa=0.7.12')
           -per/--perl Set the perl version (defaults: "5.18.2")
           -pm/perlModules Set the perl modules to be installed via cpanm (comma sep)
           -pip/--pip Set the module version of the programs that can be installed with pip (e.g. 'genmod=3.3.3')
           -sbb/sambamba Set the sambamba version (Default: "0.5.9")
           -vct/--vcfTools Set the vcftools version (Default: "0.1.14")
           -bet/--bedTools Set the bedtools version (Default: "2.25.0")
           -vt/--vt Set the vt version (Default: "0.57")
           -plk/--plink  Set the plink version (Default: "151117")
           -vep/--variantEffectPredictor Set the VEP version (Default: "83")
	   -vepc/--vepDirectoryCache Specify the cache directory to use (whole path; defaults to "~/miniconda/envs/condaEnvironment/ensembl-tools-release-variantEffectPredictorVersion/cache")
           -vepp/--variantEffectPredictorPlugin Supply a comma separated list of VEP plugins (Default: "UpDownDistance")
           -pbc/--preferBioConda Bioconda will used for overlapping shell and biconda installations (Default: 1 (=yes))
           -ppd/--printParameterDefaults Print the parameter defaults
           -u/--update Always install program (Default: 1 (=yes))
           -h/--help Display this help message   
           -v/--version Display version
        };    
}
my ($installDirectory) = (".MIP");

my %parameter; 

### Set parameter default

##Conda
$parameter{'update'} = 1;
$parameter{'preferBioConda'} = 1;
$parameter{'condaEnvironment'} = "mip";
$parameter{'condaPath'} = $ENV{HOME}."/miniconda";

$parameter{'bioConda'}{'bwa'} = "0.7.12";
$parameter{'bioConda'}{'fastqc'} = "0.11.4";
$parameter{'bioConda'}{'samtools'} = "1.2";
$parameter{'bioConda'}{'bcftools'} = "1.2";
#$parameter{'bioConda'}{'snpeff'} = "4.1";
$parameter{'bioConda'}{'picard'} = "1.141";
$parameter{'bioConda'}{'mosaik'} = "2.2.26";
$parameter{'bioConda'}{'htslib'} = "1.2.1";
$parameter{'bioConda'}{'bedtools'} = "2.25.0";
$parameter{'bioConda'}{'vt'} = "2015.11.10";
$parameter{'bioConda'}{'sambamba'} = "0.5.9";

##Perl Modules
$parameter{'perl'} = "5.18.2";
$parameter{'perlModules'} = ["YAML",
			     "Log::Log4perl",
			     "List::MoreUtils",
			     "DateTime",
			     "DateTime::Format::ISO8601",
			     "DateTime::Format::HTTP",
			     "DateTime::Format::Mail",
			     "Set::IntervalTree",  # vcfParser
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
$parameter{'pip'}{'genmod'} = "3.4";
$parameter{'pip'}{'chanjo'} = "3.1.1";
$parameter{'pip'}{'cosmid'} = "0.4.9.1";
$parameter{'pip'}{'python-Levenshtein'} = "0.12.0";

## Programs currently not supported by conda or other packet manager
$parameter{'sambamba'} = "0.5.9";
$parameter{'vcfTools'} = "0.1.14";
$parameter{'bedTools'} = "2.25.0";
$parameter{'vt'} = "gitRepo";
$parameter{'plink'} = "160110";
$parameter{'snpEff'} = "v4_2";
$parameter{'variantEffectPredictor'} = "83";
$parameter{'vepDirectoryCache'} = $parameter{'condaPath'}.q?/envs/?.$parameter{'condaEnvironment'}.q?/ensembl-tools-release-?.$parameter{'variantEffectPredictor'}.q?/cache?;  #Cache directory;
$parameter{'variantEffectPredictorPlugin'} = "UpDownDistance";

my $installVersion = "0.0.2";

###User Options
GetOptions('env|condaEnvironment:s'  => \$parameter{'condaEnvironment'},
	   'cdp|condaPath:s' => \$parameter{'condaPath'},
	   'bcv|bioConda=s'  => \%{$parameter{'bioConda'}},
	   'per|perl=s' => \$parameter{'perl'},
	   'pm|perlModules:s'  => \@{$parameter{'perlModules'}},  #Comma separated list
	   'pip|pip=s'  => \%{$parameter{'pip'}},
	   'sbb|sambamba:s'  => \$parameter{'sambamba'},
	   'vct|vcfTools:s'  => \$parameter{'vcfTools'},
	   'bet|bedTools:s' =>\$parameter{'bedTools'}, 
	   'vt|vt:s'  => \$parameter{'vt'},
	   'plk|plink:s'  => \$parameter{'plink'},
	   'vep|variantEffectPredictor:s'  => \$parameter{'variantEffectPredictor'},
	   'vepc|vepDirectoryCache:s'  => \$parameter{'vepDirectoryCache'},  #path to vep cache dir
	   'vepp|variantEffectPredictorPlugin:s'  => \$parameter{'variantEffectPredictorPlugin'},  #Comma sep string
	   'pbc|preferBioConda:s'  => \$parameter{'preferBioConda'},  # Bioconda will used for overlapping shell and biconda installlations
	   'ppd|printParameterDefaults'  => sub { &PrintParameters(\%parameter); exit;},  #Display parameter defaults
	   'u|update:n' => \$parameter{'update'},
	   'h|help' => sub { print STDOUT $USAGE, "\n"; exit;},  #Display help text
	   'v|version' => sub { print STDOUT "\ninstall.pl ".$installVersion, "\n\n"; exit;},  #Display version number
    );

###MAIN###


#my $LOGFILEHANDLE = &OpenLogFile("MIP_installation.log");

my $BASHFILEHANDLE = &CreateBashFile("mip.sh");

&CreateConda(\%parameter, $BASHFILEHANDLE);

&CreateCondaEnvironment(\%parameter, $BASHFILEHANDLE);

&Perl(\%parameter, $BASHFILEHANDLE);

&PipInstall(\%parameter, $BASHFILEHANDLE);

if ($parameter{'preferBioConda'} != 1) {

    &Sambamba(\%parameter, $BASHFILEHANDLE);

    &BedTools(\%parameter, $BASHFILEHANDLE);

    &VT(\%parameter, $BASHFILEHANDLE);

    &SnpEff(\%parameter, $BASHFILEHANDLE);
}

&VcfTools(\%parameter, $BASHFILEHANDLE);

&Plink(\%parameter, $BASHFILEHANDLE);

&VariantEffectPredictor(\%parameter, $BASHFILEHANDLE);

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

    if (! -f $parameter{'condaPath'}."/envs/".$parameter{'condaEnvironment'}."/bin/sambamba_v".${$parameterHashRef}{'bioConda'}{'sambamba'}) {
	
	&AddSoftLink({'parameterHashRef' => $parameterHashRef,
		      'FILEHANDLE' => $BASHFILEHANDLE,
		      'binary' => "sambamba",
		      'softLink' => "sambamba_v".${$parameterHashRef}{'bioConda'}{'sambamba'},
		     });
    }
    if (! -f $parameter{'condaPath'}."/envs/".$parameter{'condaEnvironment'}."/bin/picard.jar") {
	
	&AddSoftLink({'parameterHashRef' => $parameterHashRef,
		      'FILEHANDLE' => $BASHFILEHANDLE,
		      'binary' => q?../share/picard-?.${$parameterHashRef}{'bioConda'}{'picard'}.q?-1/picard.jar?,
		      'softLink' => "picard.jar",
		     });
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
	    
	    ## Removing specific Perl version
	    print $FILEHANDLE "### Removing specific Perl version\n";
	    print $FILEHANDLE q?rm -rf $HOME/perl-?.${$parameterHashRef}{'perl'};
	    print $FILEHANDLE "\n\n";
	    
	    &InstallPerlCpnam($parameterHashRef, $BASHFILEHANDLE); 

	    &PerlModules($parameterHashRef, $BASHFILEHANDLE,);
	}
    }
    else {
    	
	&InstallPerlCpnam($parameterHashRef, $BASHFILEHANDLE, "AddPath");

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
    print $FILEHANDLE "-O sambamba_v".${$parameterHashRef}{'sambamba'}."_linux.tar.bz2";  #Dowload outfile
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

    if (! -f $parameter{'condaPath'}."/envs/".$parameter{'condaEnvironment'}."/bin/sambamba") {
	
	&AddSoftLink({'parameterHashRef' => $parameterHashRef,
		      'FILEHANDLE' => $BASHFILEHANDLE,
		      'binary' => "sambamba_v".${$parameterHashRef}{'bioConda'}{'sambamba'},
		      'softLink' => "sambamba",
		     });
    }

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
    print $FILEHANDLE "-O vcftools-".${$parameterHashRef}{'vcfTools'}.".tar.gz";  #Dowload outfile
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
    print $FILEHANDLE "-O bedtools-".${$parameterHashRef}{'bedTools'}.".tar.gz";  #Dowload outfile
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

    ## Install sambamba
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
    print $FILEHANDLE "wget https://www.cog-genomics.org/static/bin/plink".${$parameterHashRef}{'plink'}."/plink_linux_x86_64.zip ";
    print $FILEHANDLE "-O plink-".${$parameterHashRef}{'plink'}."-x86_64.zip";  #Dowload outfile
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
    print $FILEHANDLE "-O snpEff_".${$parameterHashRef}{'snpEff'}."_core.zip";  #Dowload outfile
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "unzip snpEff_".${$parameterHashRef}{'snpEff'}."_core.zip";
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    print $FILEHANDLE "## Make available from conda environment\n";
    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?snpEff/snpEff.jar ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/?;
    print $FILEHANDLE "\n\n";

    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?snpEff/SnpSift.jar ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/?;
    print $FILEHANDLE "\n\n";
    
    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?snpEff/snpEff.config ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/?;
    print $FILEHANDLE "\n\n";

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
    print $FILEHANDLE "-O VariantEffectPredictor-".${$parameterHashRef}{'variantEffectPredictor'}.".zip";  #Dowload outfile
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

    print $FILEHANDLE "ln -s ";
    print $FILEHANDLE $binary.q? ?.$softLink;
    print $FILEHANDLE "\n\n";

    ## Move to back
    print $FILEHANDLE "## Move to original working directory\n";
    print $FILEHANDLE "cd ".$pwd;
    print $FILEHANDLE "\n\n";
}
