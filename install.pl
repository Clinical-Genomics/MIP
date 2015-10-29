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
           -bvc/--bioConda Set the module version of the programs that can be installed with bioConda (i.e. 'bwa=0.7.12')
           -pm/perlModules Set the perl modules to be installed via cpanm (comma sep)
           -pip/--pip Set the module version of the programs that can be installed with pip (i.e. 'genmod=3.3.3')
           -sbb/sambamba Set the sambamba version (Default: "0.5.9")
           -vct/--vcfTools Set the vcftools version (Default: "0.1.14")
           -vt/--vt Set the vt version (Default: "0.57")
           -plk/--plink  Set the plink version (Default: "1.07")
           -vep/--VariantEffectPredictor Set the VEP version (Default: "82")
	   -vepc/--vepDirectoryCache Specify the cache directory to use (whole path; defaults to "~/miniconda/envs/condaEnvironment/ensembl-tools-release-VariantEffectPredictorVersion/cache")
           -ppd/--printParameterDefaults Print the parameter defaults
           -h/--help Display this help message   
           -v/--version Display version
        };    
}
my ($installDirectory) = (".MIP");

my %parameter; 

### Set parameter default

##Conda
$parameter{'condaEnvironment'} = "mip";
$parameter{'bioConda'}{'bwa'} = "0.7.12";
$parameter{'bioConda'}{'fastqc'} = "0.11.2";
$parameter{'bioConda'}{'samtools'} = "1.2";
$parameter{'bioConda'}{'bcftools'} = "1.2";
$parameter{'bioConda'}{'snpeff'} = "4.1";
$parameter{'bioConda'}{'picard'} = "1.139";
$parameter{'bioConda'}{'bedtools'} = "2.25";
$parameter{'bioConda'}{'mosaik'} = "2.2.26";

##Perl Modules
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
    ];

## PIP
$parameter{'pip'}{'genmod'} = "3.3.3";
$parameter{'pip'}{'chanjo'} = "3.0.1";
$parameter{'pip'}{'cosmid'} = "0.4.9.1";
$parameter{'pip'}{'python-Levenshtein'} = "0.12.0";

## Programs currently not supported by conda or other packet manager
$parameter{'sambamba'} = "0.5.9";
$parameter{'vcfTools'} = "0.1.14";
$parameter{'vt'} = "0.57";
$parameter{'plink'} = "1.07";
$parameter{'VariantEffectPredictor'} = "82";
$parameter{'vepDirectoryCache'} = q?~/miniconda/envs/?.$parameter{'condaEnvironment'}.q?/ensembl-tools-release-?.$parameter{'VariantEffectPredictor'}.q?/cache?;  #Cache directory;

my $installVersion = "0.0.1";

###User Options
GetOptions('env|condaEnvironment:s'  => \$parameter{'condaEnvironment'},
	   'bcv|bioConda=s'  => \%{$parameter{'bioConda'}},
	   'pm|perlModules:s'  => \@{$parameter{'perlModules'}},  #Comma separated list
	   'pip|pip=s'  => \%{$parameter{'pip'}},
	   'sbb|sambamba:s'  => \$parameter{'sambamba'},
	   'vct|vcfTools:s'  => \$parameter{'vcfTools'},
	   'vt|vt:s'  => \$parameter{'vt'},
	   'plk|plink:s'  => \$parameter{'plink'},
	   'vep|VariantEffectPredictor:s'  => \$parameter{'VariantEffectPredictor'},
	   'vepc|vepDirectoryCache:s'  => \$parameter{'vepDirectoryCache'},  #path to vep cache dir
	   'ppd|printParameterDefaults'  => sub { &PrintParameters(\%parameter); exit;},  #Display parameter defaults
	   'h|help' => sub { print STDOUT $USAGE, "\n"; exit;},  #Display help text
	   'v|version' => sub { print STDOUT "\ninstall.pl ".$installVersion, "\n\n"; exit;},  #Display version number
    );

###MAIN###


#my $LOGFILEHANDLE = &OpenLogFile("MIP_installation.log");

my $BASHFILEHANDLE = &CreateBashFile("mip.sh");

&CreateCondaEnvironment(\%parameter, $BASHFILEHANDLE);

&InstallCpanmAndModules(\%parameter, $BASHFILEHANDLE);

&PipInstall(\%parameter, $BASHFILEHANDLE);

&Sambamba(\%parameter, $BASHFILEHANDLE);

&VcfTools(\%parameter, $BASHFILEHANDLE);

&VT(\%parameter, $BASHFILEHANDLE);

&Plink(\%parameter, $BASHFILEHANDLE);

&VariantEffectPredictor(\%parameter, $BASHFILEHANDLE);

close($BASHFILEHANDLE);
#close($LOGFILEHANDLE);

###SubRoutines###

sub CreateBashFile {

    my $fileName = $_[0];

    my $FILEHANDLE = IO::Handle->new();  #Create anonymous filehandle
    
    ## Open batch file
    open ($FILEHANDLE, ">",$fileName) or die("Can't write to '".$fileName."' :".$!."\n");

    print $FILEHANDLE "#!/usr/bin/env bash", "\n\n";
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


sub CreateCondaEnvironment {

    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];

    ## Create MIP environment
    print $FILEHANDLE "### Creating Conda Environment: ".${$parameterHashRef}{'condaEnvironment'}."\n";
    print $FILEHANDLE "conda create ";
    print $FILEHANDLE "-y ";
    print $FILEHANDLE "-n ".${$parameterHashRef}{'condaEnvironment'}." "; 
    print $FILEHANDLE "-c bioconda ";

    ## Install all bioConda packages
    foreach my $program (keys %{${$parameterHashRef}{'bioConda'}}) {

	print $FILEHANDLE $program."=".${$parameterHashRef}{'bioConda'}{$program}." ";
    }
    print $FILEHANDLE "\n\n";
}


sub InstallCpanmAndModules {

    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];

    &ActivateCondaEnvironment($parameterHashRef, $FILEHANDLE);
    
    ## Install Cpanm
    print $FILEHANDLE "### Install Cpanm in conda environment: ".${$parameterHashRef}{'condaEnvironment'}."\n";
    print $FILEHANDLE "conda install ";
    print $FILEHANDLE "-y ";
    print $FILEHANDLE "-c https://conda.anaconda.org/dan_blanchard perl-app-cpanminus ";
    print $FILEHANDLE "\n\n";

    ## Install Perl modules via cpanm
    print $FILEHANDLE "## Install Perl modules via cpanm\n";
    print $FILEHANDLE "cpanm ";
    print $FILEHANDLE join(" ", @{${$parameterHashRef}{'perlModules'}})." ";
    print $FILEHANDLE "\n\n";

    &DeactivateCondaEnvironment($FILEHANDLE);
}


sub PipInstall {

    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];

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
    print $FILEHANDLE q?sambamba_v?.${$parameterHashRef}{'sambamba'}.q? ~/miniconda/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/?;
    print $FILEHANDLE "\n\n";

    &CleanUpModuleInstall($FILEHANDLE, $pwd);
}


sub VcfTools {

    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];

    my $pwd = cwd();

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
    my $filePath = glob(q?~/miniconda/envs/?.${$parameterHashRef}{'condaEnvironment'});

    print $FILEHANDLE "## Configure\n";
    print $FILEHANDLE q?./configure --prefix=?.$filePath;
    print $FILEHANDLE "\n";

    print $FILEHANDLE "make";
    print $FILEHANDLE "\n";

    print $FILEHANDLE "make install";
    print $FILEHANDLE "\n\n";

    &CleanUpModuleInstall($FILEHANDLE, $pwd);
}


sub VT {

    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];

    my $pwd = cwd();

    ## Install sambamba
    print $FILEHANDLE "### Install VT\n\n";

    &CreateInstallDirectory($FILEHANDLE);
    
    ## Download
    print $FILEHANDLE "## Download VT\n";

    print $FILEHANDLE "wget --quiet https://github.com/atks/vt/archive/".${$parameterHashRef}{'vt'}.".tar.gz ";
    print $FILEHANDLE "-O vt-".${$parameterHashRef}{'vt'}.".tar.gz";  #Dowload outfile
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "tar xvf vt-".${$parameterHashRef}{'vt'}.".tar.gz";
    print $FILEHANDLE "\n\n";

    ## Move to vt directory
    print $FILEHANDLE "## Move to vt directory\n";
    print $FILEHANDLE "cd vt-".${$parameterHashRef}{'vt'};
    print $FILEHANDLE "\n\n";

    ## Configure

    print $FILEHANDLE "make";
    print $FILEHANDLE "\n";

    print $FILEHANDLE "make test";
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    print $FILEHANDLE "## Make available from conda environment\n";
    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?vt ~/miniconda/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/?;
    print $FILEHANDLE "\n\n";

    &CleanUpModuleInstall($FILEHANDLE, $pwd);
}


sub Plink {

    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];

    my $pwd = cwd();

    ## Install Plink
    print $FILEHANDLE "### Install Plink\n\n";

    &CreateInstallDirectory($FILEHANDLE);
    
    ## Download
    print $FILEHANDLE "## Download Plink\n";

    print $FILEHANDLE "wget --quiet http://pngu.mgh.harvard.edu/~purcell/plink/dist/plink-".${$parameterHashRef}{'plink'}."-x86_64.zip ";
    print $FILEHANDLE "-O plink-".${$parameterHashRef}{'plink'}."-x86_64.zip";  #Dowload outfile
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "unzip plink-".${$parameterHashRef}{'plink'}."-x86_64.zip";
    print $FILEHANDLE "\n\n";

    ## Move to plink directory
    print $FILEHANDLE "## Move to plink directory\n";
    print $FILEHANDLE "cd plink-".${$parameterHashRef}{'plink'}."-x86_64";
    print $FILEHANDLE "\n\n";

    ## Make available from conda environment
    print $FILEHANDLE "## Make available from conda environment\n";
    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?plink ~/miniconda/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/?;
    print $FILEHANDLE "\n\n";

    &CleanUpModuleInstall($FILEHANDLE, $pwd);
}


sub VariantEffectPredictor {

    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];

    my $pwd = cwd();

    ## Install VEP
    print $FILEHANDLE "### Install VariantEffectPredictor\n\n";

    &ActivateCondaEnvironment($parameterHashRef, $FILEHANDLE);

    ##Make sure that the cache directory exists
    print $FILEHANDLE "mkdir -p ".${$parameterHashRef}{'vepDirectoryCache'}." ";  #Cache directory
    print $FILEHANDLE "\n\n";

    ## Move to miniconda environment
    print $FILEHANDLE q?cd ~/miniconda/envs/?.${$parameterHashRef}{'condaEnvironment'};
    print $FILEHANDLE "\n\n";

    ## Download
    print $FILEHANDLE "## Download VEP\n";
    print $FILEHANDLE "wget --quiet https://github.com/Ensembl/ensembl-tools/archive/release/".${$parameterHashRef}{'VariantEffectPredictor'}.".zip ";
    print $FILEHANDLE "-O VariantEffectPredictor-".${$parameterHashRef}{'VariantEffectPredictor'}.".zip";  #Dowload outfile
    print $FILEHANDLE "\n\n";

    ## Extract
    print $FILEHANDLE "## Extract\n";
    print $FILEHANDLE "unzip VariantEffectPredictor-".${$parameterHashRef}{'VariantEffectPredictor'}.".zip";
    print $FILEHANDLE "\n\n";    

    ## Move to VariantEffectPredictor directory
    print $FILEHANDLE "## Move to VariantEffectPredictor directory\n";
    print $FILEHANDLE "cd ensembl-tools-release-".${$parameterHashRef}{'VariantEffectPredictor'}."/scripts/variant_effect_predictor/";
    print $FILEHANDLE "\n\n";

    ## Install VEP
    print $FILEHANDLE "## Install VEP\n";
    print $FILEHANDLE "perl INSTALL.pl ";
    print $FILEHANDLE "--AUTO alcf ";  #a (API), l (FAIDX/htslib), c (cache), f (FASTA), p (plugins)
    print $FILEHANDLE "-c ".${$parameterHashRef}{'vepDirectoryCache'}." ";  #Cache directory
    print $FILEHANDLE "-s homo_sapiens ";
    print $FILEHANDLE "--ASSEMBLY GRCh37 ";
    print $FILEHANDLE "\n\n";

    ## Clean up
    print $FILEHANDLE "## Clean up\n";
    print $FILEHANDLE q?cd ~/miniconda/envs/?.${$parameterHashRef}{'condaEnvironment'};
    print $FILEHANDLE "\n\n";

    print $FILEHANDLE "rm -rf VariantEffectPredictor-".${$parameterHashRef}{'VariantEffectPredictor'}.".zip";;
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
