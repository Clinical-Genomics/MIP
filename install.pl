#!/usr/bin/env perl

##Assumes you have a working conda installation

use strict;
use warnings;
use warnings qw( FATAL utf8 );
use utf8;
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );

use Getopt::Long;
use Cwd;
use FindBin qw($Bin); #Find directory of script
use vars qw($USAGE);
use IO::Handle;
use File::Spec::Functions qw(catfile), qw(catdir);

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
           -vepp/--variantEffectPredictorPlugin Supply a comma separated list of VEP plugins (Default: "UpDownDistance,LoFtool,LoF")
           -cnv/--CNVnator Set the CNVnator version (Default: "0.3.2")
           -ftr/--FindTranslocations Set the FindTranslocations version (Default: "0")

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
$parameter{'bioCondaMantaPatch'} = "-0";
$parameter{'bioConda'}{'multiqc'} = "0.4";
$parameter{'bioConda'}{'plink2'} = "1.90b3.35";
$parameter{'bioConda'}{'gcc'} = "4.8.5";
$parameter{'bioConda'}{'cmake'} = "3.3.1";
$parameter{'bioConda'}{'boost'} = "1.57.0";
$parameter{'bioCondaBoostPatch'} = "-4";


##Perl Modules
$parameter{'perlInstall'} = 0;
$parameter{'perl'} = "5.18.2";
$parameter{'perlModules'} = ["Modern::Perl",  #MIP
			     "IPC::System::Simple",  #MIP
			     "Path::Iterator::Rule",  #MIP
			     "YAML",  #MIP
			     "Log::Log4perl",  #MIP
			     "Set::IntervalTree",  # MIP/vcfParser.pl
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
			     "Bio::Root::Version",  #VEP
			     "Module::Build", #VEP
    ];

## PIP
$parameter{'pip'}{'genmod'} = "3.4.5";
$parameter{'pip'}{'chanjo'} = "3.3.1";
$parameter{'pip'}{'cosmid'} = "0.4.9.1";
$parameter{'pip'}{'python-Levenshtein'} = "0.12.0";

## Programs currently installable by SHELL
$parameter{'MIPScripts'} = "Your current MIP version";
$parameter{'picardTools'} = "2.0.1";
$parameter{'sambamba'} = "0.5.9";
$parameter{'vcfTools'} = "0.1.14";
$parameter{'bedTools'} = "2.25.0";
$parameter{'vt'} = "gitRepo";
$parameter{'plink2'} = "160316";
$parameter{'snpEff'} = "v4_2";
$parameter{'snpEffGenomeVersion'} = "GRCh37.75";
$parameter{'variantEffectPredictor'} = "83";
$parameter{'vepDirectoryCache'} = $parameter{'condaPath'}.q?/envs/?.$parameter{'condaEnvironment'}.q?/ensembl-tools-release-?.$parameter{'variantEffectPredictor'}.q?/cache?;  #Cache directory;
$parameter{'variantEffectPredictorPlugin'} = "UpDownDistance,LoFtool,LoF";
$parameter{'CNVnator'} = "0.3.2";
$parameter{'FindTranslocations'} = "0";

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
	   'plk|plink2:s' => \$parameter{'plink2'},
	   'vep|variantEffectPredictor:s' => \$parameter{'variantEffectPredictor'},
	   'vepc|vepDirectoryCache:s' => \$parameter{'vepDirectoryCache'},  #path to vep cache dir
	   'vepp|variantEffectPredictorPlugin:s' => \$parameter{'variantEffectPredictorPlugin'},  #Comma sep string
	   'cnv|CNVnator:s' => \$parameter{'CNVnator'},
	   'ftr|FindTranslocations:s' => \$parameter{'FindTranslocations'},
	   'pbc|preferBioConda:s' => \$parameter{'preferBioConda'},  # Bioconda will used for overlapping shell and biconda installlations
	   'ppd|printParameterDefaults' => sub { &PrintParameters(\%parameter); exit;},  #Display parameter defaults
	   'u|update:n' => \$parameter{'update'},
	   'sp|selectPrograms:s' => \@{$parameter{'selectPrograms'}},  #Comma sep string
	   'h|help' => sub { say STDOUT $USAGE; exit;},  #Display help text
	   'v|version' => sub { say STDOUT "\ninstall.pl ".$installVersion, "\n"; exit;},  #Display version number
    );

###MAIN###

#my $LOGFILEHANDLE = &OpenLogFile("MIP_installation.log");

my $BASHFILEHANDLE = &CreateBashFile("mip.sh");

&CreateConda(\%parameter, $BASHFILEHANDLE);

&CreateCondaEnvironment(\%parameter, $BASHFILEHANDLE);

if (scalar(@{$parameter{'selectPrograms'}}) > 0) {
    
    if ( ( grep {$_ eq "perl"} @{$parameter{'selectPrograms'}} ) ) { #If element is part of array

	&Perl(\%parameter, $BASHFILEHANDLE);
    }
}
else {
    
    &Perl(\%parameter, $BASHFILEHANDLE);
}


&PipInstall(\%parameter, $BASHFILEHANDLE);

if ($parameter{'preferBioConda'} != 1) {

    if (scalar(@{$parameter{'selectPrograms'}}) > 0) {
	
	if ( ( grep {$_ eq "picardTools"} @{$parameter{'selectPrograms'}} ) ) { #If element is part of array
	    
	    &PicardTools(\%parameter, $BASHFILEHANDLE);
	}
	if ( ( grep {$_ eq "sambamba"} @{$parameter{'selectPrograms'}} ) ) { #If element is part of array

	    &Sambamba(\%parameter, $BASHFILEHANDLE);
	}
	if ( ( grep {$_ eq "bedTools"} @{$parameter{'selectPrograms'}} ) ) { #If element is part of array

	    &BedTools(\%parameter, $BASHFILEHANDLE);
	}
	if ( ( grep {$_ eq "vt"} @{$parameter{'selectPrograms'}} ) ) { #If element is part of array

	    &VT(\%parameter, $BASHFILEHANDLE);
	}
	if ( ( grep {$_ eq "snpEff"} @{$parameter{'selectPrograms'}} ) ) { #If element is part of array

	    &SnpEff(\%parameter, $BASHFILEHANDLE);
	}
	if ( ( grep {$_ eq "plink2"} @{$parameter{'selectPrograms'}} ) ) { #If element is part of array
	
	    &Plink2(\%parameter, $BASHFILEHANDLE);
	}
    }
    else {
	
	&PicardTools(\%parameter, $BASHFILEHANDLE);
	
	&Sambamba(\%parameter, $BASHFILEHANDLE);
	
	&BedTools(\%parameter, $BASHFILEHANDLE);
	
	&VT(\%parameter, $BASHFILEHANDLE);

	&Plink2(\%parameter, $BASHFILEHANDLE);
	
	&SnpEff(\%parameter, $BASHFILEHANDLE);
    }
}

if (scalar(@{$parameter{'selectPrograms'}}) > 0) {

    if ( ( grep {$_ eq "MIPScripts"} @{$parameter{'selectPrograms'}} ) ) { #If element is part of array
	&MIPScripts(\%parameter, $BASHFILEHANDLE);
    }
    if ( ( grep {$_ eq "vcfTools"} @{$parameter{'selectPrograms'}} ) ) { #If element is part of array
	
	&VcfTools(\%parameter, $BASHFILEHANDLE);
    }
    if ( ( grep {$_ eq "variantEffectPredictor"} @{$parameter{'selectPrograms'}} ) ) { #If element is part of array
	
	&VariantEffectPredictor(\%parameter, $BASHFILEHANDLE);
    }
    if ( ( grep {$_ eq "CNVnator"} @{$parameter{'selectPrograms'}} ) ) { #If element is part of array
	
	&CNVnator(\%parameter, $BASHFILEHANDLE);
    }
    if ( ( grep {$_ eq "FindTranslocations"} @{$parameter{'selectPrograms'}} ) ) { #If element is part of array
	
	&FindTranslocations(\%parameter, $BASHFILEHANDLE);
    }
}
else {
    
    &MIPScripts(\%parameter, $BASHFILEHANDLE);

    &VcfTools(\%parameter, $BASHFILEHANDLE);
    
    &VariantEffectPredictor(\%parameter, $BASHFILEHANDLE);

    &CNVnator(\%parameter, $BASHFILEHANDLE);

    &FindTranslocations(\%parameter, $BASHFILEHANDLE);
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

    say $FILEHANDLE "#!/usr/bin/env bash", "\n";
 
    say STDOUT "Will write install instructions to '".$pwd."/".$fileName, "'";

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

	    say STDOUT $key." ".${$parameterHashRef}{$key};
	}
	elsif (ref(${$parameterHashRef}{$key})=~/HASH/) {
	    
	    foreach my $program (keys %{${$parameterHashRef}{$key}}) {
		
		say STDOUT $key." ".$program.": ".${$parameterHashRef}{$key}{$program};
	    }
	}
	elsif (ref(${$parameterHashRef}{$key})=~/ARRAY/)  {

	    say STDOUT $key.": ".join(" ", @{${$parameterHashRef}{$key}});
	}
    }
}


sub CreateConda {

    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];

    my $program = "conda";

    if ($ENV{PATH}=~/conda/) {

	say STDERR "ProgramCheck: ".$program." installed";
    }
    else {
	
	say STDERR "Could not detect ".$program." in your PATH";
	exit 1;
    }

    ## Check Conda path
    if (! -d $parameter{'condaPath'}) {

	say STDERR "Could not find miniconda directory in: ".$parameter{'condaPath'};
	exit 1;
    }

    say STDERR "Writting install instructions for Conda packages";

    ## Update Conda
    say $FILEHANDLE "### Update Conda";
    print $FILEHANDLE "conda update -y conda ";
    say $FILEHANDLE "\n";
}

sub CreateCondaEnvironment {

    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];

    ## Check Conda environment
    if (! -d $parameter{'condaPath'}."/envs/".$parameter{'condaEnvironment'}) {

	## Create conda environment
	say $FILEHANDLE "### Creating Conda Environment and install: ".${$parameterHashRef}{'condaEnvironment'};
	print $FILEHANDLE "conda create -n ".${$parameterHashRef}{'condaEnvironment'}." ";
	print $FILEHANDLE "-y ";
	print $FILEHANDLE "pip ";
	say $FILEHANDLE "\n";
    }
    
    ## Install into conda environment
    say $FILEHANDLE "### Installing into Conda Environment: ".${$parameterHashRef}{'condaEnvironment'};
    print $FILEHANDLE "conda install ";	
    print $FILEHANDLE "-n ".${$parameterHashRef}{'condaEnvironment'}." ";
    print $FILEHANDLE "-y ";
    print $FILEHANDLE "-c bioconda ";
    
    ## Install all bioConda packages
    foreach my $program (keys %{${$parameterHashRef}{'bioConda'}}) {
	    
	print $FILEHANDLE $program."=".${$parameterHashRef}{'bioConda'}{$program}." ";
    }

    say $FILEHANDLE "\n";

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
	    unless (-d $parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/share/snpeff-?.${$parameterHashRef}{'bioConda'}{'snpeff'}.${$parameterHashRef}{'bioCondaSnpeffPatch'}.q?/data/?.${$parameterHashRef}{'snpEffGenomeVersion'}) {
		
		&SnpEffDownload({'parameterHashRef' => $parameterHashRef,
				 'FILEHANDLE' => $BASHFILEHANDLE,
				});
	    }	
	}
	if ($program eq "manta") {
	    
	    &AddSoftLink({'parameterHashRef' => $parameterHashRef,
                          'FILEHANDLE' => $BASHFILEHANDLE,
                          'binary' => q?../share/manta-?.${$parameterHashRef}{'bioConda'}{'manta'}.${$parameterHashRef}{'bioCondaMantaPatch'}.q?/bin/configManta.py?,
                          'softLink' => "configManta.py",
                         });
	    
	    &EnableExecutable({'parameterHashRef' => $parameterHashRef,
			       'FILEHANDLE' => $BASHFILEHANDLE,
			       'binary' => q?configManta.py?,
			 });
	    &AddSoftLink({'parameterHashRef' => $parameterHashRef,
			  'FILEHANDLE' => $BASHFILEHANDLE,
			  'binary' => q?../share/manta-?.${$parameterHashRef}{'bioConda'}{'manta'}.${$parameterHashRef}{'bioCondaMantaPatch'}.q?/bin/configManta.py.ini?,
			  'softLink' => "configManta.py.ini",
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
	    
	    say STDERR "Found perl-".${$parameterHashRef}{'perl'}.". in your path";
	    say STDERR q?Skipping writting installation for perl-?.${$parameterHashRef}{'perl'};  
	}
	else {

	    if (${$parameterHashRef}{'perlInstall'} == 1) {
	    
		## Removing specific Perl version
		say $FILEHANDLE "### Removing specific Perl version";
		print $FILEHANDLE q?rm -rf $HOME/perl-?.${$parameterHashRef}{'perl'};
		say $FILEHANDLE "\n";
		
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

    say STDERR "Writting install instructions for Perl and Cpanm";
    
    ## Install specific Perl version
    say $FILEHANDLE "### Install specific Perl version";
    
    ## Move to Home
    say $FILEHANDLE "## Move HOME";
    print $FILEHANDLE q?cd $HOME?;
    say $FILEHANDLE "\n";
    
    ## Download
    say $FILEHANDLE "## Download Perl";
    print $FILEHANDLE "wget --quiet http://www.cpan.org/src/5.0/perl-".${$parameterHashRef}{'perl'}.".tar.gz ";
    print $FILEHANDLE "-O perl-".${$parameterHashRef}{'perl'}.".tar.gz";  #Dowload outfile
    say $FILEHANDLE "\n";
    
    ## Extract
    say $FILEHANDLE "## Extract";
    print $FILEHANDLE "tar xzf perl-".${$parameterHashRef}{'perl'}.".tar.gz";
    say $FILEHANDLE "\n";
    
    ## Move to perl directory
    say $FILEHANDLE "## Move to perl directory";
    print $FILEHANDLE "cd perl-".${$parameterHashRef}{'perl'};
    say $FILEHANDLE "\n";
    
    ## Configure
    say $FILEHANDLE "## Configure";
    say $FILEHANDLE q?./Configure -des -Dprefix=$HOME/perl-?.${$parameterHashRef}{'perl'};    
    say $FILEHANDLE "make";
    say $FILEHANDLE "make test";    
    say $FILEHANDLE "make install", "\n";
    
    if ($path) {
	
	## Export path
	say $FILEHANDLE "## Export path";
	print $FILEHANDLE q?echo 'export PATH=$HOME/perl-?.${$parameterHashRef}{'perl'}.q?/:$PATH' >> ~/.bashrc?;
	say $FILEHANDLE "\n";
	print $FILEHANDLE q?export PATH=$HOME/perl-?.${$parameterHashRef}{'perl'}.q?/:$PATH?;  #Use newly installed perl
	say $FILEHANDLE "\n";
    }

    ## Remove tar file
    say $FILEHANDLE "## Remove tar file";
    print $FILEHANDLE "cd && rm perl-".${$parameterHashRef}{'perl'}.".tar.gz";
    say $FILEHANDLE "\n";
    
    ## Move to back
    say $FILEHANDLE "## Move to original working directory";
    print $FILEHANDLE "cd ".$pwd;
    say $FILEHANDLE "\n";

    #if ($path) {

	print $FILEHANDLE q?echo 'eval `perl -I ~/perl-?.${$parameterHashRef}{'perl'}.q?/lib/perl5/ -Mlocal::lib=~/perl-?.${$parameterHashRef}{'perl'}.q?/`' >> ~/.bash_profile ?;  #Add at start-up
	say $FILEHANDLE "\n";
	print $FILEHANDLE q?echo 'export PERL_UNICODE=SAD' >> ~/.bash_profile ?;  #Add at start-up
	say $FILEHANDLE "\n";
    #}

    ## Install Perl modules via cpanm
    say $FILEHANDLE "## Install cpanm";
    print $FILEHANDLE q?wget -O- http://cpanmin.us | perl - -l $HOME/perl-?.${$parameterHashRef}{'perl'}.q?/bin App::cpanminus --local-lib=~/perl-?.${$parameterHashRef}{'perl'}.q?/ local::lib ?;
    say $FILEHANDLE "\n";

    ## Use newly installed perl
    print $FILEHANDLE q?eval `perl -I ~/perl-?.${$parameterHashRef}{'perl'}.q?/lib/perl5/ -Mlocal::lib=~/perl-?.${$parameterHashRef}{'perl'}.q?/` ?;
    say $FILEHANDLE "\n";

    ## Use newly installed perl
    print $FILEHANDLE q?PERL5LIB=~/perl-?.${$parameterHashRef}{'perl'}.q?/lib/perl5?;
    say $FILEHANDLE "\n";
}
    

sub PerlModules {
    
    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];
    
    ## Install Perl modules via cpanm
    say $FILEHANDLE "## Install Perl modules via cpanm";
    print $FILEHANDLE "cpanm ";
    print $FILEHANDLE join(" ", @{${$parameterHashRef}{'perlModules'}})." ";
    say $FILEHANDLE "\n";
}


sub PipInstall {

    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];

    say STDERR "Writting install instructions for pip packages";

    ## Install PIP packages in conda environment
    say $FILEHANDLE "### Install PIP packages in conda environment: ".${$parameterHashRef}{'condaEnvironment'};
    &ActivateCondaEnvironment($parameterHashRef, $FILEHANDLE);

    ## Install PIP packages
    say $FILEHANDLE "## Install PIP packages";
    print $FILEHANDLE "pip install ";

    ## Install all PIP packages
    foreach my $program (keys %{${$parameterHashRef}{'pip'}}) {

	print $FILEHANDLE $program."==".${$parameterHashRef}{'pip'}{$program}." ";
    }
    say $FILEHANDLE "\n";

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
    say $FILEHANDLE "### Install Picard";

    &CreateInstallDirectory($FILEHANDLE);
    
    ## Download
    say $FILEHANDLE "## Download Picard";
    print $FILEHANDLE "wget --quiet https://github.com/broadinstitute/picard/releases/download/".${$parameterHashRef}{'picardTools'}."/picard-tools-".${$parameterHashRef}{'picardTools'}.".zip ";
    print $FILEHANDLE "-O picard-tools-".${$parameterHashRef}{'picardTools'}.".zip";  #Download outfile
    say $FILEHANDLE "\n";

    ## Extract
    say $FILEHANDLE "## Extract";
    print $FILEHANDLE "unzip picard-tools-".${$parameterHashRef}{'picardTools'}.".zip";
    say $FILEHANDLE "\n";

    ## Make available from conda environment
    if (-d $parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/share/picard-tools-?.${$parameterHashRef}{'picardTools'}) {

	print $FILEHANDLE "rm -rf ".$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/share/picard-tools-?.${$parameterHashRef}{'picardTools'};
	say $FILEHANDLE "\n";
    }

    say $FILEHANDLE "## Make available from conda environment";
    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?picard-tools-?.${$parameterHashRef}{'picardTools'}.q? ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/share/?;
    say $FILEHANDLE "\n";

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
    say $FILEHANDLE "### Install sambamba";

    &CreateInstallDirectory($FILEHANDLE);

    ## Download
    say $FILEHANDLE "## Download sambamba release";
    print $FILEHANDLE q?wget --quiet https://github.com/lomereiter/sambamba/releases/download/v?.${$parameterHashRef}{'sambamba'}.q?/sambamba_v?.${$parameterHashRef}{'sambamba'}.q?_linux.tar.bz2 ?;
    print $FILEHANDLE "-O sambamba_v".${$parameterHashRef}{'sambamba'}."_linux.tar.bz2";  #Download outfile
    say $FILEHANDLE "\n";

    ## Decompress
    say $FILEHANDLE "## Decompress sambamba file";
    print $FILEHANDLE "bzip2 ";
    print $FILEHANDLE "-f ";  #Force
    print $FILEHANDLE "-d ";  #Decompress
    print $FILEHANDLE "sambamba_v".${$parameterHashRef}{'sambamba'}."_linux.tar.bz2";
    say $FILEHANDLE "\n";

    ## Extract files
    say $FILEHANDLE "## Extract files";
    print $FILEHANDLE "tar xvf sambamba_v".${$parameterHashRef}{'sambamba'}."_linux.tar";
    say $FILEHANDLE "\n";

    ## Make executable
    say $FILEHANDLE "## Make executable";
    print $FILEHANDLE "chmod 755 ";
    print $FILEHANDLE "sambamba_v".${$parameterHashRef}{'sambamba'};
    say $FILEHANDLE "\n";

    ## Make available from conda environment
    say $FILEHANDLE "## Make available from conda environment";
    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?sambamba_v?.${$parameterHashRef}{'sambamba'}.q? ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/?;
    say $FILEHANDLE "\n";

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
    say $FILEHANDLE "### Install vcfTools";

    &CreateInstallDirectory($FILEHANDLE);
    
    ## Download
    say $FILEHANDLE "## Download vcfTools";
    print $FILEHANDLE "wget --quiet https://github.com/vcftools/vcftools/releases/download/v".${$parameterHashRef}{'vcfTools'}."/vcftools-".${$parameterHashRef}{'vcfTools'}.".tar.gz ";
    print $FILEHANDLE "-O vcftools-".${$parameterHashRef}{'vcfTools'}.".tar.gz";  #Download outfile
    say $FILEHANDLE "\n";

    ## Extract
    say $FILEHANDLE "## Extract";
    print $FILEHANDLE "tar xvf vcftools-".${$parameterHashRef}{'vcfTools'}.".tar.gz";
    say $FILEHANDLE "\n";

    ## Export PERL5LIB environment variable
    say $FILEHANDLE "## Export PERL5LIB environment variable";
    print $FILEHANDLE q?export PERL5LIB=?.$Bin.q?/vcftools-?.${$parameterHashRef}{'vcfTools'}.q?/src/perl/?;
    say $FILEHANDLE "\n";

    ## Move to vcfTools directory
    say $FILEHANDLE "## Move to vcfTools directory";
    print $FILEHANDLE "cd vcftools-".${$parameterHashRef}{'vcfTools'};
    say $FILEHANDLE "\n";

    ## Configure
    my $filePath = $parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'};

    say $FILEHANDLE "## Configure";
    say $FILEHANDLE q?./configure --prefix=?.$filePath;
    say $FILEHANDLE "make";
    print $FILEHANDLE "make install";
    say $FILEHANDLE "\n";

    ## Move Perl Module
    say $FILEHANDLE "## Move Perl Module";
    print $FILEHANDLE q?cp src/perl/Vcf.pm $HOME/perl-?.${$parameterHashRef}{'perl'}.q?/lib/perl5/?;
    say $FILEHANDLE "\n";

    &CleanUpModuleInstall($FILEHANDLE, $pwd);

    ## Reset perl envionment
    print $FILEHANDLE q?PERL5LIB=~/perl-?.${$parameterHashRef}{'perl'}.q?/lib/perl5?;
    say $FILEHANDLE "\n";
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
    say $FILEHANDLE "### Install bedTools";
    
    &CreateInstallDirectory($FILEHANDLE);
    
    ## Download
    say $FILEHANDLE "## Download bedTools";
    print $FILEHANDLE "wget --quiet https://github.com/arq5x/bedtools".$bedToolsMainVersion."/releases/download/v".${$parameterHashRef}{'bedTools'}."/bedtools-".${$parameterHashRef}{'bedTools'}.".tar.gz ";
    print $FILEHANDLE "-O bedtools-".${$parameterHashRef}{'bedTools'}.".tar.gz";  #Download outfile
    say $FILEHANDLE "\n";
    
    ## Extract
    say $FILEHANDLE "## Extract";
    print $FILEHANDLE "tar xvf bedtools-".${$parameterHashRef}{'bedTools'}.".tar.gz";
    say $FILEHANDLE "\n";

    ## Move to bedtools directory
    say $FILEHANDLE "## Move to bedtools directory";
    print $FILEHANDLE "cd bedtools".$bedToolsMainVersion;
    say $FILEHANDLE "\n";

    print $FILEHANDLE "make";
    say $FILEHANDLE "\n";
       
    ## Make available from conda environment
    say $FILEHANDLE "## Make available from conda environment";
    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?./bin/* ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/?;
    say $FILEHANDLE "\n";
    
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
    say $FILEHANDLE "### Install VT";

    &CreateInstallDirectory($FILEHANDLE);
    
    ## Download
    say $FILEHANDLE "## Download VT";

    print $FILEHANDLE "git clone https://github.com/atks/vt.git ";
    say $FILEHANDLE "\n";

    ## Move to vt directory
    say $FILEHANDLE "## Move to vt directory";
    print $FILEHANDLE "cd vt ";
    say $FILEHANDLE "\n";

    ## Configure
    say $FILEHANDLE "## Configure";
    say $FILEHANDLE "make";
    print $FILEHANDLE "make test";
    say $FILEHANDLE "\n";

    ## Make available from conda environment
    say $FILEHANDLE "## Make available from conda environment";
    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?vt ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/?;
    say $FILEHANDLE "\n";

    &CleanUpModuleInstall($FILEHANDLE, $pwd);
}


sub Plink2 {

    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];

    my $pwd = cwd();

    if (&CheckCondaBinFileExists($parameterHashRef, "plink")) {

	return
    }

    ## Install Plink
    say $FILEHANDLE "### Install Plink";

    &CreateInstallDirectory($FILEHANDLE);
    
    ## Download
    say $FILEHANDLE "## Download Plink";
    print $FILEHANDLE "wget --quiet https://www.cog-genomics.org/static/bin/plink".${$parameterHashRef}{'plink2'}."/plink_linux_x86_64.zip ";
    print $FILEHANDLE "-O plink-".${$parameterHashRef}{'plink2'}."-x86_64.zip";  #Download outfile
    say $FILEHANDLE "\n";

    ## Extract
    say $FILEHANDLE "## Extract";
    print $FILEHANDLE "unzip plink-".${$parameterHashRef}{'plink2'}."-x86_64.zip";
    say $FILEHANDLE "\n";

    ## Make available from conda environment
    say $FILEHANDLE "## Make available from conda environment";
    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?plink ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/plink2?;
    say $FILEHANDLE "\n";

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
    say $FILEHANDLE "### Install SnpEff";

    &CreateInstallDirectory($FILEHANDLE);
    
    ## Download
    say $FILEHANDLE "## Download SnpEff";
    print $FILEHANDLE "wget --quiet http://sourceforge.net/projects/snpeff/files/snpEff_".${$parameterHashRef}{'snpEff'}."_core.zip/download ";
    print $FILEHANDLE "-O snpEff_".${$parameterHashRef}{'snpEff'}."_core.zip";  #Download outfile
    say $FILEHANDLE "\n";

    ## Extract
    say $FILEHANDLE "## Extract";
    print $FILEHANDLE "unzip snpEff_".${$parameterHashRef}{'snpEff'}."_core.zip";
    say $FILEHANDLE "\n";

    ## Make available from conda environment
    if (-d $parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/share/snpEff.?.${$parameterHashRef}{'snpEff'}) {
	
	print $FILEHANDLE "rm -rf ".$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/share/snpEff.?.${$parameterHashRef}{'snpEff'};
	say $FILEHANDLE "\n";
    }

    say $FILEHANDLE "## Make available from conda environment";
    print $FILEHANDLE "mkdir -p ".$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/share/snpEff.?.${$parameterHashRef}{'snpEff'};
    say $FILEHANDLE "\n";

    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?snpEff/*.jar ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/share/snpEff.?.${$parameterHashRef}{'snpEff'}.q?/?;
    say $FILEHANDLE "\n";

    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?snpEff/snpEff.config ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/share/snpEff.?.${$parameterHashRef}{'snpEff'}.q?/?;
    say $FILEHANDLE "\n";

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
    unless (-d $parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/share/snpEff.?.${$parameterHashRef}{'snpEff'}.q?/data/?.${$parameterHashRef}{'snpEffGenomeVersion'}) {
	
	&SnpEffDownload({'parameterHashRef' => $parameterHashRef,
			 'FILEHANDLE' => $BASHFILEHANDLE,
			});
    }
    &CleanUpModuleInstall($FILEHANDLE, $pwd);
}


sub VariantEffectPredictor {

    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];

    my $pwd = cwd();
    
    my $minicondaBinDirectory = $parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/ensembl-tools-release-?.${$parameterHashRef}{'variantEffectPredictor'};

    if (-d $minicondaBinDirectory) {

	say STDERR q?Found VariantEffectPredictor in miniconda directory: ?.$minicondaBinDirectory;
	
	if (${$parameterHashRef}{'update'} == 0) {

	    say STDERR "Skipping writting installation process for VariantEffectPredictor";  	    
	    return
	}
	else {

	    ## Removing VariantEffectPredictor
	    say $FILEHANDLE "### Removing VariantEffectPredictor";
	    print $FILEHANDLE q?rm -rf ?.$minicondaBinDirectory;
	    say $FILEHANDLE "\n";
	}
    }
    else {
	
	say STDERR "Writting install instructions for VariantEffectPredictor";
    }

    ## Install VEP
    say $FILEHANDLE "### Install VariantEffectPredictor";

    &ActivateCondaEnvironment($parameterHashRef, $FILEHANDLE);

    ##Make sure that the cache directory exists
    print $FILEHANDLE "mkdir -p ".${$parameterHashRef}{'vepDirectoryCache'}." ";  #Cache directory
    say $FILEHANDLE "\n";

    ## Move to miniconda environment
    print $FILEHANDLE q?cd ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'};
    say $FILEHANDLE "\n";

    ## Download
    say $FILEHANDLE "## Download VEP";
    print $FILEHANDLE "wget --quiet https://github.com/Ensembl/ensembl-tools/archive/release/".${$parameterHashRef}{'variantEffectPredictor'}.".zip ";
    print $FILEHANDLE "-O VariantEffectPredictor-".${$parameterHashRef}{'variantEffectPredictor'}.".zip";  #Download outfile
    say $FILEHANDLE "\n";

    ## Extract
    say $FILEHANDLE "## Extract";
    print $FILEHANDLE "unzip VariantEffectPredictor-".${$parameterHashRef}{'variantEffectPredictor'}.".zip";
    say $FILEHANDLE "\n";    

    ## Move to VariantEffectPredictor directory
    say $FILEHANDLE "## Move to VariantEffectPredictor directory";
    print $FILEHANDLE "cd ensembl-tools-release-".${$parameterHashRef}{'variantEffectPredictor'}."/scripts/variant_effect_predictor/";
    say $FILEHANDLE "\n";

    ## Install VEP
    say $FILEHANDLE "## Install VEP";
    print $FILEHANDLE "perl INSTALL.pl ";
    print $FILEHANDLE "--AUTO alcfp ";  #a (API), l (FAIDX/htslib), c (cache), f (FASTA), p (plugins)
    print $FILEHANDLE "-g ".$parameter{'variantEffectPredictorPlugin'}." ";  #Plugins 
    print $FILEHANDLE "-c ".${$parameterHashRef}{'vepDirectoryCache'}." ";  #Cache directory
    print $FILEHANDLE "-s homo_sapiens ";
    print $FILEHANDLE "--ASSEMBLY GRCh37 ";
    say $FILEHANDLE "\n";

    ##Add LofTool required text file
    say $FILEHANDLE "##Add LofTool required text file";
    print $FILEHANDLE "wget https://raw.githubusercontent.com/Ensembl/VEP_plugins/master/LoFtool_scores.txt ";
    print $FILEHANDLE q?-O $HOME/.vep/Plugins/LoFtool_scores.txt ?;
    say $FILEHANDLE "\n";

    ##Add Lof required perl splice script
    say $FILEHANDLE "##Add Lof required perl splice script";
    print $FILEHANDLE "wget https://raw.githubusercontent.com/konradjk/loftee/master/splice_module.pl ";
    print $FILEHANDLE q?-O $HOME/.vep/Plugins/splice_module.pl ?;
    say $FILEHANDLE "\n";

    ##Add Lof optional human_ancestor_fa
    say $FILEHANDLE "##Add Lof optional human_ancestor_fa";
    print $FILEHANDLE "wget https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz ";
    print $FILEHANDLE q?-O ?.${$parameterHashRef}{'vepDirectoryCache'}.q?/human_ancestor.fa.gz ?;
    say $FILEHANDLE "\n";

    ##Uncompress
    print $FILEHANDLE "bgzip -d ".${$parameterHashRef}{'vepDirectoryCache'}."/human_ancestor.fa.gz ";
    say $FILEHANDLE "\n";

    print $FILEHANDLE "wget https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz.fai ";
    print $FILEHANDLE q?-O ?.${$parameterHashRef}{'vepDirectoryCache'}.q?/human_ancestor.fa.fai ?;
    say $FILEHANDLE "\n";

    ## Clean up
    say $FILEHANDLE "## Clean up";
    print $FILEHANDLE q?cd ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'};
    say $FILEHANDLE "\n";

    print $FILEHANDLE "rm -rf VariantEffectPredictor-".${$parameterHashRef}{'variantEffectPredictor'}.".zip";;
    say $FILEHANDLE "\n";

    ## Moving up
    say $FILEHANDLE "## Moving back to original working directory";
    print $FILEHANDLE "cd ".$pwd;  #Go back to subroutine origin
    say $FILEHANDLE "\n";

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

	say STDERR q?Found Root in miniconda directory: ?.$minicondaBinDirectory;
	
	if (${$parameterHashRef}{'update'} == 0) {

	    say STDERR "Skipping writting installation process for Root";  	    
	    return
	}
	else {

	    ## Removing Root
	    say $FILEHANDLE "### Removing Root";
	    print $FILEHANDLE q?rm -rf ?.$minicondaBinDirectory;
	    say $FILEHANDLE "\n";
	}
    }
    else {
	
	say STDERR "Writting install instructions for Root";
    }

    ## Install Root
    say $FILEHANDLE "### Install CNVnator/Root";

    ## Move to miniconda environment
    print $FILEHANDLE q?cd ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'};
    say $FILEHANDLE "\n";

    ## Download
    say $FILEHANDLE "## Download Root";

    print $FILEHANDLE "wget --quiet https://root.cern.ch/download/root_v5.34.34.Linux-slc6-x86_64-gcc4.4.tar.gz ";  #Currently hardcoded
    print $FILEHANDLE "-O root_v5.34.34.Linux-slc6-x86_64-gcc4.4.tar.gz ";  #Download outfile
    say $FILEHANDLE "\n";

    ## Extract
    say $FILEHANDLE "## Extract";
    print $FILEHANDLE "tar xvf root_v5.34.34.Linux-slc6-x86_64-gcc4.4.tar.gz ";
    say $FILEHANDLE "\n";

    unless ($ENV{PATH}=~/$parameter{'condaPath'}\/envs\/${$parameterHashRef}{'condaEnvironment'}\/root\/bin/) {
	
	## Export path
	say $FILEHANDLE "## Export path";
	print $FILEHANDLE q?echo 'source ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/root/bin/thisroot.sh' >> ~/.bashrc?;
	say $FILEHANDLE "\n";

	## Use newly installed root
	print $FILEHANDLE q?source ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/root/bin/thisroot.sh ?;
	say $FILEHANDLE "\n";
    }
    
    ## Moving up
    say $FILEHANDLE "## Moving back to original working directory";
    print $FILEHANDLE "cd ".$pwd;  #Go back to subroutine origin
    say $FILEHANDLE "\n";

    ## Install CNVNator
    say $FILEHANDLE "### Install CNVnator";

    &CreateInstallDirectory($FILEHANDLE);
    
    ## Download
    say $FILEHANDLE "## Download CNVNator";
    print $FILEHANDLE "wget --quiet https://github.com/abyzovlab/CNVnator/releases/download/v".${$parameterHashRef}{'CNVnator'}."/CNVnator_v".${$parameterHashRef}{'CNVnator'}.".zip ";
    print $FILEHANDLE "-O CNVnator_v".${$parameterHashRef}{'CNVnator'}.".zip";  #Download outfile
    say $FILEHANDLE "\n";

    ## Extract
    say $FILEHANDLE "## Extract";
    print $FILEHANDLE "unzip CNVnator_v".${$parameterHashRef}{'CNVnator'}.".zip";
    say $FILEHANDLE "\n";

    ## Move to CNVnator directory
    say $FILEHANDLE "## Move to CNVnator directory";
    print $FILEHANDLE "cd CNVnator_v".${$parameterHashRef}{'CNVnator'}."/src/samtools/";
    say $FILEHANDLE "\n";

    ## Configure
    say $FILEHANDLE "## Configure CNVnator samTools specific version";
    print $FILEHANDLE "make";
    say $FILEHANDLE "\n";

    say $FILEHANDLE "## Move to CNVnator directory";
    print $FILEHANDLE "cd ..";
    print $FILEHANDLE "\n";

    print $FILEHANDLE "make";
    say $FILEHANDLE "\n";

    ## Make available from conda environment
    say $FILEHANDLE "## Make available from conda environment";
    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?cnvnator ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/?;
    say $FILEHANDLE "\n";

    ## Make available from conda environment
    say $FILEHANDLE "## Make available from conda environment";
    print $FILEHANDLE "mv ";
    print $FILEHANDLE q?cnvnator2VCF.pl ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/?;
    say $FILEHANDLE "\n";

    say $FILEHANDLE "## Make executable from conda environment";
    print $FILEHANDLE "chmod +x ";
    print $FILEHANDLE $parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/cnvnator2VCF.pl?;
    say $FILEHANDLE "\n";
    

    &CleanUpModuleInstall($FILEHANDLE, $pwd);
}


sub FindTranslocations {

    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];

    my $pwd = cwd();

    if (&CheckCondaBinFileExists($parameterHashRef, "FindTranslocations")) {

	return
    }

    ## Install FindTranslocations
    say $FILEHANDLE "### Install FindTranslocations";

    &ActivateCondaEnvironment($parameterHashRef, $FILEHANDLE);

    ## Add to bashrc
    unless (-d $parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/FindTranslocations/bin?) {
	
	## Export path
	say $FILEHANDLE "## Export to bashrc";
	print $FILEHANDLE q?printf '\nif [ -f ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/FindTranslocations/bin/FindTranslocations ]; then\n?;
	print $FILEHANDLE q?\t\texport LD_LIBRARY_PATH=$LD_LIBRARY_PATH:?.$parameter{'condaPath'}.q?/pkgs/boost-?.$parameter{'bioConda'}{'boost'}.$parameter{'bioCondaBoostPatch'}.q?/lib\n?;
	print $FILEHANDLE q?fi\n\n' >> ~/.bashrc?;
	say $FILEHANDLE "\n";
    }

    ## Move to miniconda environment
    print $FILEHANDLE q?cd ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'};
    say $FILEHANDLE "\n";

    ## Download
    say $FILEHANDLE "## Download FindTranslocations";
    print $FILEHANDLE "wget --quiet https://github.com/J35P312/FindTranslocations/archive/version_".${$parameterHashRef}{'FindTranslocations'}.".zip ";
    print $FILEHANDLE "-O FindTranslocations-".${$parameterHashRef}{'FindTranslocations'}.".zip";  #Download outfile
    say $FILEHANDLE "\n";

    ## Extract
    say $FILEHANDLE "## Extract";
    print $FILEHANDLE "rm -rf FindTranslocations";
    say $FILEHANDLE "\n";
    print $FILEHANDLE "unzip FindTranslocations-".${$parameterHashRef}{'FindTranslocations'}.".zip ";
    say $FILEHANDLE "\n";
    print $FILEHANDLE "mv FindTranslocations-version_".${$parameterHashRef}{'FindTranslocations'}." ";
    print $FILEHANDLE "FindTranslocations ";
    say $FILEHANDLE "\n";

    ## Move to FindTranslocations directory
    say $FILEHANDLE "## Move to FindTranslocations directory";
    print $FILEHANDLE "cd FindTranslocations";
    say $FILEHANDLE "\n";
    print $FILEHANDLE "mkdir -p build";
    say $FILEHANDLE "\n";
    print $FILEHANDLE "cd build";
    say $FILEHANDLE "\n";

    ## Configure
    say $FILEHANDLE "## Configure";
    print $FILEHANDLE "cmake .. -DBoost_NO_BOOST_CMAKE=ON";
    say $FILEHANDLE "\n";

    print $FILEHANDLE "make";
    say $FILEHANDLE "\n";

    print $FILEHANDLE "cd ../bin";
    say $FILEHANDLE "\n";
    print $FILEHANDLE "chmod a+x FindTranslocations";
    say $FILEHANDLE "\n";

    ## Make available from conda environment
    my $cwd = cwd();
    say $FILEHANDLE "## Make available from conda environment";
    print $FILEHANDLE "ln -f -s  ";
    print $FILEHANDLE $parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/FindTranslocations/bin/FindTranslocations ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/?;
    say $FILEHANDLE "\n";    

    ## Clean-up
    print $FILEHANDLE q?cd ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'};
    say $FILEHANDLE "\n";
    print $FILEHANDLE "rm -rf FindTranslocations-".${$parameterHashRef}{'FindTranslocations'}.".zip";
    say $FILEHANDLE "\n";

    &DeactivateCondaEnvironment($FILEHANDLE);
}


sub MIPScripts {

    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];

    my $pwd = cwd();
    my $mipDirectory = $Bin;  #Alias

    ## Define MIP scripts and yaml files
    my @mipScripts = ("calculateAF.pl",
		      "maxAF.pl",
		      "mip.pl",
		      "qcCollect.pl",
		      "vcfParser.pl",
	);
    my %mipSubScripts;
    $mipSubScripts{"definitions"} = ["defineParameters.yaml"];
    $mipSubScripts{"t"} = ["test.t"];
    $mipSubScripts{"templates"} = ["mip_config.yaml"];

    if (&CheckCondaBinFileExists($parameterHashRef, "mip.pl")) {  #Proxy for all 

	return
    }

    ## Install MIPScripts
    say $FILEHANDLE "### Install MIPScripts";
    
    ## Create directories
    say $FILEHANDLE "## Create directories";
    foreach my $directory (keys %mipSubScripts) {

	print $FILEHANDLE "mkdir -p ";
	print $FILEHANDLE catdir($parameter{'condaPath'}, "envs", ${$parameterHashRef}{'condaEnvironment'}, "bin", $directory);
	say $FILEHANDLE "\n";
    }

    ## Copy mip scripts and sub scripts to conda env and make executable
    say $FILEHANDLE "## Copy mip scripts and subdirectory scripts to conda env and make executable\n";
    foreach my $script (@mipScripts) {
	    
	print $FILEHANDLE "cp ";
	print $FILEHANDLE catfile($Bin, $script)." ";
	say $FILEHANDLE catdir($parameter{'condaPath'}, "envs", ${$parameterHashRef}{'condaEnvironment'}, "bin");
	print $FILEHANDLE "chmod a+x ".catfile($parameter{'condaPath'}, "envs", ${$parameterHashRef}{'condaEnvironment'}, "bin", $script);
	say $FILEHANDLE "\n";
    }

    foreach my $directory (keys %mipSubScripts) {

	foreach my $script (@{$mipSubScripts{$directory}}) {

	    print $FILEHANDLE "cp ";
	    print $FILEHANDLE catfile($Bin, $directory, $script)." ";
	    say $FILEHANDLE catdir($parameter{'condaPath'}, "envs", ${$parameterHashRef}{'condaEnvironment'}, "bin", $directory);
	    print $FILEHANDLE "chmod a+x ".catfile($parameter{'condaPath'}, "envs", ${$parameterHashRef}{'condaEnvironment'}, "bin", $directory, $script);
	    say $FILEHANDLE "\n";
	}
    }
}


sub ActivateCondaEnvironment {

    my $parameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];

    ## Activate conda environment and install cpanm and MIP modules
    say $FILEHANDLE "## Activate conda environment";
    print $FILEHANDLE "source activate ".${$parameterHashRef}{'condaEnvironment'}." ";
    say $FILEHANDLE "\n";
}


sub DeactivateCondaEnvironment {

    my $FILEHANDLE = $_[0];

    ## Deactivate conda environment
    say $FILEHANDLE "## Deactivate conda environment";
    print $FILEHANDLE "source deactivate ";
    say $FILEHANDLE "\n";
}


sub CleanUpModuleInstall {

    my $FILEHANDLE = $_[0];
    my $pwd = $_[1];

    ## Moving up
    say $FILEHANDLE "## Moving back to original working directory";
    print $FILEHANDLE "cd ".$pwd;  #Go back to subroutine origin
    say $FILEHANDLE "\n";

    ## Clean up
    say $FILEHANDLE "## Clean up";
    print $FILEHANDLE "rm -rf .MIP";
    say $FILEHANDLE "\n";
}

sub CreateInstallDirectory {

    my $FILEHANDLE = $_[0];

    ## Create temp install directory
    say $FILEHANDLE "## Create temp install directory";
    say $FILEHANDLE "mkdir -p .MIP ";
    print $FILEHANDLE "cd .MIP";
    say $FILEHANDLE "\n";
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
	    
	    say STDERR q?Found ?.$programName.q? version ?.$programVersion.q? in miniconda directory: ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/?;
	    
	    if (${$parameterHashRef}{'update'} == 0) {

		say STDERR q?Skipping writting installation process for ?.$programName.q? ?.$programVersion;  
		return 1;
	    }
	    say STDERR "Writting install instructions for ".$programName;
	}   
	else {

	    say STDERR q?Found ?.$programName.q? in miniconda directory: ?.$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/?;
	    
	    if (${$parameterHashRef}{'update'} == 0) {
		
		say STDERR q?Skipping writting installation process for ?.$programName;  	    
		return 1;
	    }
	    say STDERR "Writting install instructions for ".$programName;
	}
	return 0;
    }
    else {
	
	say STDERR "Writting install instructions for ".$programName;
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
    say $FILEHANDLE "## Add softlink";
    say $FILEHANDLE "cd ".$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/?;

    print $FILEHANDLE "ln -f -s ";
    print $FILEHANDLE $binary.q? ?.$softLink;
    say $FILEHANDLE "\n";

    ## Move to back
    say $FILEHANDLE "## Move to original working directory";
    print $FILEHANDLE "cd ".$pwd;
    say $FILEHANDLE "\n";
}


sub EnableExecutable {

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $parameterHashRef = ${$argHashRef}{'parameterHashRef'};
    my $FILEHANDLE = ${$argHashRef}{'FILEHANDLE'};
    my $binary = ${$argHashRef}{'binary'};
    
    my $pwd = cwd();

    ## Add softlink
    say $FILEHANDLE "## Enable executable";
    say $FILEHANDLE "cd ".$parameter{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/?;

    print $FILEHANDLE "chmod a+x ";
    print $FILEHANDLE $binary.q? ?;
    say $FILEHANDLE "\n";

    ## Move to back
    say $FILEHANDLE "## Move to original working directory";
    print $FILEHANDLE "cd ".$pwd;
    say $FILEHANDLE "\n";
}


sub CheckMTCodonTable {

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $parameterHashRef = ${$argHashRef}{'parameterHashRef'};
    my $FILEHANDLE = ${$argHashRef}{'FILEHANDLE'};
    my $shareDirectory = ${$argHashRef}{'shareDirectory'};
    my $binary = ${$argHashRef}{'binary'};
    
    my $pwd = cwd();

    my $detectRegExp = q?perl -nae 'if($_=~/?.${$parameterHashRef}{'snpEffGenomeVersion'}.q?.MT.codonTable/) {print 1}' ?;
    my $addRegExp = q?perl -nae 'if($_=~/?.${$parameterHashRef}{'snpEffGenomeVersion'}.q?.reference/) {print $_; print "?.${$parameterHashRef}{'snpEffGenomeVersion'}.q?.MT.codonTable : Vertebrate_Mitochondrial\n"} else {print $_;}' ?;
    my $ret;

    if (-f $shareDirectory."/".$binary) {

	$ret = `$detectRegExp $shareDirectory/$binary`;
    }
    if (!$ret) {  #No MT.codonTable in config

	say $FILEHANDLE q?## Adding ?.${$parameterHashRef}{'snpEffGenomeVersion'}.q?.MT.codonTable : Vertebrate_Mitochondrial to ?.$shareDirectory.$binary;

	## Add MT.codon Table to config
	say $FILEHANDLE $addRegExp." ".$shareDirectory.$binary." > ".$shareDirectory.$binary.".tmp";
	print $FILEHANDLE "mv ".$shareDirectory.$binary.".tmp ".$shareDirectory.$binary;
	say $FILEHANDLE "\n";
	
    }
    else {

	say STDERR  "Found MT.codonTable in ".$shareDirectory."snpEff.config. Skipping addition to snpEff config"; 
    }
}


sub SnpEffDownload {
    
    my ($argHashRef) = @_;
    
    ## Flatten argument(s)
    my $parameterHashRef = ${$argHashRef}{'parameterHashRef'};
    my $FILEHANDLE = ${$argHashRef}{'FILEHANDLE'};
    
    print $FILEHANDLE "java -Xmx2g ";
    print $FILEHANDLE q?-jar ?.${$parameterHashRef}{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/snpEff.jar ?;
    print $FILEHANDLE "download ";
    print $FILEHANDLE " -v ";
    print $FILEHANDLE ${$parameterHashRef}{'snpEffGenomeVersion'}." ";
    print $FILEHANDLE q?-c ?.${$parameterHashRef}{'condaPath'}.q?/envs/?.${$parameterHashRef}{'condaEnvironment'}.q?/bin/snpEff.config ?;
    say $FILEHANDLE "\n";
    
}
