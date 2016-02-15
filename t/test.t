#!/usr/bin/env perl


###Will test perl modules and some selected funtions as well as vcf keys both in header and body. Adjusts dynamically according to supplied config file.
###Copyright 2016 Henrik Stranneheim

use strict;
use warnings;

BEGIN {


    my @modules = ("YAML",
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
	);	

    ## Evaluate that all modules required are installed
    &EvalModules(\@modules);
    
    sub EvalModules {
	
	##EvalModules
	
	##Function : Evaluate that all modules required are installed 
	##Returns  : ""
	##Arguments: $modulesArrayRef
	##         : $modulesArrayRef => Array of module names
	
	my $modulesArrayRef = $_[0];
	
	foreach my $module (@{$modulesArrayRef}) {

	    $module =~s/::/\//g;  #Replace "::" with "/" since the automatic replacement magic only occurs for barewords.
	    $module .= ".pm";  #Add perl module ending for the same reason
	    
	    eval { 
		
		require $module; 
	    };
	    if($@) {
		
		warn("NOTE: ".$module." not installed - Please install to run MIP.\n");
		warn("NOTE: Aborting!\n");
		exit 1;
	    }
	}
    }
}

use Test::More;
use Getopt::Long;
use FindBin qw($Bin); #Find directory of script
use File::Basename qw(dirname);
use File::Spec::Functions qw(catdir);

##Cpan
use YAML;

use vars qw($USAGE);

BEGIN {
    $USAGE =
	qq{test.pl infile.vcf [VCF] configFile [YAML]
           -h/--help Display this help message   
           -v/--version Display version
        };    
}

my ($infile, $configFile);
my (%scriptParameter, %vcfParserData, %consequenceSeverity, %metaData, %siftTerm, %polyPhenTerm);

my $testVersion = "1.0.0";

if(scalar(@ARGV) == 0) {

    print STDOUT $USAGE, "\n";
    exit;
}


############
####MAIN####
############

$infile = $ARGV[0];
$configFile = $ARGV[1];

###User Options
GetOptions('h|help' => sub { print STDOUT $USAGE, "\n"; exit;},  #Display help text
	   'v|version' => sub { print STDOUT "\nvcfParser.pl ".$testVersion, "\n\n"; exit;},  #Display version number
    );

unless (defined($infile)) {
    
    print STDERR "Please supply an infile\n";
    exit;
}

unless (defined($configFile)) {
    
    print STDERR "Please supply a configFile\n";
    exit;
}

## Test perl modules and functions
&TestModules();

if (defined($configFile)) {  #Input from cmd

    ## Loads a YAML file into an arbitrary hash and returns it.
    %scriptParameter = &LoadYAML($configFile);  #Load parameters from configfile
}


if($infile =~/.selected.vcf/) {

    ## Reads a file containg features to be annotated using range queries
    &ReadRangeFile({'vcfParserDataHashRef' => \%vcfParserData,
		    'rangeCoulumnsArrayRef' => \@{$scriptParameter{'vcfParserSelectFeatureAnnotationColumns'}},
		    'infilePath' => $scriptParameter{'referencesDir'}."/".$scriptParameter{'vcfParserSelectFile'},
		    'rangeFileKey' => "SelectFile",
		   });
}
else {  #Range file

    ## Reads a file containg features to be annotated using range queries
    &ReadRangeFile({'vcfParserDataHashRef' => \%vcfParserData,
		    'rangeCoulumnsArrayRef' => \@{$scriptParameter{'vcfParserRangeFeatureAnnotationColumns'}},
		    'infilePath' => $scriptParameter{'referencesDir'}."/".$scriptParameter{'vcfParserRangeFeatureFile'},
		    'rangeFileKey' => "RangeFile",
		   });
}

## Reads infile in vcf format and parses annotations
&ReadInfileVCF({'scriptParameterHashRef' => \%scriptParameter,
		'vcfParserDataHashRef' => \%vcfParserData,
		'metaData' => \%metaData,
		'consequenceSeverityHashRef' => \%consequenceSeverity,
		'siftTermHashRef' => \%siftTerm,
		'polyPhenTermHashRef' => \%polyPhenTerm,
		'infile' => $infile,
	       });

done_testing();   # reached the end safely


######################
####SubRoutines#######
######################


sub TestModules {

##TestModules
    
##Function : Test perl modules and functions 
##Returns  : ""
##Arguments: 
##         : 

    print STDOUT "\nTesting perl modules and selected functions\n\n";

    use FindBin qw($Bin); #Find directory of script

    ok(defined($Bin),"FindBin: Locate directory of script");

    use File::Basename qw(dirname);  #Strip the last part of directory

    ok(dirname($Bin), "File::Basename qw(dirname): Strip the last part of directory");
 
    use File::Spec::Functions qw(catdir);

    ok(catdir(dirname($Bin), "t"),"File::Spec::Functions qw(catdir): Concatenate directories");
    
    use YAML;
    
    my $yamlFile = catdir(dirname($Bin), "definitions", "defineParameters.yaml");
    ok( -f $yamlFile,"YAML: File= $yamlFile in MIP directory");
    
    my $yaml = YAML::LoadFile($yamlFile);  #Create an object
    ok( defined $yaml,"YAML: Load File" );  #Check that we got something
    ok(Dump( $yaml ),"YAML: Dump file");
    
    use Log::Log4perl;
    ## Creates log
    my $logFile = catdir(dirname($Bin), "templates", "mip_config.yaml");
    ok( -f $logFile,"Log::Log4perl: File= $logFile in MIP directory");

    ## Create log4perl config file
    my $config = &CreateLog4perlCongfig(\$logFile);
    
    ok(Log::Log4perl->init(\$config), "Log::Log4perl: Initate");
    ok(Log::Log4perl->get_logger("MIPLogger"), "Log::Log4perl: Get logger");
    
    my $logger = Log::Log4perl->get_logger("MIPLogger");
    ok($logger->info("1"), "Log::Log4perl: info");
    ok($logger->warn("1"), "Log::Log4perl: warn");
    ok($logger->error("1"), "Log::Log4perl: error");
    ok($logger->fatal("1"), "Log::Log4perl: fatal");
    
    use Getopt::Long;
    push(@ARGV, ("-verbose", "2"));
    my $verbose = 1;
    ok(GetOptions("verbose:n"  => \$verbose), "Getopt::Long: Get options call");
    ok ($verbose == 2, "Getopt::Long: Get options modified");
    
    use DateTime;
    my $dateTime = DateTime->now(time_zone=>'local');
    ok (defined($dateTime), "DateTime: Now");
    my $dateTimeStamp = $dateTime->datetime();
    ok (defined($dateTimeStamp), "DateTime: Time stamp");
    my $date = $dateTime->ymd('-');  #Catches current date
    ok (defined($date), "DateTime: Current date");
    
    use List::MoreUtils qw(any);
    
    my @array = ("apples", "oranges");
    ok( (any {$_ eq "apples"} @array),"List::MoreUtils: Any");
}

sub CreateLog4perlCongfig {

##CreateLog4perlCongfig
    
##Function : Create log4perl config file. 
##Returns  : "$config"
##Arguments: $fileName
##         : $fileName => log4perl config file {REF}

    my $fileNameRef = $_[0];

    my $conf = q?
        log4perl.category.MIPLogger = TRACE, LogFile, ScreenApp
        log4perl.appender.LogFile = Log::Log4perl::Appender::File
        log4perl.appender.LogFile.filename = ?.$$fileNameRef.q?
        log4perl.appender.LogFile.layout=PatternLayout
        log4perl.appender.LogFile.layout.ConversionPattern = [%p] %d %c - %m%n
        log4perl.appender.ScreenApp = Log::Log4perl::Appender::Screen
        log4perl.appender.ScreenApp.layout = PatternLayout
        log4perl.appender.ScreenApp.layout.ConversionPattern = [%p] %d %c - %m%n
        ?;
    return $conf;
}


sub ReadInfileVCF {

##ReadInfileVCF

##Function : Reads infile in vcf format and adds and parses annotations
##Returns  : ""
##Arguments: $scriptParameterHashRef, $vcfParserDataHashRef, $metaDataHashRef, $consequenceSeverityHashRef, $siftTermHashRef, $polyPhenTermHashRef, $infile
##         : $scriptParameterHashRef     => The active parameters for this analysis hash {REF}
##         : vcfParserDataHashRef        => The keys from vcfParser i.e. range file and select file
##         : $metaDataHashRef            => Vcf meta data {REF}
##         : $consequenceSeverityHashRef => Consequence severity for SO-terms {REF}
##         : $siftTermHashRef            => Sift prediction terms {REF}
##         : $polyPhenTermHashRef        => PolyPhen prediction terms {REF}
##         : $infile                     => Infile

    my ($argHashRef) = @_;
    
    ## Flatten argument(s)
    my $scriptParameterHashRef = ${$argHashRef}{'scriptParameterHashRef'};
    my $vcfParserDataHashRef = ${$argHashRef}{'vcfParserDataHashRef'};
    my $metaDataHashRef = ${$argHashRef}{'metaDataHashRef'};
    my $consequenceSeverityHashRef = ${$argHashRef}{'consequenceSeverityHashRef'};
    my $siftTermHashRef = ${$argHashRef}{'siftTermHashRef'};
    my $polyPhenTermHashRef = ${$argHashRef}{'polyPhenTermHashRef'};
    my $infile = ${$argHashRef}{'infile'};

    my @vepFormatField;
    my %vepFormatFieldColumn;

    my %vcfHeader;
    my %vcfInfoKey;

    print STDOUT "\nTesting VCF Header:\n\n";
    
    open(VCF, "<".${$argHashRef}{'infile'}) or die "Can't open ".${$argHashRef}{'infile'}.":".$!, "\n"; 

    while (<VCF>) {

	chomp $_;  # Remove newline
	if($. > 5000) {  #Exit after parsing X lines
	    last;
	} 
	if (m/^\s+$/) {	# Avoid blank lines
	    next;
	}
	if ($_=~/^##(\S+)=/) {  # MetaData

	    &ParseMetaData(\%{$metaDataHashRef}, $_);
	    
	    if ($_=~/INFO\=\<ID\=([^,]+)/) {
	
		$vcfHeader{'INFO'}{$1} = $1; #Save to hash
	    }
	    if ($_=~/INFO\=\<ID\=CSQ/) { #Find VEP INFO Field

		ok(${$scriptParameterHashRef}{'pVariantEffectPredictor'}, "VEP: CSQ in header and VEP should have been executed");

		if ($_=~/Format:\s(\S+)"\>/) { #Locate Format within VEP INFO meta line

		    @vepFormatField = split(/\|/, $1);  
		    
		    for (my $fieldCounter=0;$fieldCounter<scalar(@vepFormatField);$fieldCounter++) {
			
			$vepFormatFieldColumn{$vepFormatField[$fieldCounter]} = $fieldCounter; #Save the order of VEP features
		    }
		}
		next;
	    }
	    next;
	}
	if ($_=~/^#CHROM/) {

	    ### Check Header now that we read all

	    ## VT
	    if (${$scriptParameterHashRef}{'pVT'} > 0) {
		
		if (${$scriptParameterHashRef}{'VTDecompose'} > 0) {
       
		    ok( defined($vcfHeader{'INFO'}{'OLD_MULTIALLELIC'}), "VTDecompose key: OLD_MULTIALLELIC");
		}
		if (${$scriptParameterHashRef}{'VTNormalize'} > 0) {
		    
		    ok( defined($vcfHeader{'INFO'}{'OLD_VARIANT'}), "VTNormalize key: OLD_VARIANT");
		}
	    }

	    ## VCFParser
	    if (${$scriptParameterHashRef}{'pVCFParser'} > 0) {
		
		for my $key (keys %{$vcfParserDataHashRef}) {
		    
		    ok( defined($vcfHeader{'INFO'}{$key}), "VCFParser key: $key");
		}

		## Keys from vcfParser that are dynamically created from parsing the data
		my @vcfParserDynamicKeys = ("HGVScp",
					    "MostSevereConsequence",
					    "GeneticRegionAnnotation",
					    "Sift",
					    "PolyPhen",
		    );

		foreach my $key (@vcfParserDynamicKeys) {

		    ok( defined($vcfHeader{'INFO'}{$key}), "VCFParser dynamic keys: $key");
		}
	    }

	    ## SnpEff
	    if (${$scriptParameterHashRef}{'pSnpEff'} > 0) {
		
		my @splittedKeys;
		for my $annotationFile (keys %{${$scriptParameterHashRef}{'snpSiftAnnotationFiles'}}) {
		    
		    if (${$scriptParameterHashRef}{'snpSiftAnnotationOutInfoKey'}{$annotationFile}) {
			
			@splittedKeys = split(',', ${$scriptParameterHashRef}{'snpSiftAnnotationOutInfoKey'}{$annotationFile});
		    }
		    else {
			
			@splittedKeys = split(',', ${$scriptParameterHashRef}{'snpSiftAnnotationFiles'}{$annotationFile});
		    }
		    foreach my $key (@splittedKeys) {		  
			
			ok( defined($vcfHeader{'INFO'}{$key}), "SnpSift Annotation key: $key");
		    }
		}
		for my $key (@{${$scriptParameterHashRef}{'snpSiftDbNSFPAnnotations'}}) {		  
		    
		    $key =~ s/\+/_/g;  #Special case due to the fact that snpEff v4.2 transforms + to _ for some reason
		    $key = "dbNSFP_".$key;
		    ok( defined($vcfHeader{'INFO'}{$key}), "SnpSift dbNSFP_key: $key");
		}
	    }
	    ## GENMOD
	    ## 1000G key from genmodFilter
	    if ( (${$scriptParameterHashRef}{'pVT'} > 0) && (${$scriptParameterHashRef}{'VTgenmodFilter'} > 0) ) {

		if (defined(${$scriptParameterHashRef}{'VTgenmodFilter1000G'})) {
		    
		    ok( defined($vcfHeader{'INFO'}{'1000GAF'}), "GENMODFilter: 1000GAF key");
		}
	    }	    
	    if (${$scriptParameterHashRef}{'pRankVariants'} > 0) {

		## Keys from genmod
		my @genmodKeys = ("Compounds",
				  "RankScore",
				  "ModelScore",
				  "GeneticModels",
		    );

		foreach my $key (@genmodKeys) {

		    ok( defined($vcfHeader{'INFO'}{$key}), "GENMOD: $key");
		}

		## Spidex key from genmodAnnotate
		if (${$scriptParameterHashRef}{'spidexFile'}) {
		    
		    ok( defined($vcfHeader{'INFO'}{'SPIDEX'}), "GENMODAnnotate: SPIDEX key");
		}
		## CADD key from genmodAnnotate
		if ( (${$scriptParameterHashRef}{'cadd1000GenomesFile'}) || (${$scriptParameterHashRef}{'cadd1000GenomesFile'}) ) {
		    
		    ok( defined($vcfHeader{'INFO'}{'CADD'}), "GENMODAnnotate: CADD key");
		}
	    }
	    next;
	}
	if ( $_ =~/^(\S+)/ ) { 

	    my @lineElements = split("\t",$_); #Loads vcf elements
	    
	    my @keyValues = split(/;/, $lineElements[7]); #Split INFO field to key=value items

	    for my $element (@keyValues) {

		my @keys = split("=", $element);  #key = 0 and value = 1

		$vcfInfoKey{$keys[0]}++;  #Increment
	    }
	}
    }
    close(VCF);

    ##Check keys found in INFO field

    print STDOUT "\nTesting VCF INFO fields and presence in header:\n\n";

    foreach my $key (keys %vcfInfoKey) {

	ok( defined($vcfHeader{'INFO'}{$key}), "Found both header and line field key for: $key. KeyCount: ".$vcfInfoKey{$key});
    }
}


sub ReadRangeFile {

##ReadRangeFile

##Function : Reads a file containg features to be annotated using range queries e.g. EnsemblGeneID.
##Returns  : ""
##Arguments: $vcfParserDataHashRef, $rangeCoulumnsArrayRef, $rangeFileKey, $infilePath, $paddingRef, $selectFeatureMatchingColumn
##         : $vcfParserDataHashRef            => Range file hash {REF}
##         : $rangeCoulumnsArrayRef       => Range columns to include {REF}
##         : $rangeFileKey                => Range file key used to seperate range file(s) i.e., select and range
##         : $infilePath                  => Infile path
##         : $selectFeatureMatchingColumn => Column in the select file to match with vcf key annotation

    my ($argHashRef) = @_;
    
    ## Flatten argument(s)
    my $vcfParserDataHashRef = ${$argHashRef}{'vcfParserDataHashRef'};
    my $rangeCoulumnsArrayRef = ${$argHashRef}{'rangeCoulumnsArrayRef'};

    my @headers; #Save headers from rangeFile

    open(RRF, "<".${$argHashRef}{'infilePath'}) or die "Can't open ".${$argHashRef}{'infilePath'}.":".$!, "\n"; 

    while (<RRF>) {
	
	chomp $_; #Remove newline

	if (m/^\s+$/) {		# Avoid blank lines
	    
	    next;
	}
	if ($_=~/^##/) {#MetaData - Avoid

	    next;
	}
	if ($_=~/^#/) {#Header/Comment
	    
	    @headers = split(/\t/, $_);

	    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar(@{$rangeCoulumnsArrayRef});$extractColumnsCounter++) { #Defines what scalar to store
	
		my $headerKeyRef = \$headers[ $$rangeCoulumnsArrayRef[$extractColumnsCounter] ];
		${$vcfParserDataHashRef}{$$headerKeyRef} = $extractColumnsCounter; #Column position in supplied range input file
	    }
	    next;
	}
    }
    close(RRF);
    print STDOUT "\nFinished Reading ".${$argHashRef}{'rangeFileKey'}." file: ".${$argHashRef}{'infilePath'},"\n";
}


sub ParseMetaData {

##ParseMetaData
    
##Function : Writes metadata to filehandle specified by order in metaDataOrders.
##Returns  : ""
##Arguments: $metaDataHashRef, $metaDataString
##         : $metaDataHashRef => Hash for metaData {REF}
##         : $metaDataString  => The metaData string from vcf header
    
    my $metaDataHashRef = $_[0];
    my $metaDataString = $_[1];
    
    if ($metaDataString=~/^##fileformat/) {  #Catch fileformat as it has to be at the top of header
	
	push(@{${$metaDataHashRef}{'fileformat'}{'fileformat'}}, $metaDataString);  #Save metadata string
    }
    elsif ($metaDataString=~/^##contig/) {  #catch contigs to not sort them later

	push(@{${$metaDataHashRef}{'contig'}{'contig'}}, $metaDataString);  #Save metadata string
    }
    elsif ($metaDataString=~/^##(\w+)=(\S+)/) {  #FILTER, FORMAT, INFO etc and more custom records
	
	push(@{${$metaDataHashRef}{$1}{$2}}, $metaDataString);  #Save metadata string
    }
    else {  #All oddities
	
	push(@{${$metaDataHashRef}{'other'}{'other'}}, $metaDataString);  #Save metadata string
    }
}


sub LoadYAML {
 
##LoadYAML
    
##Function : Loads a YAML file into an arbitrary hash and returns it. Note: Currently only supports hashreferences and hashes and no mixed entries.
##Returns  : %yamlHash
##Arguments: $yamlFile
##         : $yamlFile => The yaml file to load

    my $yamlFile = $_[0];

    my %yamlHash;

    open (my $YAML, "<", $yamlFile) or die "can't open ".$yamlFile.":".$!, "\n";  #Log4perl not initialised yet, hence no logdie

    local $YAML::QuoteNumericStrings=1;  #Force numeric values to strings in YAML representation
    if (defined ($YAML::QuoteNumericStrings) && $YAML::QuoteNumericStrings == 1) {

	%yamlHash = %{ YAML::LoadFile($yamlFile) };  #Load hashreference as hash
    }
    close($YAML);

    return %yamlHash;
}
