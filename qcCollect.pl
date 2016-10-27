#!/usr/bin/env perl

use Modern::Perl '2014';
use warnings qw( FATAL utf8 );
use autodie;
use v5.18;  #Require at least perl 5.18
use utf8;  #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );

##Collects MPS QC from MIP. Loads information on files to examine and values to extract from in YAML format and outputs exracted metrics in YAML format. 
#Copyright 2013 Henrik Stranneheim 

use Pod::Usage;
use Pod::Text;
use Getopt::Long;
use POSIX;
use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case

## Third party module(s)
use YAML;

our $USAGE;

BEGIN {
    $USAGE =
        qq{qcCollect.pl -si [sampleinfo.yaml] -r [regexp.yaml] -o [outfile]
               -si/--sampleInfoFile SampleInfo file (YAML format, Supply whole path, mandatory)
               -r/--regExpFile Regular expression file (YAML format, Supply whole path, mandatory)
               -o/--outfile The data file output (Supply whole path, defaults to "qcmetrics.yaml")
               -preg/--printRegExp Print the regexp used at CMMS switch (Defaults to "0" (=no))
               -prego/--printRegExpOutFile RegExp YAML outfile (Defaults to "qc_regExp.yaml")
               -h/--help Display this help message
               -v/--version Display version};
}
my ($sampleInfoFile, $regExpFile, $printRegExp);
my ($outfile, $printRegExpOutFile) = ("qcmetrics.yaml", "qc_regExp.yaml");
my (%qcData, %evaluateMetric);
my %qcHeader; #Save header(s) in each outfile
my %qcProgramData; #Save data in each outFile

my $qcCollectVersion = "2.0.0";

GetOptions('si|sampleInfoFile:s' => \$sampleInfoFile,
	   'r|regExpFile:s' => \$regExpFile,
	   'o|outfile:s'  => \$outfile, 
	   'preg|printRegExp:n' => \$printRegExp,
	   'prego|printRegExpOutFile:s' => \$printRegExpOutFile,
	   'h|help' => sub { say STDOUT $USAGE; exit;},  #Display help text
	   'v|version' => sub { say STDOUT "\nqcCollect.pl ".$qcCollectVersion, "\n"; exit;},  #Display version number
    )  or &Help({USAGE => $USAGE,
		 exitCode => 1,
		});

if ($printRegExp) {

    ## Write default regExp to YAML
    &RegExpToYAML({printRegExpOutFile => $printRegExpOutFile,
		  });
    print STDOUT "Wrote RegExp YAML file to: ".$printRegExpOutFile, "\n";
    exit;
}

if (! $sampleInfoFile) {

    print STDERR "\n";
    print STDOUT $USAGE, "\n";
    print STDERR "Must supply a '-sampleInfoFile' (supply whole path)", "\n\n";
    exit;
}
if (! $regExpFile) {

    print STDERR "\n";
    print STDOUT $USAGE, "\n";
    print STDERR "Must supply a '-regExpFile' (supply whole path)", "\n\n";
    exit;
}

####MAIN

## Loads a YAML file into an arbitrary hash and returns it
my %sampleInfo = &LoadYAML({yamlFile => $sampleInfoFile,
			       });

## Loads a YAML file into an arbitrary hash and returns it
my %regExp = &LoadYAML({yamlFile => $regExpFile,
			   });


## Extracts all qcdata on sampleID level using information in %sampleInfo and %regExp
&SampleQC({sampleInfoHashRef => \%sampleInfo,
	   regExpHashRef => \%regExp,
	   qcDataHashRef => \%qcData,
	   qcHeaderHashRef => \%qcHeader,
	   qcProgramDataHashRef => \%qcProgramData,
	  });


## Extracts all qcdata on family level using information in %sampleInfoFile and %regExp
&FamilyQC({sampleInfoHashRef => \%sampleInfo,
	   regExpHashRef => \%regExp,
	   qcDataHashRef => \%qcData,
	   qcHeaderHashRef => \%qcHeader,
	   qcProgramDataHashRef => \%qcProgramData,
	  }); 

##Add qcCollect version to qcData yaml file
$qcData{Program}{QCCollect}{Version} = $qcCollectVersion;
$qcData{Program}{QCCollect}{RegExpFile} = $regExpFile;

foreach my $sampleID (keys %{$sampleInfo{sample}}) {

    ## Defines programs, etrics and thresholds to evaluate
    &DefineEvaluateMetric({sampleID => $sampleID,
			  });
}


## Evaluate the metrics
&EvaluateQCParameters({qcDataHashRef => \%qcData,
		       evaluateMetricHashRef => \%evaluateMetric,
		      });

## Writes a YAML hash to file
&WriteYAML({yamlHashRef => \%qcData,
	    yamlFilePathRef => \$outfile,
	   });

####SubRoutines

sub FamilyQC {

##FamilyQC

##Function : Extracts all qcdata on family level using information in %sampleInfoFile and %regExp
##Returns  : ""
##Arguments: $sampleInfoHashRef, $regExpHashRef, $qcDataHashRef, $qcHeaderHashRef, $qcProgramDataHashRef
##         : $sampleInfoHashRef    => Info on samples and family hash {REF}
##         : $regExpHashRef        => RegExp hash {REF}
##         : $qcDataHashRef        => QCData hash {REF}
##         : $qcHeaderHashRef      => Save header(s) in each outfile {REF}
##         : $qcProgramDataHashRef => Hash to save data in each outFile {REF}
    
    my ($argHashRef) = @_;
    
    ## Flatten argument(s)
    my $sampleInfoHashRef;
    my $regExpHashRef;
    my $qcDataHashRef;
    my $qcHeaderHashRef;
    my $qcProgramDataHashRef;
    
    my $tmpl = { 
	sampleInfoHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$sampleInfoHashRef},
	regExpHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$regExpHashRef},
	qcDataHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$qcDataHashRef},
	qcHeaderHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$qcHeaderHashRef},
	qcProgramDataHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$qcProgramDataHashRef},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    for my $program ( keys %{ ${$sampleInfoHashRef}{Program} } ) { #For every program

	my $outDirectory;
	my $outFile;	
	
	if (${$sampleInfoHashRef}{Program}{$program}{Version} ) {
	    
	    ${$qcDataHashRef}{Program}{$program}{Version} = ${$sampleInfoHashRef}{Program}{$program}{Version}; #Add version to qcData
	}
	if (${$sampleInfoHashRef}{Program}{$program}{OutDirectory} ) {
	    
	    $outDirectory = ${$sampleInfoHashRef}{Program}{$program}{OutDirectory}; #Extract OutDirectory
	}
	if (${$sampleInfoHashRef}{Program}{$program}{OutFile} ) {
	    
	    $outFile = ${$sampleInfoHashRef}{Program}{$program}{OutFile}; #Extract OutFile
	}  
	
	## Parses the RegExpHash structure to identify if the info is 1) Paragraf section(s) (both header and data line(s)); 2) Seperate data line.
	&ParseRegExpHashAndCollect({regExpHashRef => $regExpHashRef,
				    qcProgramDataHashRef => $qcProgramDataHashRef,
				    qcHeaderHashRef => $qcHeaderHashRef,
				    program => $program,
				    outDirectory => $outDirectory,
				    outFile => $outFile,
				   });
	## Add extracted information to qcData
	&AddToqcData({sampleInfoHashRef => $sampleInfoHashRef,
		      regExpHashRef => $regExpHashRef,
		      qcDataHashRef => $qcDataHashRef,
		      qcHeaderHashRef => $qcHeaderHashRef,
		      qcProgramDataHashRef => $qcProgramDataHashRef,
		      program => $program,
		     });
    }
}

sub SampleQC {

##SampleQC

##Function : Collects all sample qc in files defined by sampleInfoFile and regular expressions defined by regExp.
##Returns  : ""
##Arguments: $sampleInfoHashRef, $regExpHashRef, $qcDataHashRef, $qcHeaderHashRef, $qcProgramDataHashRef
##         : $sampleInfoHashRef    => Info on samples and family hash {REF}
##         : $regExpHashRef        => RegExp hash {REF}
##         : $qcDataHashRef        => QCData hash {REF}
##         : $qcHeaderHashRef      => Save header(s) in each outfile {REF}
##         : $qcProgramDataHashRef => Hash to save data in each outFile {REF}
    
    my ($argHashRef) = @_;
    
    ## Flatten argument(s)
    my $sampleInfoHashRef;
    my $regExpHashRef;
    my $qcDataHashRef;
    my $qcHeaderHashRef;
    my $qcProgramDataHashRef;
    
    my $tmpl = { 
	sampleInfoHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$sampleInfoHashRef},
	regExpHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$regExpHashRef},
	qcDataHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$qcDataHashRef},
	qcHeaderHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$qcHeaderHashRef},
	qcProgramDataHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$qcProgramDataHashRef},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];
    
    for my $sampleID ( keys %{ ${$sampleInfoHashRef}{sample} } ) { #For every sample id
	
	for my $program ( keys %{ ${$sampleInfoHashRef}{sample}{$sampleID}{Program} } ) { #For every program  
	    
	    for my $infile ( keys %{ ${$sampleInfoHashRef}{sample}{$sampleID}{Program}{$program} } ) { #For every infile
		
		my $outDirectory = ${$sampleInfoHashRef}{sample}{$sampleID}{Program}{$program}{$infile}{OutDirectory};
		my $outFile = ${$sampleInfoHashRef}{sample}{$sampleID}{Program}{$program}{$infile}{OutFile};
		
		## Parses the RegExpHash structure to identify if the info is 1) Paragraf section(s) (both header and data line(s)); 2) Seperate data line.
		&ParseRegExpHashAndCollect({regExpHashRef => $regExpHashRef,
					    qcProgramDataHashRef => $qcProgramDataHashRef,
					    qcHeaderHashRef => $qcHeaderHashRef,
					    program => $program,
					    outDirectory => $outDirectory,
					    outFile => $outFile,
					   });

		## Add extracted information to qcData
		&AddToqcData({sampleInfoHashRef => $sampleInfoHashRef,
			      regExpHashRef => $regExpHashRef,
			      qcDataHashRef => $qcDataHashRef,
			      qcHeaderHashRef => $qcHeaderHashRef,
			      qcProgramDataHashRef => $qcProgramDataHashRef,
			      sampleID => $sampleID,
			      program => $program,
			      infile => $infile,
			     });
	    }   
	}
    }
}

sub ParseRegExpHashAndCollect {

##ParseRegExpHashAndCollect

##Function  : Parses the RegExpHash structure to identify if the info is 1) Paragraf section(s) (both header and data line(s)); 2) Seperate data line.
##Returns   : ""
##Arguments : $regExpHashRef, qcProgramDataHashRef, $qcHeaderHashRef, $program, $outDirectory, $outFile
##          : $regExpHashRef        => RegExp hash {REF}
##          : $qcHeaderHashRef      => Save header(s) in each outfile {REF}
##          : $qcProgramDataHashRef => Hash to save data in each outFile {REF}
##          : $program              => The program to examine
##          : $outDirectory         => Programs out directory
##          : $outFile              => Programs out file containing parameter to evaluate

     my ($argHashRef) = @_;

     ## Flatten argument(s)
     my $regExpHashRef;
     my $qcHeaderHashRef;
     my $qcProgramDataHashRef;
     my $program;
     my $outDirectory;
     my $outFile;
     
     my $tmpl = { 
	 regExpHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$regExpHashRef},
	 qcHeaderHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$qcHeaderHashRef},
	 qcProgramDataHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$qcProgramDataHashRef},
	 program => { required => 1, defined => 1, strict_type => 1, store => \$program},
	 outDirectory => { strict_type => 1, store => \$outDirectory},
	 outFile => { strict_type => 1, store => \$outFile},
     };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];
    
    my $regExp; #Holds the current regExp
    my @separators = ('\s+','!',','); #Covers both whitespace and tab. Add other separators if required
    
    for my $regExpKey ( keys %{ ${$regExpHashRef}{$program} } ) { #Find the actual regular expression(s) for each program that is used
	
	if ($regExpKey =~/^header|header$/i) { #Detect if the outfile contains paragrafs/header info in the outfile i.e. data is formated as a paragraf with header(s) and line(s). "regExpKey" should either start with or end with "header". This section extracts the header line(s) for the entire outdata file. Necessary to assign correct data entry to header entry later (headers and data are saved in seperate hashes). 
	    
	    ##Format outfile: Paragraf section
	    for my $regExpHeaderKey ( keys %{ ${$regExpHashRef}{$program}{$regExpKey} } ) { #Paragraf Header 
		
		$regExp = ${$regExpHashRef}{$program}{$regExpKey}{$regExpHeaderKey}; #The regular expression used to collect paragraf header info
		
		for (my $separatorElement=0;$separatorElement<scalar(@separators);$separatorElement++) { #Loop through possible separators to seperate any eventual header elements
                    
		    if ($regExpHeaderKey =~/^header|header$/i) { #Detect if the regExp key is a paragraf header and not paragraf data (Header line and data line(s))
                        
			 @{ ${$qcHeaderHashRef}{$program}{$regExpKey}{$regExpHeaderKey} } = split(/$separators[$separatorElement]/, `$regExp $outDirectory/$outFile`); #Collect paragraf header                           
			
			if ( defined(${$qcHeaderHashRef}{$program}{$regExpKey}{$regExpHeaderKey})) { #Then split should have been successful                                                                                          
			    last; #Found correct separator - do not continue                          
			}   
		    }
		    else { #For paragraf data line(s)

			@{ ${$qcProgramDataHashRef}{$program}{$regExpKey}{$regExpHeaderKey} } = split(/$separators[$separatorElement]/, `$regExp $outDirectory/$outFile`); #Collect paragraf data

			if ( defined(${$qcProgramDataHashRef}{$program}{$regExpKey}{$regExpHeaderKey}[1])) { #Then split should have been successful                                                                                                                           
			    last; #Found correct separator - do not continue                                                                                                                                      
			}
		    }
		}
	    }
	}
	else { #For info contained in Entry --> Value i.e. same line.

	    $regExp = ${$regExpHashRef}{$program}{$regExpKey}; #The regular expression used to collect info
	    
	    for (my $separatorElement=0;$separatorElement<scalar(@separators);$separatorElement++) { #Loop through possible separators 
		
		@{ ${$qcProgramDataHashRef}{$program}{$regExpKey} } = split(/$separators[$separatorElement]/, `$regExp $outDirectory/$outFile`); #Collect data. Use regExpKey as element header

		if ( defined(${$qcProgramDataHashRef}{$program}{$regExpKey}[1])) { #Then split should have been successful

		    last; #Found correct separator do not continue
		}
	    }
	}
    }
}

sub AddToqcData {

##AddToqcData

##Function  : Add to qcData hash to enable write to yaml format
##Returns   : ""
##Arguments : $sampleInfoHashRef, $regExpHashRef, $qcDataHashRef, $qcHeaderHashRef, $qcProgramDataHashRef, $sampleID, $program, $inFile
##          : $sampleInfoHashRef    => Info on samples and family hash {REF}
##          : $regExpHashRef        => RegExp hash {REF}
##          : $qcDataHashRef        => QCData hash {REF}
##          : $qcHeaderHashRef      => Save header(s) in each outfile {REF}
##          : $qcProgramDataHashRef => Hash to save data in each outFile {REF}
##          : $sampleID             => SampleID
##          : $program              => The program to examine 
##          : $inFile               => infile to program

    my ($argHashRef) = @_;
    
    ## Flatten argument(s)
    my $sampleInfoHashRef;
    my $regExpHashRef;
    my $qcDataHashRef;
    my $qcHeaderHashRef;
    my $qcProgramDataHashRef;
    my $sampleID;
    my $program;
    my $infile;

    my $tmpl = { 
	sampleInfoHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$sampleInfoHashRef},
	regExpHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$regExpHashRef},
	qcDataHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$qcDataHashRef},
	qcHeaderHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$qcHeaderHashRef},
	qcProgramDataHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$qcProgramDataHashRef},
	sampleID => { strict_type => 1, store => \$sampleID},
	program => { required => 1, defined => 1, strict_type => 1, store => \$program},
	infile => { strict_type => 1, store => \$infile},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];
    
    for my $regExpKey ( keys %{ ${$regExpHashRef}{$program} } ) { #All regExp per program 
	
	if ($regExpKey !~/^header|header$/i) { #For info contained in Entry --> Value i.e. same line 

	    if (scalar(@{ ${$qcProgramDataHashRef}{$program}{$regExpKey} }) == 1) { #Enable seperation of writing array or key-->value in qcData
		
		if ( ($sampleID) && ($infile) ) {
		    
		    ${$qcDataHashRef}{sample}{$sampleID}{$infile}{$program}{$regExpKey} = ${$qcProgramDataHashRef}{$program}{$regExpKey}[0]; #key-->value for sampleID
		}
		else {  #Family level
		    
		    ${$qcDataHashRef}{Program}{$program}{$regExpKey} = ${$qcProgramDataHashRef}{$program}{$regExpKey}[0]; #key-->value for familyID
		}
		if ($program eq "ChanjoSexCheck") {#Check gender for sampleID
		    
		    my $ChanjoSexCheck = @{${$qcProgramDataHashRef}{$program}{$regExpKey}}[0]; #ArrayRef

		    ## Check that assumed gender is supported by coverage on chrX and chrY
		    &GenderCheck({sampleInfoHashRef => $sampleInfoHashRef,
				  qcDataHashRef => $qcDataHashRef,
				  sampleIDRef => \$sampleID,
				  infileRef => \$infile,
				  chanjoSexCheckGenderRef => \$ChanjoSexCheck,
				 });
		}
	    }
	    else { #Write array to qcData

		for (my $regExpKeyCounter=0;$regExpKeyCounter<scalar(@{ ${$qcProgramDataHashRef}{$program}{$regExpKey} });$regExpKeyCounter++ ) {
		    
		    if ( ($sampleID) && ($infile) ) {
		
			${$qcDataHashRef}{sample}{$sampleID}{$infile}{$program}{$regExpKey}[$regExpKeyCounter] = ${$qcProgramDataHashRef}{$program}{$regExpKey}[$regExpKeyCounter];

		    }
		    else {

			${$qcDataHashRef}{Program}{$program}{$regExpKey}[$regExpKeyCounter] = ${$qcProgramDataHashRef}{$program}{$regExpKey}[$regExpKeyCounter];			
		    }
		    if ($program eq "SexCheck") {#Check gender for sampleID
			
			my @sexChecks = split(":", @{${$qcProgramDataHashRef}{$program}{$regExpKey}}[$regExpKeyCounter]); #ArrayRef
			## Check that assumed gender is supported by variants on chrX and chrY
			&PlinkGenderCheck({sampleInfoHashRef => $sampleInfoHashRef,
					   qcDataHashRef => $qcDataHashRef,
					   sampleIDRef => \$sexChecks[0],
					   plinkSexCheckGenderRef => \$sexChecks[1],
					  });
		    }
		}
		if (defined(${$qcDataHashRef}{Program}{RelationCheck}{Sample_RelationCheck}) && (defined(${$qcDataHashRef}{Program}{pedigreeCheck}{Sample_order}) ) ) {
		    
		    &RelationCheck({sampleInfoHashRef => $sampleInfoHashRef,
				    qcDataHashRef => $qcDataHashRef,
				    relationshipValuesRef => \@{${$qcDataHashRef}{Program}{RelationCheck}{Sample_RelationCheck}},
				    sampleOrderRef => \@{${$qcDataHashRef}{Program}{pedigreeCheck}{Sample_order}},
				   });

		    delete(${$qcDataHashRef}{Program}{RelationCheck}); #Not of any use anymore
		}
	    }
	}
	else { #Paragraf data i.e. header and subsequent data lines
	    
	    for my $regExpHeaderKey ( keys %{ ${$qcHeaderHashRef}{$program}{$regExpKey} }) { #Find Header info
		
		for my $regExpKeyHeader ( keys %{ ${$regExpHashRef}{$program}{$regExpKey} } ) { #All paragraf keys (header and data line(s))
		    
		    if ($regExpKeyHeader !~/^header|header$/i) { #Detect if the regExp id for headers and not data. 
			
			for (my $qcHeadersCounter=0;$qcHeadersCounter<scalar( @{ ${$qcHeaderHashRef}{$program}{$regExpKey}{$regExpHeaderKey} } );$qcHeadersCounter++) { #For all collected headers
                            
			    if ( ($sampleID) && ($infile)) {
				
				${$qcDataHashRef}{sample}{$sampleID}{$infile}{$program}{$regExpHeaderKey}{$regExpKeyHeader}{ ${$qcHeaderHashRef}{$program}{$regExpKey}{$regExpHeaderKey}[$qcHeadersCounter] } = ${$qcProgramDataHashRef}{$program}{$regExpKey}{$regExpKeyHeader}[$qcHeadersCounter]; #Add to qcData using header element[X] --> data[X] to correctly position elements in qcData hash 
                            } 
                            else {
				
				${$qcDataHashRef}{$program}{$regExpHeaderKey}{$regExpKeyHeader}{ ${$qcHeaderHashRef}{$program}{$regExpKey}{$regExpHeaderKey}[$qcHeadersCounter] } = ${$qcProgramDataHashRef}{$program}{$regExpKey}{$regExpKeyHeader}[$qcHeadersCounter]; #Add to qcData using header element[X] --> data[X] to correctly position elements in qcData hash
				
                            }
                        }    
		    }
		}
	    }
	}
    }
}

sub DefineEvaluateMetric {

##DefineEvaluateMetric

##Function  : Sets programs and program metrics and thresholds to be evaluated
##Returns   : ""
##Arguments : $sampleID
##          : $sampleID => SampleID

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $sampleID;
    
    my $tmpl = { 
	sampleID => { required => 1, defined => 1, strict_type => 1, store => \$sampleID},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    $evaluateMetric{$sampleID}{MosaikAligner}{Total_aligned}{threshold} = 95;
    $evaluateMetric{$sampleID}{MosaikAligner}{Uniquely_aligned_mates}{threshold} = 90;
    $evaluateMetric{$sampleID}{BamStats}{percentage_mapped_reads}{threshold} = 0.95;
    $evaluateMetric{$sampleID}{CalculateHsMetrics}{PCT_TARGET_BASES_10X}{threshold} = 0.95;
    $evaluateMetric{$sampleID}{CollectMultipleMetrics}{PCT_PF_READS_ALIGNED}{threshold} = 0.95;
    $evaluateMetric{$sampleID}{CalculateHsMetrics}{PCT_ADAPTER}{threshold} = 0.0001;

    if ($sampleInfo{sample}{$sampleID}{analysisType} eq "wes") {

	$evaluateMetric{$sampleID}{CalculateHsMetrics}{MEAN_TARGET_COVERAGE}{threshold} = 100;
	$evaluateMetric{$sampleID}{CalculateHsMetrics}{PCT_TARGET_BASES_30X}{threshold} = 0.90;
    }
    else {

	$evaluateMetric{$sampleID}{CalculateHsMetrics}{MEAN_TARGET_COVERAGE}{threshold} = 20;
    }    
}
sub EvaluateQCParameters {

##EvaluateQCParameters

##Function  : Evaluate parameters to detect parameters falling below threshold 
##Returns   : ""
##Arguments : $qcDataHashRef, $evaluateMetricHashRef
##          : $qcDataHashRef         => QCData hash {REF}
##          : $evaluateMetricHashRef => HAsh for metrics to evaluate
    
    my ($argHashRef) = @_;
    
    ## Flatten argument(s)
    my $qcDataHashRef;
    my $evaluateMetricHashRef;
    
    my $tmpl = { 
	qcDataHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$qcDataHashRef},
	evaluateMetricHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$evaluateMetricHashRef},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my $status;

    for my $sampleID ( keys %{${$qcDataHashRef}{sample}} ) { #Can be both sampleID and familyID with current structure
	
	for my $infile ( keys %{${$qcDataHashRef}{sample}{$sampleID}} ) {
	    
	    if ($infile =~/RelationCheck/) { #Special case
		
		if (${$qcDataHashRef}{sample}{$sampleID}{$infile} ne "PASS") {
		    
		    $status = "Status:".$infile.":".${$qcDataHashRef}{sample}{$sampleID}{$infile};
		    push(@{${$qcDataHashRef}{Evaluation}{$infile}}, $status); #Add to QC data at family level
		}
		next;
	    }
	    if ($infile =~/Evaluation/) { #Special case
		
		next;
	    }
	    for my $program ( keys %{${$qcDataHashRef}{sample}{$sampleID}{$infile}} ) {
		
		if (defined(${$evaluateMetricHashRef}{$sampleID}{$program})) { #Program to be evaluated
		    
		    for my $metric ( keys %{${$evaluateMetricHashRef}{$sampleID}{$program}}) { #Metric to be evaluated
			
			if (defined(${$qcDataHashRef}{sample}{$sampleID}{$infile}{$program}{$metric})) {
			    
			    if (${$qcDataHashRef}{sample}{$sampleID}{$infile}{$program}{$metric} < ${$evaluateMetricHashRef}{$sampleID}{$program}{$metric}{threshold}) { #Determine status - if below add to hash. otherwise PASS and do not include
				
				$status = "FAILED:".$sampleID."_".$program."_".$metric.":".${$qcDataHashRef}{sample}{$sampleID}{$infile}{$program}{$metric};
				push(@{${$qcDataHashRef}{Evaluation}{$program}}, $status);
			    }		
			    last;
			}
			else {
			    
			    for my $key ( keys %{${$qcDataHashRef}{sample}{$sampleID}{$infile}{$program}} ) {
				
				if ($key eq "Header") {
				    
				    for my $dataHeader ( keys %{${$qcDataHashRef}{sample}{$sampleID}{$infile}{$program}{$key}} ) {
					
					if (defined(${$qcDataHashRef}{sample}{$sampleID}{$infile}{$program}{$key}{$dataHeader}{$metric})) {
					    
					    if (${$qcDataHashRef}{sample}{$sampleID}{$infile}{$program}{$key}{$dataHeader}{$metric} < ${$evaluateMetricHashRef}{$sampleID}{$program}{$metric}{threshold}) { #Determine status - if below add to hash. otherwise PASS and do not include
						
						$status = "FAILED:".$sampleID."_".$program."_".$metric.":".${$qcDataHashRef}{sample}{$sampleID}{$infile}{$program}{$key}{$dataHeader}{$metric};
						push(@{${$qcDataHashRef}{Evaluation}{$program}}, $status);
					    }
					    next; #Metric go to next section
					}
				    }
				    last; #Metric found no need to continue
				}
			    }
			}
		    }
		}
	    }
	}
    }
}


sub RelationCheck {

##RelationCheck

##Function : Uses the .mibs file produced by PLINK to test if family members are indeed related.
##Returns  : ""
##Arguments: $sampleInfoHashRef, $qcDataHashRef, $relationshipValuesRef, $sampleOrderRef
##         : $sampleInfoHashRef     => Info on samples and family hash {REF}
##         : $qcDataHashRef        => QCData hash {REF}
##         : $sampleOrderRef        => The sample order so that correct estimation can be connected to the correct sampleIDs {REF}
##         : $relationshipValuesRef => All relationship estimations {REF}

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $sampleInfoHashRef;
    my $qcDataHashRef;
    my $relationshipValuesRef;
    my $sampleOrderRef;

    my $tmpl = { 
	sampleInfoHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$sampleInfoHashRef},
	qcDataHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$qcDataHashRef},
	relationshipValuesRef => { required => 1, defined => 1, strict_type => 1, store => \$relationshipValuesRef},
	sampleOrderRef => { required => 1, defined => 1, strict_type => 1, store => \$sampleOrderRef},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my %family; #Stores family relations and pairwise comparisons family{$sampleID}{$sampleID}["column"] -> [pairwise]
    my $sampleIDCounter = 0;
    my $incorrectRelation=0;
    my @pairwiseComparisons;

    ## Splice all relationship extimations from regExp into pairwise comparisons calculated for each sampleID
    for (my $realtionshipCounter=0;$realtionshipCounter<scalar(@{$relationshipValuesRef});$realtionshipCounter++) {
	
	my @pairwiseComparisons = splice(@{$relationshipValuesRef}, 0, scalar(@{$sampleOrderRef})); #Splices array into each sampleIDs line
	
	for (my $column=0;$column<scalar(@{$sampleOrderRef});$column++) { #All columns in .mibs file
	    
	    push ( @{ $family{ @{$sampleOrderRef}[$sampleIDCounter] }{ @{$sampleOrderRef}[$column]} }, $pairwiseComparisons[$column]); #Store sampleID, family membersID (including self) and each pairwise comparison. Uses array for to accomodate sibling info.
	}
	$sampleIDCounter++;
    }
    my $fatherID = "YYY"; #fatherID for the family
    my $motherID = "XXX"; #motherID for the family

    for my $sampleID ( keys %family ) { #For all sampleIDs

	## Currently only 1 father or Mother per pedigree is supported

	if (${$sampleInfoHashRef}{sample}{$sampleID}{father} ne 0) { #Save fatherID if not 0

	    $fatherID = ${$sampleInfoHashRef}{sample}{$sampleID}{father};
	}
	if (${$sampleInfoHashRef}{sample}{$sampleID}{mother} ne 0) { #Save motherID if not 0

	    $motherID = ${$sampleInfoHashRef}{sample}{$sampleID}{mother};
	}
    }
    
    for my $sampleID ( keys %family ) { #For all sampleIDs

	for my $members ( keys %{$family{$sampleID} } ) { #For every relation within family (mother/father/child)
	    
	    for (my $membersCount=0;$membersCount<scalar( @{$family{$sampleID}{$members} } );$membersCount++) { #@ Necessary for siblings
		
		if ($family{$sampleID}{$members}[$membersCount] == 1 ) { #Should only hit self

		    if ( $sampleID eq  $members) {
			#print "Self: ".$sampleID,"\t", $members, "\t", $family{$sampleID}{$members}[$membersCount], "\n";
		    }
		    else {
			$incorrectRelation++;
			${$qcDataHashRef}{sample}{$sampleID}{RelationCheck} = "FAIL: Duplicated sample?;";
			#print  "Incorrect should be self: ".$sampleID,"\t", $members, "\t", $family{$sampleID}{$members}[$membersCount], "\n";
		    }
		}
		elsif ($family{$sampleID}{$members}[$membersCount] >= 0.63 ) { #Should include parent to child and child to siblings unless inbreed parents

		    if ( ( ($sampleID ne $fatherID) && ($sampleID ne $motherID) ) || ( ($members ne $fatherID) && ($members ne $motherID) ) ) { #Correct
			#print "Parent-to-child or child-to-child: ".$sampleID,"\t", $members, "\t", $family{$sampleID}{$members}[$membersCount], "\n";
		    }
		    else {

			$incorrectRelation++;
			${$qcDataHashRef}{sample}{$sampleID}{RelationCheck} = "FAIL: Parents related?;";
			#print "Incorrect: ".$sampleID,"\t", $members, "\t", $family{$sampleID}{$members}[$membersCount], "\n";
		    }
		}
		elsif ($family{$sampleID}{$members}[$membersCount] < 0.63 ) { #Parents unless inbreed

		    if ( ($sampleID eq $fatherID) && ($members eq $motherID) ) {
			#print "Parents: ".$sampleID,"\t", $members, "\t", $family{$sampleID}{$members}[$membersCount], "\n";
		    }
		    elsif ( ($sampleID eq $motherID) && ($members eq $fatherID) ) {
			#print "Parents: ".$sampleID,"\t", $members, "\t", $family{$sampleID}{$members}[$membersCount], "\n";
		    }
		    else {
			$incorrectRelation++;
			${$qcDataHashRef}{sample}{$sampleID}{RelationCheck} = "FAIL:".$sampleID." not related to ".$members.";";
			#print "Incorrect: ".$sampleID,"\t", $members, "\t", $family{$sampleID}{$members}[$membersCount], "\n";
		    }
		}
	    }
	}	
	if ($incorrectRelation == 0) {

	    ${$qcDataHashRef}{sample}{$sampleID}{RelationCheck} = "PASS";  
	}
    }   
    return;
}

sub GenderCheck {

##Function : Checks that the gender predicted by ChanjoSexCheck is confirmed in the pedigee for the sample
##Returns  : ""
##Arguments: $sampleInfoHashRef, $qcDataHashRef, sampleIDRef, $infileRef, $chanjoSexCheckGenderRef
##         : $sampleInfoHashRef       => Info on samples and family hash {REF}
##         : $qcDataHashRef           => QCData hash {REF}
##         : $sampleIDRef             => SampleID {REF}
##         : $infileRef               => Infile {REF}
##         : $chanjoSexCheckGenderRef => Chanjo calculated gender {REF}

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $sampleInfoHashRef;
    my $qcDataHashRef;
    my $sampleIDRef;
    my $infileRef;
    my $chanjoSexCheckGenderRef;

    my $tmpl = { 
	sampleInfoHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$sampleInfoHashRef},
	qcDataHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$qcDataHashRef},
	sampleIDRef => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$sampleIDRef},
	infileRef => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$infileRef},
	chanjoSexCheckGenderRef => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$chanjoSexCheckGenderRef},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!]; 
    
    if ( ($$chanjoSexCheckGenderRef eq "female") && (${$sampleInfoHashRef}{sample}{$$sampleIDRef}{sex} =~/2|female/) ) { #Female
	
	${$qcDataHashRef}{sample}{$$sampleIDRef}{$$infileRef}{GenderCheck} = "PASS";
    }
    elsif ( ($$chanjoSexCheckGenderRef eq "male") && (${$sampleInfoHashRef}{sample}{$$sampleIDRef}{sex} =~/1|male/) ) { #Male
	
	${$qcDataHashRef}{sample}{$$sampleIDRef}{$$infileRef}{GenderCheck} = "PASS";
    }
    else {
	
	${$qcDataHashRef}{sample}{$$sampleIDRef}{$$infileRef}{GenderCheck} = "FAIL";
    }
    return;
}


sub PlinkGenderCheck {

##PlinkGenderCheck

##Function : Checks that the gender predicted by Plink SexCheck is confirmed in the pedigee for the sample
##Returns  : ""
##Arguments: $sampleInfoHashRef, $qcDataHashRef, $sampleIDRef, $plinkSexCheckGenderRef
##         : $sampleInfoHashRef      => Info on samples and family hash {REF}
##         : $qcDataHashRef          => QCData hash {REF}
##         : $sampleIDRef            => SampleID {REF}
##         : $plinkSexCheckGenderRef => Plink calculated gender {REF}

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $sampleInfoHashRef;
    my $qcDataHashRef;
    my $sampleIDRef;
    my $plinkSexCheckGenderRef;

    my $tmpl = {
	sampleInfoHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$sampleInfoHashRef},
	qcDataHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$qcDataHashRef},
	sampleIDRef => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$sampleIDRef},
	plinkSexCheckGenderRef => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$plinkSexCheckGenderRef},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];
    
    if ( ($$plinkSexCheckGenderRef eq "2") && (${$sampleInfoHashRef}{sample}{$$sampleIDRef}{sex} =~/2|female/) ) { #Female
	
	push(@{${$qcDataHashRef}{Program}{PlinkGenderCheck}}, $$sampleIDRef.":PASS");
    }
    elsif ( ($$plinkSexCheckGenderRef eq "1") && (${$sampleInfoHashRef}{sample}{$$sampleIDRef}{sex} =~/1|male/) ) { #Male
	
	push(@{${$qcDataHashRef}{Program}{PlinkGenderCheck}}, $$sampleIDRef.":PASS");
    }
    else {
	
	push(@{${$qcDataHashRef}{Program}{PlinkGenderCheck}}, $$sampleIDRef.":FAIL");
    }
    return;
}

sub WriteYAML {

##WriteYAML
    
##Function : Writes a YAML hash to file
##Returns  : ""
##Arguments: $yamlHashRef, $yamlFilePathRef
##         : $yamlHashRef     => The hash to dump {REF}
##         : $yamlFilePathRef => The yaml file to write to {REF}

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $yamlHashRef;
    my $yamlFilePathRef;

    my $tmpl = { 
	yamlHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$yamlHashRef},
	yamlFilePathRef => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$yamlFilePathRef},
    };

    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    open (my $YAML, ">", $$yamlFilePathRef) or die "can't open ".$$yamlFilePathRef.": $!\n";
    say $YAML Dump( $yamlHashRef );
    close($YAML);
}

sub LoadYAML {

##LoadYAML
    
##Function : Loads a YAML file into an arbitrary hash and returns it. Note: Currently only supports hashreferences and hashes and no mixed entries.
##Returns  : %yamlHash
##Arguments: $yamlFile
##         : $yamlFile => The yaml file to load

    my ($argHashRef) = @_;

    ##Flatten argument(s)
    my $yamlFile;

    my $tmpl = { 
	yamlFile => { required => 1, defined => 1, strict_type => 1, store => \$yamlFile},
    };

    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my %yamlHash;

    open (my $YAML, "<", $yamlFile) or die "can't open ".$yamlFile.":".$!, "\n"; 

    %yamlHash = %{ YAML::LoadFile($yamlFile) };  #Load hashreference as hash

    close($YAML);

    return %yamlHash;
}


sub RegExpToYAML {

##RegExpToYAML

##Function : Write default regExp to YAML
##Returns  : ""
##Arguments: $printRegExpOutFile
##         : $printRegExpOutFile => File to print regExp to

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $printRegExpOutFile;

    my $tmpl = { 
	printRegExpOutFile => { required => 1, defined => 1, strict_type => 1, store => \$printRegExpOutFile},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my %regExp;

    #Add to %regExp to enable print in YAML
    $regExp{FastQC}{Version} = q?perl -nae' if ($_=~/##FastQC\\s+(\\S+)/) {print $1;last;}' ?; #Collect FastQC version

    $regExp{FastQC}{Encoding} = q?perl -nae' if ($_=~/Encoding\s+(\S+\s\S+\s\S+\s\S+|\S+\s\S+)/) { my $encoding = $1;$encoding=~s/\s/\_/g; print $encoding;last;}' ?; #Collect Encoding
    
    $regExp{FastQC}{Sequence_length} = q?perl -nae' if ($_=~/Sequence length\s(\d+)/) {print $1;last;}' ?; #Collect Sequence length
    
    $regExp{FastQC}{Total_number_of_reads} = q?perl -nae' if ($_=~/Total Sequences\s(\d+)/) {print $1;last;}' ?; #Collect Total sequences 
    
    $regExp{FastQC}{GC} = q?perl -nae' if ($_=~/%GC\s(\d+)/) {print $1;last;}' ?; #Collect GC content 
    
    $regExp{FastQC}{Sequence_duplication} = q?perl -nae' if ($_=~/#Total Duplicate Percentage\s+(\d+.\d)/) {print $1;last;}' ?; #Collect Sequence duplication level
    
    $regExp{FastQC}{Basic_statistics} = q?perl -nae' if ($_=~/>>Basic Statistics\s+(\S+)/) {print $1;last;}' ?; #Collect Basic Statistics
    
    $regExp{FastQC}{Per_base_sequence_quality} = q?perl -nae' if ($_=~/>>Per base sequence quality\s+(\S+)/) {print $1;last;}' ?; #Collect Per base sequence quality
    
    $regExp{FastQC}{Per_sequence_quality_scores} = q?perl -nae' if ($_=~/>>Per sequence quality scores\s+(\S+)/) {print $1;last;}' ?; #Collect Per sequence quality scores
    
    $regExp{FastQC}{Per_base_sequence_content} = q?perl -nae' if ($_=~/>>Per base sequence content\s+(\S+)/) {print $1;last;}' ?; #Collect Per base sequence content
    
    $regExp{FastQC}{Per_base_GC_content} = q?perl -nae' if ($_=~/>>Per base GC content\s+(\S+)/) {print $1;last;}' ?; #Collect Per base GC content
    
    $regExp{FastQC}{Per_sequence_GC_content} = q?perl -nae' if ($_=~/>>Per sequence GC content\s+(\S+)/) {print $1;last;}' ?; #Collect Per sequence GC content
    
    $regExp{FastQC}{Per_base_N_content} = q?perl -nae' if ($_=~/>>Per base N content\s+(\S+)/) {print $1;last;}' ?; #Collect Per base N content
    
    $regExp{FastQC}{Sequence_duplication_levels} = q?perl -nae' if ($_=~/>>Sequence Duplication Levels\s+(\S+)/) {print $1;last;}' ?; #Collect Sequence Duplication Levels
    
    $regExp{FastQC}{Overrepresented_sequences} = q?perl -nae' if ($_=~/>>Overrepresented sequences\s+(\S+)/) {print $1;last;}' ?; #Collect Overrepresented sequences
    
    $regExp{FastQC}{Kmer_content} = q?perl -nae' if ($_=~/>>Kmer Content\s+(\S+)/) {print $1;last;}' ?; #Collect Kmer Content
    
    $regExp{MosaikAligner}{Version} = q?perl -nae' if ($_=~/(\d+\.\d+\.\d+)\s/) {print $1;last;}' ?; #Collect Mosaik Version 
    
    $regExp{MosaikAligner}{Unaligned_mates} = q?perl -nae' if ($_=~/# unaligned mates\S+\s+(\d+)\s\(\s+(\d+\.\d+)/) {print $2;last;}' ?; #Collect Nr of unaligned mates
    
    $regExp{MosaikAligner}{Filtered_out} = q?perl -nae' if ($_=~/# filtered out\S+\s+(\d+)\s\(\s+(\d+\.\d+)/) {print $2;last;}' ?; #Collect Nr of filtered out reads 
    
    $regExp{MosaikAligner}{Uniquely_aligned_mates} = q?perl -nae' if ($_=~/# uniquely aligned mates\S+\s+(\d+)\s\(\s+(\d+\.\d+)/) {print $2;last;}' ?; #Collect Uniquely aligned mates
    
    $regExp{MosaikAligner}{Multiply_aligned_mates} = q?perl -nae' if ($_=~/# multiply aligned mates\S+\s+(\d+)\s\(\s+(\d+\.\d+)/) {print $2;last;}' ?; #Collect Multiply aligned mates
    
    $regExp{MosaikAligner}{Total_aligned} = q?perl -nae' if ($_=~/total aligned:\s+\S+\s+(\S+)\s\(\S+\s(\d+.\d+)/ ) {print $2;last;} elsif ($_=~/total aligned:\s+(\S+)\s\(\S+\s(\d+.\d+)/ ) { print $2;last;}' ?; #Collect total aligned sequences

    $regExp{BamStats}{percentage_mapped_reads} = q?perl -nae 'if($_=~/percentage mapped reads:\s+(\S+)/) {print $1;last}' ?; #Collect % mapped reads from BAM alignment

    $regExp{BamStats}{raw_total_sequences} = q?perl -nae 'if($_=~/raw total sequences:\s+(\S+)/) {print $1;last}' ?; #Collect raw total sequences from BAM alignment

    $regExp{BamStats}{reads_mapped} = q?perl -nae 'if($_=~/reads mapped:\s+(\S+)/) {print $1;last}' ?; #Collect reads mapped from BAM alignment

    $regExp{ChanjoSexCheck}{gender} = q?perl -nae 'if( ($F[0]!~/^#/) && ($F[2] =~/\S+/) ) {print $F[2];}' ?;  #Collect gender from ChanjoSexCheck

    $regExp{pedigreeCheck}{Sample_order} = q?perl -nae 'if ($_=~/^#CHROM/) {chomp $_; my @line = split(/\t/,$_); for (my $sample=9;$sample<scalar(@line);$sample++) { print $line[$sample], "\t";}last;}' ?; #Collect sample order from vcf file used to create ".ped", ".map" and hence ".mibs".
    
    $regExp{InbreedingFactor}{Sample_InbreedingFactor}  = q?perl -nae 'my @inbreedingFactor; if ($. > 1) {my @temp = split(/\s/,$_);push(@inbreedingFactor, $F[0].":".$F[5]); print $inbreedingFactor[0], "\t"; }' ?;

    $regExp{SexCheck}{Sample_SexCheck}  = q?perl -nae 'my @sexCheckFactor; if ($. > 1) {my @temp = split(/\s+/,$_);push(@sexCheckFactor,$temp[2].":".$temp[4]); print $sexCheckFactor[0], "\t"; }' ?;

    $regExp{RelationCheck}{Sample_RelationCheck}  = q?perl -nae 'print $_;' ?; #Note will return whole file

    $regExp{MarkDuplicates}{Fraction_duplicates} = q?perl -nae 'if($_=~/Fraction Duplicates\: (\S+)/) {print $1;}' ?; #Collect fraction duplicates
    
    $regExp{CalculateHsMetrics}{Header_info}{Header} = q?perl -nae' if ($_ =~/^BAIT_SET/ ) {print $_;last;}' ?; #Note return whole line (Header) 
    
    $regExp{CalculateHsMetrics}{Header_info}{Data} = q?perl -nae' if ( ($. ==8) && ($_ =~/(\S+)/) ) {print $_;last;}' ?; #Note return whole line and only look at line 8, where the data action is
    
    $regExp{CollectMultipleMetrics}{Header_info}{Header} = q?perl -nae' if ($_ =~/^CATEGORY/ ) {print $_;last;}' ?; #Note return whole line (Header) 
    
    $regExp{CollectMultipleMetrics}{Header_info}{First_of_pair} = q?perl -nae' if ($_ =~/^FIRST_OF_PAIR/ ) {print $_;last;}' ?; #Note return whole line (FIRST_OF_PAIR)
    
    $regExp{CollectMultipleMetrics}{Header_info}{Second_of_pair} = q?perl -nae' if ($_ =~/^SECOND_OF_PAIR/ ) {print $_;last;}' ?; #Note return whole line (SECOND_OF_PAIR)  
    
    $regExp{CollectMultipleMetrics}{Header_info}{Pair} = q?perl -nae' if ($_ =~/^PAIR/ ) {print $_;last;}'  ?; #Note return whole line (PAIR)
    
    $regExp{CollectMultipleMetricsInsertSize}{Header_info}{Header} = q?perl -nae' if ($_ =~/^MEDIAN_INSERT_SIZE/ ) {print $_;last;}' ?; #Note return whole line (Header) 
    
    $regExp{CollectMultipleMetricsInsertSize}{Header_info}{Data} = q?perl -nae' if ( ($. ==8) && ($_ =~/(\S+)/) ) {print $_;last;}' ?; #Note return whole line and only look at line 8, where the data action is

    $regExp{VariantEval_All}{CompOverlap_header}{CompOverlap_header} = q?perl -nae' if ($_ =~/^CompOverlap\s+\CompRod/ ) {print $_;last;}' ?; #Note return whole line (Header)
    
    $regExp{VariantEval_All}{CompOverlap_header}{CompOverlap_data_all} = q?perl -nae' if ( ($_ =~/^CompOverlap/) && ($_ =~/all/) && ($_ =~/none/)) {print $_;last;}' ?; #Note return whole line
    
    $regExp{VariantEval_All}{CompOverlap_header}{CompOverlap_data_known} = q?perl -nae' if ( ($_ =~/^CompOverlap/) && ($_ =~/known\s/) ) {print $_;last;}' ?; #Note return whole line

    $regExp{VariantEval_All}{CompOverlap_header}{CompOverlap_data_novel} = q?perl -nae' if ( ($_ =~/^CompOverlap/) && ($_ =~/novel\s/) ) {print $_;last;}' ?; #Note return whole line

    $regExp{VariantEval_All}{CountVariants_header}{CountVariants_header} = q?perl -nae' if ($_ =~/^CountVariants\s+\CompRod/ ) {print $_;last;}' ?; #Note return whole line (Header)                                               
    $regExp{VariantEval_All}{CountVariants_header}{CountVariants_data_all} = q?perl -nae' if ( ($_ =~/^CountVariants/) && ($_ =~/all\s/) ) {print $_;last;}' ?; #Note return whole line                                              
    $regExp{VariantEval_All}{CountVariants_header}{CountVariants_data_known} = q?perl -nae' if ( ($_ =~/^CountVariants/) && ($_ =~/known\s/) ) {print $_;last;}' ?; #Note return whole line                                          
    $regExp{VariantEval_All}{CountVariants_header}{CountVariants_data_novel} = q?perl -nae' if ( ($_ =~/^CountVariants/) && ($_ =~/novel\s/) ) {print $_;last;}' ?; #Note return whole line

    $regExp{VariantEval_All}{IndelSummary_header}{IndelSummary_header} = q?perl -nae' if ($_ =~/^IndelSummary\s+\CompRod/ ) {print $_;last;}' ?; #Note return whole line (Header)                                                                                                                                                                                                                         
    $regExp{VariantEval_All}{IndelSummary_header}{IndelSummary_data_all} = q?perl -nae' if ( ($_ =~/^IndelSummary/) && ($_ =~/all\s/) ) {print $_;last;}' ?; #Note return whole line                                                                                                                                                                                                        
    $regExp{VariantEval_All}{IndelSummary_header}{IndelSummary_data_known} = q?perl -nae' if ( ($_ =~/^IndelSummary/) && ($_ =~/known\s/) ) {print $_;last;}' ?; #Note return whole line                                     
    $regExp{VariantEval_All}{IndelSummary_header}{IndelSummary_data_novel} = q?perl -nae' if ( ($_ =~/^IndelSummary/) && ($_ =~/novel\s/) ) {print $_;last;}' ?; #Note return whole line 

    $regExp{VariantEval_All}{MultiallelicSummary_header}{MultiallelicSummary_header} = q?perl -nae' if ($_ =~/^MultiallelicSummary\s+\CompRod/ ) {print $_;last;}' ?; #Note return whole line (Header)                         
    $regExp{VariantEval_All}{MultiallelicSummary_header}{MultiallelicSummary_data_all} = q?perl -nae' if ( ($_ =~/^MultiallelicSummary/) && ($_ =~/all\s/) ) {print $_;last;}' ?; #Note return whole line                                          
    $regExp{VariantEval_All}{MultiallelicSummary_header}{MultiallelicSummary_data_known} = q?perl -nae' if ( ($_ =~/^MultiallelicSummary/) && ($_ =~/known\s/) ) {print $_;last;}' ?; #Note return whole line                                        
    $regExp{VariantEval_All}{MultiallelicSummary_header}{MultiallelicSummary_data_novel} = q?perl -nae' if ( ($_ =~/^MultiallelicSummary/) && ($_ =~/novel\s/) ) {print $_;last;}' ?; #Note return whole line

    $regExp{VariantEval_All}{TiTvVariantEvaluator_header}{TiTvVariantEvaluator_header} = q?perl -nae' if ($_ =~/^TiTvVariantEvaluator\s+\CompRod/ ) {print $_;last;}' ?; #Note return whole line (Header)                                                                                                                                                                                        
    $regExp{VariantEval_All}{TiTvVariantEvaluator_header}{TiTvVariantEvaluator_data_all} = q?perl -nae' if ( ($_ =~/^TiTvVariantEvaluator/) && ($_ =~/all\s/) ) {print $_;last;}' ?; #Note return whole line                                                                                                                                                                                                       
    $regExp{VariantEval_All}{TiTvVariantEvaluator_header}{TiTvVariantEvaluator_data_known} = q?perl -nae' if ( ($_ =~/^TiTvVariantEvaluator/) && ($_ =~/known\s/) ) {print $_;last;}' ?; #Note return whole line                                                                                                                                                                                               
    $regExp{VariantEval_All}{TiTvVariantEvaluator_header}{TiTvVariantEvaluator_data_novel} = q?perl -nae' if ( ($_ =~/^TiTvVariantEvaluator/) && ($_ =~/novel\s/) ) {print $_;last;}' ?; #Note return whole line 

    $regExp{VariantEval_All}{ValidationReport_header}{ValidationReport_header} = q?perl -nae' if ($_ =~/^ValidationReport\s+\CompRod/ ) {print $_;last;}' ?; #Note return whole line (Header)                                                                                                                                                                                                                          
    $regExp{VariantEval_All}{ValidationReport_header}{ValidationReport_data_all} = q?perl -nae' if ( ($_ =~/^ValidationReport/) && ($_ =~/all\s/) && ($_ =~/none\s/)) {print $_;last;}' ?; #Note return whole line                                                                                                                                                                                                                         
    $regExp{VariantEval_All}{ValidationReport_header}{ValidationReport_data_known} = q?perl -nae' if ( ($_ =~/^ValidationReport/) && ($_ =~/known\s/) ) {print $_;last;}' ?; #Note return whole line                                                                                                                                                                                                        
    $regExp{VariantEval_All}{ValidationReport_header}{ValidationReport_data_novel} = q?perl -nae' if ( ($_ =~/^ValidationReport/) && ($_ =~/novel\s/) ) {print $_;last;}' ?; #Note return whole line

    $regExp{VariantEval_All}{VariantSummary_header}{VariantSummary_header} = q?perl -nae' if ($_ =~/^VariantSummary\s+\CompRod/ ) {print $_;last;}' ?; #Note return whole line (Header)                                                                                                                                                                                                                     
    $regExp{VariantEval_All}{VariantSummary_header}{VariantSummary_data_all} = q?perl -nae' if ( ($_ =~/^VariantSummary/) && ($_ =~/all\s/) ) {print $_;last;}' ?; #Note return whole line                                                                                                                                                                                                                   
    $regExp{VariantEval_All}{VariantSummary_header}{VariantSummary_data_known} = q?perl -nae' if ( ($_ =~/^VariantSummary/) && ($_ =~/known\s/) ) {print $_;last;}' ?; #Note return whole line                                                                                                                                                                                                                
    $regExp{VariantEval_All}{VariantSummary_header}{VariantSummary_data_novel} = q?perl -nae' if ( ($_ =~/^VariantSummary/) && ($_ =~/novel\s/) ) {print $_;last;}' ?; #Note return whole line

    $regExp{VariantEval_Exome} = $regExp{VariantEval_All};

    $regExp{Genmod}{Version} = q?perl -nae 'if($_=~/##Software=<ID=genmod\w+,Version=(\d+.\d+.\d+)/) {print $1;last;}' ?; #Collect Genmod version

    $regExp{SnpEff}{Version} = q?perl -nae 'if($_=~/##SnpSiftVersion=\"(.+),/) {my $ret=$1; $ret=~s/\s/_/g;print $ret;last;}' ?; #Collect SnpEff version

    $regExp{VariantEffectPredictor}{Version} = q?perl -nae 'if($_=~/##VEP=(\w+)/) {print $1;last;}' ?; #Collect VariantEffectPredictor version

    $regExp{VariantEffectPredictor}{Cache} = q?perl -nae 'if($_=~/##VEP=\w+\s+cache=(\S+)/) {print $1;last;}' ?; #Collect VariantEffectPredictor cache directory

    $regExp{VariantEffectPredictor}{PolyPhen} = q?perl -nae 'if($_=~/##VEP=/ && $_=~/polyphen=(\S+)/) {print $1;last;}' ?; #Collect VariantEffectPredictor polyPhen version

    $regExp{VariantEffectPredictor}{Sift} = q?perl -nae 'if($_=~/##VEP=/ && $_=~/sift=sift(\S+)/) {print $1;last;}' ?; #Collect VariantEffectPredictor sift version

    $regExp{VariantEffectPredictor}{GeneBuild} = q?perl -nae 'if($_=~/##VEP=/ && $_=~/genebuild=(\S+)/) {print $1;last;}' ?; #Collect VariantEffectPredictor geneBuild

    $regExp{VariantEffectPredictor}{Assembly} = q?perl -nae 'if($_=~/##VEP=/ && $_=~/assembly=(\S+)/) {print $1;last;}' ?; #Collect VariantEffectPredictor assembly

    $regExp{VariantEffectPredictor}{HGMDPublic} = q?perl -nae 'if($_=~/##VEP=/ && $_=~/HGMD-PUBLIC=(\S+)/) {print $1;last;}' ?; #Collect VariantEffectPredictor HGMD-PUBLIC version

    $regExp{VariantEffectPredictor}{RegBuild} = q?perl -nae 'if($_=~/##VEP=/ && $_=~/regbuild=(\S+)/) {print $1;last;}' ?; #Collect VariantEffectPredictor regbuild version

    $regExp{VariantEffectPredictor}{Gencode} = q?perl -nae 'if($_=~/##VEP=/ && $_=~/gencode=\S+\s+(\d+)/) {print $1;last;}' ?; #Collect VariantEffectPredictor gencode version

    $regExp{VCFParser}{Version} = q?perl -nae 'if($_=~/##Software=<ID=vcfParser.pl,Version=(\d+.\d+.\d+)/) {print $1;last;}' ?; #Collect VCFParser version

    $regExp{Bwa}{Version} = q?perl -nae 'if($_=~/\[main\]\sVersion:\s(\S+)/) {print $1;last;}' ?; #Collect Bwa version

    $regExp{Chanjo}{Version} = q?perl -nae 'if($_=~/version\s(\d+.\d+.\d+)/) {print $1;last;}' ?; #Collect Chanjo version

    $regExp{vt}{Version} = q?perl -nae 'if($_=~/decompose\sv(\S+)/) {print $1;last;}' ?; #Collect vt version

    $regExp{Samtools}{Version} = q?perl -nae 'if($_=~/samtoolsVersion=(\S+)/) {print $1;last;}' ?; #Collect Samtools version

    $regExp{Bcftools}{Version} = q?perl -nae 'if($_=~/bcftools_\w+Version=(\S+)/) {print $1;last;}' ?; #Collect Bcftools version

    $regExp{Freebayes}{Version} = q?perl -nae 'if($_=~/source=freeBayes\s(\S+)/) {print $1;last;}' ?; #Collect Freebayes version

    $regExp{Delly}{Version} = q?perl -nae 'if($_=~/SVMETHOD=EMBL\.DELLY(v\d+\.\d+\.\d+)/) {print $1;last }' ?; #Collect Delly version

    $regExp{Manta}{Version} = q?perl -nae 'if($_=~/GenerateSVCandidates\s+(\S+)/) {print $1;last}' ?; #Collect Manta version

    $regExp{SVCombineVariantCallSets}{VcfAnno} = q?perl -nae 'if($_=~/vcfanno\sversion\s(\S+)/) {print $1;last;}' ?; #Collect SVVCFAnno version

    $regExp{SVVariantEffectPredictor}{Version} = q?perl -nae 'if($_=~/##VEP=(\w+)/) {print $1;last;}' ?; #Collect SVVariantEffectPredictor version

    $regExp{SVVariantEffectPredictor}{Cache} = q?perl -nae 'if($_=~/##VEP=\w+\s+cache=(\S+)/) {print $1;last;}' ?; #Collect SVVariantEffectPredictor cache directory

    $regExp{SVVariantEffectPredictor}{PolyPhen} = q?perl -nae 'if($_=~/##VEP=/ && $_=~/polyphen=(\S+)/) {print $1;last;}' ?; #Collect SVVariantEffectPredictor polyPhen version

    $regExp{SVVariantEffectPredictor}{Sift} = q?perl -nae 'if($_=~/##VEP=/ && $_=~/sift=sift(\S+)/) {print $1;last;}' ?; #Collect SVVariantEffectPredictor sift version

    $regExp{SVVariantEffectPredictor}{GeneBuild} = q?perl -nae 'if($_=~/##VEP=/ && $_=~/genebuild=(\S+)/) {print $1;last;}' ?; #Collect SVVariantEffectPredictor geneBuild

    $regExp{SVVariantEffectPredictor}{Assembly} = q?perl -nae 'if($_=~/##VEP=/ && $_=~/assembly=(\S+)/) {print $1;last;}' ?; #Collect SVVariantEffectPredictor assembly

    $regExp{SVVariantEffectPredictor}{HGMDPublic} = q?perl -nae 'if($_=~/##VEP=/ && $_=~/HGMD-PUBLIC=(\S+)/) {print $1;last;}' ?; #Collect SVVariantEffectPredictor HGMD-PUBLIC version

    $regExp{SVVariantEffectPredictor}{RegBuild} = q?perl -nae 'if($_=~/##VEP=/ && $_=~/regbuild=(\S+)/) {print $1;last;}' ?; #Collect SVVariantEffectPredictor regbuild version

    $regExp{SVVariantEffectPredictor}{Gencode} = q?perl -nae 'if($_=~/##VEP=/ && $_=~/gencode=\S+\s+(\d+)/) {print $1;last;}' ?; #Collect SVVariantEffectPredictor gencode version

    $regExp{SVVCFParser}{Version} = q?perl -nae 'if($_=~/##Software=<ID=vcfParser.pl,Version=(\d+.\d+.\d+)/) {print $1;last;}' ?; #Collect SVVCFParser version

    $regExp{SVGenmod}{Version} = q?perl -nae 'if($_=~/##Software=<ID=genmod\w+,Version=(\d+.\d+.\d+)/) {print $1;last;}' ?; #Collect SVGenmod version

    $regExp{Vcftools}{Version} = q?perl -nae 'if($_=~/VCFtools\s-\s(\d+.\d+.\d+)/) {print $1;last;}' ?; #Collect VCFTools version

    $regExp{Plink2}{Version} = q?perl -nae 'if($_=~/PLINK\s(\S+\s\S+\s\S+\s\S+\s\S+)/) {my $ret = $1;$ret =~s/\s/_/g;print $ret;last;}' ?; #Collect Plink2 version

    $regExp{VariantIntegrityMendel}{FractionofErrors}  = q?perl -nae 'unless ($_=~/^#/) {print $F[1];last;}' ?;

    $regExp{VariantIntegrityMendel}{MendelianErrors}  = q?perl -nae 'unless ($_=~/^#/) {print $F[2];last;}' ?;

    $regExp{VariantIntegrityFather}{FractionofCommonVariants}  = q?perl -nae 'unless ($_=~/^#/) {print $F[1];last;}' ?;

    $regExp{VariantIntegrityFather}{CommonVariants}  = q?perl -nae 'unless ($_=~/^#/) {print $F[2];last;}' ?;

#$regExp{}{} = ;

    ## Writes a YAML hash to file
    &WriteYAML({yamlHashRef => \%regExp,
		yamlFilePathRef => \$printRegExpOutFile,
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

    ## Flatten argument(s)
    my $USAGE;
    my $exitCode;

    my $tmpl = { 
	USAGE => {required => 1, defined => 1, strict_type => 1, store => \$USAGE},
	exitCode => { default => 0, strict_type => 1, store => \$exitCode},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];    
    
    say STDOUT $USAGE;
    exit $exitCode;
}
