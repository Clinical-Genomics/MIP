#!/usr/bin/env perl

use strict;
use warnings;

##Collects MPS QC from MIP. Loads information on files to examine and values to extract from in YAML format and outputs exracted metrics in YAML format. 
#Copyright 2013 Henrik Stranneheim 

use Pod::Usage;
use Pod::Text;
use Getopt::Long;
use POSIX;
use YAML;

use vars qw($USAGE);

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

my ($sampleInfoFile, $regExpFile, $outfile, $printRegExp, $printRegExpOutFile, $version, $help) = (0,0,"qcmetrics.yaml", 0, "qc_regExp.yaml", 0, 0);
my (%sampleInfo, %regExp, %qcData, %evaluateMetric);
my %qcHeader; #Save header(s) in each outfile
my %qcProgramData; #Save data in each outFile

GetOptions('si|sampleInfoFile:s' => \$sampleInfoFile,
	   'r|regExpFile:s' => \$regExpFile,
	   'o|outfile:s'  => \$outfile, 
	   'preg|printRegExp:n' => \$printRegExp,
	   'prego|printRegExpOutFile:s' => \$printRegExpOutFile,
	   'h|help' => \$help,
	   'v|version' => \$version,
    );

if($help) {

    print STDOUT "\n".$USAGE, "\n";
    exit;
}

my $qcCollectVersion = "1.0.2";
if($version) {

    print STDOUT "\nqcCollect.pl v".$qcCollectVersion,"\n\n";
    exit;
}

if ($printRegExp) {

    RegExpToYAML($printRegExpOutFile); #Write regExp used @ CMMS to YAML
    print STDOUT "Wrote RegExp YAML file to: ".$printRegExpOutFile, "\n";
    exit;
}

if ($sampleInfoFile eq 0) {

    print STDERR "\n";
    print STDOUT $USAGE, "\n";
    print STDERR "Must supply a '-sampleInfoFile' (supply whole path)", "\n\n";
    exit;
}
if ($regExpFile eq 0) {

    print STDERR "\n";
    print STDOUT $USAGE, "\n";
    print STDERR "Must supply a '-regExpFile' (supply whole path)", "\n\n";
    exit;
}

####MAIN

my %sampleInfoFile = &LoadYAML($sampleInfoFile); #Load sampleInfoFile (YAML) and transfer to sampleInfoFile (Hash)

my %regExpFile = &LoadYAML($regExpFile); #Load regExpFile (YAML) and transfer to regExpFile (Hash)

&SampleQC(); #Extracts all qcdata on sampleID level using information in %sampleInfoFile and %regExpFile

&FamilyQC(); #Extracts all qcdata on family level using information in %sampleInfoFile and %regExpFile

##Add qcCollect version to yaml file
for my $familyID ( keys %sampleInfoFile ) { #For every family id

    $qcData{$familyID}{$familyID}{'Program'}{'QCCollect'}{'Version'} = $qcCollectVersion;
    $qcData{$familyID}{$familyID}{'Program'}{'QCCollect'}{'RegExpFile'} = $regExpFile;

    &DefineEvaluateMetric($familyID); #Defines programs, etrics and thresholds to evaluate
}

&EvaluateQCParameters(); #Evaluate the metrics

&WriteYAML($outfile, \%qcData ); #Writes to YAML file

####SubRoutines

sub FamilyQC {

##FamilyQC

##Function  : Extracts all qcdata on family level using information in %sampleInfoFile and %regExpFile
##Returns   : ""
##Arguments : 

    for my $familyID ( keys %sampleInfoFile ) { #For every family id 

        if ($sampleInfoFile{$familyID}{$familyID}{'Program'}) { #Only examine programs     

            for my $program ( keys %{ $sampleInfoFile{$familyID}{$familyID}{'Program'} } ) { #For every programs           

		my $outDirectory;
		my $outFile;	

                if ($sampleInfoFile{$familyID}{$familyID}{'Program'}{$program}{'Version'} ) {
                    
                    $qcData{$familyID}{$familyID}{'Program'}{$program}{'Version'} = $sampleInfoFile{$familyID}{$familyID}{'Program'}{$program}{'Version'}; #Add version to qcData
                }
                if ($sampleInfoFile{$familyID}{$familyID}{'Program'}{$program}{'OutDirectory'} ) {

                    $outDirectory = $sampleInfoFile{$familyID}{$familyID}{'Program'}{$program}{'OutDirectory'}; #Extract OutDirectory
                }
                if ($sampleInfoFile{$familyID}{$familyID}{'Program'}{$program}{'OutFile'} ) {

                    $outFile = $sampleInfoFile{$familyID}{$familyID}{'Program'}{$program}{'OutFile'}; #Extract OutFile
                }  
		
                &ParseRegExpHashAndCollect($program, $outDirectory, $outFile); #Loads qcHeader and qcProgramData

                &AddToqcData($familyID, "", $program, ""); #Add extracted information to qcData
            }
        }
    }
}

sub SampleQC {

##SampleQC

##Function  : Collects all sample qc in files defined by sampleInfoFile and regular expressions defined by regExpFile.
##Returns   : ""
##Arguments : 
    
    for my $familyID ( keys %sampleInfoFile ) { #For every family id
	
	for my $sampleID ( keys %{ $sampleInfoFile{$familyID} } ) { #For every sample id

	    unless ($sampleID eq $familyID) { #Family data is on the same level sa sampleIDs

		if ($sampleInfoFile{$familyID}{$sampleID}{'Program'}) { #Only examine programs
		    
		    for my $program ( keys %{ $sampleInfoFile{$familyID}{$sampleID}{'Program'} } ) { #For every program  
			
			for my $infile ( keys %{ $sampleInfoFile{$familyID}{$sampleID}{'Program'}{$program} } ) { #For every infile
			    
			    my $outDirectory = $sampleInfoFile{$familyID}{$sampleID}{'Program'}{$program}{$infile}{'OutDirectory'};
			    my $outFile = $sampleInfoFile{$familyID}{$sampleID}{'Program'}{$program}{$infile}{'OutFile'};
			    
			    &ParseRegExpHashAndCollect($program, $outDirectory, $outFile); #Loads qcHeader and qcProgramData
			    
			    &AddToqcData($familyID, $sampleID, $program, $infile); #Add extracted information to qcData
			    
			}   
		    }
		}
	    }
	}
    }   
}

sub ParseRegExpHashAndCollect {

##ParseRegExpHashAndCollect

##Function  : Parses the RegExpHash structure to identify if the info is 1) Paragraf section(s) (both header and data line(s)); 2) Seperate data line.
##Returns   : ""
##Arguments : $program, $outDirectory, $outFile
##          : $program      => The program to examine
##          : $outDirectory => Programs out directory
##          : $outFile      => Programs out file containing parameter to evaluate

    my $program = $_[0];
    my $outDirectory = $_[1];
    my $outFile = $_[2]; 
    
    my $regExp; #Holds the current regExp
    my @separators = ('\s+','!',','); #Covers both whitespace and tab. Add other separators if required
    
    for my $regExpKey ( keys %{ $regExpFile{$program} } ) { #Find the actual regular expression(s) for each program that is used
	
	if ($regExpKey =~/^header|header$/i) { #Detect if the outfile contains paragrafs/header info in the outfile i.e. data is formated as a paragraf with header(s) and line(s). "regExpKey" should either start with or end with "header". This section extracts the header line(s) for the entire outdata file. Necessary to assign correct data entry to header entry later (headers and data are saved in seperate hashes). 
	    
	    ##Format outfile: Paragraf section
	    for my $regExpHeaderKey ( keys %{ $regExpFile{$program}{$regExpKey} } ) { #Paragraf Header 
		
		$regExp = $regExpFile{$program}{$regExpKey}{$regExpHeaderKey}; #The regular expression used to collect paragraf header info
		
		for (my $separatorElement=0;$separatorElement<scalar(@separators);$separatorElement++) { #Loop through possible separators to seperate any eventual header elements
                    
		    if ($regExpHeaderKey =~/^header|header$/i) { #Detect if the regExp key is a paragraf header and not paragraf data (Header line and data line(s))
                        
			@{ $qcHeader{$program}{$regExpKey}{$regExpHeaderKey} } = split(/$separators[$separatorElement]/, `$regExp $outDirectory/$outFile`); #Collect paragraf header                           
			if ( defined($qcHeader{$program}{$regExpKey}{$regExpHeaderKey})) { #Then split should have been successful                                                                                          
			    last; #Found correct separator - do not continue                                                                                         
			}   
		    }
		    else { #For paragraf data line(s)                                                                                                                        
			@{ $qcProgramData{$program}{$regExpKey}{$regExpHeaderKey} } = split(/$separators[$separatorElement]/, `$regExp $outDirectory/$outFile`); #Collect paragraf data
			if ( defined($qcProgramData{$program}{$regExpKey}{$regExpHeaderKey}[1])) { #Then split should have been successful                                                                                                                           
			    last; #Found correct separator - do not continue                                                                                                                                      
			}
		    }
		}
	    }
	}
	else { #For info contained in Entry --> Value i.e. same line.

	    $regExp = $regExpFile{$program}{$regExpKey}; #The regular expression used to collect info
	    
	    for (my $separatorElement=0;$separatorElement<scalar(@separators);$separatorElement++) { #Loop through possible separators 
		
		@{ $qcProgramData{$program}{$regExpKey} } = split(/$separators[$separatorElement]/, `$regExp $outDirectory/$outFile`); #Collect data. Use regExpKey as element header

		if ( defined($qcProgramData{$program}{$regExpKey}[1])) { #Then split should have been successful                                                     

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
##Arguments : $familyID, $sampleID, $program, $inFile
##          : $familyID => FamilyID
##          : $sampleID => SampleID
##          : $program  => The program to examine 
##          : $inFile   => infile to program
    
    my $familyID = $_[0]; 
    my $sampleID = $_[1];  
    my $program = $_[2]; 
    my $infile = $_[3];
    
    for my $regExpKey ( keys %{ $regExpFile{$program} } ) { #All regExp per program 
	
	if ($regExpKey !~/^header|header$/i) { #For info contained in Entry --> Value i.e. same line 

	    if (scalar(@{ $qcProgramData{$program}{$regExpKey} }) == 1) { #Enable seperation of writing array or key-->value in qcData
		
		if ( ($familyID) && ($sampleID) && ($infile) ) {
		    
		    $qcData{$familyID}{$sampleID}{$infile}{$program}{$regExpKey} = $qcProgramData{$program}{$regExpKey}[0]; #key-->value for sampleID
		}
		elsif ($familyID) {
		    
		    $qcData{$familyID}{$familyID}{'Program'}{$program}{$regExpKey} = $qcProgramData{$program}{$regExpKey}[0]; #key-->value for familyID
		}
		if ($program eq "ChanjoSexCheck") {#Check gender for sampleID
		    
		    my $ChanjoSexCheck = @{$qcProgramData{$program}{$regExpKey}}[0]; #ArrayRef
		    &GenderCheck(\$familyID, \$sampleID, \$infile, \$ChanjoSexCheck); #Check that assumed gender is supported by coverage on chrX and chrY
		}
	    }
	    else { #Write array to qcData

		for (my $regExpKeyCounter=0;$regExpKeyCounter<scalar(@{ $qcProgramData{$program}{$regExpKey} });$regExpKeyCounter++ ) {
		    
		    if ( ($familyID) && ($sampleID) && ($infile) ) {
		
			$qcData{$familyID}{$sampleID}{$infile}{$program}{$regExpKey}[$regExpKeyCounter] = $qcProgramData{$program}{$regExpKey}[$regExpKeyCounter];

		    }
		    elsif ($familyID) {

			$qcData{$familyID}{$familyID}{'Program'}{$program}{$regExpKey}[$regExpKeyCounter] = $qcProgramData{$program}{$regExpKey}[$regExpKeyCounter];			
		    }
		    if ($program eq "SexCheck") {#Check gender for sampleID
			
			my @sexChecks = split(":", @{$qcProgramData{$program}{$regExpKey}}[$regExpKeyCounter]); #ArrayRef
			&PlinkGenderCheck(\$familyID, \$sexChecks[0], \$sexChecks[1]); #Check that assumed gender is supported by variants on chrX and chrY
		    }
		}
		if (defined($qcData{$familyID}{$familyID}{'Program'}{'RelationCheck'}{'Sample_RelationCheck'}) && (defined($qcData{$familyID}{$familyID}{'Program'}{'pedigreeCheck'}{'Sample_order'}) ) ) {
		    
		    &RelationCheck(\$familyID, \@{$qcData{$familyID}{$familyID}{'Program'}{'RelationCheck'}{'Sample_RelationCheck'}}, \@{$qcData{$familyID}{$familyID}{'Program'}{'pedigreeCheck'}{'Sample_order'}});
		    delete($qcData{$familyID}{$familyID}{'Program'}{'RelationCheck'}); #Not of any use anymore
		}
	    }
	}
	else { #Paragraf data i.e. header and subsequent data lines
	    
	    for my $regExpHeaderKey ( keys %{ $qcHeader{$program}{$regExpKey} }) { #Find Header info
		
		for my $regExpKeyHeader ( keys %{ $regExpFile{$program}{$regExpKey} } ) { #All paragraf keys (header and data line(s))
		    
		    if ($regExpKeyHeader !~/^header|header$/i) { #Detect if the regExp id for headers and not data. 
			
			for (my $qcHeadersCounter=0;$qcHeadersCounter<scalar( @{ $qcHeader{$program}{$regExpKey}{$regExpHeaderKey} } );$qcHeadersCounter++) { #For all collected headers
                            
			    if ( ($familyID) && ($sampleID) && ($infile)) {
				
				$qcData{$familyID}{$sampleID}{$infile}{$program}{$regExpHeaderKey}{$regExpKeyHeader}{ $qcHeader{$program}{$regExpKey}{$regExpHeaderKey}[$qcHeadersCounter] } = $qcProgramData{$program}{$regExpKey}{$regExpKeyHeader}[$qcHeadersCounter]; #Add to qcData using header element[X] --> data[X] to correctly position elements in qcData hash 
                            } 
                            elsif ($familyID) {
				
				$qcData{$familyID}{$familyID}{$program}{$regExpHeaderKey}{$regExpKeyHeader}{ $qcHeader{$program}{$regExpKey}{$regExpHeaderKey}[$qcHeadersCounter] } = $qcProgramData{$program}{$regExpKey}{$regExpKeyHeader}[$qcHeadersCounter]; #Add to qcData using header element[X] --> data[X] to correctly position elements in qcData hash
				
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
##Arguments : $familyID
##          : $familyID => FamilyID

    my $familyID = $_[0];

    $evaluateMetric{"MosaikAligner"}{"Total_aligned"}{'threshold'} = 95;
    $evaluateMetric{"MosaikAligner"}{"Uniquely_aligned_mates"}{'threshold'} = 90;
    $evaluateMetric{"BamStats"}{"percentag_mapped_reads"}{'threshold'} = 95;
    $evaluateMetric{"CalculateHsMetrics"}{"PCT_TARGET_BASES_10X"}{'threshold'} = 0.95;
    $evaluateMetric{"CollectMultipleMetrics"}{"PCT_PF_READS_ALIGNED"}{'threshold'} = 0.95;
    $evaluateMetric{"CalculateHsMetrics"}{"PCT_ADAPTER"}{'threshold'} = 0.0001;

    if ($sampleInfoFile{$familyID}{$familyID}{'AnalysisType'} eq "exomes") {

	$evaluateMetric{"CalculateHsMetrics"}{"MEAN_TARGET_COVERAGE"}{'threshold'} = 100;
	$evaluateMetric{"CalculateHsMetrics"}{"PCT_TARGET_BASES_30X"}{'threshold'} = 0.90;
    }
    else {

	$evaluateMetric{"CalculateHsMetrics"}{"MEAN_TARGET_COVERAGE"}{'threshold'} = 20;
    }    
}
sub EvaluateQCParameters {

##EvaluateQCParameters

##Function  : Evaluate parameters to detect parameters falling below threshold 
##Returns   : ""
##Arguments : 

    my $status;

    for my $familyID ( keys %qcData ) {

	for my $ID ( keys %{$qcData{$familyID}} ) { #Can be both sampleID and familyID with current structure

	    for my $infile ( keys %{$qcData{$familyID}{$ID}} ) {
		
		if ($infile =~/RelationCheck/) { #Special case
		  
		    if ($qcData{$familyID}{$ID}{$infile} ne "PASS") {

			$status = "Status:".$infile.":".$qcData{$familyID}{$ID}{$infile};
			push(@{$qcData{$familyID}{$familyID}{'Evaluation'}{$infile}}, $status); #Add to QC data at family level
		    }
		    next;
		}
		if ($infile =~/Evaluation/) { #Special case
		    
		    next;
		}
		for my $program ( keys %{$qcData{$familyID}{$ID}{$infile}} ) {

		    if (defined($evaluateMetric{$program})) { #Program to be evaluated
	
			for my $metric ( keys %{$evaluateMetric{$program}}) { #Metric to be evaluated

			    if (defined($qcData{$familyID}{$ID}{$infile}{$program}{$metric})) {

				if ($qcData{$familyID}{$ID}{$infile}{$program}{$metric} < $evaluateMetric{$program}{$metric}{'threshold'}) { #Determine status - if below add to hash. otherwise PASS and do not include
				    
				    $status = "FAILED:".$ID."_".$program."_".$metric.":".$qcData{$familyID}{$ID}{$infile}{$program}{$metric};
				    push(@{$qcData{$familyID}{$familyID}{'Evaluation'}{$program}}, $status);
				}		
				last;
			    }
			    else {

				for my $key ( keys %{$qcData{$familyID}{$ID}{$infile}{$program}} ) {
				    
				    if ($key eq "Header") {
					
					for my $dataHeader ( keys %{$qcData{$familyID}{$ID}{$infile}{$program}{$key}} ) {
				
					    if (defined($qcData{$familyID}{$ID}{$infile}{$program}{$key}{$dataHeader}{$metric})) {
						
						if ($qcData{$familyID}{$ID}{$infile}{$program}{$key}{$dataHeader}{$metric} < $evaluateMetric{$program}{$metric}{'threshold'}) { #Determine status - if below add to hash. otherwise PASS and do not include
	
						    $status = "FAILED:".$ID."_".$program."_".$metric.":".$qcData{$familyID}{$ID}{$infile}{$program}{$key}{$dataHeader}{$metric};
						    push(@{$qcData{$familyID}{$familyID}{'Evaluation'}{$program}}, $status);
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
}


sub RelationCheck {
##Uses the .mibs file produced by PLINK to test if family members are indeed related.

    my $familyIDRef = $_[0]; #From SampleInfo
    my $relationshipValuesRef = $_[1]; #All relationship estimations
    my $sampleOrderRef = $_[2]; #The sample order so that correct estimation can be connected to the correct sampleIDs

    my %family; #Stores family relations and pairwise comparisons family{$sampleID}{$sampleID}["column"] -> [pairwise]
    my $sampleIDCounter = 0;
    my $incorrectRelation=0;
    my @pairwiseComparisons;

    ## Splice all relationship extimations from regExp into pairwise comparisons calculated for each sampleID
    for (my $realtionshipCounter=0;$realtionshipCounter<scalar(@{$relationshipValuesRef});$realtionshipCounter++) {
	
	my @pairwiseComparisons = splice(@{$relationshipValuesRef},0,scalar(@{$sampleOrderRef})); #Splices array into each sampleIDs line
	
	for (my $column=0;$column<scalar(@{$sampleOrderRef});$column++) { #All columns in .mibs file
	    
	    push ( @{ $family{ @{$sampleOrderRef}[$sampleIDCounter] }{ @{$sampleOrderRef}[$column]} }, $pairwiseComparisons[$column]); #Store sampleID, family membersID (including self) and each pairwise comparison. Uses array for to accomodate sibling info.
	}
	$sampleIDCounter++;
    }
    my $fatherID = "YYY"; #fatherID for the family
    my $motherID = "XXX"; #motherID for the family

    for my $sampleID ( keys %family ) { #For all sampleIDs

	## Currently only 1 father or Mother per pedigree is supported

	if ($sampleInfoFile{$$familyIDRef}{$sampleID}{'Father'} ne 0) { #Save fatherID if not 0

	    $fatherID = $sampleInfoFile{$$familyIDRef}{$sampleID}{'Father'};
	}
	if ($sampleInfoFile{$$familyIDRef}{$sampleID}{'Mother'} ne 0) { #Save motherID if not 0

	    $motherID = $sampleInfoFile{$$familyIDRef}{$sampleID}{'Mother'};
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
			$qcData{$$familyIDRef}{$sampleID}{'RelationCheck'} = "FAIL: Duplicated sample?;";
			#print  "Incorrect should be self: ".$sampleID,"\t", $members, "\t", $family{$sampleID}{$members}[$membersCount], "\n";
		    }
		}
		elsif ($family{$sampleID}{$members}[$membersCount] >= 0.63 ) { #Should include parent to child and child to siblings unless inbreed parents

		    if ( ( ($sampleID ne $fatherID) && ($sampleID ne $motherID) ) || ( ($members ne $fatherID) && ($members ne $motherID) ) ) { #Correct
			#print "Parent-to-child or child-to-child: ".$sampleID,"\t", $members, "\t", $family{$sampleID}{$members}[$membersCount], "\n";
		    }
		    else {

			$incorrectRelation++;
			$qcData{$$familyIDRef}{$sampleID}{'RelationCheck'} = "FAIL: Parents related?;";
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
			$qcData{$$familyIDRef}{$sampleID}{'RelationCheck'} = "FAIL:".$sampleID." not related to ".$members.";";
			#print "Incorrect: ".$sampleID,"\t", $members, "\t", $family{$sampleID}{$members}[$membersCount], "\n";
		    }
		}
	    }
	}	
	if ($incorrectRelation == 0) {
	    $qcData{$$familyIDRef}{$sampleID}{'RelationCheck'} = "PASS";  
	}
    }   
    return;
}

sub GenderCheck {
#Checks that the gender predicted by ChanjoSexCheck is confirmed in the pedigee for the sample
    
    my $familyIDRef = $_[0]; #From SampleInfo
    my $sampleIDRef = $_[1]; #From SampleInfo 
    my $infileRef = $_[2]; #From SampleInfo
    my $chanjoSexCheckGenderRef = $_[3]; #From ChanjoSexCheck 
    
    if ( ($$chanjoSexCheckGenderRef eq "female") && ($sampleInfoFile{$$familyIDRef}{$$sampleIDRef}{'Sex'} eq 2) ) { #Female
	
	$qcData{$$familyIDRef}{$$sampleIDRef}{$$infileRef}{'GenderCheck'} = "PASS";
    }
    elsif ( ($$chanjoSexCheckGenderRef eq "male") && ($sampleInfoFile{$$familyIDRef}{$$sampleIDRef}{'Sex'} eq 1) ) { #Male
	
	$qcData{$$familyIDRef}{$$sampleIDRef}{$$infileRef}{'GenderCheck'} = "PASS";
    }
    else {
	
	$qcData{$$familyIDRef}{$$sampleIDRef}{$$infileRef}{'GenderCheck'} = "FAIL";
    }
    return;
}


sub PlinkGenderCheck {
    
    #Checks that the gender predicted by Plink SexCheck is confirmed in the pedigee for the sample
    
    my $familyIDRef = $_[0]; #From SampleInfo
    my $sampleIDRef = $_[1]; #From plink sex check output 
    my $plinkSexCheckGenderRef = $_[2]; #From Plink SexCheck 
    
    if ( ($$plinkSexCheckGenderRef eq "2") && ($sampleInfoFile{$$familyIDRef}{$$sampleIDRef}{'Sex'} eq 2) ) { #Female
	
	push(@{$qcData{$$familyIDRef}{$$familyIDRef}{'Program'}{'PlinkGenderCheck'}}, $$sampleIDRef.":PASS");
    }
    elsif ( ($$plinkSexCheckGenderRef eq "1") && ($sampleInfoFile{$$familyIDRef}{$$sampleIDRef}{'Sex'} eq 1) ) { #Male
	
	push(@{$qcData{$$familyIDRef}{$$familyIDRef}{'Program'}{'PlinkGenderCheck'}}, $$sampleIDRef.":PASS");
    }
    else {
	
	push(@{$qcData{$$familyIDRef}{$$familyIDRef}{'Program'}{'PlinkGenderCheck'}}, $$sampleIDRef.":FAIL");
    }
    return;
}

sub WriteYAML {
###Writes a YAML hash to file. 
###Note: 2nd argument should be a hash reference

    my $yamlFile = $_[0]; #Filename
    my $yamlHashRef = $_[1]; #Hash reference to write to file

    open (YAML, ">". $yamlFile) or die "can't open ".$yamlFile.": $!\n";
    print YAML Dump( $yamlHashRef ), "\n";
    close(YAML);
}

sub LoadYAML {
###Loads a YAML file into an arbitrary hash and returns it. Note: Currently only supports hashreferences and hashes and no mixed entries 

    my $yamlFile = $_[0];
    my %yamlHash;

    my $fileType = &DetectYamlContentType($yamlFile);

    open (YAML, "<". $yamlFile) or die "can't open ".$yamlFile.": $!\n";    
        
        if ($fileType eq "reference") {

	    %yamlHash = %{ YAML::LoadFile($yamlFile) }; #Load hashreference as hash
        }
        if ($fileType eq "hash") {
        
	    %yamlHash = YAML::LoadFile($yamlFile); #File contained a hash = no workup
        }
    close(YAML);

    return %yamlHash;
}

sub DetectYamlContentType {
###Check the content of the YAML file for seperating hashreferences and hash. Return the content type.

    my $yamlFile = $_[0];
    my $fileType;

    open (YAML, "<". $yamlFile) or die "can't open ".$yamlFile.": $!\n";
        
        while (<YAML>) {

            if ($. == 1 && $_=~/^---$/) { #YAML file contains a hashreference
                $fileType = "reference";
                last;
            }
            else {
                $fileType = "hash";
                last;
            }
        }
    close(YAML);
    return $fileType;
}

sub RegExpToYAML {

    my $printRegExpOutFile = $_[0]; #File to print regExpe to

    my %regExp;
    #Add to %regExp to enable print in YAML

    $regExp{'FastQC'}{'Version'} = q?perl -nae' if ($_=~/##FastQC\\s+(\\S+)/) {print $1;last;}' ?; #Collect FastQC version

    $regExp{'FastQC'}{'Encoding'} = q?perl -nae' if ($_=~/Encoding\s+(\S+\s\S+\s\S+\s\S+|\S+\s\S+)/) { my $encoding = $1;$encoding=~s/\s/\_/g; print $encoding;last;}' ?; #Collect Encoding
    
    $regExp{'FastQC'}{'Sequence_length'} = q?perl -nae' if ($_=~/Sequence length\s(\d+)/) {print $1;last;}' ?; #Collect Sequence length
    
    $regExp{'FastQC'}{'Total_number_of_reads'} = q?perl -nae' if ($_=~/Total Sequences\s(\d+)/) {print $1;last;}' ?; #Collect Total sequences 
    
    $regExp{'FastQC'}{'GC'} = q?perl -nae' if ($_=~/%GC\s(\d+)/) {print $1;last;}' ?; #Collect GC content 
    
    $regExp{'FastQC'}{'Sequence_duplication'} = q?perl -nae' if ($_=~/#Total Duplicate Percentage\s+(\d+.\d)/) {print $1;last;}' ?; #Collect Sequence duplication level
    
    $regExp{'FastQC'}{'Basic_statistics'} = q?perl -nae' if ($_=~/>>Basic Statistics\s+(\S+)/) {print $1;last;}' ?; #Collect Basic Statistics
    
    $regExp{'FastQC'}{'Per_base_sequence_quality'} = q?perl -nae' if ($_=~/>>Per base sequence quality\s+(\S+)/) {print $1;last;}' ?; #Collect Per base sequence quality
    
    $regExp{'FastQC'}{'Per_sequence_quality_scores'} = q?perl -nae' if ($_=~/>>Per sequence quality scores\s+(\S+)/) {print $1;last;}' ?; #Collect Per sequence quality scores
    
    $regExp{'FastQC'}{'Per_base_sequence_content'} = q?perl -nae' if ($_=~/>>Per base sequence content\s+(\S+)/) {print $1;last;}' ?; #Collect Per base sequence content
    
    $regExp{'FastQC'}{'Per_base_GC_content'} = q?perl -nae' if ($_=~/>>Per base GC content\s+(\S+)/) {print $1;last;}' ?; #Collect Per base GC content
    
    $regExp{'FastQC'}{'Per_sequence_GC_content'} = q?perl -nae' if ($_=~/>>Per sequence GC content\s+(\S+)/) {print $1;last;}' ?; #Collect Per sequence GC content
    
    $regExp{'FastQC'}{'Per_base_N_content'} = q?perl -nae' if ($_=~/>>Per base N content\s+(\S+)/) {print $1;last;}' ?; #Collect Per base N content
    
    $regExp{'FastQC'}{'Sequence_duplication_levels'} = q?perl -nae' if ($_=~/>>Sequence Duplication Levels\s+(\S+)/) {print $1;last;}' ?; #Collect Sequence Duplication Levels
    
    $regExp{'FastQC'}{'Overrepresented_sequences'} = q?perl -nae' if ($_=~/>>Overrepresented sequences\s+(\S+)/) {print $1;last;}' ?; #Collect Overrepresented sequences
    
    $regExp{'FastQC'}{'Kmer_content'} = q?perl -nae' if ($_=~/>>Kmer Content\s+(\S+)/) {print $1;last;}' ?; #Collect Kmer Content
    
    $regExp{'MosaikAligner'}{'Version'} = q?perl -nae' if ($_=~/(\d+\.\d+\.\d+)\s/) {print $1;last;}' ?; #Collect Mosaik Version 
    
    $regExp{'MosaikAligner'}{'Unaligned_mates'} = q?perl -nae' if ($_=~/# unaligned mates\S+\s+(\d+)\s\(\s+(\d+\.\d+)/) {print $2;last;}' ?; #Collect Nr of unaligned mates
    
    $regExp{'MosaikAligner'}{'Filtered_out'} = q?perl -nae' if ($_=~/# filtered out\S+\s+(\d+)\s\(\s+(\d+\.\d+)/) {print $2;last;}' ?; #Collect Nr of filtered out reads 
    
    $regExp{'MosaikAligner'}{'Uniquely_aligned_mates'} = q?perl -nae' if ($_=~/# uniquely aligned mates\S+\s+(\d+)\s\(\s+(\d+\.\d+)/) {print $2;last;}' ?; #Collect Uniquely aligned mates
    
    $regExp{'MosaikAligner'}{'Multiply_aligned_mates'} = q?perl -nae' if ($_=~/# multiply aligned mates\S+\s+(\d+)\s\(\s+(\d+\.\d+)/) {print $2;last;}' ?; #Collect Multiply aligned mates
    
    $regExp{'MosaikAligner'}{'Total_aligned'} = q?perl -nae' if ($_=~/total aligned:\s+\S+\s+(\S+)\s\(\S+\s(\d+.\d+)/ ) {print $2;last;} elsif ($_=~/total aligned:\s+(\S+)\s\(\S+\s(\d+.\d+)/ ) { print $2;last;}' ?; #Collect total aligned sequences

    $regExp{'BamStats'}{'percentag_mapped_reads'} = q?perl -nae 'if($_=~/percentag mapped reads:\s+(\S+)/) {print $1;last}' ?; #Collect % mapped reads from BAm alignment

    $regExp{'ChanjoSexCheck'}{'gender'} = q?perl -nae 'if( ($F[0]!~/^#/) && ($F[2] =~/\S+/) ) {print $F[2];}' ?;  #Collect gender from ChanjoSexCheck

    $regExp{'pedigreeCheck'}{'Sample_order'} = q?perl -nae 'if ($_=~/^#CHROM/) {chomp $_; my @line = split(/\t/,$_); for (my $sample=9;$sample<scalar(@line);$sample++) { print $line[$sample], "\t";}last;}' ?; #Collect sample order from vcf file used to create ".ped", ".map" and hence ".mibs".
    
    $regExp{'InbreedingFactor'}{'Sample_InbreedingFactor'}  = q?perl -nae 'my @inbreedingFactor; if ($. > 1) {my @temp = split(/\s/,$_);push(@inbreedingFactor,$temp[0].":".$temp[4]); print $inbreedingFactor[0], "\t"; }' ?;

    $regExp{'SexCheck'}{'Sample_SexCheck'}  = q?perl -nae 'my @sexCheckFactor; if ($. > 1) {my @temp = split(/\s+/,$_);push(@sexCheckFactor,$temp[2].":".$temp[4]); print $sexCheckFactor[0], "\t"; }' ?;

    $regExp{'RelationCheck'}{'Sample_RelationCheck'}  = q?perl -nae 'print $_;' ?; #Note will return whole file

    $regExp{'MarkDuplicates'}{'Fraction_duplicates'} = q?perl -nae 'if($_=~/Fraction Duplicates\: (\S+)/) {print $1;}' ?; #Collect fraction duplicates
    
    $regExp{'CalculateHsMetrics'}{'Header_info'}{'Header'} = q?perl -nae' if ($_ =~/^BAIT_SET/ ) {print $_;last;}' ?; #Note return whole line (Header) 
    
    $regExp{'CalculateHsMetrics'}{'Header_info'}{'Data'} = q?perl -nae' if ( ($. ==8) && ($_ =~/(\S+)/) ) {print $_;last;}' ?; #Note return whole line and only look at line 8, where the data action is
    
    $regExp{'CollectMultipleMetrics'}{'Header_info'}{'Header'} = q?perl -nae' if ($_ =~/^CATEGORY/ ) {print $_;last;}' ?; #Note return whole line (Header) 
    
    $regExp{'CollectMultipleMetrics'}{'Header_info'}{'First_of_pair'} = q?perl -nae' if ($_ =~/^FIRST_OF_PAIR/ ) {print $_;last;}' ?; #Note return whole line (FIRST_OF_PAIR)
    
    $regExp{'CollectMultipleMetrics'}{'Header_info'}{'Second_of_pair'} = q?perl -nae' if ($_ =~/^SECOND_OF_PAIR/ ) {print $_;last;}' ?; #Note return whole line (SECOND_OF_PAIR)  
    
    $regExp{'CollectMultipleMetrics'}{'Header_info'}{'Pair'} = q?perl -nae' if ($_ =~/^PAIR/ ) {print $_;last;}'  ?; #Note return whole line (PAIR)
    
    $regExp{'VariantEval_All'}{'CompOverlap_header'}{'CompOverlap_header'} = q?perl -nae' if ($_ =~/^CompOverlap\s+\CompRod/ ) {print $_;last;}' ?; #Note return whole line (Header)
    
    $regExp{'VariantEval_All'}{'CompOverlap_header'}{'CompOverlap_data_all'} = q?perl -nae' if ( ($_ =~/^CompOverlap/) && ($_ =~/all/) && ($_ =~/none/)) {print $_;last;}' ?; #Note return whole line
    
    $regExp{'VariantEval_All'}{'CompOverlap_header'}{'CompOverlap_data_known'} = q?perl -nae' if ( ($_ =~/^CompOverlap/) && ($_ =~/known\s/) ) {print $_;last;}' ?; #Note return whole line

    $regExp{'VariantEval_All'}{'CompOverlap_header'}{'CompOverlap_data_novel'} = q?perl -nae' if ( ($_ =~/^CompOverlap/) && ($_ =~/novel\s/) ) {print $_;last;}' ?; #Note return whole line

    $regExp{'VariantEval_All'}{'CountVariants_header'}{'CountVariants_header'} = q?perl -nae' if ($_ =~/^CountVariants\s+\CompRod/ ) {print $_;last;}' ?; #Note return whole line (Header)                                               
    $regExp{'VariantEval_All'}{'CountVariants_header'}{'CountVariants_data_all'} = q?perl -nae' if ( ($_ =~/^CountVariants/) && ($_ =~/all\s/) ) {print $_;last;}' ?; #Note return whole line                                              
    $regExp{'VariantEval_All'}{'CountVariants_header'}{'CountVariants_data_known'} = q?perl -nae' if ( ($_ =~/^CountVariants/) && ($_ =~/known\s/) ) {print $_;last;}' ?; #Note return whole line                                          
    $regExp{'VariantEval_All'}{'CountVariants_header'}{'CountVariants_data_novel'} = q?perl -nae' if ( ($_ =~/^CountVariants/) && ($_ =~/novel\s/) ) {print $_;last;}' ?; #Note return whole line

    $regExp{'VariantEval_All'}{'IndelSummary_header'}{'IndelSummary_header'} = q?perl -nae' if ($_ =~/^IndelSummary\s+\CompRod/ ) {print $_;last;}' ?; #Note return whole line (Header)                                                                                                                                                                                                                         
    $regExp{'VariantEval_All'}{'IndelSummary_header'}{'IndelSummary_data_all'} = q?perl -nae' if ( ($_ =~/^IndelSummary/) && ($_ =~/all\s/) ) {print $_;last;}' ?; #Note return whole line                                                                                                                                                                                                        
    $regExp{'VariantEval_All'}{'IndelSummary_header'}{'IndelSummary_data_known'} = q?perl -nae' if ( ($_ =~/^IndelSummary/) && ($_ =~/known\s/) ) {print $_;last;}' ?; #Note return whole line                                     
    $regExp{'VariantEval_All'}{'IndelSummary_header'}{'IndelSummary_data_novel'} = q?perl -nae' if ( ($_ =~/^IndelSummary/) && ($_ =~/novel\s/) ) {print $_;last;}' ?; #Note return whole line 

    $regExp{'VariantEval_All'}{'MultiallelicSummary_header'}{'MultiallelicSummary_header'} = q?perl -nae' if ($_ =~/^MultiallelicSummary\s+\CompRod/ ) {print $_;last;}' ?; #Note return whole line (Header)                         
    $regExp{'VariantEval_All'}{'MultiallelicSummary_header'}{'MultiallelicSummary_data_all'} = q?perl -nae' if ( ($_ =~/^MultiallelicSummary/) && ($_ =~/all\s/) ) {print $_;last;}' ?; #Note return whole line                                          
    $regExp{'VariantEval_All'}{'MultiallelicSummary_header'}{'MultiallelicSummary_data_known'} = q?perl -nae' if ( ($_ =~/^MultiallelicSummary/) && ($_ =~/known\s/) ) {print $_;last;}' ?; #Note return whole line                                        
    $regExp{'VariantEval_All'}{'MultiallelicSummary_header'}{'MultiallelicSummary_data_novel'} = q?perl -nae' if ( ($_ =~/^MultiallelicSummary/) && ($_ =~/novel\s/) ) {print $_;last;}' ?; #Note return whole line

    $regExp{'VariantEval_All'}{'TiTvVariantEvaluator_header'}{'TiTvVariantEvaluator_header'} = q?perl -nae' if ($_ =~/^TiTvVariantEvaluator\s+\CompRod/ ) {print $_;last;}' ?; #Note return whole line (Header)                                                                                                                                                                                        
    $regExp{'VariantEval_All'}{'TiTvVariantEvaluator_header'}{'TiTvVariantEvaluator_data_all'} = q?perl -nae' if ( ($_ =~/^TiTvVariantEvaluator/) && ($_ =~/all\s/) ) {print $_;last;}' ?; #Note return whole line                                                                                                                                                                                                       
    $regExp{'VariantEval_All'}{'TiTvVariantEvaluator_header'}{'TiTvVariantEvaluator_data_known'} = q?perl -nae' if ( ($_ =~/^TiTvVariantEvaluator/) && ($_ =~/known\s/) ) {print $_;last;}' ?; #Note return whole line                                                                                                                                                                                               
    $regExp{'VariantEval_All'}{'TiTvVariantEvaluator_header'}{'TiTvVariantEvaluator_data_novel'} = q?perl -nae' if ( ($_ =~/^TiTvVariantEvaluator/) && ($_ =~/novel\s/) ) {print $_;last;}' ?; #Note return whole line 

    $regExp{'VariantEval_All'}{'ValidationReport_header'}{'ValidationReport_header'} = q?perl -nae' if ($_ =~/^ValidationReport\s+\CompRod/ ) {print $_;last;}' ?; #Note return whole line (Header)                                                                                                                                                                                                                          
    $regExp{'VariantEval_All'}{'ValidationReport_header'}{'ValidationReport_data_all'} = q?perl -nae' if ( ($_ =~/^ValidationReport/) && ($_ =~/all\s/) && ($_ =~/none\s/)) {print $_;last;}' ?; #Note return whole line                                                                                                                                                                                                                         
    $regExp{'VariantEval_All'}{'ValidationReport_header'}{'ValidationReport_data_known'} = q?perl -nae' if ( ($_ =~/^ValidationReport/) && ($_ =~/known\s/) ) {print $_;last;}' ?; #Note return whole line                                                                                                                                                                                                        
    $regExp{'VariantEval_All'}{'ValidationReport_header'}{'ValidationReport_data_novel'} = q?perl -nae' if ( ($_ =~/^ValidationReport/) && ($_ =~/novel\s/) ) {print $_;last;}' ?; #Note return whole line

    $regExp{'VariantEval_All'}{'VariantSummary_header'}{'VariantSummary_header'} = q?perl -nae' if ($_ =~/^VariantSummary\s+\CompRod/ ) {print $_;last;}' ?; #Note return whole line (Header)                                                                                                                                                                                                                     
    $regExp{'VariantEval_All'}{'VariantSummary_header'}{'VariantSummary_data_all'} = q?perl -nae' if ( ($_ =~/^VariantSummary/) && ($_ =~/all\s/) ) {print $_;last;}' ?; #Note return whole line                                                                                                                                                                                                                   
    $regExp{'VariantEval_All'}{'VariantSummary_header'}{'VariantSummary_data_known'} = q?perl -nae' if ( ($_ =~/^VariantSummary/) && ($_ =~/known\s/) ) {print $_;last;}' ?; #Note return whole line                                                                                                                                                                                                                
    $regExp{'VariantEval_All'}{'VariantSummary_header'}{'VariantSummary_data_novel'} = q?perl -nae' if ( ($_ =~/^VariantSummary/) && ($_ =~/novel\s/) ) {print $_;last;}' ?; #Note return whole line

    $regExp{'VariantEval_Exome'} = $regExp{'VariantEval_All'};

    $regExp{'Genmod'}{'Version'} = q?perl -nae 'if($_=~/##Software=<ID=genmod\w+,Version=(\d+.\d+.\d+)/) {print $1;last;}' ?; #Collect Genmod version

    $regExp{'SnpEff'}{'Version'} = q?perl -nae 'if($_=~/##SnpSiftVersion=\"(.+),/) {my $ret=$1; $ret=~s/\s/_/g;print $ret;last;}' ?; #Collect SnpEff version

    $regExp{'VariantEffectPredictor'}{'Version'} = q?perl -nae 'if($_=~/##VEP=(\w+)/) {print $1;last;}' ?; #Collect VariantEffectPredictor version

    $regExp{'VariantEffectPredictor'}{'Cache:'} = q?perl -nae 'if($_=~/##VEP=\w+\s+cache=(\S+)/) {print $1;last;}' ?; #Collect VariantEffectPredictor cache directory 

    $regExp{'VCFParser'}{'Version'} = q?perl -nae 'if($_=~/##Software=<ID=vcfParser.pl,Version=(\d+.\d+.\d+)/) {print $1;last;}' ?; #Collect VCFParser version

    $regExp{'Bwa'}{'Version'} = q?perl -nae 'if($_=~/\[main\]\sVersion:\s(\S+)/) {print $1;last;}' ?; #Collect Bwa version

    $regExp{'Chanjo'}{'Version'} = q?perl -nae 'if($_=~/version\s(\d+.\d+.\d+)/) {print $1;last;}' ?; #Collect Chanjo version

    $regExp{'vt'}{'Version'} = q?perl -nae 'if($_=~/decompose\sv(\S+)/) {print $1;last;}' ?; #Collect vt version

    $regExp{'Samtools'}{'Version'} = q?perl -nae 'if($_=~/samtoolsVersion=(\S+)/) {print $1;last;}' ?; #Collect Samtools version

    $regExp{'Bcftools'}{'Version'} = q?perl -nae 'if($_=~/bcftools_\w+Version=(\S+)/) {print $1;last;}' ?; #Collect Bcftools version

    $regExp{'Freebayes'}{'Version'} = q?perl -nae 'if($_=~/source=freeBayes\s(\S+)/) {print $1;last;}' ?; #Collect Freebayes version

#$regExp{''}{''} = ;

    
&WriteYAML($printRegExpOutFile, \%regExp);

}
