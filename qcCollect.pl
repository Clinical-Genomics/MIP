#!/usr/bin/perl - w                                                                                                                                                                                       
use strict;
use warnings;

#Collects MPS QC from MIP. Loads information on files to examine and values to extract from in YAML format and outputs exracted metrics in YAML format. 
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

my ($sampleInfoFile, $regExpFile, $outfile, $printRegExp, $printRegExpOutFile, $version, $help) = (0,0,"qcmetrics.yaml", 0, "qc_regExp.yaml");
my (%sampleInfo, %regExp, %qcData);
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

die $USAGE if( $help );

die "\nqcCollect.pl v1.0\n\n" if($version);

if ($printRegExp) {

    RegExpToYAML($printRegExpOutFile); #Write regExp used @ CMMS to YAML
    print STDOUT "Wrote RegExp YAML file to: ".$printRegExpOutFile, "\n";
    exit;
}

if ($sampleInfoFile eq 0) {

    print STDERR "\n";
    print STDERR "Must supply a '-sampleInfoFile' (supply whole path)", "\n\n";
    die $USAGE;
}
if ($regExpFile eq 0) {

    print STDERR "\n";
    print STDERR "Must supply a '-regExpFile' (supply whole path)", "\n\n";
    die $USAGE;
}

####MAIN

my %sampleInfoFile = LoadYAML($sampleInfoFile); #Load sampleInfoFile (YAML) and transfer to sampleInfoFile (Hash)

my %regExpFile = LoadYAML($regExpFile); #Load regExpFile (YAML) and transfer to regExpFile (Hash)

SampleQC(); #Extracts all qcdata on sampleID level using information in %sampleInfoFile and %regExpFile

FamilyQC(); #Extracts all qcdata on family level using information in %sampleInfoFile and %regExpFile

WriteYAML($outfile, \%qcData ); #Writes to YAML file

####SubRoutines

sub FamilyQC {

    for my $familyID ( keys %sampleInfoFile ) { #For every family id 

        if ($sampleInfoFile{$familyID}{'program'}) { #Only examine programs     

            for my $program ( keys %{ $sampleInfoFile{$familyID}{'program'} } ) { #For every programs           

		my $outDirectory;
		my $outFile;	

                if ($sampleInfoFile{$familyID}{'program'}{$program}{'Version'} ) {
                    
                    $qcData{$familyID}{'program'}{$program}{'Version'} = $sampleInfoFile{$familyID}{'program'}{$program}{'Version'}; #Add version to qcData
                }
                if ($sampleInfoFile{$familyID}{'program'}{$program}{'OutDirectory'} ) {

                    $outDirectory = $sampleInfoFile{$familyID}{'program'}{$program}{'OutDirectory'}; #Extract OutDirectory
                }
                if ($sampleInfoFile{$familyID}{'program'}{$program}{'OutFile'} ) {

                    $outFile = $sampleInfoFile{$familyID}{'program'}{$program}{'OutFile'}; #Extract OutFile
                }  
		
                ParseRegExpHashAndCollect($program, $outDirectory, $outFile); #Loads qcHeader and qcProgramData
		
                AddToqcData($familyID, "", $program, ""); #Add extracted information to qcData
            }
        }
    }
}

sub SampleQC {
###Collects all sample qc in files defined by sampleInfoFile and regular expressions defined by regExpFile.
    
    for my $familyID ( keys %sampleInfoFile ) { #For every family id
	
	for my $sampleID ( keys %{ $sampleInfoFile{$familyID} } ) { #For every sample id
	    
	    if ($sampleInfoFile{$familyID}{$sampleID}{'program'}) { #Only examine programs
		
		for my $program ( keys %{ $sampleInfoFile{$familyID}{$sampleID}{'program'} } ) { #For every program  
		    
		    for my $infile ( keys %{ $sampleInfoFile{$familyID}{$sampleID}{'program'}{$program} } ) { #For every infile
			
			my $outDirectory = $sampleInfoFile{$familyID}{$sampleID}{'program'}{$program}{$infile}{'OutDirectory'};
			my $outFile = $sampleInfoFile{$familyID}{$sampleID}{'program'}{$program}{$infile}{'OutFile'};
			
			ParseRegExpHashAndCollect($program, $outDirectory, $outFile); #Loads qcHeader and qcProgramData
			
			AddToqcData($familyID, $sampleID, $program, $infile); #Add extracted information to qcData
			
		    }   
		}
	    }
	}
    }
}   

sub ParseRegExpHashAndCollect {
###Parses the RegExpHash structure to identify if the info is 1) Paragraf section(s) (both header and data line(s)); 2) Seperate data line

    my $program = $_[0]; #From SampleInfo
    my $outDirectory = $_[1]; #From SampleInfo
    my $outFile = $_[2]; #From SampleInfo
    
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
			    last; #Found correct separator do not continue                                                                                         
			}   
		    }
		    else { #For paragraf data line(s)                                                                                                                        
			@{ $qcProgramData{$program}{$regExpKey}{$regExpHeaderKey} } = split(/$separators[$separatorElement]/, `$regExp $outDirectory/$outFile`); #Collect paragraf data
			if ( defined($qcProgramData{$program}{$regExpKey}{$regExpHeaderKey}[1])) { #Then split should have been successful                                                                                                                           
			    last; #Found correct separator do not continue                                                                                                                                      
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
    return;
}

sub AddToqcData {
    ##Add to qcData hash to enable write to yaml format
    
    my $familyID = $_[0]; #From SampleInfo
    my $sampleID = $_[1]; #From SampleInfo 
    my $program = $_[2]; #From SampleInfo
    my $infile = $_[3]; #From SampleInfo
    
    for my $regExpKey ( keys %{ $regExpFile{$program} } ) { #All regExp per program 
	
	if ($regExpKey !~/^header|header$/i) { #For info contained in Entry --> Value i.e. same line 
	    
	    if (scalar(@{ $qcProgramData{$program}{$regExpKey} }) == 1) { #Enable seperation of writing array or key-->value in qcData
		
		if ( ($familyID) && ($sampleID) && ($infile) ) {
		    $qcData{$familyID}{$sampleID}{$infile}{$program}{$regExpKey} = $qcProgramData{$program}{$regExpKey}[0]; #key-->value for sampleID
		}
		elsif ($familyID) {
		    $qcData{$familyID}{'program'}{$program}{$regExpKey} = $qcProgramData{$program}{$regExpKey}[0]; #key-->value for familyID
		}
	    }
	    else { #Write array to qcData
		
		for (my $regExpKeyCounter=0;$regExpKeyCounter<scalar(@{ $qcProgramData{$program}{$regExpKey} });$regExpKeyCounter++ ) {
		    
		    if ( ($familyID) && ($sampleID) && ($infile) ) {
			$qcData{$familyID}{$sampleID}{$infile}{$program}{$regExpKey}[$regExpKeyCounter] = $qcProgramData{$program}{$regExpKey}[$regExpKeyCounter];
		    }
		    elsif ($familyID) {
			$qcData{$familyID}{'program'}{$program}{$regExpKey}[$regExpKeyCounter] = $qcProgramData{$program}{$regExpKey}[$regExpKeyCounter];
		    }
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
				$qcData{$familyID}{$program}{$regExpHeaderKey}{$regExpKeyHeader}{ $qcHeader{$program}{$regExpKey}{$regExpHeaderKey}[$qcHeadersCounter] } = $qcProgramData{$program}{$regExpKey}{$regExpKeyHeader}[$qcHeadersCounter]; #Add to qcData using header element[X] --> data[X] to correctly position elements in qcData hashe
				
                            }
                        }    
		    }
		}
	    }
	}
    }
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

    my $fileType = DetectYamlContentType($yamlFile);

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

    $regExp{'FastQC'}{'Encoding'} = q?perl -nae' if ($_=~/Encoding\s+(\S+\s\S+\s\S+\s\S+|\S+\s\S+)/) { my $encoding = $1;$encoding=~s/\s/\_/g; print $encoding;}' ?; #Collect Encoding
    
    $regExp{'FastQC'}{'Sequence_length'} = q?perl -nae' if ($_=~/Sequence length\s(\d+)/) {print $1;}' ?; #Collect Sequence length
    
    $regExp{'FastQC'}{'Total_number_of_reads'} = q?perl -nae' if ($_=~/Total Sequences\s(\d+)/) {print $1;}' ?; #Collect Total sequences 
    
    $regExp{'FastQC'}{'GC'} = q?perl -nae' if ($_=~/%GC\s(\d+)/) {print $1;}' ?; #Collect GC content 
    
    $regExp{'FastQC'}{'Sequence_duplication'} = q?perl -nae' if ($_=~/#Total Duplicate Percentage\s+(\d+.\d)/) {print $1;}' ?; #Collect Sequence duplication level
    
    $regExp{'FastQC'}{'Basic_statistics'} = q?perl -nae' if ($_=~/>>Basic Statistics\s+(\S+)/) {print $1;}' ?; #Collect Basic Statistics
    
    $regExp{'FastQC'}{'Per_base_sequence_quality'} = q?perl -nae' if ($_=~/>>Per base sequence quality\s+(\S+)/) {print $1;}' ?; #Collect Per base sequence quality
    
    $regExp{'FastQC'}{'Per_sequence_quality_scores'} = q?perl -nae' if ($_=~/>>Per sequence quality scores\s+(\S+)/) {print $1;}' ?; #Collect Per sequence quality scores
    
    $regExp{'FastQC'}{'Per_base_sequence_content'} = q?perl -nae' if ($_=~/>>Per base sequence content\s+(\S+)/) {print $1;}' ?; #Collect Per base sequence content
    
    $regExp{'FastQC'}{'Per_base_GC_content'} = q?perl -nae' if ($_=~/>>Per base GC content\s+(\S+)/) {print $1;}' ?; #Collect Per base GC content
    
    $regExp{'FastQC'}{'Per_sequence_GC_content'} = q?perl -nae' if ($_=~/>>Per sequence GC content\s+(\S+)/) {print $1;}' ?; #Collect Per sequence GC content
    
    $regExp{'FastQC'}{'Per_base_N_content'} = q?perl -nae' if ($_=~/>>Per base N content\s+(\S+)/) {print $1;}' ?; #Collect Per base N content
    
    $regExp{'FastQC'}{'Sequence_duplication_levels'} = q?perl -nae' if ($_=~/>>Sequence Duplication Levels\s+(\S+)/) {print $1;}' ?; #Collect Sequence Duplication Levels
    
    $regExp{'FastQC'}{'Overrepresented_sequences'} = q?perl -nae' if ($_=~/>>Overrepresented sequences\s+(\S+)/) {print $1;}' ?; #Collect Overrepresented sequences
    
    $regExp{'FastQC'}{'Kmer_content'} = q?perl -nae' if ($_=~/>>Kmer Content\s+(\S+)/) {print $1;}' ?; #Collect Kmer Content
    
    $regExp{'MosaikAligner'}{'Version'} = q?perl -nae' if ($_=~/(\d+\.\d+\.\d+)\s/) {print $1;}' ?; #Collect Mosaik Version 
    
    $regExp{'MosaikAligner'}{'Unaligned_mates'} = q?perl -nae' if ($_=~/# unaligned mates\S+\s+(\d+)\s\(\s+(\d+\.\d+)/) {print $2;}' ?; #Collect Nr of unaligned mates
    
    $regExp{'MosaikAligner'}{'Filtered_out'} = q?perl -nae' if ($_=~/# filtered out\S+\s+(\d+)\s\(\s+(\d+\.\d+)/) {print $2;}' ?; #Collect Nr of filtered out reads 
    
    $regExp{'MosaikAligner'}{'Uniquely_aligned_mates'} = q?perl -nae' if ($_=~/# uniquely aligned mates\S+\s+(\d+)\s\(\s+(\d+\.\d+)/) {print $2;}' ?; #Collect Uniquely aligned mates
    
    $regExp{'MosaikAligner'}{'Multiply_aligned_mates'} = q?perl -nae' if ($_=~/# multiply aligned mates\S+\s+(\d+)\s\(\s+(\d+\.\d+)/) {print $2;}' ?; #Collect Multiply aligned mates
    
    $regExp{'MosaikAligner'}{'Total_aligned'} = q?perl -nae' if ($_=~/total aligned:\s+\S+\s+(\S+)\s\(\S+\s(\d+.\d+)/ ) {print $2;} elsif ($_=~/total aligned:\s+(\S+)\s\(\S+\s(\d+.\d+)/ ) { print $2}' ?; #Collect total aligned sequences
    
    $regExp{'QaCompute'}{'X_Y_coverage'} = q?perl -nae' if ($F[0]=~/^X/ || $F[0]=~/^Y/ ) {print "$F[2],";}' ?; #Collect X and Y coverage. "," required for later split 
    
    $regExp{'pedigreeCheck'}{'Sample_order'} = q?perl -nae 'if ($_=~/^#CHROM/) {chomp $_; my @line = split(/\t/,$_); for (my $sample=9;$sample<scalar(@line);$sample++) { print $line[$sample], "\t";}last;}' ?; #Collect sample order from vcf file used to create ".ped", ".map" and hence ".mibs".
    
    $regExp{'InbreedingFactor'}{'Sample_InbreedingFactor'}  = q?perl -nae 'my @inbreedingFactor; if ($. > 1) {my @temp = split(/\s/,$_);push(@inbreedingFactor,$temp[0].":".$temp[4]); print $inbreedingFactor[0], "\t"; }' ?;

    $regExp{'MarkDuplicates'}{'Header_info'}{'Header'} = q?perl -nae' if ($_ =~/^LIBRARY/ ) {print $_;last;}' ?; #Note return whole line (Header) 
    
    $regExp{'MarkDuplicates'}{'Header_info'}{'Data'} = q?perl -nae' if ( ($. ==8) && ($_ =~/(\S+)/) ) {print $_;last;}' ?; #Note return whole line and only look at line 8, where the data action is               

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
#$regExp{''}{''} = ;

    
WriteYAML($printRegExpOutFile, \%regExp);

}
