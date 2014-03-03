#!/usr/bin/perl - w

use strict;
use warnings;
use File::Basename;
use IO::File;

=for comment
Intersects and collects information based on 1-4 keys present (mandatory) in
each file to be investigated. The set of elements to be interrogated are decided
by the first db file elements unless the merge option is used. The db files are
supplied using the -db flag, which should point to a db_master file (tab-sep)
with the format:
(DbPath\tSeparator\tColumn_Keys\tChr_Column\tMatching\tColumns_to_Extract\tFile_Size\t).

N.B. matching should be either range or exact. Currently the range option only
supports 3-4 keys i.e. only 3 keys are used to define the range look up
(preferbly chr, start, stop). Range db file should be sorted -k1,1 -k2,2n if it
contains chr information. If the merge option is used then all overlapping and
unique elements are added to the final list. Beware that this option is memory
demanding.
=cut

# Copyright 2012 Henrik Stranneheim

=head1 SYNOPSIS
    
intersectCollect.pl -db db_master.txt -o outFile.txt
    
=head2 COMMANDS AND OPTIONS

-db/--dbFile A tab-sep file containing 1 db per line with format (DbPath\tSeparator\tColumn_Keys\tChr_Column\tMatching\tColumns_to_Extract\tFile_Size\t). NOTE: db file and col nr are 0-based. 

-s/--sampleIDs The sample ID(s)

-o/--outFile The output file (defaults to intersectCollect.txt)

-oinfo/--outInfos The headers and order for each column in output file. The information can also be recorded in the db master file as "outinfo:header1=0_2,,headerN=N_N". NOTE: db file and col nr are 0-based. (if col 2,3 in 1st db file & col 0,1 in 0th db file is desired then ocol should be 1_2,1_3,0_0,0_1. You are free to include, exclude or rearrange the order as you like. Precedence 1. command line 2. Recorded in db master file 3. Order of appearance in db master file.

-m/--merge Merge all entries found within db files. Unique entries (not found in first db file) will be included. Do not support range matching. (Defaults to "0")

-sl/--select Select all entries in first infile matching keys in subsequent db files. Do not support range matching. (Defaults to "0") 

-sofs/--selectOutFiles Selected variants and orphan db files out data directory. Comma sep (Defaults to ".";Supply whole path(s) and in the same order as the '-db' db file) 

-prechr/--prefixChromosomes "chrX" or just "X" (defaults to "X")

-h/--help Display this help message    

-v/--version Display version

=head3 I/O

Input format db master file (tab-sep)

Output format (tab separate list)

=cut

use Pod::Usage;
use Pod::Text;
use Getopt::Long;
use Tabix;

use vars qw($USAGE);

BEGIN {
    $USAGE =
	qq{intersectCollect.pl -db db_master.txt -o outFile.txt
               -db/--dbFile A tab-sep file containing 1 db per line with format (DbPath\tSeparator\tColumn_Keys\tChr_Column\tMatching\tColumns_to_Extract\tFile_Size\t). NOTE: db file and col nr are 0-based.
                  1. DbPath = Complete path to db file. First file determines the nr of elements that are subsequently used to collects information. [string]
                  2. Separator = Anything that can be inserted into perls split function. [string]
                  3. Column_Keys = The column number(s) for the key values that are to be matched (1-4 keys supported). [Number]
                  4. Chr_Column = The column number(s) for the chr information. Set to "Na" if not applicable and the program will not try to use chr information. [Number, Na]
                  5. Columns_to_Extract = The column number(s) for the information to extract from the db file(s). [Number]
                  6. Matching = The type of matching to apply to the keys. Range db file should be sorted -k1,1 -k2,2n if it contains chr information. Currently range has only been tested using chromosomal coordinates. ["range", "exact"].
                  7. File_Size = The size of the db file. If it is large another sub routine is used to collect the information to ensure speed and proper memory handling. ["small", "large"]
               -s/--sampleIDs The sample ID(s),comma sep
               -o/--outFile The output file (defaults to intersectCollect.txt)
               -oinfo/--outInfos The headers and order for each column in output file. The information can also be recorded in the db master file as "outinfo:header1=0_2,,headerN=N_N". NOTE: db file and col nr are 0-based. (if col 2,3 in 1st db file & col 0,1 in 0th db file is desired then ocol should be 1_2,1_3,0_0,0_1. You are free to include, exclude or rearrange the order as you like. Precedence 1. command line 2. Recorded in db master file 3. Order of appearance in db master file.
               -m/--merge Merge all entries found within db files. Unique entries (not found in first db file) will be included. Do not support range matching. (Defaults to "0")
               -sl/--select Select all entries in first infile matching keys in subsequent db files. Do not support range matching. (Defaults to "0")            
               -sofs/--selectOutFiles Selected variants and orphan db files out data directory. Comma sep (Defaults to ".";Supply whole path(s) and in the same order as the '-db' db file)
               -prechr/--prefixChromosomes "chrX" or just "X" (defaults to "X")
               -h/--help Display this help message    
               -v/--version Display version
	   };    
}

my ($dbFile) = (0);
my ($outFile, $outInfo) = ("intersectCollect.txt", 0);
my ($merge, $select) = (0, 0);
my ($prefixChromosomes, $help, $version) = (0);

my (@sampleIDs, @chromosomes, @outColumns, @outHeaders, @outInfos, @selectOutFiles, @metaData);


###User Options
GetOptions('db|dbFile:s'  => \$dbFile,
	   's|sampleIDs:s'  => \@sampleIDs, #Comma separated list
	   'o|outFile:s'  => \$outFile,
	   'oinfo|outInfos:s'  => \@outInfos, #comma separated
	   'm|merge:n'  => \$merge,
	   'sl|select:n'  => \$select,
	   'sofs|selectOutFiles:s'  => \@selectOutFiles, #Comma separated list
	   'prechr|prefixChromosomes:n'  => \$prefixChromosomes,
	   'h|help' => \$help,  #Display help text
	   'v|version' => \$version, #Display version number
    );

if($help) {
    
    print STDOUT "\n".$USAGE, "\n";
    exit;
}

if($version) {
    
    print STDOUT "\nintersectCollect.pl v1.1", "\n\n";
    exit
}

if ($dbFile eq 0) {
    
    print STDOUT "\n".$USAGE, "\n";
    print STDERR "\n", "Need to specify db file by using flag -db. 1 db file per line. Format: DbPath\tSeparator\tColumn_Keys\tChr_Column\tMatching\tColumns_to_Extract\tFile_Size\t","\n\n";
    exit;
}

if (@outInfos) {

    @outInfos = split(/,/,join(',',@outInfos)); #Enables comma separated list
    print STDOUT "Order of output header and columns as supplied by user:", "\n";
    
    for (my $outInfoCounter=0;$outInfoCounter<scalar(@outInfos);$outInfoCounter++) {
	
	print STDOUT $outInfos[$outInfoCounter], "\t";
	
	if ($outInfos[$outInfoCounter] =~/IDN_GT_Call\=/) { #Handle IDN exception
	    
	    push(@outColumns, $'); #' #Collect output columns
	    push(@outHeaders, $`); #Collect output headers
	}
	elsif ($outInfos[$outInfoCounter] =~/\=/) {

	    push(@outColumns, $'); #' #Collect output columns
	    push(@outHeaders, $`); #Collect output headers
	}
    } 
    print STDOUT "\n";
    $outInfo = 1; #To not rewrite order supplied by user with the order in the Db master file
}

if ($prefixChromosomes == 0) { #Ensembl - no prefix and MT

    @chromosomes = ("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT"); #Chromosomes for enhanced speed in collecting information and reducing memory consumption
}
else { #Refseq - prefix and MT

    @chromosomes = ("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrMT");
}
	
my $dbFileCounter=0; #Global count of the nr of files in $dbFile
my (@allVariants, @allVariantsUnique, @allVariantsSorted); #temporary arrays for sorting
my (%dbFile, %RangeVariants, %RangeBins, %dbFilePos); #Db file structure, range/bin intervalls, binary positions in db files
my (%selectFilehandles, %selectVariants);
my (%allVariants, %allVariantsChromosome, %allVariantsChromosomeUnique, %allVariantsChromosomeSorted); #hash for collected information and temporary hashes for sorting
my %unSorted; #For later sorting using ST 

@sampleIDs = split(/,/,join(',',@sampleIDs)); #Enables comma separated sampleID(s)
@selectOutFiles = split(/,/,join(',',@selectOutFiles)); #Enables comma separated selectOutFiles(s)

###
#MAIN
###

if ($dbFile) {

    &ReadDbMaster($dbFile,$outInfo); #Collect information on db files from master file supplied with -db
}

#Read all range Db file first to enable check against first db file keys as it is read.
for (my $dbFileNr=0;$dbFileNr<$dbFileCounter;$dbFileNr++) {

    if ( ($dbFile{$dbFileNr}{'Matching'} eq "range") && ($dbFile{$dbFileNr}{'Size'} ne "tabix") ) {

	if ( ($merge == 0) && ( $select == 0) ) { #Only include elements found in first db file. Not supported by merge option or select option

	    &ReadDbRange($dbFile{$dbFileNr}{'File'},$dbFileNr);
	}
    }
}
	
for (my $dbFileNr=0;$dbFileNr<$dbFileCounter;$dbFileNr++) {

    if ($dbFileNr ==0) {#Add first keys and columns to extract and determines valid keys for subsequent Db files
	
	if ( ($merge == 0) && ( $select == 0) ) { #Only include elements found in first db file 
	
	    &ReadInfile($dbFile{$dbFileNr}{'File'}, $dbFileNr);
	    %RangeVariants = (); #All done with range queries.
	}
	if ( $select == 1 ) {
	    
	    if ($dbFile{0}{'Chr_column'} eq "Na") { #No chromosome queries
		
		&ReadDbFilesNoChrSelect();
		&ReadInfileSelect($dbFile{$dbFileNr}{'File'}, $dbFileNr);
	    }
	}
    }
    elsif ($dbFile{$dbFileNr}{'Size'} eq "large") {#Reads large Db file(s). Large (long) Db file(s) will slow down the performance due to the fact that the file using the ReadDbFiles sub routine would be scanned linearly for every key (i.e. if chr coordinates are used). To circumvent this the large db file(s) are completely read and all entries matching the infile are saved using the keys supplied by the user. This is done for all large db file(s), and the extracted information is then keept in memory. This ensures that only entries matching the keys in the infile are keept, keeping the memory use low compared to reading the large db and keeping the whole large db file(s) in memory. As each chr in is processed in ReadDbFiles the first key is erased reducing the search space (somewhat) and memory consumption (more).
	
	if ( $merge == 0 ) { #Only include elements found in first db file 
	    
	    &ReadDbLarge($dbFile{$dbFileNr}{'File'},$dbFileNr);
	}
    }
}

#if chromsome Number are used
if ($dbFile{0}{'Chr_column'} ne "Na") { #For chromosome queries
 
    if ($merge == 0) { #Only include elements found in first db file
	
	for (my $chromosomeNumber=0;$chromosomeNumber<scalar(@chromosomes);$chromosomeNumber++) {

	    &ReadDbFilesTabix($chromosomes[$chromosomeNumber]);#Check that tabix is specified in sub routine, otherwise leave untouched

	    if ($chromosomes[$chromosomeNumber+1]) {
		
		&ReadDbFiles($chromosomes[$chromosomeNumber],$chromosomes[$chromosomeNumber+1]); #Scans each file for chromosomes entries only processing those i.e. will scan the db file once for each chr in @chromosomes.   
	    }
	    else {
		
		&ReadDbFiles($chromosomes[$chromosomeNumber]); #Last chromosome
	    }
	    &SortAllVariantsST($chromosomes[$chromosomeNumber]); #Sort all variants per chromsome
	    &WriteChrVariants($outFile, $chromosomes[$chromosomeNumber]); #Write all variants to file per chrosome
	    
	    #Reset for next chromosome
	    $allVariants{$chromosomes[$chromosomeNumber]} = (); %allVariantsChromosome = (); %allVariantsChromosomeUnique = (); %allVariantsChromosomeSorted = ();
	    @allVariants = (); @allVariantsUnique = (); @allVariantsSorted = ();
	    
	}
    }
    else  { #Only for db files that are to be merged

	for (my $chromosomeNumber=0;$chromosomeNumber<scalar(@chromosomes);$chromosomeNumber++) {
	    
	    if ($chromosomes[$chromosomeNumber+1]) {

		&ReadInfileMerge($dbFile{0}{'File'}, 0,$chromosomes[$chromosomeNumber],$chromosomes[$chromosomeNumber+1]);
		&ReadDbFilesMerge($chromosomes[$chromosomeNumber], $chromosomes[$chromosomeNumber+1]); #Scans each file for chromsome entries only processing those i.e. will scan the db file once for each chromsome in @chromosomes.   
	    }
	    else {

		&ReadInfileMerge($dbFile{0}{'File'}, 0,$chromosomes[$chromosomeNumber]);
		&ReadDbFilesMerge($chromosomes[$chromosomeNumber]); #Scans each file for chromosome entries only processing those i.e. will scan the db file once for each chromosome in @chromosomes.
	    }
	    &SortAllVariantsMergeST($chromosomes[$chromosomeNumber]); #Sort all variants per chrosome
	    &WriteAllVariantsMerge($outFile, $chromosomes[$chromosomeNumber]); #Write all variants to file per chromosome
	    
	    #Reset for next chromosome
	    %allVariants = (); %allVariantsChromosome = (); %allVariantsChromosomeUnique = (); %allVariantsChromosomeSorted = ();
	    @allVariants = (); @allVariantsUnique = (); @allVariantsSorted = ();
	}
    }
}
else { #Other type of keys

    if ( $select == 0) {

	#&ReadDbFilesNoChr();
	#&WriteAll($outFile);
    }
    else {
	#Do nothing because if select mode is on the db files have already beeen read
    }
}

###
#Sub Routines
###

sub ReadDbFilesTabix {
#Reads all db files collected from Db master file, except first db file and db files with features "large", "range". These db files are handled by different subroutines.   

    my $chromosomeNumber = $_[0];
    
    for (my $dbFileNr=1;$dbFileNr<$dbFileCounter;$dbFileNr++) { #All db files (in order of appearance in $dbFile) except first db which has already been handled
	
	if ($dbFile{$dbFileNr}{'Size'} eq "tabix") { #Only for files with tabix index, other files are handled downstream
	    
	    my $tabix = Tabix->new('-data' => $dbFile{$dbFileNr}{'File'});
	    my @tabixContigs = $tabix->getnames; #Locate all contigs in file

	    if ( grep( /^$chromosomeNumber$/, @tabixContigs ) ) { #Only collect from contigs present in file

		for my $secondKey (keys % {$allVariants{$chromosomeNumber} }) {
		    
		    for my $thirdKey (keys % {$allVariants{$chromosomeNumber}{$secondKey} }) {
			
			my $iteration = $tabix->query( $chromosomeNumber, ($secondKey-1), $thirdKey+1); #Create overlapping region, creates a bit of search overhead, but should not affect accuracy since return is check for association before collecting
			    
			while (my $variantLine = $tabix->read($iteration)){ #Iterate over all postions within region
			    
			    if (defined($variantLine)) {
				
				my @tabixReturnArray = split('\t', $variantLine);
				
				if (scalar( @{$dbFile{$dbFileNr}{'Column_Keys'}}) == 3) {
				    
				    if ( $allVariants{$chromosomeNumber}{$secondKey}{$tabixReturnArray[ $dbFile{$dbFileNr}{'Column_Keys'}[2] ]}) {
					
					for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$dbFileNr}{'Column_To_Extract'}});$extractColumnsCounter++) {
					    
					    my $columnIdRef = \($dbFileNr."_".$dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter]);
					    $allVariants{$chromosomeNumber }{$secondKey}{$tabixReturnArray[ $dbFile{$dbFileNr}{'Column_Keys'}[2] ]}{$$columnIdRef} = $tabixReturnArray[ $dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter] ];
					}
				    }
				}
				if (scalar( @{$dbFile{$dbFileNr}{'Column_Keys'}}) == 4) {
				    
				    if ( $allVariants{$chromosomeNumber}{$secondKey}{$thirdKey}{$tabixReturnArray[ $dbFile{$dbFileNr}{'Column_Keys'}[3] ]}) {
					
					for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$dbFileNr}{'Column_To_Extract'}});$extractColumnsCounter++) {
					    
					    my $columnIdRef = \($dbFileNr."_".$dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter]);
					    $allVariants{$chromosomeNumber }{$secondKey}{$thirdKey }{$tabixReturnArray[ $dbFile{$dbFileNr}{'Column_Keys'}[3] ]}{$$columnIdRef} = $tabixReturnArray[ $dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter] ];
					}
				    }
				}
			    }
			}
		    } 
		}
	    }
	    print STDOUT "Finished Reading chromosome".$chromosomeNumber." in Infile: ".$dbFile{$dbFileNr}{'File'},"\n";   
	}
    }
}

sub ReadDbMaster {
#Reads DbMaster file
#DbPath\tSeparator\Column_Keys\tColumns_to_Extract\t. First 4 columns are mandatory.

    my $fileName = $_[0]; #File name for databse master file
    my $UserSuppliedOutInfoSwitch = $_[1]; #"0" || "1", depending on the user supplied an output headerN=N_N (1) or not (0)

    my (@dbColumnKeys, @dbColumnKeysNrs, @dbExtractColumns);
    
    open(DBM, "<".$fileName) or die "Can't open ".$fileName.":".$!, "\n";    
    
    while (<DBM>) {
	
	chomp $_; #Remove newline
	
	if (m/^\s+$/) {		# Avoid blank lines
	    next;
	}
	if ($_=~/^##/) {#MetaData

	    push(@metaData, $_); #Save metadata string
	    next;
	}
	if ($_=~/^outInfo:/i) { #Locate order of out columns if recorded in db master (precedence: 1. Command line, 2. Recorded in db master file 3. Order of appearance in db master file)
	   
	    if ($UserSuppliedOutInfoSwitch == 0) { #Create output headers and order determined by entry in db master file
		
		$outInfo = 1; #Ensure that precedence is kept
		@outInfos = split(",", $'); #'

		my $IDNCounter = 0;

		for (my $outInfoCounter=0;$outInfoCounter<scalar(@outInfos);$outInfoCounter++) {
		    
		    if ($outInfos[$outInfoCounter] =~/\=>/) {
			
			push(@outColumns, $'); #'
			
			if ( ($sampleIDs[$IDNCounter]) && ($outInfos[$outInfoCounter] =~/IDN_GT_Call\=/) ) { #Handle IDN exception
			
			    push(@outHeaders, "IDN:".$sampleIDs[$IDNCounter]);
			    $IDNCounter++;
			}
			else {
			
			    push(@outHeaders, $`);
			}
		    }
		} 
	    }
	    next;
	}		
	if ( $_ =~/^(\S+)/ ) {	

	    my @dbElements = split("\t",$_); #Loads databse elements description

	    #Set mandatory values
	    $dbFile{$dbFileCounter}{'File'}=$dbElements[0]; #Add dbFile Name
	    $dbFile{$dbFileCounter}{'Separator'}=$dbElements[1]; #Add dbFile Seperator
	    @dbColumnKeys = split(/,/, join(',',$dbElements[2]) ); #Enable comma separeted entry for column keys
	    push (@dbColumnKeysNrs, scalar(@dbColumnKeys)); #Check that correct Nr of keys are supplied for all files (dictated by first file)
	    $dbFile{$dbFileCounter}{'Column_Keys'}= [@dbColumnKeys]; #Add dbFile column keys
	    $dbFile{$dbFileCounter}{'Chr_column'}=$dbElements[3]; #Add dbFile Chromosome column (For skipping lines in files, should be replaced by general indexing of db files in the future
	    $dbFile{$dbFileCounter}{'Matching'}=$dbElements[4]; # (exact or range). Range only valid for chromosome coordinates or similiar. 

	    if ($dbFile{$dbFileCounter}{'Matching'} eq "range") { #Check range queries 

		if ($dbFileCounter == 0) { #Cannot handle range queries for the determining elements
		    
		    print STDERR "\n", "First Db file should be used with exact mathing \n";
		    print STDOUT "\n".$USAGE, "\n";
		    exit;
		}
		if ( (scalar(@dbColumnKeys)  <=2 ) ) { #Must have at least three keys presently to perform range queries

		    print STDERR "\n", "To few keys to use for range comparison, currently only three keys are supported for the range file, but for are allowed within same run for exact Db files \n";
		    print STDOUT "\n".$USAGE, "\n";
		    exit;
		}
	    }
	    if ( scalar(@dbColumnKeysNrs) > 1 ) { #Check that same number of keys are consistently added. 

		if ( ($dbColumnKeysNrs[0] ne $dbColumnKeysNrs[1]) && ($dbFile{$dbFileCounter}{'Matching'} ne "range") ) {

		    print STDERR "Not the same number of keys as in previous db files in db file: ".$dbFile{$dbFileCounter}{'File'}, "\n";
		    print STDOUT "\n".$USAGE, "\n";
		    exit;
		}
		pop(@dbColumnKeysNrs); #Remove database keys entry 
	    }

	    @dbExtractColumns = split(/,/, join(',',$dbElements[5]) ); #Enable comma separeted entry for columns to extract

	    if ($outInfo eq 0) { #Create output order determined by appearance in db master file

		for (my $outColumnCounter=0;$outColumnCounter<scalar(@dbExtractColumns);$outColumnCounter++) {

		    push (@outColumns, $dbFileCounter."_".$dbExtractColumns[$outColumnCounter]);
		}
	    }
	    $dbFile{$dbFileCounter}{'Column_To_Extract'}= [@dbExtractColumns]; #Add dbFile columns to extract
	    $dbFile{$dbFileCounter}{'Size'}=$dbElements[6]; #Determine the way to parse the db file
###
#Validation	    
###
	    print STDOUT $dbFile{$dbFileCounter}{'File'}."\t".$dbFile{$dbFileCounter}{'Separator'},"\t";

	    for (my $columnKeysCounter=0;$columnKeysCounter<scalar( @ {$dbFile{$dbFileCounter}{'Column_Keys'} });$columnKeysCounter++) {
		
		print STDOUT $dbFile{$dbFileCounter}{'Column_Keys'}[$columnKeysCounter], ",";
	    }
	    print STDOUT "\t";
	    print STDOUT $dbFile{$dbFileCounter}{'Chr_column'},"\t";

	    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @ {$dbFile{$dbFileCounter}{'Column_To_Extract'} });$extractColumnsCounter++) {

		print STDOUT $dbFile{$dbFileCounter}{'Column_To_Extract'}[$extractColumnsCounter], ",";
	    }
	    print STDOUT "\t";
	    print STDOUT $dbFile{$dbFileCounter}{'Matching'}, "\t", $dbFile{$dbFileCounter}{'Size'}, "\n";
	    $dbFileCounter++;
	}
    }
    if ($UserSuppliedOutInfoSwitch eq 0) { #No order of output headers and columns supplied by user

	if ($outInfo == 1) {

	    print STDOUT "Order of output headers and columns determined by db master file: ","\n";

	    for (my $outInfoCounter=0;$outInfoCounter<scalar(@outInfos);$outInfoCounter++) {
	
		print STDOUT $outInfos[$outInfoCounter], "\t";
	    }
	    print STDOUT "\n";
	}
 	else {
	    print STDOUT "No users supplied order of output header and columns. Will order the columns according to appearance in db master file\n";
	}
    }
    if ( scalar(@selectOutFiles) eq 0 ) { #Add relative path if none were specified. Need to know the number of db files to link the output paths to each db file
	
	for (my $dbFileNr=0;$dbFileNr<$dbFileCounter;$dbFileNr++) {
	
	    my $dbfileBaseName = basename($dbFile{$dbFileNr}{'File'});
	    $selectOutFiles[$dbFileNr] = $dbfileBaseName."selectVariants";
	}
    }
    close(DBM);
    print STDOUT "Finished Reading Db file: ".$fileName,"\n";
    print STDOUT "Found ".$dbFileCounter." Db files","\n";
    return;
}

sub ReadDbFiles {
#Reads all db files collected from Db master file, except first db file and db files with features "large", "range". These db files are handled by different subroutines.

    my $chromosome = $_[0];
    my $nextChromosome = $_[1];

    for (my $dbFileNr=1;$dbFileNr<$dbFileCounter;$dbFileNr++) { #All db files (in order of appearance in $dbFile) except first db which has already been handled     
	
	if ( ($dbFile{$dbFileNr}{'Size'} eq "tabix") || ($dbFile{$dbFileNr}{'Size'} eq "large") || ($dbFile{$dbFileNr}{'Matching'} eq "range") ) { #Already handled

	    next;
	}
	else { #Read files per chrosome (i.e. multiple times)
	    
	    open(DBF, "<".$dbFile{$dbFileNr}{'File'}) or die "Can't open ".$dbFile{$dbFileNr}{'File'}.": ".$!, "\n";    
	    
	    if ( defined($dbFilePos{$dbFileNr}) ) { #if file has been searched previously
	
		seek(DBF, $dbFilePos{$dbFileNr},0) or die "Couldn't seek to ".$dbFilePos{$dbFileNr}." in ".$dbFile{$dbFileNr}{'File'}.": ".$!,"\n"; #Seek to binary position in file where we left off
	    }
	    while (<DBF>) {
		
		chomp $_;
		###VALIDATION
		#if ($.==1) {
		  #  print STDERR "Started at line $.\n";
		   # my $pos = tell(DBF);
		   # if ( defined ($dbFilePos{$dbFileNr}) ) {print STDERR "Started at pos ", $dbFilePos{$dbFileNr}, "\n";}
		   # else {print STDERR "Started at pos ",  $pos, "\n";}
		#}
		if (m/^\s+$/) {		# Avoid blank lines
		    next;
		}
		if (m/^#/) {		# Avoid #
		    next;
		}		
		if ( $_ =~/^(\S+)/) {
		    
		    my @lineElements = split($dbFile{$dbFileNr}{'Separator'},$_); #Splits columns on separator and loads line
		    
		    if ($nextChromosome && $lineElements[$dbFile{$dbFileNr}{'Chr_column'}] eq $nextChromosome) { #If next chrosome is found return (Since all numerically infiles are sorted this is ok)
			
			###VALIDATION
			#print STDOUT "Finished Reading chromosomes$chromosome in Infile $dbFile{$dbFileNr}{'File'}","\n";
			#$dbFilePos{$dbFileNr} = tell(DBF); # Save  binary position in file to enable seek when revisiting e.g. next chromsome
			#print "Ended at pos $dbFilePos{$dbFileNr} in $dbFile{$dbFileNr}{'File'}\n";
			close(DBF);
			last;
		    }
		    if ( $lineElements[$dbFile{$dbFileNr}{'Chr_column'}] eq $chromosome) { #If chromosome number and chromosome in line match - go ahead and process
			
                        #Depending on the number of column keys supplied by user in db master file . 
			if (scalar( @{$dbFile{$dbFileNr}{'Column_Keys'}}) == 1) {
			    
			    if ( $allVariants{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[0] ] }) { #If first key match to already processed first db file 

				for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$dbFileNr}{'Column_To_Extract'}});$extractColumnsCounter++) {
				    
				    my $columnIdRef = \($dbFileNr."_".$dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter]);
				    $allVariants{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[0] ] }{$$columnIdRef} = $lineElements[ $dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter] ];
				}
			    }
			}
			if (scalar( @{$dbFile{$dbFileNr}{'Column_Keys'}}) == 2) {

			    if ( $allVariants{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[1] ] }) {

				for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$dbFileNr}{'Column_To_Extract'}});$extractColumnsCounter++) {
				    
				    my $columnIdRef = \($dbFileNr."_".$dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter]);
				    $allVariants{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[1] ] }{$$columnIdRef} = $lineElements[ $dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter] ];
				}
			    }
			}
			if (scalar( @{$dbFile{$dbFileNr}{'Column_Keys'}}) == 3) {

			    if ( $allVariants{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[1] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[2] ]}) {

				for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$dbFileNr}{'Column_To_Extract'}});$extractColumnsCounter++) {
				    
				    my $columnIdRef = \($dbFileNr."_".$dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter]);
				    $allVariants{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[0] ]}{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[1] ]}{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[2] ]}{$$columnIdRef} = $lineElements[ $dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter] ];
				}
			    }
			}
			if (scalar( @{$dbFile{$dbFileNr}{'Column_Keys'}}) == 4) {
			    
			    if ( $allVariants{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[1] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[2] ]}{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[3] ]}) {
				for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$dbFileNr}{'Column_To_Extract'}});$extractColumnsCounter++) {
				    
				    my $columnIdRef = \($dbFileNr."_".$dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter]);
				    $allVariants{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[1] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[2] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[3] ]}{$$columnIdRef} = $lineElements[ $dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter] ];
				}
			    }
			}
		    }
		} 	
	    }
	}
	close(DBF);
	print STDOUT "Finished Reading chromosome".$chromosome." in Infile: ".$dbFile{$dbFileNr}{'File'},"\n";
    }
    return;
}

sub ReadInfileSelect {
#Reads the first db file, which is the file that all subsequent elements will matched to i.e. only information for elements present in the first file will be added)

    my $dbFileName = $_[0];
    my $DbFileNumber = $_[1];

    for (my $dbFileNr=0;$dbFileNr<$dbFileCounter;$dbFileNr++) { #All db files (in order of appearance in $db) except first file

	$selectFilehandles{ $dbFile{$dbFileNr}{'File'} }=IO::Handle->new(); #Create anonymous filehandle
	
	open($selectFilehandles{ $dbFile{$dbFileNr}{'File'} }, ">".$selectOutFiles[$dbFileNr]) or die "Can't open ".$selectOutFiles[$dbFileNr].":".$!, "\n"; #open file(s) for db output
	print STDOUT "Select Mode: Writing Selected Variants to: ".$selectOutFiles[$dbFileNr], "\n";  
    }

    my %selectedSwithc; #For printing not selected records to orphan file
    my %writeTracker; #For tracking the number of prints to each db file and wether to print to orphan as well

    open(RIFS, "<".$dbFileName) or die "Can't open ".$dbFileName.":".$!, "\n"; #Open first infile
    
    while (<RIFS>) {
	
	chomp $_; #Remove newline
	
	if (m/^\s+$/) {		# Avoid blank lines
	    next;
	}
	if ($_=~/^##/) {#MetaData

	    for (my $dbFileNr=0;$dbFileNr<$dbFileCounter;$dbFileNr++) { #All db files (in order of appearance in $db) except first file

		print { $selectFilehandles{ $dbFile{$dbFileNr}{'File'} } } $_,"\n";
	    }
	    next;
        }
	if (m/^#/) {		#Header info

	    print { $selectFilehandles{ $dbFile{0}{'File'} } } $_, "\n"; #Header in original file, print to unselected variants

	    for (my $dbFileNr=0;$dbFileNr<$dbFileCounter;$dbFileNr++) { #All db files (in order of appearance in $db) except first file

		if ( ($outInfo == 1) && ($dbFileNr>0) ) { #Print header if supplied, but not for the unselected variants then keep original header (if any)

		    print { $selectFilehandles{ $dbFile{$dbFileNr}{'File'} } } "#";

		    my $IDNCounter = 0;

		    for (my $outHeaderCounter=0;$outHeaderCounter<scalar(@outHeaders);$outHeaderCounter++) {
			
			if ( ($sampleIDs[$IDNCounter]) && ($outHeaders[$outHeaderCounter] =~/IDN_GT_Call\=/) ) { #Handle IDN exception
			    
			    print { $selectFilehandles{ $dbFile{$dbFileNr}{'File'} } } "IDN:".$sampleIDs[$IDNCounter],"\t";
			    $IDNCounter++;
			}
			else {
			   
			    print { $selectFilehandles{ $dbFile{$dbFileNr}{'File'} } } $outHeaders[$outHeaderCounter],"\t";
			}
			
		    }
		    print { $selectFilehandles{ $dbFile{$dbFileNr}{'File'} } } "\n";
		}
	    }
	    next;
	}		
	if ( $_ =~/^(\S+)/ ) {

	    my @lineElements = split($dbFile{$DbFileNumber}{'Separator'},$_); #Loads line

	    if (scalar( @{$dbFile{0}{'Column_Keys'}}) == 1) {
		
		my @parsedColumns = split(';',$lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ]); #For entries with X;Y
		
		for (my $dbFileNr=1;$dbFileNr<$dbFileCounter;$dbFileNr++) { #All db files (in order of appearance in $db) except first file

		    $selectedSwithc{$dbFileNr}=0;
		    $writeTracker{$dbFileNr}=0;
		    my $dbWroteSwitch=0; #make sure that there are no duplictate entries

		    for (my $parsedColumnsCounter=0;$parsedColumnsCounter<scalar(@parsedColumns);$parsedColumnsCounter++) { #Loop through all
		
			if ( $selectVariants{$dbFileNr}{$parsedColumns[ $parsedColumnsCounter ]} ) { #If key exists in db file

   
			    $selectedSwithc{$dbFileNr}++; #Increment switch for correct Db file

			    my $filehandle = $selectFilehandles{ $dbFile{$dbFileNr}{'File'} }; #Collect correct anonymous filehandle
			    		
			    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$DbFileNumber}{'Column_To_Extract'}});$extractColumnsCounter++) {
				
				my $columnIdRef = \($DbFileNumber."_".$dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter]);
				$allVariants{ $parsedColumns[ $parsedColumnsCounter ] }{$$columnIdRef}= $lineElements[ $dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter] ]; #Collect all columns to enable print later
			    }
			    if ( ($selectedSwithc{$dbFileNr} == 1) &&  ($dbWroteSwitch == 0) ) { #Print record only once to avoid duplicates

				for (my $outColumnCounter=0;$outColumnCounter<scalar(@outColumns);$outColumnCounter++ ) { #Print all outInfo from both files

				    print { $filehandle } $allVariants{ $parsedColumns[ $parsedColumnsCounter ] }{$outColumns[$outColumnCounter]}, "\t";
				}
				print { $filehandle } "\n";
				$dbWroteSwitch++; #Do not print an entry more than once to the outFile i.e. no variant duplications
				#last; #Do not print an entry more than once to the outFile i.e. no variant duplications
			    }
			    if ( ($selectedSwithc{$dbFileNr} > 0) && ($selectedSwithc{$dbFileNr} == scalar(@parsedColumns)) ) { #Only Hit in Db file and no other genes outside the Db.
				
				$writeTracker{$dbFileNr}++;

			    }
			}
		    }    
		}
		for (my $dbFileNr=1;$dbFileNr<$dbFileCounter;$dbFileNr++) { #All db files (in order of appearance in $db) except first file

		    if ( $writeTracker{$dbFileNr} == 0 ) { #Hit both Db and with another geneID outside of the Db
		
			print { $selectFilehandles{ $dbFile{0}{'File'} } } $_, "\n"; #write original record to orphan file
			last; #Do not print an entry more than once to the outFile i.e. no variant duplications
		    }
		    else { #Reset switch and tracker

			$selectedSwithc{$dbFileNr}=0;
			$writeTracker{$dbFileNr}=0;
		    }
		}
	    }
	    if (scalar( @{$dbFile{0}{'Column_Keys'}}) == 2) {

		for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$DbFileNumber}{'Column_To_Extract'}});$extractColumnsCounter++) {
		    
		    my $columnIdRef = \($DbFileNumber."_".$dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter]);
		    $allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] }{$$columnIdRef}= $lineElements[ $dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter] ];
		}
		push (@{$unSorted{$lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ]}}, $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ]);
	    }
	    if (scalar( @{$dbFile{0}{'Column_Keys'}}) == 3) {

		for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$DbFileNumber}{'Column_To_Extract'}});$extractColumnsCounter++) {
		    
		    my $columnIdRef = \($DbFileNumber."_".$dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter]);
		    $allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ] }{$$columnIdRef}= $lineElements[ $dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter] ];	    
		}
		push (@{$unSorted{$lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ]}}, $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ]);

		for my $rangeDbFiles (keys %RangeVariants) { #All Range Db Files
		    
		    for my $firstKey (keys % { $RangeVariants{$rangeDbFiles} })  { #For all first key e.g. chromosomeNr
			
			if ($firstKey eq $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ]) { #Find correct key e.g. chromosomeNr (enhance speed)
			   
			    for my $binNr (keys %{ $RangeBins{$rangeDbFiles}{$firstKey} }) { #For all bins within range Db file and first key (enhance speed)
					
				if ($lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] >= $RangeBins{$rangeDbFiles}{$firstKey}{$binNr}{'lowerboundary'} && $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] <= $RangeBins{$rangeDbFiles}{$firstKey}{$binNr}{'upperboundary'} ) { #Found in Range Bin Db file, correct firstKey, within bin and second key overlapps bin range.
				    
				    for my $startPosition (keys % { $RangeVariants{$rangeDbFiles}{$firstKey}{$binNr} })  { #For all Db file range start positions within bin
				 	
					if ( $startPosition <= $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ]) { #second key is downstream of start position
					    
					    for my $stopPosition (keys % { $RangeVariants{$rangeDbFiles}{$firstKey}{$binNr}{$startPosition} })  { #For all range stop positions 
						
						if ( $stopPosition >= $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] ) { #second key is upstream of stop position i.e overlapping
						    
						    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$rangeDbFiles}{'Column_To_Extract'}});$extractColumnsCounter++) {

							my $columnIdRef = \($rangeDbFiles."_".$dbFile{$rangeDbFiles}{'Column_To_Extract'}[$extractColumnsCounter]); #RangeDb file is not random key but ordered Db file from initial for loop.		
							$allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ] }{$$columnIdRef}.= $RangeVariants{$rangeDbFiles}{$firstKey}{$binNr}{$startPosition}{$stopPosition}{$$columnIdRef}.";";
						    }
						}
					    }
					}
					#Catches events where second key is upstream of range start, but third key is within range i.e. we do not require complete feature overlap
					elsif ( $startPosition <= $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ]) { #third key is downstream of start position
					    
					    for my $stopPosition (keys % { $RangeVariants{$rangeDbFiles}{$firstKey}{$binNr}{$startPosition} })  { #For all range stop positions 
						
						if ( $stopPosition >= $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ] ) { #third key is upstream of stop position i.e overlapping
						    
						    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$rangeDbFiles}{'Column_To_Extract'}});$extractColumnsCounter++) {

							my $columnIdRef = \($rangeDbFiles."_".$dbFile{$rangeDbFiles}{'Column_To_Extract'}[$extractColumnsCounter]); #RangeDb file is not random key but ordered Db file from initial for loop.		
							$allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ] }{$$columnIdRef}.= $RangeVariants{$rangeDbFiles}{$firstKey}{$binNr}{$startPosition}{$stopPosition}{$$columnIdRef}.";";
						    }
						}
					    }
					}
				    }
				    last; #No overlapping bins - move to next line
				}
			    }
			}
		    }
		}
	    } 
	    if (scalar( @{$dbFile{0}{'Column_Keys'}}) == 4) {
		
		for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$DbFileNumber}{'Column_To_Extract'}});$extractColumnsCounter++) {
		    
		    my $columnIdRef = \($DbFileNumber."_".$dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter]);
		    $allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[3] ]}{$$columnIdRef}= $lineElements[ $dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter] ];
		    
		}
		push (@{$unSorted{$lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ]}}, $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ]);

		for my $rangeDbFiles (keys %RangeVariants) { #All Range Db Files
		    
		    for my $firstKey (keys % { $RangeVariants{$rangeDbFiles} })  { #For all first key e.g. chromosomeNr
			
			if ($firstKey eq $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ]) { #Find correct key e.g. chromosomeNr (enhance speed)
			    
			    for my $binNr (keys %{ $RangeBins{$rangeDbFiles}{$firstKey} }) { #For all bins within range Db file and first key (enhance speed)
				
				if ($lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] >= $RangeBins{$rangeDbFiles}{$firstKey}{$binNr}{'lowerboundary'} && $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] <= $RangeBins{$rangeDbFiles}{$firstKey}{$binNr}{'upperboundary'} ) { #Found in Range Bin Db file, correct firstKey, within bin and second key overlapps bin range.
				    
				    for my $startPosition (keys % { $RangeVariants{$rangeDbFiles}{$firstKey}{$binNr} })  { #For all Db file range start positions within bin

					if ( $startPosition <= $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ]) { #second key is downstream of start position
					    
					    for my $stopPosition (keys % { $RangeVariants{$rangeDbFiles}{$firstKey}{$binNr}{$startPosition} })  { #For all range stop positions 
						
						if ( $stopPosition >= $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] ) { #second key is upstream of stop position i.e overlapping
						    
						    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$rangeDbFiles}{'Column_To_Extract'}});$extractColumnsCounter++) {

							my $columnIdRef = \($rangeDbFiles."_".$dbFile{$rangeDbFiles}{'Column_To_Extract'}[$extractColumnsCounter]); #RangeDb file is not random key but ordered Db file from initial for loop.		
							$allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[3] ]}{$$columnIdRef}.= $RangeVariants{$rangeDbFiles}{$firstKey}{$binNr}{$startPosition}{$stopPosition}{$$columnIdRef}.";";
						    }
						}
					    }
					}
					#Catches events where second key is upstream of range start, but third key is within range i.e. we do not require complete feature overlap
					elsif ( $startPosition <= $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ]) { #third key is downstream of start position
					    
					    for my $stopPosition (keys % { $RangeVariants{$rangeDbFiles}{$firstKey}{$binNr}{$startPosition} })  { #For all range stop positions 
						
						if ( $stopPosition >= $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ] ) { #third key is upstream of stop position i.e overlapping
						    
						    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$rangeDbFiles}{'Column_To_Extract'}});$extractColumnsCounter++) {

							my $columnIdRef = \($rangeDbFiles."_".$dbFile{$rangeDbFiles}{'Column_To_Extract'}[$extractColumnsCounter]); #RangeDb file is not random key but ordered Db file from initial for loop.		
							$allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[3] ]}{$$columnIdRef}.= $RangeVariants{$rangeDbFiles}{$firstKey}{$binNr}{$startPosition}{$stopPosition}{$$columnIdRef}.";";
						    }
						}
					    }
					}
				    }
				    last; #No overlapping bins since border spanning events are added to all Bin that the event is spanned - move to next line
				}
			    }
			}
		    }
		}
	    } 	
	}
    }
    close(RIFS);
    print STDOUT "Select Mode: Finished Reading key Db file: ".$dbFileName,"\n";

    for (my $dbFileNr=1;$dbFileNr<$dbFileCounter;$dbFileNr++) { #All db files (in order of appearance in $db) except first db which has already been handled
	
	close($selectFilehandles{ $dbFile{$dbFileNr}{'File'} }); 
    }
    return;
}

sub ReadDbFilesNoChrSelect {
#Reads all db files collected from Db master file, except first db file and db files with features "large", "range". These db files are handled by different subroutines.
    
    for (my $dbFileNr=1;$dbFileNr<$dbFileCounter;$dbFileNr++) { #All db files (in order of appearance in $db) except first db which has already been handled 
	
	
	if ( ($dbFile{$dbFileNr}{'Size'} eq "large") || ($dbFile{$dbFileNr}{'Matching'} eq "range") ) { #Already handled

	    next;
	}
	else { #Read files

	    open(DBF, "<".$dbFile{$dbFileNr}{'File'}) or die "Can't open ".$dbFile{$dbFileNr}{'File'}.":".$!,"\n";    
	    
	    while (<DBF>) {
		
		chomp $_; #remove newline
		
		if (m/^\s+$/) {		# Avoid blank lines
		    next;
		}
		if (m/^#/) {		# Avoid #
		    next;
		}		
		if ( $_ =~/^(\S+)/) {
		    
		    my @lineElements = split($dbFile{$dbFileNr}{'Separator'},$_); #Splits columns on separator and loads line

		    if (scalar( @{$dbFile{$dbFileNr}{'Column_Keys'}}) == 1) {

			my @parsed_column = split(';',$lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[0] ]); #For entries with X;Y
			
			for (my $parsedColumnCounter=0;$parsedColumnCounter<scalar(@parsed_column);$parsedColumnCounter++) { #Loop through all
			    
			    $selectVariants{$dbFileNr}{$parsed_column[$parsedColumnCounter]} = $parsed_column[$parsedColumnCounter]; #Add key entrie(s)
			    
			    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$dbFileNr}{'Column_To_Extract'}});$extractColumnsCounter++) { #Enable collecton of columns from db file
				
				my $columnIdRef = \($dbFileNr."_".$dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter]);
			    
				$allVariants{ $parsed_column[$parsedColumnCounter] }{$$columnIdRef} = $lineElements[ $dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter] ];
			    }
			}
		    }
###
#NOTE Only 1 key supported so far 130204
###		    
		    if (scalar( @{$dbFile{$dbFileNr}{'Column_Keys'}}) == 2) {
		
			if ( $allVariants{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[1] ] }) {
			
			    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$dbFileNr}{'Column_To_Extract'}});$extractColumnsCounter++) {
				
				my $columnIdRef = \($dbFileNr."_".$dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter]);
				$allVariants{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[1] ] }{$$columnIdRef} = $lineElements[ $dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter] ];
			    }
			}
		    }
		    if (scalar( @{$dbFile{$dbFileNr}{'Column_Keys'}}) == 3) {

			if ( $allVariants{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[1] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[2] ]}) {

			    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$dbFileNr}{'Column_To_Extract'}});$extractColumnsCounter++) {
				
				my $columnIdRef = \($dbFileNr."_".$dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter]);
				$allVariants{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[0] ]}{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[1] ]}{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[2] ]}{$$columnIdRef} = $lineElements[ $dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter] ];
				
			    }
			}
		    }
		    if (scalar( @{$dbFile{$dbFileNr}{'Column_Keys'}}) == 4) {

			if ( $allVariants{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[1] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[2] ]}{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[3] ]}) {
			    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$dbFileNr}{'Column_To_Extract'}});$extractColumnsCounter++) {
				
				my $columnIdRef = \($dbFileNr."_".$dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter]);
				$allVariants{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[1] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[2] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[3] ]}{$$columnIdRef} = $lineElements[ $dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter] ];
			    }
			}
		    }
		}
	    } 	
	}
	close(DBF);
	print STDOUT "Finished Reading Db Infile:".$dbFile{$dbFileNr}{'File'},"\n";
    }
    return;
}

sub ReadDbFilesNoChr {
#Reads all db files collected from Db master file, except first db file and db files with features "large", "range". These db files are handled by different subroutines.
    
    for (my $dbFileNr=1;$dbFileNr<$dbFileCounter;$dbFileNr++) { #All db files (in order of appearance in $db) except first db which has already been handled 
	
	
	if ( ($dbFile{$dbFileNr}{'Size'} eq "large") || ($dbFile{$dbFileNr}{'Matching'} eq "range") ) { #Already handled

	    next;
	}
	else { #Read files

	    open(DBF, "<".$dbFile{$dbFileNr}{'File'}) or die "Can't open ".$dbFile{$dbFileNr}{'File'}.":".$!, "\n";    
	    
	    while (<DBF>) {
		
		chomp $_; #Remove newline
		
		if (m/^\s+$/) {		# Avoid blank lines
		    next;
		}
		if (m/^#/) {		# Avoid #
		    next;
		}		
		if ( $_ =~/^(\S+)/) {
		    
		    my @lineElements = split($dbFile{$dbFileNr}{'Separator'},$_); #Splits columns on separator and loads line
		    
		    ##Depending on the number of column keys supplied by user. NOTE: Must always be the same nr of columns containing the same keys 
		    if (scalar( @{$dbFile{$dbFileNr}{'Column_Keys'}}) == 1) {
		
			if ( $allVariants{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[0] ] } ) { #If first key match
			    
			    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$dbFileNr}{'Column_To_Extract'}});$extractColumnsCounter++) {
				
				my $columnIdRef = \($dbFileNr."_".$dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter]);
				
				$allVariants{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[0] ] }{$$columnIdRef}=$lineElements[ $dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter] ];
##Code to collapse some entries and make serial additions of others. CMMS_External specific and not part of original programe. To enable comment previous line and remove comments from subsequent lines.
				#if ($lineElements[ $dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter] ]) {
				#   if ($extractColumnsCounter>=3) {
				#	$allVariants{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[0] ] }{$$columnIdRef}.="$lineElements[ $dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter] ];";
				#   }
				#  else {
				#	$allVariants{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[0] ] }{$$columnIdRef}=$lineElements[ $dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter] ];
				#   }
				#}
			    }
			}
		    }
		    if (scalar( @{$dbFile{$dbFileNr}{'Column_Keys'}}) == 2) {
			
			if ( $allVariants{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[1] ] }) {
			
			    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$dbFileNr}{'Column_To_Extract'}});$extractColumnsCounter++) {
				
				my $columnIdRef = \($dbFileNr."_".$dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter]);
				$allVariants{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[1] ] }{$$columnIdRef} = $lineElements[ $dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter] ];
			    }
			}
		    }
		    if (scalar( @{$dbFile{$dbFileNr}{'Column_Keys'}}) == 3) {

			if ( $allVariants{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[1] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[2] ]}) {

			    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$dbFileNr}{'Column_To_Extract'}});$extractColumnsCounter++) {
				
				my $columnIdRef = \($dbFileNr."_".$dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter]);
				$allVariants{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[0] ]}{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[1] ]}{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[2] ]}{$$columnIdRef} = $lineElements[ $dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter] ];
				
			    }
			}
		    }
		    if (scalar( @{$dbFile{$dbFileNr}{'Column_Keys'}}) == 4) {
			
			if ( $allVariants{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[1] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[2] ]}{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[3] ]}) {
			    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$dbFileNr}{'Column_To_Extract'}});$extractColumnsCounter++) {
				
				my $columnIdRef = \($dbFileNr."_".$dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter]);
				$allVariants{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[1] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[2] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[3] ]}{$$columnIdRef}= $lineElements[ $dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter] ];
			    }
			}
		    }
		}
	    } 	
	}
	close(DBF);
	print STDOUT "Finished Reading Db Infile: ".$dbFile{$dbFileNr}{'File'},"\n";
    }
    return;
}

sub ReadDbLarge {
#Reads large Db file(s). Large (long) Db file(s) will slow down the performance due to the fact that the file is scanned linearly for every key (i.e. if chromosome coordinates are used). To circumvent this the large db file(s) are completely read and all entries matching the infile are saved using the keys supplied by the user. This is done for all large db file(s), and the extracted information is then keept in memory. This ensures that only entries matching the keys in the infile are keept, keeping the memory use low compared to reading the large db and keeping the whole large db file(s) in memory. As each chromosome in is processed in ReadDbFiles the first key is erased reducing the search space and memory consumption.
    
    my $fileName = $_[0];
    my $DbFileNumber = $_[1];

    open(RDBL, "<".$fileName) or die "Can't open ".$fileName.":".$!, "\n";    
    
    while (<RDBL>) {
	
	chomp $_; #Remove newline
	
	if (m/^\s+$/) {		# Avoid blank lines
	    next;
	}
	if (m/^#/) {		# Avoid #
	    next;
	}		
	if ( $_ =~/^(\S+)/ ) {	

	    my @lineElements = split($dbFile{$DbFileNumber}{'Separator'},$_); #Loads line

	    if (scalar( @{$dbFile{$DbFileNumber}{'Column_Keys'}}) == 1) {

		if ( $allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] } ) { #Check full entry
	
		    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$DbFileNumber}{'Column_To_Extract'}});$extractColumnsCounter++) { #Extract info
			
			my $columnIdRef = \($DbFileNumber."_".$dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter]);			
			$allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] }{$$columnIdRef}= $lineElements[ $dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter] ];
		    }		    
		}
	    }
	    if (scalar( @{$dbFile{$DbFileNumber}{'Column_Keys'}}) == 2) {	    
		
		if ( $allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] } ) { #Needed to not run out of memory in hash look-up
		    
		    if ( $allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] } ){ #Check full entry

			for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$DbFileNumber}{'Column_To_Extract'}});$extractColumnsCounter++) { #Extract info
			    
			    my $columnIdRef = \($DbFileNumber."_".$dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter]);			
			    $allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] }{$$columnIdRef}= $lineElements[ $dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter] ];
			}		    
		    }
		}
	    }
	    if (scalar( @{$dbFile{$DbFileNumber}{'Column_Keys'}}) == 3) {

		if ( $allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] } ) { #Needed to not run out of memory in hash look-up
		    
		    if ( $allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] } ){ #Needed to not run out of memory in hash look-up
			
			if ( $allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ] }) { #Check full entry
			  
			    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$DbFileNumber}{'Column_To_Extract'}});$extractColumnsCounter++) { #Extract info
				
				my $columnIdRef = \($DbFileNumber."_".$dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter]);			
				$allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ] }{$$columnIdRef}= $lineElements[ $dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter] ];
			    }		    
			}
		    }		
		    
		}
	    }
	    if (scalar( @{$dbFile{$DbFileNumber}{'Column_Keys'}}) == 4) {

		if ( $allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] } ) { #Needed to not run out of memory in hash look-up
		    
		    if ( $allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] } ){ #Needed to not run out of memory in hash look-up
			
			if ( $allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[3] ]} ) { #Check full entry

			    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$DbFileNumber}{'Column_To_Extract'}});$extractColumnsCounter++) { #Extract info
				
				my $columnIdRef = \($DbFileNumber."_".$dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter]);			
				$allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[3] ]}{$$columnIdRef}= $lineElements[ $dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter] ];
			    }		    
			}
		    }
		}
	    }
	}
    } 	
    close(RDBL);
    print STDOUT "Finished Reading Large Db file: ".$fileName,"\n";
    return;
}

sub ReadInfile {
#Reads the first db file, which is the file that all subsequent elements will matched to i.e. only information for elements present in the first file will be added)

    my $DbFileName = $_[0];
    my $DbFileNumber = $_[1];
    
    open(RIF, "<".$DbFileName) or die "Can't open ".$DbFileName.":".$!, "\n";    
    
    while (<RIF>) {
	
	chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
	    next;
	}
	if (m/^#/) {		# Avoid #
	    next;
	}		
	if ( $_ =~/^(\S+)/ ) {

	    my @lineElements = split($dbFile{$DbFileNumber}{'Separator'},$_); #Loads line

	    if (scalar( @{$dbFile{0}{'Column_Keys'}}) == 1) {
		
		for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$DbFileNumber}{'Column_To_Extract'}});$extractColumnsCounter++) {
			
		    my $columnIdRef = \($DbFileNumber."_".$dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter]);
		    $allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] }{$$columnIdRef}= $lineElements[ $dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter] ];
		}
	    }
	    if (scalar( @{$dbFile{0}{'Column_Keys'}}) == 2) {

		for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$DbFileNumber}{'Column_To_Extract'}});$extractColumnsCounter++) {
		    
		    my $columnIdRef = \($DbFileNumber."_".$dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter]);
		    $allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] }{$$columnIdRef}= $lineElements[ $dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter] ];
		}
		push (@{$unSorted{$lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ]}}, $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ]);
	    }
	    if (scalar( @{$dbFile{0}{'Column_Keys'}}) == 3) {

		for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$DbFileNumber}{'Column_To_Extract'}});$extractColumnsCounter++) {
		    
		    my $columnIdRef = \($DbFileNumber."_".$dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter]);
		    $allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ] }{$$columnIdRef} = $lineElements[ $dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter] ];	    
		}
		push (@{$unSorted{$lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ]}}, $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ]);

		for my $rangeDbFiles (keys %RangeVariants) { #All Range Db Files
		    
		    for my $firstKey (keys % { $RangeVariants{$rangeDbFiles} })  { #For all first key e.g. chromosomeNr
			
			if ($firstKey eq $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ]) { #Find correct key e.g. chromosomeNr (enhance speed)
			   
			    for my $binNr (keys %{ $RangeBins{$rangeDbFiles}{$firstKey} }) { #For all bins within range Db file and first key (enhance speed)
					
				if ($lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] >= $RangeBins{$rangeDbFiles}{$firstKey}{$binNr}{'lowerboundary'} && $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] <= $RangeBins{$rangeDbFiles}{$firstKey}{$binNr}{'upperboundary'} ) { #Found in Range Bin Db file, correct firstKey, within bin and second key overlapps bin range.
				    
				    for my $startPosition (keys % { $RangeVariants{$rangeDbFiles}{$firstKey}{$binNr} })  { #For all Db file range start positions within bin
				 	
					if ( $startPosition <= $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ]) { #second key is downstream of start position
					    
					    for my $stopPosition (keys % { $RangeVariants{$rangeDbFiles}{$firstKey}{$binNr}{$startPosition} })  { #For all range stop positions 
						
						if ( $stopPosition >= $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] ) { #second key is upstream of stop position i.e overlapping
						    
						    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$rangeDbFiles}{'Column_To_Extract'}});$extractColumnsCounter++) {
							
							my $columnIdRef = \($rangeDbFiles."_".$dbFile{$rangeDbFiles}{'Column_To_Extract'}[$extractColumnsCounter]); #RangeDb file is not random key but ordered Db file from initial for loop.		
							$allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ] }{$$columnIdRef}.= $RangeVariants{$rangeDbFiles}{$firstKey}{$binNr}{$startPosition}{$stopPosition}{$$columnIdRef}.";";
						    }
						}
					    }
					}
					#Catches events where second key is upstream of range start, but third key is within range i.e. we do not require complete feature overlap
					elsif ( $startPosition <= $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ]) { #third key is downstream of start position
					    
					    for my $stopPosition (keys % { $RangeVariants{$rangeDbFiles}{$firstKey}{$binNr}{$startPosition} })  { #For all range stop positions 
						
						if ( $stopPosition >= $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ] ) { #third key is upstream of stop position i.e overlapping
						    
						    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$rangeDbFiles}{'Column_To_Extract'}});$extractColumnsCounter++) {

							my $columnIdRef = \($rangeDbFiles."_".$dbFile{$rangeDbFiles}{'Column_To_Extract'}[$extractColumnsCounter]); #RangeDb file is not random key but ordered Db file from initial for loop.
							$allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ] }{$$columnIdRef}.= $RangeVariants{$rangeDbFiles}{$firstKey}{$binNr}{$startPosition}{$stopPosition}{$$columnIdRef}.";";
						    }
						}
					    }
					}
				    }
				    last; #No overlapping bins - move to next line
				}
			    }
			}
		    }
		}
	    } 
	    if (scalar( @{$dbFile{0}{'Column_Keys'}}) == 4) {
		
		for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$DbFileNumber}{'Column_To_Extract'}});$extractColumnsCounter++) {
		    
		    my $columnIdRef = \($DbFileNumber."_".$dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter]);
		    $allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[3] ]}{$$columnIdRef} = $lineElements[ $dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter] ];
		    
		}
		push (@{$unSorted{$lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ]}}, $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ]);

		for my $rangeDbFiles (keys %RangeVariants) { #All Range Db Files
		    
		    for my $firstKey (keys % { $RangeVariants{$rangeDbFiles} })  { #For all first key e.g. chromosomeNr
			
			if ($firstKey eq $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ]) { #Find correct key e.g. chromosomeNr (enhance speed)
			    
			    for my $binNr (keys %{ $RangeBins{$rangeDbFiles}{$firstKey} }) { #For all bins within range Db file and first key (enhance speed)
				
				if ($lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] >= $RangeBins{$rangeDbFiles}{$firstKey}{$binNr}{'lowerboundary'} && $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] <= $RangeBins{$rangeDbFiles}{$firstKey}{$binNr}{'upperboundary'} ) { #Found in Range Bin Db file, correct firstKey, within bin and second key overlapps bin range.
				    
				    for my $startPosition (keys % { $RangeVariants{$rangeDbFiles}{$firstKey}{$binNr} })  { #For all Db file range start positions within bin

					if ( $startPosition <= $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ]) { #second key is downstream of start position
					    
					    for my $stopPosition (keys % { $RangeVariants{$rangeDbFiles}{$firstKey}{$binNr}{$startPosition} })  { #For all range stop positions 
						
						if ( $stopPosition >= $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] ) { #second key is upstream of stop position i.e overlapping
						    
						    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$rangeDbFiles}{'Column_To_Extract'}});$extractColumnsCounter++) {
							
							my $columnIdRef = \($rangeDbFiles."_".$dbFile{$rangeDbFiles}{'Column_To_Extract'}[$extractColumnsCounter]); #RangeDb file is not random key but ordered Db file from initial for loop.
							$allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[3] ]}{$$columnIdRef}.="$RangeVariants{$rangeDbFiles}{$firstKey}{$binNr}{$startPosition}{$stopPosition}{$$columnIdRef};";
						    }
						}
					    }
					}
					#Catches events where second key is upstream of range start, but third key is within range i.e. we do not require complete feature overlap
					elsif ( $startPosition <= $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ]) { #third key is downstream of start position
					    
					    for my $stopPosition (keys % { $RangeVariants{$rangeDbFiles}{$firstKey}{$binNr}{$startPosition} })  { #For all range stop positions 
						
						if ( $stopPosition >= $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ] ) { #third key is upstream of stop position i.e overlapping
						    
						    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$rangeDbFiles}{'Column_To_Extract'}});$extractColumnsCounter++) {

							my $columnIdRef = \($rangeDbFiles."_".$dbFile{$rangeDbFiles}{'Column_To_Extract'}[$extractColumnsCounter]); #RangeDb file is not random key but ordered Db file from initial for loop.		
							$allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[3] ]}{$$columnIdRef}.="$RangeVariants{$rangeDbFiles}{$firstKey}{$binNr}{$startPosition}{$stopPosition}{$$columnIdRef};";
						    }
						}
					    }
					}
				    }
				    last; #No overlapping bins since border spanning events are added to all Bin that the event is spanned - move to next line
				}
			    }
			}
		    }
		}
	    } 	
	}
    }
    close(RIF);
    print STDOUT "Finished Reading key Db file: ".$DbFileName,"\n";
    return;
}



sub ReadDbRange {
#Reads db file for allowing overlapping feature look-up
    
    my $DbFileName = $_[0];
    my $DbFileNumber = $_[1];

#Loop through file to assign bins based on third key highest value
    my %EstimateBinsize;

    open(DBR, "<".$DbFileName) or die "Can't open ".$DbFileName.":".$!, "\n";    

    while (<DBR>) {
	
	chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
	    next;
	}
	if (m/^#/) {		# Avoid #
	    next;
	}		
	if ( $_ =~/^(\S+)/ ) {	

	    my @lineElements = split($dbFile{$DbFileNumber}{'Separator'},$_);	    #Loads line
	 
	    if ( (scalar( @{$dbFile{0}{'Column_Keys'}}) == 3) || (scalar( @{$dbFile{0}{'Column_Keys'}} ) == 4) ) {

		if ( $EstimateBinsize{$DbFileNumber}{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ]}{'highest'} && $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ] > $EstimateBinsize{$DbFileNumber}{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ]}{'highest'} ) { #Within FirstKey and highest value of thirdkey e.g. chromosomeNr and stopPosition

		$EstimateBinsize{$DbFileNumber}{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ]}{'highest'} = $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ];
		}
		else { #Initiate binsize
		 
		    $EstimateBinsize{$DbFileNumber}{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ]}{'highest'} = $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ];
		}
	    }
	}
    }
    
    for my $firstKey (keys %{$EstimateBinsize{$DbFileNumber} }) { #Set size of bins and add to %RangeBins. Used later to add values in Rnage Db file to correct bin and when infile is read to reduce the number of hash look-ups that needs to be performed for range Db files. Assumes random distribution of ranges across bins.
    my $lowerboundary=0;
    my $upperboundary;
    my $binsize = $EstimateBinsize{$DbFileNumber}{$firstKey}{'highest'}/100;
    $upperboundary = $binsize; 

    for (my $binNr=0;$binNr<100;$binNr++) {
	
	$RangeBins{$DbFileNumber}{$firstKey}{$binNr}{'lowerboundary'}=$lowerboundary;
	$RangeBins{$DbFileNumber}{$firstKey}{$binNr}{'upperboundary'}= $upperboundary;
	$lowerboundary = $lowerboundary + $binsize; #Starting from 0
	$upperboundary= $upperboundary + $binsize; #Starting from binsize
    }
}
    %EstimateBinsize=();
    close(DBR);

##Collect values within Range Db file, now that the bin values have been set. A range can end up in several bins if the range spans bin borders. Then as a exackt match is done a bin the range within that bin can actually expand outside of the bin. This is to catch elements that independent of bin/range size.

    open(DBR, "<".$DbFileName) or die "Can't open ".$DbFileName.":".$!, "\n";    

    while (<DBR>) {
	
	chomp $_; #Remove newline
	
	if (m/^\s+$/) {		# Avoid blank lines
	    next;
	}
	if (m/^#/) {		# Avoid #
	    next;
	}		
	if ( $_ =~/^(\S+)/ ) {	

	    my @lineElements = split($dbFile{$DbFileNumber}{'Separator'},$_); #Loads line

	    if ( (scalar( @{$dbFile{0}{'Column_Keys'}}) == 3) || (scalar( @{$dbFile{0}{'Column_Keys'}} ) == 4) ) {

		for my $firstKey (keys %{$RangeBins{$DbFileNumber} }) { #All firstKeys
		    
		    if ( $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] eq $firstKey ) { #Correct firstKey - check bins
		
			for (my $binNr=0;$binNr<100;$binNr++) { #For all bins
			  
			    if ( $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] >= $RangeBins{$DbFileNumber}{$firstKey}{$binNr}{'lowerboundary'} && $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ] <= $RangeBins{$DbFileNumber}{$firstKey}{$binNr}{'upperboundary'}) { #Whole range within bin - go ahead and add
				
				for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$DbFileNumber}{'Column_To_Extract'}});$extractColumnsCounter++) {
				    
				    my $columnIdRef = \($DbFileNumber."_".$dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter]);
				    $RangeVariants{$DbFileNumber}{$lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ]}{$binNr}{$lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ]}{$lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ]}{$$columnIdRef} = $lineElements[ $dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter] ]; #RangeVariants{dbNr}{key1}{BinNr}{range_start}{range_stop}{columnId}
				}
				last;
			    }
			    elsif ( ($lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] >= $RangeBins{$DbFileNumber}{$firstKey}{$binNr}{'lowerboundary'}) && ($lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ] >= $RangeBins{$DbFileNumber}{$firstKey}{$binNr}{'upperboundary'}) && ($lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] <= $RangeBins{$DbFileNumber}{$firstKey}{$binNr}{'upperboundary'}) ) { #Range spans bins

				for (my $spanBinNr=$binNr;$spanBinNr<100;$spanBinNr++) {

				    if ( $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ] >= $RangeBins{$DbFileNumber}{$firstKey}{$spanBinNr}{'upperboundary'} ) { #Not found upper limit yet

					for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$DbFileNumber}{'Column_To_Extract'}});$extractColumnsCounter++) {
					    
					    my $columnIdRef = \($DbFileNumber."_".$dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter]);
					    $RangeVariants{$DbFileNumber}{$lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ]}{$spanBinNr}{$lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ]}{$lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ]}{$$columnIdRef}= $lineElements[ $dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter] ]; #RangeVariants{dbNr}{key1}{BinNr}{range_start}{range_stop}{columnId}
					}
				    }
				    elsif ( $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ] >= $RangeBins{$DbFileNumber}{$firstKey}{$spanBinNr}{'lowerboundary'} ) { #Found upper limit and check to see if upper limit is before start of bin - if so then add to bin and break
					for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$DbFileNumber}{'Column_To_Extract'}});$extractColumnsCounter++) {
					    
					    my $columnIdRef = \($DbFileNumber."_".$dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter]);
					    $RangeVariants{$DbFileNumber}{$lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ]}{$spanBinNr}{$lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ]}{$lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ]}{$$columnIdRef}= $lineElements[ $dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter] ]; #RangeVariants{dbNr}{key1}{BinNr}{range_start}{range_stop}{columnId}
					}
					last;
				    }
				    else {
					
					$binNr=$spanBinNr;
					last;
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
    close(DBR);
    print STDOUT "Finished Reading Range key Db file: ".$DbFileName,"\n";
    return;
}

sub SortAllVariantsST {
#Creates an array of all position which are unique and in sorted ascending order   
    
    my $firstKey = $_[0];
    my %seen = (); 

#Unique entries only
    @{$allVariantsChromosomeUnique{$firstKey} } = grep { ! $seen{$_} ++ } @{$unSorted{$_[0]} };
#Sort using the Schwartzian Transform 
    @{$allVariantsChromosomeSorted{$firstKey} } =  map { $_->[0] }
    sort { $a->[0] <=> $b->[0] }
    map { [$_] }
    @{$allVariantsChromosomeUnique{$firstKey} };

    $unSorted{$_[0]}=(); #Reset for this chromosome
}

sub SortAllVariants {
#Creates an array of all position which are unique and in sorted ascending order  
    
    my $firstKey = $_[0];
    my %seen = (); 

    if (scalar( @{$dbFile{0}{'Column_Keys'}}) > 1) {
	
	for my $secondKey (keys %{ $allVariants{$firstKey} } )  { #For all secondkeys
	    
	    push ( @{$allVariantsChromosome{$firstKey} }, $secondKey );
	}
	@{$allVariantsChromosomeUnique{$firstKey} } = grep { ! $seen{$_} ++ } @{$allVariantsChromosome{$firstKey} }; #Unique entries only 
	@{$allVariantsChromosomeSorted{$firstKey} } = sort { $a <=> $b } @{ $allVariantsChromosomeUnique{$firstKey} }; #Sorts keys to be able to print sorted table later 
	print STDOUT "Sorted all non overlapping entries\n";
    }
}

sub WriteChrVariants {
#Prints tab separated file of all collected db file info in ascending order dictaded by %$allVariantsChromosome
    
    my $filename = $_[0];
    my $chromosomeNumber = $_[1];
    
    if ( ($chromosomeNumber eq 1) || ($chromosomeNumber eq "chr1") ) {

	open (WAV, ">".$filename) or die "Can't write to ".$filename.":".$!, "\n";

	if (@metaData) { #Print metaData if supplied

	    for (my $metaDataCounter=0;$metaDataCounter<scalar(@metaData);$metaDataCounter++) {

		print WAV $metaData[$metaDataCounter],"\n";
	    }
	}
	if ($outInfo == 1) { #Print header if supplied

	    print WAV "#";

	    for (my $outHeaderCounter=0;$outHeaderCounter<scalar(@outHeaders);$outHeaderCounter++) {

		print WAV $outHeaders[$outHeaderCounter]."\t";
	    }
	    print WAV "\n";
	}
    }
    else {

	open (WAV, ">>".$filename) or die "Can't write to ".$filename.":".$!,"\n";
    }

    if (scalar( @{$dbFile{'0'}{'Column_Keys'}}) == 1) { #Any db file should be fine since all have to have the same nr of keys
	
	for (my $outColumnCounter=0;$outColumnCounter<scalar(@outColumns);$outColumnCounter++ ) {
	    
	    if ( defined($allVariants{$chromosomeNumber}{$outColumns[$outColumnCounter]}) ) {

		print WAV $allVariants{$chromosomeNumber}{$outColumns[$outColumnCounter]}, "\t";
	    }
	    else {

		print WAV "-", "\t";
	    }
	}
	print WAV "\n";
    }
    if (scalar( @{$dbFile{'0'}{'Column_Keys'}}) == 2) { #Any db file should be fine since all have to have the same nr of keys

	for (my $position=0;$position<scalar( @{$allVariantsChromosomeSorted{$chromosomeNumber} } );$position++)  { #For all pos per chromosome	
	    
	    my $secondKey = $allVariantsChromosomeSorted{$chromosomeNumber}[$position]; #pos keys to hash from sorted arrray
	    
	    for (my $outColumnCounter=0;$outColumnCounter<scalar(@outColumns);$outColumnCounter++ ) {
		
		if ( defined($allVariants{$chromosomeNumber}{$secondKey}{$outColumns[$outColumnCounter]}) ) {

		    print WAV $allVariants{$chromosomeNumber}{$secondKey}{$outColumns[$outColumnCounter]}, "\t";
		}
		else {

		    print WAV "-", "\t";
		}
	    }
	    print WAV "\n";
	}
    }
    if (scalar( @{$dbFile{'0'}{'Column_Keys'}}) == 3) { #Any db file should be fine since all have to have the same nr of keys
	
	for (my $position=0;$position<scalar( @{$allVariantsChromosomeSorted{$chromosomeNumber} } );$position++)  { #For all pos per chr	
	    
	    my $secondKey = $allVariantsChromosomeSorted{$chromosomeNumber}[$position]; #pos keys to hash from sorted arrray
	    
	    for my $thirdKey (keys % {$allVariants{$chromosomeNumber}{$secondKey} }) {
		
		for (my $outColumnCounter=0;$outColumnCounter<scalar(@outColumns);$outColumnCounter++ ) {
		    
		    if ( defined($allVariants{$chromosomeNumber}{$secondKey}{$thirdKey}{$outColumns[$outColumnCounter]}) ) {
	
			print WAV $allVariants{$chromosomeNumber}{$secondKey}{$thirdKey}{$outColumns[$outColumnCounter]}, "\t";
		    }
		    else {
			
			print WAV "-", "\t";
		    }
		}
		print WAV "\n";   
	    }
	}
    }
    if (scalar( @{$dbFile{'0'}{'Column_Keys'}}) == 4) { #Any db file should be fine since all have to have the same nr of keys
	
	for (my $position=0;$position<scalar( @{$allVariantsChromosomeSorted{$chromosomeNumber} } );$position++)  { #For all pos per chr	
	    
	    my $secondKey = $allVariantsChromosomeSorted{$chromosomeNumber}[$position]; #pos keys to hash from sorted arrray

	    for my $thirdKey (keys % {$allVariants{$chromosomeNumber}{$secondKey} }) {
		
		for my $fourthKey (keys % {$allVariants{$chromosomeNumber}{$secondKey}{$thirdKey} }) {
		    
		    for (my $outColumnCounter=0;$outColumnCounter<scalar(@outColumns);$outColumnCounter++ ) {
			
			if ( defined($allVariants{$chromosomeNumber}{$secondKey}{$thirdKey}{$fourthKey}{$outColumns[$outColumnCounter]}) ) { #Seperate zero from undef
			
			    print WAV $allVariants{$chromosomeNumber}{$secondKey}{$thirdKey}{$fourthKey}{$outColumns[$outColumnCounter]}, "\t";
			}
			else {
			    print WAV "-", "\t";
			}
		    }
		    print WAV "\n";
		}
	    }
	}
    }
    close (WAV);

    if ($prefixChromosomes == 0) {

	print STDOUT "Finished Writing Master file for: chromosome".$chromosomeNumber,"\n";
    }
    else {

	print STDOUT "Finished Writing Master file for: ".$chromosomeNumber,"\n";
    }
}

sub WriteAll {
#Prints tab separated file of all collected db file info 
    
    my $fileName = $_[0];

    open (WAV, ">".$fileName) or die "Can't write to ".$fileName.":".$!, "\n";

    if (@metaData) { #Print metaData if supplied

	for (my $metaDataCounter=0;$metaDataCounter<scalar(@metaData);$metaDataCounter++) {

	    print WAV $metaData[$metaDataCounter],"\n";
	}
    }
    if ($outInfo == 1 ) { #Print header if supplied

	print WAV "#";

	for (my $outHeaderCounter=0;$outHeaderCounter<scalar(@outHeaders);$outHeaderCounter++) {

	    print WAV $outHeaders[$outHeaderCounter], "\t";
	}
	print WAV "\n";
    }
    if (scalar( @{$dbFile{'0'}{'Column_Keys'}}) == 1) { #Any db file should be fine since all have to have the same nr of keys
	
	for my $firstKey (keys %allVariants) { #All firstKeys
	    
	    for (my $outColumnCounter=0;$outColumnCounter<scalar(@outColumns);$outColumnCounter++ ) {
		
		if ( defined($allVariants{$firstKey}{$outColumns[$outColumnCounter]}) ) {

		    print WAV $allVariants{$firstKey}{$outColumns[$outColumnCounter]}, "\t";
		}
		else {

		    print WAV "-", "\t";
		}
	    }
	    print WAV "\n";
	}
    }
    if (scalar( @{$dbFile{'0'}{'Column_Keys'}}) == 2) { #Any db file should be fine since all have to have the same nr of keys
	
	for my $firstKey (keys %allVariants) { #All firstKeys
	    
	    for my $secondKey (keys % {$allVariants{$firstKey} }) {
		
		for (my $outColumnCounter=0;$outColumnCounter<scalar(@outColumns);$outColumnCounter++ ) {
		    
		    if ( defined($allVariants{$firstKey}{$secondKey}{$outColumns[$outColumnCounter]}) ) {
		
			print WAV $allVariants{$firstKey}{$secondKey}{$outColumns[$outColumnCounter]}, "\t";
		    }
		    else {
			
			print WAV "-", "\t";
		    }
		}
		print WAV "\n";
	    }
	}
    }
    if (scalar( @{$dbFile{'0'}{'Column_Keys'}}) == 3) { #Any db file should be fine since all have to have the same nr of keys

	for my $firstKey (keys %allVariants) { #All firstKeys
	    
	    for my $secondKey (keys % {$allVariants{$firstKey} }) {
		
		for my $thirdKey (keys % {$allVariants{$firstKey}{$secondKey} }) {
		    
		    for (my $outColumnCounter=0;$outColumnCounter<scalar(@outColumns);$outColumnCounter++ ) {
			
			if ( defined($allVariants{$firstKey}{$secondKey}{$thirdKey}{$outColumns[$outColumnCounter]}) ) {
			
			    print WAV $allVariants{$firstKey}{$secondKey}{$thirdKey}{$outColumns[$outColumnCounter]}, "\t";
			}
			else {
			    
			    print WAV "-", "\t";
			}
		    }
		    print WAV "\n";   
		}
	    }
	}
    }
    if (scalar( @{$dbFile{'0'}{'Column_Keys'}}) == 4) { #Any db file should be fine since all have to have the same nr of keys
	
	for my $firstKey (keys %allVariants) { #All firstKeys
	    
	    for my $secondKey (keys % {$allVariants{$firstKey} }) {
		
		for my $thirdKey (keys % {$allVariants{$firstKey}{$secondKey} }) {
		    
		    for my $fourthKey (keys % {$allVariants{$firstKey}{$secondKey}{$thirdKey} }) {
			
			for (my $outColumnCounter=0;$outColumnCounter<scalar(@outColumns);$outColumnCounter++ ) {
			    
			    if ( defined($allVariants{$firstKey}{$secondKey}{$thirdKey}{$fourthKey}{$outColumns[$outColumnCounter]}) ) {
		
				print WAV $allVariants{$firstKey}{$secondKey}{$thirdKey}{$fourthKey}{$outColumns[$outColumnCounter]}, "\t";
			    }
			    else {
				
				print WAV "-", "\t";
			    }
			}
			print WAV "\n";
		    }
		}
	    }
	}
    }
    close (WAV);
    print STDOUT "Finished Writing Master file for: ".$fileName,"\n";
} 

####
#Merged Sub Routines
####


sub ReadInfileMerge {
#Reads the first db file per chromosome.

    my $DbFileName = $_[0];
    my $DbFileNumber = $_[1];
    my $chromosomeNumber = $_[2];
    my $nextChromosome = $_[3];

    open(RIFM, "<".$DbFileName) or die "Can't open ".$DbFileName.":".$!, "\n";
    
    if ( defined($dbFilePos{$DbFileNumber}) ) { #if file has been searched previously

	seek(RIFM, $dbFilePos{$DbFileNumber},0) or die "Couldn't seek to ".$dbFilePos{$DbFileNumber}." in ".$DbFileName.":".$!, "\n"; #Seek to binary position in file where we left off
    }    

    while (<RIFM>) {
	
	chomp $_;

	if ($.==1) {
	    
	    print STDERR "Started at line $.\n";
	    my $position = tell(RIFM);
	    if ( defined ($dbFilePos{$DbFileNumber}) ) {
		
		print STDERR "Started at position ", $dbFilePos{$DbFileNumber}, "\n";
	    }
	    else {
	
		print STDERR "Started at position ",  $position, "\n";
	    }
	}
	if (m/^\s+$/) {		# Avoid blank lines
	    next;
	}
	if (m/^#/) {		# Avoid #
	    next;
	}		
	if ( $_ =~/^(\S+)/ ) {

	    my @lineElements = split($dbFile{$DbFileNumber}{'Separator'},$_); #Loads line

	    if ($nextChromosome && $lineElements[$dbFile{$DbFileNumber}{'Chr_column'}] eq $nextChromosome) { #If next chromosome is found return (Since all numerically infiles are sorted this is ok)

		print STDOUT "Finished Reading chromosome".$chromosomeNumber." in Infile ".$DbFileName,"\n";
		$dbFilePos{$DbFileNumber} = tell(RIFM); # Save  binary position in file to enable seek when revisiting e.g. next chromosome
		print STDERR "Ended at position ".$dbFilePos{$DbFileNumber}." in ".$dbFile{$DbFileNumber}{'File'}, "\n";
		close(RIFM);
		return;
	    }
	    if ($lineElements[$dbFile{$DbFileNumber}{'Chr_column'}] eq $chromosomeNumber) { #Only correct chromosome

		if (scalar( @{$dbFile{0}{'Column_Keys'}}) == 1) {
		    
		    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$DbFileNumber}{'Column_To_Extract'}});$extractColumnsCounter++) {
			
			my $columnIdRef = \($DbFileNumber."_".$dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter]);
			$allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] }{$$columnIdRef} = $lineElements[ $dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter] ];
		    }
		    push (@{$unSorted{$lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ]}}, $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ]);
		}
		if (scalar( @{$dbFile{0}{'Column_Keys'}}) == 2) {
		    
		    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$DbFileNumber}{'Column_To_Extract'}});$extractColumnsCounter++) {
			
			my $columnIdRef = \($DbFileNumber."_".$dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter]);
			$allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] }{$$columnIdRef} = $lineElements[ $dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter] ];
		    }
		    push (@{$unSorted{$lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ]}}, $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ]);
		}
		if (scalar( @{$dbFile{0}{'Column_Keys'}}) == 3) {

		    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$DbFileNumber}{'Column_To_Extract'}});$extractColumnsCounter++) {
			
			my $columnIdRef = \($DbFileNumber."_".$dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter]);
			$allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ] }{$$columnIdRef} = $lineElements[ $dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter] ];	    
		    }
		    push (@{$unSorted{$lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ]}}, $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ]);
		} 
		if (scalar( @{$dbFile{0}{'Column_Keys'}}) == 4) {
		    
		    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$DbFileNumber}{'Column_To_Extract'}});$extractColumnsCounter++) {
			
			my $columnIdRef = \($DbFileNumber."_".$dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter]);
			$allVariants{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[2] ] }{ $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[3] ]}{$$columnIdRef} = $lineElements[ $dbFile{$DbFileNumber}{'Column_To_Extract'}[$extractColumnsCounter] ];
		    }
		    push (@{$unSorted{$lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[0] ]}}, $lineElements[ $dbFile{$DbFileNumber}{'Column_Keys'}[1] ]); 	
		}
	    }
	}
    }
    close(RIFM);
    print STDOUT "Finished Reading key Db file: ".$DbFileName." in Merged Mode","\n";
    return;
}

sub ReadDbFilesMerge {
#Reads all db files collected from Db master file, except first db file and db files with features "large", "range". These db files are handled by different subroutines. Requires numerically sorted infiles.

    my $chromosomeNumber = $_[0];
    my $nextChromosome = $_[1];
    
    for (my $dbFileNr=1;$dbFileNr<$dbFileCounter;$dbFileNr++) { #All db files (in order of appearance in $db) except first db which has already been handled	
	
	#Read files per chromosome (i.e. multiple times)
	open(DBFM, "<".$dbFile{$dbFileNr}{'File'}) or die "Can't open ".$dbFile{$dbFileNr}{'File'}.":".$!, "\n";    

	if ( defined($dbFilePos{$dbFileNr}) ) { #if file has been searched previously

	    seek(DBFM, $dbFilePos{$dbFileNr},0) or die "Couldn't seek to ".$dbFilePos{$dbFileNr}.":".$!, "\n"; #Seek to binary position in file where we left off
	}

	while (<DBFM>) {
	    
	    chomp $_; #Remove newline
	    
	    if ($.==1) {

		print STDERR "Started at line ".$., "\n";
		my $position = tell(DBFM);
		if ( defined ($dbFilePos{$dbFileNr}) ) {
		
		    print STDERR "Started at position ".$dbFilePos{$dbFileNr}, "\n";
		}
		else {
		
		    print STDERR "Started at position ".$position, "\n";
		}
	    }
	    if (m/^\s+$/) {		# Avoid blank lines
		next;
	    }
	    if (m/^#/) {		# Avoid #
		next;
	    }		
	    if ( $_ =~/^(\S+)/) {
		
		my @lineElements = split($dbFile{$dbFileNr}{'Separator'},$_); #Splits columns on separator and loads line
		
		if ($nextChromosome && $lineElements[$dbFile{$dbFileNr}{'Chr_column'}] eq $nextChromosome) { #If next chromosome is found return (Since all numerically infiles are sorted this is ok)
		    
		    $dbFilePos{$dbFileNr} = tell(DBFM); # Save  binary position in file to enable seek when revisiting e.g. next chromsome
		    print STDERR "Ended at position ".$dbFilePos{$dbFileNr}." in ".$dbFile{$dbFileNr}{'File'}, "\n";		    
		    close(DBFM);
		    last;
		}
		if ( $lineElements[$dbFile{$dbFileNr}{'Chr_column'}] eq $chromosomeNumber) { #If chromosome number and chromosome in line match - go ahead and process
		    
		    #Depending on the number of column keys supplied by user in db master file . 
		    if (scalar( @{$dbFile{$dbFileNr}{'Column_Keys'}}) == 1) {
		       
			for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$dbFileNr}{'Column_To_Extract'}});$extractColumnsCounter++) {
			    
			    my $columnIdRef = \($dbFileNr."_".$dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter]);
			    $allVariants{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[0] ] }{$$columnIdRef} = $lineElements[ $dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter] ];
			}
		    }
		    if (scalar( @{$dbFile{$dbFileNr}{'Column_Keys'}}) == 2) {

			for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$dbFileNr}{'Column_To_Extract'}});$extractColumnsCounter++) {
			    
			    my $columnIdRef = \($dbFileNr."_".$dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter]);
			    $allVariants{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[1] ] }{$$columnIdRef}= $lineElements[ $dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter] ];
			}
		    }
		    if (scalar( @{$dbFile{$dbFileNr}{'Column_Keys'}}) == 3) {
			
			for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$dbFileNr}{'Column_To_Extract'}});$extractColumnsCounter++) {
			    
			    my $columnIdRef = \($dbFileNr."_".$dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter]);
			    $allVariants{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[0] ]}{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[1] ]}{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[2] ]}{$$columnIdRef}= $lineElements[ $dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter] ];
			}
		    }
		    if (scalar( @{$dbFile{$dbFileNr}{'Column_Keys'}}) == 4) {
			
			for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$dbFileNr}{'Column_To_Extract'}});$extractColumnsCounter++) {
			    
			    my $columnIdRef = \($dbFileNr."_".$dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter]);
			    $allVariants{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[0] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[1] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[2] ] }{ $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[3] ]}{$$columnIdRef}= $lineElements[ $dbFile{$dbFileNr}{'Column_To_Extract'}[$extractColumnsCounter] ];
			    
			}
			push (@{$unSorted{$lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[0] ]}}, $lineElements[ $dbFile{$dbFileNr}{'Column_Keys'}[1] ]);
		    }
		}
	    } 	
	}
	close(DBFM);
	print STDOUT "Finished Reading chromsome".$chromosomeNumber." in Infile: ".$dbFile{$dbFileNr}{'File'},"\n";
    }
}

sub SortAllVariantsMergeST {
#Creates an array of all position which are unique and in sorted ascending order  
    
    my $firstKey = $_[0];
    my %seen = (); 

#Unique entries only
    @{$allVariantsChromosomeUnique{$firstKey} } = grep { ! $seen{$_} ++ } @{$unSorted{$_[0]} };
#Sort using the Schwartzian Transform 
    @{$allVariantsChromosomeSorted{$firstKey} } =  map { $_->[0] }
    sort { $a->[0] <=> $b->[0] }
    map { [$_] }
    @{$allVariantsChromosomeUnique{$firstKey} };

    @{$unSorted{$_[0] } }=();  #Reset for this chromosome
}

sub SortAllVariantsMerge {
#Creates an array of all position   
    
    my $firstKey = $_[0];
    my %seen = (); 

    if (scalar( @{$dbFile{0}{'Column_Keys'}}) > 1) {
	
	for my $secondKey (keys %{ $allVariants{$firstKey} } )  { #For all secondkeys
	    
	    push ( @{$allVariantsChromosome{$firstKey} },$secondKey );
	}
	@{$allVariantsChromosomeUnique{$firstKey} } = grep { ! $seen{$_} ++ } @{$allVariantsChromosome{$firstKey} }; #Unique entries only 
	@{$allVariantsChromosomeSorted{$firstKey} } = sort { $a <=> $b } @{ $allVariantsChromosomeUnique{$firstKey} }; #Sorts keys to be able to print sorted table later 
	print STDOUT "Sorted all non overlapping entries\n";
    }
}

sub WriteAllVariantsMerge {
#Prints tab separated file of all collected db file info in ascending order.
#$_[0] = filename
#$_[1] = chromsome number
    
    my $fileName = $_[0];
    my $chromosome = $_[1];

    if ( ($chromosome eq 1) || ($chromosome eq "chr1") ) {

	open (WAVM, ">".$fileName) or die "Can't write to ".$fileName.":".$!,"\n";

	if (@metaData) { #Print metaData if supplied

	    for (my $metaDataCounter=0;$metaDataCounter<scalar(@metaData);$metaDataCounter++) {

		print WAV $metaData[$metaDataCounter],"\n";
	    }
	}
	if ($outInfo == 1 ) { #Print header if supplied

	    print WAVM "#";

	    for (my $outHeaderCounter=0;$outHeaderCounter<scalar(@outHeaders);$outHeaderCounter++) {

		print WAVM $outHeaders[$outHeaderCounter], "\t";
	    }
	    print WAVM "\n";
	}
    }
    else {

	open (WAVM, ">>".$fileName) or die "Can't write to ".$fileName.":".$!, "\n";
    }

    if (scalar( @{$dbFile{'0'}{'Column_Keys'}}) == 1) { #Any db file should be fine since all have to have the same nr of keys
	
	for (my $outColumnCounter=0;$outColumnCounter<scalar(@outColumns);$outColumnCounter++ ) {
	    
	    if ( defined($allVariants{$chromosome}{$outColumns[$outColumnCounter]}) ) {

		print WAVM $allVariants{$chromosome}{$outColumns[$outColumnCounter]}, "\t";
	    }
	    else {

		print WAVM "-", "\t";
	    }
	}
	print WAVM "\n";
    }
    if (scalar( @{$dbFile{'0'}{'Column_Keys'}}) == 2) { #Any db file should be fine since all have to have the same nr of keys

	for (my $position=0;$position<scalar( @{$allVariantsChromosomeSorted{$chromosome} } );$position++)  { #For all position per chromosome	
	    
	    my $secondKey = $allVariantsChromosomeSorted{$chromosome}[$position]; #position keys to hash from sorted arrray
	    
	    for (my $outColumnCounter=0;$outColumnCounter<scalar(@outColumns);$outColumnCounter++ ) {
		
		if ( defined($allVariants{$chromosome}{$secondKey}{$outColumns[$outColumnCounter]}) ) {
	
		    print WAVM $allVariants{$chromosome}{$secondKey}{$outColumns[$outColumnCounter]}, "\t";
		}
		else {
		    
		    print WAVM "-", "\t";
		}
	    }
	    print WAVM "\n";
	}
    }
    if (scalar( @{$dbFile{'0'}{'Column_Keys'}}) == 3) { #Any db file should be fine since all have to have the same nr of keys
	
	for (my $position=0;$position<scalar( @{$allVariantsChromosomeSorted{$chromosome} } );$position++)  { #For all position per chromosome	
	    
	    my $secondKey = $allVariantsChromosomeSorted{$chromosome}[$position]; #position keys to hash from sorted arrray
	    
	    for my $thirdKey (keys % {$allVariants{$chromosome}{$secondKey} }) {
		
		for (my $outColumnCounter=0;$outColumnCounter<scalar(@outColumns);$outColumnCounter++ ) {
		    
		    if ( defined($allVariants{$chromosome}{$secondKey}{$thirdKey}{$outColumns[$outColumnCounter]}) ) {
		
			print WAVM $allVariants{$chromosome}{$secondKey}{$thirdKey}{$outColumns[$outColumnCounter]}, "\t";
		    }
		    else {
			
			print WAVM "-", "\t";
		    }
		}
		print WAVM "\n";   
	    }
	}
    }
    if (scalar( @{$dbFile{'0'}{'Column_Keys'}}) == 4) { #Any db file should be fine since all have to have the same nr of keys
	
	for (my $position=0;$position<scalar( @{$allVariantsChromosomeSorted{$chromosome} } );$position++)  { #For all position per chromsome	
	    
	    my $secondKey = $allVariantsChromosomeSorted{$chromosome}[$position]; #position keys to hash from sorted arrray

	    for my $thirdKey (keys % {$allVariants{$chromosome}{$secondKey} }) {
		
		for my $fourthKey (keys % {$allVariants{$chromosome}{$secondKey}{$thirdKey} }) {
		    
		    print WAVM $chromosome."\t".$secondKey."\t".$thirdKey."\t".$fourthKey."\t";

		    for (my $outColumnCounter=0;$outColumnCounter<scalar(@outColumns);$outColumnCounter++ ) {
			
			if ( defined($allVariants{$chromosome}{$secondKey}{$thirdKey}{$fourthKey}{$outColumns[$outColumnCounter]}) ) {

			    print WAVM $allVariants{$chromosome}{$secondKey}{$thirdKey}{$fourthKey}{$outColumns[$outColumnCounter]}, "\t";
			}
			else {

			    print WAVM "-", "\t";
			}
		    }
		    print WAVM "\n";
		}
	    }
	}
    }
    close (WAVM);
    print STDOUT "Finished Writing Master file for: chromsome".$chromosome,"\n";
}


###Decommissoned###
sub TabixQuery {
##Returns 

    my $chromosomeRef = $_[0];
    my $chromosomeStartRef = $_[1];
    my $chromosomeStopRef = $_[2];
    my $arrayColumnsRef = $_[3];    
    my $tabix = $_[4];
    my $dbFileNrRef = $_[5];

print $$chromosomeRef, "\t", $$chromosomeStartRef, "\t", $$chromosomeStopRef, "\n";
my $var = $tabix->read(

  $tabix->query( $$chromosomeRef, $$chromosomeStartRef, $$chromosomeStopRef) 

); #Starts from position but not with the actual position so need to subtrackt -1 from start
    if (defined($var)) {
	print $var, "\n";
	my @tabixReturnArray = split('\t', $var);
	my @array;
	if ( $allVariants{$$chromosomeRef}{$$chromosomeStartRef}{$$chromosomeStopRef}{$tabixReturnArray[ $dbFile{$$dbFileNrRef}{'Column_Keys'}[3] ]}) {
	
	    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar( @{$dbFile{$$dbFileNrRef}{'Column_To_Extract'}});$extractColumnsCounter++) {
		
		my $columnIdRef = \($$dbFileNrRef."_".$dbFile{$$dbFileNrRef}{'Column_To_Extract'}[$extractColumnsCounter]);
		$allVariants{$$chromosomeRef }{$$chromosomeStartRef}{$$chromosomeStopRef }{$tabixReturnArray[ $dbFile{$$dbFileNrRef}{'Column_Keys'}[3] ]}{$$columnIdRef} = $tabixReturnArray[ $dbFile{$$dbFileNrRef}{'Column_To_Extract'}[$extractColumnsCounter] ];
	    }
	}
    }
}
