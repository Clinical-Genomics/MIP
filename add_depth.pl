#!/usr/bin/perl - w
  
###Copyright 2012 Henrik Stranneheim

##Add the DP for samples lacking a genotype call where at least 1 other sampleID has a GT call.

use strict;
use warnings;
use Pod::Usage;
use Pod::Text;
use Getopt::Long;

use vars qw($USAGE);

BEGIN {
    $USAGE =
	qq{add_depth.pl -i [inFile1] -infnv [inFilesnv...n]  -s [sampleID...n] -o [outFile.txt]
               -i/--inFile inFile
               -infnv/--inFilesNoVariants Novariant inFile(s), comma sep (Same order as sampleIDs)
               -s/--sampleIDs Sampleid(s), comma sep (Same order as infnv)
	       -o/--outFile The output file (defaults to add_depth_out.txt)
               -prechr/--prefixChromosomes "chrX" or just "X" (defaults to "X" i.e. no prefix)
               -h/--help Display this help message
               -v/--version Display version
	   };
}

my ($outFile, $inFile, $prefixChromosome, $help, $version) = ("add_depth_out.txt", 0, 0);
my (@inFilesNoVariants, @sampleIDs, @contigs, @metaData);

GetOptions('i|inFile:s'  => \$inFile,
	   'infnv|inFilesNoVariants:s'  => \@inFilesNoVariants, #Comma separated list
	   's|sampleIDs:s'  => \@sampleIDs, #Comma separated list
	   'o|outFile:s'  => \$outFile,
	   'prechr|prefixChromosomes:n'  => \$prefixChromosome,
	   'h|help' => \$help,
	   'v|version' => \$version,
	   );

if($help) {

    print STDOUT $USAGE, "\n";
    exit;
}
if($version) {

    print STDOUT "\nAdd_depth.pl v1.0\n\n";
    exit;
}

if ($inFile eq 0) {

    print STDERR "\n";
    print STDOUT $USAGE, "\n";
    print STDERR "Must supply an annovar_all.txt infile", "\n\n";
    exit;
}

if (@inFilesNoVariants == 0) {

    print STDERR "\n";
    print STDOUT $USAGE, "\n";
    print STDERR "Must supply infile(s) as comma separeted list", "\n\n";
    exit;
}

@inFilesNoVariants = split(/,/,join(',',@inFilesNoVariants));  # Enables comma separated indir(s)
@sampleIDs = split(/,/,join(',',@sampleIDs));      # Enables comma separated indir(s)

## Set contigs prefix and chromosome names depending on reference used
if ($prefixChromosome == 0) {
    # Ensembl - no prefix and MT
    # Contigs for enhanced speed in collecting information and reducing memory
    # consumption
    @contigs = ("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15",
            "16","17","18","19","20","21","22","X","Y","MT");
}
else {
    # Refseq - prefix and M
    @contigs = ("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
            "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
            "chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM");
}

my $header;  # To preserve header information
my (@allVariants);
my (%allVariants, %allVariantsContig, %allVariantsContigUnique,
    %allVariantsContigSorted);
my (%sampleVariants, %column);

###
#Main
###

for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {

    # Read all positions per sampleID that lacks variation
    &ReadNonVariant($inFilesNoVariants[$sampleIDCounter],$sampleIDs[$sampleIDCounter]);
    &FindColumn($inFile,$sampleIDs[$sampleIDCounter]);  # Collect column of sampleID
}

&ReadInFile($inFile);  # Read annovar_all.txt master file
&SortAllVariants();     # Sorts all variants
&WriteAddedDepth($outFile);  # Write all variants

###
#Sub routines
###

sub FindColumn {
#Finds sampleID columns number in annovar_all.txt file and then breaks
    
    my $fileName = $_[0];
    my $sampleID = $_[1];

    open(FC, "<".$fileName) or die "Can't open ".$fileName.":".$!, "\n";    
    
    while (<FC>) {
	
	chomp $_; #Remove newline
	
	if (m/^\s+$/) {		# Avoid blank lines
            next;
        }	
	if (m/^\#/) {		# Avoid #
            next;
        }			
	if ( $_ =~/$sampleID/ ) {
	    
	    my @lineElements = split("\t",$_);	    #Loads variant calls
	    
	    for (my $columnCounter=0;$columnCounter<scalar(@lineElements);$columnCounter++) { 
		
		if ( $lineElements[$columnCounter]=~/$sampleID/ ) {

		    $column{$sampleID} = $columnCounter;
		    last; #Break
		}
	    }
	    last;
	} 	
    }
    close(FC);
    print STDOUT "SampleID: ".$sampleID." has column: ".$column{$sampleID},"\n";
}

sub ReadNonVariant {
#Reads sampleID_mpileup_nonvariants.txt file to be able to add depth later

    my $fileName = $_[0];
    my $sampleID = $_[1];
    
    open(NV, "<".$fileName) or die "Can't open ".$fileName.":".$!, "\n";    
    
    while (<NV>) {
	
	chomp $_; #Remove newline
	
	if (m/^\s+$/) {		# Avoid blank lines
            next;
        }	
	if (m/^\#/) {		# Avoid #
            next;
        }			
	if ( $_ =~/\S+/ ) {
	    
	    my @lineElements = split("\t",$_);	    #Loads nonvariant calls
	    $sampleVariants{$sampleID}{$lineElements[0]}{$lineElements[1]} = $lineElements[2]; # Hash{sampleID}{contigs}{pos} = [depth]
	}
    } 	
    close(NV);
    print STDOUT "Finished Reading Infile ".$fileName,"\n";
}

sub ReadInFile {
#Reads infile
    
    my $fileName = $_[0];

    open(RIF, "<".$fileName) or die "Can't open ".$fileName.":".$!, "\n";    
    
    while (<RIF>) {
	
	chomp $_; #Remove newline
	
	if (m/^\s+$/) {		# Avoid blank lines
            next;
        }
	if ($_=~/^##/) {#MetaData

	    push(@metaData, $_); #Save metadata string
	    next;
        }	
	if (m/^\#/) { #Catch header if present
           
	    $header = $_;
	    next;
        }			
	if ( $_ =~/\S+/ ) {
	    
	    my @lineElements = split("\t",$_);	    #Loads variant calls
	    $allVariants{$lineElements[0]}{$lineElements[1]}{$lineElements[4]} = [@lineElements]; # Hash{contigs}{pos}{variant}, all variants non overlapping and array [contigs->unknown] All info starting from contigs

	    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) { #Check if any nonvariants
		
		if ($sampleVariants{$sampleIDs[$sampleIDCounter]}{$lineElements[0]}{$lineElements[1]}) {
		    
		    $allVariants{$lineElements[0]}{$lineElements[1]}{$lineElements[4]}[ $column{$sampleIDs[$sampleIDCounter]} ] .= ":DP=".$sampleVariants{$sampleIDs[$sampleIDCounter]}{$lineElements[0]}{$lineElements[1]}; #Add DP info the original GT call
		}
	    }
	}
    } 	
    close(RIF);
    print STDOUT "Finished Reading Infile ".$fileName,"\n";
}

sub SortAllVariants {
#Creates an array of all position which are unique and in sorted ascending order
    
    for my $contigs (keys %allVariants)  { #For all contigs
	
	for my $position (keys %{ $allVariants{$contigs} } )  { #For all positions
    
	    for my $variant (keys % { $allVariants{$contigs}{$position} })  { #For all variants
		push ( @{$allVariantsContig{$contigs} },$position );
	    }
	}
	my %seen = (); @{$allVariantsContigUnique{$contigs} } = grep { ! $seen{$_} ++ } @{$allVariantsContig{$contigs} }; #Unique entries only
	@{$allVariantsContigSorted{$contigs} } = sort { $a <=> $b } @{ $allVariantsContigUnique{$contigs} }; #Sorts keys to be able to print sorted table later
    }
    print STDOUT "Sorted all non overlapping entries per contig and position\n";
}

sub WriteAddedDepth {
#Write added depth variants to new masterfile
    
    my $fileName = $_[0];

    open (WAD, ">".$fileName) or die "Can't write to ".$fileName.":".$!, "\n";
    
    if (@metaData) { #Print metaData if supplied

	for (my $metaDataCounter=0;$metaDataCounter<scalar(@metaData);$metaDataCounter++) {

	    print WAD $metaData[$metaDataCounter],"\n";
	}
    }
    if ($header) { #Print original header

	print WAD $header, "\n";
    } 
    for (my $contigCounter=0;$contigCounter<scalar(@contigs);$contigCounter++)  { #For all contigs

	if ($allVariantsContigSorted{$contigs[$contigCounter]}) { #If present in infile

	    for (my $i=0;$i<scalar( @{$allVariantsContigSorted{$contigs[$contigCounter]} } );$i++)  { #For all pos per contigs	
		
		my $positionRef = \($allVariantsContigSorted{$contigs[$contigCounter]}[$i]); #pos keys to hash from sorted arrray
		
		for my $variant (keys % { $allVariants{$contigs[$contigCounter]}{$$positionRef} })  { #For all variants
		    
		    for (my $variants=0;$variants<scalar( @{ $allVariants{$contigs[$contigCounter]}{$$positionRef}{$variant} } );$variants++)  {

			print WAD $allVariants{$contigs[$contigCounter]}{$$positionRef}{$variant}[$variants], "\t";
		    }
		    print WAD "\n";
		}
	    }
	}	
    }
    close(WAD);
    print STDOUT "Wrote all variants with added depth (DP) to ".$fileName,"\n";
    return;
}


