#!/usr/bin/perl - w
  
###Copyright 2012 Henrik Stranneheim

=head1 SYNOPSIS

add_depth.pl -i [infile1] -infnv [infilesnv...n]  -sid [sampleID...n] -o [outfile.txt]

=head2 COMMANDS AND OPTIONS

-i/--infile Infile, comma sep

-infnv/--infileNoVariant Novariant infile(s)

-sid/--sampleid Sampleid(s), comma sep

-o/--outfile The output file (defaults to annovar_master.txt)

-prechr/--prefix_chromosomes "chrX" or just "X" (defaults to "X" i.e. no prefix)

=head3 I/O

Input format (annovar_all_variants(tab sep) - added depth )

Output format (tab separate list)

=cut

use strict;
use warnings;
use Pod::Usage;
use Pod::Text;
use Getopt::Long;

use vars qw($USAGE);

BEGIN {
    $USAGE =
	qq{add_depth.pl -i [infile1] -infnv [infilesnv...n]  -sid [sampleID...n] -o [outfile.txt]
               -i/--infile infile
               -infnv/--infileNoVariant Novariant infile(s), comma sep (Same order as sid)
               -sid/--sampleid Sampleid(s), comma sep (Same order as infnv)
	           -o/--outfile The output file (defaults to annovar_master.txt)
               -prechr/--prefix_chromosomes "chrX" or just "X" (defaults to "X" i.e. no prefix)
	   };
}

my ($of, $inf, $prechr, $help) = ("annovar_master.txt", 0, 0);
my (@infnv, @sid, @chr, @metaData);

GetOptions('i|infile:s'  => \$inf,
	   'infnv|infileNoVariant:s'  => \@infnv, #Comma separated list
	   'sid|sampleid:s'  => \@sid, #Comma separated list
	   'o|outfile:s'  => \$of,
	   'prechr|prefix_chromosomes:n'  => \$prechr,
	   'h|help' => \$help,
	   );

die $USAGE if( $help );

if ($inf eq 0) {
   my $verbosity = 2;
   print"\n";
   pod2usage({-message => "Must supply an annovar_all.txt infile.\n",
     -verbose => $verbosity
   });
}

if (@infnv == 0) {
   my $verbosity = 2;
   print"\n";
   pod2usage({-message => "Must supply infile(s) as comma separeted list.\n",
     -verbose => $verbosity
   });
}

@infnv = split(/,/,join(',',@infnv));  # Enables comma separated indir(s)
@sid = split(/,/,join(',',@sid));      # Enables comma separated indir(s)

## Set chr prefix and chromosome names depending on reference used
if ($prechr == 0) {
    # Ensembl - no prefix and MT
    # Chr for enhanced speed in collecting information and reducing memory
    # consumption
    @chr = ("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15",
            "16","17","18","19","20","21","22","X","Y","MT");
}
else {
    # Refseq - prefix and M
    @chr = ("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
            "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
            "chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM");
}

my $header;  # To preserve header information
my (@allVariants, @allVariants_unique, @allVariants_sorted);
my (%allVariants, %allVariants_chr, %allVariants_chr_unique,
    %allVariants_chr_sorted);
my (%sampleVariants, %col);

###
#Main
###

for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {
    # Read all positions per sampleID that lacks variation
    ReadNonVariant($infnv[$sampleid],$sid[$sampleid]);
    FindCol($inf,$sid[$sampleid]);  # Collect column of sampleID
}

ReadAnnovarAll($inf);  # Read annovar_all.txt master file
SortAllVariants();     # Sorts all variants
WriteAddedDepth($of);  # Write all variants

###
#Sub routines
###

sub FindCol {
#Finds sampleID columns number in annovar_all.txt file and then breaks
#$_[0] = filename
#$_[1] = sampleid
    
    open(FC, "<$_[0]") or die "Can't open $_[0]:$!, \n";    
    
    while (<FC>) {
	
	chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
            next;
        }	
	if (m/^\#/) {		# Avoid #
            next;
        }			
	if ( $_ =~/$_[1]/ ) {
	    
	    my @temp = split("\t",$_);	    #Loads variant calls
	    
	    for (my $colnr=0;$colnr<scalar(@temp);$colnr++) { 
		
		if ( $temp[$colnr]=~/$_[1]/ ) {
		    $col{$_[1]} = $colnr;
		    last; #Break
		}
	    }
	    last;
	} 	
    }
    close(FC);
    print STDOUT "SampleID: $_[1] has column: $col{$_[1]}","\n";
    return;
}

sub ReadNonVariant {
#Reads sampleID_mpileup_nonvariants.txt file to be able to add depth later
#$_[0] = filename
#$_[1] = sampleid
    
    open(NV, "<$_[0]") or die "Can't open $_[0]:$!, \n";    
    
    while (<NV>) {
	
	chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
            next;
        }	
	if (m/^\#/) {		# Avoid #
            next;
        }			
	if ( $_ =~/\S+/ ) {
	    
	    my @temp = split("\t",$_);	    #Loads nonvariant calls
	    $sampleVariants{$_[1]}{$temp[0]}{$temp[1]} = $temp[2]; # Hash{sampleID}{chr}{pos}, [depth]
	}
    } 	
    close(NV);
    print STDOUT "Finished Reading Infile $_[0]","\n";
    return;
}

sub ReadAnnovarAll {
#Reads annovar_all.txt file
#$_[0] = filename
    
    open(AALL, "<$_[0]") or die "Can't open $_[0]:$!, \n";    
    
    while (<AALL>) {
	
	chomp $_;
	
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
	    
	    my @temp = split("\t",$_);	    #Loads variant calls
	    $allVariants{$temp[0]}{$temp[1]}{$temp[4]} = [@temp]; # Hash{chr}{pos}{variant}, all variants non overlapping and array [chr->unknown] All info starting from chr
	    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) { #Check if any nonvariants
		
		if ($sampleVariants{$sid[$sampleid]}{$temp[0]}{$temp[1]}) {
		    
		    $allVariants{$temp[0]}{$temp[1]}{$temp[4]}[ $col{$sid[$sampleid]} ] .= ":DP=$sampleVariants{$sid[$sampleid]}{$temp[0]}{$temp[1]}"; #Add DP info the original GT call
		}
	    }
	}
    } 	
    close(AALL);
    print STDOUT "Finished Reading Infile $_[0]","\n";
    return;
}

sub SortAllVariants {
#Creates an array of all position which are unique and in sorted ascending order
    
    for my $chr (keys %allVariants)  { #For all chr
	
	for my $pos (keys %{ $allVariants{$chr} } )  { #For all pos
    
	    for my $variant (keys % { $allVariants{$chr}{$pos} })  { #For all variants
		push ( @{$allVariants_chr{$chr} },$pos );
	    }
	}
	my %seen = (); @{$allVariants_chr_unique{$chr} } = grep { ! $seen{$_} ++ } @{$allVariants_chr{$chr} }; #Unique entries only
	@{$allVariants_chr_sorted{$chr} } = sort { $a <=> $b } @{ $allVariants_chr_unique{$chr} }; #Sorts keys to be able to print sorted table later
    }
    print STDOUT "Sorted all non overlapping entries per chr and position\n";
}

sub WriteAddedDepth {
#Write added depth variants to new masterfile
#$_[0] = filename
    
    open (WAD, ">$_[0]") or die "Can't write to $_[0]: $!\n";
    
    if (@metaData) { #Print metaData if supplied
	for (my $metaDataCounter=0;$metaDataCounter<scalar(@metaData);$metaDataCounter++) {
	    print WAD $metaData[$metaDataCounter],"\n";
	}
    }
    if ($header) { #Print original header
	print WAD $header, "\n";
    } 
    for (my $chrc=0;$chrc<scalar(@chr);$chrc++)  { #For all chr
	
	my $chr = $chr[$chrc];

	if ($allVariants_chr_sorted{$chr}) { #If present in infile

	    for (my $i=0;$i<scalar( @{$allVariants_chr_sorted{$chr} } );$i++)  { #For all pos per chr	
		
		my $pos = $allVariants_chr_sorted{$chr}[$i]; #pos keys to hash from sorted arrray
		
		for my $variant (keys % { $allVariants{$chr}{$pos} })  { #For all variants
		    
		    for (my $variants=0;$variants<scalar( @{ $allVariants{$chr}{$pos}{$variant} } );$variants++)  {
			print WAD $allVariants{$chr}{$pos}{$variant}[$variants], "\t";
		    }
		    print WAD "\n";
		}
	    }
	}	
    }
    close(WAD);
    print STDOUT "Wrote all variants with added depth (DP) to $_[0]","\n";
    return;
}


