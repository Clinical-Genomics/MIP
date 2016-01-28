#!/usr/bin/env perl

use strict;
use warnings;

use Pod::Usage;
use Pod::Text;
use Getopt::Long;


use vars qw($USAGE);

BEGIN {
    $USAGE =
        qq{afMax.pl - [vcf]
               -h/--help Display this help message
               -v/--version Display version};
}

my ($infile, $version, $help) = ("");

my $maxAFVersion = "0.0.1";

## Enables cmd "maxAF.pl" to print usage help 
if(scalar(@ARGV) == 0) {

    print STDOUT $USAGE, "\n";
    exit;
}
elsif ( (defined($ARGV)) && ($ARGV[0]!~/^-/) ) { #Collect potential infile - otherwise read from STDIN
    
    $infile = $ARGV[0];
}


GetOptions('h|help' => sub { print STDOUT $USAGE, "\n"; exit;},  #Display help text
	   'v|version' => sub { print STDOUT "\nmaxAF.pl ".$maxAFVersion, "\n\n"; exit;},  #Display version number
    );


###
#MAIN
###

&ReadInfile($infile);



###
#Sub Routines
###


sub ReadInfile {

    my $infile = $_[0];
    my $reference = $_[1];

    while (<>) {
	
	my $variantLine;
	my $sampleIDInfo;

	chomp $_;  # Remove newline
	
	if (m/^\s+$/) {	# Avoid blank lines
	    next;
	}
	if ($_=~/^##(\S+)=/) { # MetaData

	    print $_, "\n";  #Write to STDOUT
	    next;
	}
	elsif ($_=~/^#CHROM/) {

	    ## Write vcf Header
	    print '##INFO=<ID=MAX_AF,Number=A,Type=Float,Description="Maximum Allele Frequency for variant, across populations ">', "\n";   
	    print $_, "\n";  #Write to STDOUT
	    next;
	}
	else {

	    my $maxAF = 0;
	    my @lineElements = split("\t",$_);
	    
	    for (my $lineElementsCounter=0;$lineElementsCounter<scalar(@lineElements);$lineElementsCounter++) { #Add until INFO field
			
		if ($lineElementsCounter < 7) { #Save fields until INFO field

		    $variantLine .= $lineElements[$lineElementsCounter]."\t";
		}
		elsif ($lineElementsCounter > 7) { #Save GT:PL: and sample(s) GT Call fields and add to proper line last
		    
		    if ($lineElementsCounter == (scalar(@lineElements) - 1)) {

			$sampleIDInfo .= $lineElements[$lineElementsCounter];
		    }
		    else {

			$sampleIDInfo .= $lineElements[$lineElementsCounter]."\t";
		    }
		}
	    }

	    my @keyValues = split(/;/, $lineElements[7]); #Split INFO field to key=value items
	    
	    for my $element (@keyValues) {
		
		my @keys = split("=", $element);  #Key = 0 and value = 1
		
		if ($keys[0]=~/_AF$/) {  #Allele frequency key
		    
		    if ($keys[1] > $maxAF) {
			
			$maxAF = $keys[1];
		    }
		}
	    }
	    if ($maxAF != 0) {

		$variantLine .= join(";", "MAX_AF=".$maxAF, $lineElements[7]);
	    }
	    if (defined($sampleIDInfo)) {

		$variantLine .= "\t".$sampleIDInfo;
	    }
	    print $variantLine, "\n";  #Write to STDOUT
	}
    }
    close()
}
