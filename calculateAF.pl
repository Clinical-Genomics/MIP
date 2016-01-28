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

my $calculateAFVersion = "0.0.1";

## Enables cmd "calculateAF.pl" to print usage help 
if(scalar(@ARGV) == 0) {

    print STDOUT $USAGE, "\n";
    exit;
}
elsif ( (defined($ARGV)) && ($ARGV[0]!~/^-/) ) { #Collect potential infile - otherwise read from STDIN
    
    $infile = $ARGV[0];
}


GetOptions('h|help' => sub { print STDOUT $USAGE, "\n"; exit;},  #Display help text
	   'v|version' => sub { print STDOUT "\ncalculateAF.pl ".$calculateAFVersion, "\n\n"; exit;},  #Display version number
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

    my %afKey;

    while (<>) {
	
	my $variantLine;
	my $sampleIDInfo;

	chomp $_;  # Remove newline
	
	if (m/^\s+$/) {	# Avoid blank lines
	    next;
	}
	if ($_=~/^##(\S+)=/) { # MetaData

	    my $key = $1;

	    if ($1=~/AN_(\w+)/) {  #Collect AN keys

		$afKey{'AN'}{$1} = $key;
	    }
	    if ($1=~/AC_(\w+)/) {  #Collect AC keys

		$afKey{'AC'}{$1} = $key;
	    }
	    print $_, "\n";  #Write to STDOUT
	    next;
	}
	elsif ($_=~/^#CHROM/) {

	    foreach my $key (keys %{$afKey{'AN'}}) {  #Either "AN" or "AC" keys will do

		## Write vcf Header
		print '##INFO=<ID='.$key.'_AF,Number=A,Type=Float,Description="Allele frequency in the '.$key.' populations calculated from AC and AN, in the range (0,1)">', "\n"
	    }

	    print $_, "\n";  #Write to STDOUT
	    next;
	}
	else {

	    my $calculateAF = 0;
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
		
		if ($keys[0]=~/AN_(\w+)/) {  #"AN" key
		    
		    $afKey{'AN'}{$1} = $keys[1]; #Save value for each population
		}
		if ($keys[0]=~/AC_(\w+)/) {  #"AC" key
		    
		    $afKey{'AC'}{$1} = $keys[1];  #Save value for each population
		}
	    }
	    my @afs;

	    foreach my $key (keys %{$afKey{'AN'}}) {  #To generate AF keys
	
		if( ($afKey{'AC'}{$key}) && ($afKey{'AN'}{$key} != 0)) {

		    $afKey{'AF'}{$key} = $afKey{'AC'}{$key} / $afKey{'AN'}{$key};
		    push(@afs, "AF_".$key."=".$afKey{'AF'}{$key});
		}
	    }

	    $variantLine .= join(";", @afs, $lineElements[7]);

	    if (defined($sampleIDInfo)) {

		$variantLine .= "\t".$sampleIDInfo;
	    }
	    print $variantLine, "\n";  #Write to STDOUT
	}
    }
    close()
}
