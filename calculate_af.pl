#!/usr/bin/env perl

use strict;
use warnings;
use Modern::Perl '2014';
use warnings qw( FATAL utf8 );
use autodie;
use v5.18;  #Require at least perl 5.18
use utf8;  #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );

use Pod::Usage;
use Pod::Text;
use Getopt::Long;
use Params::Check qw[check allow last_error];


use vars qw($USAGE);

BEGIN {
    $USAGE =
        qq{calculate_af.pl - [vcf]
               -h/--help Display this help message
               -v/--version Display version};
}

my ($infile, $version, $help) = ("");

my $calculate_af_version = "0.0.2";

## Enables cmd "calculate_af.pl" to print usage help 
if(scalar(@ARGV) == 0) {

    print STDOUT $USAGE, "\n";
    exit;
}
elsif ( (defined($ARGV)) && ($ARGV[0]!~/^-/) ) { #Collect potential infile - otherwise read from STDIN
    
    $infile = $ARGV[0];
}


GetOptions('h|help' => sub { print STDOUT $USAGE, "\n"; exit;},  #Display help text
	   'v|version' => sub { print STDOUT "\ncalculate_af.pl ".$calculate_af_version, "\n\n"; exit;},  #Display version number
    );


###
#MAIN
###

&read_infile({calculate_af_version => $calculate_af_version,
	     });



###
#Sub Routines
###


sub read_infile {

##read_infile

##Function : Read infile and calculate allele frequency
##Returns  : ""
##Arguments: $calculate_af_version 
##         : $calculate_af_version => Calculate af version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $calculate_af_version;

    my $tmpl = { 
	calculate_af_version => { required => 1, defined => 1, strict_type => 1, store => \$calculate_af_version},
    };
    
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];


    my %af_key;

    while (<>) {
	
	my $variant_line;
	my $sample_id_info;

	chomp $_;  # Remove newline
	
	if (m/^\s+$/) {	# Avoid blank lines
	    next;
	}
	if ($_=~/^##(\S+)=/) { # MetaData

	    my $key = $1;

	    if ($1=~/AN_(\w+)/) {  #Collect AN keys

		$af_key{AN}{$1} = $key;
	    }
	    if ($1=~/AC_(\w+)/) {  #Collect AC keys

		$af_key{AC}{$1} = $key;
	    }
	    say STDOUT $_;
	    next;
	}
	elsif ($_=~/^#CHROM/) {

	    foreach my $key (keys %{$af_key{AN}}) {  #Either "AN" or "AC" keys will do

		## Write vcf Header
		say '##INFO=<ID='.$key.'_AF,Number=A,Type=Float,Description="Allele frequency in the '.$key.' populations calculated from AC and AN, in the range (0,1)">';
	    }

	    ## Write program cmd
	    my ($base, $script) = (`date +%Y%m%d`,`basename $0`);  #Catches current date and script name
	    chomp($base,$script);  #Remove \n;
	    say "##software=<ID=".$script.",Version=".$calculate_af_version.",Date=".$base;

	    say STDOUT $_;
	    next;
	}
	else {

	    my $calculate_af = 0;
	    my @line_elements = split("\t",$_);
	    
	    for (my $line_elements_counter=0;$line_elements_counter<scalar(@line_elements);$line_elements_counter++) { #Add until INFO field
			
		if ($line_elements_counter < 7) { #Save fields until INFO field

		    $variant_line .= $line_elements[$line_elements_counter]."\t";
		}
		elsif ($line_elements_counter > 7) { #Save GT:PL: and sample(s) GT Call fields and add to proper line last
		    
		    if ($line_elements_counter == (scalar(@line_elements) - 1)) {

			$sample_id_info .= $line_elements[$line_elements_counter];
		    }
		    else {

			$sample_id_info .= $line_elements[$line_elements_counter]."\t";
		    }
		}
	    }

	    my @key_values = split(/;/, $line_elements[7]); #Split INFO field to key=value items

	    for my $element (@key_values) {
		
		my @keys = split("=", $element);  #Key = 0 and value = 1
		
		if ($keys[0]=~/AN_(\w+)/) {  #"AN" key
		    
		    $af_key{AN}{$1} = $keys[1]; #Save value for each population
		}
		if ($keys[0]=~/AC_(\w+)/) {  #"AC" key
		    
		    $af_key{AC}{$1} = $keys[1];  #Save value for each population
		}
	    }
	    my @afs;

	    foreach my $key (keys %{$af_key{AN}}) {  #To generate AF keys
	
		if( ($af_key{AC}{$key}) && ($af_key{AN}{$key} != 0)) {

		    $af_key{AF}{$key} = $af_key{AC}{$key} / $af_key{AN}{$key};
		    push(@afs, "AF_".$key."=".$af_key{AF}{$key});
		}
	    }

	    $variant_line .= join(";", @afs, $line_elements[7]);

	    if (defined($sample_id_info)) {

		$variant_line .= "\t".$sample_id_info;
	    }
	    say STDOUT $variant_line;
	}
    }
    close()
}
