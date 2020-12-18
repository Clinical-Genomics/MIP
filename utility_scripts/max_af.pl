#!/usr/bin/env perl

use warnings;
use Modern::Perl qw{ 2018 };
use warnings qw( FATAL utf8 );
use autodie;
use v5.18;    #Require at least perl 5.18
use utf8;     #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );

use Pod::Usage;
use Pod::Text;
use Getopt::Long;

use vars qw($USAGE);

BEGIN {
    $USAGE = qq{af_max.pl - [vcf]
               -h/--help Display this help message
               -v/--version Display version};
}

my ( $infile, $version, $help ) = ("");

my $max_af_version = "0.0.2";

## Enables cmd "max_af.pl" to print usage help
if ( scalar(@ARGV) == 0 ) {

    print STDOUT $USAGE, "\n";
    exit;
}
elsif ( ( defined($ARGV) ) && ( $ARGV[0] !~ /^-/ ) )
{    #Collect potential infile - otherwise read from STDIN

    $infile = $ARGV[0];
}

GetOptions(
    'h|help' => sub { print STDOUT $USAGE, "\n"; exit; },    #Display help text
    'v|version' => sub { print STDOUT "\nmax_af.pl " . $max_af_version, "\n\n"; exit; }
    ,                                                        #Display version number
);

###
#MAIN
###

read_infile();

###
#Sub Routines
###

sub read_infile {

##read_infile

##Function : Read infile and calculate max allele frequency
##Returns  : ""
##Arguments:

    while (<>) {

        my $variant_line;
        my $sample_id_info;

        chomp $_;    # Remove newline

        if (m/^\s+$/) {    # Avoid blank lines
            next;
        }
        if ( $_ =~ /^##(\S+)=/ ) {    # MetaData

            say STDOUT $_;
            next;
        }
        elsif ( $_ =~ /^#CHROM/ ) {

            ## Write vcf Header
            say
'##INFO=<ID=MAX_AF,Number=A,Type=Float,Description="Maximum Allele Frequency for variant, across populations ">';
            say STDOUT $_;
            next;
        }
        else {

            my $max_af        = 0;
            my @line_elements = split( "\t", $_ );

            for (
                my $line_elements_counter = 0 ;
                $line_elements_counter < scalar(@line_elements) ;
                $line_elements_counter++
              )
            {    #Add until INFO field

                if ( $line_elements_counter < 7 ) {    #Save fields until INFO field

                    $variant_line .= $line_elements[$line_elements_counter] . "\t";
                }
                elsif ( $line_elements_counter > 7 )
                {    #Save GT:PL: and sample(s) GT Call fields and add to proper line last

                    if ( $line_elements_counter == ( scalar(@line_elements) - 1 ) ) {

                        $sample_id_info .= $line_elements[$line_elements_counter];
                    }
                    else {

                        $sample_id_info .= $line_elements[$line_elements_counter] . "\t";
                    }
                }
            }

            my @key_values =
              split( /;/, $line_elements[7] );    #Split INFO field to key=value items

            for my $element (@key_values) {

                my @keys = split( "=", $element );    #Key = 0 and value = 1

                if ( $keys[0] =~ /_AF$/ ) {           #Allele frequency key

                    if ( $keys[1] > $max_af ) {

                        $max_af = $keys[1];
                    }
                }
                elsif ( $keys[0] =~ /^AF_/ ) {        #Allele frequency key

                    if ( $keys[1] > $max_af ) {

                        $max_af = $keys[1];
                    }
                }
            }
            if ( $max_af != 0 ) {

                $variant_line .= join( ";", "MAX_AF=" . $max_af, $line_elements[7] );
            }
            else {

                $variant_line .= join( "", $line_elements[7] );
            }
            if ( defined($sample_id_info) ) {

                $variant_line .= "\t" . $sample_id_info;
            }
            say $variant_line;
        }
    }
    close();
}
