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
use Params::Check qw[check allow last_error];

## Third party module(s)
use List::Util qw(any);
use Scalar::Util::Numeric qw(isnum);

use vars qw($USAGE);

BEGIN {
    $USAGE = qq{calculate_af.pl - [vcf]
               -cak/--af_key_suffixes Key(s) to calculate AF for (AC_af_key_suffix / AN_af_key_suffix)
               -h/--help Display this help message
               -v/--version Display version};
}

my ( $infile, $version, $help ) = ("");
my @af_key_suffixes;

my $calculate_af_version = "0.0.3";

## Enables cmd "calculate_af.pl" to print usage help
if ( scalar(@ARGV) == 0 ) {

    print STDOUT $USAGE, "\n";
    exit;
}
elsif ( ( defined($ARGV) ) && ( $ARGV[0] !~ /^-/ ) )
{    #Collect potential infile - otherwise read from STDIN

    $infile = $ARGV[0];
}

GetOptions(
    'cak|af_key_suffixes:s' => \@af_key_suffixes,
    'h|help'    => sub { print STDOUT $USAGE, "\n"; exit; },    #Display help text
    'v|version' => sub { print STDOUT "\ncalculate_af.pl " . $calculate_af_version, "\n\n"; exit; }
    ,                                                           #Display version number
);

#####
#MAIN
#####

read_infile(
    {
        af_key_suffixes_ref  => \@af_key_suffixes,
        calculate_af_version => $calculate_af_version,
    }
);

###
#Sub Routines
###

sub read_infile {

##read_infile

##Function : Read infile and calculate allele frequency
##Returns  : ""
##Arguments: $af_key_suffixes_ref, $calculate_af_version
##         : $af_key_suffixes_ref => AN and AC key to calculate frequency for
##         : $calculate_af_version => Calculate af version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $calculate_af_version;
    my $af_key_suffixes_ref;

    my $tmpl = {
        af_key_suffixes_ref  => { default => [], strict_type => 1, store => \$af_key_suffixes_ref },
        calculate_af_version => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$calculate_af_version
        },
    };

    check( $tmpl, $arg_href, 1 ) or die qw[Could not parse arguments!];

    my %af_header;
    my %af_key;
    my @af_prefixes = ( "AN", "AC" );

    while (<>) {

        my $variant_line;
        my $sample_id_info;

        chomp $_;    # Remove newline

        if (m/^\s+$/) {    # Avoid blank lines
            next;
        }
        if ( $_ =~ /^##(\S+)=/ ) {    # MetaData

            my $key = $1;

            foreach my $af_prefix (@af_prefixes) {

                if ( $key =~ /$af_prefix\_(\w+)/ ) {    #Collect af prefix keys

                    my $af_key_suffix = $1;

                    if (@$af_key_suffixes_ref) {

                        if ( any { $_ eq $af_key_suffix } @$af_key_suffixes_ref )
                        {                               #If element is part of array

                            $af_header{$af_prefix}{$af_key_suffix}       = $key; #For writing header
                            $af_key{$af_prefix}{$af_key_suffix}          = undef;
                            $af_key{ $af_prefix . "_" . $af_key_suffix } = undef;
                        }
                    }
                    else {

                        $af_header{$af_prefix}{$af_key_suffix}       = $key;     #For writing header
                        $af_key{$af_prefix}{$af_key_suffix}          = undef;
                        $af_key{ $af_prefix . "_" . $af_key_suffix } = undef;
                    }
                }
            }
            say STDOUT $_;
            next;
        }
        elsif ( $_ =~ /^#CHROM/ ) {

            if (@$af_key_suffixes_ref) {

                foreach my $af_key_suffix (@$af_key_suffixes_ref) {

                    if ( !exists( $af_header{AN}{$af_key_suffix} ) ) {

                        say STDERR "Could not find af_key_suffix: " . $af_key_suffix . " in header";
                        exit 1;
                    }

                    ## Write vcf Header
                    say '##INFO=<ID=AF_'
                      . $af_key_suffix
                      . ',Number=A,Type=Float,Description="Allele frequency in the '
                      . $af_key_suffix
                      . ' populations calculated from AC and AN, in the range (0,1)">';
                }
            }
            else {

                foreach my $key ( keys %{ $af_header{AN} } ) {    #Either "AN" or "AC" keys will do

                    ## Write vcf Header
                    say '##INFO=<ID=AF_'
                      . $key
                      . ',Number=A,Type=Float,Description="Allele frequency in the '
                      . $key
                      . ' populations calculated from AC and AN, in the range (0,1)">';
                }
            }

            ## Write program cmd
            my ( $base, $script ) =
              ( `date +%Y%m%d`, `basename $0` );    #Catches current date and script name
            chomp( $base, $script );                #Remove \n;
            say "##software=<ID="
              . $script
              . ",Version="
              . $calculate_af_version
              . ",Date="
              . $base;

            say STDOUT $_;
            next;
        }
        else {

            my $calculate_af  = 0;
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

            my $af_counter = 0;
            for my $element (@key_values) {

                my @keys = split( "=", $element );    #Key = 0 and value = 1

                if ( exists( $af_key{ $keys[0] } ) ) {

                    if ( $keys[0] =~ /AN_(\w+)/ ) {    #"AN" key

                        if ( exists( $af_key{AN}{$1} ) ) {

                            $af_key{AN}{$1} = $keys[1];    #Save value for each population
                            $af_counter++;
                        }
                    }
                    elsif ( $keys[0] =~ /AC_(\w+)/ ) {     #"AC" key

                        if ( exists( $af_key{AC}{$1} ) ) {

                            $af_key{AC}{$1} = $keys[1];    #Save value for each population
                            $af_counter++;
                        }
                    }
                    if (   (@$af_key_suffixes_ref)
                        && ( $af_counter eq ( scalar(@$af_key_suffixes_ref) * 2 ) ) )
                    {

                        last;                              #Found all AC and AN
                    }
                }
            }
            my @afs;

            if (@$af_key_suffixes_ref) {                   #For selected AF keys

                foreach my $af_key_suffix (@$af_key_suffixes_ref) {

                    calculate_af(
                        {
                            af_key_href   => \%af_key,
                            afs_ref       => \@afs,
                            af_key_suffix => $af_key_suffix,
                        }
                    );

                }
            }
            else {

                foreach my $key ( keys %{ $af_key{AN} } ) {    #To generate AF keys

                    calculate_af(
                        {
                            af_key_href   => \%af_key,
                            afs_ref       => \@afs,
                            af_key_suffix => $key,
                        }
                    );
                }
            }

            $variant_line .= join( ";", @afs, $line_elements[7] );

            if ( defined($sample_id_info) ) {

                $variant_line .= "\t" . $sample_id_info;
            }
            say STDOUT $variant_line;
        }
    }
    close();
}

sub calculate_af {

##calculate_af

##Function : Calculate af for af key
##Returns  : ""
##Arguments: $af_key_href, $afs_ref, $af_key_suffix
##         : $af_key_href   => Holds allele frequency keys
##         : $afs_ref       => Calculated allele frequencies
##         : $af_key_suffix => Key to calcluate allele frequency for

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $af_key_href;
    my $afs_ref;
    my $af_key_suffix;

    my $tmpl = {
        af_key_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$af_key_href
        },
        afs_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$afs_ref
        },
        af_key_suffix =>
          { required => 1, defined => 1, strict_type => 1, store => \$af_key_suffix },
    };

    check( $tmpl, $arg_href, 1 ) or die qw[Could not parse arguments!];

    if (   ( $af_key_href->{AC}{$af_key_suffix} )
        && ( isnum( $af_key_href->{AC}{$af_key_suffix} ) ) )
    {    #For floats and integers

        if ( $af_key_href->{AN}{$af_key_suffix} != 0 ) {    #Do not divide by 0

            $af_key_href->{AF}{$af_key_suffix} =
              $af_key_href->{AC}{$af_key_suffix} / $af_key_href->{AN}{$af_key_suffix};
            push( @$afs_ref, "AF_" . $af_key_suffix . "=" . $af_key_href->{AF}{$af_key_suffix} );
        }
    }
}
