#!/usr/bin/perl - w

use strict;
use warnings;

use Getopt::Long;

use vars qw($USAGE);

#########
##Tests##
#########
# - The file should contain a header line starting with "#".
# - The number of headers should match the number of field elements for each entry.
# - Do not use whitespace in headers.
# - Do not use ";" in file.
# - Separate elements in fields with ",". Do not use ", ".
# - No whitespace in the beginning or end within fields.
# - No entries should be duplicated within database.
# - Length of gene coordinates should be greater than 0
# - Only digits in gene coordinate entries

BEGIN {
    $USAGE =
    q?
    dbParser.pl dbfile.txt
    -s/--separator The db fields separator (defaults to "\t")
    -c/--contigColumn Contig column number (defaults to "0")
    -fs/--featureStartColumn Feature start coordinate column number (defaults to "1")
    -fe/--featureEndColumn Feature end coordinate column number (defaults to "2")
    -h/--help Display this help message
    -v/--version Display version
    ?;    
}

my ($separator, $contigColumn, $featureStartColumn, $featureEndColumn) = ("\t", 0, 1, 2);
my %validateMetric;
my ($help, $version) = (0, 0);

my $dbParserVersion = "0.0.1";

###User Options
GetOptions('s|separator:s' => \$separator,
    'c|contigColumn:s' => \$contigColumn,
    'fs|featureStartColumn:n' => \$featureStartColumn,
    'fe|featureEndColumn:n' => \$featureEndColumn,
    'h|help' => \$help,  #Display help text
    'v|version' => \$version, #Display version number
);
if (scalar(@ARGV) eq 0) {

    print STDOUT "\ndbParser.pl v".$dbParserVersion;
    print STDOUT $USAGE, "\n";
    exit;

}
else {

    foreach my $parameter (@ARGV) {

        if( ($help) || ($parameter=~/-h/) ) {

            print STDOUT $USAGE, "\n";
            exit;
        }
        if( ($version) || ($parameter=~/-v/) ) {

            print STDOUT "\ndbParser.pl v".$dbParserVersion, "\n\n";
            exit
        }
    }
}

# print messages to STDERR

{

    my $had_errors = 0; # flag set on printing to STDERR

    sub pr {
        ## pr
        ## Function: prints to STDERR and sets a flag
        ## Returns: ""
        ## Arguments: @_ will be printed to STDERR

        print STDERR @_;
        $had_errors = 1;
    }

    sub HadErrors {
        ## HadErrors
        ## Function: 
        return $had_errors;
    }

}

###
#MAIN
###

&ReadInfileDB();

###
#Sub Routines
###

sub ReadInfileDB {

##ReadInfileDB

##Function : Reads the database
##Returns  : ""
##Arguments: $infileName
##         : $infileName => The database file

    my @characterCheck = (";", ", ", "^\\s+", "\\s+\$");  

    while (<>) {

        chomp $_; # Remove newline

        if ($. eq 1) {  # Header line

            if ($_=~/^#/) {  #Header present

                my @headerElements = &SplitLine(\$_, \$separator);

                &CheckHeaderBlanks(\@headerElements, \$.);

                $validateMetric{'Header'}{'NrofElements'} = scalar(@headerElements);
            }
            else {

                pr("Could not detect any header!\n");
                pr("Aborting\n");
                exit;
            }
            next;
        }
        else {

            my @lineElements = &SplitLine(\$_, \$separator);

            &CheckNumberofElements(\@lineElements, \$., \$validateMetric{'Header'}{'NrofElements'});
            &CheckFeatureLength(\@lineElements, \$., \$featureStartColumn, \$featureEndColumn);
            &CheckElementChar(\@lineElements, \$., \@characterCheck);

            if (@lineElements[$contigColumn, $featureStartColumn, $featureEndColumn]) {  #Check that we have chromosomal and feature coordinates

                my $queryKey = @lineElements[$contigColumn, $featureStartColumn, $featureEndColumn];  #Create key for comparison
                &CheckDuplicatedEntries(\$queryKey, \$validateMetric{'line'}{$queryKey}, \$.)

            }	    
        }
    }
    my $exit_code = &HadErrors();
    close();
    exit $exit_code;
}

sub CheckElementChar {

##CheckElementChar

##Function : Check that the elements does not contain char within the elements
##Returns  : ""
##Arguments: $arrayRef, $lineCounterRef, $characterArrayRef
##         : $arrayRef          => The line elements array {REF}
##         : $lineCounterRef    => The line number {REF}
##         : $characterArrayRef => The characters to detect {REF}

    my $arrayRef = $_[0];
    my $lineCounterRef = $_[1];
    my $characterArrayRef = $_[2];

    foreach my $element (@{$arrayRef}) {

        foreach my $character (@{$characterArrayRef}) {

            if ($element =~/$character/) {

                pr("LINE: ".$$lineCounterRef." Element contains '".$character."': ".$element, "\n");
            }
        }
    }
}

sub CheckFeatureLength {

##CheckFeatureLength

##Function : Check that the genomic coordinates only contain digits and has a length gt zero
##Returns  : ""
##Arguments: $arrayRef, $lineCounterRef, $featureStartRef, $featureEndRef
##         : $arrayRef        => The line elements array {REF}
##         : $lineCounterRef  => The line number {REF}
##         : $featureStartRef => The feature start coordinate
##         : $featureEndRef => The feature end coordinate

    my $arrayRef = $_[0];
    my $lineCounterRef = $_[1];
    my $featureStartRef = $_[2];
    my $featureEndRef = $_[3];

    my $featureStartisDigit = 0;  #1 equals valid entry
    my $featureEndisDigit = 0;  #1 equals valid entry

    if (defined(@{$arrayRef}[$$featureStartRef])) {

        $featureStartisDigit = &CheckCoordinates($arrayRef, $lineCounterRef, $featureStartRef);
    }
    if (defined(@{$arrayRef}[$$featureEndRef])) {

        $featureEndisDigit = &CheckCoordinates($arrayRef, $lineCounterRef, $featureEndRef);
    }
    if ( ($featureStartisDigit == 1) && ($featureEndisDigit == 1) ) {

        my $featureLength = @{$arrayRef}[$$featureEndRef] - @{$arrayRef}[$$featureStartRef];

        if ($featureLength <= 0) {

            pr("Line: ".$$lineCounterRef." Feature length is: ".$featureLength."\n");
        }
    }
}

sub CheckCoordinates {

##CheckCoordinates

##Function : Check that the genomic coordinates only contain digits
##Returns  : "0|1"
##Arguments: $arrayRef, $lineCounterRef, $coordinateRef
##         : $arrayRef       => The line elements array {REF}
##         : $lineCounterRef => The line number {REF}
##         : $coordinateRef  => The genomic cordinate to check {REF}

    my $arrayRef = $_[0];
    my $lineCounterRef = $_[1];
    my $coordinateRef = $_[2];

    if (defined(@{$arrayRef}[$$coordinateRef])) {

        unless (@{$arrayRef}[$$coordinateRef] =~/^\d+$/) {

            pr("Line: ".$$lineCounterRef." Gene coordinate (column ".$$coordinateRef.") contains invalid characters at: ".@{$arrayRef}[$$coordinateRef]."\n");
            return 0;
        }
    }
    return 1;
}

sub CheckDuplicatedEntries {

##CheckDuplicatedEntries

##Function : Checks if there are any duplicated entries by comparing col 1-3 (chr,start,stop)
##Returns  : ""
##Arguments: $queryKeyRef, $keyRef, $separatorRef
##         : $queryKeyRef => The key to be used as query {REF}
##         : $keyRef      => The key that is stored {REF}
##         : $lineCounterRef => The line number {REF}

    my $queryKeyRef = $_[0];
    my $keyRef = $_[1];
    my $lineCounterRef = $_[2];

    if (defined($$keyRef)) {  #Key already present, hence duplicated entry

        pr("LINE: ".$$lineCounterRef." Contains a duplicated entry:", "\n");
        pr($_, "\n");
    }
    else {  #No duplicated entry

        $validateMetric{'line'}{$$queryKeyRef} = $$queryKeyRef;
    }
}

sub CheckNumberofElements {

##CheckNumberofElements

##Function : Check that the header elements does not contain whitespace within the elements
##Returns  : ""
##Arguments: $arrayRef, $lineCounterRef, $headerElementCountRef
##         : $arrayRef              => The line elements array {REF}
##         : $lineCounterRef        => The line number {REF}
##         : $headerElementCountRef => The number of header elements {REF}

    my $arrayRef = $_[0];
    my $lineCounterRef = $_[1];
    my $headerElementCountRef = $_[2];


    if (scalar(@{$arrayRef}) ne $$headerElementCountRef) {

        pr("LINE: ".$$lineCounterRef." Does not contain as many fields as header: ".scalar(@{$arrayRef}). " vs ".$$headerElementCountRef, "\n");
    }


}

sub SplitLine {

##SplitLine

##Function : Split the db line depending on the separator
##Returns  : "array"
##Arguments: $headerLineRef, $separatorRef
##         : $headerLineRef => The header line {REF}
##         : $separatorRef => The separator used to split the fields {REF}

    my $headerLineRef = $_[0];
    my $separatorRef = $_[1];

    my @lineElements = split(/$$separatorRef/, $$headerLineRef);
    return @lineElements;
}

sub CheckHeaderBlanks {

##CheckHeaderBlanks

##Function : Check that the header elements does not contain whitespace within the elements
##Returns  : ""
##Arguments: $arrayRef, $lineCounterRef
##         : $arrayRef => The line elements array {REF}
##         : $lineCounterRef => The line number {REF}

    my $arrayRef = $_[0];
    my $lineCounterRef = $_[1];

    foreach my $element (@{$arrayRef}) {

        if ($element =~/\s/) {

            pr("LINE: ".$$lineCounterRef." Element contains whitespace: ".$element, "\n");
        }
    }
}

