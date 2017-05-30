package Language::Java;

use strict;
use warnings;
use warnings qw( FATAL utf8 );
use utf8;  #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );

BEGIN {
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.00;

    # Inherit from Exporter to export functions and variables
    our @ISA = qw(Exporter);

    # Functions and variables which are exported by default
    our @EXPORT = qw();

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(core);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case

sub core {

##core

##Function : Perl wrapper for writing GATK realignertargetcreator recipe to $FILEHANDLE. Based on java openjdk version "1.8.0_92".
##Returns  : ""
##Arguments: $memory_allocation, $FILEHANDLE, $temp_directory, $java_jar, $java_use_large_pages
##         : $memory_allocation    => Memory allocation for java
##         : $FILEHANDLE           => Filehandle to write to
##         : $temp_directory       => Redirect tmp files to java temp {Optional}
##         : $java_jar             => The JAR
##         : $java_use_large_pages => Use java large pages

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $memory_allocation;
    my $FILEHANDLE;
    my $java_use_large_pages;
    my $java_jar;
    my $temp_directory;

    my $tmpl = {
	memory_allocation => { required => 1, defined => 1, strict_type => 1, store => \$memory_allocation},
	FILEHANDLE => { store => \$FILEHANDLE},
	temp_directory => { strict_type => 1, store => \$temp_directory},
	java_jar => { strict_type => 1, store => \$java_jar},
	java_use_large_pages => { default => 0,
				  allow => [0, 1],
				  strict_type => 1, store => \$java_use_large_pages},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Java core
    my @commands = qw(java);  #Stores commands depending on input parameters

    push(@commands, "-".$memory_allocation);

    if ($java_use_large_pages) {

	push(@commands, "-XX:-UseLargePages");  #UseLargePages for requiring large memory pages (cross-platform flag)
    }
    if (defined($temp_directory)) {

	push(@commands, "-Djava.io.tmpdir=".$temp_directory);  #Temporary directory
    }
    if (defined($java_jar)) {

	push(@commands, "-jar ".$java_jar);
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}
