package Program::Qc::Fastqc;

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
    our @EXPORT_OK = qw(fastqc);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case

sub fastqc {

##fastqc

##Function : Perl wrapper for writing fastqc recipe to already open $FILEHANDLE or return commands array. Based on cp 0.11.5
##Returns  : "@commands"
##Arguments: $FILEHANDLE, $infile_path, $outdirectory_path, $extract, $quiet
##         : $FILEHANDLE        => Filehandle to write to
##         : $infile_path       => Infile path
##         : $outdirectory_path => Outdirectory path
##         : $extract           => If set then the zipped output file will be uncompressed in the same directory after it has been created
##         : $quiet             => Supress all progress messages on stdout and only report errors

    my ($arg_href) = @_;

    ## Default(s)
    my $extract;
    my $quiet;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $outdirectory_path;

    my $tmpl = { 
	FILEHANDLE => { store => \$FILEHANDLE},
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	outdirectory_path => { strict_type => 1, store => \$outdirectory_path},
	extract => { default => 1,
		     strict_type => 1, store => \$extract},
	quiet => { default => 0,
		   strict_type => 1, store => \$quiet},
    };
    
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Fastqc
    my @commands = qw(fastqc);  #Stores commands depending on input parameters

    if($quiet) {

	push(@commands, "--quiet");  #Supress all progress messages on stdout and only report errors
    }
    if($extract) {

	push(@commands, "--extract");  #Zipped output file will be uncompressed in the same directory after it has been created
    }
    if ($outdirectory_path) {

	push(@commands, "--outdir ".$outdirectory_path);  #Create all output files in the specified output directory
    }

    ## Infile
    push(@commands, $infile_path);

    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


1;
