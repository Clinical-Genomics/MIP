package Program::Gnu::Software::Less;

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
    our @EXPORT_OK = qw(less);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub less {

##less

##Function : Perl wrapper for writing less recipe to already open $FILEHANDLE or return commands array. Based on less 436
##Returns  : "@commands"
##Arguments: $FILEHANDLE, $infile_path, $outfile_path, $stderrfile_path
##         : $FILEHANDLE       => Filehandle to write to
##         : $infile_path      => Infile path
##         : $outfile_path     => Outfile path
##         : $stderrfile_path  => Stderrfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;

    my $tmpl = { 
	FILEHANDLE => { store => \$FILEHANDLE},
	infile_path => { strict_type => 1, store => \$infile_path},
	outfile_path => { strict_type => 1, store => \$outfile_path},
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
    };
    
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## less
    my @commands = qw(less);  #Stores commands depending on input parameters

    ## Options

    ## Infile
    push(@commands, $infile_path);

    ## Outfile
    if ($outfile_path) {

	push(@commands, "> ".$outfile_path);
    }
    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


1;
