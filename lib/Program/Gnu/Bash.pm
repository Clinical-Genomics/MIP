package Program::Gnu::Bash;

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
    our @EXPORT_OK = qw(cd);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub cd {

##cd

##Function : Perl wrapper for writing cd recipe to already open $FILEHANDLE or return commands array. Based on cd 4.0
##Returns  : "@commands"
##Arguments: $FILEHANDLE, $directory_path, $stderrfile_path
##         : $FILEHANDLE      => Filehandle to write to
##         : $directory_path  => Directory path
##         : $stderrfile_path => Stderrfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $directory_path;
    my $stderrfile_path;

    my $tmpl = { 
	FILEHANDLE => { store => \$FILEHANDLE},
	directory_path => { strict_type => 1, store => \$directory_path},
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
    };
    
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## cd
    my @commands = qw(cd);  #Stores commands depending on input parameters

    ## Indirectory
    if ($directory_path) {

	push(@commands, $directory_path);
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
