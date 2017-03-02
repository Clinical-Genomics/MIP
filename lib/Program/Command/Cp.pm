package Program::Command::Cp;

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
    our @EXPORT_OK = qw(cp);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case

sub cp {

##cp

##Function : Perl wrapper for writing cp recipe to already open $FILEHANDLE or return commands array. Based on cp 8.4
##Returns  : "@commands"
##Arguments: $preserve_attributes_ref, $FILEHANDLE, $infile_path, $outfile_path, $stderrfile_path, $preserve, $$recursive, $force, $verbose
##         : $preserve_attributes_ref => Preserve the specified attributes (default:mode,ownership,timestamps), if possible additional attributes: context, links, xattr, all
##         : $FILEHANDLE              => Filehandle to write to
##         : $infile_path             => Infile path
##         : $outfile_path            => Outfile path
##         : $stderrfile_path         => Stderrfile path
##         : $preserve                => Same as --preserve=mode,ownership,timestamps
##         : $recursive               => Copy directories recursively
##         : $force                   => If an existing destination file cannot be opened, remove it and try again
##         : $verbose                 => Verbosity

    my ($arg_href) = @_;

    ## Default(s)
    my $preserve;
    my $recursive;
    my $force;
    my $verbose;

    ## Flatten argument(s)
    my $preserve_attributes_ref;
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;

    my $tmpl = { 
	preserve_attributes_ref => { default => [], strict_type => 1, store => \$preserve_attributes_ref},
	FILEHANDLE => { store => \$FILEHANDLE},
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	outfile_path => { required => 1, defined => 1, strict_type => 1, store => \$outfile_path},
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	recursive => { default => 0,
		       strict_type => 1, store => \$recursive},
	force => { default => 0,
		   strict_type => 1, store => \$force},
	preserve => { default => 0,
		      strict_type => 1, store => \$preserve},
	verbose => { default => 0,
		     strict_type => 1, store => \$verbose},
    };
    
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## cp
    my @commands = qw(cp);  #Stores commands depending on input parameters

    if(@$preserve_attributes_ref) {

	push(@commands, "--preserve=".join(",", @$preserve_attributes_ref));  #Preserve the specified attributes
    }
    elsif($preserve) {

	push(@commands, "-p");
    }
    if ($recursive) {

	push(@commands, "--recursive");
    }
    if ($force) {

	push(@commands, "--force");
    }
    if ($verbose) {

	push(@commands, "--verbose");  #Explain what is being done
    }
    push(@commands, $infile_path);
    push(@commands, $outfile_path);

    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}

1;
