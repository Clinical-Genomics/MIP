package Program::Variantcalling::Rhocall;

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
    our @EXPORT_OK = qw(aggregate annotate);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub aggregate {

##aggregate

##Function : Perl wrapper for writing rhocall aggregate recipe to $FILEHANDLE or return commands array. Based on rhocall 0.3.
##Returns  : "@commands"
##Arguments: $infile_path, $outfile_path, $stderrfile_path, $stdoutfile_path, $FILEHANDLE
##         : $infile_path     => Infile path to read from
##         : $outfile_path    => Outfile path to write to
##         : $stderrfile_path => Stderr file path to write to
##         : $stdoutfile_path => Stdoutfile file path to write to
##         : $FILEHANDLE      => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stdoutfile_path;
    my $FILEHANDLE;

    my $tmpl = {
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path },
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## rhocall
    my @commands = qw(rhocall aggregate);  #Stores commands depending on input parameters

    ## Options
    if ($outfile_path) {

	push(@commands, "--output ".$outfile_path);  #Specify output filename
    }

    ## Infile
    if ($infile_path) {

	push(@commands, $infile_path);
    }
    if ($stdoutfile_path) {

	push(@commands, "1> ".$stdoutfile_path);  #Redirect stdout to program specific stderr file
    }
    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


sub annotate {

##annotate

##Function : Perl wrapper for writing rhocall annotate recipe to $FILEHANDLE or return commands array. Based on rhocall 0.3.
##Returns  : "@commands"
##Arguments: $infile_path, $outfile_path, $stderrfile_path, $stdoutfile_path, $FILEHANDLE, $bedfile_path
##         : $infile_path     => Infile path to read from
##         : $outfile_path    => Outfile path to write to
##         : $stderrfile_path => Stderr file path to write to
##         : $stdoutfile_path => Stdoutfile file path to write to
##         : $FILEHANDLE      => Filehandle to write to
##         : $bedfile_path    => BED file with AZ windows

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stdoutfile_path;
    my $FILEHANDLE;
    my $bedfile_path;

    my $tmpl = {
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path },
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	bedfile_path => { strict_type => 1, store => \$bedfile_path },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## rhocall
    my @commands = qw(rhocall annotate);  #Stores commands depending on input parameters

    ## Options
    if ($bedfile_path) {

	push(@commands, "-b ".$bedfile_path);
    }
    if ($outfile_path) {

	push(@commands, "--output ".$outfile_path);  #Specify output filename
    }

    ## Infile
    if ($infile_path) {

	push(@commands, $infile_path);
    }
    if ($stdoutfile_path) {

	push(@commands, "1> ".$stdoutfile_path);  #Redirect stdout to program specific stderr file
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
