package Program::Interval::Picardtools;

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
    our @EXPORT_OK = qw(intervallisttools);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub intervallisttools {

##intervallisttools

##Function : Perl wrapper for writing picardtools intervallisttools recipe to $FILEHANDLE. Based on picardtools 2.5.0.
##Returns  : "@commands"
##Arguments: $infile_paths_ref, $outfile_path, $FILEHANDLE, $stderrfile_path, $padding
##         : $infile_paths_ref    => Infile paths {REF}
##         : $outfile_path        => Outfile path
##         : $FILEHANDLE          => Sbatch filehandle to write to
##         : $stderrfile_path     => Stderrfile path
##         : $padding             => The amount to pad each end of the intervals by before other operations are undertaken

    my ($arg_href) = @_;

    ## Default(s)
    my $padding;

    ## Flatten argument(s)
    my $infile_paths_ref;
    my $outfile_path;
    my $FILEHANDLE;
    my $stderrfile_path;
    
    my $tmpl = {
	infile_paths_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$infile_paths_ref},
	outfile_path => { required => 1, defined => 1, strict_type => 1, store => \$outfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	padding => { default => 0,
		     allow => qr/^\d+$/,
		     strict_type => 1, store => \$padding},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Picardtools intervallisttools
    my @commands = qw(IntervalListTools);  #Stores commands depending on input parameters

    ##Options
    if ($padding) {

	push(@commands, "PADDING=".$padding);
    }

    ## Infile
    push(@commands, "INPUT=".join(" INPUT=", @$infile_paths_ref));

    ## Output
    push(@commands, "OUTPUT=".$outfile_path);  #Specify output filename

    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


1;

