package Program::Variantcalling::Svdb;

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
    our @EXPORT_OK = qw(merge);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub merge {

##merge

##Function : Perl wrapper for writing svdb merge recipe to $FILEHANDLE or return commands array. Based on svdb 0.1.2.
##Returns  : "@commands"
##Arguments: $infile_paths_ref, $outfile_path, $FILEHANDLE, $priority
##         : $infile_paths_ref => Infile path
##         : $outfile_path     => Outfile path
##         : $FILEHANDLE       => Filehandle to write to
##         : $priority         => Priority order of structural variant calls

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_paths_ref;
    my $outfile_path;
    my $FILEHANDLE;
    my $priority;

    my $tmpl = {
	infile_paths_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$infile_paths_ref },
	outfile_path => { strict_type => 1, store => \$outfile_path },
	priority => { strict_type => 1, store => \$priority },
	FILEHANDLE => { store => \$FILEHANDLE },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## svdb
    my @commands = qw(svdb --merge);  #Stores commands depending on input parameters

    ## Options
    if ($priority) {

	push(@commands, "-priority ".$priority);  #Priority order of structural variant calls
    }

    ## Infile
    push(@commands, "-vcf ".join(" ", @{ $infile_paths_ref }));

    if ($outfile_path) {

	push(@commands, "> ".$outfile_path);  #Outfile prefix
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


1;
