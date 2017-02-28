package Program::Alignment::Sambamba;

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
    our @EXPORT_OK = qw(view index);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub view {

##view

##Function : Perl wrapper for writing sambamba view recipe to $FILEHANDLE. Based on sambamba 0.6.5
##Returns  : ""
##Arguments: $regions_ref, $FILEHANDLE, $infile_path, $outfile_path, stderrfile_path, $with_header, $show_progress, $output_format
##         : $regions_ref     => The regions to process {REF}
##         : $FILEHANDLE      => Sbatch filehandle to write to
##         : $infile_path     => Infile path
##         : $outfile_path    => Outfile path
##         : $stderrfile_path => Stderrfile path
##         : $with_header     => Include header
##         : $show_progress   => Show progress
##         : $output_format   => Output format

    my ($arg_href) = @_;

    ## Default(s)
    my $with_header;
    my $show_progress;
    my $output_format;

    ## Flatten argument(s)
    my $regions_ref;
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    
    my $tmpl = {
	regions_ref => { default => [], strict_type => 1, store => \$regions_ref},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	outfile_path => { strict_type => 1, store => \$outfile_path},
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	with_header => { default => 0,
			 allow => [0, 1],
			 strict_type => 1, store => \$with_header},
	show_progress => { default => 0,
			   allow => [0, 1],
			   strict_type => 1, store => \$show_progress},
	output_format => { default => "bam",
			   allow => ["sam", "bam", "cram", "json"],
			   strict_type => 1, store => \$output_format},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    print $FILEHANDLE "sambamba view ";  #Command

    if ($with_header) {  #Include header

	print $FILEHANDLE "--with-header ";
    }
    if ($output_format) {

	print $FILEHANDLE "--format ".$output_format." ";  #Output format
    }
    if ($show_progress) {

	print $FILEHANDLE "--show-progress ";  #Show progressbar in STDERR (works only for BAM files with no regions specified)
    }
    if ($outfile_path) {
	
	print $FILEHANDLE "--output-filename=".$outfile_path." ";  #Specify output filename
    }
    print $FILEHANDLE $infile_path." ";  #InFile

    if(@$regions_ref) {  #Limit output to regions

	print $FILEHANDLE join(" ", @{ $regions_ref })." ";
    }
    if ($stderrfile_path) {

	print $FILEHANDLE "2> ".$stderrfile_path;  #Redirect stderr output to program specific stderr file
    }
}


sub index {

##index

##Function : Perl wrapper for writing sambamba index recipe to $FILEHANDLE. Based on sambamba 0.6.5
##Returns  : ""
##Arguments: $FILEHANDLE, $infile_path, stderrfile_path, $show_progress
##         : $FILEHANDLE      => Sbatch filehandle to write to
##         : $infile_path     => Infile path
##         : $stderrfile_path => Stderrfile path
##         : $show_progress   => Show progress

    my ($arg_href) = @_;

    ## Default(s)
    my $show_progress;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $stderrfile_path;
    
    my $tmpl = {
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	show_progress => { default => 0,
			   allow => [0, 1],
			   strict_type => 1, store => \$show_progress},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    print $FILEHANDLE "sambamba index ";  #Command

    if ($show_progress) {

	print $FILEHANDLE "--show-progress ";  #Show progressbar in STDERR (works only for BAM files with no regions specified)
    }
    print $FILEHANDLE $infile_path." ";  #InFile

    if ($stderrfile_path) {

	print $FILEHANDLE "2> ".$stderrfile_path;  #Redirect stderr output to program specific stderr file
    }
}


1;
