package Program::Alignment::Picardtools;

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
    our @EXPORT_OK = qw(mergesamfiles markduplicates);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub mergesamfiles {

##mergesamfiles

##Function : Perl wrapper for writing picardtools mergesamfiles recipe to $FILEHANDLE. Based on picardtools 2.5.0.
##Returns  : "@commands"
##Arguments: $infile_paths_ref, $outfile_path, $stderrfile_path, $FILEHANDLE, $create_index, $threading, $regionsfile_path
##         : $infile_paths_ref => Infile paths {REF}
##         : $outfile_path     => Outfile path
##         : $stderrfile_path  => Stderrfile path
##         : $FILEHANDLE       => Sbatch filehandle to write to
##         : $create_index     => Create index
##         : $threading        => Create a background thread to encode, compress and write to disk the output file.
##         : $regionsfile_path => The regions to process {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_paths_ref;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $create_index;
    my $threading;
    my $regionsfile_path;
    
    my $tmpl = {
	infile_paths_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$infile_paths_ref},
	outfile_path => { required => 1, defined => 1, strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	create_index => { allow => ["true", "false"],
			   strict_type => 1, store => \$threading },
	threading => { allow => ["true", "false"],
		       strict_type => 1, store => \$threading },
	regionsfile_path => { strict_type => 1, store => \$regionsfile_path },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Picardtools mergesamfiles
    my @commands = qw(MergeSamFiles);  #Stores commands depending on input parameters

    if ($threading) {

	push(@commands, "USE_THREADING=".$threading);  #Create a background thread to encode, compress and write to disk the output file
    }
    if ($create_index) {

	push(@commands, "CREATE_INDEX=".$create_index);
    }

    ## Infile
    push(@commands, "INPUT=".join(" INPUT=", @$infile_paths_ref));

    ## Output
    push(@commands, "OUTPUT=".$outfile_path);  #Specify output filename

    if($regionsfile_path) {  #Limit output to regions

	push(@commands, "INTERVALS=".$regionsfile_path);
    }
    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


sub markduplicates {

##markduplicates

##Function : Perl wrapper for writing picardtools markduplicates recipe to $FILEHANDLE. Based on picardtools 2.5.0.
##Returns  : "@commands"
##Arguments: $infile_paths_ref, $outfile_path, $metrics_file, $stderrfile_path, $FILEHANDLE, $create_index
##         : $infile_paths_ref => Infile paths {REF}
##         : $outfile_path     => Outfile path
##         : $metrics_file     => File to write duplication metrics to
##         : $stderrfile_path  => Stderrfile path
##         : $create_index     => Create index
##         : $FILEHANDLE       => Sbatch filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_paths_ref;
    my $outfile_path;
    my $metrics_file;
    my $stderrfile_path;
    my $create_index;
    my $FILEHANDLE;
    my $threading;
    my $regionsfile_path;
    
    my $tmpl = {
	infile_paths_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$infile_paths_ref},
	outfile_path => { required => 1, defined => 1, strict_type => 1, store => \$outfile_path },
	metrics_file => { required => 1, defined => 1, strict_type => 1, store => \$metrics_file },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	create_index => { allow => ["true", "false"],
			   strict_type => 1, store => \$threading },
	FILEHANDLE => { store => \$FILEHANDLE },
	threading => { allow => ["true", "false"],
		       strict_type => 1, store => \$threading },
	regionsfile_path => { strict_type => 1, store => \$regionsfile_path },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Picardtools markduplicates
    my @commands = qw(MarkDuplicates);  #Stores commands depending on input parameters

    if ($metrics_file) {

	push(@commands, "METRICS_FILE=".$metrics_file);  #File to write duplication metrics to
    }
    if ($create_index) {

	push(@commands, "CREATE_INDEX=".$create_index);
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
