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
    our @EXPORT_OK = qw(view index sort flagstat);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub view {

##view

##Function : Perl wrapper for writing sambamba view recipe to $FILEHANDLE. Based on sambamba 0.6.5
##Returns  : "@commands"
##Arguments: $regions_ref, $FILEHANDLE, $infile_path, $outfile_path, $stderrfile_path, $with_header, $show_progress, $output_format
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

    ## Sambamba
    my @commands = qw(sambamba view);  #Stores commands depending on input parameters

    if ($with_header) {  #Include header

	push(@commands, "--with-header");
    }
    if ($output_format) {

	push(@commands, "--format ".$output_format);  #Output format
    }
    if ($show_progress) {

	push(@commands, "--show-progress");  #Show progressbar in STDERR (works only for BAM files with no regions specified)
    }
    if ($outfile_path) {
	
	push(@commands, "--output-filename=".$outfile_path);  #Specify output filename
    }

    ## Infile
    push(@commands, $infile_path);

    if(@$regions_ref) {  #Limit output to regions

	push(@commands, join(" ", @{ $regions_ref }));
    }
    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


sub index {

##index

##Function : Perl wrapper for writing sambamba index recipe to $FILEHANDLE. Based on sambamba 0.6.5
##Returns  : "@commands"
##Arguments: $FILEHANDLE, $infile_path, $stderrfile_path, $show_progress
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

    ## Sambamba
    my @commands = qw(sambamba index);  #Stores commands depending on input parameters

    if ($show_progress) {

	push(@commands, "--show-progress");  #Show progressbar in STDERR (works only for BAM files with no regions specified)
    }

    ##Infile
    push(@commands, $infile_path);

    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


sub sort {

##sort

##Function : Perl wrapper for writing sambamba sort recipe to $FILEHANDLE. Based on sambamba 0.6.5
##Returns  : "@commands"
##Arguments: $FILEHANDLE, $infile_path, $outfile_path, $stderrfile_path, $show_progress, $memory_limit, $temp_directory
##         : $FILEHANDLE      => Sbatch filehandle to write to
##         : $infile_path     => Infile path
##         : $outfile_path    => Outfile path
##         : $stderrfile_path => Stderrfile path
##         : $show_progress   => Show progress
##         : $memory_limit    => Approximate total memory limit for all threads
##         : $temp_directory  => Directory for storing intermediate files; default is system directory for temporary files

    my ($arg_href) = @_;

    ## Default(s)
    my $show_progress;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $memory_limit;
    my $temp_directory;
    
    my $tmpl = {
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	outfile_path => { strict_type => 1, store => \$outfile_path},
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	memory_limit => { allow => qr/^\d+G/,
			  strict_type => 1, store => \$memory_limit},
	temp_directory => { strict_type => 1, store => \$temp_directory},
	show_progress => { default => 0,
			   allow => [0, 1],
			   strict_type => 1, store => \$show_progress},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Sambamba
    my @commands = qw(sambamba sort);  #Stores commands depending on input parameters

    if ($show_progress) {

	push(@commands, "--show-progress");  #Show progressbar in STDERR (works only for BAM files with no regions specified)
    }
    if ($memory_limit) {

	push(@commands, "--memory-limit=".$memory_limit);  #Approximate total memory limit for all threads
    }
    if ($temp_directory) {
	
	push(@commands, "--tmpdir=".$temp_directory);  #Directory for storing intermediate files
    }

    ## Outfile
    if ($outfile_path) {
	
	push(@commands, "--out=".$outfile_path);
    }

    ##Infile
    if ($infile_path) {

	push(@commands, $infile_path);
    }
    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


sub flagstat {

##flagstat

##Function : Perl wrapper for writing sambamba flagstat recipe to $FILEHANDLE. Based on sambamba 0.6.5
##Returns  : "@commands"
##Arguments: $FILEHANDLE, $infile_path, $outfile_path, $stderrfile_path, $show_progress, $memory_limit, $temp_directory
##         : $FILEHANDLE      => Sbatch filehandle to write to
##         : $infile_path     => Infile path
##         : $outfile_path    => Outfile path
##         : $stderrfile_path => Stderrfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;

    my $tmpl = {
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	outfile_path => { strict_type => 1, store => \$outfile_path},
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	FILEHANDLE => { store => \$FILEHANDLE},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Sambamba
    my @commands = qw(sambamba flagstat);  #Stores commands depending on input parameters

    ##Infile
    if ($infile_path) {

	push(@commands, $infile_path);
    }
    ## Outfile
    if ($outfile_path) {
	
	push(@commands, "> ".$outfile_path);
    }
    if ($stderrfile_path) {

	push(@commands, "2>> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


1;
