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
    our @EXPORT_OK = qw(view index sort markdup flagstat depth);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub view {

##view

##Function : Perl wrapper for writing sambamba view recipe to $FILEHANDLE. Based on sambamba 0.6.5
##Returns  : "@commands"
##Arguments: $regions_ref, $FILEHANDLE, $infile_path, $outfile_path, $stderrfile_path, $with_header, $show_progress, $output_format, $referencefile_path
##         : $regions_ref        => The regions to process {REF}
##         : $FILEHANDLE         => Sbatch filehandle to write to
##         : $infile_path        => Infile path
##         : $outfile_path       => Outfile path
##         : $stderrfile_path    => Stderrfile path
##         : $referencefile_path => Reference for writing CRAM
##         : $with_header        => Include header
##         : $show_progress      => Show progress
##         : $output_format      => Output format

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
    my $referencefile_path;
    
    my $tmpl = {
	regions_ref => { default => [], strict_type => 1, store => \$regions_ref },
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path },
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	referencefile_path => { strict_type => 1, store => \$referencefile_path },
	with_header => { default => 0,
			 allow => [0, 1],
			 strict_type => 1, store => \$with_header },
	show_progress => { default => 0,
			   allow => [0, 1],
			   strict_type => 1, store => \$show_progress },
	output_format => { default => "bam",
			   allow => ["sam", "bam", "cram", "json"],
			   strict_type => 1, store => \$output_format },
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
    if ($referencefile_path) {

	push(@commands, "--ref-filename=".$referencefile_path);  #Reference for writing CRAM
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
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	show_progress => { default => 0,
			   allow => [0, 1],
			   strict_type => 1, store => \$show_progress },
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
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path },
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	memory_limit => { allow => qr/^\d+G$/,
			  strict_type => 1, store => \$memory_limit },
	temp_directory => { strict_type => 1, store => \$temp_directory },
	show_progress => { default => 0,
			   allow => [0, 1],
			   strict_type => 1, store => \$show_progress },
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


sub markdup {

##markdup

##Function : Perl wrapper for writing sambamba markdup recipe to $FILEHANDLE. Based on sambamba 0.6.5
##Returns  : "@commands"
##Arguments: $infile_path, $outfile_path, $stderrfile_path, $FILEHANDLE, temp_directory, $show_progress, $hash_table_size, $overflow_list_size, $io_buffer_size
##         : $infile_path        => Infile path
##         : $outfile_path       => Outfile path
##         : $stderrfile_path    => Stderrfile path
##         : $FILEHANDLE         => Sbatch filehandle to write to
##         : $temp_directory     => Specify directory for temporary files
##         : $show_progress      => Show progress
##         : $hash_table_size    => Size of hash table for finding read pairs
##         : $overflow_list_size => Size of the overflow list where reads, thrown from the hash table, get a second chance to meet their pairs
##         : $io_buffer_size      => Two buffers of BUFFER_SIZE *megabytes* each are used for reading and writing BAM during the second pass

    my ($arg_href) = @_;

    ## Default(s)
    my $show_progress;

    ## Flatten argument(s)
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $temp_directory;
    my $hash_table_size;
    my $overflow_list_size;
    my $io_buffer_size;

    my $tmpl = {
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path },
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
	temp_directory => { strict_type => 1, store => \$temp_directory },
	hash_table_size => { allow => qr/^\d+$/,
			     strict_type => 1, store => \$hash_table_size },
	overflow_list_size => { allow => qr/^\d+$/,
				strict_type => 1, store => \$overflow_list_size },
	io_buffer_size => { allow => qr/^\d+$/,
			    strict_type => 1, store => \$io_buffer_size },
	show_progress => { default => 0,
			   allow => [0, 1],
			   strict_type => 1, store => \$show_progress },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Sambamba
    my @commands = qw(sambamba markdup);  #Stores commands depending on input parameters

    if ($temp_directory) {

	push(@commands, "--tmpdir=".$temp_directory);  #Specify directory for temporary files
    }
    if ($hash_table_size) {

	push(@commands, "--hash-table-size=".$hash_table_size);  #Size of hash table for finding read pairs
    }
    if ($overflow_list_size) {

	push(@commands, "--overflow-list-size=".$overflow_list_size);  #Size of the overflow list
    }
    if ($io_buffer_size) {

	push(@commands, "--io-buffer-size=".$io_buffer_size);  #Two buffers of BUFFER_SIZE *megabytes*
    }
    if ($show_progress) {

	push(@commands, "--show-progress");  #Show progressbar in STDERR (works only for BAM files with no regions specified)
    }

    ## Infile
    push(@commands, $infile_path);

    ## Outfile
    if ($outfile_path) {
	
	push(@commands, $outfile_path);
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
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path },
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
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


sub depth {

##depth

##Function : Perl wrapper for writing sambamba depth recipe to $FILEHANDLE. Based on sambamba 0.6.5
##Returns  : "@commands"
##Arguments: $depth_cutoffs_ref, $FILEHANDLE, $infile_path, $outfile_path, $stderrfile_path, $region, $filter, $min_base_quality, $mode, $fix_mate_overlap
##         : $depth_cutoffs_ref => Multiple thresholds can be provided, for each one an extra column will be added, the percentage of bases in the region where coverage is more than this value {REF}
##         : $FILEHANDLE        => Sbatch filehandle to write to
##         : $infile_path       => Infile path
##         : $outfile_path      => Outfile path
##         : $stderrfile_path   => Stderrfile path
##         : $region            => List or regions of interest or a single region in form chr:beg-end
##         : $filter            => Set custom filter for alignments
##         : $fix_mate_overlap  => Detect overlaps of mate reads and handle them on per-base basis
##         : $min_base_quality  => don't count bases with lower base quality
##         : $mode              => Mode unit to print the statistics on

    my ($arg_href) = @_;

    ## Default(s)
    my $min_base_quality;
    my $mode;
    my $fix_mate_overlap;

    ## Flatten argument(s)
    my $depth_cutoffs_ref;
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $region;
    my $filter;
    
    my $tmpl = {
	depth_cutoffs_ref => { default => [], strict_type => 1, store => \$depth_cutoffs_ref },
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path },
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	region => { strict_type => 1, store => \$region },
	filter => { strict_type => 1, store => \$filter },
	fix_mate_overlap => { default => 0,
			      allow => [0, 1],
			      strict_type => 1, store => \$fix_mate_overlap },
	min_base_quality => { default => 0,
			      allow => qr/^\d+$/,,
			      strict_type => 1, store => \$min_base_quality },
	mode => { default => "region",
			   allow => ["base", "region", "window"],
			   strict_type => 1, store => \$mode },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Sambamba
    my @commands = qw(sambamba depth);  #Stores commands depending on input parameters

    if ($mode) {

	push(@commands, $mode);
    }
    if($region) {  #Limit output to regions

	push(@commands, "--regions ".$region);
    }
    if(@$depth_cutoffs_ref) {

	push(@commands, "--cov-threshold ".join(" --cov-threshold ", @{ $depth_cutoffs_ref }));
    }
    if ($min_base_quality) {

	push(@commands, "--min-base-quality ".$min_base_quality);
    }
    if ($fix_mate_overlap) {

	push(@commands, "--fix-mate-overlaps");
    }
    if ($filter) {

	push(@commands, "--filter ".$filter);
    }
    if ($outfile_path) {
	
	push(@commands, "--output-filename=".$outfile_path);  #Specify output filename
    }

    ## Infile
    push(@commands, $infile_path);

    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


1;
