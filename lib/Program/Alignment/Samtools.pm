package Program::Alignment::Samtools;

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
    our @EXPORT_OK = qw(view index stats mpileup faidx);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub view {

##view

##Function : Perl wrapper for writing samtools view recipe to $FILEHANDLE. Based on samtools 1.3.1 (using htslib 1.3.1).
##Returns  : "@commands"
##Arguments: $regions_ref, $infile_path, $outfile_path, $stderrfile_path, $FILEHANDLE, $thread_number, $with_header, $output_format, $auto_detect_input_format, $uncompressed_bam_output
##         : $regions_ref              => The regions to process {REF}
##         : $infile_path              => Infile path
##         : $outfile_path             => Outfile path
##         : $stderrfile_path          => Stderrfile path
##         : $FILEHANDLE               => Sbatch filehandle to write to
##         : $thread_number            => Number of BAM/CRAM compression threads
##         : $with_header              => Include header
##         : $output_format            => Output format
##         : $auto_detect_input_format => Ignored (input format is auto-detected)
##         : $uncompressed_bam_output  => Uncompressed bam output

    my ($arg_href) = @_;

    ## Default(s)
    my $with_header;
    my $output_format;
    my $auto_detect_input_format;
    my $uncompressed_bam_output;

    ## Flatten argument(s)
    my $regions_ref;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $thread_number;
    
    my $tmpl = {
	regions_ref => { default => [], strict_type => 1, store => \$regions_ref },
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path },
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	thread_number => { allow => qr/^\d+$/,
			   strict_type => 1, store => \$thread_number },
	FILEHANDLE => { store => \$FILEHANDLE },
	with_header => { default => 0,
			 allow => [0, 1],
			 strict_type => 1, store => \$with_header },
	output_format => { default => "bam",
			   allow => ["sam", "bam", "cram", "json"],
			   strict_type => 1, store => \$output_format },
	auto_detect_input_format => { default => 0,
				      allow => [0, 1],
				      strict_type => 1, store => \$auto_detect_input_format },
	uncompressed_bam_output => { default => 0,
				     allow => [undef, 0, 1],
				     strict_type => 1, store => \$uncompressed_bam_output },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Samtools view
    my @commands = qw(samtools view);  #Stores commands depending on input parameters

    ## Options
    if ($thread_number) {

	push(@commands, "--threads ".$thread_number);  #Number of threads
    }
    if ($with_header) {  #Include header

	push(@commands, "-h");
    }
    if ($output_format) {

	push(@commands, "--output-fmt ".uc($output_format));  #Output format
    }
    if ($auto_detect_input_format) {

	push(@commands, "-S");
    }
    if ($outfile_path) {
	
	push(@commands, "-o".$outfile_path);  #Specify output filename
    }
    if ($uncompressed_bam_output) {

	push(@commands, "-u");
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

##Function : Perl wrapper for writing samtools index recipe to $FILEHANDLE. Based on samtools 1.3.1 (using htslib 1.3.1).
##Returns  : "@commands"
##Arguments: $infile_path, $stderrfile_path, $FILEHANDLE, $bai_format
##         : $infile_path     => Infile path
##         : $stderrfile_path => Stderrfile path
##         : $FILEHANDLE      => Sbatch filehandle to write to
##         : $bai_format      => Generate BAI-format index for BAM files

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $bai_format;
    
    my $tmpl = {
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path },
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	bai_format => { strict_type => 1, store => \$bai_format },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Samtools index
    my @commands = qw(samtools index);  #Stores commands depending on input parameters

    ## Options
    if ($bai_format) {

	push(@commands, "-b");  #Generate BAI-format index for BAM files
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


sub stats {

##stats

##Function : Perl wrapper for writing samtools stats recipe to $FILEHANDLE. Based on samtools 1.3.1 (using htslib 1.3.1).
##Returns  : "@commands"
##Arguments: $regions_ref, $infile_path, $outfile_path, $stderrfile_path, $FILEHANDLE, $auto_detect_input_format
##         : $regions_ref              => The regions to process {REF}
##         : $infile_path              => Infile path
##         : $outfile_path             => Outfile path
##         : $stderrfile_path          => Stderrfile path
##         : $FILEHANDLE               => Sbatch filehandle to write to
##         : $auto_detect_input_format => Ignored (input format is auto-detected)

    my ($arg_href) = @_;

    ## Default(s)
    my $auto_detect_input_format;

    ## Flatten argument(s)
    my $regions_ref;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    
    my $tmpl = {
	regions_ref => { default => [], strict_type => 1, store => \$regions_ref },
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path },
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	auto_detect_input_format => { default => 0,
				      allow => [0, 1],
				      strict_type => 1, store => \$auto_detect_input_format },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Samtools stats
    my @commands = qw(samtools stats);  #Stores commands depending on input parameters

    if ($auto_detect_input_format) {

	push(@commands, "-s");
    }

    ## Infile
    push(@commands, $infile_path);

    if(@$regions_ref) {  #Limit output to regions

	push(@commands, join(" ", @{ $regions_ref }));
    }
    if ($outfile_path) {
	
	push(@commands, "> ".$outfile_path);  #Specify output filename
    }
    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


sub mpileup {

##mpileup

##Function : Perl wrapper for writing samtools mpileup recipe to $FILEHANDLE. Based on samtools 1.3.1 (using htslib 1.3.1).
##Returns  : "@commands"
##Arguments: $infile_paths_ref, $output_tags_ref, $outfile_path, $referencefile_path, $stderrfile_path, $FILEHANDLE, $output_bcf, $adjust_mq
##         : $infile_paths_ref                 => Infile paths {REF}
##         : $output_tags_ref                  => Optional tags to output {REF}
##         : $outfile_path                     => Outfile path
##         : $referencefile_path               => Reference sequence file
##         : $stderrfile_path                  => Stderrfile path
##         : $FILEHANDLE                       => Sbatch filehandle to write to
##         : $region                           => The regions to process {REF}
##         : $output_bcf                       => Generate genotype likelihoods in BCF format
##         : $per_sample_increased_sensitivity => Apply -m and -F per-sample for increased sensitivity
##         : $adjust_mq                        => Adjust mapping quality

    my ($arg_href) = @_;

    ## Default(s)
    my $per_sample_increased_sensitivity;
    my $adjust_mq;

    ## Flatten argument(s)
    my $infile_paths_ref;
    my $output_tags_ref;
    my $outfile_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $region;
    my $output_bcf;
    
    my $tmpl = {
	infile_paths_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$infile_paths_ref},
	output_tags_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$output_tags_ref},
	outfile_path => { strict_type => 1, store => \$outfile_path },
	referencefile_path => { required => 1, defined => 1, strict_type => 1, store => \$referencefile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	region => { strict_type => 1, store => \$region },
	output_bcf => { strict_type => 1, store => \$output_bcf },
	per_sample_increased_sensitivity => { default => 0,
					      allow => [undef, 0, 1],
					      strict_type => 1, store => \$per_sample_increased_sensitivity },
	adjust_mq => { default => 50,
		       allow => qr/^\d+$/,
		       strict_type => 1, store => \$adjust_mq },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Samtools mpileup
    my @commands = qw(samtools mpileup);  #Stores commands depending on input parameters

    ## Options
    push(@commands, "--adjust-MQ ".$adjust_mq);

    if($per_sample_increased_sensitivity) {
	
	push(@commands, "--per-sample-mF");
    }
    if(@$output_tags_ref) {
	
	push(@commands, "--output-tags ".join(",", @$output_tags_ref));
    }
    if($region) {  #Limit output to region
	
	push(@commands, "--region ".$region);
    }
    if ($referencefile_path) {

	push(@commands, "--fasta-ref ".$referencefile_path);  #Reference sequence file
    }
    if($output_bcf) {
	
	push(@commands, "--BCF");
    }
    if ($outfile_path) {
	
	push(@commands, "--output ".$outfile_path);  #Specify output filename
    }

    ## Infile
    push(@commands, join(" ", @$infile_paths_ref));

    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


sub faidx {

##faidx

##Function : Perl wrapper for writing samtools faidx recipe to $FILEHANDLE. Based on samtools 1.3.1 (using htslib 1.3.1).
##Returns  : "@commands"
##Arguments: $regions_ref, $infile_path, $outfile_path, $stderrfile_path, $FILEHANDLE
##         : $regions_ref     => The regions to process {REF}
##         : $infile_path     => Infile path
##         : $outfile_path    => Outfile path
##         : $stderrfile_path => Stderrfile path
##         : $FILEHANDLE      => Sbatch filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $regions_ref;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    
    my $tmpl = {
	regions_ref => { default => [], strict_type => 1, store => \$regions_ref },
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path },
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Samtools faidx
    my @commands = qw(samtools faidx);  #Stores commands depending on input parameters

    ## Infile
    push(@commands, $infile_path);

    if(@$regions_ref) {  #Limit output to regions

	push(@commands, join(" ", @{ $regions_ref }));
    }
    if ($outfile_path) {
	
	push(@commands, "> ".$outfile_path);  #Specify output filename
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
