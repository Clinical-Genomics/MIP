package Program::Alignment::Bwa;

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
    our @EXPORT_OK = qw(mem run_bwamem);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub mem {

##mem

##Function : Perl wrapper for writing bwa mem recipe to $FILEHANDLE. Based on bwa 0.7.15-r1140.
##Returns  : "@commands"
##Arguments: $infile_path, $idxbase, $second_infile_path, $outfile_path, $stderrfile_path, $FILEHANDLE, $thread_number, $read_group_header, $mark_split_as_secondary, $interleaved_fastq_file
##         : $infile_path             => Infile path (read 1 or interleaved i.e. read 1 and 2)
##         : $idxbase                 => Idxbase (human genome references and bwa mem idx files)
##         : $second_infile_path      => Second infile path (read 2)
##         : $outfile_path            => Outfile path
##         : $stderrfile_path         => Stderrfile path
##         : $FILEHANDLE              => Sbatch filehandle to write to
##         : $thread_number           => Number of threads
##         : $read_group_header       => Read group header line, such as '@RG\tID:foo\tSM:bar'
##         : $mark_split_as_secondary => Mark shorter split hits as secondary
##         : $interleaved_fastq_file  => Smart pairing

    my ($arg_href) = @_;

    ## Default(s)
    my $mark_split_as_secondary;
    my $interleaved_fastq_file;

    ## Flatten argument(s)
    my $infile_path;
    my $idxbase;
    my $second_infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $thread_number;
    my $read_group_header;

    my $tmpl = {
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path },
	idxbase => { required => 1, defined => 1, strict_type => 1, store => \$idxbase },
	second_infile_path => { strict_type => 1, store => \$second_infile_path },
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	thread_number => { allow => qr/^\d+$/,
			   strict_type => 1, store => \$thread_number },
	read_group_header => { strict_type => 1, store => \$read_group_header },
	mark_split_as_secondary => { default => 0,
				     allow => [0, 1],
				     strict_type => 1, store => \$mark_split_as_secondary },
	interleaved_fastq_file => { default => 0,
				    allow => [undef, 0, 1],
				    strict_type => 1, store => \$interleaved_fastq_file },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Bwa
    my @commands = qw(bwa mem);  #Stores commands depending on input parameters

    ##Options
    if ($thread_number) {

	push(@commands, "-t ".$thread_number);  #Number of threads
    }
    if ($interleaved_fastq_file) {

	push(@commands, "-p");  #Interleaved fastq file
    }
    if ($mark_split_as_secondary) {

	push(@commands, "-M");  #Mark shorter split hits as secondary
    }
    if ($read_group_header) {

	push(@commands, "-R ".$read_group_header);  #Read group header line
    }

    ##Human reference genome and bwa mem files
    push(@commands, $idxbase);

    ## Infile
    push(@commands, $infile_path);

    ## Read 2
    if ($second_infile_path) {
	
	push(@commands, $second_infile_path);
    }
    if ($outfile_path) {
	
	push(@commands, "> ".$outfile_path);  #Outfile
    }
    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


sub run_bwamem {

##run_bwamem

##Function : Perl wrapper for writing run_bwamem recipe to $FILEHANDLE. Based on bwakit 0.7.12.
##Returns  : "@commands"
##Arguments: $infile_path, $idxbase, $outfiles_prefix_path, $second_infile_path, $stderrfile_path, $FILEHANDLE, $thread_number, $read_group_header, $hla_typing
##         : $infile_path          => Infile path (read 1 or interleaved i.e. read 1 and 2)
##         : $idxbase              => Idxbase (human genome references and bwa mem idx files)
##         : $outfiles_prefix_path => Prefix for output files
##         : $second_infile_path   => Second infile path (read 2)
##         : $stderrfile_path      => Stderrfile path
##         : $FILEHANDLE           => Sbatch filehandle to write to
##         : $thread_number        => Number of threads
##         : $read_group_header    => Read group header line, such as '@RG\tID:foo\tSM:bar'
##         : $hla_typing           => Apply HLA typing

    my ($arg_href) = @_;

    ## Default(s)
    my $hla_typing;

    ## Flatten argument(s)
    my $infile_path;
    my $idxbase;
    my $outfiles_prefix_path;
    my $second_infile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $thread_number;
    my $read_group_header;

    my $tmpl = {
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path },
	idxbase => { required => 1, defined => 1, strict_type => 1, store => \$idxbase },
	outfiles_prefix_path => { required => 1, defined => 1, strict_type => 1, store => \$outfiles_prefix_path },
	second_infile_path => { strict_type => 1, store => \$second_infile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	thread_number => { allow => qr/^\d+$/,
			   strict_type => 1, store => \$thread_number },
	read_group_header => { strict_type => 1, store => \$read_group_header },
	hla_typing => { default => 0,
			allow => [undef, 0, 1],
			strict_type => 1, store => \$hla_typing },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Run_bwamem
    my @commands = qw(run-bwamem);  #Stores commands depending on input parameters

    ##Options
    if ($thread_number) {

	push(@commands, "-t ".$thread_number);  #Number of threads
    }
    if ($hla_typing) {

	push(@commands, "-H");  #Apply HLA typing
    }
    if ($read_group_header) {
	
	push(@commands, "-R ".$read_group_header);  #Read group header line
    }
    if ($outfiles_prefix_path) {
	
	push(@commands, "-o ".$outfiles_prefix_path);  #Outfiles prefix
    }

    ##Human reference genome and bwa mem files
    push(@commands, $idxbase);

    ## Infile
    push(@commands, $infile_path);

    ## Read 2
    if ($second_infile_path) {
	
	push(@commands, $second_infile_path);
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
