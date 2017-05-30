package Program::Htslib;

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
    our @EXPORT_OK = qw(bgzip tabix);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub bgzip {

##bgzip

##Function : Perl wrapper for writing bgzip recipe to $FILEHANDLE or return commands array. Based on htslib 1.3.1.
##Returns  : "@commands"
##Arguments: $infile_path, $stderrfile_path, $FILEHANDLE, $output_type, $decompress, $write_to_stdout
##         : $infile_path     => Infile path to read from
##         : $outfile_path    => Outfile path to write to
##         : $stderrfile_path => Stderr file path to write to {OPTIONAL}
##         : $FILEHANDLE      => Filehandle to write to
##         : $decompress      => Decompress file
##         : $write_to_stdout => Write on standard output, keep original files unchanged

    my ($arg_href) = @_;

    ## Default(s)
    my $decompress;
    my $write_to_stdout;

    ## Flatten argument(s)
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;

    my $tmpl = {
	infile_path => { strict_type => 1, store => \$infile_path},
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	FILEHANDLE => { store => \$FILEHANDLE},
	decompress => { default => 0,
		     allow => [0, 1],
		     strict_type => 1, store => \$decompress},
	write_to_stdout => { default => 0,
			     allow => [0, 1],
			     strict_type => 1, store => \$write_to_stdout},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## bgzip
    my @commands = qw(bgzip);  #Stores commands depending on input parameters

    ## Options
    if ($decompress) {

	push(@commands, "--decompress");
    }
    if ($write_to_stdout) {

	push(@commands, "--stdout");
    }

    ## Infile
    if ($infile_path) {

	push(@commands, $infile_path);
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


sub tabix {

##tabix

##Function : Perl wrapper for writing tabix recipe to $FILEHANDLE or return commands array. Based on htslib 1.3.1.
##Returns  : "@commands"
##Arguments: $regions_ref, $infile_path, $stderrfile_path, $FILEHANDLE, $output_type, $force, $preset, $with_header
##         : $regions_ref     => The regions to process {REF}
##         : $infile_path     => Infile path to read from
##         : $outfile_path    => Outfile path to write to
##         : $stderrfile_path => Stderr file path to write to {OPTIONAL}
##         : $FILEHANDLE      => Filehandle to write to
##         : $force           => Overwrite existing index without asking
##         : $preset          => Preset
##         : $with_header     => Include header

    my ($arg_href) = @_;

    ## Default(s)
    my $with_header;
    my $force;
    my $preset;

    ## Flatten argument(s)
    my $regions_ref;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;

    my $tmpl = {
	regions_ref => { default => [], strict_type => 1, store => \$regions_ref },
	infile_path => { strict_type => 1, store => \$infile_path},
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	FILEHANDLE => { store => \$FILEHANDLE},
	force => { default => 0,
		   allow => [0, 1],
		   strict_type => 1, store => \$force},
	preset => { allow => [undef, "gff", "bed", "sam", "vcf"],
		    strict_type => 1, store => \$preset},
	with_header => { default => 0,
			 allow => [0, 1],
			 strict_type => 1, store => \$with_header },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## tabix
    my @commands = qw(tabix);  #Stores commands depending on input parameters

    ## Options
    if ($force) {

	push(@commands, "--force");
    }
    if ($preset) {

	push(@commands, "--preset ".$preset);
    }
    if ($with_header) {  #Include header

	push(@commands, "--print-header");
    }

    ## Infile
    if ($infile_path) {

	push(@commands, $infile_path);
    }
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
