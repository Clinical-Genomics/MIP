package Program::Gzip;

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
    our @EXPORT_OK = qw(gzip);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub gzip {

##gzip

##Function : Perl wrapper for writing gzip recipe to $FILEHANDLE. Based on gzip 1.3.12.
##Returns  : ""
##Arguments: $FILEHANDLE, $stdout, $decompress, $infile_path, $outfile_path, $quiet, $verbose
##         : $FILEHANDLE   => Filehandle to write to
##         : $stdout       => Write on standard output, keep original files unchanged
##         : $decompress   => Decompress
##         : $infile_path  => Infile path
##         : $outfile_path => Outfile path. Write documents to FILE 
##         : $quiet        => Suppress all warnings
##         : $verbose      => Verbosity

    my ($arg_href) = @_;

    ## Default(s)
    my $quiet;
    my $verbose;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $stdout;
    my $decompress;
    my $infile_path;
    my $outfile_path;

    my $tmpl = {
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
	stdout => { strict_type => 1, store => \$stdout},
	decompress => { strict_type => 1, store => \$decompress},
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	outfile_path => { strict_type => 1, store => \$outfile_path},
	quiet => { default => 0,
		   strict_type => 1, store => \$quiet},
	verbose => { default => 0,
		     strict_type => 1, store => \$verbose},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Gzip
    print $FILEHANDLE "gzip ";

    ## Options
    if ($quiet) {

	print $FILEHANDLE "--quiet ";
    }
    if ($verbose) {

	print $FILEHANDLE "--verbose ";
    }
    if ($decompress) {

	print $FILEHANDLE "--decompress ";
    }
    if ($stdout) {  #Write to stdout stream

	print $FILEHANDLE "--stdout ";
    }

    ## FILE
    print $FILEHANDLE $infile_path." ";

    ## Outfile
    if ($outfile_path) {

	print $FILEHANDLE "> ".$outfile_path;  #Outfile
    }
}


1;
