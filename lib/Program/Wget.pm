package Program::Wget;

use strict;
use warnings;
use warnings qw( FATAL utf8 );
use utf8;  #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case

use Exporter qw(import);
 
our @EXPORT_OK = ("wget",
    );

sub wget {

##wget

##Function : Perl wrapper for writing wget recipe to $FILEHANDLE. Based on GNU Wget 1.12, a non-interactive network retriever.
##Returns  : ""
##Arguments: $FILEHANDLE, $outfile_path, $url, $quiet, $verbose
##         : $FILEHANDLE   => Filehandle to write to
##         : $outfile_path => Outfile path. Write documents to FILE 
##         : $url          => Url to use for download
##         : $quiet        => Quiet (no output)
##         : $verbose      => Verbosity

    my ($arg_href) = @_;

    ## Default(s)
    my $quiet;
    my $verbose;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $outfile_path;
    my $url;

    my $tmpl = {
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
	outfile_path => { strict_type => 1, store => \$outfile_path},
	url => { required => 1, defined => 1, strict_type => 1, store => \$url},
	quiet => { default => 0,
		   strict_type => 1, store => \$quiet},
	verbose => { default => 1,
		     strict_type => 1, store => \$verbose},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Wget
    print $FILEHANDLE "wget ";

    if ($quiet) {

	print $FILEHANDLE "--quiet ";
    }
    if ($verbose) {

	print $FILEHANDLE "--verbose ";
    }
    else {

	print $FILEHANDLE "--no-verbose ";
    }
    print $FILEHANDLE $url." ";

    if ($outfile_path) {

	print $FILEHANDLE "-O ".$outfile_path;  #Outfile
    }
    if ($outfile_path ne "-") { #Write to stdout stream

	print $FILEHANDLE "\n\n";
    }
}

return 1;
