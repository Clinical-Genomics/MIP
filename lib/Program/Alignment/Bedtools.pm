package Program::Alignment::Bedtools;

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
    our @EXPORT_OK = qw(genomecov);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub genomecov {

##genomecov

##Function : Perl wrapper for writing bedtools genomecov recipe to $FILEHANDLE. Based on bedtools 2.26.0.
##Returns  : "@commands"
##Arguments: $infile_path, $outfile_path, $referencefile_path, $stderrfile_path, $FILEHANDLE, $max_coverage
##         : $infile_path        => Infile paths
##         : $outfile_path       => Outfile path
##         : $referencefile_path => Genome reference file
##         : $stderrfile_path    => Stderrfile path
##         : $FILEHANDLE         => Sbatch filehandle to write to
##         : $max_coverage       => Combine all positions with a depth >= max into a single bin in the histogram

    my ($arg_href) = @_;

    ## Default(s)
    my $max_coverage;

    ## Flatten argument(s)
    my $infile_path;
    my $outfile_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    
    my $tmpl = {
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	outfile_path => { strict_type => 1, store => \$outfile_path },
	referencefile_path => { required => 1, defined => 1, strict_type => 1, store => \$referencefile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	max_coverage => { allow => qr/^\d+$/,
			  strict_type => 1, store => \$max_coverage},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Bedtools genomecov
    my @commands = qw(bedtools genomecov);  #Stores commands depending on input parameters

    ## Options
    if (defined($max_coverage)) {

	push(@commands, "-max ".$max_coverage);
    }

    ## Infile
    push(@commands, "-ibam ".$infile_path);

    if ($referencefile_path) {

	push(@commands, "-g ".$referencefile_path);
    }

    ## Output
    if($outfile_path) {

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
