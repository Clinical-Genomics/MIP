package Program::Variantcalling::Delly;

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
    our @EXPORT_OK = qw(call);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub call {

##call

##Function : Perl wrapper for writing Delly call recipe to $FILEHANDLE or return commands array. Based on Delly 0.7.6.
##Returns  : "@commands"
##Arguments: $infile_path, $referencefile_path, $sv_type, $outfile_path, $stderrfile_path, $stdoutfile_path, $FILEHANDLE, $exclude_file_path
##         : $infile_path        => Infile paths {REF}
##         : $referencefile_path => Reference sequence file
##         : $sv_type            => Calculate the marginal probability of genotypes and report as GQ in each sample field in the VCF output
##         : $outfile_path       => Outfile path
##         : $stderrfile_path    => Stderrfile path
##         : $stdoutfile_path    => Stdoutfile path
##         : $FILEHANDLE         => Filehandle to write to
##         : $exclude_file_path  => File with regions to exclude

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $referencefile_path;
    my $sv_type;
    my $outfile_path;
    my $stderrfile_path;
    my $stdoutfile_path;
    my $FILEHANDLE;
    my $exclude_file_path;

    my $tmpl = {
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	referencefile_path => { required => 1, defined => 1, strict_type => 1, store => \$referencefile_path },
	sv_type => { required => 1, defined => 1,
		     allow => ["DEL", "DUP", "INV", "INS", "TRA"],
		     strict_type => 1, store => \$sv_type },	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	exclude_file_path => { strict_type => 1, store => \$exclude_file_path },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## delly
    my @commands = qw(delly call);  #Stores commands depending on input parameters

    ## Options
    if ($sv_type) {

	push(@commands, "--type ".$sv_type);
    }
    if ($exclude_file_path) {

	push(@commands, "--exclude ".$exclude_file_path); 
    }
    if ($referencefile_path) {

	push(@commands, "--genome ".$referencefile_path);  #Reference sequence file
    }
    if ($outfile_path) {

	push(@commands, "--outfile ".$outfile_path);  #Specify output filename
    }

    ## Infile
    push(@commands, $infile_path);

    if ($stdoutfile_path) {

	push(@commands, "1> ".$stdoutfile_path);  #Redirect stdout to program specific stderr file
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
