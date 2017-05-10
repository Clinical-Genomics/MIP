package Program::Qc::Multiqc;

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
    our @EXPORT_OK = qw(multiqc);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case

sub multiqc {

##multiqc

##Function : Perl wrapper for writing multiqc recipe to already open $FILEHANDLE or return commands array. Based on multiqc 0.8.dev0.
##Returns  : "@commands"
##Arguments: $indir_path, $outdir_path, $stderrfile_path, $stdoutfile_path, $FILEHANDLE, $append_stderr_info, $force
##         : $indir_path         => Indir path
##         : $outdir_path        => Outdir path
##         : $stderrfile_path    => Stderrfile path
##         : $stdoutfile_path    => Stdoutfile path
##         : $FILEHANDLE         => Filehandle to write to
##         : $append_stderr_info => Append stderr info to file
##         : $force              => Force overwrite of output files

    my ($arg_href) = @_;

    ## Default(s)
    my $append_stderr_info;
    my $force;
	
    ## Flatten argument(s)
    my $indir_path;
    my $outdir_path;
    my $stderrfile_path;
    my $stdoutfile_path;
    my $FILEHANDLE;

    my $tmpl = { 
	indir_path => { required => 1, defined => 1, strict_type => 1, store => \$indir_path},
	outdir_path => { strict_type => 1, store => \$outdir_path},
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
	FILEHANDLE => { store => \$FILEHANDLE},
	append_stderr_info => { default => 0,
				allow => [0, 1],
				strict_type => 1, store => \$append_stderr_info},
	force => { default => 0,
		   allow => [0, 1],
		   strict_type => 1, store => \$force},
    };
    
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## multiqc
    my @commands = qw(multiqc);  #Stores commands depending on input parameters

    ## Options
    if ($force) {

	push(@commands, "--force");
    }
    ## Outdir
    if ($outdir_path) {

	push(@commands, "--outdir ".$outdir_path);
    }
    ## Indir
    if ($indir_path) {

	push(@commands, $indir_path);
    }
    if ($stdoutfile_path) {

	push(@commands, "1> ".$stdoutfile_path);  #Redirect stdout to program specific stdout file
    }
    if ($stderrfile_path) {

	if ($append_stderr_info) {

	    push(@commands, "2>> ".$stderrfile_path);  #Redirect and append stderr output to program specific stderr file
	}
	else {
	    
	    push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
	}
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


1;
