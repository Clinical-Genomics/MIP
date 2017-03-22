package Program::Variantcalling::Manta;

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
    our @EXPORT_OK = qw(config workflow);

}
use File::Spec::Functions qw(catfile);
use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub config {

##config

##Function : Perl wrapper for writing Manta config recipe to $FILEHANDLE or return commands array. Based on Manta 1.0.0.
##Returns  : "@commands"
##Arguments: $infile_paths_ref, $referencefile_path, $outdirectory_path, $stderrfile_path, $FILEHANDLE, $exome_analysis
##         : $infile_paths_ref   => Infile paths {REF}
##         : $referencefile_path => Reference sequence file
##         : $outdirectory_path  => Outfile path
##         : $stderrfile_path    => Stderrfile path
##         : $FILEHANDLE         => Filehandle to write to
##         : $exome_analysis     => Set options for WES input: turn off depth filters

    my ($arg_href) = @_;

    ## Default(s)
    my $exome_analysis;

    ## Flatten argument(s)
    my $infile_paths_ref;
    my $referencefile_path;
    my $outdirectory_path;
    my $stderrfile_path;
    my $FILEHANDLE;

    my $tmpl = {
	infile_paths_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$infile_paths_ref},
	referencefile_path => { required => 1, defined => 1, strict_type => 1, store => \$referencefile_path },
	outdirectory_path => { strict_type => 1, store => \$outdirectory_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	exome_analysis => { default => 0,
			    allow => [undef, 0, 1],
			    strict_type => 1, store => \$exome_analysis},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Manta
    my @commands = qw(configManta.py);  #Stores commands depending on input parameters

    ## Options
    if ($referencefile_path) {

	push(@commands, "--referenceFasta ".$referencefile_path);  #Reference sequence file
    }
    if ($exome_analysis) {

	push(@commands, "--exome");
    }

    ## Infile
    push(@commands, "--bam ".join(" --bam ", @$infile_paths_ref));

    if ($outdirectory_path) {

	push(@commands, "--runDir ".$outdirectory_path);
    }
    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


sub workflow {

##workflow

##Function : Perl wrapper for writing Manta workflow recipe to $FILEHANDLE or return commands array. Based on Manta 1.0.0.
##Returns  : "@commands"
##Arguments: $outdirectory_path, $stderrfile_path, $FILEHANDLE, $mode
##         : $outdirectory_path  => Outfile path
##         : $stderrfile_path    => Stderrfile path
##         : $FILEHANDLE         => Filehandle to write to
##         : $mode               => Mode of parallel 

    my ($arg_href) = @_;

    ## Default(s)
    my $mode;

    ## Flatten argument(s)
    my $outdirectory_path;
    my $stderrfile_path;
    my $FILEHANDLE;

    my $tmpl = {
	outdirectory_path => { required => 1, defined => 1, strict_type => 1, store => \$outdirectory_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	mode => { default => "local",
		  allow => [undef, "local", "sge"],
		  strict_type => 1, store => \$mode},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Manta
    my @commands = (catfile($outdirectory_path, "runWorkflow.py"));  #Stores commands depending on input parameters

    ## Options
    if ($mode) {

	push(@commands, "--mode ".$mode);
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
