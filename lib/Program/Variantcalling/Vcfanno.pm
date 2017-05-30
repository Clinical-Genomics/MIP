package Program::Variantcalling::Vcfanno;

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
    our @EXPORT_OK = qw(vcfanno);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub vcfanno {

##vcfanno

##Function : Perl wrapper for writing vcfanno recipe to $FILEHANDLE or return commands array. Based on vcfanno 0.1.0.
##Returns  : "@commands"
##Arguments: $infile_path, $toml_configfile_path, $outfile_path, $stderrfile_path, $FILEHANDLE, $lua, $ends
##         : $infile_path          => Infile path to read from
##         : $toml_configfile_path => Toml config file 
##         : $outfile_path         => Outfile path to write to
##         : $stderrfile_path      => Stderr file path to write to {OPTIONAL}
##         : $FILEHANDLE           => Filehandle to write to
##         : $lua                  => optional path to a file containing custom javascript functions to be used as ops
##         : $ends                 => annotate the start and end as well as the interval itself

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $toml_configfile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $ends;
    my $lua;

    my $tmpl = {
	infile_path => { strict_type => 1, store => \$infile_path},
	toml_configfile_path => { required => 1, defined => 1, strict_type => 1, store => \$toml_configfile_path },
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	FILEHANDLE => { store => \$FILEHANDLE},
	ends => { default => 0,
		  allow => [0, 1],
		  strict_type => 1, store => \$ends },
	lua => { strict_type => 1, store => \$lua},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Vcfanno
    my @commands = qw(vcfanno);  #Stores commands depending on input parameters

    ## Options
    if ($lua) {

	push(@commands, "-lua ".$lua);
    }
    if ($ends) {

	push(@commands, "-ends");
    }
    if ($toml_configfile_path) {

	push(@commands, $toml_configfile_path);
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


1;
