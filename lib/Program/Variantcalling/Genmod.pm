package Program::Variantcalling::Genmod;

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
    our @EXPORT_OK = qw(annotate filter);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub annotate {

##annotate

##Function : Perl wrapper for writing Genmod annotate recipe to $FILEHANDLE or return commands array. Based on genmod 3.7.0.
##Returns  : "@commands"
##Arguments: $infile_path, $stderrfile_path, $FILEHANDLE, $output_type, $verbosity, $temp_directory_path, $thousand_g_file_path
##         : $infile_path          => Infile path to read from
##         : $outfile_path         => Outfile path to write to
##         : $stderrfile_path      => Stderr file path to write to {OPTIONAL}
##         : $FILEHANDLE           => Filehandle to write to
##         : $verbosity            => Increase output verbosity
##         : $temp_directory_path  => Directory for storing intermediate files
##         : $thousand_g_file_path => Specify the path to a bgzipped vcf file (with index) with 1000g variants

    my ($arg_href) = @_;

    ## Default(s)
    my $verbosity;

    ## Flatten argument(s)
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $temp_directory_path;
    my $thousand_g_file_path;

    my $tmpl = {
	infile_path => { strict_type => 1, store => \$infile_path},
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	FILEHANDLE => { store => \$FILEHANDLE},
	temp_directory_path => { strict_type => 1, store => \$temp_directory_path },
	thousand_g_file_path => { strict_type => 1, store => \$thousand_g_file_path },
	verbosity => { allow => qr/^\w+$/,
		       strict_type => 1, store => \$verbosity},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Genmod annotate
    my @commands = qw(genmod);  #Stores commands depending on input parameters

    ## Options
    if ($verbosity) {

	push(@commands, "-".$verbosity);
    }
    push(@commands, "annotate");

    if ($temp_directory_path) {

	push(@commands, "--temp_dir ".$temp_directory_path);
    }
    if ($thousand_g_file_path) {

	push(@commands, "--thousand-g ".$thousand_g_file_path);
    }

    if ($outfile_path) {

	push(@commands, "-o ".$outfile_path);  #Specify output filename
    }

    ## Infile
    if ($infile_path) {

	push(@commands, $infile_path);
    }
    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


sub filter {

##filter

##Function : Perl wrapper for writing Genmod filter recipe to $FILEHANDLE or return commands array. Based on genmod 3.7.0.
##Returns  : "@commands"
##Arguments: $infile_path, $stderrfile_path, $FILEHANDLE, $output_type, $verbosity, $threshold
##         : $infile_path     => Infile path to read from
##         : $outfile_path    => Outfile path to write to
##         : $stderrfile_path => Stderr file path to write to {OPTIONAL}
##         : $FILEHANDLE      => Filehandle to write to
##         : $verbosity       => Increase output verbosity
##         : $threshold       => Directory for storing intermediate files

    my ($arg_href) = @_;

    ## Default(s)
    my $verbosity;

    ## Flatten argument(s)
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $threshold;
    my $thousand_g_file_path;

    my $tmpl = {
	infile_path => { strict_type => 1, store => \$infile_path},
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	FILEHANDLE => { store => \$FILEHANDLE},
	threshold => { allow => qr/^\d+$|^\d+.\d+$/,
		       strict_type => 1, store => \$threshold },
	verbosity => { allow => qr/^\w+$/,
		       strict_type => 1, store => \$verbosity},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Genmod filter
    my @commands = qw(genmod);  #Stores commands depending on input parameters

    ## Options
    if ($verbosity) {

	push(@commands, "-".$verbosity);
    }
    push(@commands, "filter");

    if ($threshold) {

	push(@commands, "--threshold ".$threshold);
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
