package Program::Variantcalling::Variant_integrity;

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
    our @EXPORT_OK = qw(mendel father);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub mendel {

##mendel

##Function : Perl wrapper for writing Variant_integrity mendel recipe to $FILEHANDLE or return commands array. Based on Variant_integrity 3.7.0.
##Returns  : "@commands"
##Arguments: $infile_path, $family_file, $outfile_path, $stderrfile_path, $FILEHANDLE, $verbosity, $family_type
##         : $infile_path          => Infile path to read from
##         : $family_file          => Family file
##         : $outfile_path         => Outfile path to write to
##         : $stderrfile_path      => Stderr file path to write to {OPTIONAL}
##         : $FILEHANDLE           => Filehandle to write to
##         : $verbosity            => Increase output verbosity
##         : $family_type          => Setup of family file

    my ($arg_href) = @_;

    ## Default(s)
    my $verbosity;

    ## Flatten argument(s)
    my $infile_path;
    my $family_file;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $family_type;

    my $tmpl = {
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	family_file => { required => 1, defined => 1, strict_type => 1, store => \$family_file},
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	FILEHANDLE => { store => \$FILEHANDLE},
	family_type => { allow => ["ped", "alt", "cmms", "mip"],
			 strict_type => 1, store => \$family_type },
	verbosity => { allow => qr/^\w+$/,
		       strict_type => 1, store => \$verbosity},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Variant_integrity mendel
    my @commands = qw(variant_integrity);  #Stores commands depending on input parameters

    ## Options
    if ($verbosity) {

	push(@commands, "-".$verbosity);
    }
    if ($family_file) {

	push(@commands, "--family_file ".$family_file);
    }
    if ($family_type) {

	push(@commands, "--family_type ".$family_type);
    }
    if ($outfile_path) {

	push(@commands, "--outfile ".$outfile_path);  #Specify output filename
    }

    ## Infile
    if ($infile_path) {

	push(@commands, $infile_path);
    }

    push(@commands, "mendel");

    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


sub father {

##father

##Function : Perl wrapper for writing Variant_integrity father recipe to $FILEHANDLE or return commands array. Based on Variant_integrity 3.7.0.
##Returns  : "@commands"
##Arguments: $infile_path, $family_file, $outfile_path, $stderrfile_path, $FILEHANDLE, $verbosity, $family_type
##         : $infile_path          => Infile path to read from
##         : $family_file          => Family file
##         : $outfile_path         => Outfile path to write to
##         : $stderrfile_path      => Stderr file path to write to {OPTIONAL}
##         : $FILEHANDLE           => Filehandle to write to
##         : $verbosity            => Increase output verbosity
##         : $family_type          => Setup of family file

    my ($arg_href) = @_;

    ## Default(s)
    my $verbosity;

    ## Flatten argument(s)
    my $infile_path;
    my $family_file;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $family_type;

    my $tmpl = {
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	family_file => { required => 1, defined => 1, strict_type => 1, store => \$family_file},
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	FILEHANDLE => { store => \$FILEHANDLE},
	family_type => { allow => ["ped", "alt", "cmms", "mip"],
			 strict_type => 1, store => \$family_type },
	verbosity => { allow => qr/^\w+$/,
		       strict_type => 1, store => \$verbosity},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Variant_integrity father
    my @commands = qw(variant_integrity);  #Stores commands depending on input parameters

    ## Options
    if ($verbosity) {

	push(@commands, "-".$verbosity);
    }
    if ($family_file) {

	push(@commands, "--family_file ".$family_file);
    }
    if ($family_type) {

	push(@commands, "--family_type ".$family_type);
    }
    if ($outfile_path) {

	push(@commands, "--outfile ".$outfile_path);  #Specify output filename
    }

    ## Infile
    if ($infile_path) {

	push(@commands, $infile_path);
    }

    push(@commands, "father");

    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


1;
