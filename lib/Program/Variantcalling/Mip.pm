package Program::Variantcalling::Mip;

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
    our @EXPORT_OK = qw(vcfparser calculate_af max_af);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub calculate_af {

##calculate_af

##Function : Perl wrapper for writing calculate_af recipe to already open $FILEHANDLE or return commands array. Based on calculate_af 0.0.2.
##Returns  : "@commands"
##Arguments: $FILEHANDLE, $infile_path, $outfile_path, $stderrfile_path, $append_stderr_info
##         : $FILEHANDLE         => Filehandle to write to
##         : $infile_path        => Infile path
##         : $outfile_path       => Outfile path
##         : $stderrfile_path    => Stderrfile path
##         : $append_stderr_info => Append stderr info to file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $append_stderr_info;

    my $tmpl = {
	FILEHANDLE => { store => \$FILEHANDLE},
	infile_path => { strict_type => 1, store => \$infile_path},
	outfile_path => { strict_type => 1, store => \$outfile_path},
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	append_stderr_info => { default => 0,
				allow => [0, 1],
				strict_type => 1, store => \$append_stderr_info},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## calculate_af
    my @commands = qw(calculate_af);  #Stores commands depending on input parameters

    ## Options

    ## Infile
    if ($infile_path) {

	push(@commands, $infile_path);
    }
    ## Outfile
    if ($outfile_path) {

	push(@commands, "> ".$outfile_path);
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


sub max_af {

##max_af

##Function : Perl wrapper for writing max_af recipe to already open $FILEHANDLE or return commands array. Based on max_af 0.0.2.
##Returns  : "@commands"
##Arguments: $FILEHANDLE, $infile_path, $outfile_path, $stderrfile_path, $append_stderr_info
##         : $FILEHANDLE         => Filehandle to write to
##         : $infile_path        => Infile path
##         : $outfile_path       => Outfile path
##         : $stderrfile_path    => Stderrfile path
##         : $append_stderr_info => Append stderr info to file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $append_stderr_info;

    my $tmpl = {
	FILEHANDLE => { store => \$FILEHANDLE},
	infile_path => { strict_type => 1, store => \$infile_path},
	outfile_path => { strict_type => 1, store => \$outfile_path},
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	append_stderr_info => { default => 0,
				allow => [0, 1],
				strict_type => 1, store => \$append_stderr_info},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## max_af
    my @commands = qw(max_af);  #Stores commands depending on input parameters

    ## Options

    ## Infile
    if ($infile_path) {

	push(@commands, $infile_path);
    }
    ## Outfile
    if ($outfile_path) {

	push(@commands, "> ".$outfile_path);
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
