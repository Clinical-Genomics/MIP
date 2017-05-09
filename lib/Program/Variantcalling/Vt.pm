package Program::Variantcalling::Vt;

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
    our @EXPORT_OK = qw(decompose normalize vt_uniq);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub decompose {

##decompose

##Function : Perl wrapper for writing Vt decompose recipe to $FILEHANDLE or return commands array. Based on Vt v0.5.
##Returns  : "@commands"
##Arguments: $infile_path, $outfile_path, $stderrfile_path, $FILEHANDLE, $smart_decomposition
##         : $infile_path         => Infile path to read from
##         : $outfile_path        => Outfile path to write to
##         : $stderrfile_path     => Stderr file path to write to {OPTIONAL}
##         : $FILEHANDLE          => Filehandle to write to
##         : $smart_decomposition => Smart decomposition

    my ($arg_href) = @_;

    ## Default(s)
    my $smart_decomposition;

    ## Flatten argument(s)
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;

    my $tmpl = {
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	FILEHANDLE => { store => \$FILEHANDLE},
	smart_decomposition => { default => 0,
		     allow => [0, 1],
		     strict_type => 1, store => \$smart_decomposition},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Vt decompose
    my @commands = qw(vt decompose);  #Stores commands depending on input parameters

    ## Options
    if ($smart_decomposition) {

	push(@commands, "-s");
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


sub normalize {

##normalize

##Function : Perl wrapper for writing Vt normalize recipe to $FILEHANDLE or return commands array. Based on Vt v0.5.
##Returns  : "@commands"
##Arguments: $infile_path, $referencefile_path, $outfile_path, $stderrfile_path, $FILEHANDLE, $append_stderr_info, $no_fail_inconsistent_reference
##         : $infile_path                    => Infile path to read from
##         : $referencefile_path             => Reference sequence fasta file
##         : $outfile_path                   => Outfile path to write to
##         : $stderrfile_path                => Stderr file path to write to {OPTIONAL}
##         : $FILEHANDLE                     => Filehandle to write to
##         : $append_stderr_info             => Append stderr info to file
##         : $no_fail_inconsistent_reference => Do not fail when REF is inconsistent with reference sequence for non SNPs

    my ($arg_href) = @_;

    ## Default(s)
    my $no_fail_inconsistent_reference;

    ## Flatten argument(s)
    my $infile_path;
    my $referencefile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $append_stderr_info;

    my $tmpl = {
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	outfile_path => { strict_type => 1, store => \$outfile_path },
	referencefile_path => { required => 1, defined => 1, strict_type => 1, store => \$referencefile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	FILEHANDLE => { store => \$FILEHANDLE},
	append_stderr_info => { default => 0,
				allow => [0, 1],
				strict_type => 1, store => \$append_stderr_info},
	no_fail_inconsistent_reference => { default => 0,
					    allow => [0, 1],
					    strict_type => 1, store => \$no_fail_inconsistent_reference},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Vt normalize
    my @commands = qw(vt normalize);  #Stores commands depending on input parameters

    ## Options
    if ($no_fail_inconsistent_reference) {

	push(@commands, "-n");
    }
    if ($referencefile_path) {

	push(@commands, "-r ".$referencefile_path);
    }
    if ($outfile_path) {

	push(@commands, "-o ".$outfile_path);  #Specify output filename
    }

    ## Infile
    if ($infile_path) {

	push(@commands, $infile_path);
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


sub vt_uniq {

##vt_uniq

##Function : Perl wrapper for writing Vt uniq recipe to $FILEHANDLE or return commands array. Based on Vt v0.5.
##Returns  : "@commands"
##Arguments: $infile_path, $outfile_path, $stderrfile_path, $FILEHANDLE
##         : $infile_path         => Infile path to read from
##         : $outfile_path        => Outfile path to write to
##         : $stderrfile_path     => Stderr file path to write to {OPTIONAL}
##         : $FILEHANDLE          => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;

    my $tmpl = {
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	FILEHANDLE => { store => \$FILEHANDLE},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Vt uniq
    my @commands = qw(vt uniq);  #Stores commands depending on input parameters

    ## Options
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


1;
