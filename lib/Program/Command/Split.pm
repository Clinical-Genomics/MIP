package Program::Command::Split;

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
    our @EXPORT_OK = qw(split);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub split {

##split

##Function : Perl wrapper for writing split recipe to $FILEHANDLE or return commands array. Based on split 8.4.
##Returns  : "@commands"
##Arguments: $FILEHANDLE, $infile_path, $prefix, $lines, $suffix_length, $numeric_suffixes, $quiet, $verbose
##         : $FILEHANDLE       => Filehandle to write to
##         : $infile_path      => Infile path
##         : $prefix           => Prefix of output files
##         : $lines            => Put number lines per output file
##         : $suffix_length    => Use suffixes of length N
##         : $numeric_suffixes => Use numeric suffixes instead of alphabetic
##         : $quiet            => Suppress all warnings
##         : $verbose          => Verbosity

    my ($arg_href) = @_;

    ## Default(s)
    my $quiet;
    my $verbose;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $prefix;
    my $lines;
    my $suffix_length;
    my $numeric_suffixes;

    my $tmpl = {
	FILEHANDLE => { store => \$FILEHANDLE},
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	prefix => { strict_type => 1, store => \$prefix},
	lines => { allow => qr/^\d+$/,
		   strict_type => 1, store => \$lines},
	suffix_length => { allow => qr/^\d+$/,
			   strict_type => 1, store => \$suffix_length},
	numeric_suffixes => { default => 0,
			      allow => [0, 1],
			      strict_type => 1, store => \$numeric_suffixes},
	quiet => { default => 0,
		   allow => [0, 1],
		   strict_type => 1, store => \$quiet},
	verbose => { default => 0,
		     allow => [0, 1],
		     strict_type => 1, store => \$verbose},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Split
    my @commands = qw(split);  #Stores commands depending on input parameters

    ## Options
    if ($lines) {

	push(@commands, "--lines=".$lines);
    }
    if ($numeric_suffixes) {

	push(@commands, "--numeric-suffixes");
    }
    if ($numeric_suffixes) {
	
	push(@commands, "--suffix-length=".$numeric_suffixes);
    }
    if ($quiet) {

	push(@commands, "--quiet");
    }
    if ($verbose) {

	push(@commands, "--verbose");
    }

    ## FILE
    push(@commands, $infile_path);

    if ($prefix) {

	push(@commands, $prefix);
    }

    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


1;
