package Program::Variantcalling::Svdb;

use strict;
use warnings;
use warnings qw( FATAL utf8 );
use utf8;  #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );

BEGIN {
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.01;

    # Inherit from Exporter to export functions and variables
    our @ISA = qw(Exporter);

    # Functions and variables which are exported by default
    our @EXPORT = qw();

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(merge query);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub merge {

##merge

##Function : Perl wrapper for writing svdb merge recipe to $FILEHANDLE or return commands array. Based on svdb 1.0.7.
##Returns  : "@commands"
##Arguments: $infile_paths_ref, $outfile_path, $stderrfile_path, $FILEHANDLE, $priority, $notag
##         : $infile_paths_ref => Infile path {REF}
##         : $outfile_path     => Outfile path
##         : $stderrfile_path  => Stderr file path to write to {OPTIONAL}
##         : $FILEHANDLE       => Filehandle to write to
##         : $priority         => Priority order of structural variant calls
##         : $notag            => Do not add the the VARID and set entries to the info field

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_paths_ref;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $priority;
    my $notag;

    my $tmpl = {
	infile_paths_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$infile_paths_ref },
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	priority => { strict_type => 1, store => \$priority },
		notag => { strict_type => 1, store => \$notag },
	FILEHANDLE => { store => \$FILEHANDLE },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## svdb
    my @commands = qw(svdb --merge);  #Stores commands depending on input parameters

    ## Options
    if ($priority) {

	push(@commands, "--priority ".$priority);  #Priority order of structural variant calls
    }
    if ($notag) {

        # Do not tag variant with origin file
	push @commands, q{--notag };
    }

    ## Infile
    push(@commands, "--vcf ".join(" ", @{ $infile_paths_ref }));

    if ($outfile_path) {

	push(@commands, "> ".$outfile_path);  #Outfile prefix
    }
    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {

	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


sub query {

##query

##Function : Perl wrapper for writing svdb query recipe to $FILEHANDLE or return commands array. Based on svdb 0.1.2.
##Returns  : "@commands"
##Arguments: $infile_path, $dbfile_path, $outfile_path, $stderrfile_path, $FILEHANDLE, $hit_tag, $frequency_tag, $bnd_distance, $overlap
##         : $infile_path     => Infile path
##         : $outfile_path    => Outfile path
##         : $stderrfile_path => Stderr file path to write to {OPTIONAL}
##         : $FILEHANDLE      => Filehandle to write to
##         : $hit_tag         => The tag used to describe the number of hits within the info field of the output vcf
##         : $bnd_distance    => The maximum distance between two similar precise breakpoints
##         : $overlap         => The overlap required to merge two events

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $dbfile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $hit_tag;
    my $frequency_tag;
    my $bnd_distance;
    my $overlap;

    my $tmpl = {
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	dbfile_path => { required => 1, defined => 1, strict_type => 1, store => \$dbfile_path},
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	hit_tag => { strict_type => 1, store => \$hit_tag },
	frequency_tag => { strict_type => 1, store => \$frequency_tag },
	bnd_distance => { allow => qr/^\d+$/,
			  strict_type => 1, store => \$bnd_distance },
	overlap => { allow => qr/^\d+|d+.d+$/,
		     strict_type => 1, store => \$overlap },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## svdb
    my @commands = qw(svdb --query);  #Stores commands depending on input parameters

    ## Options
    if ($bnd_distance) {

	push(@commands, "--bnd_distance ".$bnd_distance);
    }
    if ($overlap) {

	push(@commands, "--overlap ".$overlap);
    }
    if ($hit_tag) {

	push(@commands, "--hit_tag ".$hit_tag);
    }
    if ($frequency_tag) {

	push(@commands, "--frequency_tag ".$frequency_tag);
    }

    push(@commands, "--db ".$dbfile_path);

    ## Infile
    push(@commands, "--query_vcf ".$infile_path);

    if ($outfile_path) {

	push(@commands, "> ".$outfile_path);  #Outfile prefix
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
