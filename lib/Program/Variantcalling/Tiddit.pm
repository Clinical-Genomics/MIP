package Program::Variantcalling::Tiddit;

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
    our @EXPORT_OK = qw(sv);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub sv {

##sv

##Function : Perl wrapper for writing tiddit sv recipe to $FILEHANDLE or return commands array. Based on tiddit 1.0.2.
##Returns  : "@commands"
##Arguments: $infile_path, $outfile_path_prefix, $FILEHANDLE, $minimum_number_supporting_pairs
##         : $infile_path                     => Infile path
##         : $outfile_path_prefix             => Outfile path. Write documents to FILE
##         : $FILEHANDLE                      => Filehandle to write to
##         : $minimum_number_supporting_pairs => Minimum number of supporting pairs in order to call a variation event

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $outfile_path_prefix;
    my $FILEHANDLE;
    my $minimum_number_supporting_pairs;

    my $tmpl = {
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path },
	outfile_path_prefix => { strict_type => 1, store => \$outfile_path_prefix },
	FILEHANDLE => { store => \$FILEHANDLE },
	minimum_number_supporting_pairs => { allow => qr/^\d+$/,
					 strict_type => 1, store => \$minimum_number_supporting_pairs },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## tiddit
    my @commands = qw(TIDDIT --sv);  #Stores commands depending on input parameters

    ## Options
    if ($minimum_number_supporting_pairs) {

	push(@commands, "-p ".$minimum_number_supporting_pairs);  #Minimum number of supporting pairs in order to call a variation event
    }
    if ($outfile_path_prefix) {

	push(@commands, "-o ".$outfile_path_prefix);  #Outfile prefix
    }

    ## Infile
    push(@commands, "-b ".$infile_path);

    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


1;
