package Program::Variantcalling::Freebayes;

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
    our @EXPORT_OK = qw(calling);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub calling {

##calling

##Function : Perl wrapper for writing freebayes recipe to $FILEHANDLE or return commands array. Based on freebayes 1.0.2.
##Returns  : "@commands"
##Arguments: $infile_paths_ref, $referencefile_path, $outfile_path, $stderrfile_path, $FILEHANDLE, $apply_standard_filter, $calculate_genotype_quality
##         : $infile_paths_ref           => Infile paths {REF}
##         : $referencefile_path         => Reference sequence file
##         : $outfile_path               => Outfile path
##         : $stderrfile_path            => Stderrfile path
##         : $FILEHANDLE                 => Filehandle to write to
##         : $apply_standard_filter      => Use stringent input base and mapping quality filters. Equivalent to -m 30 -q 20 -R 0 -S 0
##         : $calculate_genotype_quality => Calculate the marginal probability of genotypes and report as GQ in each sample field in the VCF output

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_paths_ref;
    my $referencefile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $apply_standard_filter;
    my $calculate_genotype_quality;

    my $tmpl = {
	infile_paths_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$infile_paths_ref},
	referencefile_path => { required => 1, defined => 1, strict_type => 1, store => \$referencefile_path },
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	apply_standard_filter => { default => 0,
			     allow => [0, 1],
			     strict_type => 1, store => \$apply_standard_filter },
	calculate_genotype_quality => { default => 0,
					allow => [0, 1],
					strict_type => 1, store => \$calculate_genotype_quality },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## freebayes
    my @commands = qw(freebayes);  #Stores commands depending on input parameters

    ## Options
    if ($apply_standard_filter) {

	push(@commands, "--standard-filters");  #Equivalent to -m 30 -q 20 -R 0 -S 0
    }
    if ($calculate_genotype_quality) {

	push(@commands, "--genotype-qualities");
    }
    if ($referencefile_path) {

	push(@commands, "--fasta-reference ".$referencefile_path);  #Reference sequence file
    }

    ## Infile
    push(@commands, join(" ", @$infile_paths_ref));

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
