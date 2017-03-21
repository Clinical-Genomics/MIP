package Program::Variantcalling::Gatk;

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
    our @EXPORT_OK = qw(catvariants);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub catvariants {

##catvariants

##Function : Perl wrapper for writing GATK catvariants recipe to $FILEHANDLE. Based on GATK 3.7.0.
##Returns  : "@commands"
##Arguments: $gatk_path, $intervals_ref, $infile_paths_ref, $outfile_path, $referencefile_path, $stderrfile_path, $FILEHANDLE, $downsample_to_coverage, $gatk_disable_auto_index_and_file_lock, $logging_level, $assume_sorted
##         : $gatk_path => Path to java jar and analysis to run
##         : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##         : $infile_paths_ref                      => Infile paths {REF}
##         : $outfile_path                          => Outfile path
##         : $referencefile_path                    => Reference sequence file
##         : $stderrfile_path                       => Stderrfile path
##         : $FILEHANDLE                            => Sbatch filehandle to write to
##         : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##         : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##         : $logging_level                         => Set the minimum level of logging
##         : $assume_sorted                         =>Assume_sorted should be true if the input files are already sorted

    my ($arg_href) = @_;

    ## Default(s)
    my $gatk_disable_auto_index_and_file_lock;
    my $logging_level;
    my $assume_sorted;

    ## Flatten argument(s)
    my $gatk_path;
    my $intervals_ref;
    my $infile_paths_ref;
    my $outfile_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $downsample_to_coverage;
    
    my $tmpl = {
	gatk_path => { required => 1, defined => 1, strict_type => 1, store => \$gatk_path},
	intervals_ref => { default => [], strict_type => 1, store => \$intervals_ref},
	infile_paths_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$infile_paths_ref},
	outfile_path => { required => 1, defined => 1, strict_type => 1, store => \$outfile_path },
	referencefile_path => { required => 1, defined => 1, strict_type => 1, store => \$referencefile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	downsample_to_coverage => { allow => qr/^\d+$/,
				    strict_type => 1, store => \$downsample_to_coverage },
	gatk_disable_auto_index_and_file_lock => { default => 0,
						   allow => [0, 1],
						   strict_type => 1, store => \$gatk_disable_auto_index_and_file_lock },
	logging_level => { default => "INFO",
			   allow => ["INFO", "ERROR", "FATAL"],
			   strict_type => 1, store => \$logging_level },
	assume_sorted => { default => 1,
			   allow => [undef, 0, 1],
			   strict_type => 1, store => \$assume_sorted },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Gatk catvariants
    my @commands = qw(-cp);  #Stores commands depending on input parameters

    push(@commands, $gatk_path);

    ## Options
    push(@commands, "--logging_level ".$logging_level);

    if ($downsample_to_coverage) {

	push(@commands, "--downsample_to_coverage ".$downsample_to_coverage);
    }
    if ($gatk_disable_auto_index_and_file_lock) {

	push(@commands, "--disable_auto_index_creation_and_locking_when_reading_rods");
    }
    if (@$intervals_ref) {

	push(@commands, "--intervals ".join(" --intervals ", @$intervals_ref));
    }
    if ($referencefile_path) {

	push(@commands, "--reference ".$referencefile_path);  #Reference sequence file
    }
    if ($assume_sorted) {

	push(@commands, "--assumeSorted");  #AssumeSorted should be true if the input files are already sorted
    }

    ## Infile
    if (@$infile_paths_ref) {

	push(@commands, "--variant: ".join(" --variant: ", @$infile_paths_ref));
    }

    ## Output
    push(@commands, "--outputFile ".$outfile_path);  #Specify output filename

    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


1;
