package Program::Variantcalling::Picardtools;

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
    our @EXPORT_OK = qw(sortvcf genotypeconcordance);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub sortvcf {

##sortvcf

##Function : Perl wrapper for writing picardtools sortvcf recipe to $FILEHANDLE. Based on picardtools 2.5.0.
##Returns  : "@commands"
##Arguments: $infile_paths_ref, $outfile_path, $sequence_dictionary, $FILEHANDLE, $stderrfile_path, $create_index
##         : $infile_paths_ref    => Infile paths {REF}
##         : $outfile_path        => Outfile path
##         : $sequence_dictionary => Sequence dictionary
##         : $FILEHANDLE          => Sbatch filehandle to write to
##         : $stderrfile_path     => Stderrfile path
##         : $create_index        => Create index

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_paths_ref;
    my $outfile_path;
    my $sequence_dictionary;
    my $FILEHANDLE;
    my $stderrfile_path;
    my $create_index;
    
    my $tmpl = {
	infile_paths_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$infile_paths_ref},
	outfile_path => { required => 1, defined => 1, strict_type => 1, store => \$outfile_path },
	sequence_dictionary => { required => 1, defined => 1, strict_type => 1, store => \$sequence_dictionary },
	FILEHANDLE => { store => \$FILEHANDLE },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	create_index => { allow => ["true", "false"],
			  strict_type => 1, store => \$create_index },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Picardtools sortvcf
    my @commands = qw(SortVcf);  #Stores commands depending on input parameters

    if ($sequence_dictionary) {

	push(@commands, "SEQUENCE_DICTIONARY=".$sequence_dictionary);
    }
    if ($create_index) {

	push(@commands, "CREATE_INDEX=".$create_index);
    }

    ## Infile
    push(@commands, "INPUT=".join(" INPUT=", @$infile_paths_ref));

    ## Output
    push(@commands, "OUTPUT=".$outfile_path);  #Specify output filename

    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


sub genotypeconcordance {

##genotypeconcordance

##Function : Perl wrapper for writing picardtools genotypeconcordance recipe to $FILEHANDLE. Based on picardtools 2.5.0.
##Returns  : "@commands"
##Arguments: $intervals_ref, $infile_path, $truth_file_path, $outfile_prefix_path, $truth_sample, $call_sample, $FILEHANDLE, $stderrfile_path, $min_genotype_quality, $min_depth
##         : $intervals_ref        => One or more genomic intervals over which to operate {REF}
##         : $infile_path          => Infile paths
##         : $truth_file_path      => VCF containing the truth sample
##         : $outfile_prefix_path  => Outfile path
##         : $truth_sample         => Name of the truth sample within the truth VCF
##         : $call_sample          => Name of the call sample within the call VCF
##         : $FILEHANDLE           => Sbatch filehandle to write to
##         : $stderrfile_path      => Stderrfile path
##         : $min_genotype_quality => Genotypes below this genotype quality will have genotypes classified as LowGq
##         : $min_depth            => Genotypes below this depth will have genotypes classified as LowDp

    my ($arg_href) = @_;

    ## Default(s)
    my $min_genotype_quality;
    my $min_depth;

    ## Flatten argument(s)
    my $intervals_ref;
    my $infile_path;
    my $truth_file_path;
    my $outfile_prefix_path;
    my $truth_sample;
    my $call_sample;
    my $FILEHANDLE;
    my $stderrfile_path;
    
    my $tmpl = {
	intervals_ref => { default => [], strict_type => 1, store => \$intervals_ref},
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	truth_file_path  => { required => 1, defined => 1, strict_type => 1, store => \$truth_file_path},
	outfile_prefix_path => { required => 1, defined => 1, strict_type => 1, store => \$outfile_prefix_path },
	truth_sample => { required => 1, defined => 1, strict_type => 1, store => \$truth_sample },
	call_sample => { required => 1, defined => 1, strict_type => 1, store => \$call_sample },
	FILEHANDLE => { store => \$FILEHANDLE },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	min_genotype_quality => { default => 0,
				  allow => qr/^\d+$/,
				  strict_type => 1, store => \$min_genotype_quality},
	min_depth  => { default => 0,
			allow => qr/^\d+$/,
			strict_type => 1, store => \$min_depth},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Picardtools genotypeconcordance
    my @commands = qw(GenotypeConcordance);  #Stores commands depending on input parameters

    ##Options
    if ($min_genotype_quality) {

	push(@commands, "MIN_GQ=".$min_genotype_quality);
    }
    if ($min_depth) {

	push(@commands, "MIN_DP=".$min_depth);
    }
    if (@$intervals_ref) {

	push(@commands, "INTERVALS=".join(" INTERVALS=", @$intervals_ref));
    }

    ## Infile
    push(@commands, "CALL_VCF=".$infile_path);

    push(@commands, "TRUTH_VCF=".$truth_file_path);

    ## Output
    push(@commands, "OUTPUT=".$outfile_prefix_path);

    push(@commands, "TRUTH_SAMPLE=".$truth_sample);

    push(@commands, "CALL_SAMPLE=".$call_sample);

    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


1;
