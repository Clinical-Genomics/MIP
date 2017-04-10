package Program::Variantcalling::Vep;

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
    our @EXPORT_OK = qw(variant_effect_predictor);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub variant_effect_predictor {

##variant_effect_predictor

##Function : Perl wrapper for writing variant_effect_predictor recipe to $FILEHANDLE or return commands array. Based on VEP 87.
##Returns  : "@commands"
##Arguments: $plugins_ref, $regions_ref, $outfile_path, $infile_path, $stderrfile_path, $stdoutfile_path, $FILEHANDLE, $reference_path, $script_path, $assembly, $buffer_size, $infile_format, $outfile_format, $fork, $no_progress, $offline
##         : $plugins_ref      => Use named plugin {REF}
##         : $regions_ref      => The regions to process {REF}
##         : $vep_features_ref => Features to add to VEP
##         : $outfile_path     => Outfile path to write to
##         : $infile_path      => Infile path to read from
##         : $stderrfile_path  => Stderr file path to write to {OPTIONAL}
##         : $stdoutfile_path  => Stdoutfile path
##         : $FILEHANDLE       => Filehandle to write to
##         : $reference_path   => Reference sequence file
##         : $script_path      => Path to variant_effect_predictor script
##         : $assembly         => Assembly version to use
##         : $cache_directory  => VEP chache directory 
##         : $buffer_size      => Sets the internal buffer size, corresponding to the number of variants that are read in to memory simultaneously
##         : $infile_format    => Input file format - one of "ensembl", "vcf", "hgvs", "id"
##         : $outfile_format   => Output file format

    my ($arg_href) = @_;

    ## Default(s)
    my $infile_format;
    my $outfile_format;
    my $fork;

    ## Flatten argument(s)
    my $plugins_ref;
    my $regions_ref;
    my $vep_features_ref;
    my $outfile_path;
    my $infile_path;
    my $stderrfile_path;
    my $stdoutfile_path;
    my $FILEHANDLE;
    my $reference_path;
    my $script_path;
    my $assembly;
    my $cache_directory;
    my $buffer_size;

    my $tmpl = {
	plugins_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$plugins_ref },
	regions_ref => { default => [], strict_type => 1, store => \$regions_ref },
	vep_features_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$vep_features_ref },
	outfile_path => { strict_type => 1, store => \$outfile_path },
	infile_path => { strict_type => 1, store => \$infile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	reference_path => { strict_type => 1, store => \$reference_path },
	script_path => { strict_type => 1, store => \$script_path },
	assembly => { strict_type => 1, store => \$assembly },
	cache_directory => { strict_type => 1, store => \$cache_directory },
	buffer_size => { allow => qr/^\d+$/,
			 strict_type => 1, store => \$buffer_size },
	infile_format => { default => "vcf",
			   allow => ["ensembl", "vcf", "hgvs", "id"],
			   strict_type => 1, store => \$infile_format },
	outfile_format => { default => "vcf",
			    allow => ["vcf", "tab", "json"],
			    strict_type => 1, store => \$outfile_format },
	fork  => { default => 0,
		   allow => qr/^\d+$/,
		   strict_type => 1, store => \$fork },
    };
    
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## vep
    my @commands = qw(perl);  #Stores commands depending on input parameters

    if ($script_path) {

	push(@commands, $script_path);
    }
    else {

	push(@commands, "variant_effect_predictor.pl");
    }

    ## Options
    if ($fork) {

	push(@commands, "--fork ".$fork);
    }
    if ($buffer_size) {

	push(@commands, "--buffer_size ".$buffer_size);
    }
    if ($assembly) {

	push(@commands, "--assembly ".$assembly);  
    }
    if ($reference_path) {

	push(@commands, "--fasta ".$reference_path);  
    }
    if ($cache_directory) {

	push(@commands, "--dir_cache ".$cache_directory);
    }
    if ($infile_format) {

	push(@commands, "--format ".$infile_format);
    }
    if ($outfile_format) {

	push(@commands, "--".$outfile_format);
    }
    if(@$regions_ref) {  #Limit output to regions

        push(@commands, "--chr ".join(",", @{ $regions_ref }));
    }
    if (@$plugins_ref) {

	push(@commands, "--plugin ".join(" --plugin ", @$plugins_ref));
    }
    if (@$vep_features_ref) {

	push(@commands, "--".join(" --", @$vep_features_ref));
    }

    ## Infile
    if ($infile_path) {

	push(@commands, "--input_file ".$infile_path);
    }
    if ($outfile_path) {

	push(@commands, "--output_file ".$outfile_path);  #Specify output filename
    }
    if ($stdoutfile_path) {

        push(@commands, "1> ".$stdoutfile_path);  #Redirect stdout to program specific stdout file
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
