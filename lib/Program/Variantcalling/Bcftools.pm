package Program::Variantcalling::Bcftools;

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
    our @EXPORT_OK = qw(call filter norm);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub call {

##call

##Function : Perl wrapper for writing bcftools call recipe to $FILEHANDLE or return commands array. Based on bcftools 1.3.1.
##Returns  : "@commands"
##Arguments: $form_fields_ref, $outfile_path, $infile_path, $stderrfile_path, $FILEHANDLE, $samples_file, $multiallelic_caller, $output_type, $variants_only
##         : $form_fields_ref     => Output format fields {REF}
##         : $outfile_path        => Outfile path to write to
##         : $infile_path         => Infile path to read from
##         : $stderrfile_path     => Stderr file path to write to {OPTIONAL}
##         : $FILEHANDLE          => Filehandle to write to
##         : $samples_file        => PED file or a file with an optional column with sex
##         : $constrain           => One of: alleles, trio
##         : $multiallelic_caller => Alternative model for multiallelic and rare-variant calling
##         : $output_type         => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##         : $variants_only       => Output variant sites only

    my ($arg_href) = @_;

    ## Default(s)
    my $multiallelic_caller;
    my $output_type;
    my $variants_only;

    ## Flatten argument(s)
    my $form_fields_ref;
    my $outfile_path;
    my $infile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $samples_file;
    my $constrain;

    my $tmpl = {
	form_fields_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$form_fields_ref },
	outfile_path => { strict_type => 1, store => \$outfile_path },
	infile_path => { strict_type => 1, store => \$infile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	samples_file => { strict_type => 1, store => \$samples_file },
	constrain => { allow => [undef, "alleles", "trio"],
		       strict_type => 1, store => \$constrain },
	multiallelic_caller => { default => 1,
				 allow => [undef, 0, 1],
				 strict_type => 1, store => \$multiallelic_caller },
	output_type => { default => "v",
			 allow => ["b", "u", "z", "v"],
			 strict_type => 1, store => \$output_type },
	variants_only => { default => 1,
			   allow => [undef, 0, 1],
			   strict_type => 1, store => \$variants_only },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## bcftools
    my @commands = qw(bcftools call);  #Stores commands depending on input parameters

    ## Options
    if ($multiallelic_caller) {

	push(@commands, "--multiallelic-caller");  
    }
    if (@$form_fields_ref) {

	push(@commands, "--format-fields ".join(",", @$form_fields_ref));
    }
    if ($variants_only) {

	push(@commands, "--variants-only");  
    }
    if ($samples_file) {

	push(@commands, "--samples-file ".$samples_file);
    }
    if ($constrain) {

	push(@commands, "--constrain ".$constrain);
    }
    if ($output_type) {

	push(@commands, "--output-type ".$output_type);  #Specify output type
    }
    if ($outfile_path) {

	push(@commands, "--output ".$outfile_path);  #Specify output filename
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


sub filter {

##filter

##Function : Perl wrapper for writing bcftools filter recipe to $FILEHANDLE or return commands array. Based on bcftools 1.3.1.
##Returns  : "@commands"
##Arguments: $infile_path, $outfile_path, $stderrfile_path, $FILEHANDLE, $exclude, $soft_filter, $snp_gap, $indel_gap
##         : $infile_path      => Infile paths
##         : $outfile_path     => Outfile path
##         : $stderrfile_path  => Stderrfile path
##         : $FILEHANDLE       => Filehandle to write to
##         : $exclude          => Exclude sites for which the expression is true
##         : $soft_filter      => Annotate FILTER column with <string> or unique filter name
##         : $snp_gap          => Filter SNPs within <int> base pairs of an indel
##         : $indel_gap        => Filter clusters of indels separated by <int> or fewer base pairs allowing only one to pass

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $exclude;
    my $soft_filter;
    my $snp_gap;
    my $indel_gap;

    my $tmpl = {
	infile_path => { strict_type => 1, store => \$infile_path },
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	exclude => { strict_type => 1, store => \$exclude },
	soft_filter => { strict_type => 1, store => \$soft_filter },
	snp_gap => { allow => qr/^\d+$/,
		       strict_type => 1, store => \$snp_gap },
	indel_gap => { allow => qr/^\d+$/,
		       strict_type => 1, store => \$indel_gap },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## bcftools
    my @commands = qw(bcftools filter);  #Stores commands depending on input parameters

    ## Options
    if ($exclude) {

	push(@commands, "--exclude ".$exclude);  
    }
    if ($soft_filter) {

	push(@commands, "--soft-filter ".$soft_filter);  
    }
    if ($snp_gap) {

	push(@commands, "--SnpGap ".$snp_gap);
    }
    if ($indel_gap) {

	push(@commands, "--IndelGap ".$indel_gap);
    }

    ## Infile
    if ($infile_path) {

	push(@commands, $infile_path);
    }
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


sub norm {

##norm

##Function : Perl wrapper for writing bcftools norm recipe to $FILEHANDLE or return commands array. Based on bcftools 1.3.1.
##Returns  : "@commands"
##Arguments: $outfile_path, $reference_path, $infile_path, $stderrfile_path, $FILEHANDLE, $multiallelic, $output_type, $multiallelic_type
##         : $outfile_path      => Outfile path to write to
##         : $reference_path    => Human genome reference path
##         : $infile_path       => Infile path to read from
##         : $stderrfile_path   => Stderr file path to write to {OPTIONAL}
##         : $FILEHANDLE        => Filehandle to write to
##         : $multiallelic      => To split/join multiallelic calls or not
##         : $output_type       => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##         : $multiallelic_type => Type of multiallelic to split/join {OPTIONAL}

    my ($arg_href) = @_;

    ## Default(s)
    my $output_type;
    my $multiallelic_type;

    ## Flatten argument(s)
    my $outfile_path;
    my $reference_path;
    my $infile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $multiallelic;

    my $tmpl = {
	outfile_path => { required => 1, defined => 1, strict_type => 1, store => \$outfile_path},
	reference_path => { required => 1, defined => 1, strict_type => 1, store => \$reference_path},
	infile_path => { strict_type => 1, store => \$infile_path},
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	FILEHANDLE => { store => \$FILEHANDLE},
	multiallelic => { allow => ["+", "-"],
			  strict_type => 1, store => \$multiallelic},
	output_type => { default => "v",
			 allow => ["b", "u", "z", "v"],
			 strict_type => 1, store => \$output_type},
	multiallelic_type => { default => "both",
			       allow => ["snps", "indels", "both", "any"],
			       strict_type => 1, store => \$multiallelic_type},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## bcftools
    my @commands = qw(bcftools norm);  #Stores commands depending on input parameters

    ## Options
    if ($multiallelic) {

	push(@commands, "--multiallelics ".$multiallelic.$multiallelic_type);  
    }
    if ($reference_path) {

	push(@commands, "--fasta-ref ".$reference_path);  
    }
    if ($output_type) {

	push(@commands, "--output-type ".$output_type);  #Specify output type
    }
    if ($outfile_path) {

	push(@commands, "--output ".$outfile_path);  #Specify output filename
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
