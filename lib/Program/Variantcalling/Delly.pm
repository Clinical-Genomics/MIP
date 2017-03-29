package Program::Variantcalling::Delly;

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
    our @EXPORT_OK = qw(call merge filter);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub call {

##call

##Function : Perl wrapper for writing Delly call recipe to $FILEHANDLE or return commands array. Based on Delly 0.7.6.
##Returns  : "@commands"
##Arguments: $infile_path, $referencefile_path, $sv_type, $outfile_path, $stderrfile_path, $stdoutfile_path, $genotypefile_path, $FILEHANDLE, $exclude_file_path
##         : $infile_path        => Infile paths {REF}
##         : $referencefile_path => Reference sequence file
##         : $sv_type            => Calculate the marginal probability of genotypes and report as GQ in each sample field in the VCF output
##         : $outfile_path       => Outfile path
##         : $stderrfile_path    => Stderrfile path
##         : $stdoutfile_path    => Stdoutfile path
##         : $genotypefile_path  => Input VCF/BCF file for re-genotyping
##         : $FILEHANDLE         => Filehandle to write to
##         : $exclude_file_path  => File with regions to exclude

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $referencefile_path;
    my $sv_type;
    my $outfile_path;
    my $stderrfile_path;
    my $stdoutfile_path;
    my $genotypefile_path;
    my $FILEHANDLE;
    my $exclude_file_path;

    my $tmpl = {
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	referencefile_path => { required => 1, defined => 1, strict_type => 1, store => \$referencefile_path },
	sv_type => { required => 1, defined => 1,
		     allow => ["DEL", "DUP", "INV", "INS", "TRA"],
		     strict_type => 1, store => \$sv_type },
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
	genotypefile_path => { strict_type => 1, store => \$genotypefile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	exclude_file_path => { strict_type => 1, store => \$exclude_file_path },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## delly
    my @commands = qw(delly call);  #Stores commands depending on input parameters

    ## Options
    if ($sv_type) {

	push(@commands, "--type ".$sv_type);
    }
    if ($exclude_file_path) {

	push(@commands, "--exclude ".$exclude_file_path); 
    }
    if ($referencefile_path) {

	push(@commands, "--genome ".$referencefile_path);  #Reference sequence file
    }
    if ($genotypefile_path) {

	push(@commands, "--vcffile ".$genotypefile_path);  #Specify output filename
    }
    if ($outfile_path) {

	push(@commands, "--outfile ".$outfile_path);  #Specify output filename
    }

    ## Infile
    push(@commands, $infile_path);

    if ($stdoutfile_path) {

	push(@commands, "1> ".$stdoutfile_path);  #Redirect stdout to program specific stderr file
    }
    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


sub merge {

##merge

##Function : Perl wrapper for writing Delly merge recipe to $FILEHANDLE or return commands array. Based on Delly 0.7.6.
##Returns  : "@commands"
##Arguments: $infile_paths_ref, $sv_type, $outfile_path, $stderrfile_path, $stdoutfile_path, $FILEHANDLE, $min_size, $max_size
##         : $infile_paths_ref => Infile paths {REF}
##         : $sv_type          => Calculate the marginal probability of genotypes and report as GQ in each sample field in the VCF output
##         : $outfile_path     => Outfile path
##         : $stderrfile_path  => Stderrfile path
##         : $stdoutfile_path  => Stdoutfile path
##         : $FILEHANDLE       => Filehandle to write to
##         : $min_size         => Min. SV size
##         : $max_size         => Max. SV size

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_paths_ref;
    my $sv_type;
    my $outfile_path;
    my $stderrfile_path;
    my $stdoutfile_path;
    my $FILEHANDLE;
    my $min_size;
    my $max_size;

    my $tmpl = {
	infile_paths_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$infile_paths_ref},
	sv_type => { required => 1, defined => 1,
		     allow => ["DEL", "DUP", "INV", "INS", "TRA"],
		     strict_type => 1, store => \$sv_type },
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	min_size => { allow => qr/^\d+$/,
		      strict_type => 1, store => \$min_size },
	max_size => { allow => qr/^\d+$/,
		      strict_type => 1, store => \$max_size },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## delly
    my @commands = qw(delly merge);  #Stores commands depending on input parameters

    ## Options
    if ($sv_type) {

	push(@commands, "--type ".$sv_type);
    }
    if (defined($min_size)) {

	push(@commands, "--minsize ".$min_size); 
    }
    if (defined($max_size)) {

	push(@commands, "--maxsize ".$max_size); 
    }
    if ($outfile_path) {

	push(@commands, "--outfile ".$outfile_path);  #Specify output filename
    }

    ## Infile
    push(@commands, join(" ", @$infile_paths_ref));

    if ($stdoutfile_path) {

	push(@commands, "1> ".$stdoutfile_path);  #Redirect stdout to program specific stderr file
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

##Function : Perl wrapper for writing Delly filter recipe to $FILEHANDLE or return commands array. Based on Delly 0.7.6.
##Returns  : "@commands"
##Arguments: $infile_path, $sv_type, $filter_mode, $outfile_path, $stderrfile_path, $stdoutfile_path, $FILEHANDLE, $min_size, $max_size
##         : $infile_path      => Infile paths
##         : $sv_type          => Calculate the marginal probability of genotypes and report as GQ in each sample field in the VCF output
##         : $filter_mode      => Filter mode
##         : $outfile_path     => Outfile path
##         : $stderrfile_path  => Stderrfile path
##         : $stdoutfile_path  => Stdoutfile path
##         : $FILEHANDLE       => Filehandle to write to
##         : $min_size         => Min. SV size
##         : $max_size         => Max. SV size

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $sv_type;
    my $filter_mode;
    my $outfile_path;
    my $stderrfile_path;
    my $stdoutfile_path;
    my $FILEHANDLE;
    my $min_size;
    my $max_size;

    my $tmpl = {
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	sv_type => { required => 1, defined => 1,
		     allow => ["DEL", "DUP", "INV", "INS", "TRA"],
		     strict_type => 1, store => \$sv_type },
	filter_mode => { required => 1, defined => 1,
			 allow => ["somatic", "germline"],
			 strict_type => 1, store => \$filter_mode },
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	min_size => { allow => qr/^\d+$/,
		      strict_type => 1, store => \$min_size },
	max_size => { allow => qr/^\d+$/,
		      strict_type => 1, store => \$max_size },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## delly
    my @commands = qw(delly filter);  #Stores commands depending on input parameters

    ## Options
    if ($sv_type) {

	push(@commands, "--type ".$sv_type);
    }
    if ($filter_mode) {

	push(@commands, "--filter ".$filter_mode);
    }
    if (defined($min_size)) {

	push(@commands, "--minsize ".$min_size); 
    }
    if (defined($max_size)) {

	push(@commands, "--maxsize ".$max_size); 
    }
    if ($outfile_path) {

	push(@commands, "--outfile ".$outfile_path);  #Specify output filename
    }

    ## Infile
    push(@commands, $infile_path);

    if ($stdoutfile_path) {

	push(@commands, "1> ".$stdoutfile_path);  #Redirect stdout to program specific stderr file
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
