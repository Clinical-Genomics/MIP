package Program::Variantcalling::Cnvnator;

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
    our @EXPORT_OK = qw(read_extraction histogram statistics partition calling convert_to_vcf);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub read_extraction {

##read_extraction

##Function : Perl wrapper for writing cnvnator recipe to $FILEHANDLE or return commands array. Based on cnvnator 0.3.3.
##Returns  : "@commands"
##Arguments: $regions_ref, $infile_paths_ref, $outfile_path, $stderrfile_path, $stdoutfile_path, $FILEHANDLE, $unique
##         : $regions_ref      => The regions to process {REF}
##         : $infile_paths_ref => Infile paths {REF}
##         : $outfile_path     => Outfile path
##         : $stderrfile_path  => Stderrfile path
##         : $stdoutfile_path  => Stdoutfile path
##         : $FILEHANDLE       => Filehandle to write to
##         : $unique           => Ensure correct q0 field for CNV calls

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $regions_ref;
    my $infile_paths_ref;
    my $outfile_path;
    my $stderrfile_path;
    my $stdoutfile_path;
    my $FILEHANDLE;
    my $unique;

    my $tmpl = {
	regions_ref => { default => [], strict_type => 1, store => \$regions_ref },
	infile_paths_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$infile_paths_ref},
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	unique => { default => 0,
		    allow => [0, 1],
		    strict_type => 1, store => \$unique },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## cnvnator
    my @commands = qw(cnvnator);  #Stores commands depending on input parameters

    ## Options
    if(@$regions_ref) {  #Limit output to regions

	push(@commands, "-chrom", join(" ", @{ $regions_ref }));
    }
    if ($unique) {

	push(@commands, "-unique");
    }
    if ($outfile_path) {

	push(@commands, "-root ".$outfile_path);  #Specify output filename
    }

    ## Infile
    push(@commands, "-tree", join(" ", @$infile_paths_ref));

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


sub histogram {

##histogram

##Function : Perl wrapper for writing cnvnator recipe to $FILEHANDLE or return commands array. Based on cnvnator 0.3.3.
##Returns  : "@commands"
##Arguments: $regions_ref, $infile_path, $referencedirectory_path, $stderrfile_path, $stdoutfile_path, $FILEHANDLE, $cnv_bin_size 
##         : $regions_ref             => The regions to process {REF}
##         : $infile_path             => Infile paths
##         : $referencedirectory_path => Reference sequence file
##         : $stderrfile_path         => Stderrfile path
##         : $stdoutfile_path         => Stdoutfile path
##         : $FILEHANDLE              => Filehandle to write to
##         : $cnv_bin_size            => Copy number variant bin size

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $regions_ref;
    my $infile_path;
    my $referencedirectory_path;
    my $stderrfile_path;
    my $stdoutfile_path;
    my $FILEHANDLE;
    my $cnv_bin_size;

    my $tmpl = {
	regions_ref => { default => [], strict_type => 1, store => \$regions_ref },
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	referencedirectory_path => { required => 1, defined => 1, strict_type => 1, store => \$referencedirectory_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	cnv_bin_size => { allow => qr/^\d+$/,
			  strict_type => 1, store => \$cnv_bin_size },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## cnvnator
    my @commands = qw(cnvnator);  #Stores commands depending on input parameters

    ## Options
    if(@$regions_ref) {  #Limit output to regions

	push(@commands, "-chrom", join(" ", @{ $regions_ref }));
    }
    if ($referencedirectory_path) {

	push(@commands, "-d ".$referencedirectory_path);
    }
    if ($cnv_bin_size) {

	push(@commands, "-his ".$cnv_bin_size);
    }

    ## Infile
    push(@commands, "-root", $infile_path);

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


sub statistics {

##statistics

##Function : Perl wrapper for writing cnvnator recipe to $FILEHANDLE or return commands array. Based on cnvnator 0.3.3.
##Returns  : "@commands"
##Arguments: $regions_ref, $infile_path, $stderrfile_path, $stdoutfile_path, $FILEHANDLE, $cnv_bin_size 
##         : $regions_ref     => The regions to process {REF}
##         : $infile_path     => Infile paths
##         : $stderrfile_path => Stderrfile path
##         : $stdoutfile_path => Stdoutfile path
##         : $FILEHANDLE      => Filehandle to write to
##         : $cnv_bin_size    => Copy number variant bin size

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $regions_ref;
    my $infile_path;
    my $stderrfile_path;
    my $stdoutfile_path;
    my $FILEHANDLE;
    my $cnv_bin_size;

    my $tmpl = {
	regions_ref => { default => [], strict_type => 1, store => \$regions_ref },
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	cnv_bin_size => { allow => qr/^\d+$/,
			  strict_type => 1, store => \$cnv_bin_size },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## cnvnator
    my @commands = qw(cnvnator);  #Stores commands depending on input parameters

    ## Options
    if(@$regions_ref) {  #Limit output to regions

	push(@commands, "-chrom", join(" ", @{ $regions_ref }));
    }
    if ($cnv_bin_size) {

	push(@commands, "-stat ".$cnv_bin_size);
    }

    ## Infile
    push(@commands, "-root", $infile_path);

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


sub partition {

##statistics

##Function : Perl wrapper for writing cnvnator recipe to $FILEHANDLE or return commands array. Based on cnvnator 0.3.3.
##Returns  : "@commands"
##Arguments: $regions_ref, $infile_path, $stderrfile_path, $stdoutfile_path, $FILEHANDLE, $cnv_bin_size 
##         : $regions_ref     => The regions to process {REF}
##         : $infile_path     => Infile paths
##         : $stderrfile_path => Stderrfile path
##         : $stdoutfile_path => Stdoutfile path
##         : $FILEHANDLE      => Filehandle to write to
##         : $cnv_bin_size    => Copy number variant bin size

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $regions_ref;
    my $infile_path;
    my $stderrfile_path;
    my $stdoutfile_path;
    my $FILEHANDLE;
    my $cnv_bin_size;

    my $tmpl = {
	regions_ref => { default => [], strict_type => 1, store => \$regions_ref },
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	cnv_bin_size => { allow => qr/^\d+$/,
			  strict_type => 1, store => \$cnv_bin_size },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## cnvnator
    my @commands = qw(cnvnator);  #Stores commands depending on input parameters

    ## Options
    if(@$regions_ref) {  #Limit output to regions

	push(@commands, "-chrom", join(" ", @{ $regions_ref }));
    }
    if ($cnv_bin_size) {

	push(@commands, "-partition ".$cnv_bin_size);
    }

    ## Infile
    push(@commands, "-root", $infile_path);

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


sub calling {

##calling

##Function : Perl wrapper for writing cnvnator recipe to $FILEHANDLE or return commands array. Based on cnvnator 0.3.3.
##Returns  : "@commands"
##Arguments: $regions_ref, $infile_path, $outfile_path, $stderrfile_path, $stdoutfile_path, $FILEHANDLE, $cnv_bin_size 
##         : $regions_ref             => The regions to process {REF}
##         : $infile_path             => Infile paths
##         : $outfile_path            => Outfile path
##         : $stderrfile_path         => Stderrfile path
##         : $stdoutfile_path         => Stdoutfile path
##         : $FILEHANDLE              => Filehandle to write to
##         : $cnv_bin_size            => Copy number variant bin size

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $regions_ref;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stdoutfile_path;
    my $FILEHANDLE;
    my $cnv_bin_size;

    my $tmpl = {
	regions_ref => { default => [], strict_type => 1, store => \$regions_ref },
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	cnv_bin_size => { allow => qr/^\d+$/,
			  strict_type => 1, store => \$cnv_bin_size },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## cnvnator
    my @commands = qw(cnvnator);  #Stores commands depending on input parameters

    ## Options
    if(@$regions_ref) {  #Limit output to regions

	push(@commands, "-chrom", join(" ", @{ $regions_ref }));
    }
    if ($cnv_bin_size) {

	push(@commands, "-call ".$cnv_bin_size);
    }

    ## Infile
    push(@commands, "-root", $infile_path);

    if ($outfile_path) {

	push(@commands, "> ".$outfile_path);  #Specify output filename
    }
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


sub convert_to_vcf {

##convert_to_vcf

##Function : Perl wrapper for writing cnvnator recipe to $FILEHANDLE or return commands array. Based on cnvnator 0.3.3.
##Returns  : "@commands"
##Arguments: $infile_path, $referencedirectory_path, $outfile_path, $stderrfile_path, $FILEHANDLE
##         : $infile_path             => Infile paths
##         : $referencedirectory_path => Reference sequence file
##         : $outfile_path            => Stdoutfile path
##         : $stderrfile_path         => Stderrfile path
##         : $FILEHANDLE              => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $referencedirectory_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;

    my $tmpl = {
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	referencedirectory_path => { required => 1, defined => 1, strict_type => 1, store => \$referencedirectory_path },
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## cnvnator
    my @commands = qw(cnvnator2VCF.pl);  #Stores commands depending on input parameters

    ## Infile
    push(@commands, $infile_path);

    ## Options
    if ($referencedirectory_path) {

	push(@commands, $referencedirectory_path);
    }    
    if ($outfile_path) {

	push(@commands, "> ".$outfile_path);
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
