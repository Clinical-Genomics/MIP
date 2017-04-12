package Program::Variantcalling::Genmod;

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
    our @EXPORT_OK = qw(annotate models score compound filter);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub annotate {

##annotate

##Function : Perl wrapper for writing Genmod annotate recipe to $FILEHANDLE or return commands array. Based on genmod 3.7.0.
##Returns  : "@commands"
##Arguments: $cadd_file_paths_ref, $infile_path, $outfile_path, $stderrfile_path, $FILEHANDLE, $verbosity, $temp_directory_path, $annotate_region, $thousand_g_file_path, $spidex_file_path
##         : $cadd_file_paths_ref       => Specify the path to a bgzipped cadd file (with index) with variant scores
##         : $infile_path          => Infile path to read from
##         : $outfile_path         => Outfile path to write to
##         : $stderrfile_path      => Stderr file path to write to {OPTIONAL}
##         : $FILEHANDLE           => Filehandle to write to
##         : $verbosity            => Increase output verbosity
##         : $temp_directory_path  => Directory for storing intermediate files
##         : $annotate_region      => Annotate what regions a variant belongs to
##         : $thousand_g_file_path => Specify the path to a bgzipped vcf file (with index) with 1000g variants
##         : $spidex_file_path     => Specify the path to a bgzipped tsv file (with index) with spidex information

    my ($arg_href) = @_;

    ## Default(s)
    my $verbosity;
    my $annotate_region;

    ## Flatten argument(s)
    my $cadd_file_paths_ref;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $temp_directory_path;
    my $thousand_g_file_path;
    my $spidex_file_path;

    my $tmpl = {
	cadd_file_paths_ref => { default => [], strict_type => 1, store => \$cadd_file_paths_ref},
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	FILEHANDLE => { store => \$FILEHANDLE},
	temp_directory_path => { strict_type => 1, store => \$temp_directory_path },
	thousand_g_file_path => { strict_type => 1, store => \$thousand_g_file_path },
	spidex_file_path => { strict_type => 1, store => \$spidex_file_path },
	verbosity => { allow => qr/^\w+$/,
		       strict_type => 1, store => \$verbosity},
	annotate_region => { default => 0,
			     allow => [0, 1],
			     strict_type => 1, store => \$annotate_region},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Genmod annotate
    my @commands = qw(genmod);  #Stores commands depending on input parameters

    ## Options
    if ($verbosity) {

	push(@commands, "-".$verbosity);
    }
    push(@commands, "annotate");

    if ($temp_directory_path) {

	push(@commands, "--temp_dir ".$temp_directory_path);
    }
    if ($annotate_region) {

	push(@commands, "--annotate_regions");
    }
    if ($thousand_g_file_path) {

	push(@commands, "--thousand-g ".$thousand_g_file_path);
    }
    if ($spidex_file_path) {

	push(@commands, "--spidex ".$spidex_file_path);
    }
    if (@$cadd_file_paths_ref) {

	push(@commands, "--cadd_file ".join(" --cadd_file ", @{ $cadd_file_paths_ref }));
    }

    if ($outfile_path) {

	push(@commands, "--outfile ".$outfile_path);  #Specify output filename
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


sub models {

##models

##Function : Perl wrapper for writing Genmod models recipe to $FILEHANDLE or return commands array. Based on genmod 3.7.0.
##Returns  : "@commands"
##Arguments: $infile_path, $family_file, $outfile_path, $stderrfile_path, $FILEHANDLE, $verbosity, $temp_directory_path, $reduced_penetrance_file_path, $family_type, $thread_number, $vep, $whole_gene
##         : $infile_path                  => Infile path to read from
##         : $family_file                  => Family file
##         : $outfile_path                 => Outfile path to write to
##         : $stderrfile_path              => Stderr file path to write to {OPTIONAL}
##         : $FILEHANDLE                   => Filehandle to write to
##         : $verbosity                    => Increase output verbosity
##         : $temp_directory_path          => Directory for storing intermediate files
##         : $reduced_penetrance_file_path => Reduced penetrance file path
##         : $family_type                  => Setup of family file
##         : $thread_number                => Define how many processes that should be use for annotation
##         : $vep                          => If variants are annotated with the Variant Effect Predictor
##         : $whole_gene                   => If compounds should be checked for all variants in a gene

    my ($arg_href) = @_;

    ## Default(s)
    my $verbosity;
    my $thread_number;
    my $vep;
    my $whole_gene;

    ## Flatten argument(s)
    my $infile_path;
    my $family_file;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $temp_directory_path;
    my $reduced_penetrance_file_path;
    my $family_type;

    my $tmpl = {
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	family_file => { required => 1, defined => 1, strict_type => 1, store => \$family_file},
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	FILEHANDLE => { store => \$FILEHANDLE},
	temp_directory_path => { strict_type => 1, store => \$temp_directory_path },
	reduced_penetrance_file_path => { strict_type => 1, store => \$reduced_penetrance_file_path },
	family_type => { allow => ["ped", "alt", "cmms", "mip"],
			 strict_type => 1, store => \$family_type },
	verbosity => { allow => qr/^\w+$/,
		       strict_type => 1, store => \$verbosity},
	thread_number => { default => 1,
			   allow => qr/^\d+$/,
			   strict_type => 1, store => \$thread_number},
	vep => { default => 0,
		 allow => [undef, 0, 1],
		 strict_type => 1, store => \$vep},
	whole_gene => { default => 0,
			allow => [undef, 0, 1],
			strict_type => 1, store => \$whole_gene},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Genmod annotate
    my @commands = qw(genmod);  #Stores commands depending on input parameters

    ## Options
    if ($verbosity) {

	push(@commands, "-".$verbosity);
    }
    push(@commands, "models");

    if ($temp_directory_path) {

	push(@commands, "--temp_dir ".$temp_directory_path);
    }
    if ($family_file) {

	push(@commands, "--family_file ".$family_file);
    }
    if ($family_type) {

	push(@commands, "--family_type ".$family_type);
    }
    if ($reduced_penetrance_file_path) {

	push(@commands, "--reduced_penetrance ".$reduced_penetrance_file_path);
    }
    if ($thread_number) {

	push(@commands, "--processes ".$thread_number);
    }
    if ($vep) {

	push(@commands, "--vep ".$vep);
    }
    if ($whole_gene) {

	push(@commands, "--whole_gene");
    }

    if ($outfile_path) {

	push(@commands, "--outfile ".$outfile_path);  #Specify output filename
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


sub score {

##score

##Function : Perl wrapper for writing Genmod score recipe to $FILEHANDLE or return commands array. Based on genmod 3.7.0.
##Returns  : "@commands"
##Arguments: $infile_path, $family_file, $rank_model_file_path, $outfile_path, $stderrfile_path, $FILEHANDLE, $verbosity, $temp_directory_path, $rank_result, $family_type
##         : $infile_path          => Infile path to read from
##         : $family_file          => Family file
##         : $rank_model_file_path => The plug-in config file
##         : $outfile_path         => Outfile path to write to
##         : $stderrfile_path      => Stderr file path to write to {OPTIONAL}
##         : $FILEHANDLE           => Filehandle to write to
##         : $verbosity            => Increase output verbosity
##         : $temp_directory_path  => Directory for storing intermediate files
##         : $family_type          => Setup of family file
##         : $rank_result          => Add a info field that shows how the different categories contribute to the rank score

    my ($arg_href) = @_;

    ## Default(s)
    my $verbosity;
    my $rank_result;

    ## Flatten argument(s)
    my $infile_path;
    my $family_file;
    my $rank_model_file_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $temp_directory_path;
    my $family_type;

    my $tmpl = {
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	family_file => { required => 1, defined => 1, strict_type => 1, store => \$family_file},
	rank_model_file_path => { required => 1, defined => 1, strict_type => 1, store => \$rank_model_file_path},
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	FILEHANDLE => { store => \$FILEHANDLE},
	temp_directory_path => { strict_type => 1, store => \$temp_directory_path },
	family_type => { allow => ["ped", "alt", "cmms", "mip"],
			 strict_type => 1, store => \$family_type },
	verbosity => { allow => qr/^\w+$/,
		       strict_type => 1, store => \$verbosity},
	rank_result => { default => 0,
			 allow => [undef, 0, 1],
			 strict_type => 1, store => \$rank_result },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Genmod annotate
    my @commands = qw(genmod);  #Stores commands depending on input parameters

    ## Options
    if ($verbosity) {

	push(@commands, "-".$verbosity);
    }
    push(@commands, "score");

    if ($temp_directory_path) {

	push(@commands, "--temp_dir ".$temp_directory_path);
    }
    if ($family_file) {

	push(@commands, "--family_file ".$family_file);
    }
    if ($family_type) {

	push(@commands, "--family_type ".$family_type);
    }
    if ($rank_result) {

	push(@commands, "--rank_results");
    }
    if ($rank_model_file_path) {

	push(@commands, "--score_config ".$rank_model_file_path);
    }
    if ($outfile_path) {

	push(@commands, "--outfile ".$outfile_path);  #Specify output filename
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


sub compound {

##compound

##Function : Perl wrapper for writing Genmod compound recipe to $FILEHANDLE or return commands array. Based on genmod 3.7.0.
##Returns  : "@commands"
##Arguments: $infile_path, $outfile_path, $stderrfile_path, $FILEHANDLE, $temp_directory_path, $thread_number, $verbosity, $vep
##         : $infile_path         => Infile path to read from
##         : $outfile_path        => Outfile path to write to
##         : $stderrfile_path     => Stderr file path to write to {OPTIONAL}
##         : $FILEHANDLE          => Filehandle to write to
##         : $temp_directory_path => Directory for storing intermediate files
##         : $thread_number       => Define how many processes that should be use for annotation
##         : $verbosity           => Increase output verbosity
##         : $vep                 => If variants are annotated with the Variant Effect Predictor


    my ($arg_href) = @_;

    ## Default(s)
    my $thread_number;
    my $verbosity;
    my $vep;

    ## Flatten argument(s)
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $temp_directory_path;

    my $tmpl = {
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	FILEHANDLE => { store => \$FILEHANDLE},
	temp_directory_path => { strict_type => 1, store => \$temp_directory_path },
	thread_number => { default => 0,
			   allow => qr/^\d+$/,
			   strict_type => 1, store => \$thread_number},
	verbosity => { allow => qr/^\w+$/,
		       strict_type => 1, store => \$verbosity},
	vep => { default => 0,
		 allow => [undef, 0, 1],
		 strict_type => 1, store => \$vep},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Genmod compound
    my @commands = qw(genmod);  #Stores commands depending on input parameters

    ## Options
    if ($verbosity) {

	push(@commands, "-".$verbosity);
    }
    push(@commands, "compound");

    if ($temp_directory_path) {

	push(@commands, "--temp_dir ".$temp_directory_path);
    }
    if ($thread_number) {

	push(@commands, "--processes ".$thread_number);
    }
    if ($vep) {

	push(@commands, "--vep ".$vep);
    }
    if ($outfile_path) {

	push(@commands, "--outfile ".$outfile_path);  #Specify output filename
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

##Function : Perl wrapper for writing Genmod filter recipe to $FILEHANDLE or return commands array. Based on genmod 3.7.0.
##Returns  : "@commands"
##Arguments: $infile_path, $outfile_path, $stderrfile_path, $FILEHANDLE, $verbosity, $threshold
##         : $infile_path     => Infile path to read from
##         : $outfile_path    => Outfile path to write to
##         : $stderrfile_path => Stderr file path to write to {OPTIONAL}
##         : $FILEHANDLE      => Filehandle to write to
##         : $verbosity       => Increase output verbosity
##         : $threshold       => Directory for storing intermediate files

    my ($arg_href) = @_;

    ## Default(s)
    my $verbosity;

    ## Flatten argument(s)
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $threshold;

    my $tmpl = {
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	FILEHANDLE => { store => \$FILEHANDLE},
	threshold => { allow => qr/^\d+$|^\d+.\d+$/,
		       strict_type => 1, store => \$threshold },
	verbosity => { allow => qr/^\w+$/,
		       strict_type => 1, store => \$verbosity},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Genmod filter
    my @commands = qw(genmod);  #Stores commands depending on input parameters

    ## Options
    if ($verbosity) {

	push(@commands, "-".$verbosity);
    }
    push(@commands, "filter");

    if ($threshold) {

	push(@commands, "--threshold ".$threshold);
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


1;
