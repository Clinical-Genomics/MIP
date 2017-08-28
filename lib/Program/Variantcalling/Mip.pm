package Program::Variantcalling::Mip;

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
    our @EXPORT_OK = qw(vcfparser calculate_af max_af);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub vcfparser {

##vcfparser

##Function : Perl wrapper for writing MIP vcfparser recipe to $FILEHANDLE or return commands array. Based on vcfparser 1.2.10.
##Returns  : "@commands"
##Arguments: $range_feature_annotation_columns_ref, $select_feature_annotation_columns_ref $infile_path, $outfile_path, $stderrfile_path, $range_feature_file_path, $select_feature_file_path, $select_feature_matching_column, $select_outfile, $FILEHANDLE, $append_stderr_info, $parse_vep, $per_gene, $padding
##         : $range_feature_annotation_columns_ref  => Range file annotation columns
##         : $select_feature_annotation_columns_ref => Range file annotation columns
##         : $infile_path                           => Infile path
##         : $outfile_path                          => Outfile path
##         : $stderrfile_path                       => Stderr file path to write to {OPTIONAL}
##         : $range_feature_file_path               => Path to range feature file
##         : $select_feature_file_path              => Path to range feature file
##         : $select_feature_matching_column        => Select feature matching column
##         : $select_outfile                        => Select outfile
##         : $FILEHANDLE                            => Filehandle to write to
##         : $append_stderr_info                    => Append stderr info to file
##         : $parse_vep                             => Parse VEP transcript specific entries
##         : $per_gene                              => Output most severe consequence transcript
##         : $padding                               => Pad each gene with X number of nucleotides

    my ($arg_href) = @_;

    ## Default(s)
    my $parse_vep;
    my $per_gene;
    my $padding;
    my $append_stderr_info;

    ## Flatten argument(s)
    my $range_feature_annotation_columns_ref;
    my $select_feature_annotation_columns_ref;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $range_feature_file_path;
    my $select_feature_file_path;
    my $select_feature_matching_column;
    my $select_outfile;
    my $FILEHANDLE;

    my $tmpl = {
	range_feature_annotation_columns_ref => { default => [], strict_type => 1, store => \$range_feature_annotation_columns_ref },
	select_feature_annotation_columns_ref => { default => [], strict_type => 1, store => \$select_feature_annotation_columns_ref },
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path },
	outfile_path => { strict_type => 1, store => \$outfile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	range_feature_file_path => { strict_type => 1, store => \$range_feature_file_path },
	select_feature_file_path => { strict_type => 1, store => \$select_feature_file_path },
	select_feature_matching_column => { allow => [undef, qr/^\d+$/],
					    strict_type => 1, store => \$select_feature_matching_column },
	select_outfile => { strict_type => 1, store => \$select_outfile },
	FILEHANDLE => { store => \$FILEHANDLE },
	append_stderr_info => { default => 0,
				allow => [0, 1],
				strict_type => 1, store => \$append_stderr_info},
	parse_vep => { default => 0,
		       allow => [undef, 0, 1, 2],
		       strict_type => 1, store => \$parse_vep},
	per_gene => { default => 0,
		      allow => [undef, 0, 1],
		      strict_type => 1, store => \$per_gene},
	padding => { allow => [undef, qr/^\d+$/],
		     strict_type => 1, store => \$padding},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## vcfparser
    my @commands = qw(vcfparser);  #Stores commands depending on input parameters
    ## Infile
    push(@commands, $infile_path);

    ## Options
    if ($parse_vep) {

	push(@commands, "--parse_vep");
    }
    if ($per_gene) {

	push(@commands, "--per_gene");
    }
    if (defined($padding)) {

	push(@commands, "--padding ".$padding);
    }
    if ($range_feature_file_path) {

	push(@commands, "--range_feature_file ".$range_feature_file_path);
    }
    if(@$range_feature_annotation_columns_ref) {  #Limit output to regions

        push(@commands, "--range_feature_annotation_columns ".join(",", @{ $range_feature_annotation_columns_ref }));
    }
    if ($select_feature_file_path) {

	push(@commands, "--select_feature_file ".$select_feature_file_path);
    }
    if ($select_feature_matching_column) {

	push(@commands, "--select_feature_matching_column ".$select_feature_matching_column);
    }
    if(@$select_feature_annotation_columns_ref) {  #Limit output to regions

        push(@commands, "--select_feature_annotation_columns ".join(",", @{ $select_feature_annotation_columns_ref }));
    }
    if ($select_outfile) {

	push(@commands, "--select_outfile ".$select_outfile);
    }
    if ($outfile_path) {

	push(@commands, "> ".$outfile_path);  #Outfile prefix
    }
    if ($stderrfile_path) {
	
	if ($append_stderr_info) {

	    push(@commands, "2>> ".$stderrfile_path);  #Redirect and append stderr output to program specific stderr file
	}
	else {
	    
	    push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
	}
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


sub calculate_af {

##calculate_af

##Function : Perl wrapper for writing calculate_af recipe to already open $FILEHANDLE or return commands array. Based on calculate_af 0.0.2.
##Returns  : "@commands"
##Arguments: $FILEHANDLE, $infile_path, $outfile_path, $stderrfile_path, $append_stderr_info
##         : $FILEHANDLE         => Filehandle to write to
##         : $infile_path        => Infile path
##         : $outfile_path       => Outfile path
##         : $stderrfile_path    => Stderrfile path
##         : $append_stderr_info => Append stderr info to file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $append_stderr_info;

    my $tmpl = { 
	FILEHANDLE => { store => \$FILEHANDLE},
	infile_path => { strict_type => 1, store => \$infile_path},
	outfile_path => { strict_type => 1, store => \$outfile_path},
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	append_stderr_info => { default => 0,
				allow => [0, 1],
				strict_type => 1, store => \$append_stderr_info},
    };
    
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## calculate_af
    my @commands = qw(calculate_af);  #Stores commands depending on input parameters

    ## Options

    ## Infile
    if ($infile_path) {

	push(@commands, $infile_path);
    }
    ## Outfile
    if ($outfile_path) {

	push(@commands, "> ".$outfile_path);
    }
    if ($stderrfile_path) {

	if ($append_stderr_info) {

	    push(@commands, "2>> ".$stderrfile_path);  #Redirect and append stderr output to program specific stderr file
	}
	else {

	    push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
	}
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


sub max_af {

##max_af

##Function : Perl wrapper for writing max_af recipe to already open $FILEHANDLE or return commands array. Based on max_af 0.0.2.
##Returns  : "@commands"
##Arguments: $FILEHANDLE, $infile_path, $outfile_path, $stderrfile_path, $append_stderr_info
##         : $FILEHANDLE         => Filehandle to write to
##         : $infile_path        => Infile path
##         : $outfile_path       => Outfile path
##         : $stderrfile_path    => Stderrfile path
##         : $append_stderr_info => Append stderr info to file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $append_stderr_info;

    my $tmpl = { 
	FILEHANDLE => { store => \$FILEHANDLE},
	infile_path => { strict_type => 1, store => \$infile_path},
	outfile_path => { strict_type => 1, store => \$outfile_path},
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path},
	append_stderr_info => { default => 0,
				allow => [0, 1],
				strict_type => 1, store => \$append_stderr_info},
    };
    
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## max_af
    my @commands = qw(max_af);  #Stores commands depending on input parameters

    ## Options

    ## Infile
    if ($infile_path) {

	push(@commands, $infile_path);
    }
    ## Outfile
    if ($outfile_path) {

	push(@commands, "> ".$outfile_path);
    }
    if ($stderrfile_path) {

	if ($append_stderr_info) {

	    push(@commands, "2>> ".$stderrfile_path);  #Redirect and append stderr output to program specific stderr file
	}
	else {

	    push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
	}
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


1;
