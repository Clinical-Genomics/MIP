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
    our @EXPORT_OK = qw(vcfparser);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub vcfparser {

##vcfparser

##Function : Perl wrapper for writing MIP vcfparser recipe to $FILEHANDLE or return commands array. Based on vcfparser 1.2.8.
##Returns  : "@commands"
##Arguments: $range_feature_annotation_columns_ref, $select_feature_annotation_columns_ref $infile_path, $outfile_path, $stderrfile_path, $range_feature_file_path, $select_feature_file_path, $select_feature_matching_column, $select_outfile, $FILEHANDLE, $parse_vep, $per_gene, $padding
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
##         : $parse_vep                             => Parse VEP transcript specific entries
##         : $per_gene                              => Output most severe consequence transcript
##         : $padding                               => Pad each gene with X number of nucleotides

    my ($arg_href) = @_;

    ## Default(s)
    my $parse_vep;
    my $per_gene;
    my $padding;

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
	select_feature_matching_column => { allow => qr/^\d+$/,
					    strict_type => 1, store => \$select_feature_matching_column },
	select_outfile => { strict_type => 1, store => \$select_outfile },
	FILEHANDLE => { store => \$FILEHANDLE },
	parse_vep => { default => 0,
		       allow => [0, 1],
		       strict_type => 1, store => \$parse_vep},
	per_gene => { default => 0,
		      allow => [0, 1],
		      strict_type => 1, store => \$per_gene},
	padding => { allow => [undef, qr/^\d+$/],
		     strict_type => 1, store => \$padding},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## svdb
    my @commands = qw(vcfparser);  #Stores commands depending on input parameters
    ## Infile
    push(@commands, $infile_path);

    ## Options
    if ($parse_vep) {

	push(@commands, "--parse_vep ".$parse_vep);
    }
    if ($per_gene) {

	push(@commands, "--per_gene ".$per_gene);
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

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


1;
