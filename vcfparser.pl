#!/usr/bin/env perl

use Modern::Perl '2014';  #CPAN
use warnings qw( FATAL utf8 );
use autodie qw(open close :all);    #CPAN
use v5.18;  #Require at least perl 5.18
use utf8;  #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );

use Cwd;
use FindBin qw($Bin);  #Find directory of script
use File::Basename qw(basename);
use File::Spec::Functions qw(catdir catfile devnull);
use Getopt::Long;
use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case
use IO::File;

## Third party module(s)
use Set::IntervalTree; #CPAN

##MIPs lib/
use lib catdir($Bin, "lib");
use MIP_log::Log4perl qw(initiate_logger);
use Check::Check_modules qw(check_modules);
use Script::Utils qw(help);

our $USAGE;

BEGIN {

    my @modules = ("Modern::Perl",
		   "autodie",
		   "Set::IntervalTree",
		   "Log::Log4perl",
		   "MIP_log::Log4perl",
	);

    ## Evaluate that all modules required are installed
    Check::Check_modules::check_modules({modules_ref => \@modules,
					 program_name => $0,
					});

    $USAGE =
	basename($0).qq{ infile.vcf [OPTIONS] > outfile.vcf
           -pvep/--parse_vep Parse VEP transcript specific entries (Supply flag to enable)
           -rf/--range_feature_file (tsv)
           -rf_ac/--range_feature_annotation_columns
           -sf/--select_feature_file (tsv)
           -sf_mc/--select_feature_matching_column
           -sf_ac/--select_feature_annotation_columns
           -sof/--select_outfile (vcf)
           -pad/--padding (Default: "5000" nucleotides)
           -peg/--per_gene Output most severe consequence transcript (Supply flag to enable)
           -wst/--write_software_tag (Default: "1")
           -l/--log_file Log file (Default: "vcfparser.log")
           -h/--help Display this help message
           -v/--version Display version
        };
}

my ($infile, $select_feature_file, $select_feature_matching_column, $range_feature_file, $select_outfile);

##Scalar parameters with defaults
my ($write_software_tag, $padding, $log_file) = (1, 5000, 0, catfile(cwd(), "vcfparser.log"));

##Boolean
my ($parse_vep, $per_gene);

my (@range_feature_annotation_columns, @select_feature_annotation_columns);
my (%consequence_severity, %range_data, %select_data, %snpeff_cmd, %tree, %meta_data);

my $vcfparser_version = "1.2.10";

## Enables cmd "vcfparser.pl" to print usage help
if(!@ARGV) {

    help({USAGE => $USAGE,
	  exit_code => 0,
	 });
}
elsif ( (defined($ARGV)) && ($ARGV[0]!~/^-/) ) { #Collect potential infile - otherwise read from STDIN

    $infile = $ARGV[0];
}

###User Options
GetOptions('pvep|parse_vep' => \$parse_vep,
	   'rf|range_feature_file:s' => \$range_feature_file,
	   'rf_ac|range_feature_annotation_columns:s'  => \@range_feature_annotation_columns, #Comma separated list
	   'sf|select_feature_file:s' => \$select_feature_file,
	   'sf_mc|select_feature_matching_column:n' => \$select_feature_matching_column,
	   'sf_ac|select_feature_annotation_columns:s'  => \@select_feature_annotation_columns, #Comma separated list
	   'sof|select_outfile:s' => \$select_outfile,
	   'wst|write_software_tag:n' => \$write_software_tag,
	   'pad|padding:n' => \$padding,
	   'peg|per_gene' => \$per_gene,
	   'l|log_file:s' => \$log_file,
	   'h|help' => sub { say STDOUT $USAGE; exit;},  #Display help text
	   'v|version' => sub { say STDOUT "\n".basename($0)." ".$vcfparser_version, "\n"; exit;},  #Display version number
    )  or Script::Utils::help({USAGE => $USAGE,
			       exit_code => 1,
			      });

## Creates log object
my $log = MIP_log::Log4perl::initiate_logger({file_path_ref => \$log_file,
					      log_name => "Vcfparser",
					     });


## Basic flag option check
if ( (!@range_feature_annotation_columns) && ($range_feature_file) ) {

    $log->info($USAGE);
    $log->fatal("Need to specify which feature column(s) to use with range feature file: ".$range_feature_file." when annotating variants by using flag -rf_ac","\n");
    exit 1;
}
if ( (!$select_feature_matching_column) && ($select_feature_file) ) {

    $log->info($USAGE);
    $log->fatal("Need to specify which feature column to use with select feature file: ".$select_feature_file." when selecting variants by using flag -sf_mc","\n");
    exit 1;
}
if ( (!$select_outfile) && ($select_feature_file) ) {

    $log->info($USAGE);
    $log->fatal("Need to specify which a select outfile to use when selecting variants by using flag -sof","\n");
    exit 1;
}

@range_feature_annotation_columns = split(/,/,join(',',@range_feature_annotation_columns)); #Enables comma separated annotation columns on cmd
@select_feature_annotation_columns = split(/,/,join(',',@select_feature_annotation_columns)); #Enables comma separated annotation columns on cmd

###
#MAIN
###

define_select_data();

if ($range_feature_file) {

    read_feature_file({tree_href => \%tree,
		       feature_data_href => \%range_data,
		       feature_columns_ref => \@range_feature_annotation_columns,
		       range_file_key => "range_feature",
		       infile_path => $range_feature_file,
		       padding_ref => \$padding,
		      });
}

if ($select_feature_file) {

    read_feature_file({tree_href => \%tree,
		       feature_data_href => \%select_data,
		       feature_columns_ref => \@select_feature_annotation_columns,
		       range_file_key => "select_feature",
		       infile_path => $select_feature_file,
		       padding_ref => \$padding,
		       select_feature_matching_column => $select_feature_matching_column,
		      });
}

define_snpeff_annotations();

define_consequence_severity();

read_infile_vcf({meta_data_href => \%meta_data,
		 snpeff_cmd_href => \%snpeff_cmd,
		 range_data_href => \%range_data,
		 select_data_href => \%select_data,
		 consequence_severity_href => \%consequence_severity,
		 tree_href => \%tree,
		 range_feature_annotation_columns_ref => \@range_feature_annotation_columns,
		 select_feature_annotation_columns_ref => \@select_feature_annotation_columns,
		 select_outfile_path => $select_outfile,
		 vcfparser_version => $vcfparser_version,
		 select_feature_file => $select_feature_file,
		 per_gene => $per_gene,
		 parse_vep => $parse_vep,
		 write_software_tag => $write_software_tag,
		});

###
#Sub Routines
###

sub define_select_data {

##define_select_data

##Function : Defines arbitrary INFO fields based on headers in selectFile
##Returns  : ""
##Arguments: None

    $select_data{select_file}{HGNC_symbol}{info} = q?##INFO=<ID=HGNC_symbol,Number=.,Type=String,Description="The HGNC gene symbol">?;
    $select_data{select_file}{Ensembl_gene_id}{info} = q?##INFO=<ID=Ensembl_gene_id,Number=.,Type=String,Description="Ensembl gene identifier">?;
    $select_data{select_file}{OMIM_morbid}{info} = q?##INFO=<ID=OMIM_morbid,Number=.,Type=String,Description="OMIM morbid ID associated with gene(s)">?;
    $select_data{select_file}{Phenotypic_disease_model}{info} = q?##INFO=<ID=Phenotypic_disease_model,Number=.,Type=String,Description="Known disease gene(s) phenotype inheritance model">?;
    $select_data{select_file}{Clinical_db_gene_annotation}{info} = q?##INFO=<ID=Clinical_db_gene_annotation,Number=.,Type=String,Description="Gene disease group association">?;
    $select_data{select_file}{Reduced_penetrance}{info} = q?##INFO=<ID=Reduced_penetrance,Number=.,Type=String,Description="Pathogenic gene which can exhibit reduced penetrance">?;
    $select_data{select_file}{Disease_associated_transcript}{info} = q?##INFO=<ID=Disease_associated_transcript,Number=.,Type=String,Description="Known pathogenic transcript(s) for gene">?;
    $select_data{select_file}{Ensembl_transcript_to_refseq_transcript}{info} = q?##INFO=<ID=Ensembl_transcript_to_refseq_transcript,Number=.,Type=String,Description="The link between ensembl transcript and refSeq transcript IDs">?;
    $select_data{select_file}{Gene_description}{info} = q?##INFO=<ID=Gene_description,Number=.,Type=String,Description="The HGNC gene description">?;
    $select_data{select_file}{Genetic_disease_model}{info} = q?##INFO=<ID=Genetic_disease_model,Number=.,Type=String,Description="Known disease gene(s) inheritance model">?;
    $select_data{select_file}{No_hgnc_symbol}{info} = q?##INFO=<ID=No_hgnc_symbol,Number=.,Type=String,Description="Clinically relevant genetic regions lacking a HGNC_symbol or Ensembl gene ">?;
}

sub define_snpeff_annotations {

##define_snpeff_annotations

##Function : Defines the snpEff annotations that can be parsed and modified
##Returns  : ""
##Arguments: None

    $snpeff_cmd{snpeff}{phastCons100way_vertebrate_prediction_term}{File} = q?SnpSift dbnsfp?;
    $snpeff_cmd{snpeff}{phastCons100way_vertebrate_prediction_term}{vcf_key} = q?dbNSFP_phastCons100way_vertebrate?;
    $snpeff_cmd{snpeff}{phastCons100way_vertebrate_prediction_term}{info} = q?##INFO=<ID=phastCons100way_vertebrate_prediction_term,Number=A,Type=String,Description="PhastCons conservation prediction term">?;

    $snpeff_cmd{snpeff}{phyloP100way_vertebrate_prediction_term}{File} = q?SnpSift dbnsfp?;
    $snpeff_cmd{snpeff}{phyloP100way_vertebrate_prediction_term}{vcf_key} = q?dbNSFP_phyloP100way_vertebrate?;
    $snpeff_cmd{snpeff}{phyloP100way_vertebrate_prediction_term}{info} = q?##INFO=<ID=phyloP100way_vertebrate_prediction_term,Number=A,Type=String,Description="PhyloP conservation prediction term">?;

    $snpeff_cmd{snpeff}{'GERP++_RS_prediction_term'}{File} = q?SnpSift dbnsfp?;
    $snpeff_cmd{snpeff}{'GERP++_RS_prediction_term'}{vcf_key} = q?dbNSFP_GERP___RS?;
    $snpeff_cmd{snpeff}{'GERP++_RS_prediction_term'}{info} = q?##INFO=<ID=GERP++_RS_prediction_term,Number=A,Type=String,Description="GERP RS conservation prediction term">?;

}

sub define_consequence_severity {

##define_consequence_severity

##Function : Defines the precedence of consequences for SO-terms
##Returns  : ""
##Arguments: None

    $consequence_severity{transcript_ablation}{rank} = 1;
    $consequence_severity{transcript_ablation}{genetic_region_annotation} = "exonic";
    $consequence_severity{splice_donor_variant}{rank} = 2;
    $consequence_severity{splice_donor_variant}{genetic_region_annotation} = "splicing";
    $consequence_severity{splice_acceptor_variant}{rank} = 2;
    $consequence_severity{splice_acceptor_variant}{genetic_region_annotation} = "splicing";
    $consequence_severity{stop_gained}{rank} = 3;
    $consequence_severity{stop_gained}{genetic_region_annotation} = "exonic";
    $consequence_severity{frameshift_variant}{rank} = 4;
    $consequence_severity{frameshift_variant}{genetic_region_annotation} = "exonic";
    $consequence_severity{stop_lost}{rank} = 5;
    $consequence_severity{stop_lost}{genetic_region_annotation} = "exonic";
    $consequence_severity{start_lost}{rank} = 5;
    $consequence_severity{start_lost}{genetic_region_annotation} = "exonic";
    $consequence_severity{initiator_codon_variant}{rank} = 6;
    $consequence_severity{initiator_codon_variant}{genetic_region_annotation} = "exonic";
    $consequence_severity{inframe_insertion}{rank} = 6;
    $consequence_severity{inframe_insertion}{genetic_region_annotation} = "exonic";
    $consequence_severity{inframe_deletion}{rank} = 6;
    $consequence_severity{inframe_deletion}{genetic_region_annotation} = "exonic";
    $consequence_severity{missense_variant}{rank} = 6;
    $consequence_severity{missense_variant}{genetic_region_annotation} = "exonic";
    $consequence_severity{protein_altering_variant}{rank} = 6;
    $consequence_severity{protein_altering_variant}{genetic_region_annotation} = "exonic";
    $consequence_severity{transcript_amplification}{rank} = 7;
    $consequence_severity{transcript_amplification}{genetic_region_annotation} = "exonic";
    $consequence_severity{splice_region_variant}{rank} = 8;
    $consequence_severity{splice_region_variant}{genetic_region_annotation} = "splicing";
    $consequence_severity{incomplete_terminal_codon_variant}{rank} = 9;
    $consequence_severity{incomplete_terminal_codon_variant}{genetic_region_annotation} = "exonic";
    $consequence_severity{synonymous_variant}{rank} = 10;
    $consequence_severity{synonymous_variant}{genetic_region_annotation} = "exonic";
    $consequence_severity{stop_retained_variant}{rank} = 10;
    $consequence_severity{stop_retained_variant}{genetic_region_annotation} = "exonic";
    $consequence_severity{coding_sequence_variant}{rank} = 11;
    $consequence_severity{coding_sequence_variant}{genetic_region_annotation} = "exonic";
    $consequence_severity{mature_miRNA_variant}{rank} = 12;
    $consequence_severity{mature_miRNA_variant}{genetic_region_annotation} = "ncRNA_exonic";
    $consequence_severity{'5_prime_UTR_variant'}{rank} = 13;
    $consequence_severity{'5_prime_UTR_variant'}{genetic_region_annotation} = "5UTR";
    $consequence_severity{'3_prime_UTR_variant'}{rank} = 14;
    $consequence_severity{'3_prime_UTR_variant'}{genetic_region_annotation} = "3UTR";
    $consequence_severity{non_coding_transcript_exon_variant}{rank} = 15;
    $consequence_severity{non_coding_transcript_exon_variant}{genetic_region_annotation} = "ncRNA_exonic";
    $consequence_severity{non_coding_transcript_variant}{rank} = 15;
    $consequence_severity{non_coding_transcript_variant}{genetic_region_annotation} = "ncRNA";
    $consequence_severity{intron_variant}{rank} = 16;
    $consequence_severity{intron_variant}{genetic_region_annotation} = "intronic";
    $consequence_severity{NMD_transcript_variant}{rank} = 17;
    $consequence_severity{NMD_transcript_variant}{genetic_region_annotation} = "ncRNA";
    $consequence_severity{upstream_gene_variant}{rank} = 18;
    $consequence_severity{upstream_gene_variant}{genetic_region_annotation} = "upstream";
    $consequence_severity{downstream_gene_variant}{rank} = 19;
    $consequence_severity{downstream_gene_variant}{genetic_region_annotation} = "downstream";
    $consequence_severity{TFBS_ablation}{rank} = 20;
    $consequence_severity{TFBS_ablation}{genetic_region_annotation} = "TFBS";
    $consequence_severity{TFBS_amplification}{rank} = 21;
    $consequence_severity{TFBS_amplification}{genetic_region_annotation} = "TFBS";
    $consequence_severity{TF_binding_site_variant}{rank} = 22;
    $consequence_severity{TF_binding_site_variant}{genetic_region_annotation} = "TFBS";
    $consequence_severity{regulatory_region_variant}{rank} = 22;
    $consequence_severity{regulatory_region_variant}{genetic_region_annotation} = "regulatory_region";
    $consequence_severity{regulatory_region_ablation}{rank} = 23;
    $consequence_severity{regulatory_region_ablation}{genetic_region_annotation} = "regulatory_region";
    $consequence_severity{regulatory_region_amplification}{rank} = 24;
    $consequence_severity{regulatory_region_amplification}{genetic_region_annotation} = "regulatory_region";
    $consequence_severity{feature_elongation}{rank} = 25;
    $consequence_severity{feature_elongation}{genetic_region_annotation} = "genomic_feature";
    $consequence_severity{feature_truncation}{rank} = 26;
    $consequence_severity{feature_truncation}{genetic_region_annotation} = "genomic_feature";
    $consequence_severity{intergenic_variant}{rank} = 27;
    $consequence_severity{intergenic_variant}{genetic_region_annotation} = "intergenic"

}


sub read_feature_file {

##read_feature_file

##Function : Reads a file containg features to be annotated using range queries e.g. EnsemblGeneID. Adds to Metadata hash and creates Interval tree for feature.
##Returns  : ""
##Arguments: $tree_href, $feature_data_href, $feature_columns_ref, $range_file_key, $infile_path, $padding_ref, $select_feature_matching_column
##         : $tree_href                      => Interval tree hash {REF}
##         : $feature_data_href              => Range file hash {REF}
##         : $feature_columns_ref            => Range columns to include {REF}
##         : $range_file_key                 => Range file key used to seperate range file(s) i.e., select and range
##         : $infile_path                    => Infile path
##         : $padding_ref                    => Padding distance {REF}
##         : $select_feature_matching_column => Column in the select file to match with vcf key annotation {Optional}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $tree_href;
    my $feature_data_href;
    my $feature_columns_ref;
    my $range_file_key;
    my $infile_path;
    my $padding_ref;
    my $select_feature_matching_column;

    my $tmpl = {
	tree_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$tree_href},
	feature_data_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$feature_data_href},
	feature_columns_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$feature_columns_ref},
	range_file_key => { required => 1, defined => 1, strict_type => 1, store => \$range_file_key},
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	padding_ref => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$padding_ref},
	select_feature_matching_column => { strict_type => 1, store => \$select_feature_matching_column},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Retrieve logger object now that log_file has been set
    my $log = Log::Log4perl->get_logger("Vcfparser");

    my @headers; #Save headers from rangeFile

    my $FILEHANDLE = IO::Handle->new();  #Create anonymous filehandle
    open($FILEHANDLE, "<".$infile_path) or $log->logdie("Cannot open ".$infile_path.":".$!, "\n");

    while (<$FILEHANDLE>) {

	chomp $_; #Remove newline

	if (m/^\s+$/) {		# Avoid blank lines

	    next;
	}
	if ($_=~/^##/) {#MetaData - Avoid

	    next;
	}
	if ($_=~/^#/) {#Header/Comment

	    @headers = split(/\t/, $_);

	    for (my $extract_columns_counter=0;$extract_columns_counter<scalar(@$feature_columns_ref);$extract_columns_counter++) { #Defines what scalar to store

		my $header_ref = \$headers[ $$feature_columns_ref[$extract_columns_counter] ];  #Alias

		add_meta_data_info({meta_data_href => $feature_data_href,
				    range_file_key => $range_file_key,
				    header_ref => $header_ref,
				    position_ref => \$extract_columns_counter,
				    range_file_path_ref => \$infile_path,
				   });
	    }
	    next;
	}
	if ( $_ =~/^(\S+)/ ) {

	    my @line_elements = split("\t", $_); #Loads range file line elements

	    if (defined($select_feature_matching_column) ) {

		$line_elements[ $select_feature_matching_column ] =~ s/\s/_/g; # Replace whitespace with "_"
		$select_data{ $line_elements[ $select_feature_matching_column ] } = $line_elements[ $select_feature_matching_column ];
	    }

	    ## Create Interval Tree
	    if (@$feature_columns_ref) {#Annotate vcf with features from range file

		feature_annotations({tree_href => $tree_href,
				     feature_columns_ref => $feature_columns_ref,
				     line_elements_ref => \@line_elements,
				     range_file_key => $range_file_key,
				     padding_ref => $padding_ref,
				    });
	    }
	}
    }
    close($FILEHANDLE);
    $log->info("Finished reading ".$range_file_key." file: ".$infile_path);
}


sub read_infile_vcf {

##read_infile_vcf

##Function : Reads infile in vcf format and adds and parses annotations
##Returns  : ""
##Arguments: $meta_data_href, $snpeff_cmd_href, $range_data_href, $select_data_href, $consequence_severity_href, $tree_href, $range_feature_annotation_columns_ref, $select_feature_annotation_columns_ref, $select_outfile_path, $vcfparser_version, $select_feature_file, $parse_vep, $write_software_tag, $per_gene_ref
##         : $meta_data_href                        => Vcf meta data {REF}
##         : $snpeff_cmd_href                       => SnpEff meta data {REF}
##         : $range_data_href                       => Range file data {REF}
##         : $select_data_href                      => Select file data {REF}
##         : $consequence_severity_href             => Consequence severity for SO-terms {REF}
##         : $tree_href                             => Interval tree hash {REF}
##         : $range_feature_annotation_columns_ref  => Range feature columns {REF}
##         : $select_feature_annotation_columns_ref => Select feature columns {REF}
##         : $select_outfile_path                   => The select file path
##         : $vcfparser_version                     => vcfParser version
##         : $select_feature_file                   => The select feature file
##         : $parse_vep                             => Parse VEP output
##         : $write_software_tag                    => Write software tag to vcf header switch
##         : $per_gene                              => Only collect most severe transcript per gene

    my ($arg_href) = @_;

    ## Default(s)
    my $select_feature_file;
    my $parse_vep;
    my $write_software_tag;

    ## Flatten argument(s)
    my $meta_data_href;
    my $snpeff_cmd_href;
    my $range_data_href;
    my $select_data_href;
    my $consequence_severity_href;
    my $tree_href;
    my $range_feature_annotation_columns_ref;
    my $select_feature_annotation_columns_ref;
    my $select_outfile_path;
    my $vcfparser_version;
    my $per_gene;

    my $tmpl = {
	meta_data_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$meta_data_href},
	snpeff_cmd_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$snpeff_cmd_href},
	range_data_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$range_data_href},
	select_data_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$select_data_href},
	consequence_severity_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$consequence_severity_href},
	tree_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$tree_href},
	range_feature_annotation_columns_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$range_feature_annotation_columns_ref},
	select_feature_annotation_columns_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$select_feature_annotation_columns_ref},
	select_outfile_path => { strict_type => 1, store => \$select_outfile_path},
	vcfparser_version => { required => 1, defined => 1, strict_type => 1, store => \$vcfparser_version},
	select_feature_file => { default => 0,
				 strict_type => 1, store => \$select_feature_file},
	per_gene => { default => 0,
		      allow => [undef, 0, 1],
		      strict_type => 1, store => \$per_gene},
	parse_vep => { default => 0,
		       allow => [undef, 0, 1],
		       strict_type => 1, store => \$parse_vep},
	write_software_tag => { default => 1,
				allow => [0, 1],
				strict_type => 1, store => \$write_software_tag},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Retrieve logger object now that log_file has been set
    my $log = Log::Log4perl->get_logger("Vcfparser");

    my $FILEHANDLE = IO::Handle->new();  #Create anonymous filehandle for select file

    my @vep_format_fields;
    my %vep_format_field_column;
    my %vcf_header;

    my @vcf_format_columns;  #Catch #vcf header #CHROM line
    
    if ($select_feature_file) {
	
	open($FILEHANDLE, ">".$select_outfile_path) or $log->logdie("Cannot open ".$select_outfile_path.":".$!, "\n");
    }
    
    while (<>) {
	
	chomp $_;  # Remove newline
	
	if (m/^\s+$/) {	# Avoid blank lines
	    next;
	}
	if ($_=~/^##(\S+)=/) {  # MetaData
	    
	    parse_meta_data({meta_data_href => $meta_data_href,
			     meta_data_string => $_,
			    });

	    if ($_=~/INFO\=\<ID\=(\w+)/) { # Collect all INFO keys

		$vcf_header{info}{$1} = $1; #Save to hash
	    }
	    if ($_=~/SnpSiftCmd\=/) { #Find SnpEff command meta line
		
		for my $database (keys %{ $snpeff_cmd_href->{snpeff} }) {
		    
		    if ($_=~/$snpeff_cmd_href->{snpeff}{$database}{File}/) { #SnpEff/Sift has been used to annotate input vcf
			
			unless (defined($vcf_header{info}{$database})) { #Unless INFO header is already present add to meta_dataHeader
			    
			    $snpeff_cmd_href->{present}{database}{$database} = $database; #Save which frequency db has been used for later
			    push(@{ $meta_data_href->{info}{$database} }, $snpeff_cmd_href->{snpeff}{$database}{info});
			    
			    if (defined($snpeff_cmd_href->{snpeff}{$database}{fix_info})) { #If FIX_INFO flag is present add to meta_dataHeader
				
				push(@{ $meta_data_href->{fix_info}{$database} }, $snpeff_cmd_href->{snpeff}{$database}{fix_info});
			    }
			}
		    }
		}
		next;
	    }
	    if ($_=~/INFO\=\<ID\=CSQ/) { #Find VEP INFO Field
		
		if ($_=~/Format:\s(\S+)"\>/) { #Locate Format within VEP INFO meta line
		    
		    @vep_format_fields = split(/\|/, $1);
		    
		    while (my ($field_index, $field) = each (@vep_format_fields) ) {

			$vep_format_field_column{$field} = $field_index; #Save the order of VEP features
		    }
		}
		if ($parse_vep) {
		    
		    if ( ($vep_format_field_column{HGNC_ID}) && ($vep_format_field_column{Consequence}) ) {
			
			push(@{ $meta_data_href->{info}{most_severe_consequence} }, '##INFO=<ID=most_severe_consequence,Number=.,Type=String,Description="Most severe genomic consequence.">');
		    }
		}
		next;
	    }
	    next;
	}
	if ($_=~/^#CHROM/) {
	    
	    @vcf_format_columns = split(/\t/, $_);  #Split vcf format line
	    
	    add_feature_file_meta_data_to_vcf({meta_data_href => $meta_data_href,
					       vcf_header_href => \%vcf_header,
					       feature_annotation_columns_ref => $range_feature_annotation_columns_ref,
					       data_href => $range_data_href,
					       file_key => "Range",
					      });
	    add_feature_file_meta_data_to_vcf({meta_data_href => $meta_data_href,
					       vcf_header_href => \%vcf_header,
					       feature_annotation_columns_ref => $select_feature_annotation_columns_ref,
					       data_href => $select_data_href,
					       file_key => "Select",
					      });
	    
	    if ($write_software_tag) {
		
		add_program_to_meta_data_header({meta_data_href => $meta_data_href,
						 vcfparser_version => $vcfparser_version,
						});
	    }
	    if ($select_feature_file) { #SelectFile annotations
		
		write_meta_data({meta_data_href => $meta_data_href,
				 FILEHANDLE => *STDOUT,
				 SELECTFILEHANDLE => $FILEHANDLE,
				});
		say STDOUT $_;  #Write #CHROM header line
		say $FILEHANDLE $_;  #Write #CHROM header line
	    }
	    else {
		
		write_meta_data({meta_data_href => $meta_data_href,
				 FILEHANDLE => *STDOUT,
				});
		say STDOUT $_;  #Write #CHROM header line
	    }
	    next;
	}
	if ( $_ =~/^(\S+)/ ) {
	    
	    my %record;
	    my %consequence;
	    my %noid_region;
	    my $variant_line;
	    my $selected_variant_line;	    
	    
	    my @line_elements = split("\t", $_);  #Loads vcf file elements
	    
	    ##Check that we have an INFO field
	    unless ($line_elements[7]) {
		
		$log->fatal("No INFO field at line number: ".$.);
		$log->fatal("Displaying malformed line: ".$_);
		exit 1;
	    }
	    
	    ##Add line elements to record hash
	    while (my ($element_index, $element) = each (@line_elements) ) {

		$record{ $vcf_format_columns[$element_index] } = $element;  #Link vcf format headers to the line elements
	    }
	    
	    my @info_elements = split(/;/, $record{INFO});  #Add INFO elements
	    
	    ## Collect key value pairs in INFO field
	    foreach my $element (@info_elements) {
		
		my @key_value_pairs = split("=", $element);  #key index = 0 and value index = 1
		
		$record{INFO_key_value}{$key_value_pairs[0]} = $key_value_pairs[1];
	    }

	    for my $database (keys % {$snpeff_cmd_href->{present}{database}}) { #Note that the vcf should only contain 1 database entry
				    
		my $vcf_key = $snpeff_cmd_href->{snpeff}{$database}{vcf_key};

		if ($record{INFO_key_value}{$vcf_key}) {
			
		    my @allele_scores = split(",", $record{INFO_key_value}{$vcf_key}); #Split on ","
		    my $conservation_term;
		
		    if($database eq "phastCons100way_vertebrate_prediction_term") {
	
			$conservation_term = find_conserved({elements_ref => \@allele_scores,
							     score_cutoff => 0.8,
							    });
		    }
		    if($database eq "phyloP100way_vertebrate_prediction_term") {
			
			$conservation_term = find_conserved({elements_ref => \@allele_scores,
							     score_cutoff => 2.5,
							    });
		    }
		    if($database eq "GERP++_RS_prediction_term") {

			$conservation_term = find_conserved({elements_ref => \@allele_scores,
							     score_cutoff => 2,
							    });
		    }
		    
		    if (defined($conservation_term)) {
			
			## Save database info
			$record{INFO_addition}{$database} = $conservation_term;
		    }
		}
	    }
	    
	    ## Checks if an interval tree exists (per chr) and collects features from input array and adds annotations to line
	    %noid_region = tree_annotations({tree_href => $tree_href,
					     data_href => $select_data_href,
					     record_href => \%record,
					     line_elements_ref => \@line_elements,
					     range_file_key => "select_feature",
					    });  #noid_region is only for selectfile since all variants are passed to research file
	    
	    ## Checks if an interval tree exists (per chr) and collects features from input array and adds annotations to line
	    tree_annotations({tree_href => $tree_href,
			      data_href => $range_data_href,
			      record_href => \%record,
			      range_file_key => "range_feature",
			      line_elements_ref => \@line_elements,
			     });	   
	    
	    if ($parse_vep) {
		
		parse_vep_csq({record_href => \%record,
			       vep_format_field_column_href => \%vep_format_field_column,
			       select_data_href => $select_data_href,
			       consequence_href => \%consequence,
			       consequence_severity_href => $consequence_severity_href,
			       per_gene => $per_gene,
			      });
	    }
	    
	    for (my $line_elements_counter=0;$line_elements_counter<scalar(@line_elements);$line_elements_counter++) { #Add until INFO field
		
		if ($line_elements_counter < 7) { #Save fields until INFO field
		  
		    if ($record{select_transcripts}) {

			print $FILEHANDLE $record{ $vcf_format_columns[$line_elements_counter] }."\t";
		    }
		    print STDOUT $record{ $vcf_format_columns[$line_elements_counter] }."\t";
		}

		if ($line_elements_counter == 7) {

		    if ( ! $parse_vep) {

			if ($record{select_transcripts}) {
			    
			    print $FILEHANDLE $record{ $vcf_format_columns[$line_elements_counter] };
			}
			print STDOUT $record{ $vcf_format_columns[$line_elements_counter] };
		    }
		    else {

			my $counter = 0;
			foreach my $key (keys %{ $record{INFO_key_value} }) { 
			    
			    if (! $counter) {
				
				if (defined($record{INFO_key_value}{$key})) {
				    
				    if($key eq "CSQ") {
					
					if ($record{range_transcripts}) {
					    
					    print STDOUT $key."=".join(",", @{ $record{range_transcripts} });
					}
					if ($record{select_transcripts}) {
					    
					    print $FILEHANDLE $key."=".join(",", @{ $record{select_transcripts} });
					}
				    }
				    else {
					
					if ($record{select_transcripts}) {
					    
					    print $FILEHANDLE $key."=".$record{INFO_key_value}{$key};
					}
					print STDOUT $key."=".$record{INFO_key_value}{$key};
				    }
				}
				else {
				    
				    if ($record{select_transcripts}) {
					
					print $FILEHANDLE $key;
				    }
				    print STDOUT $key;
				}
			    }
			    else {
				
				if (defined($record{INFO_key_value}{$key})) {
				    
				    if($key eq "CSQ") {
					
					if ($record{range_transcripts}) {
					    
					    print STDOUT ";".$key."=".join(",", @{ $record{range_transcripts} });
					}
					if ($record{select_transcripts}) {
					    
					    print $FILEHANDLE ";".$key."=".join(",", @{ $record{select_transcripts} });
					}
				    }
				    else {
					
					if ($record{select_transcripts}) {
					    
					    print $FILEHANDLE ";".$key."=".$record{INFO_key_value}{$key};
					}
					print STDOUT ";".$key."=".$record{INFO_key_value}{$key};
				    }
				}
				else {
				    
				    if ($record{select_transcripts}) {			    
					
					print $FILEHANDLE ";".$key;
				    }
				    print STDOUT ";".$key;
				}
			    }
			    $counter++;
			}
		    }
		    
		    foreach my $key (keys %{ $record{INFO_addition} }) {
			
			if ($record{select_transcripts}) {
			    
			    print $FILEHANDLE ";".$key."=".$record{INFO_addition}{$key};
			}
			print STDOUT ";".$key."=".$record{INFO_addition}{$key};
		    }
		    if ($record{select_transcripts}) {
			
			foreach my $key (keys %{ $record{INFO_addition_select_feature} }) {
			    
			    print $FILEHANDLE ";".$key."=".$record{INFO_addition_select_feature}{$key};
			}
		    }
		    foreach my $key (keys %{ $record{INFO_addition_range_feature} }) {
			
			print STDOUT ";".$key."=".$record{INFO_addition_range_feature}{$key};
		    }
		    
		    if ($record{select_transcripts}) {

			print $FILEHANDLE "\t";
		    }
		    print STDOUT "\t";
		}
		if ($line_elements_counter > 7) {
		    
		    if ($record{select_transcripts}) {
			
			print $FILEHANDLE $record{ $vcf_format_columns[$line_elements_counter] }."\t";
		    }
		    print STDOUT $record{ $vcf_format_columns[$line_elements_counter] }."\t";
		}
	    }
	    if ($record{select_transcripts}) {
		
		print $FILEHANDLE "\n";
	    }
	    print STDOUT "\n";
	}
    }
    if ($select_feature_file) {
	
	close($FILEHANDLE);
    }
    $log->info("Finished Processing VCF\n");
}

sub parse_vep_csq {

##parse_vep_csq

##Function : 
##Returns  : ""
##Arguments: $record_href, $vep_format_field_column_href, select_data_href, consequence_href, consequence_severity_href, $per_gene
##         : $record_href                  => VCF record {REF}
##         : $vep_format_field_column_href => VEP format columns {REF}
##         : $select_data_href             => Select file data {REF}
##         : $consequence_href             => Variant consequence {REF}
##         : $consequence_severity_href    => Consequence severity for SO-terms {REF}
##         : $per_gene                     => Only collect most severe transcript per gene

    my ($arg_href) = @_;

    ## Default(s)

    ## Flatten argument(s)
    my $record_href;
    my $vep_format_field_column_href;
    my $select_data_href;
    my $consequence_href;
    my $consequence_severity_href;
    my $per_gene;

    my $tmpl = { 
	record_href  => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$record_href},
	vep_format_field_column_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$vep_format_field_column_href},
	select_data_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$select_data_href},
	consequence_href  => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$consequence_href},
	consequence_severity_href  => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$consequence_severity_href},
	per_gene => { default => 0,
		      allow => [undef, 0, 1],
		      strict_type => 1, store => \$per_gene},
    };
    
   
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    if($record_href->{INFO_key_value}{CSQ}) {

	my %most_severe_transcript;  #Collect most severe transcript per gene
	my @transcripts = split(/,/, $record_href->{INFO_key_value}{CSQ});  #Split into transcripts

	foreach my $transcript (@transcripts) {

	    my @transcript_effects = split(/\|/, $transcript); #Split in "|"

	    my $hgnc_id_column_ref = \$vep_format_field_column_href->{HGNC_ID};  #Alias hgnc_id column number

	    if ( (defined($transcript_effects[$$hgnc_id_column_ref]))
		 && $transcript_effects[$$hgnc_id_column_ref] ne "") {  #If gene
		
		my $hgnc_id_ref = \$transcript_effects[$$hgnc_id_column_ref];  #Alias HGNC Symbol
		my $transcript_id_ref = \$transcript_effects[ $vep_format_field_column_href->{Feature} ];  #Alias transcript_id
		my $allele_ref = \$transcript_effects[ $vep_format_field_column_href->{Allele} ];  #Alias allele

		my @consequences = split(/\&/, $transcript_effects[ $vep_format_field_column_href->{Consequence} ]); #Split consequence
		
		foreach my $consequence_term (@consequences) {
		    
		    check_terms({data_href => $consequence_severity_href,
				 term_ref => \$consequence_term,
				 data_category_name => "SO",
				});
		    
		    my $most_severe_consequence = $$hgnc_id_ref.":".$$allele_ref."|".$consequence_term;

		    if(defined($consequence_href->{$$hgnc_id_ref}{$$allele_ref}{score})) { #Compare to previous record
	
			if ($consequence_severity_href->{$consequence_term}{rank} < $consequence_href->{$$hgnc_id_ref}{$$allele_ref}{score}) { #Collect most severe consequence
			    
			    add_to_consequence_hash({key_ref => \$consequence_href->{$$hgnc_id_ref}{$$allele_ref}{score},
						     value_ref => \$consequence_severity_href->{$consequence_term}{rank},
						    });

			    add_to_consequence_hash({key_ref => \$consequence_href->{$$hgnc_id_ref}{$$allele_ref}{most_severe_consequence},
						     value_ref => \$most_severe_consequence,
						    });
			    add_to_consequence_hash({key_ref => \$consequence_href->{$$hgnc_id_ref}{$$allele_ref}{most_severe_transcript},
						     value_ref => \$transcript,
						    });
			}
		    }
		    else { #First pass
			
			add_to_consequence_hash({key_ref => \$consequence_href->{$$hgnc_id_ref}{$$allele_ref}{score},
						 value_ref => \$consequence_severity_href->{$consequence_term}{rank},
						});
			
			add_to_consequence_hash({key_ref => \$consequence_href->{$$hgnc_id_ref}{$$allele_ref}{most_severe_consequence},
						 value_ref => \$most_severe_consequence,
						});
			 add_to_consequence_hash({key_ref => \$consequence_href->{$$hgnc_id_ref}{$$allele_ref}{most_severe_transcript},
						  value_ref => \$transcript,
						 });
		    }
		}
		
		if(! $per_gene) {

		    if ($select_data_href->{$$hgnc_id_ref}) { #Exists in selected Features

			push(@{ $record_href->{select_transcripts} }, $transcript);  #Add all transcripts to selected transcripts
		    }
		    push(@{ $record_href->{range_transcripts} }, $transcript);  #Add all transcripts to range transcripts
		}
	    }
	    else {

	    ## Intergenic
	    push(@{ $record_href->{range_transcripts} }, $transcript);  #Add all transcripts to range transcripts
	    }
	}
	my @most_severe_select_consequences;
	my @most_severe_range_consequences;

      GENES:
	for my $gene (keys %$consequence_href) {
	    
	  ALLELS:
	    for my $allele (keys %{ $consequence_href->{$gene} }) {  #All alleles
		
		if ($select_data_href->{$gene}) { #Exists in selected Features
		    
		    push(@most_severe_select_consequences, $consequence_href->{$gene}{$allele}{most_severe_consequence});

		    if($per_gene) {

			push(@{ $record_href->{select_transcripts} }, $consequence_href->{$gene}{$allele}{most_severe_transcript});
		    }
		}
		push(@most_severe_range_consequences, $consequence_href->{$gene}{$allele}{most_severe_consequence});

		if($per_gene) {
		    
		    push(@{ $record_href->{range_transcripts} }, $consequence_href->{$gene}{$allele}{most_severe_transcript});
		}
	    }
	}
	if(@most_severe_select_consequences) {

	    $record_href->{INFO_addition_select_feature}{most_severe_consequence} = join(",", @most_severe_select_consequences);
	}
	if(@most_severe_range_consequences) {

	    $record_href->{INFO_addition_range_feature}{most_severe_consequence} = join(",", @most_severe_range_consequences);
	}
    }
}


sub add_to_consequence_hash {

##add_to_consequence_hash

##Function : Adds the most severe consequence or prediction to gene.
##Returns  : ""
##Arguments: $key_ref, $value_ref
##         : $key_ref   => The hash key to update {REF}
##         : $value_ref => The value to add to hash key

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $key_ref;
    my $value_ref;

    my $tmpl = {
	key_ref => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$key_ref},
	value_ref => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$value_ref},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    $$key_ref = $$value_ref;
}


sub add_field_to_element {

##add_field_to_element

##Function : Adds adds a field to an element.
##Returns  : ""
##Arguments: $selected_transcript_tracker_ref, $selected_line_ref, $line_ref, $value_ref, $separator
##         : $selected_transcript_tracker_ref => The selected transcript tracker for belonging to select file {REF}
##         : $selected_line_ref               => Selected line to add annotations to {REF}
##         : $line_ref                        => Variant line to add annotations to {REF}
##         : $value_ref                       => Field value {REF}
##         : $separator                       => Separator for field

    my ($arg_href) = @_;

    ## Default(s)
    my $separator;

    ## Flatten argument(s)
    my $selected_transcript_tracker_ref;
    my $selected_line_ref;
    my $line_ref;
    my $value_ref;

    my $tmpl = {
	selected_transcript_tracker_ref => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$selected_transcript_tracker_ref},
	selected_line_ref => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$selected_line_ref},
	line_ref => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$line_ref},
	separator => { default => ":", strict_type => 1, store => \$separator},
	value_ref => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$value_ref},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    if ($$selected_transcript_tracker_ref) {

	$$selected_line_ref .= $separator.$$value_ref;
    }

    ## Always include all transcripts in orphan list
    $$line_ref .= $separator.$$value_ref;
}


sub collect_consequence_genes {

##collect_consequence_genes

##Function : Collects all consequence and predictors per gene and adds info to line to be written.
##Returns  : ""
##Arguments: $consequence_href, $select_data_href, $fields_ref, $selected_variant_line_ref, $variantVariantLineRef
##         : $consequence_href          => Consequence(s) for each gene {REF}
##         : $select_data_href          => Select file data {REF}
##         : $fields_ref                => Features to be processed as determined by CSQ {REF}
##         : $selected_variant_line_ref => Selected line to add annotations to {REF}
##         : $variantVariantLineRef     => Variant line to add annotations to {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $consequence_href;
    my $select_data_href;
    my $fields_ref;
    my $selected_variant_line_ref;
    my $variant_line_ref;

    my $tmpl = {
	consequence_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$consequence_href},
	select_data_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$select_data_href},
	fields_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$fields_ref},
	selected_variant_line_ref => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$selected_variant_line_ref},
	variant_line_ref => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$variant_line_ref},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my %gene_counter;
    my %selected_gene_counter;
    my @temp_fields;
    my @selected_temp_fields;

    for (my $field_counter=0;$field_counter<scalar(@$fields_ref);$field_counter++) { #Set transcript counter to "0"

	$selected_gene_counter{$field_counter} = 0;
	$gene_counter{$field_counter} = 0;
    }

  GENES:
    for my $gene (keys %$consequence_href) {

      FIELDS:
	for (my $field_counter=0;$field_counter<scalar(@$fields_ref);$field_counter++) {

	    if ($select_data_href->{$gene}) { #Exists in selected Features

		collect_consequence_field({gene_counter_href => \%selected_gene_counter,
					   consequence_href => $consequence_href,
					   fields_ref => \@{$fields_ref},
					   selected_fields_ref => \@selected_temp_fields,
					   field_counter_ref => \$field_counter,
					   gene_ref => \$gene,
					  });
	    }

	    ## Always include all transcripts in research list
	    collect_consequence_field({gene_counter_href => \%gene_counter,
				       consequence_href => $consequence_href,
				       fields_ref => \@{$fields_ref},
				       selected_fields_ref => \@temp_fields,
				       field_counter_ref => \$field_counter,
				       gene_ref => \$gene,
				      });
	}
    }

    ## Adds to present line
    add_to_line({fields_ref => \@{$fields_ref},
		 temp_fields_ref => \@selected_temp_fields,
		 line_ref => \$$selected_variant_line_ref,
		});

    ## Adds to present line
    add_to_line({fields_ref => \@{$fields_ref},
		 temp_fields_ref => \@temp_fields,
		 line_ref => \$$variant_line_ref,
		});
}

sub collect_consequence_field {

##collect_consequence_field

##Function : Collects consequences for features in @feature_fields to temporary array for adding to line once all information are collected.
##Returns  : ""
##Arguments: $gene_counter_href, $consequence_href, $fields_ref, $selected_fields_ref, $field_counter_ref, $gene_ref
##         : $gene_counter_href   => Counts the number of transcripts per gene {REF}
##         : $consequence_href    => Consequence(s) for each gene {REF}
##         : $fields_ref          => Features to be processed as determined by CSQ {REF}
##         : $selected_fields_ref => Selected array {REF}
##         : $field_counter_ref   => Field number in feature {REF}
##         : $gene_ref            => The gene symbol {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $gene_counter_href;
    my $consequence_href;
    my $fields_ref;
    my $selected_fields_ref;
    my $field_counter_ref;
    my $gene_ref;

    my $tmpl = {
	gene_counter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$gene_counter_href},
	consequence_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$consequence_href},
	fields_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$fields_ref},
	selected_fields_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$selected_fields_ref},
	field_counter_ref => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$field_counter_ref},
	gene_ref => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$gene_ref},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $field_ref = \${$fields_ref}[$$field_counter_ref];  #Alias

    if (! $$gene_counter_href{$$field_counter_ref}) { #First time

	my $allel_counter = 0;

      ALLELS:
	for my $allele (keys %{ $consequence_href->{$$gene_ref} }) {  #All alleles

	    if (defined($$consequence_href{$$gene_ref}{$allele}{$$field_ref})
		&& ($$consequence_href{$$gene_ref}{$allele}{$$field_ref} ne "")) { #If feature exists - else do nothing

		if (! $allel_counter) {

		    $$selected_fields_ref[$$field_counter_ref] .= ";".$$field_ref."=".$$gene_ref.":".$allele."|".$$consequence_href{$$gene_ref}{$allele}{$$field_ref};
		}
		else {

		    $$selected_fields_ref[$$field_counter_ref] .= ",".$$gene_ref.":".$allele."|".$$consequence_href{$$gene_ref}{$allele}{$$field_ref};
		}
		$$gene_counter_href{$$field_counter_ref}++;
		$allel_counter++;
	    }
	}
    }
    else { #Subsequent passes

      ALLELS:
	for my $allele (keys %{ $consequence_href->{$$gene_ref} }) {  #All alleles

	    if (defined($$consequence_href{$$gene_ref}{$allele}{$$field_ref})
		&& ($$consequence_href{$$gene_ref}{$allele}{$$field_ref} ne "")) { #If feature exists - else do nothing

		$$selected_fields_ref[$$field_counter_ref] .= ",".$$gene_ref.":".$allele."|".$$consequence_href{$$gene_ref}{$allele}{$$field_ref};
	    }
	}
    }
}

sub add_to_line {

##add_to_line

##Function : Adds to present line.
##Returns  : ""
##Arguments: $fields_ref, $temp_fields_ref, $line_ref
##         : $fields_ref      => Features to be processed as determined by CSQ {REF}
##         : $temp_fields_ref => Annotations for feature {REF}
##         : $line_ref        => Line to add to {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $fields_ref;
    my $temp_fields_ref;
    my $line_ref;

    my $tmpl = {
	fields_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$fields_ref},
	temp_fields_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$temp_fields_ref},
	line_ref => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$line_ref},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    for (my $field_counter=0;$field_counter<scalar(@$fields_ref);$field_counter++) {

	if (defined($$temp_fields_ref[$field_counter])) {

	    $$line_ref .= $$temp_fields_ref[$field_counter];
	}
    }
}


sub convert_to_range {

##convert_to_range

##Function : Converts vcf sv to corresponding range coordinates.
##Returns  : "$final_start_position, $final_stop_position"
##Arguments: $fields_ref
##         : $fields_ref => Holds the chromosomal coordinates and allel data

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $fields_ref;

    my $tmpl = {
	fields_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$fields_ref},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $chromosome = $fields_ref->[0];
    my $start_position = $fields_ref->[1];
    my $reference_allele = $fields_ref->[3];
    my $alternative_allele = $fields_ref->[4];

    my $final_start_position = $start_position; #The most "uppstream" position per variant
    my $final_stop_position = 0; #The most "downstream" position per variant

    ## Convert to upper case
    ($reference_allele, $alternative_allele) = (uc $reference_allele, uc $alternative_allele);

    if ($alternative_allele eq ".") { #No Variant Call

	next;
    }
    my @alternative_alleles = split(/,/, $fields_ref->[4]);

    for (my $allel_counter=0;$allel_counter<scalar(@alternative_alleles);$allel_counter++) {

	my ($head, $newstart, $newend, $newref, $newalt);

	if (length ($reference_allele) == 1 and length ($alternative_alleles[$allel_counter]) == 1) { #SNV

	    ($newstart, $newend) = ($start_position, $start_position+length($reference_allele)-1);
	    ($newref, $newalt) = ($reference_allele, $alternative_alleles[$allel_counter]);
	}
	elsif (length($reference_allele) >= length($alternative_alleles[$allel_counter])) { #deletion or block substitution

	    $head = substr ($reference_allele, 0, length($alternative_alleles[$allel_counter]));

	    if ($head eq $alternative_alleles[$allel_counter]) {

		($newstart, $newend) = ($start_position+length($head), $start_position + length($reference_allele)-1);
		($newref, $newalt) = (substr($reference_allele, length($alternative_alleles[$allel_counter])), '-');
	    }
	    else {

		($newstart, $newend) = ($start_position, $start_position+length($reference_allele)-1);
		($newref, $newalt) = ($reference_allele, $alternative_alleles[$allel_counter]);
	    }
	}
	elsif (length($reference_allele) < length($alternative_alleles[$allel_counter])) { #insertion or block substitution

	    $head = substr($alternative_alleles[$allel_counter], 0, length($reference_allele));

	    if ($head eq $reference_allele) {

		($newstart, $newend) = ($start_position+length($reference_allele)-1, $start_position+length($reference_allele)-1);
		($newref, $newalt) = ('-', substr($alternative_alleles[$allel_counter], length($reference_allele)));
	    }
	    else {

		($newstart, $newend) = ($start_position, $start_position+length($reference_allele)-1);
		($newref, $newalt) = ($reference_allele, $alternative_alleles[$allel_counter]);
	    }
	}

	## Collect largest range per variant based on all alternative_alleles
	if ($final_start_position < $newstart) { #New start is upstream of old

	    $final_start_position = $newstart;
	}
	if ($final_stop_position < $newend) { #New end is downstream of old

	    $final_stop_position = $newend;
	}
    }
    return $final_start_position, $final_stop_position;
}


sub add_meta_data_info {

##add_meta_data_info

##Function : Adds arbitrary INFO fields to hash based on supplied headers unless header is already defined.
##Returns  : ""
##Arguments: $meta_data_href, $range_file_key, $header_ref, $position_ref, $range_file_path_ref
##         : $meta_data_href      => Hash to store meta_data in {REF}
##         : $range_file_key      => Range file key
##         : $header_ref          => Header from range file {REF}
##         : $position_ref        => Column position in supplied range file {REF}
##         : $range_file_path_ref => Range file path {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $meta_data_href;
    my $range_file_key;
    my $header_ref;
    my $position_ref;
    my $range_file_path_ref;

    my $tmpl = {
	meta_data_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$meta_data_href},
	range_file_key => { required => 1, defined => 1, strict_type => 1, store => \$range_file_key},
	header_ref => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$header_ref},
	position_ref => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$position_ref},
	range_file_path_ref => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$range_file_path_ref},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    if (defined($meta_data_href->{$range_file_key}{$$header_ref})) { #Add INFO from predefined entries

	$meta_data_href->{present}{$$header_ref}{info} = $meta_data_href->{$range_file_key}{$$header_ref}{info};
	$meta_data_href->{present}{$$header_ref}{column_order} = $$position_ref; #Column position in supplied range input file
    }
    else { #Add arbitrary INFO field using input header

	$meta_data_href->{present}{$$header_ref}{info} = q?##INFO=<ID=?.$$header_ref.q?,Number=.,Type=String,Description="String taken from ?.$$range_file_path_ref.q?">?;
	$meta_data_href->{present}{$$header_ref}{column_order} = $$position_ref; #Column position in supplied -sf_ac
    }
}

sub feature_annotations {

##feature_annotations

##Function : Creates the interval tree(s) for range and select feature files. Adds corresponding INFO fields to metadata.
##Returns  : ""
##Arguments: $tree_href, $feature_columns_ref, $line_elements_ref, $range_file_key, $padding_ref
##         : $tree_href           => Interval tree hash {REF}
##         : $feature_columns_ref => Feature annotations columns {REF}
##         : $line_elements_ref   => Feature file line {REF}
##         : $range_file_key      => Range file key
##         : $padding_ref         => The nucleotide distance to pad the range with {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s
    my $tree_href;
    my $feature_columns_ref;
    my $line_elements_ref;
    my $range_file_key;
    my $padding_ref;

    my $tmpl = {
	tree_href  => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$tree_href},
	feature_columns_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$feature_columns_ref},
	line_elements_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$line_elements_ref},
	range_file_key => { required => 1, defined => 1, strict_type => 1, store => \$range_file_key},
	padding_ref => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$padding_ref},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $features; #Features to collect (Format: ";" separated elements)

    for (my $extract_columns_counter=0;$extract_columns_counter<scalar(@$feature_columns_ref);$extract_columns_counter++) { #Defines what scalar to store

	my $feature_column_ref = \${$feature_columns_ref}[$extract_columns_counter];  #Alias

	if(defined($line_elements_ref->[ $$feature_column_ref ])) {

	    $line_elements_ref->[ $$feature_column_ref ] =~ tr/ /_/; #Remove whitespace since this is not allowed in vcf INFO field

	    if (!$extract_columns_counter) {

		$features .= $line_elements_ref->[ $$feature_column_ref ];
	    }
	    else {

		$features .= ";".$line_elements_ref->[ $$feature_column_ref ];
	    }
	}
    }
    unless(defined($tree_href->{$range_file_key}{ $line_elements_ref->[0] })) { #Only create once per firstKey

	$tree_href->{$range_file_key}{ $line_elements_ref->[0] } = Set::IntervalTree->new(); #Create tree
    }
    my $padded_start = $line_elements_ref->[1] - $$padding_ref;
    my $padded_stop = $line_elements_ref->[2] + $$padding_ref;
    $tree_href->{$range_file_key}{ $line_elements_ref->[0] }->insert($features, $padded_start, $padded_stop); #Store range and ";" sep string
}


sub tree_annotations {

##tree_annotations

##Function : Checks if an interval tree exists (per chr) and collects features from input array and adds annotations to line.
##Returns  : ""
##Arguments: $tree_href, $data_href, $record_href, $line_elements_ref, $range_file_key,
##         : $tree_href         => Interval tree hash {REF}
##         : $data_href         => Range file hash {REF}
##         : $record_href       => Record hash info {REF}
##         : $line_elements_ref => Infile vcf line elements array {REF}
##         : $range_file_key    => Range file key

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $tree_href;
    my $data_href;
    my $record_href;
    my $line_elements_ref;
    my $range_file_key;

    my $tmpl = {
	tree_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$tree_href},
	data_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$data_href},
	record_href  => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$record_href},
	line_elements_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$line_elements_ref},
	range_file_key => { required => 1, defined => 1, strict_type => 1, store => \$range_file_key},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my %noid_region;  #No HGNC_symbol or ensemblGeneID, but still clinically releveant e.g. mtD-loop

    if(defined($tree_href->{$range_file_key}{ $line_elements_ref->[0] }) ) { #Range annotations

	my $feature; #Features to be collected

	## #Convert SVs to range coordinates from vcf coordinates
	my ($start, $stop) = convert_to_range({fields_ref => $line_elements_ref,
					      });

	if ($start eq $stop) { #SNV

	    $feature = $tree_href->{$range_file_key}{ $line_elements_ref->[0] }->fetch($start, $stop+1); #Add 1 to SNV to create range input.
	}
	else {#Range input

	    $feature = $tree_href->{$range_file_key}{ $line_elements_ref->[0] }->fetch($start, $stop);
	}
	if (@$feature) { #Features found in tree

	    my %collected_annotation; #Collect all features before adding to line

	  FEATURE:
	    for (my $feature_counter=0;$feature_counter<scalar(@$feature);$feature_counter++) { #All features

		my @annotations = split(/;/, @$feature[$feature_counter]); #Split feature array ref into annotations

	      ANNOTATION:
		for (my $annotations_counter=0;$annotations_counter<scalar(@annotations);$annotations_counter++) { #All annotations

		    if( (defined($annotations[$annotations_counter])) && ($annotations[$annotations_counter] ne "") ) {

			push(@{ $collected_annotation{$annotations_counter} }, $annotations[$annotations_counter]);
		    }
		    if ($feature_counter == (scalar(@$feature - 1)) ) { #Last for this feature tuple

		      SELECTED_ANNOTATION:
			for my $range_annotation (keys % {$$data_href{present}}) { #All selected annotations

			    if ($$data_href{present}{$range_annotation}{column_order} eq $annotations_counter) { #Correct feature

				if ($range_annotation eq "Clinical_db_gene_annotation") {  #Special case, which is global and not gene centric

				    ## Collect unique elements from array reference and return array reference with unique elements
				    my $unique_ref = uniq_elements({elements_ref => \@{ $collected_annotation{$annotations_counter} },
								   });

				    @{ $collected_annotation{$annotations_counter} } = @{$unique_ref};
				}
				if ($range_annotation eq "No_hgnc_symbol") {  #Special case, where there is no HGNC or Ensembl gene ID but the region should be included in the select file anyway

				    my $id_key = join ("_", @$line_elements_ref[0..1, 3..4]);
				    $noid_region{$id_key}++;
				}
				if ( (defined($collected_annotation{$annotations_counter})) && (@{ $collected_annotation{$annotations_counter} }) ) {

				    $record_href->{"INFO_addition_".$range_file_key}{$range_annotation} = join(",", @{ $collected_annotation{$annotations_counter} });
				}
			    }
			}
		    }
		}
	    }
	}
    }
    return %noid_region;
}


sub add_program_to_meta_data_header {

##add_program_to_meta_data_header

##Function : Adds the program version and run date to the vcf meta-information section
##Returns  : ""
##Arguments: $meta_data_href, $vcfparser_version
##         : $meta_data_href    => Vcf meta data {REF}
##         : $vcfparser_version => vcfParser version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $meta_data_href;
    my $vcfparser_version;

    my $tmpl = {
	meta_data_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$meta_data_href},
	vcfparser_version => { required => 1, defined => 1, strict_type => 1, store => \$vcfparser_version},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my ($base, $script) = (`date +%Y%m%d`,`basename $0`);  #Catches current date and script name
    chomp($base,$script);  #Remove \n;
    push(@{ $meta_data_href->{software}{$script} }, "##Software=<ID=".$script.",Version=".$vcfparser_version.",Date=".$base);
}


sub find_af {

##find_af

##Function : Adds the least alternative allele(s) frequency to each line
##Returns  : ""
##Arguments: $elements_ref, $regexp
##         : $elements_ref => The INFO array {REF}
##         : $regexp       => The regexp to used to locate correct ID field

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $elements_ref;
    my $regexp;

    my $tmpl = {
	elements_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$elements_ref},
	regexp => { required => 1, defined => 1, strict_type => 1, store => \$regexp},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $temp_maf;

    for my $element (@$elements_ref) {

	if ($element =~/$regexp/) {  #Find the key=value field

	    my @value = split(/=/, $element);  #Split key=value pair

	    $temp_maf = $value[1];  #Collect whole string to represent all possible alleles
	}
    }
    return $temp_maf
}


sub find_conserved {

##find_conserved

##Function : Adds the least common alternative allele frequency to each line
##Returns  : "" or string prediction terms
##Arguments: $elements_ref, $score_cutoff
##         : $elements_ref => The INFO array {REF}
##         : $score_cutoff => Cut-off for conserved or not

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $elements_ref;
    my $score_cutoff;

    my $tmpl = {
	elements_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$elements_ref},
	score_cutoff => { required => 1, defined => 1, strict_type => 1, store => \$score_cutoff},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my @allel_values;

    for my $allel_value (@$elements_ref) {
	
	if ($allel_value >= $score_cutoff) {
	    
	    push(@allel_values, "Conserved");
	}
	else {
	    
	    push(@allel_values, "NotConserved");
	}
    }
    if (scalar(@allel_values) > 0) {
	
	return join(',', @allel_values);
    }
    else {
	
	return
    }
}


sub FindLCAF {

##FindLCAF

##Function : Adds the least common alternative allele frequency to each line
##Returns  : ""
##Arguments: $elements_ref, $regexp, $frequencyPosition
##         : $elements_ref      => The INFO array {REF}
##         : $regexp            => The regexp to used to locate correct ID field
##         : $frequencyPosition => The position to extract frequency from

    my ($arg_href) = @_;

    ## Default(s)
    my $frequencyPosition;

    ## Flatten argument(s)
    my $elements_ref;
    my $regexp;

    my $tmpl = {
	elements_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$elements_ref},
	regexp => { required => 1, defined => 1, strict_type => 1, store => \$regexp},
	frequencyPosition => { default => 0,
			       allow => qr/^\d+$/,
			       strict_type => 1, store => \$frequencyPosition},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my $temp_maf;

    for my $element (@$elements_ref) {

	if ($element =~/$regexp/) {  #Find the key=value field

	    my @value = split(/=/, $element);  #Split key=value pair

	    my @temp_mafs = sort {$a <=> $b} grep { $_ ne "." } split(",", $value[1]); #Split on ",", remove entries containing only "." and sort remaining entries numerically ascending order

	    if (scalar(@temp_mafs) > 0) {

		## We are interested in the least common allele listed for this position. We cannot connect the frequency position in the list and the multiple alternative alleles. So the best we can do is report the least common allele frequency for multiple alternative allels. Unless the least common frequency is lower than the frequency defined as pathogenic for rare disease (usually 0.01) then this will work. In that case this will be a false positive, but it is better than taking the actual MAF which would be a false negative if the pathogenic variant found in the patient(s) has a lower frequency than the MAF.
		$temp_maf = $temp_mafs[$frequencyPosition];
	    }
	}
    }
    return $temp_maf
}


sub parse_meta_data {

##parse_meta_data

##Function : Writes metadata to filehandle specified by order in meta_data_sections.
##Returns  : ""
##Arguments: $meta_data_href, $meta_data_string
##         : $meta_data_href   => Hash for meta_data {REF}
##         : $meta_data_string => The meta_data string from vcf header

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $meta_data_href;
    my $meta_data_string;

    my $tmpl = {
	meta_data_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$meta_data_href},
	meta_data_string => { required => 1, defined => 1, strict_type => 1, store => \$meta_data_string},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    if ($meta_data_string=~/^##fileformat/) {  #Catch fileformat as it has to be at the top of header

	push(@{ $meta_data_href->{fileformat}{fileformat} }, $meta_data_string);  #Save metadata string
    }
    elsif ($meta_data_string=~/^##contig/) {  #catch contigs to not sort them later

	push(@{ $meta_data_href->{contig}{contig} }, $meta_data_string);  #Save metadata string
    }
    elsif ($meta_data_string=~/^##(\w+)=(\S+)/) {  #FILTER, FORMAT, INFO etc and more custom records

	push(@{ $meta_data_href->{$1}{$2} }, $meta_data_string);  #Save metadata string
    }
    else {  #All oddities

	push(@{ $meta_data_href->{other}{other} }, $meta_data_string);  #Save metadata string
    }
}


sub write_meta_data {

##write_meta_data

##Function : Writes metadata to filehandle specified by order in meta_data_sections.
##Returns  : ""
##Arguments: $meta_data_href, $FILEHANDLE, $SELECTFILEHANDLE
##         : $meta_data_href   => Hash for meta_data {REF}
##         : $FILEHANDLE       => The filehandle to write to
##         : $SELECTFILEHANDLE => The filehandle to write to {Optional}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $meta_data_href;
    my $FILEHANDLE;
    my $SELECTFILEHANDLE;

    my $tmpl = {
	meta_data_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$meta_data_href},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
	SELECTFILEHANDLE => { store => \$SELECTFILEHANDLE},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my @meta_data_sections = ("fileformat", "FILTER", "FORMAT", "INFO", "FIX_INFO", "contig", "Software");  #Determine order to print for standard records
    my @lines;

    for (my $line_counter=0;$line_counter<scalar(@meta_data_sections);$line_counter++) {

	my $meta_data_record_ref = \$meta_data_sections[$line_counter];  #Alias

	if ($meta_data_href->{ $$meta_data_record_ref }) {  #MetaDataRecordExists

	    if ($$meta_data_record_ref eq "contig") {  #Should not be sorted, but printed "as is"

		foreach my $line (@{ $meta_data_href->{$$meta_data_record_ref}{ $$meta_data_record_ref } }) {

		    say $FILEHANDLE $line;

		    if (defined($SELECTFILEHANDLE)) {

			say $SELECTFILEHANDLE $line;
		    }
		}
		delete($meta_data_href->{ $$meta_data_record_ref });  #Enable print of rest later
	    }
	    else {

		foreach my $line (sort( keys %{ $meta_data_href->{ $$meta_data_record_ref } })) {

		    say $FILEHANDLE @{ $meta_data_href->{ $$meta_data_record_ref }{$line} };

		    if (defined($SELECTFILEHANDLE)) {

			say $SELECTFILEHANDLE @{ $meta_data_href->{ $$meta_data_record_ref }{$line} };
		    }
		}
		if (defined($SELECTFILEHANDLE)) {

		    foreach my $line (sort( keys %{ $meta_data_href->{select}{ $$meta_data_record_ref } })) {

			say $SELECTFILEHANDLE @{ $meta_data_href->{select}{ $$meta_data_record_ref }{$line} };
		    }
		}
		foreach my $line (sort( keys %{ $meta_data_href->{range}{ $$meta_data_record_ref } })) {

		    say $FILEHANDLE @{ $meta_data_href->{range}{ $$meta_data_record_ref }{$line} };
		}
		delete($meta_data_href->{ $$meta_data_record_ref });  #Enable print of rest later

		if ($meta_data_href->{select}{ $$meta_data_record_ref }) {

		    delete($meta_data_href->{select}{ $$meta_data_record_ref });
		}
		if ($meta_data_href->{range}{ $$meta_data_record_ref }) {

		    delete($meta_data_href->{range}{ $$meta_data_record_ref });
		}
	    }
	}
    }
    for my $keys (keys %$meta_data_href) {

	for my $second_key (keys %{ $meta_data_href->{$keys} }) {

	    if (ref($meta_data_href->{$keys}{$second_key}) eq "HASH") {

		for my $line ( sort(keys %{ $meta_data_href->{$keys}{$second_key} }) ) {

		    say $FILEHANDLE @{ $meta_data_href->{$keys}{$second_key}{$line} };

		    if (defined($SELECTFILEHANDLE)) {

			say $SELECTFILEHANDLE @{ $meta_data_href->{$keys}{$second_key}{$line} };
		    }
		    delete $meta_data_href->{$keys}{$second_key}{$line};
		}
	    }
	    elsif (ref($meta_data_href->{$keys}{$second_key}) eq "ARRAY") {

		foreach my $element (@{ $meta_data_href->{$keys}{$second_key} }) {

		    say $element;

		    if (defined($SELECTFILEHANDLE)) {

			say $SELECTFILEHANDLE $element;
		    }
		}
		delete $meta_data_href->{$keys}{$second_key};
	    }
	    else {

		say $FILEHANDLE @{ $meta_data_href->{$keys}{$second_key} };

		if (defined($SELECTFILEHANDLE)) {

		    say $SELECTFILEHANDLE @{ $meta_data_href->{$keys}{$second_key} };
		}
		delete $meta_data_href->{$keys}{$second_key};
	    }
	}
    }
}


sub uniq_elements {

##uniq_elements

##Function : Collect unique elements from array reference and return array reference with unique elements
##Returns  : "array reference"
##Arguments: $elements_ref
##         : $elements_ref => The array whose elements are to be made distinct {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $elements_ref;

    my $tmpl = {
	elements_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$elements_ref},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my %seen;

    return [ grep { !$seen{$_}++ } @$elements_ref ];  #For each element in array, see if seen before and only return list distinct elements
}

sub add_feature_file_meta_data_to_vcf {

##add_feature_file_meta_data_to_vcf

##Function : Adds feature file meta data annotation headers to meta data hash.
##Returns  : ""
##Arguments: $meta_data_href, $vcf_header_href, $data_href, $feature_annotation_columns_ref, $file_key
##         : $meta_data_href                 => Vcf meta data {REF}
##         : $vcf_header_href                => Vcf header meta data
##         : $data_href                      => Range file hash {REF}
##         : $feature_annotation_columns_ref => Range columns to include {REF}
##         : $file_key                       => Range file key used to seperate range file(s) i.e., select and range

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $meta_data_href;
    my $vcf_header_href;
    my $data_href;
    my $feature_annotation_columns_ref;
    my $file_key;

    my $tmpl = {
	meta_data_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$meta_data_href},
	vcf_header_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$vcf_header_href},
	data_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$data_href},
	feature_annotation_columns_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$feature_annotation_columns_ref},
	file_key => { required => 1, defined => 1, strict_type => 1, store => \$file_key},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    if (@$feature_annotation_columns_ref) { #Feature file annotations

	for my $annotation (keys %{ $data_href->{present} }) {

	    unless (defined($vcf_header_href->{info}{$annotation})) { #Unless INFO header is already present add to file

		push(@{ $meta_data_href->{ $file_key }{info}{$annotation} }, ${$data_href}{present}{$annotation}{info});  #Save specific featureFile INFO
	    }
	}
    }
}


sub check_terms {

##check_terms

##Function : Check that the found terms in the vcf corresond to known terms - otherwise croak and exit.
##Returns  : ""
##Arguments: $data_href, $term_ref, $data_category_name
##         : $data_href          => The term hash {REF}
##         : $term_ref           => The found term {REF}
##         : $data_category_name => The origin of the term i.e SO

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $data_href;
    my $term_ref;
    my $data_category_name;

    my $tmpl = {
	data_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$data_href},
	term_ref => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$term_ref},
	data_category_name => { required => 1, defined => 1, strict_type => 1, store => \$data_category_name},
    };

    if ( (defined($$term_ref)) && ($$term_ref ne "") ) {

	unless (exists($data_href->{ $$term_ref })) {

	    warn("Could not find ".$data_category_name." term from vcf in corresponding hash. Update hash to contain term: '".$$term_ref."'\n");
	    exit 1;
	}
    }
}
