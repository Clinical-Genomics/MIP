#!/usr/bin/env perl

use Modern::Perl '2014';
use warnings qw( FATAL utf8 );
use autodie;
use v5.18;  #Require at least perl 5.18
use utf8;  #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );

use Getopt::Long;
use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case
use IO::File;
use Set::IntervalTree; #CPAN

use vars qw($USAGE);

BEGIN {
    $USAGE =
	qq{vcfParser.pl infile.vcf > outfile.vcf
           -pVEP/--parseVEP Parse VEP transcript specific entries (Default: 0 (=no))
           -rf/--rangeFeatureFile (tsv)
           -rf_ac/--rangeFeatureAnnotationColumns
           -sf/--selectFeatureFile (tsv)
           -sf_mc/--selectFeatureMatchingColumn
           -sf_ac/--selectFeatureAnnotationColumns
           -sof/--selectOutfile selectOutfile (vcf)
           -pad/--padding (Default: "5000" nucleotides)
           -wst/--writeSoftwareTag (Default: "1")
           -h/--help Display this help message   
           -v/--version Display version
        };    
}

my ($infile, $parseVEP, $rangeFeatureFile, $selectFeatureFile, $selectFeatureMatchingColumn, $selectOutfile, $writeSoftwareTag, $padding) = ("", 0, "", "", "", "", 1, 5000);
my (@metaData, @selectMetaData, @rangeFeatureAnnotationColumns, @selectFeatureAnnotationColumns); 
my (%geneAnnotation, %consequenceSeverity, %rangeData, %selectData, %snpEffCmd, %tree, %metaData, %siftTerm, %polyPhenTerm);

my $vcfParserVersion = "1.2.8";

## Enables cmd "vcfParser.pl" to print usage help 
if(!@ARGV) {

    &Help({USAGE => $USAGE,
	   exitCode => 0,
	  });
}
elsif ( (defined($ARGV)) && ($ARGV[0]!~/^-/) ) { #Collect potential infile - otherwise read from STDIN
    
    $infile = $ARGV[0];
}

###User Options
GetOptions('pVEP|parseVEP:s' => \$parseVEP,
	   'rf|rangeFeatures:s' => \$rangeFeatureFile,
	   'rf_ac|rangeFeatureAnnotationColumns:s'  => \@rangeFeatureAnnotationColumns, #Comma separated list
	   'sf|selectFeatures:s' => \$selectFeatureFile,
	   'sf_mc|selectFeatureMatchingColumn:n' => \$selectFeatureMatchingColumn,
	   'sf_ac|selectFeatureAnnotationColumns:s'  => \@selectFeatureAnnotationColumns, #Comma separated list
	   'sof|selectOutfile:s' => \$selectOutfile,
	   'wst|writeSoftwareTag:n' => \$writeSoftwareTag,
	   'pad|padding:n' => \$padding,
	   'h|help' => sub { say STDOUT $USAGE; exit;},  #Display help text
	   'v|version' => sub { say STDOUT "\nvcfParser.pl ".$vcfParserVersion, "\n"; exit;},  #Display version number
    )  or &Help({USAGE => $USAGE,
		 exitCode => 1,
		});

## Basic flag option check
if ( (!@rangeFeatureAnnotationColumns) && ($rangeFeatureFile) ) {
    
    say STDOUT "\n".$USAGE;
    say STDERR "\n", "Need to specify which feature column(s) to use with range file: ".$rangeFeatureFile." when annotating variants by using flag -rf_ac","\n";
    exit;
}
if ( (!$selectFeatureMatchingColumn) && ($selectFeatureFile) ) {
    
    say STDOUT "\n".$USAGE;
    say STDERR "\n", "Need to specify which feature column to use with select file: ".$selectFeatureFile." when selecting variants by using flag -sf_mc","\n";
    exit;
}
if ( (!$selectOutfile) && ($selectFeatureFile) ) {
    
    say STDOUT "\n".$USAGE;
    say STDERR "\n", "Need to specify which a select outfile to use when selecting variants by using flag -sof","\n";
    exit;
}

@rangeFeatureAnnotationColumns = split(/,/,join(',',@rangeFeatureAnnotationColumns)); #Enables comma separated annotation columns on cmd
@selectFeatureAnnotationColumns = split(/,/,join(',',@selectFeatureAnnotationColumns)); #Enables comma separated annotation columns on cmd

###
#MAIN
###

&DefineSelectData();

if ($rangeFeatureFile) {

    &ReadFeatureFile({featureDataHashRef => \%rangeData,
		      featureColumnsArrayRef => \@rangeFeatureAnnotationColumns,
		      rangeFileKey => "RangeFile",
		      infilePath => $rangeFeatureFile,
		      paddingRef => \$padding,
		     });
}

if ($selectFeatureFile) {

    &ReadFeatureFile({featureDataHashRef => \%selectData,
		      featureColumnsArrayRef => \@selectFeatureAnnotationColumns,
		      rangeFileKey => "SelectFile",
		      infilePath => $selectFeatureFile,
		      paddingRef => \$padding,
		      selectFeatureMatchingColumn => $selectFeatureMatchingColumn,
		     });
}

&DefineSnpEffAnnotations();
&DefineConsequenceSeverity();
&DefineSiftTerms();
&DefinePolyPhenTerms();

&ReadInfileVCF({metaDataHashRef => \%metaData,
		snpEffCmdHashRef => \%snpEffCmd,
		rangeDataHashRef => \%rangeData,
		selectDataHashRef => \%selectData,
		consequenceSeverityHashRef => \%consequenceSeverity,
		siftTermHashRef => \%siftTerm,
		polyPhenTermHashRef => \%polyPhenTerm,
		rangeFeatureAnnotationColumnsArrayRef => \@rangeFeatureAnnotationColumns,
		selectFeatureAnnotationColumnsArrayRef => \@selectFeatureAnnotationColumns,
		selectOutFilePath => $selectOutfile,
		vcfParserVersion => $vcfParserVersion,
		selectFeatureFile => $selectFeatureFile,
		parseVEP => $parseVEP,
		writeSoftwareTag => $writeSoftwareTag,
	       });

###
#Sub Routines
###

sub DefineSelectData {

##DefineSelectData

##Function : Defines arbitrary INFO fields based on headers in selectFile
##Returns  : ""
##Arguments: None

    $selectData{SelectFile}{HGNC_symbol}{INFO} = q?##INFO=<ID=HGNC_symbol,Number=.,Type=String,Description="The HGNC gene symbol">?;
    $selectData{SelectFile}{Ensembl_gene_id}{INFO} = q?##INFO=<ID=Ensembl_gene_id,Number=.,Type=String,Description="Ensembl gene identifier">?;
    $selectData{SelectFile}{OMIM_morbid}{INFO} = q?##INFO=<ID=OMIM_morbid,Number=.,Type=String,Description="OMIM morbid ID associated with gene(s)">?;
    $selectData{SelectFile}{Phenotypic_disease_model}{INFO} = q?##INFO=<ID=Phenotypic_disease_model,Number=.,Type=String,Description="Known disease gene(s) phenotype inheritance model">?;
    $selectData{SelectFile}{Clinical_db_gene_annotation}{INFO} = q?##INFO=<ID=Clinical_db_gene_annotation,Number=.,Type=String,Description="Gene disease group association">?;
    $selectData{SelectFile}{Reduced_penetrance}{INFO} = q?##INFO=<ID=Reduced_penetrance,Number=.,Type=String,Description="Pathogenic gene which can exhibit reduced penetrance">?;
    $selectData{SelectFile}{Disease_associated_transcript}{INFO} = q?##INFO=<ID=Disease_associated_transcript,Number=.,Type=String,Description="Known pathogenic transcript(s) for gene">?;
    $selectData{SelectFile}{Ensembl_transcript_to_refseq_transcript}{INFO} = q?##INFO=<ID=Ensembl_transcript_to_refseq_transcript,Number=.,Type=String,Description="The link between ensembl transcript and refSeq transcript IDs">?;
    $selectData{SelectFile}{Gene_description}{INFO} = q?##INFO=<ID=Gene_description,Number=.,Type=String,Description="The HGNC gene description">?;
    $selectData{SelectFile}{Genetic_disease_model}{INFO} = q?##INFO=<ID=Genetic_disease_model,Number=.,Type=String,Description="Known disease gene(s) inheritance model">?;
    $selectData{SelectFile}{No_hgnc_symbol}{INFO} = q?##INFO=<ID=No_hgnc_symbol,Number=.,Type=String,Description="Clinically relevant genetic regions lacking a HGNC_symbol or Ensembl gene ">?;
}

sub DefineSnpEffAnnotations {

##DefineSnpEffAnnotations

##Function : Defines the snpEff annotations that can be parsed and modified 
##Returns  : ""
##Arguments: None
    
    $snpEffCmd{SnpEff}{Dbsnp129LCAF}{File} = q?dbsnp_\S+.excluding_sites_after_129.vcf?;
    $snpEffCmd{SnpEff}{Dbsnp129LCAF}{INFO} = q?##INFO=<ID=Dbsnp129LCAF,Number=1,Type=Float,Description="Least common AF in dbSNP excluding sites after 129.">?;
    $snpEffCmd{SnpEff}{Dbsnp129LCAF}{FIX_INFO} = q?##INFO=<ID=SnpSift_CAF,Number=.,Type=String,Description="An ordered, comma delimited list of allele frequencies based on 1000Genomes, starting with the reference allele followed by alternate alleles as ordered in the ALT column. Where a 1000Genomes alternate allele is not in the dbSNPs alternate allele set, the allele is added to the ALT column.  The minor allele is the second largest value in the list, and was previuosly reported in VCF as the GMAF.  This is the GMAF reported on the RefSNP and EntrezSNP pages and VariationReporter">?;

    $snpEffCmd{SnpEff}{DbsnpLCAF}{File} = q?dbsnp_\d+.\w\d+.vcf?;
    $snpEffCmd{SnpEff}{DbsnpLCAF}{INFO} = q?##INFO=<ID=DbsnpLCAF,Number=1,Type=Float,Description="Least common AF in the DbSNP database.">?;
    $snpEffCmd{SnpEff}{DbsnpLCAF}{FIX_INFO} = q?##INFO=<ID=SnpSift_CAF,Number=.,Type=String,Description="An ordered, comma delimited list of allele frequencies based on 1000Genomes, starting with the reference allele followed by alternate alleles as ordered in the ALT column. Where a 1000Genomes alternate allele is not in the dbSNPs alternate allele set, the allele is added to the ALT column.  The minor allele is the second largest value in the list, and was previuosly reported in VCF as the GMAF.  This is the GMAF reported on the RefSNP and EntrezSNP pages and VariationReporter">?;

    $snpEffCmd{SnpEff}{ESPMAF}{File} = q?ESP\d+SI-V\d+-\w+.updatedProteinHgvs.snps_indels.vcf?;
    $snpEffCmd{SnpEff}{ESPMAF}{INFO} = q?##INFO=<ID=ESPMAF,Number=1,Type=Float,Description="Global Minor Allele Frequency in the ESP database.">?;
    $snpEffCmd{SnpEff}{ESPMAF}{FIX_INFO} = q?##INFO=<ID=SnpSift_MAF,Number=.,Type=String,Description="Minor Allele Frequency in percent in the order of EA,AA,All">?;

    $snpEffCmd{SnpEff}{EXACMAXAF}{File} = q?ExAC.r\d+.\d+.sites.vep.vcf?;
    $snpEffCmd{SnpEff}{EXACMAXAF}{INFO} = q?##INFO=<ID=EXACMAXAF,Number=A,Type=Float,Description="Estimated max allele frequency in the range (0,1) in Exac">?;
    $snpEffCmd{SnpEff}{EXACMAXAF}{FIX_INFO} = q?##INFO=<ID=SnpSift_AF,Number=.,Type=String,Description="Estimated max allele frequency in the range (0,1)">?;

    $snpEffCmd{SnpEff}{phastCons100way_vertebrate_prediction_term}{File} = q?SnpSift dbnsfp?;
    $snpEffCmd{SnpEff}{phastCons100way_vertebrate_prediction_term}{INFO} = q?##INFO=<ID=phastCons100way_vertebrate_prediction_term,Number=A,Type=String,Description="PhastCons conservation prediction term">?;

    $snpEffCmd{SnpEff}{phyloP100way_vertebrate_prediction_term}{File} = q?SnpSift dbnsfp?;
    $snpEffCmd{SnpEff}{phyloP100way_vertebrate_prediction_term}{INFO} = q?##INFO=<ID=phyloP100way_vertebrate_prediction_term,Number=A,Type=String,Description="PhyloP conservation prediction term">?;

    $snpEffCmd{SnpEff}{'GERP++_RS_prediction_term'}{File} = q?SnpSift dbnsfp?;
    $snpEffCmd{SnpEff}{'GERP++_RS_prediction_term'}{INFO} = q?##INFO=<ID=GERP++_RS_prediction_term,Number=A,Type=String,Description="GERP RS conservation prediction term">?;

}

sub DefineConsequenceSeverity {

##DefineConsequenceSeverity

##Function : Defines the precedence of consequences for SO-terms
##Returns  : ""
##Arguments: None

    $consequenceSeverity{transcript_ablation}{Rank} = 1;
    $consequenceSeverity{transcript_ablation}{GeneticRegionAnnotation} = "exonic";
    $consequenceSeverity{splice_donor_variant}{Rank} = 2;
    $consequenceSeverity{splice_donor_variant}{GeneticRegionAnnotation} = "splicing";
    $consequenceSeverity{splice_acceptor_variant}{Rank} = 2;
    $consequenceSeverity{splice_acceptor_variant}{GeneticRegionAnnotation} = "splicing";
    $consequenceSeverity{stop_gained}{Rank} = 3;
    $consequenceSeverity{stop_gained}{GeneticRegionAnnotation} = "exonic";
    $consequenceSeverity{frameshift_variant}{Rank} = 4;
    $consequenceSeverity{frameshift_variant}{GeneticRegionAnnotation} = "exonic";
    $consequenceSeverity{stop_lost}{Rank} = 5;
    $consequenceSeverity{stop_lost}{GeneticRegionAnnotation} = "exonic";
    $consequenceSeverity{start_lost}{Rank} = 5;
    $consequenceSeverity{start_lost}{GeneticRegionAnnotation} = "exonic";
    $consequenceSeverity{initiator_codon_variant}{Rank} = 6;
    $consequenceSeverity{initiator_codon_variant}{GeneticRegionAnnotation} = "exonic";
    $consequenceSeverity{inframe_insertion}{Rank} = 6;
    $consequenceSeverity{inframe_insertion}{GeneticRegionAnnotation} = "exonic";
    $consequenceSeverity{inframe_deletion}{Rank} = 6;
    $consequenceSeverity{inframe_deletion}{GeneticRegionAnnotation} = "exonic";
    $consequenceSeverity{missense_variant}{Rank} = 6;
    $consequenceSeverity{missense_variant}{GeneticRegionAnnotation} = "exonic";
    $consequenceSeverity{protein_altering_variant}{Rank} = 6;
    $consequenceSeverity{protein_altering_variant}{GeneticRegionAnnotation} = "exonic";
    $consequenceSeverity{transcript_amplification}{Rank} = 7;
    $consequenceSeverity{transcript_amplification}{GeneticRegionAnnotation} = "exonic";
    $consequenceSeverity{splice_region_variant}{Rank} = 8;
    $consequenceSeverity{splice_region_variant}{GeneticRegionAnnotation} = "splicing";
    $consequenceSeverity{incomplete_terminal_codon_variant}{Rank} = 9;
    $consequenceSeverity{incomplete_terminal_codon_variant}{GeneticRegionAnnotation} = "exonic";
    $consequenceSeverity{synonymous_variant}{Rank} = 10;
    $consequenceSeverity{synonymous_variant}{GeneticRegionAnnotation} = "exonic";
    $consequenceSeverity{stop_retained_variant}{Rank} = 10;
    $consequenceSeverity{stop_retained_variant}{GeneticRegionAnnotation} = "exonic";
    $consequenceSeverity{coding_sequence_variant}{Rank} = 11;
    $consequenceSeverity{coding_sequence_variant}{GeneticRegionAnnotation} = "exonic";
    $consequenceSeverity{mature_miRNA_variant}{Rank} = 12;
    $consequenceSeverity{mature_miRNA_variant}{GeneticRegionAnnotation} = "ncRNA_exonic";
    $consequenceSeverity{'5_prime_UTR_variant'}{Rank} = 13;
    $consequenceSeverity{'5_prime_UTR_variant'}{GeneticRegionAnnotation} = "5UTR";
    $consequenceSeverity{'3_prime_UTR_variant'}{Rank} = 14;
    $consequenceSeverity{'3_prime_UTR_variant'}{GeneticRegionAnnotation} = "3UTR";
    $consequenceSeverity{non_coding_transcript_exon_variant}{Rank} = 15;
    $consequenceSeverity{non_coding_transcript_exon_variant}{GeneticRegionAnnotation} = "ncRNA_exonic";
    $consequenceSeverity{non_coding_transcript_variant}{Rank} = 15;
    $consequenceSeverity{non_coding_transcript_variant}{GeneticRegionAnnotation} = "ncRNA";
    $consequenceSeverity{intron_variant}{Rank} = 16;
    $consequenceSeverity{intron_variant}{GeneticRegionAnnotation} = "intronic";
    $consequenceSeverity{NMD_transcript_variant}{Rank} = 17;
    $consequenceSeverity{NMD_transcript_variant}{GeneticRegionAnnotation} = "ncRNA";
    $consequenceSeverity{upstream_gene_variant}{Rank} = 18;
    $consequenceSeverity{upstream_gene_variant}{GeneticRegionAnnotation} = "upstream";
    $consequenceSeverity{downstream_gene_variant}{Rank} = 19;
    $consequenceSeverity{downstream_gene_variant}{GeneticRegionAnnotation} = "downstream";
    $consequenceSeverity{TFBS_ablation}{Rank} = 20;
    $consequenceSeverity{TFBS_ablation}{GeneticRegionAnnotation} = "TFBS";
    $consequenceSeverity{TFBS_amplification}{Rank} = 21;
    $consequenceSeverity{TFBS_amplification}{GeneticRegionAnnotation} = "TFBS";
    $consequenceSeverity{TF_binding_site_variant}{Rank} = 22;
    $consequenceSeverity{TF_binding_site_variant}{GeneticRegionAnnotation} = "TFBS";
    $consequenceSeverity{regulatory_region_variant}{Rank} = 22;
    $consequenceSeverity{regulatory_region_variant}{GeneticRegionAnnotation} = "regulatory_region";
    $consequenceSeverity{regulatory_region_ablation}{Rank} = 23;
    $consequenceSeverity{regulatory_region_ablation}{GeneticRegionAnnotation} = "regulatory_region";
    $consequenceSeverity{regulatory_region_amplification}{Rank} = 24;
    $consequenceSeverity{regulatory_region_amplification}{GeneticRegionAnnotation} = "regulatory_region";
    $consequenceSeverity{feature_elongation}{Rank} = 25;
    $consequenceSeverity{feature_elongation}{GeneticRegionAnnotation} = "genomic_feature";
    $consequenceSeverity{feature_truncation}{Rank} = 26;
    $consequenceSeverity{feature_truncation}{GeneticRegionAnnotation} = "genomic_feature";
    $consequenceSeverity{intergenic_variant}{Rank} = 27;
    $consequenceSeverity{intergenic_variant}{GeneticRegionAnnotation} = "intergenic" 

}


sub DefineSiftTerms {

##DefineSiftTerms

##Function : Defines the allowed sift terms
##Returns  : ""
##Arguments: None

    $siftTerm{deleterious} = "";
    $siftTerm{tolerated} = "";
    $siftTerm{tolerated_low_confidence} = "";
    $siftTerm{deleterious_low_confidence} = "";

}

sub DefinePolyPhenTerms {

##DefinePolyPhenTerms

##Function : Defines the allowed polyPhen terms
##Returns  : ""
##Arguments: None

    $polyPhenTerm{probably_damaging} = "";
    $polyPhenTerm{possibly_damaging} = "";
    $polyPhenTerm{benign} = "";
    $polyPhenTerm{unknown} = "";

}

sub ReadFeatureFile {

##ReadFeatureFile

##Function : Reads a file containg features to be annotated using range queries e.g. EnsemblGeneID. Adds to Metadata hash and creates Interval tree for feature.
##Returns  : ""
##Arguments: $featureDataHashRef, $featureColumnsArrayRef, $rangeFileKey, $infilePath, $paddingRef, $selectFeatureMatchingColumn
##         : $featureDataHashRef          => Range file hash {REF}
##         : $featureColumnsArrayRef      => Range columns to include {REF}
##         : $rangeFileKey                => Range file key used to seperate range file(s) i.e., select and range
##         : $infilePath                  => Infile path
##         : $paddingRef                  => Padding distance {REF}
##         : $selectFeatureMatchingColumn => Column in the select file to match with vcf key annotation {Optional}

    my ($argHashRef) = @_;
    
    ## Flatten argument(s)
    my $featureDataHashRef;
    my $featureColumnsArrayRef;
    my $rangeFileKey;
    my $infilePath;
    my $paddingRef;
    my $selectFeatureMatchingColumn;

    my $tmpl = { 
	featureDataHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$featureDataHashRef},
	featureColumnsArrayRef => { required => 1, defined => 1, default => [], strict_type => 1, store => \$featureColumnsArrayRef},
	rangeFileKey => { required => 1, defined => 1, strict_type => 1, store => \$rangeFileKey},
	infilePath => { required => 1, defined => 1, strict_type => 1, store => \$infilePath},
	paddingRef => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$paddingRef},
	selectFeatureMatchingColumn => { strict_type => 1, store => \$selectFeatureMatchingColumn},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my @headers; #Save headers from rangeFile

    open(RRF, "<".$infilePath) or die "Cannot open ".$infilePath.":".$!, "\n"; 

    while (<RRF>) {
	
	chomp $_; #Remove newline

	if (m/^\s+$/) {		# Avoid blank lines
	    
	    next;
	}
	if ($_=~/^##/) {#MetaData - Avoid

	    next;
	}
	if ($_=~/^#/) {#Header/Comment
	    
	    @headers = split(/\t/, $_);

	    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar(@{$featureColumnsArrayRef});$extractColumnsCounter++) { #Defines what scalar to store
	
		my $headerRef = \$headers[ $$featureColumnsArrayRef[$extractColumnsCounter] ];  #Alias
		
		&AddMetaDataINFO({metaDataHashRef => $featureDataHashRef,
				  rangeFileKey => $rangeFileKey,
				  headerRef => $headerRef,
				  positionRef => \$extractColumnsCounter,
				  rangeFilePathRef => \$infilePath,
				 });
	    }
	    next;
	}
	if ( $_ =~/^(\S+)/ ) {	
	    
	    my @lineElements = split("\t",$_); #Loads range file line elements

	    if (defined($selectFeatureMatchingColumn) ) {

		$lineElements[ $selectFeatureMatchingColumn ] =~ s/\s/_/g; # Replace whitespace with "_"
		$selectData{$lineElements[ $selectFeatureMatchingColumn ]} = $lineElements[ $selectFeatureMatchingColumn ];
	    }

	    ## Create Interval Tree
	    if (@{$featureColumnsArrayRef}) {#Annotate vcf with features from range file

		&FeatureAnnotations({featureColumnsArrayRef => $featureColumnsArrayRef,
				     lineElementsArrayRef => \@lineElements,
				     rangeFileKey => $rangeFileKey,
				     paddingRef => $paddingRef,
				    });
	    }
	}
    }
    close(RRF);
    say STDERR "Finished Reading ".$rangeFileKey." file: ".$infilePath;
}


sub ReadInfileVCF {

##ReadInfileVCF

##Function : Reads infile in vcf format and adds and parses annotations
##Returns  : ""
##Arguments: $metaDataHashRef, $snpEffCmdHashRef, $rangeDataHashRef, $selectDataHashRef, $consequenceSeverityHashRef, $siftTermHashRef, $polyPhenTermHashRef, $rangeFeatureAnnotationColumnsArrayRef, $selectFeatureAnnotationColumnsArrayRef, $selectOutFilePath, $vcfParserVersion, $selectFeatureFile, $parseVEP, $writeSoftwareTag
##         : $metaDataHashRef                        => Vcf meta data {REF}
##         : $snpEffCmdHashRef                       => SnpEff meta data {REF}
##         : $rangeDataHashRef                       => Range file data {REF}
##         : $selectDataHashRef                      => Select file data {REF}
##         : $consequenceSeverityHashRef             => Consequence severity for SO-terms {REF}
##         : $siftTermHashRef                        => Sift prediction terms {REF}
##         : $polyPhenTermHashRef                    => PolyPhen prediction terms {REF}
##         : $rangeFeatureAnnotationColumnsArrayRef  => Range feature columns {REF}
##         : $selectFeatureAnnotationColumnsArrayRef => Select feature columns {REF}
##         : $selectOutFilePath                      => The select file path
##         : $vcfParserVersion                       => vcfParser version
##         : $selectFeatureFile                      => The select feature file
##         : $parseVEP                               => Parse VEP output
##         : $writeSoftwareTag                       => Write software tag to vcf header switch

    my ($argHashRef) = @_;

    ## Default(s)
    my $selectFeatureFile;
    my $parseVEP;
    my $writeSoftwareTag;

    ## Flatten argument(s)
    my $metaDataHashRef;
    my $snpEffCmdHashRef;
    my $rangeDataHashRef;
    my $selectDataHashRef;
    my $consequenceSeverityHashRef;
    my $siftTermHashRef;
    my $polyPhenTermHashRef;
    my $rangeFeatureAnnotationColumnsArrayRef;
    my $selectFeatureAnnotationColumnsArrayRef;
    my $selectOutFilePath;
    my $vcfParserVersion;

    my $tmpl = { 
	metaDataHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$metaDataHashRef},
	snpEffCmdHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$snpEffCmdHashRef},
	rangeDataHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$rangeDataHashRef},
	selectDataHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$selectDataHashRef},
	consequenceSeverityHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$consequenceSeverityHashRef},
	siftTermHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$siftTermHashRef},
	polyPhenTermHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$polyPhenTermHashRef},
	rangeFeatureAnnotationColumnsArrayRef => { required => 1, defined => 1, default => [], strict_type => 1, store => \$rangeFeatureAnnotationColumnsArrayRef},
	selectFeatureAnnotationColumnsArrayRef => { required => 1, defined => 1, default => [], strict_type => 1, store => \$selectFeatureAnnotationColumnsArrayRef},
	selectOutFilePath => { required => 1, defined => 1, strict_type => 1, store => \$selectOutFilePath},
	vcfParserVersion => { required => 1, defined => 1, strict_type => 1, store => \$vcfParserVersion},
	selectFeatureFile => { default => 0,
			       strict_type => 1, store => \$selectFeatureFile},
	parseVEP => { default => 0,
		      allow => [0, 1],
		      strict_type => 1, store => \$parseVEP},
	writeSoftwareTag => { default => 1,
			      allow => [0, 1],
			      strict_type => 1, store => \$writeSoftwareTag},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my @vepFormatField;
    my %vepFormatFieldColumn;

    my @printTSVFields = ("FeatureType");
    my @featureFields;

    my %vcfHeader;

    if ($selectFeatureFile) {

	open(WOSFTSV, ">".$selectOutFilePath) or die "Cannot open ".$selectOutFilePath.":".$!, "\n";
    }
    
    while (<>) {
	
	chomp $_;  # Remove newline
	
	if (m/^\s+$/) {	# Avoid blank lines
	    next;
	}
	if ($_=~/^##(\S+)=/) {  # MetaData

	    &ParseMetaData({metaDataHashRef => $metaDataHashRef,
			    metaDataString => $_,
			   });

	    if ($_=~/INFO\=\<ID\=(\w+)/) { # Collect all INFO keys
	    
		$vcfHeader{INFO}{$1} = $1; #Save to hash
	    }
	    if ($_=~/SnpSiftCmd\=/) { #Find SnpEff command meta line

		for my $database (keys %{${$snpEffCmdHashRef}{SnpEff}}) {

		    if ($_=~/${$snpEffCmdHashRef}{SnpEff}{$database}{File}/) { #SnpEff/Sift has been used to annotate input vcf
			
			unless (defined($vcfHeader{INFO}{$database})) { #Unless INFO header is already present add to metaDataHeader
			    
			    ${$snpEffCmdHashRef}{Present}{database}{$database} = $database; #Save which frequency db has been used for later
			    push(@{${$metaDataHashRef}{INFO}{$database}}, ${$snpEffCmdHashRef}{SnpEff}{$database}{INFO});

			    if (defined(${$snpEffCmdHashRef}{SnpEff}{$database}{FIX_INFO})) { #If FIX_INFO flag is present add to metaDataHeader

				push(@{${$metaDataHashRef}{FIX_INFO}{$database}}, ${$snpEffCmdHashRef}{SnpEff}{$database}{FIX_INFO});
			    }
			}
		    }
		}
		next;
	    }
	    if ($_=~/INFO\=\<ID\=CSQ/) { #Find VEP INFO Field

		if ($_=~/Format:\s(\S+)"\>/) { #Locate Format within VEP INFO meta line
		
		    @vepFormatField = split(/\|/, $1);  
		    
		    for (my $fieldCounter=0;$fieldCounter<scalar(@vepFormatField);$fieldCounter++) {
			
			$vepFormatFieldColumn{$vepFormatField[$fieldCounter]} = $fieldCounter; #Save the order of VEP features
		    }
		}
		if ($parseVEP) {
			
		    if ( ($vepFormatFieldColumn{SYMBOL}) && ($vepFormatFieldColumn{HGVSc}) && ($vepFormatFieldColumn{HGVSp})) {
			
			push(@{${$metaDataHashRef}{INFO}{HGVScp}}, '##INFO=<ID=HGVScp,Number=.,Type=String,Description="Transcript and protein functional annotation.">');
			push(@{${$metaDataHashRef}{INFO}{MostSevereConsequence}}, '##INFO=<ID=MostSevereConsequence,Number=.,Type=String,Description="Most severe genomic consequence.">');
			push(@{${$metaDataHashRef}{INFO}{GeneticRegionAnnotation}}, '##INFO=<ID=GeneticRegionAnnotation,Number=.,Type=String,Description="Genetic region that variant falls into.">');
			
		    }
		    if ($vepFormatFieldColumn{SIFT}) {
			
			push(@{${$metaDataHashRef}{INFO}{Sift}}, '##INFO=<ID=Sift,Number=.,Type=String,Description="Sift protein function prediction term">');
		    }
		    if ($vepFormatFieldColumn{PolyPhen}) {
			
			push(@{${$metaDataHashRef}{INFO}{PolyPhen}}, '##INFO=<ID=PolyPhen,Number=.,Type=String,Description="PolyPhen protein function prediction term">');
		    }
		}
		next;
	    }
	    next;
	}
	if ($_=~/^#CHROM/) {

	    &AddFeatureFileMetaDataToVcf({metaDataHashRef => $metaDataHashRef,
					  vcfHeaderHashRef => \%vcfHeader,
					  featureAnnotationColumnsArrayRef => $rangeFeatureAnnotationColumnsArrayRef,
					  dataHashRef => $rangeDataHashRef,
					  fileKey => "Range",
					 });
	    &AddFeatureFileMetaDataToVcf({metaDataHashRef => $metaDataHashRef,
					  vcfHeaderHashRef => \%vcfHeader,
					  featureAnnotationColumnsArrayRef => $selectFeatureAnnotationColumnsArrayRef,
					  dataHashRef => $selectDataHashRef,
					  fileKey => "Select",
					 });

	    if ($writeSoftwareTag) {

		&AddProgramToMeta({metaDataHashRef => $metaDataHashRef,
				   vcfParserVersion => $vcfParserVersion,
				  });
	    }
	    if (@{$selectFeatureAnnotationColumnsArrayRef}) { #SelectFile annotations

		&WriteMetaData({metaDataHashRef => $metaDataHashRef,
				FILEHANDLE => *STDOUT,
				SELECTFILEHANDLE => *WOSFTSV,
			       });
		say STDOUT $_;  #Write #CHROM header line
		say WOSFTSV $_;  #Write #CHROM header line
	    }
	    else {

		&WriteMetaData({metaDataHashRef => $metaDataHashRef,
				FILEHANDLE => *STDOUT,
			       });
		say STDOUT $_;  #Write #CHROM header line
	    }
	    
	    if ($parseVEP) {

		@featureFields = ("MostSevereConsequence", "GeneticRegionAnnotation");
		
		if ($vepFormatFieldColumn{SIFT}) {
		    
		    push(@featureFields, "Sift");
		}
		if ($vepFormatFieldColumn{PolyPhen}) {
		    
		    push(@featureFields, "PolyPhen");
		}
	    }
	    next;
	}
	if ( $_ =~/^(\S+)/ ) {	
	    
	    my %variantData;
	    my %selectedVariantData;
	    my %consequence;
	    my %noIDRegion;
	    my $variantLine;
	    my $selectedVariantLine;
	    my $sampleIDInfo;
	    my $foundSelectedFeatureCounter = 0;
	    my $transcriptsCounter = 0; #Tracks the number of transcripts to enable print of ", and ;" at correct position
	    my $selectedTranscriptCounter = 0; #Tracks the number of selected transcripts to enable print of ", and ;" at correct position

	    my @lineElements = split("\t",$_);  #Loads database elements description

	    for (my $lineElementsCounter=0;$lineElementsCounter<scalar(@lineElements);$lineElementsCounter++) { #Add until INFO field
			
		if ($lineElementsCounter < 7) { #Save fields until INFO field

		    $variantLine .= $lineElements[$lineElementsCounter]."\t";
		    $selectedVariantLine .= $lineElements[$lineElementsCounter]."\t"; #Copy vcf info to enable selective print downstream; 
		}
		elsif ($lineElementsCounter > 7) { #Save GT:PL: and sample(s) GT Call fields and add to proper line last
		    
		    if ($lineElementsCounter == (scalar(@lineElements) - 1)) {

			$sampleIDInfo .= $lineElements[$lineElementsCounter];
		    }
		    else {

			$sampleIDInfo .= $lineElements[$lineElementsCounter]."\t";
		    }
		}
	    }
	    for my $database (keys % {${$snpEffCmdHashRef}{Present}{database}}) { #Note that the vcf should only contain 1 database entry

		if ( ($database eq "Dbsnp129LCAF") || ($database eq "DbsnpLCAF") ) {

		    my @tempArray = split(/;/, $lineElements[7]);  #Split INFO field to key=value items

		    my $tempMafList = &FindAF({arrayRef => \@tempArray,
					       regexp => "Dbsnp\\S+CAF=",
					      });

		    if (defined($tempMafList)) {

			my $tempMaf;

			if ($tempMafList =~/^\[(.+)\],/) {  #Multiple entries in database for variant

			    $tempMafList = $';  #Pick last block as they should be identical. '			    
			}
			if ($tempMafList =~/\[(.+)\]/) {  #Single entry in database for variant

			    $tempMaf = $1;
			}
			
			my @tempMafs = sort {$a <=> $b} grep { $_ ne "." } split(",", $tempMaf); #Split on ",", remove entries containing only "." and sort remaining entries numerically
			if (scalar(@tempMafs) > 0) {
			    
			    ## Save LCAF frequency info
			    $variantLine .= $database."=".$tempMafs[0].";";
			    $selectedVariantLine .= $database."=".$tempMafs[0].";";
			}
		    }
		}
		elsif($database eq "ESPMAF") {
		    
		    my @tempArray = split(/;/, $lineElements[7]);  #Split INFO field to key=value items

		    my $tempMaf = &FindLCAF({arrayRef => \@tempArray,
					     regexp => "\\S+_MAF=",
					     frequencyPosition => "2",
					    });
	
		    if (defined($tempMaf)) {
			
			$tempMaf = $tempMaf / 100; #fraction for consistent representation
			
			## Save MAF frequency info  
			$variantLine .= $database."=".$tempMaf.";";
			$selectedVariantLine .= $database."=".$tempMaf.";";
		    }   
		}
		elsif($database eq "EXACMAXAF") {
		    
		    my @tempArray = split(/;/, $lineElements[7]);  #Split INFO field to key=value items
		    
		    my $tempMaf = &FindAF({arrayRef => \@tempArray,
					   regexp => "\\S+_EXACAFMAX_AF=",
					  });
		    
		    if (defined($tempMaf)) {
			
			## Save Alternative Allele frequency info  
			$variantLine .= $database."=".$tempMaf.";";
			$selectedVariantLine .= $database."=".$tempMaf.";";
		    }		    
		}
		elsif($database eq "phastCons100way_vertebrate_prediction_term") {
		    
		    my @tempArray = split(/;/, $lineElements[7]);  #Split INFO field to key=value items
		    
		    my $conservationTerm = &FindConserved({arrayRef => \@tempArray,
							   regexp => "\\S+_phastCons100way_vertebrate=",
							   scoreCutOff => 0.8,
							  });

		    if (defined($conservationTerm)) {
			
			## Save database info  
			$variantLine .= $database."=".$conservationTerm.";";
			$selectedVariantLine .= $database."=".$conservationTerm.";";
		    }
		}
		elsif($database eq "phyloP100way_vertebrate_prediction_term") {
		    
		    my @tempArray = split(/;/, $lineElements[7]);  #Split INFO field to key=value items
		    
		    my $conservationTerm = &FindConserved({arrayRef => \@tempArray,
							   regexp => "\\S+_phyloP100way_vertebrate=",
							   scoreCutOff => 2.5,
							  });

		    if (defined($conservationTerm)) {
			
			## Save database info
			$variantLine .= $database."=".$conservationTerm.";";
			$selectedVariantLine .= $database."=".$conservationTerm.";";
		    }
		}
		elsif($database eq "GERP++_RS_prediction_term") {
		    
		    
		    my @tempArray = split(/;/, $lineElements[7]);  #Split INFO field to key=value items
		    
		    my $conservationTerm = &FindConserved({arrayRef => \@tempArray,
							   regexp => "\\S+_GERP\\+\\+_RS=",
							   scoreCutOff => 2,
							  });

		    if (defined($conservationTerm)) {
			
			## Save database info
			$variantLine .= $database."=".$conservationTerm.";";
			$selectedVariantLine .= $database."=".$conservationTerm.";";
		    }
		}
	    }

	    ## Checks if an interval tree exists (per chr) and collects features from input array and adds annotations to line
	    %noIDRegion = &TreeAnnotations({hashRef => $selectDataHashRef,
					    lineElementsArrayRef => \@lineElements,
					    rangeFileKey => "SelectFile",
					    printLineRef => \$selectedVariantLine,
					   });  #noIDRegion is only for selectfile since all variants are passed to research file

	    ## Checks if an interval tree exists (per chr) and collects features from input array and adds annotations to line
	    &TreeAnnotations({rangeFileKey => "RangeFile",
			      lineElementsArrayRef => \@lineElements,
			      hashRef => $rangeDataHashRef,
			      printLineRef => \$variantLine,
			     });
	    
	    my @variantEffects = split(/;/, $lineElements[7]); #Split INFO field

	    ##Check that we have an INFO field
	    unless ($lineElements[7]) {

		say STDERR "No INFO field at line number: ".$.;
		say STDERR "Displaying malformed line: ".$_;
		exit 1;
	    }
	    my $CSQTranscripts;
	    
	    for (my $variantEffectCounter=0;$variantEffectCounter<scalar(@variantEffects);$variantEffectCounter++) {
		
		if ($parseVEP) {
		    
		    ## Find CSQ field and extract transripts that belong to select genes
		    if ($variantEffects[$variantEffectCounter]=~/CSQ\=(\S+)/) { #Find CSQ
			
			$CSQTranscripts = $1;
			my @transcripts = split(/,/, $1); #Split into transcripts elements
						    
			for (my $fieldCounter=0;$fieldCounter<scalar(@transcripts);$fieldCounter++) { #CSQ field
				
			    my @transcriptsEffects = split(/\|/, $transcripts[$fieldCounter]); #Split in "|"
	
			    if ( (defined($transcriptsEffects[ $vepFormatFieldColumn{Feature} ])) && ($transcriptsEffects[ $vepFormatFieldColumn{Feature} ] ne "") ) {

				if (defined($transcriptsEffects[ $vepFormatFieldColumn{SYMBOL} ])) { 
					
				    $variantData{Symbol} = $transcriptsEffects[ $vepFormatFieldColumn{SYMBOL} ];  #Save HGNC Symbol
					
				    if ($selectData{ $variantData{Symbol} }) { #Exists in selected Features

					if ($selectedTranscriptCounter > 0) {
						
					    $selectedVariantLine .= ",".$transcripts[$fieldCounter];
					}
					else {
						
					    $selectedVariantLine .= "CSQ=".$transcripts[$fieldCounter];
					}
					$selectedTranscriptCounter++;
				    }
				}
			    }

			    ## Always include all transcripts in research list
			    if ($transcriptsCounter > 0) {
				
				$variantLine .= ",".$transcripts[$fieldCounter];
			    }
			    else {
				
				$variantLine .= "CSQ=".$transcripts[$fieldCounter];
			    }
			    $transcriptsCounter++;
			}
		    }
		    else { #Not CSQ field from VEP 		    		    
			
			if ($variantEffectCounter == scalar(@variantEffects)-1) {
			    
			    $variantLine .= $variantEffects[$variantEffectCounter];
			    $selectedVariantLine .= $variantEffects[$variantEffectCounter];
			}
			else {

			    $variantLine .= $variantEffects[$variantEffectCounter].";";
			    $selectedVariantLine .= $variantEffects[$variantEffectCounter].";";
			}
		    }
		}
		else { #No CSQ field from VEP
		    
		    if ($variantEffectCounter == scalar(@variantEffects)-1) {
			
			$variantLine .= $variantEffects[$variantEffectCounter];
			$selectedVariantLine .= $variantEffects[$variantEffectCounter];
		    }
		    else {

			$variantLine .= $variantEffects[$variantEffectCounter].";";
			$selectedVariantLine .= $variantEffects[$variantEffectCounter].";";
		    }
		}
	    }
	    if ($parseVEP) {
		
		$transcriptsCounter = 0;
		$selectedTranscriptCounter = 0;

		if (defined($CSQTranscripts)) {

		    my @transcripts = split(/,/, $CSQTranscripts);
		    
		    for (my $fieldCounter=0;$fieldCounter<scalar(@transcripts);$fieldCounter++) { #CSQ field
			
			my @transcriptsEffects = split(/\|/, $transcripts[$fieldCounter]); #Split in "|"
			
			if ( (defined($transcriptsEffects[ $vepFormatFieldColumn{Feature} ]) ) && ($transcriptsEffects[ $vepFormatFieldColumn{Feature} ] ne "") && ($transcriptsEffects[ $vepFormatFieldColumn{Feature} ] !~/ENSR\d+/) ) { #All but regulatory regions, since there currently is no HGNC_Symbol annotated for them
			    
			    my $selectedTranscriptTracker = 0; #Track if any transcripts belong to selected features
			    
			    if (defined($transcriptsEffects[ $vepFormatFieldColumn{SYMBOL} ])) { #Save HGNC Symbol
				
				$variantData{Symbol} = $transcriptsEffects[ $vepFormatFieldColumn{SYMBOL} ];
				
				if (${$selectDataHashRef}{ $variantData{Symbol} }) { #Exists in selected Features
				    
				    $selectedTranscriptTracker = 1; #Record belongs to selected Features
				    my $alleleGeneEntry = $transcriptsEffects[ $vepFormatFieldColumn{Allele} ].":".$variantData{Symbol};
				   
				    &AddFieldToElementCounter({transcriptCounterRef => \$selectedTranscriptCounter,
							       lineRef => \$selectedVariantLine,
							       valueRef => \$alleleGeneEntry,
							       fieldID => "HGVScp=",
							      });
				}
		
				## Always include all transcripts in research list
				my $alleleGeneEntry = $transcriptsEffects[ $vepFormatFieldColumn{Allele} ].":".$variantData{Symbol};
				    
				&AddFieldToElementCounter({transcriptCounterRef => \$transcriptsCounter,
							   lineRef => \$variantLine,
							   valueRef => \$alleleGeneEntry,
							   fieldID => "HGVScp=",
							  });
			    }
			    &AddFieldToElement({selectedTranscriptTrackerRef => \$selectedTranscriptTracker,
						selectedLineRef => \$selectedVariantLine,
						lineRef => \$variantLine,
						valueRef => \$transcriptsEffects[ $vepFormatFieldColumn{Feature} ],
					       }); #Save transcript

			    if (defined($transcriptsEffects[ $vepFormatFieldColumn{Feature_type} ])) { #Save Feature_type
				
				if ($selectedTranscriptTracker == 1) {
				    
				    if ($selectedTranscriptCounter == 0) { #First selected gene
					
					$selectedVariantData{FeatureType} = $transcriptsEffects[ $vepFormatFieldColumn{Feature_type} ];
				    }
				    else {
					
					$selectedVariantData{FeatureType} .= ",".$transcriptsEffects[ $vepFormatFieldColumn{Feature_type} ];
				    }
				}
				## Always include all transcripts in research list				    
				if ($transcriptsCounter == 0) { #First Gene
				    
				    $variantData{FeatureType} = $transcriptsEffects[ $vepFormatFieldColumn{Feature_type} ];
				}
				else {
				    
				    $variantData{FeatureType} .= ",".$transcriptsEffects[ $vepFormatFieldColumn{Feature_type} ];
				}
			    }
			    if (defined($transcriptsEffects[ $vepFormatFieldColumn{Consequence} ])) { #Save Consequence
				
				&AddFieldToElement({selectedTranscriptTrackerRef => \$selectedTranscriptTracker,
						    selectedLineRef => \$selectedVariantLine,
						    lineRef => \$variantLine,
						    valueRef => \$transcriptsEffects[ $vepFormatFieldColumn{Consequence} ],
						   });

				my @consequences = split(/\&/, $transcriptsEffects[ $vepFormatFieldColumn{Consequence} ]); #Find "MostSevereConsequence
				my $allele = $transcriptsEffects[ $vepFormatFieldColumn{Allele} ];

				foreach my $consequenceTerm (@consequences) {
				    
				    &CheckTerms({termHashRef => $consequenceSeverityHashRef,
						 termRef => \$consequenceTerm,
						 termName => "SO",
						});
				    
				    if (defined($variantData{Symbol})) {

					if(defined($consequence{ $variantData{Symbol} }{$allele}{Score})) { #Compare to previous record
					
					    if (${$consequenceSeverityHashRef}{$consequenceTerm}{Rank} < $consequence{ $variantData{Symbol} }{$allele}{Score}) { #Collect most severe consequence
						
						&AddToConsequenceHash({hashKeyRef => \$consequence{ $variantData{Symbol} }{$allele}{Score},
								       valueRef => \${$consequenceSeverityHashRef}{$consequenceTerm}{Rank},
								      });
						
						&AddToConsequenceHash({hashKeyRef => \$consequence{ $variantData{Symbol} }{$allele}{GeneticRegionAnnotation},
								       valueRef => \${$consequenceSeverityHashRef}{$consequenceTerm}{GeneticRegionAnnotation},
								      });
						&AddToConsequenceHash({hashKeyRef => \$consequence{ $variantData{Symbol} }{$allele}{MostSevereConsequence},
								       valueRef => \$consequenceTerm,
								      });
						
						if (defined($vepFormatFieldColumn{SIFT}) ) {

						    &CheckTerms({termHashRef => $siftTermHashRef,
								 termRef => \$transcriptsEffects[ $vepFormatFieldColumn{SIFT} ],
								 termName => "Sift",
								});
						    &AddToConsequenceHash({hashKeyRef => \$consequence{ $variantData{Symbol} }{$allele}{Sift},
									   valueRef => \$transcriptsEffects[ $vepFormatFieldColumn{SIFT} ],
									  });
						}
						if (defined($vepFormatFieldColumn{PolyPhen}) ) {
						    
						    &CheckTerms({termHashRef => $polyPhenTermHashRef,
								 termRef => \$transcriptsEffects[ $vepFormatFieldColumn{PolyPhen} ],
								 termName => "polyPhen",
								});
						    &AddToConsequenceHash({hashKeyRef => \$consequence{ $variantData{Symbol} }{$allele}{PolyPhen},
									   valueRef => \$transcriptsEffects[ $vepFormatFieldColumn{PolyPhen} ],
									  });
						}
					    }
					}
					else { #First pass
					
					    &AddToConsequenceHash({hashKeyRef => \$consequence{ $variantData{Symbol} }{$allele}{Score},
								   valueRef => \${$consequenceSeverityHashRef}{$consequenceTerm}{Rank},
								  });
					    &AddToConsequenceHash({hashKeyRef => \$consequence{ $variantData{Symbol} }{$allele}{GeneticRegionAnnotation},
								   valueRef => \${$consequenceSeverityHashRef}{$consequenceTerm}{GeneticRegionAnnotation},
								  });
					    &AddToConsequenceHash({hashKeyRef => \$consequence{ $variantData{Symbol} }{$allele}{MostSevereConsequence},
								   valueRef => \$consequenceTerm,
								  });    
					    
					    if (defined($vepFormatFieldColumn{SIFT}) ) {
						
						&AddToConsequenceHash({hashKeyRef => \$consequence{ $variantData{Symbol} }{$allele}{Sift},
								       valueRef => \$transcriptsEffects[ $vepFormatFieldColumn{SIFT} ],
								      });
					    }
					    if (defined($vepFormatFieldColumn{PolyPhen}) ) {
						
						&AddToConsequenceHash({hashKeyRef => \$consequence{ $variantData{Symbol} }{$allele}{PolyPhen},
								       valueRef => \$transcriptsEffects[ $vepFormatFieldColumn{PolyPhen} ],
								      });
					    }
					}
				    }
				}
			    }
			    if ( (defined($transcriptsEffects[ $vepFormatFieldColumn{STRAND} ])) && ($transcriptsEffects[ $vepFormatFieldColumn{STRAND} ] ne "") ) { #Save strand 

				if ($transcriptsEffects[ $vepFormatFieldColumn{STRAND} ] == 1) { #Remap to "+" or "-"

				    $variantData{Strand} = "s.+";
				}
				else {
				    
				    $variantData{Strand} = "s.-";
				}
				&AddFieldToElement({selectedTranscriptTrackerRef => \$selectedTranscriptTracker,
						    selectedLineRef => \$selectedVariantLine,
						    lineRef => \$variantLine,
						    valueRef => \$variantData{Strand},
						   });
			    }
			    if (defined($transcriptsEffects[ $vepFormatFieldColumn{EXON} ]) && $transcriptsEffects[ $vepFormatFieldColumn{EXON}] =~/^\d/) {
				
				$variantData{Exon} = "e.".$transcriptsEffects[ $vepFormatFieldColumn{EXON} ]; #Save exon number "X/Y"

				&AddFieldToElement({selectedTranscriptTrackerRef => \$selectedTranscriptTracker,
						    selectedLineRef => \$selectedVariantLine,
						    lineRef => \$variantLine,
						    valueRef => \$variantData{Exon},
						   });
			    }
			    if (defined($transcriptsEffects[ $vepFormatFieldColumn{INTRON} ]) && $transcriptsEffects[ $vepFormatFieldColumn{INTRON} ] =~/^\d/) {
				
				$variantData{Intron} = "i.".$transcriptsEffects[ $vepFormatFieldColumn{INTRON} ]; #Save intron number "X/Y"

				&AddFieldToElement({selectedTranscriptTrackerRef => \$selectedTranscriptTracker,
						    selectedLineRef => \$selectedVariantLine,
						    lineRef => \$variantLine,
						    valueRef => \$variantData{Intron},
						   });
			    }
			    if ( (defined($transcriptsEffects[ $vepFormatFieldColumn{HGVSc} ])) && ($transcriptsEffects[ $vepFormatFieldColumn{HGVSc} ] ne "")) { #Save HGVS cDNA change
				
				my @cDNAChanges = split(/:/, $transcriptsEffects[ $vepFormatFieldColumn{HGVSc} ]);

				&AddFieldToElement({selectedTranscriptTrackerRef => \$selectedTranscriptTracker,
						    selectedLineRef => \$selectedVariantLine,
						    lineRef => \$variantLine,
						    valueRef => \$cDNAChanges[1],
						   });
			    }
			    if ( (defined($transcriptsEffects[ $vepFormatFieldColumn{HGVSp} ])) && ($transcriptsEffects[ $vepFormatFieldColumn{HGVSp} ] ne "") ) {
				
				my @pAAChanges = split(/:/, $transcriptsEffects[ $vepFormatFieldColumn{HGVSp} ]);

				&AddFieldToElement({selectedTranscriptTrackerRef => \$selectedTranscriptTracker,
						    selectedLineRef => \$selectedVariantLine,
						    lineRef => \$variantLine,
						    valueRef => \$pAAChanges[1],
						   });
			    }
			    if ($selectedTranscriptTracker == 1) {
				
				$selectedTranscriptCounter++;
			    }
			    ## Always include all transcripts in research list	
			    $transcriptsCounter++;
			}
		    }
		    &CollectConsequenceGenes({consequenceHashRef => \%consequence,
					      selectDataHashRef => $selectDataHashRef,
					      fieldsArrayRef => \@featureFields,
					      selectedVariantLineRef => \$selectedVariantLine,
					      variantLineRef => \$variantLine,
					     });
		}
	    }

	    $selectedVariantLine .= "\t".$sampleIDInfo;
	    $variantLine .= "\t".$sampleIDInfo;	

	    if ($parseVEP) {
		
		if ($selectedTranscriptCounter > 0) { #Write to selected file
		    
		    say WOSFTSV $selectedVariantLine;
		}
		if ($transcriptsCounter > 0) { #Write to transcript file
		    
		    if (%noIDRegion) {  #No HGNC_symbol or EnsemblGeneID, but clinically relevant
		
			say WOSFTSV $selectedVariantLine;
		    }
		    say STDOUT $variantLine;
		}
		elsif ( (! $selectedTranscriptCounter) && (! $transcriptsCounter) ) {

		    if (%noIDRegion) {  #No HGNC_symbol or EnsemblGeneID, but clinically relevant

			say WOSFTSV $selectedVariantLine;
		    }
		    say STDOUT $variantLine;
		}
	    }
	    else {

		say STDOUT $variantLine;	
	    }
	}
    }
    if ($selectFeatureFile) {

	close(WOSFTSV);
    }
    say STDERR "Finished Processing VCF";
}


sub PrintTSVField {
##Prints a TSV field if defined

    my $arrayRef = $_[0]; #FILEHANDLES
    my $valueRef = $_[1];

    for (my $filehandleCounter=0;$filehandleCounter<scalar(@{$arrayRef});$filehandleCounter++) {

	if (defined($$valueRef)) {
	    
	    print {${$arrayRef}[$filehandleCounter]} "\t".$$valueRef;
	}
	else {
	    
	    print {${$arrayRef}[$filehandleCounter]} "\t";
	}
    }
}


sub AddToConsequenceHash {

##AddToConsequenceHash
    
##Function : Adds the most severe consequence or prediction to gene.
##Returns  : ""
##Arguments: $hashKeyRef, $valueRef
##         : $hashKeyRef => The hash key to update {REF}
##         : $valueRef   => The value to add to hash key

    my ($argHashRef) = @_;
	
    ## Flatten argument(s)
    my $hashKeyRef;
    my $valueRef;

    my $tmpl = { 
	hashKeyRef => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$hashKeyRef},
	valueRef => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$valueRef},    
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];
    
    $$hashKeyRef = $$valueRef;
}

sub AddFieldToElementCounter {

##AddFieldToElementCounter
    
##Function : Adds a field to an element. Adjust addition depending on if field has been seen before.
##Returns  : ""
##Arguments: $transcriptCounterRef, $lineRef, $valueRef, $fieldID, $separator
##         : $transcriptCounterRef => The transcript counter {REF}
##         : $lineRef              => Variant line to add annotations to {REF}
##         : $valueRef             => Field value {REF}
##         : $fieldID              => Field key ID
##         : $separator            => Separator for field

    my ($argHashRef) = @_;

    ## Default(s)
    my $separator;

    ## Flatten argument(s)
    my $transcriptCounterRef;
    my $lineRef;
    my $valueRef;
    my $fieldID;

    my $tmpl = { 
	transcriptCounterRef => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$transcriptCounterRef},
	lineRef => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$lineRef},
	separator => { default => ",", strict_type => 1, store => \$separator},
	valueRef => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$valueRef},
	fieldID => { required => 1, defined => 1, strict_type => 1, store => \$fieldID},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];
    
    if (! $$transcriptCounterRef) {
	
	$$lineRef .= ";".$fieldID.$$valueRef; #First selected feature
    }
    else {

	$$lineRef .= $separator.$$valueRef; #Following features
    }
}

sub AddFieldToElement {

##AddFieldToElement
    
##Function : Adds adds a field to an element.
##Returns  : ""
##Arguments: $selectedTranscriptTrackerRef, $selectedLineRef, $lineRef, $valueRef, $separator
##         : $selectedTranscriptTrackerRef => The selected transcript tracker for belonging to select file {REF}
##         : $selectedLineRef              => Selected line to add annotations to {REF}
##         : $lineRef                      => Variant line to add annotations to {REF}
##         : $valueRef                     => Field value {REF}
##         : $separator                    => Separator for field

    my ($argHashRef) = @_;

    ## Default(s)
    my $separator;

    ## Flatten argument(s)
    my $selectedTranscriptTrackerRef;
    my $selectedLineRef;
    my $lineRef;
    my $valueRef;

    my $tmpl = { 
	selectedTranscriptTrackerRef => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$selectedTranscriptTrackerRef},
	selectedLineRef => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$selectedLineRef},
	lineRef => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$lineRef},
	separator => { default => ":", strict_type => 1, store => \$separator},
	valueRef => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$valueRef},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];
    
    if ($$selectedTranscriptTrackerRef) {
	
	$$selectedLineRef .= $separator.$$valueRef;
    }
    ## Always include all transcripts in orphan list
    $$lineRef .= $separator.$$valueRef;
}

sub CollectConsequenceGenes {
 
##CollectConsequenceGenes
    
##Function : Collects all consequence and predictors per gene and adds info to line to be written.
##Returns  : ""
##Arguments: $consequenceHashRef, $selectDataHashRef, $fieldsArrayRef, $selectedVariantLineRef, $variantVariantLineRef
##         : $consequenceHashRef     => Consequence(s) for each gene {REF}
##         : $selectDataHashRef      => Select file data {REF}
##         : $fieldsArrayRef         => Features to be processed as determined by CSQ {REF}
##         : $selectedVariantLineRef => Selected line to add annotations to {REF}
##         : $variantVariantLineRef  => Variant line to add annotations to {REF}

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $consequenceHashRef;
    my $selectDataHashRef;
    my $fieldsArrayRef;
    my $selectedVariantLineRef;
    my $variantLineRef;

    my $tmpl = { 
	consequenceHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$consequenceHashRef},
	selectDataHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$selectDataHashRef},
	fieldsArrayRef => { required => 1, defined => 1, default => [], strict_type => 1, store => \$fieldsArrayRef},
	selectedVariantLineRef => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$selectedVariantLineRef},
	variantLineRef => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$variantLineRef},
    };

    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];
    
    my %geneCounter;
    my %selectedGeneCounter;
    my @tempFields;
    my @selectedTempFields;
    
    for (my $fieldCounter=0;$fieldCounter<scalar(@{$fieldsArrayRef});$fieldCounter++) { #Set transcript counter to "0"
	
	$selectedGeneCounter{$fieldCounter} = 0;
	$geneCounter{$fieldCounter} = 0;
    }
    for my $genes (keys %{$consequenceHashRef}) {
	
	for (my $fieldCounter=0;$fieldCounter<scalar(@{$fieldsArrayRef});$fieldCounter++) {
	    
	    if (${$selectDataHashRef}{$genes}) { #Exists in selected Features
		
		&CollectConsequenceField({geneCounterHashRef => \%selectedGeneCounter,
					  consequenceHashRef => \%{$consequenceHashRef},
					  fieldsArrayRef => \@{$fieldsArrayRef},
					  selectedArrayRef => \@selectedTempFields,
					  fieldCounterRef => \$fieldCounter,
					  geneRef => \$genes,
					 });
	    }

	    ## Always include all transcripts in research list
	    &CollectConsequenceField({geneCounterHashRef => \%geneCounter,
				      consequenceHashRef => \%{$consequenceHashRef},
				      fieldsArrayRef => \@{$fieldsArrayRef},
				      selectedArrayRef => \@tempFields,
				      fieldCounterRef => \$fieldCounter,
				      geneRef => \$genes,
				     });
	}
    }
    
    ## Adds to present line
    &AddToLine({fieldsArrayRef => \@{$fieldsArrayRef},
		tempArrayRef => \@selectedTempFields,
		lineRef => \$$selectedVariantLineRef,
	       });

    ## Adds to present line
    &AddToLine({fieldsArrayRef => \@{$fieldsArrayRef},
		tempArrayRef => \@tempFields,
		lineRef => \$$variantLineRef,
	       });
}

sub CollectConsequenceField {
   
##CollectConsequenceField
    
##Function : Collects consequences for features in @featureFields to temporary array for adding to line once all information are collected.
##Returns  : ""
##Arguments: $geneCounterHashRef, $consequenceHashRef, $fieldsArrayRef, $selectedArrayRef, $fieldCounterRef, $geneRef
##         : $geneCounterHashRef => Counts the number of transcripts per gene {REF}
##         : $consequenceHashRef => Consequence(s) for each gene {REF}
##         : $fieldsArrayRef     => Features to be processed as determined by CSQ {REF}
##         : $selectedArrayRef   => Selected array {REF}
##         : $fieldCounterRef    => Field number in feature {REF}
##         : $geneRef            => The gene symbol {REF}

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $geneCounterHashRef;
    my $consequenceHashRef;
    my $fieldsArrayRef;
    my $selectedArrayRef;
    my $fieldCounterRef;
    my $geneRef;

     my $tmpl = { 
	geneCounterHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$geneCounterHashRef},
	consequenceHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$consequenceHashRef},
	fieldsArrayRef => { required => 1, defined => 1, default => [], strict_type => 1, store => \$fieldsArrayRef},
	selectedArrayRef => { required => 1, defined => 1, default => [], strict_type => 1, store => \$selectedArrayRef},
	fieldCounterRef => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$fieldCounterRef},
	geneRef => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$geneRef},
    };

    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my $alleleNumber = scalar(keys %{${$consequenceHashRef}{$$geneRef}});  #Enable tracking of multiple alleles
    my $fieldRef = \${$fieldsArrayRef}[$$fieldCounterRef];  #Alias

    if (! $$geneCounterHashRef{$$fieldCounterRef}) { #First time
	
	my $allelCounter = 0;

	for my $allele (keys %{${$consequenceHashRef}{$$geneRef}}) {  #All alleles
	    
	    if (defined($$consequenceHashRef{$$geneRef}{$allele}{$$fieldRef}) && ($$consequenceHashRef{$$geneRef}{$allele}{$$fieldRef} ne "")) { #If feature exists - else do nothing
		
		if (! $allelCounter) {
		    
		    $$selectedArrayRef[$$fieldCounterRef] .= ";".$$fieldRef."=".$$geneRef.":".$allele."|".$$consequenceHashRef{$$geneRef}{$allele}{$$fieldRef};
		}
		else {
		    
		    $$selectedArrayRef[$$fieldCounterRef] .= ",".$$geneRef.":".$allele."|".$$consequenceHashRef{$$geneRef}{$allele}{$$fieldRef};
		}
		$$geneCounterHashRef{$$fieldCounterRef}++;
		$allelCounter++;
	    }
	}
    }
    else { #Subsequent passes
	
	for my $allele (keys %{${$consequenceHashRef}{$$geneRef}}) {  #All alleles
	    
	    if (defined($$consequenceHashRef{$$geneRef}{$allele}{$$fieldRef}) && ($$consequenceHashRef{$$geneRef}{$allele}{$$fieldRef} ne "")) { #If feature exists - else do nothing
		
		$$selectedArrayRef[$$fieldCounterRef] .= ",".$$geneRef.":".$allele."|".$$consequenceHashRef{$$geneRef}{$allele}{$$fieldRef};
	    }
	}
    }    
}

sub AddToLine {

##AddToLine
    
##Function : Adds to present line.
##Returns  : ""
##Arguments: $fieldsArrayRef, $tempArrayRef, $lineRef
##         : $fieldsArrayRef => Features to be processed as determined by CSQ {REF}
##         : $tempArrayRef   => Annotations for feature {REF}
##         : $lineRef        => Line to add to {REF}

    my ($argHashRef) = @_;
    
    ## Flatten argument(s)
    my $fieldsArrayRef;
    my $tempArrayRef;
    my $lineRef;

    my $tmpl = { 
	fieldsArrayRef => { required => 1, defined => 1, default => [], strict_type => 1, store => \$fieldsArrayRef},
	tempArrayRef => { required => 1, defined => 1, default => [], strict_type => 1, store => \$tempArrayRef},
	lineRef => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$lineRef},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];
    
    for (my $arrayFieldCounter=0;$arrayFieldCounter<scalar(@{$fieldsArrayRef});$arrayFieldCounter++) {
	
	if (defined($$tempArrayRef[$arrayFieldCounter])) {
	    
	    $$lineRef .= $$tempArrayRef[$arrayFieldCounter];
	}
    }
}

sub ConvertToRange {

##ConvertToRange
    
##Function : Converts vcf sv to corresponding range coordinates.
##Returns  : "$finalStartPosition, $finalStopPosition"
##Arguments: $fieldsArrayRef
##         : $fieldsArrayRef => Holds the chromosomal coordinates and allel data

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $fieldsArrayRef;

    my $tmpl = { 
	fieldsArrayRef => { required => 1, defined => 1, default => [], strict_type => 1, store => \$fieldsArrayRef},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];
   
    my $chromosome = $$fieldsArrayRef[0];
    my $startPosition = $$fieldsArrayRef[1];
    my $referenceAllele = $$fieldsArrayRef[3];
    my $alternativeAllele = $$fieldsArrayRef[4];

    my $finalStartPosition = $startPosition; #The most "uppstream" position per variant
    my $finalStopPosition = 0; #The most "downstream" position per variant
	
    ## Convert to upper case
    ($referenceAllele, $alternativeAllele) = (uc $referenceAllele, uc $alternativeAllele);
    
    if ($alternativeAllele eq ".") { #No Variant Call
	
	next;
    }
    my @alternativeAlleles = split(/,/, $$fieldsArrayRef[4]);
    
    for (my $allelCounter=0;$allelCounter<scalar(@alternativeAlleles);$allelCounter++) {
	
	my ($head, $newstart, $newend, $newref, $newalt);
	
	if (length ($referenceAllele) == 1 and length ($alternativeAlleles[$allelCounter]) == 1) { #SNV
	    
	    ($newstart, $newend) = ($startPosition, $startPosition+length($referenceAllele)-1);
	    ($newref, $newalt) = ($referenceAllele, $alternativeAlleles[$allelCounter]);
	}
	elsif (length($referenceAllele) >= length($alternativeAlleles[$allelCounter])) { #deletion or block substitution
	    
	    $head = substr ($referenceAllele, 0, length($alternativeAlleles[$allelCounter]));
	    
	    if ($head eq $alternativeAlleles[$allelCounter]) {
  
		($newstart, $newend) = ($startPosition+length($head), $startPosition + length($referenceAllele)-1);
		($newref, $newalt) = (substr($referenceAllele, length($alternativeAlleles[$allelCounter])), '-');
	    } 
	    else {

		($newstart, $newend) = ($startPosition, $startPosition+length($referenceAllele)-1);
		($newref, $newalt) = ($referenceAllele, $alternativeAlleles[$allelCounter]);
	    }
	}
	elsif (length($referenceAllele) < length($alternativeAlleles[$allelCounter])) { #insertion or block substitution

	    $head = substr($alternativeAlleles[$allelCounter], 0, length($referenceAllele));
	    
	    if ($head eq $referenceAllele) {
	
		($newstart, $newend) = ($startPosition+length($referenceAllele)-1, $startPosition+length($referenceAllele)-1);
		($newref, $newalt) = ('-', substr($alternativeAlleles[$allelCounter], length($referenceAllele)));
	    } 
	    else {

		($newstart, $newend) = ($startPosition, $startPosition+length($referenceAllele)-1);
		($newref, $newalt) = ($referenceAllele, $alternativeAlleles[$allelCounter]);
	    }
	}

	## Collect largest range per variant based on all alternativeAlleles
	if ($finalStartPosition < $newstart) { #New start is upstream of old
	    
	    $finalStartPosition = $newstart;
	}
	if ($finalStopPosition < $newend) { #New end is downstream of old
	    
	    $finalStopPosition = $newend;
	}
    }
    return $finalStartPosition, $finalStopPosition;
}


sub AddMetaDataINFO {

##AddMetaDataINFO
    
##Function : Adds arbitrary INFO fields to hash based on supplied headers unless header is already defined.
##Returns  : ""
##Arguments: $metaDataHashRef, $rangeFileKey, $headerRef, $positionRef, $rangeFilePathRef
##         : $metaDataHashRef  => Hash to store metaData in {REF}
##         : $rangeFileKey     => Range file key
##         : $headerRef        => Header from range file {REF}
##         : $positionRef      => Column position in supplied range file {REF}
##         : $rangeFilePathRef => Range file path {REF}

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $metaDataHashRef;
    my $rangeFileKey;
    my $headerRef;
    my $positionRef;
    my $rangeFilePathRef;

    my $tmpl = { 
	metaDataHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$metaDataHashRef},
	rangeFileKey => { required => 1, defined => 1, strict_type => 1, store => \$rangeFileKey},
	headerRef => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$headerRef},
	positionRef => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$positionRef},
	rangeFilePathRef => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$rangeFilePathRef},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    if (defined($$metaDataHashRef{$rangeFileKey}{$$headerRef})) { #Add INFO from predefined entries

	$$metaDataHashRef{Present}{$$headerRef}{INFO} = $$metaDataHashRef{$rangeFileKey}{$$headerRef}{INFO};
	$$metaDataHashRef{Present}{$$headerRef}{ColumnOrder} = $$positionRef; #Column position in supplied range input file
    }
    else { #Add arbitrary INFO field using input header

	$$metaDataHashRef{Present}{$$headerRef}{INFO} = q?##INFO=<ID=?.$$headerRef.q?,Number=.,Type=String,Description="String taken from ?.$$rangeFilePathRef.q?">?;
	$$metaDataHashRef{Present}{$$headerRef}{ColumnOrder} = $$positionRef; #Column position in supplied -sf_ac
    }
}

sub FeatureAnnotations {

##FeatureAnnotations
    
##Function : Creates the interval tree(s) for range and select feature files. Adds corresponding INFO fields to metadata.
##Returns  : ""
##Arguments: $featureColumnsArrayRef, $lineElementsArrayRef, $rangeFileKey, $paddingRef
##         : $featureColumnsArrayRef => Feature annotations columns {REF}
##         : $lineElementsArrayRef   => Feature file line {REF} 
##         : $rangeFileKey           => Range file key
##         : $paddingRef             => The nucleotide distance to pad the range with {REF}

    my ($argHashRef) = @_;

    ## Flatten argument(s
    my $featureColumnsArrayRef;
    my $lineElementsArrayRef;
    my $rangeFileKey;
    my $paddingRef;

    my $tmpl = { 
	featureColumnsArrayRef => { required => 1, defined => 1, default => [], strict_type => 1, store => \$featureColumnsArrayRef},
	lineElementsArrayRef => { required => 1, defined => 1, default => [], strict_type => 1, store => \$lineElementsArrayRef},
	rangeFileKey => { required => 1, defined => 1, strict_type => 1, store => \$rangeFileKey},
	paddingRef => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$paddingRef},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];
    
    my $features; #Features to collect (Format: ";" separated elements)
    
    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar(@{$featureColumnsArrayRef});$extractColumnsCounter++) { #Defines what scalar to store

	my $featureColumnRef = \${$featureColumnsArrayRef}[$extractColumnsCounter];  #Alias

	if(defined($$lineElementsArrayRef[ $$featureColumnRef ])) {

	    $$lineElementsArrayRef[ $$featureColumnRef ] =~ tr/ /_/; #Remove whitespace since this is not allowed in vcf INFO field
	
	    if (!$extractColumnsCounter) {
		
		$features .= $$lineElementsArrayRef[ $$featureColumnRef ];
	    }		    
	    else {
		
		$features .= ";".$$lineElementsArrayRef[ $$featureColumnRef ];
	    }
	}
    }
    unless(defined($tree{$rangeFileKey}{ $$lineElementsArrayRef[0] })) { #Only create once per firstKey
	
	$tree{$rangeFileKey}{ $$lineElementsArrayRef[0] } = Set::IntervalTree->new(); #Create tree
    }
    my $paddedStart = $$lineElementsArrayRef[1] - $$paddingRef;
    my $paddedStop = $$lineElementsArrayRef[2] + $$paddingRef;
    $tree{$rangeFileKey}{ $$lineElementsArrayRef[0] }->insert($features, $paddedStart, $paddedStop); #Store range and ";" sep string
}


sub TreeAnnotations {

##TreeAnnotations
    
##Function : Checks if an interval tree exists (per chr) and collects features from input array and adds annotations to line.
##Returns  : ""
##Arguments: $hashRef, $lineElementsArrayRef, $rangeFileKey, $printLineRef
##         : $hashRef              => Range file hash {REF}
##         : $lineElementsArrayRef => Infile vcf line elements array {REF}
##         : $rangeFileKey         => Range file key
##         : $printLineRef         => Line to add annotations to {REF}

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $hashRef;
    my $lineElementsArrayRef;
    my $rangeFileKey;
    my $printLineRef;

    my $tmpl = { 
	    hashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$hashRef},
	    lineElementsArrayRef => { required => 1, defined => 1, default => [], strict_type => 1, store => \$lineElementsArrayRef},
	    rangeFileKey => { required => 1, defined => 1, strict_type => 1, store => \$rangeFileKey},
	    printLineRef => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$printLineRef},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];
    
    my %noIDRegion;  #No HGNC_symbol or ensemblGeneID, but still clinically releveant e.g. mtD-loop

    if(defined($tree{$rangeFileKey}{ $$lineElementsArrayRef[0] }) ) { #Range annotations
	
	my $feature; #Features to be collected
	
	## #Convert SVs to range coordinates from vcf coordinates
	my ($start, $stop) = &ConvertToRange({fieldsArrayRef => $lineElementsArrayRef,
					     });
	
	if ($start eq $stop) { #SNV
	    
	    $feature = $tree{$rangeFileKey}{ $$lineElementsArrayRef[0] }->fetch($start, $stop+1); #Add 1 to SNV to create range input.
	}
	else {#Range input
	    
	    $feature = $tree{$rangeFileKey}{ $$lineElementsArrayRef[0] }->fetch($start, $stop);
	}
	if (@{$feature}) { #Features found in tree
	    
	    my %collectedAnnotations; #Collect all features before adding to line
	    
	    for (my $featureCounter=0;$featureCounter<scalar(@{$feature});$featureCounter++) { #All features

		my @annotations = split(/;/, @{$feature}[$featureCounter]); #Split feature array ref into annotations
		
		for (my $annotationsCounter=0;$annotationsCounter<scalar(@annotations);$annotationsCounter++) { #All annotations
		    
		    if( (defined($annotations[$annotationsCounter])) && ($annotations[$annotationsCounter] ne "") ) {
			
			push(@{$collectedAnnotations{$annotationsCounter}}, $annotations[$annotationsCounter]);
		    }
		    if ($featureCounter == (scalar(@{$feature}-1)) ) { #Last for this feature tuple

			for my $rangeAnnotation (keys % {$$hashRef{Present}}) { #All selected annotations
			
			    if ($$hashRef{Present}{$rangeAnnotation}{ColumnOrder} eq $annotationsCounter) { #Correct feature
			
				if ($rangeAnnotation eq "Clinical_db_gene_annotation") {  #Special case, which is global and not gene centric

				    ## Collect unique elements from array reference and return array reference with unique elements
				    my $uniqueRef = &UniqElements({arrayRef => \@{$collectedAnnotations{$annotationsCounter}},
								  });

				    @{$collectedAnnotations{$annotationsCounter}} = @{$uniqueRef};
				}
				if ($rangeAnnotation eq "No_hgnc_symbol") {  #Special case, where there is no HGNC or Ensembl gene ID but the region should be included in the select file anyway
				    my $idKey = join ("_", @{$lineElementsArrayRef}[0..1, 3..4]);
				    $noIDRegion{$idKey}++;
				}
				if ( (defined($collectedAnnotations{$annotationsCounter})) && (@{$collectedAnnotations{$annotationsCounter}}) ) {
				
				    $$printLineRef .= $rangeAnnotation."=".join(",", @{$collectedAnnotations{$annotationsCounter}}).";"; #Add to corresponding line
				}
			    }
			}
		    }
		}
	    }
	}
    }
    return %noIDRegion;
}


sub AddProgramToMeta {
    
##AddProgramToMeta
    
##Function : Adds the program version and run date to the vcf meta-information section
##Returns  : ""
##Arguments: $metaDataHashRef, $vcfParserVersion
##         : $metaDataHashRef  => Vcf meta data {REF}
##         : $vcfParserVersion => vcfParser version
 
    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $metaDataHashRef;
    my $vcfParserVersion;
    
    my $tmpl = { 
	metaDataHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$metaDataHashRef},
	vcfParserVersion => { required => 1, defined => 1, strict_type => 1, store => \$vcfParserVersion},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my ($base, $script) = (`date +%Y%m%d`,`basename $0`);  #Catches current date and script name
    chomp($base,$script);  #Remove \n;
    push(@{${$metaDataHashRef}{Software}{$script}}, "##Software=<ID=".$script.",Version=".$vcfParserVersion.",Date=".$base);
}


sub FindAF {

##FindAF
    
##Function : Adds the least alternative allele(s) frequency to each line
##Returns  : ""
##Arguments: $arrayRef, $regexp
##         : $arrayRef => The INFO array {REF}
##         : $regexp   => The regexp to used to locate correct ID field

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $arrayRef;
    my $regexp;

    my $tmpl = { 
	arrayRef => { required => 1, defined => 1, default => [], strict_type => 1, store => \$arrayRef},
	regexp => { required => 1, defined => 1, strict_type => 1, store => \$regexp},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my $tempMaf;
    
    for my $element (@{$arrayRef}) {
	
	if ($element =~/$regexp/) {  #Find the key=value field
	    
	    my @value = split(/=/, $element);  #Split key=value pair

	    $tempMaf = $value[1];  #Collect whole string to represent all possible alleles
	}
    }
    return $tempMaf
}


sub FindConserved {

##FindConserved
    
##Function : Adds the least common alternative allele frequency to each line
##Returns  : "" or string prediction terms
##Arguments: $arrayRef, $regexp, $scoreCutOff
##         : $arrayRef    => The INFO array {REF}
##         : $regexp      => The regexp to used to locate correct ID field
##         : $scoreCutOff => Cut-off for conserved or not

    my ($argHashRef) = @_;
    
    ## Flatten argument(s)
    my $arrayRef;
    my $regexp;
    my $scoreCutOff;
    
    my $tmpl = { 
	arrayRef => { required => 1, defined => 1, default => [], strict_type => 1, store => \$arrayRef},
	regexp => { required => 1, defined => 1, strict_type => 1, store => \$regexp},
	scoreCutOff => { required => 1, defined => 1, strict_type => 1, store => \$scoreCutOff},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my @allelValues;
    
    for my $element (@{$arrayRef}) {
	
	if ($element =~/$regexp/) {  #Find the key=value field
	    	
	    my @value = split(/=/, $element);  #Split key=value pair

	    my @tempArray = split(",", $value[1]); #Split on ","
	    
	    for my $allelValue (@tempArray) {

		if ($allelValue >= $scoreCutOff) {
		    
		    push(@allelValues, "Conserved");
		}
		else {

		    push(@allelValues, "NotConserved");
		}
	    }
	}
    }
    if (scalar(@allelValues) > 0) {
	
	return join(',', @allelValues);
    }
    else {

	return
    }
}


sub FindLCAF {

##FindLCAF
    
##Function : Adds the least common alternative allele frequency to each line
##Returns  : ""
##Arguments: $arrayRef, $regexp, $frequencyPosition
##         : $arrayRef          => The INFO array {REF}
##         : $regexp            => The regexp to used to locate correct ID field
##         : $frequencyPosition => The position to extract frequency from

    my ($argHashRef) = @_;

    ## Default(s)
    my $frequencyPosition;

    ## Flatten argument(s)
    my $arrayRef;
    my $regexp;
    
    my $tmpl = { 
	arrayRef => { required => 1, defined => 1, default => [], strict_type => 1, store => \$arrayRef},
	regexp => { required => 1, defined => 1, strict_type => 1, store => \$regexp},
	frequencyPosition => { default => 0,
			       allow => qr/^\d+$/,
			       strict_type => 1, store => \$frequencyPosition},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my $tempMaf;
    
    for my $element (@{$arrayRef}) {
	
	if ($element =~/$regexp/) {  #Find the key=value field
	    
	    my @value = split(/=/, $element);  #Split key=value pair

	    my @tempMafs = sort {$a <=> $b} grep { $_ ne "." } split(",", $value[1]); #Split on ",", remove entries containing only "." and sort remaining entries numerically ascending order
	    
	    if (scalar(@tempMafs) > 0) {
		
		## We are interested in the least common allele listed for this position. We cannot connect the frequency position in the list and the multiple alternative alleles. So the best we can do is report the least common allele frequency for multiple alternative allels. Unless the least common frequency is lower than the frequency defined as pathogenic for rare disease (usually 0.01) then this will work. In that case this will be a false positive, but it is better than taking the actual MAF which would be a false negative if the pathogenic variant found in the patient(s) has a lower frequency than the MAF.
		$tempMaf = $tempMafs[$frequencyPosition];   
	    }
	}
    }
    return $tempMaf
}


sub ParseMetaData {

##ParseMetaData
    
##Function : Writes metadata to filehandle specified by order in metaDataOrders.
##Returns  : ""
##Arguments: $metaDataHashRef, $metaDataString
##         : $metaDataHashRef => Hash for metaData {REF}
##         : $metaDataString  => The metaData string from vcf header

    my ($argHashRef) = @_;
    
    ## Flatten argument(s)
    my $metaDataHashRef;
    my $metaDataString;
    
    my $tmpl = { 
	metaDataHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$metaDataHashRef},
	metaDataString => { required => 1, defined => 1, strict_type => 1, store => \$metaDataString},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];
    
    if ($metaDataString=~/^##fileformat/) {  #Catch fileformat as it has to be at the top of header
	
	push(@{${$metaDataHashRef}{fileformat}{fileformat}}, $metaDataString);  #Save metadata string
    }
    elsif ($metaDataString=~/^##contig/) {  #catch contigs to not sort them later

	push(@{${$metaDataHashRef}{contig}{contig}}, $metaDataString);  #Save metadata string
    }
    elsif ($metaDataString=~/^##(\w+)=(\S+)/) {  #FILTER, FORMAT, INFO etc and more custom records
	
	push(@{${$metaDataHashRef}{$1}{$2}}, $metaDataString);  #Save metadata string
    }
    else {  #All oddities

	push(@{${$metaDataHashRef}{other}{other}}, $metaDataString);  #Save metadata string
    }
}


sub WriteMetaData {

##WriteMetaData
    
##Function : Writes metadata to filehandle specified by order in metaDataOrders.
##Returns  : ""
##Arguments: $metaDataHashRef, $FILEHANDLE, $SELECTFILEHANDLE
##         : $metaDataHashRef  => Hash for metaData {REF}
##         : $FILEHANDLE       => The filehandle to write to
##         : $SELECTFILEHANDLE => The filehandle to write to {Optional}

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $metaDataHashRef;
    my $FILEHANDLE;
    my $SELECTFILEHANDLE;

    my $tmpl = { 
	metaDataHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$metaDataHashRef},
	FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE},
	SELECTFILEHANDLE => { store => \$SELECTFILEHANDLE},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my @metaDataOrders = ("fileformat", "FILTER", "FORMAT", "INFO", "FIX_INFO", "contig", "Software");  #Determine order to print for standard records
    my @lines;

    for (my $lineCounter=0;$lineCounter<scalar(@metaDataOrders);$lineCounter++) {
	
	my $metaDataRecordRef = \$metaDataOrders[$lineCounter];  #Alias

	if (${$metaDataHashRef}{ $$metaDataRecordRef }) {  #MetaDataRecordExists
	    
	    if ($$metaDataRecordRef eq "contig") {  #Should not be sorted, but printed "as is"
		
		foreach my $line (@{${$metaDataHashRef}{$$metaDataRecordRef}{ $$metaDataRecordRef }}) {
		    
		    say $FILEHANDLE $line;
		    
		    if (defined($SELECTFILEHANDLE)) {
			
			say $SELECTFILEHANDLE $line;
		    }
		}
		delete(${$metaDataHashRef}{ $$metaDataRecordRef });  #Enable print of rest later
	    }	    
	    else {
		
		foreach my $line (sort( keys %{${$metaDataHashRef}{ $$metaDataRecordRef }})) {
		    
		    say $FILEHANDLE @{${$metaDataHashRef}{ $$metaDataRecordRef }{$line}};
		    
		    if (defined($SELECTFILEHANDLE)) {
			
			say $SELECTFILEHANDLE @{${$metaDataHashRef}{ $$metaDataRecordRef }{$line}};
		    }
		}
		if (defined($SELECTFILEHANDLE)) {
		    
		    foreach my $line (sort( keys %{${$metaDataHashRef}{Select}{ $$metaDataRecordRef }})) {
			
			say $SELECTFILEHANDLE @{${$metaDataHashRef}{Select}{ $$metaDataRecordRef }{$line}};
		    }
		}
		foreach my $line (sort( keys %{${$metaDataHashRef}{Range}{ $$metaDataRecordRef }})) {
		    
		    say $FILEHANDLE @{${$metaDataHashRef}{Range}{ $$metaDataRecordRef }{$line}};
		}
		delete(${$metaDataHashRef}{ $$metaDataRecordRef });  #Enable print of rest later

		if (${$metaDataHashRef}{Select}{ $$metaDataRecordRef }) {
		    
		    delete(${$metaDataHashRef}{Select}{ $$metaDataRecordRef });
		}
		if (${$metaDataHashRef}{Range}{ $$metaDataRecordRef }) {
		    
		    delete(${$metaDataHashRef}{Range}{ $$metaDataRecordRef });
		}
	    }
	}
    }
    for my $keys (keys %{$metaDataHashRef}) {
	
	for my $secondKey (keys %{${$metaDataHashRef}{$keys}}) {

	    if (ref(${$metaDataHashRef}{$keys}{$secondKey}) eq "HASH") {
		
		for my $line ( sort(keys %{${$metaDataHashRef}{$keys}{$secondKey}}) ) {

		    say $FILEHANDLE @{${$metaDataHashRef}{$keys}{$secondKey}{$line}};
		    
		    if (defined($SELECTFILEHANDLE)) {
			
			say $SELECTFILEHANDLE @{${$metaDataHashRef}{$keys}{$secondKey}{$line}};
		    }
		    delete ${$metaDataHashRef}{$keys}{$secondKey}{$line};
		}
	    }
	    elsif (ref(${$metaDataHashRef}{$keys}{$secondKey}) eq "ARRAY") {

		foreach my $element (@{${$metaDataHashRef}{$keys}{$secondKey}}) {

		    say $element;

		    if (defined($SELECTFILEHANDLE)) {
			
			say $SELECTFILEHANDLE $element;
		    }
		}
		delete ${$metaDataHashRef}{$keys}{$secondKey};
	    }
	    else {

		say $FILEHANDLE @{${$metaDataHashRef}{$keys}{$secondKey}};
		
		if (defined($SELECTFILEHANDLE)) {
		    
		    say $SELECTFILEHANDLE @{${$metaDataHashRef}{$keys}{$secondKey}};
		}
		delete ${$metaDataHashRef}{$keys}{$secondKey};
	    }
	}
    }
}


sub UniqElements {
    
##UniqElements
    
##Function : Collect unique elements from array reference and return array reference with unique elements
##Returns  : "array reference"
##Arguments: $arrayRef
##         : $arrayRef => The array whose elements are to be made distinct {REF}

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $arrayRef;
    
    my $tmpl = { 
	arrayRef => { required => 1, defined => 1, default => [], strict_type => 1, store => \$arrayRef},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    my %seen;

    return [ grep { !$seen{$_}++ } @{$arrayRef} ];  #For each element in array, see if seen before and only return list distinct elements  
}

sub AddFeatureFileMetaDataToVcf {
		
##AddFeatureFileMetaDataToVcf

##Function : Adds feature file meta data annotation headers to meta data hash.
##Returns  : ""
##Arguments: $metaDataHashRef, $vcfHeaderHashRef, $dataHashRef, $featureAnnotationCoulumnsArrayRef, $fileKey
##         : $metaDataHashRef                   => Vcf meta data {REF}
##         : $vcfHeaderHashRef                  => Vcf header meta data
##         : $dataHashRef                       => Range file hash {REF}
##         : $featureAnnotationCoulumnsArrayRef => Range columns to include {REF}
##         : $fileKey                           => Range file key used to seperate range file(s) i.e., select and range
    
    my ($argHashRef) = @_;
    
    ## Flatten argument(s)
    my $metaDataHashRef;
    my $vcfHeaderHashRef;
    my $dataHashRef;
    my $featureAnnotationColumnsArrayRef;
    my $fileKey;
    
    my $tmpl = { 
	metaDataHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$metaDataHashRef},
	vcfHeaderHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$vcfHeaderHashRef},
	dataHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$dataHashRef},
	featureAnnotationColumnsArrayRef => { required => 1, defined => 1, default => [], strict_type => 1, store => \$featureAnnotationColumnsArrayRef},
	fileKey => { required => 1, defined => 1, strict_type => 1, store => \$fileKey},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];

    if (@{$featureAnnotationColumnsArrayRef}) { #FeatureFile annotations
	
	for my $annotation (keys %{${$dataHashRef}{Present}}) {
	    
	    unless (defined(${$vcfHeaderHashRef}{INFO}{$annotation})) { #Unless INFO header is already present add to file 
		
		push(@{${$metaDataHashRef}{ $fileKey }{INFO}{$annotation}}, ${$dataHashRef}{Present}{$annotation}{INFO});  #Save specific featureFile INFO
	    }
	}
    }
}


sub CheckTerms {

##CheckTerms
    
##Function : Check that the found terms in the vcf corresond to known terms - otherwise croak and exit.
##Returns  : ""
##Arguments: $termHashRef, $termRef, $termName
##         : $termHashRef => The term hash {REF}
##         : $termRef     => The found term {REF}
##         : $termName    => The origin of the term i.e Sift, PolyPhen, SO

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $termHashRef;
    my $termRef;
    my $termName;

    my $tmpl = { 
	termHashRef => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$termHashRef},
	termRef => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$termRef},
	termName => { required => 1, defined => 1, strict_type => 1, store => \$termName},
    };
    
    if ( (defined($$termRef)) && ($$termRef ne "") ) {
	
	unless (exists(${$termHashRef}{ $$termRef })) {
	    
	    warn("Could not find ".$termName." term from vcf in corresponding hash. Update hash to contain term: '".$$termRef."'\n");
	    exit 1;
	}
    }   
}

sub Help {

##Help
    
##Function : Print help text and exit with supplied exit code
##Returns  : ""
##Arguments: $USAGE, $exitCode
##         : $USAGE    => Help text
##         : $exitCode => Exit code

    my ($argHashRef) = @_;

    ## Flatten argument(s)
    my $USAGE;
    my $exitCode;

    my $tmpl = { 
	USAGE => {required => 1, defined => 1, strict_type => 1, store => \$USAGE},
	exitCode => { default => 0, strict_type => 1, store => \$exitCode},
    };
    
    check($tmpl, $argHashRef, 1) or die qw[Could not parse arguments!];    
    
    say STDOUT $USAGE;
    exit $exitCode;
}
