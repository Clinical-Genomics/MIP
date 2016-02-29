#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
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

my ($infile, $parseVEP, $rangeFeatureFile, $selectFeatureFile, $selectFeatureMatchingColumn, $selectOutfile, $writeSoftwareTag, $padding) = ("", 0, 0, 0, "nocmdinput", "nocmdinput", 1, 5000);
my (@metaData, @selectMetaData, @rangeFeatureAnnotationColumns, @selectFeatureAnnotationColumns); 
my (%geneAnnotation, %consequenceSeverity, %rangeData, %selectData, %snpEffCmd, %tree, %metaData, %siftTerm, %polyPhenTerm);

my $vcfParserVersion = "1.2.6";

## Enables cmd "vcfParser.pl" to print usage help 
if(scalar(@ARGV) == 0) {

    print STDOUT $USAGE, "\n";
    exit;
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
	   'h|help' => sub { print STDOUT $USAGE, "\n"; exit;},  #Display help text
	   'v|version' => sub { print STDOUT "\nvcfParser.pl ".$vcfParserVersion, "\n\n"; exit;},  #Display version number
    );


## Basic flag option check
if ( (scalar(@rangeFeatureAnnotationColumns) == 0) && ($rangeFeatureFile ne 0) ) {
    
    print STDOUT "\n".$USAGE, "\n";
    print STDERR "\n", "Need to specify which feature column(s) to use with range file: ".$rangeFeatureFile." when annotating variants by using flag -rf_ac","\n\n";
    exit;
}
if ( ($selectFeatureMatchingColumn eq "nocmdinput") && ($selectFeatureFile ne 0) ) {
    
    print STDOUT "\n".$USAGE, "\n";
    print STDERR "\n", "Need to specify which feature column to use with select file: ".$selectFeatureFile." when selecting variants by using flag -sf_mc","\n\n";
    exit;
}
if ( ($selectOutfile eq "nocmdinput") && ($selectFeatureFile ne 0) ) {
    
    print STDOUT "\n".$USAGE, "\n";
    print STDERR "\n", "Need to specify which a select outfile to use when selecting variants by using flag -sof","\n\n";
    exit;
}

@rangeFeatureAnnotationColumns = split(/,/,join(',',@rangeFeatureAnnotationColumns)); #Enables comma separated annotation columns on cmd
@selectFeatureAnnotationColumns = split(/,/,join(',',@selectFeatureAnnotationColumns)); #Enables comma separated annotation columns on cmd

###
#MAIN
###

&DefineSelectData();

if ($rangeFeatureFile ne 0) {

    &ReadRangeFile({'infilePath' => $rangeFeatureFile,
		    'rangeCoulumnsArrayRef' => \@rangeFeatureAnnotationColumns,
		    'rangeDataHashRef' => \%rangeData,
		    'paddingRef' => \$padding,
		    'rangeFileKey' => "RangeFile",
		   });
}

if ($selectFeatureFile ne 0) {

    &ReadRangeFile({'infilePath' => $selectFeatureFile,
		    'selectFeatureMatchingColumn' => $selectFeatureMatchingColumn,
		    'rangeCoulumnsArrayRef' => \@selectFeatureAnnotationColumns,
		    'rangeDataHashRef' => \%selectData,
		    'paddingRef' => \$padding,
		    'rangeFileKey' => "SelectFile",
		   });
}

&DefineSnpEffAnnotations();
&DefineConsequenceSeverity();
&DefineSiftTerms();
&DefinePolyPhenTerms();

&ReadInfileVCF({'metaData' => \%metaData,
		'snpEffCmdHashRef' => \%snpEffCmd,
		'rangeDataHashRef' => \%rangeData,
		'selectDataHashRef' => \%selectData,
		'consequenceSeverityHashRef' => \%consequenceSeverity,
		'siftTermHashRef' => \%siftTerm,
		'polyPhenTermHashRef' => \%polyPhenTerm,
		'rangeFeatureAnnotationColumnsArrayRef' => \@rangeFeatureAnnotationColumns,
		'selectFeatureAnnotationColumnsArrayRef' => \@selectFeatureAnnotationColumns,
		'selectOutFilePath' => $selectOutfile,
		'writeSoftwareTag' => $writeSoftwareTag,
	       });

###
#Sub Routines
###

sub DefineSelectData {

##DefineSelectData

##Function : Defines arbitrary INFO fields based on headers in selectFile
##Returns  : ""
##Arguments: None

    $selectData{'SelectFile'}{'HGNC_symbol'}{'INFO'} = q?##INFO=<ID=HGNC_symbol,Number=.,Type=String,Description="The HGNC gene symbol">?;
    $selectData{'SelectFile'}{'Ensembl_gene_id'}{'INFO'} = q?##INFO=<ID=Ensembl_gene_id,Number=.,Type=String,Description="Ensembl gene identifier">?;
    $selectData{'SelectFile'}{'OMIM_morbid'}{'INFO'} = q?##INFO=<ID=OMIM_morbid,Number=.,Type=String,Description="OMIM morbid ID associated with gene(s)">?;
    $selectData{'SelectFile'}{'Phenotypic_disease_model'}{'INFO'} = q?##INFO=<ID=Phenotypic_disease_model,Number=.,Type=String,Description="Known disease gene(s) phenotype inheritance model">?;
    $selectData{'SelectFile'}{'Clinical_db_gene_annotation'}{'INFO'} = q?##INFO=<ID=Clinical_db_gene_annotation,Number=.,Type=String,Description="Gene disease group association">?;
    $selectData{'SelectFile'}{'Reduced_penetrance'}{'INFO'} = q?##INFO=<ID=Reduced_penetrance,Number=.,Type=String,Description="Pathogenic gene which can exhibit reduced penetrance">?;
    $selectData{'SelectFile'}{'Disease_associated_transcript'}{'INFO'} = q?##INFO=<ID=Disease_associated_transcript,Number=.,Type=String,Description="Known pathogenic transcript(s) for gene">?;
    $selectData{'SelectFile'}{'Ensembl_transcript_to_refseq_transcript'}{'INFO'} = q?##INFO=<ID=Ensembl_transcript_to_refseq_transcript,Number=.,Type=String,Description="The link between ensembl transcript and refSeq transcript IDs">?;
    $selectData{'SelectFile'}{'Gene_description'}{'INFO'} = q?##INFO=<ID=Gene_description,Number=.,Type=String,Description="The HGNC gene description">?;
    $selectData{'SelectFile'}{'Genetic_disease_model'}{'INFO'} = q?##INFO=<ID=Genetic_disease_model,Number=.,Type=String,Description="Known disease gene(s) inheritance model">?;
    $selectData{'SelectFile'}{'No_hgnc_symbol'}{'INFO'} = q?##INFO=<ID=No_hgnc_symbol,Number=.,Type=String,Description="Clinically relevant genetic regions lacking a HGNC_symbol or Ensembl gene ">?;
}

sub DefineSnpEffAnnotations {

##DefineSnpEffAnnotations

##Function : Defines the snpEff annotations that can be parsed and modified 
##Returns  : ""
##Arguments: None
    
    $snpEffCmd{'SnpEff'}{'Dbsnp129LCAF'}{'File'} = q?dbsnp_\S+.excluding_sites_after_129.vcf?;
    $snpEffCmd{'SnpEff'}{'Dbsnp129LCAF'}{'INFO'} = q?##INFO=<ID=Dbsnp129LCAF,Number=1,Type=Float,Description="Least common AF in dbSNP excluding sites after 129.">?;
    $snpEffCmd{'SnpEff'}{'Dbsnp129LCAF'}{'FIX_INFO'} = q?##INFO=<ID=SnpSift_CAF,Number=.,Type=String,Description="An ordered, comma delimited list of allele frequencies based on 1000Genomes, starting with the reference allele followed by alternate alleles as ordered in the ALT column. Where a 1000Genomes alternate allele is not in the dbSNPs alternate allele set, the allele is added to the ALT column.  The minor allele is the second largest value in the list, and was previuosly reported in VCF as the GMAF.  This is the GMAF reported on the RefSNP and EntrezSNP pages and VariationReporter">?;

    $snpEffCmd{'SnpEff'}{'DbsnpLCAF'}{'File'} = q?dbsnp_\d+.\w\d+.vcf?;
    $snpEffCmd{'SnpEff'}{'DbsnpLCAF'}{'INFO'} = q?##INFO=<ID=DbsnpLCAF,Number=1,Type=Float,Description="Least common AF in the DbSNP database.">?;
    $snpEffCmd{'SnpEff'}{'DbsnpLCAF'}{'FIX_INFO'} = q?##INFO=<ID=SnpSift_CAF,Number=.,Type=String,Description="An ordered, comma delimited list of allele frequencies based on 1000Genomes, starting with the reference allele followed by alternate alleles as ordered in the ALT column. Where a 1000Genomes alternate allele is not in the dbSNPs alternate allele set, the allele is added to the ALT column.  The minor allele is the second largest value in the list, and was previuosly reported in VCF as the GMAF.  This is the GMAF reported on the RefSNP and EntrezSNP pages and VariationReporter">?;

    $snpEffCmd{'SnpEff'}{'1000GAF'}{'File'} = q?1000G_phase\d+.\S+.vcf|ALL.wgs.phase\d+.\S+.vcf?;
    $snpEffCmd{'SnpEff'}{'1000GAF'}{'INFO'} = q?##INFO=<ID=1000GAF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1) in the 1000G database.">?;
    $snpEffCmd{'SnpEff'}{'1000GAF'}{'FIX_INFO'} = q?##INFO=<ID=SnpSift_AF,Number=.,Type=String,Description="Estimated allele frequency in the range (0,1)">?;

    $snpEffCmd{'SnpEff'}{'ESPMAF'}{'File'} = q?ESP\d+SI-V\d+-\w+.updatedProteinHgvs.snps_indels.vcf?;
    $snpEffCmd{'SnpEff'}{'ESPMAF'}{'INFO'} = q?##INFO=<ID=ESPMAF,Number=1,Type=Float,Description="Global Minor Allele Frequency in the ESP database.">?;
    $snpEffCmd{'SnpEff'}{'ESPMAF'}{'FIX_INFO'} = q?##INFO=<ID=SnpSift_MAF,Number=.,Type=String,Description="Minor Allele Frequency in percent in the order of EA,AA,All">?;

    $snpEffCmd{'SnpEff'}{'EXACAF'}{'File'} = q?ExAC.r\d+.\d+.sites.vep.vcf?;
    $snpEffCmd{'SnpEff'}{'EXACAF'}{'INFO'} = q?##INFO=<ID=EXACAF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1) in Exac">?;
    $snpEffCmd{'SnpEff'}{'EXACAF'}{'FIX_INFO'} = q?##INFO=<ID=SnpSift_AF,Number=.,Type=String,Description="Estimated allele frequency in the range (0,1)">?;

    $snpEffCmd{'SnpEff'}{'EXACMAXAF'}{'File'} = q?ExAC.r\d+.\d+.sites.vep.vcf?;
    $snpEffCmd{'SnpEff'}{'EXACMAXAF'}{'INFO'} = q?##INFO=<ID=EXACMAXAF,Number=A,Type=Float,Description="Estimated max allele frequency in the range (0,1) in Exac">?;
    $snpEffCmd{'SnpEff'}{'EXACMAXAF'}{'FIX_INFO'} = q?##INFO=<ID=SnpSift_AF,Number=.,Type=String,Description="Estimated max allele frequency in the range (0,1)">?;

    $snpEffCmd{'SnpEff'}{'CLNSIG'}{'File'} = q?clinvar_\d+.vcf?;
    $snpEffCmd{'SnpEff'}{'CLNSIG'}{'INFO'} = q?##INFO=<ID=CLNSIG,Number=A,Type=String,Description="Variant Clinical Significance, 0 - Uncertain significance, 1 - not provided, 2 - Benign, 3 - Likely benign, 4 - Likely pathogenic, 5 - Pathogenic, 6 - drug response, 7 - histocompatibility, 255 - other">?;
    $snpEffCmd{'SnpEff'}{'CLNSIG'}{'FIX_INFO'} = q?##INFO=<ID=SnpSift_CLNSIG,Number=A,Type=String,Description="Variant Clinical Significance, 0 - Uncertain significance, 1 - not provided, 2 - Benign, 3 - Likely benign, 4 - Likely pathogenic, 5 - Pathogenic, 6 - drug response, 7 - histocompatibility, 255 - other">?;

    $snpEffCmd{'SnpEff'}{'CLNACC'}{'File'} = q?clinvar_\d+.vcf?;
    $snpEffCmd{'SnpEff'}{'CLNACC'}{'INFO'} = q?##INFO=<ID=CLNACC,Number=.,Type=String,Description="Variant Accession and Versions">?;
    $snpEffCmd{'SnpEff'}{'CLNACC'}{'FIX_INFO'} = q?##INFO=<ID=SnpSift_CLNACC,Number=.,Type=String,Description="Variant Accession and Versions">?;

    $snpEffCmd{'SnpEff'}{'phastCons100way_vertebrate_prediction_term'}{'File'} = q?SnpSift dbnsfp?;
    $snpEffCmd{'SnpEff'}{'phastCons100way_vertebrate_prediction_term'}{'INFO'} = q?##INFO=<ID=phastCons100way_vertebrate_prediction_term,Number=A,Type=String,Description="PhastCons conservation prediction term">?;

    $snpEffCmd{'SnpEff'}{'phyloP100way_vertebrate_prediction_term'}{'File'} = q?SnpSift dbnsfp?;
    $snpEffCmd{'SnpEff'}{'phyloP100way_vertebrate_prediction_term'}{'INFO'} = q?##INFO=<ID=phyloP100way_vertebrate_prediction_term,Number=A,Type=String,Description="PhyloP conservation prediction term">?;

    $snpEffCmd{'SnpEff'}{'GERP++_RS_prediction_term'}{'File'} = q?SnpSift dbnsfp?;
    $snpEffCmd{'SnpEff'}{'GERP++_RS_prediction_term'}{'INFO'} = q?##INFO=<ID=GERP++_RS_prediction_term,Number=A,Type=String,Description="GERP RS conservation prediction term">?;

    $snpEffCmd{'SnpEff'}{'MTAF'}{'File'} = q?genbank_haplogroup_\d+_\S+.vcf?;
    $snpEffCmd{'SnpEff'}{'MTAF'}{'INFO'} = q?##INFO=<ID=MTAF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">?;
    $snpEffCmd{'SnpEff'}{'MTAF'}{'FIX_INFO'} = q?##INFO=<ID=SnpSift_MTAF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">?;

}

sub DefineConsequenceSeverity {

##DefineConsequenceSeverity

##Function : Defines the precedence of consequences for SO-terms
##Returns  : ""
##Arguments: None

    $consequenceSeverity{'transcript_ablation'}{'Rank'} = 1;
    $consequenceSeverity{'transcript_ablation'}{'GeneticRegionAnnotation'} = "exonic";
    $consequenceSeverity{'splice_donor_variant'}{'Rank'} = 2;
    $consequenceSeverity{'splice_donor_variant'}{'GeneticRegionAnnotation'} = "splicing";
    $consequenceSeverity{'splice_acceptor_variant'}{'Rank'} = 2;
    $consequenceSeverity{'splice_acceptor_variant'}{'GeneticRegionAnnotation'} = "splicing";
    $consequenceSeverity{'stop_gained'}{'Rank'} = 3;
    $consequenceSeverity{'stop_gained'}{'GeneticRegionAnnotation'} = "exonic";
    $consequenceSeverity{'frameshift_variant'}{'Rank'} = 4;
    $consequenceSeverity{'frameshift_variant'}{'GeneticRegionAnnotation'} = "exonic";
    $consequenceSeverity{'stop_lost'}{'Rank'} = 5;
    $consequenceSeverity{'stop_lost'}{'GeneticRegionAnnotation'} = "exonic";
    $consequenceSeverity{'start_lost'}{'Rank'} = 5;
    $consequenceSeverity{'start_lost'}{'GeneticRegionAnnotation'} = "exonic";
    $consequenceSeverity{'initiator_codon_variant'}{'Rank'} = 6;
    $consequenceSeverity{'initiator_codon_variant'}{'GeneticRegionAnnotation'} = "exonic";
    $consequenceSeverity{'inframe_insertion'}{'Rank'} = 6;
    $consequenceSeverity{'inframe_insertion'}{'GeneticRegionAnnotation'} = "exonic";
    $consequenceSeverity{'inframe_deletion'}{'Rank'} = 6;
    $consequenceSeverity{'inframe_deletion'}{'GeneticRegionAnnotation'} = "exonic";
    $consequenceSeverity{'missense_variant'}{'Rank'} = 6;
    $consequenceSeverity{'missense_variant'}{'GeneticRegionAnnotation'} = "exonic";
    $consequenceSeverity{'protein_altering_variant'}{'Rank'} = 6;
    $consequenceSeverity{'protein_altering_variant'}{'GeneticRegionAnnotation'} = "exonic";
    $consequenceSeverity{'transcript_amplification'}{'Rank'} = 7;
    $consequenceSeverity{'transcript_amplification'}{'GeneticRegionAnnotation'} = "exonic";
    $consequenceSeverity{'splice_region_variant'}{'Rank'} = 8;
    $consequenceSeverity{'splice_region_variant'}{'GeneticRegionAnnotation'} = "splicing";
    $consequenceSeverity{'incomplete_terminal_codon_variant'}{'Rank'} = 9;
    $consequenceSeverity{'incomplete_terminal_codon_variant'}{'GeneticRegionAnnotation'} = "exonic";
    $consequenceSeverity{'synonymous_variant'}{'Rank'} = 10;
    $consequenceSeverity{'synonymous_variant'}{'GeneticRegionAnnotation'} = "exonic";
    $consequenceSeverity{'stop_retained_variant'}{'Rank'} = 10;
    $consequenceSeverity{'stop_retained_variant'}{'GeneticRegionAnnotation'} = "exonic";
    $consequenceSeverity{'coding_sequence_variant'}{'Rank'} = 11;
    $consequenceSeverity{'coding_sequence_variant'}{'GeneticRegionAnnotation'} = "exonic";
    $consequenceSeverity{'mature_miRNA_variant'}{'Rank'} = 12;
    $consequenceSeverity{'mature_miRNA_variant'}{'GeneticRegionAnnotation'} = "ncRNA_exonic";
    $consequenceSeverity{'5_prime_UTR_variant'}{'Rank'} = 13;
    $consequenceSeverity{'5_prime_UTR_variant'}{'GeneticRegionAnnotation'} = "5UTR";
    $consequenceSeverity{'3_prime_UTR_variant'}{'Rank'} = 14;
    $consequenceSeverity{'3_prime_UTR_variant'}{'GeneticRegionAnnotation'} = "3UTR";
    $consequenceSeverity{'non_coding_transcript_exon_variant'}{'Rank'} = 15;
    $consequenceSeverity{'non_coding_transcript_exon_variant'}{'GeneticRegionAnnotation'} = "ncRNA_exonic";
    $consequenceSeverity{'non_coding_transcript_variant'}{'Rank'} = 15;
    $consequenceSeverity{'non_coding_transcript_variant'}{'GeneticRegionAnnotation'} = "ncRNA";
    $consequenceSeverity{'intron_variant'}{'Rank'} = 16;
    $consequenceSeverity{'intron_variant'}{'GeneticRegionAnnotation'} = "intronic";
    $consequenceSeverity{'NMD_transcript_variant'}{'Rank'} = 17;
    $consequenceSeverity{'NMD_transcript_variant'}{'GeneticRegionAnnotation'} = "ncRNA";
    $consequenceSeverity{'upstream_gene_variant'}{'Rank'} = 18;
    $consequenceSeverity{'upstream_gene_variant'}{'GeneticRegionAnnotation'} = "upstream";
    $consequenceSeverity{'downstream_gene_variant'}{'Rank'} = 19;
    $consequenceSeverity{'downstream_gene_variant'}{'GeneticRegionAnnotation'} = "downstream";
    $consequenceSeverity{'TFBS_ablation'}{'Rank'} = 20;
    $consequenceSeverity{'TFBS_ablation'}{'GeneticRegionAnnotation'} = "TFBS";
    $consequenceSeverity{'TFBS_amplification'}{'Rank'} = 21;
    $consequenceSeverity{'TFBS_amplification'}{'GeneticRegionAnnotation'} = "TFBS";
    $consequenceSeverity{'TF_binding_site_variant'}{'Rank'} = 22;
    $consequenceSeverity{'TF_binding_site_variant'}{'GeneticRegionAnnotation'} = "TFBS";
    $consequenceSeverity{'regulatory_region_variant'}{'Rank'} = 22;
    $consequenceSeverity{'regulatory_region_variant'}{'GeneticRegionAnnotation'} = "regulatory_region";
    $consequenceSeverity{'regulatory_region_ablation'}{'Rank'} = 23;
    $consequenceSeverity{'regulatory_region_ablation'}{'GeneticRegionAnnotation'} = "regulatory_region";
    $consequenceSeverity{'regulatory_region_amplification'}{'Rank'} = 24;
    $consequenceSeverity{'regulatory_region_amplification'}{'GeneticRegionAnnotation'} = "regulatory_region";
    $consequenceSeverity{'feature_elongation'}{'Rank'} = 25;
    $consequenceSeverity{'feature_elongation'}{'GeneticRegionAnnotation'} = "genomic_feature";
    $consequenceSeverity{'feature_truncation'}{'Rank'} = 26;
    $consequenceSeverity{'feature_truncation'}{'GeneticRegionAnnotation'} = "genomic_feature";
    $consequenceSeverity{'intergenic_variant'}{'Rank'} = 27;
    $consequenceSeverity{'intergenic_variant'}{'GeneticRegionAnnotation'} = "intergenic" 

}


sub DefineSiftTerms {

##DefineSiftTerms

##Function : Defines the allowed sift terms
##Returns  : ""
##Arguments: None

    $siftTerm{'deleterious'} = "";
    $siftTerm{'tolerated'} = "";
    $siftTerm{'tolerated_low_confidence'} = "";
    $siftTerm{'deleterious_low_confidence'} = "";

}

sub DefinePolyPhenTerms {

##DefinePolyPhenTerms

##Function : Defines the allowed polyPhen terms
##Returns  : ""
##Arguments: None

    $polyPhenTerm{'probably_damaging'} = "";
    $polyPhenTerm{'possibly_damaging'} = "";
    $polyPhenTerm{'benign'} = "";
    $polyPhenTerm{'unknown'} = "";

}

sub ReadRangeFile {

##ReadRangeFile

##Function : Reads a file containg features to be annotated using range queries e.g. EnsemblGeneID. Adds to Metadata hash and creates Interval tree for feature.
##Returns  : ""
##Arguments: $rangeDataHashRef, $rangeCoulumnsArrayRef, $rangeFileKey, $infilePath, $paddingRef, $selectFeatureMatchingColumn
##         : $rangeDataHashRef            => Range file hash {REF}
##         : $rangeCoulumnsArrayRef       => Range columns to include {REF}
##         : $rangeFileKey                => Range file key used to seperate range file(s) i.e., select and range
##         : $infilePath                  => Infile path
##         : $paddingRef                  => Padding distance {REF}
##         : $selectFeatureMatchingColumn => Column in the select file to match with vcf key annotation

    my ($argHashRef) = @_;
    
    ## Flatten argument(s)
    my $rangeDataHashRef = ${$argHashRef}{'rangeDataHashRef'};
    my $rangeCoulumnsArrayRef = ${$argHashRef}{'rangeCoulumnsArrayRef'};

    my @headers; #Save headers from rangeFile

    open(RRF, "<".${$argHashRef}{'infilePath'}) or die "Can't open ".${$argHashRef}{'infilePath'}.":".$!, "\n"; 

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

	    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar(@{$rangeCoulumnsArrayRef});$extractColumnsCounter++) { #Defines what scalar to store
	
		&AddMetaDataINFO($rangeDataHashRef, ${$argHashRef}{'rangeFileKey'}, \$headers[ $$rangeCoulumnsArrayRef[$extractColumnsCounter] ], \$extractColumnsCounter, \${$argHashRef}{'infilePath'});
	    }
	    next;
	}
	if ( $_ =~/^(\S+)/ ) {	
	    
	    my @lineElements = split("\t",$_); #Loads range file line elements

	    if (defined(${$argHashRef}{'selectFeatureMatchingColumn'}) ) {

		$lineElements[ ${$argHashRef}{'selectFeatureMatchingColumn'} ] =~ s/\s/_/g; # Replace whitespace with "_"
		$selectData{$lineElements[ ${$argHashRef}{'selectFeatureMatchingColumn'} ]} = $lineElements[ ${$argHashRef}{'selectFeatureMatchingColumn'} ];
	    }

	    ## Create Interval Tree
	    if (scalar(@{$rangeCoulumnsArrayRef}) > 0) {#Annotate vcf with features from range file
		
		&RangeAnnotations(\@{$rangeCoulumnsArrayRef}, \@lineElements, $rangeDataHashRef, ${$argHashRef}{'rangeFileKey'}, \@headers, ${$argHashRef}{'paddingRef'});
	    }
	}
    }
    close(RRF);
    print STDERR "Finished Reading ".${$argHashRef}{'rangeFileKey'}." file: ".${$argHashRef}{'infilePath'},"\n";
}


sub ReadInfileVCF {

##ReadInfileVCF

##Function : Reads infile in vcf format and adds and parses annotations
##Returns  : ""
##Arguments: $metaDataHashRef, $snpEffCmdHashRef, $rangeDataHashRef, $selectDataHashRef, $consequenceSeverityHashRef, $siftTermHashRef, $polyPhenTermHashRef, $rangeFeatureAnnotationColumnsArrayRef, $selectOutFilePath 
##         : $metaDataHashRef                        => Vcf meta data {REF}
##         : $snpEffCmdHashRef                       => SnpEff meta data {REF}
##         : $rangeDataHashRef                       => Range file data {REF}
##         : $selectDataHashRef                      => Select file data {REF}
##         : $consequenceSeverityHashRef             => Consequence severity for SO-terms {REF}
##         : $siftTermHashRef                        => Sift prediction terms {REF}
##         : $polyPhenTermHashRef                    => PolyPhen prediction terms {REF}
##         : $rangeFeatureAnnotationColumnsArrayRef  => Range feature columns {REF}
##         : $selectFeatureAnnotationColumnsArrayRef => Select feature columns {REF}
##         : writeSoftwareTag                        => Write software tag to vcf header switch
##         : $selectOutFilePath                      => The select file path

    my ($argHashRef) = @_;
    
    my %default = ('writeSoftwareTag' => 1,
	);

    &SetDefaultArg(\%{$argHashRef}, \%default);

    ## Flatten argument(s)
    my $metaDataHashRef = ${$argHashRef}{'metaDataHashRef'};
    my $snpEffCmdHashRef = ${$argHashRef}{'snpEffCmdHashRef'};
    my $rangeDataHashRef = ${$argHashRef}{'rangeDataHashRef'};
    my $selectDataHashRef = ${$argHashRef}{'selectDataHashRef'};
    my $consequenceSeverityHashRef = ${$argHashRef}{'consequenceSeverityHashRef'};
    my $rangeFeatureAnnotationColumnsArrayRef = ${$argHashRef}{'rangeFeatureAnnotationColumnsArrayRef'};
    my $selectFeatureAnnotationColumnsArrayRef = ${$argHashRef}{'selectFeatureAnnotationColumnsArrayRef'};
    my $siftTermHashRef = ${$argHashRef}{'siftTermHashRef'};
    my $polyPhenTermHashRef = ${$argHashRef}{'polyPhenTermHashRef'};

    my @vepFormatField;
    my %vepFormatFieldColumn;

    my @printTSVFields = ("FeatureType");
    my @featureFields;

    my %vcfHeader;

    if ($selectFeatureFile ne 0) {

	open(WOSFTSV, ">".${$argHashRef}{'selectOutFilePath'}) or die "Can't open ".${$argHashRef}{'selectOutFilePath'}.":".$!, "\n";
    }
    
    while (<>) {
	
	chomp $_;  # Remove newline
	
	if (m/^\s+$/) {	# Avoid blank lines
	    next;
	}
	if ($_=~/^##(\S+)=/) {  # MetaData

	    &ParseMetaData(\%{$metaDataHashRef}, $_);

	    if ($_=~/INFO\=\<ID\=(\w+)/) { # Collect all INFO keys
	    
		$vcfHeader{'INFO'}{$1} = $1; #Save to hash
	    }
	    if ($_=~/SnpSiftCmd\=/) { #Find SnpEff command meta line

		for my $database (keys %{${$snpEffCmdHashRef}{'SnpEff'}}) {

		    if ($_=~/${$snpEffCmdHashRef}{'SnpEff'}{$database}{'File'}/) { #SnpEff/Sift has been used to annotate input vcf
			
			unless (defined($vcfHeader{'INFO'}{$database})) { #Unless INFO header is already present add to metaDataHeader
			    
			    ${$snpEffCmdHashRef}{'Present'}{'database'}{$database} = $database; #Save which frequency db has been used for later
			    push(@{${$metaDataHashRef}{'INFO'}{$database}}, ${$snpEffCmdHashRef}{'SnpEff'}{$database}{'INFO'});

			    if (defined(${$snpEffCmdHashRef}{'SnpEff'}{$database}{'FIX_INFO'})) { #If 'FIX_INFO' flag is present add to metaDataHeader

				push(@{${$metaDataHashRef}{'FIX_INFO'}{$database}}, ${$snpEffCmdHashRef}{'SnpEff'}{$database}{'FIX_INFO'});
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
		if ($parseVEP == 1) {
			
		    if ( ($vepFormatFieldColumn{'SYMBOL'}) && ($vepFormatFieldColumn{'HGVSc'}) && ($vepFormatFieldColumn{'HGVSp'})) {
			
			push(@{${$metaDataHashRef}{'INFO'}{'HGVScp'}}, '##INFO=<ID=HGVScp,Number=.,Type=String,Description="Transcript and protein functional annotation.">');
			push(@{${$metaDataHashRef}{'INFO'}{'MostSevereConsequence'}}, '##INFO=<ID=MostSevereConsequence,Number=.,Type=String,Description="Most severe genomic consequence.">');
			push(@{${$metaDataHashRef}{'INFO'}{'GeneticRegionAnnotation'}}, '##INFO=<ID=GeneticRegionAnnotation,Number=.,Type=String,Description="Genetic region that variant falls into.">');
			
		    }
		    if ($vepFormatFieldColumn{'SIFT'}) {
			
			push(@{${$metaDataHashRef}{'INFO'}{'Sift'}}, '##INFO=<ID=Sift,Number=.,Type=String,Description="Sift protein function prediction term">');
		    }
		    if ($vepFormatFieldColumn{'PolyPhen'}) {
			
			push(@{${$metaDataHashRef}{'INFO'}{'PolyPhen'}}, '##INFO=<ID=PolyPhen,Number=.,Type=String,Description="PolyPhen protein function prediction term">');
		    }
		}
		next;
	    }
	    next;
	}
	if ($_=~/^#CHROM/) {

	    &AddRangeMetaDataToVcf({'metaDataHashRef' => $metaDataHashRef,
				    'vcfHeaderHashRef' => \%vcfHeader,
				    'featureAnnotationColumnsArrayRef' => $rangeFeatureAnnotationColumnsArrayRef,
				    'dataHashRef' => $rangeDataHashRef,
				    'fileKey' => "Range",
				   });
	    &AddRangeMetaDataToVcf({'metaDataHashRef' => $metaDataHashRef,
				    'vcfHeaderHashRef' => \%vcfHeader,
				    'featureAnnotationColumnsArrayRef' => $selectFeatureAnnotationColumnsArrayRef,
				    'dataHashRef' => $selectDataHashRef,
				    'fileKey' => "Select",
				   });

	    if (${$argHashRef}{'writeSoftwareTag'} == 1) {

		&AddProgramToMeta(\%{$metaDataHashRef});
	    }
	    if (scalar(@{$selectFeatureAnnotationColumnsArrayRef}) > 0) { #SelectFile annotations

		&WriteMetaData(\%{$metaDataHashRef}, *STDOUT, *WOSFTSV);
		print STDOUT $_, "\n";  #Write header line
		print WOSFTSV $_, "\n";  #Write header line
	    }
	    else {

		&WriteMetaData(\%{$metaDataHashRef}, *STDOUT);
		print STDOUT $_, "\n";  #Write header line
	    }
	    
	    if ($parseVEP == 1) {

		@featureFields = ("MostSevereConsequence", "GeneticRegionAnnotation");
		
		if ($vepFormatFieldColumn{'SIFT'}) {
		    
		    push(@featureFields, "Sift");
		}
		if ($vepFormatFieldColumn{'PolyPhen'}) {
		    
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
	    for my $database (keys % {${$snpEffCmdHashRef}{'Present'}{'database'}}) { #Note that the vcf should only contain 1 database entry

		if ( ($database eq "Dbsnp129LCAF") || ($database eq "DbsnpLCAF") ) {

		    my @tempArray = split(/;/, $lineElements[7]);  #Split INFO field to key=value items

		     my $tempMafList = &FindAF(\@tempArray, "\\S+_CAF=");

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
		elsif($database eq "1000GAF") {
			
		    my @tempArray = split(/;/, $lineElements[7]);  #Split INFO field to key=value items

		    my $tempMaf = &FindAF(\@tempArray, "\\S+_1000GAF_AF=");

		    if (defined($tempMaf)) {
			
			## Save Alternative Allele frequency info
			$variantLine .= $database."=".$tempMaf.";";
			$selectedVariantLine .= $database."=".$tempMaf.";";
		    }
		}
		elsif($database eq "ESPMAF") {
		    
		    my @tempArray = split(/;/, $lineElements[7]);  #Split INFO field to key=value items

		    my $tempMaf = &FindLCAF(\@tempArray, "\\S+_MAF=", "2");
	
		    if (defined($tempMaf)) {
			
			$tempMaf = $tempMaf / 100; #fraction for consistent representation
			
			## Save MAF frequency info  
			$variantLine .= $database."=".$tempMaf.";";
			$selectedVariantLine .= $database."=".$tempMaf.";";
		    }   
		}
		elsif($database eq "EXACAF") {
		    
		    my @tempArray = split(/;/, $lineElements[7]);  #Split INFO field to key=value items
		    
		    my $tempMaf = &FindAF(\@tempArray, "\\S+_EXACAF_AF=");
		    
		    if (defined($tempMaf)) {
			
			## Save Alternative Allele frequency info  
			$variantLine .= $database."=".$tempMaf.";";
			$selectedVariantLine .= $database."=".$tempMaf.";";
		    }		    
		}
		elsif($database eq "EXACMAXAF") {
		    
		    my @tempArray = split(/;/, $lineElements[7]);  #Split INFO field to key=value items
		    
		    my $tempMaf = &FindAF(\@tempArray, "\\S+_EXACAF_MAX_AF=");
		    
		    if (defined($tempMaf)) {
			
			## Save Alternative Allele frequency info  
			$variantLine .= $database."=".$tempMaf.";";
			$selectedVariantLine .= $database."=".$tempMaf.";";
		    }		    
		}
		elsif($database eq "MTAF") {
			
		    my @tempArray = split(/;/, $lineElements[7]);  #Split INFO field to key=value items

		    my $tempMaf = &FindAF(\@tempArray, "\\S+_MTAF=");

		    if (defined($tempMaf)) {
			
			## Save Alternative Allele frequency info
			$variantLine .= $database."=".$tempMaf.";";
			$selectedVariantLine .= $database."=".$tempMaf.";";
		    }
		}
		elsif($database eq "CLNSIG") {
		    
		    my @tempArray = split(/;/, $lineElements[7]);  #Split INFO field to key=value items
		    
		    my $tempMaf = &FindAF(\@tempArray, "\\S+_CLNSIG=");
		    
		    if (defined($tempMaf)) {
			
			## Save Alternative Allele frequency info  
			$variantLine .= $database."=".$tempMaf.";";
			$selectedVariantLine .= $database."=".$tempMaf.";";
		    }
		}
		elsif($database eq "CLNACC") {
		    
		    my @tempArray = split(/;/, $lineElements[7]);  #Split INFO field to key=value items
		    
		    my $tempMaf = &FindAF(\@tempArray, "\\S+_CLNACC=");
		    
		    if (defined($tempMaf)) {
			
			## Save Alternative Allele frequency info  
			$variantLine .= $database."=".$tempMaf.";";
			$selectedVariantLine .= $database."=".$tempMaf.";";
		    }
		}
		elsif($database eq "phastCons100way_vertebrate_prediction_term") {
		    
		    my @tempArray = split(/;/, $lineElements[7]);  #Split INFO field to key=value items
		    
		    my $conservationTerm = &FindConserved(\@tempArray, "\\S+_phastCons100way_vertebrate=", 0.8);

		    if (defined($conservationTerm)) {
			
			## Save database info  
			$variantLine .= $database."=".$conservationTerm.";";
			$selectedVariantLine .= $database."=".$conservationTerm.";";
		    }
		}
		elsif($database eq "phyloP100way_vertebrate_prediction_term") {
		    
		    my @tempArray = split(/;/, $lineElements[7]);  #Split INFO field to key=value items
		    
		    my $conservationTerm = &FindConserved(\@tempArray, "\\S+_phyloP100way_vertebrate=", 2.5);

		    if (defined($conservationTerm)) {
			
			## Save database info
			$variantLine .= $database."=".$conservationTerm.";";
			$selectedVariantLine .= $database."=".$conservationTerm.";";
		    }
		}
		elsif($database eq "GERP++_RS_prediction_term") {
		    
		    
		    my @tempArray = split(/;/, $lineElements[7]);  #Split INFO field to key=value items
		    
		    my $conservationTerm = &FindConserved(\@tempArray, "\\S+_GERP\\+\\+_RS=", 2);

		    if (defined($conservationTerm)) {
			
			## Save database info
			$variantLine .= $database."=".$conservationTerm.";";
			$selectedVariantLine .= $database."=".$conservationTerm.";";
		    }
		}
	    }
	    %noIDRegion = &TreeAnnotations("SelectFile", \@lineElements, $selectDataHashRef, \$selectedVariantLine);  #Only for selectfile since all variants are passed to research file
	    &TreeAnnotations("RangeFile", \@lineElements, $rangeDataHashRef, \$variantLine);
	    
	    my @variantEffects = split(/;/, $lineElements[7]); #Split INFO field

	    ##Check that we have an INFO field
	    unless ($lineElements[7]) {

		print STDERR "No INFO field at line number: ".$.."\n";
		print STDERR "Displaying malformed line: ".$_, "\n";
		exit 1;
	    }
	    my $CSQTranscripts;
	    
	    for (my $variantEffectCounter=0;$variantEffectCounter<scalar(@variantEffects);$variantEffectCounter++) {
		
		if ($parseVEP == 1) {
		    
		    ## Find CSQ field and extract transripts that belong to select genes
		    if ($variantEffects[$variantEffectCounter]=~/CSQ\=(\S+)/) { #Find CSQ
			
			$CSQTranscripts = $1;
			my @transcripts = split(/,/, $1); #Split into transcripts elements
						    
			for (my $fieldCounter=0;$fieldCounter<scalar(@transcripts);$fieldCounter++) { #CSQ field
				
			    my @transcriptsEffects = split(/\|/, $transcripts[$fieldCounter]); #Split in "|"
				
			    if ( (defined($transcriptsEffects[ $vepFormatFieldColumn{'Feature'} ])) && ($transcriptsEffects[ $vepFormatFieldColumn{'Feature'} ] ne "") ) {
				    
				if (defined($transcriptsEffects[ $vepFormatFieldColumn{'SYMBOL'} ])) { 
					
				    $variantData{'Symbol'} = $transcriptsEffects[ $vepFormatFieldColumn{'SYMBOL'} ];  #Save HGNC Symbol
					
				    if ($selectData{ $variantData{'Symbol'} }) { #Exists in selected Features

					if ($selectedTranscriptCounter > 0) {
						
					    $selectedVariantLine .= ",".$transcripts[$fieldCounter];
					}
					else {
						
					    $selectedVariantLine .= "CSQ=".$transcripts[$fieldCounter];
					}
					$selectedTranscriptCounter++;
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
			}
			if ( ($selectedTranscriptCounter == 0) && ($transcriptsCounter == 0) ) {

			    $variantLine .= "CSQ=".$CSQTranscripts;  #No transcript info
			}
			unless ($variantEffectCounter == scalar(@variantEffects)-1) {  #Unless this is the last feature
				
			    $variantLine .= ";";
			    $selectedVariantLine .= ";";
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
	    if ($parseVEP == 1) {
		
		unless ( ($selectedTranscriptCounter == 0) && ($transcriptsCounter == 0) ) {  #Info to add from CSQ field
		    
		    $selectedVariantLine .= ";";
		    $variantLine .= ";";
		}
		$transcriptsCounter = 0;
		$selectedTranscriptCounter = 0;

		if (defined($CSQTranscripts)) {

		    my @transcripts = split(/,/, $CSQTranscripts);
		    
		    for (my $fieldCounter=0;$fieldCounter<scalar(@transcripts);$fieldCounter++) { #CSQ field
			
			my @transcriptsEffects = split(/\|/, $transcripts[$fieldCounter]); #Split in "|"
			
			if ( (defined($transcriptsEffects[ $vepFormatFieldColumn{'Feature'} ]) ) && ($transcriptsEffects[ $vepFormatFieldColumn{'Feature'} ] ne "") && ($transcriptsEffects[ $vepFormatFieldColumn{'Feature'} ] !~/ENSR\d+/) ) { #All but regulatory regions, since there currently is no HGNC_Symbol annotated for them
			    
			    my $selectedTranscriptTracker = 0; #Track if any transcripts belong to selected features
			    
			    if (defined($transcriptsEffects[ $vepFormatFieldColumn{'SYMBOL'} ])) { #Save HGNC Symbol
				
				$variantData{'Symbol'} = $transcriptsEffects[ $vepFormatFieldColumn{'SYMBOL'} ];
				
				if (${$selectDataHashRef}{ $variantData{'Symbol'} }) { #Exists in selected Features
				    
				    $selectedTranscriptTracker = 1; #Record belongs to selected Features
				    my $alleleGeneEntry = $transcriptsEffects[ $vepFormatFieldColumn{'Allele'} ].":".$variantData{'Symbol'};
				   
				    &AddFieldToElementCounter(\$selectedTranscriptCounter, \$selectedVariantLine, ",", \$alleleGeneEntry, "HGVScp=");
				}
				## Always include all transcripts in research list
				my $alleleGeneEntry = $transcriptsEffects[ $vepFormatFieldColumn{'Allele'} ].":".$variantData{'Symbol'};
				    
				&AddFieldToElementCounter(\$transcriptsCounter, \$variantLine, ",", \$alleleGeneEntry, "HGVScp=");
			    }
			    &AddFieldToElement(\$selectedTranscriptTracker, \$selectedVariantLine, \$variantLine, ":", \$transcriptsEffects[ $vepFormatFieldColumn{'Feature'} ]); #Save transcript

			    if (defined($transcriptsEffects[ $vepFormatFieldColumn{'Feature_type'} ])) { #Save Feature_type
				
				if ($selectedTranscriptTracker == 1) {
				    
				    if ($selectedTranscriptCounter == 0) { #First selected gene
					
					$selectedVariantData{'FeatureType'} = $transcriptsEffects[ $vepFormatFieldColumn{'Feature_type'} ];
				    }
				    else {
					
					$selectedVariantData{'FeatureType'} .= ",".$transcriptsEffects[ $vepFormatFieldColumn{'Feature_type'} ];
				    }
				}
				## Always include all transcripts in research list				    
				if ($transcriptsCounter == 0) { #First Gene
				    
				    $variantData{'FeatureType'} = $transcriptsEffects[ $vepFormatFieldColumn{'Feature_type'} ];
				}
				else {
				    
				    $variantData{'FeatureType'} .= ",".$transcriptsEffects[ $vepFormatFieldColumn{'Feature_type'} ];
				}
			    }
			    if (defined($transcriptsEffects[ $vepFormatFieldColumn{'Consequence'} ])) { #Save Consequence
				
				&AddFieldToElement(\$selectedTranscriptTracker, \$selectedVariantLine, \$variantLine, ":", \$transcriptsEffects[ $vepFormatFieldColumn{'Consequence'} ]);
				my @consequences = split(/\&/, $transcriptsEffects[ $vepFormatFieldColumn{'Consequence'} ]); #Find "MostSevereConsequence
				my $allele = $transcriptsEffects[ $vepFormatFieldColumn{'Allele'} ];

				for (my $consequencesCounter=0;$consequencesCounter<scalar(@consequences);$consequencesCounter++) {
				    
				    &CheckTerms($consequenceSeverityHashRef, \$consequences[$consequencesCounter], "SO");
				    
				    if (defined($variantData{'Symbol'})) {

					if(defined($consequence{ $variantData{'Symbol'} }{$allele}{'Score'})) { #Compare to previous record
					
					    if (${$consequenceSeverityHashRef}{$consequences[$consequencesCounter]}{'Rank'} < $consequence{ $variantData{'Symbol'} }{$allele}{'Score'}) { #Collect most severe consequence
						
						&AddToConsequenceHash(\$consequence{ $variantData{'Symbol'} }{$allele}{'Score'}, \${$consequenceSeverityHashRef}{$consequences[$consequencesCounter]}{'Rank'});
						&AddToConsequenceHash(\$consequence{ $variantData{'Symbol'} }{$allele}{'GeneticRegionAnnotation'}, \${$consequenceSeverityHashRef}{$consequences[$consequencesCounter]}{'GeneticRegionAnnotation'});
						&AddToConsequenceHash(\$consequence{ $variantData{'Symbol'} }{$allele}{'MostSevereConsequence'}, \$consequences[$consequencesCounter]);
						
						if (defined($vepFormatFieldColumn{'SIFT'}) ) {

						    &CheckTerms($siftTermHashRef, \$transcriptsEffects[ $vepFormatFieldColumn{'SIFT'} ], "Sift");
						    &AddToConsequenceHash(\$consequence{ $variantData{'Symbol'} }{$allele}{'Sift'}, \$transcriptsEffects[ $vepFormatFieldColumn{'SIFT'} ]);
						}
						if (defined($vepFormatFieldColumn{'PolyPhen'}) ) {
						    
						    &CheckTerms($polyPhenTermHashRef, \$transcriptsEffects[ $vepFormatFieldColumn{'PolyPhen'} ], "polyPhen");
						    &AddToConsequenceHash(\$consequence{ $variantData{'Symbol'} }{$allele}{'PolyPhen'}, \$transcriptsEffects[ $vepFormatFieldColumn{'PolyPhen'} ]);
						}
					    }
					}
					else { #First pass
					
					    &AddToConsequenceHash(\$consequence{ $variantData{'Symbol'} }{$allele}{'Score'}, \${$consequenceSeverityHashRef}{$consequences[$consequencesCounter]}{'Rank'});
					    &AddToConsequenceHash(\$consequence{ $variantData{'Symbol'} }{$allele}{'GeneticRegionAnnotation'}, \${$consequenceSeverityHashRef}{$consequences[$consequencesCounter]}{'GeneticRegionAnnotation'});
					    &AddToConsequenceHash(\$consequence{ $variantData{'Symbol'} }{$allele}{'MostSevereConsequence'}, \$consequences[$consequencesCounter]);    
					    
					    if (defined($vepFormatFieldColumn{'SIFT'}) ) {
						
						&AddToConsequenceHash(\$consequence{ $variantData{'Symbol'} }{$allele}{'Sift'}, \$transcriptsEffects[ $vepFormatFieldColumn{'SIFT'} ]);
					    }
					    if (defined($vepFormatFieldColumn{'PolyPhen'}) ) {
						
						&AddToConsequenceHash(\$consequence{ $variantData{'Symbol'} }{$allele}{'PolyPhen'}, \$transcriptsEffects[ $vepFormatFieldColumn{'PolyPhen'} ]);
					    }
					}
				    }
				}
			    }
			    if ( (defined($transcriptsEffects[ $vepFormatFieldColumn{'STRAND'} ])) && ($transcriptsEffects[ $vepFormatFieldColumn{'STRAND'} ] ne "") ) { #Save strand 

				if ($transcriptsEffects[ $vepFormatFieldColumn{'STRAND'} ] == 1) { #Remap to "+" or "-"

				    $variantData{'Strand'} = "s.+";
				}
				else {
				    
				    $variantData{'Strand'} = "s.-";
				}
				&AddFieldToElement(\$selectedTranscriptTracker, \$selectedVariantLine, \$variantLine, ":", \$variantData{'Strand'});
			    }
			    if (defined($transcriptsEffects[ $vepFormatFieldColumn{'EXON'} ]) && $transcriptsEffects[ $vepFormatFieldColumn{'EXON'}] =~/^\d/) {
				
				$variantData{'Exon'} = "e.".$transcriptsEffects[ $vepFormatFieldColumn{'EXON'} ]; #Save exon number "X/Y"
				&AddFieldToElement(\$selectedTranscriptTracker, \$selectedVariantLine, \$variantLine, ":", \$variantData{'Exon'});
			    }
			    if (defined($transcriptsEffects[ $vepFormatFieldColumn{'INTRON'} ]) && $transcriptsEffects[ $vepFormatFieldColumn{'INTRON'} ] =~/^\d/) {
				
				$variantData{'Intron'} = "i.".$transcriptsEffects[ $vepFormatFieldColumn{'INTRON'} ]; #Save intron number "X/Y"
				&AddFieldToElement(\$selectedTranscriptTracker, \$selectedVariantLine, \$variantLine, ":", \$variantData{'Intron'});
			    }
			    if ( (defined($transcriptsEffects[ $vepFormatFieldColumn{'HGVSc'} ])) && ($transcriptsEffects[ $vepFormatFieldColumn{'HGVSc'} ] ne "")) { #Save HGVS cDNA change
				
				my @cDNAChanges = split(/:/, $transcriptsEffects[ $vepFormatFieldColumn{'HGVSc'} ]);
				&AddFieldToElement(\$selectedTranscriptTracker, \$selectedVariantLine, \$variantLine, ":", \$cDNAChanges[1]);
			    }
			    if ( (defined($transcriptsEffects[ $vepFormatFieldColumn{'HGVSp'} ])) && ($transcriptsEffects[ $vepFormatFieldColumn{'HGVSp'} ] ne "") ) {
				
				my @pAAChanges = split(/:/, $transcriptsEffects[ $vepFormatFieldColumn{'HGVSp'} ]);
				&AddFieldToElement(\$selectedTranscriptTracker, \$selectedVariantLine, \$variantLine, ":", \$pAAChanges[1]);
			    }
			    if ($selectedTranscriptTracker == 1) {
				
				$selectedTranscriptCounter++;
			    }
			    ## Always include all transcripts in research list	
			    $transcriptsCounter++;
			}
		    }
		    &CollectConsequenceGenes(\%consequence, $selectDataHashRef, \@featureFields, \$selectedVariantLine, \$variantLine);
		}
	    }

	    $selectedVariantLine .= "\t".$sampleIDInfo;
	    $variantLine .= "\t".$sampleIDInfo;	

	    if ($parseVEP == 1) {
		
		if ($selectedTranscriptCounter > 0) { #Write to selected file
		    
		    print WOSFTSV $selectedVariantLine, "\n";
		}
		if ($transcriptsCounter > 0) { #Write to transcript file
		    
		    if (%noIDRegion) {  #No HGNC_symbol or EnsemblGeneID, but clinically relevant
		
			print WOSFTSV $selectedVariantLine, "\n";
		    }
		    print STDOUT $variantLine, "\n";
		}
		elsif ( ($selectedTranscriptCounter == 0) && ($transcriptsCounter == 0) ) {

		    if (%noIDRegion) {  #No HGNC_symbol or EnsemblGeneID, but clinically relevant

			print WOSFTSV $selectedVariantLine, "\n";
		    }
		    print STDOUT $variantLine, "\n";
		}
	    }
	    else {

		print STDOUT $variantLine, "\n";	
	    }
	}
    }
    close(WOSFTSV);
    print STDERR "Finished Processing VCF", "\n";
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

    my $hashKeyRef = $_[0];
    my $valueRef = $_[1];
    
    if (defined($valueRef)) {
	
	$$hashKeyRef = $$valueRef;
    }
}

sub AddFieldToElementCounter {

##AddFieldToElementCounter
    
##Function : Adds a field to an element. Adjust addition depending on if field has been seen before.
##Returns  : ""
##Arguments: $transcriptCounterRef, $lineRef, $separator, $valueRef
##         : $transcriptCounterRef => The transcript counter {REF}
##         : $lineRef              => Variant line to add annotations to {REF}
##         : $separator            => Separator for field
##         : $valueRef             => Field value {REF}
##         : $fieldID              => Field key ID

    my $transcriptCounterRef = $_[0];
    my $lineRef = $_[1];
    my $separator = $_[2];
    my $valueRef = $_[3];
    my $fieldID = $_[4];
    
    if ($$transcriptCounterRef == 0) {
	
	$$lineRef .= $fieldID.$$valueRef; #First selected feature
    }
    else {

	$$lineRef .= $separator.$$valueRef; #Following features
    }
}

sub AddFieldToElement {

##AddFieldToElement
    
##Function : Adds adds a field to an element.
##Returns  : ""
##Arguments: $selectedTranscriptTrackerRef, $selectedLineRef, $lineRef, $separator, $valueRef
##         : $selectedTranscriptTrackerRef => The selected transcript tracker for belonging to select file {REF}
##         : $selectedLineRef              => Selected line to add annotations to {REF}
##         : $lineRef                      => Variant line to add annotations to {REF}
##         : $separator                    => Separator for field
##         : $valueRef                     => Field value {REF}

    my $selectedTranscriptTrackerRef = $_[0];
    my $selectedLineRef = $_[1];
    my $lineRef = $_[2];
    my $separator = $_[3];
    my $valueRef = $_[4];
    
    if ($$selectedTranscriptTrackerRef == 1) {
	
	$$selectedLineRef .= $separator.$$valueRef;
    }
    ## Always include all transcripts in orphan list
    $$lineRef .= $separator.$$valueRef;
}

sub CollectConsequenceGenes {
 
##CollectConsequenceGenes
    
##Function : Collects all consequence and predictors per gene and adds info to line to be written.
##Returns  : ""
##Arguments: $consequenceHashRef, $fieldArrayRef, $selectedVariantLineRef, $variantVariantLineRef
##         : $consequenceHashRef     => Consequence(s) for each gene {REF}
##         : $selectDataHashRef      => Select file data {REF}
##         : $fieldArrayRef          => Features to be processed as determined by CSQ {REF}
##         : $selectedVariantLineRef => Selected line to add annotations to {REF}
##         : $variantVariantLineRef  => Variant line to add annotations to {REF}
			    
    my $consequenceHashRef = $_[0];
    my $selectDataHashRef = $_[1];
    my $fieldArrayRef = $_[2];
    my $selectedVariantLineRef = $_[3];
    my $variantLineRef = $_[4];
    
    my %geneCounter;
    my %selectedGeneCounter;
    my @tempFields;
    my @selectedTempFields;
    
    for (my $fieldCounter=0;$fieldCounter<scalar(@{$fieldArrayRef});$fieldCounter++) { #Set transcript counter to "0"
	
	$selectedGeneCounter{$fieldCounter} = 0;
	$geneCounter{$fieldCounter} = 0;
    }
    for my $genes (keys %{$consequenceHashRef}) {
	
	for (my $fieldCounter=0;$fieldCounter<scalar(@{$fieldArrayRef});$fieldCounter++) {
	    
	    if (${$selectDataHashRef}{$genes}) { #Exists in selected Features
		
		&CollectConsequenceField(\$fieldCounter, \%selectedGeneCounter, \%{$consequenceHashRef}, \$genes, \@{$fieldArrayRef}, \@selectedTempFields);
	    }

	    ## Always include all transcripts in research list
	    &CollectConsequenceField(\$fieldCounter, \%geneCounter, \%{$consequenceHashRef}, \$genes, \@{$fieldArrayRef}, \@tempFields);
	}
    }
    
    &AddToLine(\@{$fieldArrayRef}, \@selectedTempFields, \$$selectedVariantLineRef);
    &AddToLine(\@{$fieldArrayRef}, \@tempFields, \$$variantLineRef);
}

sub CollectConsequenceField {
   
##CollectConsequenceField
    
##Function : Collects consequences for features in @featureFields to temporary array for adding to line once all information are collected.
##Returns  : ""
##Arguments: $fieldCounterRef, $geneCounterHashRef, $consequenceHashRef, $geneRef, $fieldArrayRef, $selectedArrayRef
##         : $fieldCounterRef    => Field number in feature {REF}
##         : $geneCounterHashRef => Counts the number of transcripts per gene {REF}
##         : $consequenceHashRef => Consequence(s) for each gene {REF}
##         : $geneRef            => The gene symbol {REF}
##         : $fieldArrayRef      => Features to be processed as determined by CSQ {REF}
##         : $selectedArrayRef   => Selected array {REF}

    my $fieldCounterRef = $_[0];
    my $geneCounterHashRef = $_[1];
    my $consequenceHashRef = $_[2];
    my $geneRef = $_[3];
    my $fieldArrayRef = $_[4];
    my $selectedArrayRef = $_[5];

    my $alleleNumber = scalar(keys %{${$consequenceHashRef}{$$geneRef}});  #Enable tracking of multiple alleles

    if ($$geneCounterHashRef{$$fieldCounterRef} == 0) { #First time
	
	my $allelCounter = 1;

	for my $allele (keys %{${$consequenceHashRef}{$$geneRef}}) {  #All alleles
	    
	    if (defined($$consequenceHashRef{$$geneRef}{$allele}{ $$fieldArrayRef[$$fieldCounterRef] }) && ($$consequenceHashRef{$$geneRef}{$allele}{ $$fieldArrayRef[$$fieldCounterRef] } ne "")) { #If feature exists - else do nothing
		
		if ($allelCounter == 1) {
		    
		    $$selectedArrayRef[$$fieldCounterRef] .= ";".$$fieldArrayRef[$$fieldCounterRef]."=".$$geneRef.":".$allele."|".$$consequenceHashRef{$$geneRef}{$allele}{$$fieldArrayRef[$$fieldCounterRef]};
		}
		else {
		    
		    $$selectedArrayRef[$$fieldCounterRef] .= ",".$$geneRef.":".$allele."|".$$consequenceHashRef{$$geneRef}{$allele}{$$fieldArrayRef[$$fieldCounterRef]};
		}
		$$geneCounterHashRef{$$fieldCounterRef}++;
		$allelCounter++;
	    }
	}
    }
    else { #Subsequent passes
	
	for my $allele (keys %{${$consequenceHashRef}{$$geneRef}}) {  #All alleles
	    
	    if (defined($$consequenceHashRef{$$geneRef}{$allele}{$$fieldArrayRef[$$fieldCounterRef]}) && ($$consequenceHashRef{$$geneRef}{$allele}{ $$fieldArrayRef[$$fieldCounterRef] } ne "")) { #If feature exists - else do nothing
				    
		$$selectedArrayRef[$$fieldCounterRef] .= ",".$$geneRef.":".$allele."|".$$consequenceHashRef{$$geneRef}{$allele}{ $$fieldArrayRef[$$fieldCounterRef] };
	    }
	}
    }    
}

sub AddToLine {

##AddToLine
    
##Function : Adds to present line.
##Returns  : ""
##Arguments: $fieldArrayRef, $tempArrayRef, $lineRef
##         : $fieldArrayRef => Features to be processed as determined by CSQ {REF}
##         : $tempArrayRef  => Annotations for feature {REF}
##         : $lineRef       => Line to add to {REF}

    my $fieldsArrayRef = $_[0];
    my $tempArrayRef = $_[1];
    my $lineRef = $_[2];
    
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
##Arguments: $fieldArrayRef
##         : $fieldArrayRef => Holds the chromosomal coordinates and allel data

    my $fieldArrayRef = $_[0];
   
    my $chromosome = $$fieldArrayRef[0];
    my $startPosition = $$fieldArrayRef[1];
    my $referenceAllele = $$fieldArrayRef[3];
    my $alternativeAllele = $$fieldArrayRef[4];

    my $finalStartPosition = $startPosition; #The most "uppstream" position per variant
    my $finalStopPosition = 0; #The most "downstream" position per variant
	
    ## Convert to upper case
    ($referenceAllele, $alternativeAllele) = (uc $referenceAllele, uc $alternativeAllele);
    
    if ($alternativeAllele eq ".") { #No Variant Call
	
	next;
    }
    my @alternativeAlleles = split(/,/, $$fieldArrayRef[4]);
    
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
	if ($finalStartPosition < $newstart) { #New start is uppstream of old
	    
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
##Arguments: $hashRef, $rangeFileKey, $headerRef, $positionRef, $rangeFilePathRef
##         : $hashRef          => Hash to store metaData in {REF}
##         : $rangeFileKey     => Range file key
##         : $headerRef        => Header from range file {REF}
##         : $positionRef      => Column position in supplied range file
##         : $rangeFilePathRef => Range file path

    my $hashRef = $_[0];
    my $rangeFileKey = $_[1];
    my $headerRef = $_[2];
    my $positionRef = $_[3];
    my $rangeFilePathRef = $_[4];

    if (defined($$hashRef{$rangeFileKey}{$$headerRef})) { #Add INFO from predefined entries

	$$hashRef{'Present'}{$$headerRef}{'INFO'} = $$hashRef{$rangeFileKey}{$$headerRef}{'INFO'};
	$$hashRef{'Present'}{$$headerRef}{'ColumnOrder'} = $$positionRef; #Column position in supplied range input file
    }
    else { #Add arbitrary INFO field using input header

	$$hashRef{'Present'}{$$headerRef}{'INFO'} = q?##INFO=<ID=?.$$headerRef.q?,Number=.,Type=String,Description="String taken from ?.$$rangeFilePathRef.q?">?;
	$$hashRef{'Present'}{$$headerRef}{'ColumnOrder'} = $$positionRef; #Column position in supplied -sf_ac
    }
}

sub RangeAnnotations {

##RangeAnnotations
    
##Function : Creates the interval tree(s) for range and select files. Adds corresponding INFO fields to metadata.
##Returns  : ""
##Arguments: $rangeColumnsArrayRef, $lineElementsArrayRef, $hashRef, $headersArrayRef, $rangeFileKey
##         : $rangeColumnsArrayRef => Range feature annotations columns {REF}
##         : $lineElementsArrayRef => Range file line {REF} 
##         : $hashRef              => Hash to store metaData in {REF}
##         : $headersArrayRef      => Headers from range file {REF}
##         : $rangeFileKey         => Range file key
##         : $paddingRef           => The nucleotide distance to pad the range with {REF}

    my $rangeCoulumnsArrayRef = $_[0];
    my $lineElementsArrayRef = $_[1];
    my $hashRef = $_[2];
    my $rangeFileKey = $_[3];
    my $headersArrayRef = $_[4];
    my $paddingRef = $_[5];
    
    my $features; #Features to collect (Format: ";" separated elements)
    
    for (my $extractColumnsCounter=0;$extractColumnsCounter<scalar(@{$rangeCoulumnsArrayRef});$extractColumnsCounter++) { #Defines what scalar to store

	if(defined($$lineElementsArrayRef[ $$rangeCoulumnsArrayRef[$extractColumnsCounter] ])) {

	    $$lineElementsArrayRef[ $$rangeCoulumnsArrayRef[$extractColumnsCounter] ] =~ tr/ /_/; #Remove whitespace since this is not allowed in vcf INFO field
	
	    if ($extractColumnsCounter == 0) {
		
		$features .= $$lineElementsArrayRef[ $$rangeCoulumnsArrayRef[$extractColumnsCounter] ];
	    }		    
	    else {
		
		$features .= ";".$$lineElementsArrayRef[ $$rangeCoulumnsArrayRef[$extractColumnsCounter] ];
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
    
##Function : Checks if an interval tree exists (per chr) and collects features from input file and adds annotations to line.
##Returns  : ""
##Arguments: $rangeFileKey, $lineElementsArrayRef, $hashRef
##         : $rangeFileKey         => Range file key
##         : $lineElementsArrayRef => Infile vcf line elements array {REF}
##         : $hashRef              => Range file hash {REF}
##         : $printLineRef         => Line to add annotations to {REF}
  
    my $rangeFileKey = $_[0];
    my $lineElementsArrayRef = $_[1];
    my $hashRef = $_[2];
    my $printLineRef = $_[3];
    
    my %noIDRegion;  #No HGNC_symbol or ensemblGeneID, but still clinically releveant e.g. mtD-loop

    if(defined($tree{$rangeFileKey}{ $$lineElementsArrayRef[0] }) ) { #Range annotations
	
	my $feature; #Features to be collected
	
	my ($start, $stop) = &ConvertToRange($lineElementsArrayRef); #Convert SVs to range coordinates from vcf coordinates
	
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

			for my $rangeAnnotation (keys % {$$hashRef{'Present'}}) { #All selected annotations
			
			    if ($$hashRef{'Present'}{$rangeAnnotation}{'ColumnOrder'} eq $annotationsCounter) { #Correct feature
			
				if ($rangeAnnotation eq "Clinical_db_gene_annotation") {  #Special case, which is global and not gene centric

				    ## Collect unique elements from array reference and return array reference with unique elements
				    my $uniqueRef = &UniqElements(\@{$collectedAnnotations{$annotationsCounter}});
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
##Arguments: $hashRef
##         : $hashRef => The hash to store the meta data {REF}
    
    my $hashRef = $_[0];
    
    my ($base, $script) = (`date +%Y%m%d`,`basename $0`);  #Catches current date and script name
    chomp($base,$script);  #Remove \n;
    push(@{${$hashRef}{'Software'}{$script}}, "##Software=<ID=".$script.",Version=".$vcfParserVersion.",Date=".$base);
}


sub FindAF {

##FindAF
    
##Function : Adds the least alternative allele(s) frequency to each line
##Returns  : ""
##Arguments: $arrayRef, $regexp
##         : $arrayRef => The INFO array {REF}
##         : $regexp   => The regexp to used to locate correct ID field

    my $arrayRef = $_[0];
    my $regexp = $_[1];

    my $frequencyPosition = 0;

    if ($_[2]) {

	$frequencyPosition = $_[2];
    }
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
##Arguments: $arrayRef, $regexp
##         : $arrayRef => The INFO array {REF}
##         : $regexp   => The regexp to used to locate correct ID field

    my $arrayRef = $_[0];
    my $regexp = $_[1];
    my $scoreCutOff = $_[2];

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
##Arguments: $arrayRef, $regexp
##         : $arrayRef => The INFO array {REF}
##         : $regexp   => The regexp to used to locate correct ID field

    my $arrayRef = $_[0];
    my $regexp = $_[1];

    my $frequencyPosition = 0;

    if ($_[2]) {

	$frequencyPosition = $_[2];
    }
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
    
    my $metaDataHashRef = $_[0];
    my $metaDataString = $_[1];
    
    if ($metaDataString=~/^##fileformat/) {  #Catch fileformat as it has to be at the top of header
	
	push(@{${$metaDataHashRef}{'fileformat'}{'fileformat'}}, $metaDataString);  #Save metadata string
    }
    elsif ($metaDataString=~/^##contig/) {  #catch contigs to not sort them later

	push(@{${$metaDataHashRef}{'contig'}{'contig'}}, $metaDataString);  #Save metadata string
    }
    elsif ($metaDataString=~/^##(\w+)=(\S+)/) {  #FILTER, FORMAT, INFO etc and more custom records
	
	push(@{${$metaDataHashRef}{$1}{$2}}, $metaDataString);  #Save metadata string
    }
    else {  #All oddities
	
	push(@{${$metaDataHashRef}{'other'}{'other'}}, $metaDataString);  #Save metadata string
    }
}


sub WriteMetaData {

##WriteMetaData
    
##Function : Writes metadata to filehandle speciied by order in metaDataOrders.
##Returns  : ""
##Arguments: $hashRef, $FILEHANDLE, $SELECTFILEHANDLE
##         : $hashRef          => Hash for metaData {REF}
##         : $FILEHANDLE       => The filehandle to write to
##         : $SELECTFILEHANDLE => The filehandle to write to

    my $hashRef = $_[0];
    my $FILEHANDLE = $_[1];
    my $SELECTFILEHANDLE = $_[2];

    my @metaDataOrders = ("fileformat", "FILTER", "FORMAT", "INFO", "FIX_INFO", "contig", "Software");  #Determine order to print for standard records
    my @lines;

    for (my $lineCounter=0;$lineCounter<scalar(@metaDataOrders);$lineCounter++) {
	
	if (${$hashRef}{ $metaDataOrders[$lineCounter] }) {  #MetaDataRecordExists
	    
	    if ($metaDataOrders[$lineCounter] eq "contig") {  #Should not be sorted, but printed "as is"
		
		foreach my $line (@{${$hashRef}{$metaDataOrders[$lineCounter]}{ $metaDataOrders[$lineCounter] }}) {
		    
		    print $FILEHANDLE $line, "\n";
		    
		    if (defined($_[2])) {
			
			print $SELECTFILEHANDLE $line, "\n";
		    }
		}
		delete(${$hashRef}{ $metaDataOrders[$lineCounter] });  #Enable print of rest later
	    }	    
	    else {
		
		foreach my $line (sort( keys %{${$hashRef}{ $metaDataOrders[$lineCounter] }})) {
		    
		    print $FILEHANDLE @{${$hashRef}{ $metaDataOrders[$lineCounter] }{$line}}, "\n";
		    
		    if (defined($_[2])) {
			
			print $SELECTFILEHANDLE @{${$hashRef}{ $metaDataOrders[$lineCounter] }{$line}}, "\n";
		    }
		}
		if (defined($_[2])) {
		    
		    foreach my $line (sort( keys %{${$hashRef}{'Select'}{ $metaDataOrders[$lineCounter] }})) {
			
			print $SELECTFILEHANDLE @{${$hashRef}{'Select'}{ $metaDataOrders[$lineCounter] }{$line}}, "\n";
		    }
		}
		foreach my $line (sort( keys %{${$hashRef}{'Range'}{ $metaDataOrders[$lineCounter] }})) {
		    
		    print $FILEHANDLE @{${$hashRef}{'Range'}{ $metaDataOrders[$lineCounter] }{$line}}, "\n";
		}
		delete(${$hashRef}{ $metaDataOrders[$lineCounter] });  #Enable print of rest later
		if (${$hashRef}{'Select'}{ $metaDataOrders[$lineCounter] }) {
		    
		    delete(${$hashRef}{'Select'}{ $metaDataOrders[$lineCounter] });
		}
		if (${$hashRef}{'Range'}{ $metaDataOrders[$lineCounter] }) {
		    
		    delete(${$hashRef}{'Range'}{ $metaDataOrders[$lineCounter] });
		}
	    }
	}
    }
    for my $keys (keys %{$hashRef}) {
	
	for my $line ( sort(keys %{${$hashRef}{$keys}}) ) {
	    
	    print $FILEHANDLE @{${$hashRef}{$keys}{$line}}, "\n";
	    if (defined($_[2])) {
		
		print $SELECTFILEHANDLE @{${$hashRef}{$keys}{$line}}, "\n";
	    }
	}
    }
}


sub UniqElements {
    
##UniqElements
    
##Function : Collect unique elements from array reference and return array reference with unique elements
##Returns  : "array reference"
##Arguments: $arrayRef, $familyIDRef, $selectFile
##         : $arrayRef => The array whose elements are to be made distinct {REF}
    
    my $arrayRef = $_[0];
    my %seen;
    return [ grep { !$seen{$_}++ } @{$arrayRef} ];  #For each element in array, see if seen before and only return list distinct elements  
}

sub AddRangeMetaDataToVcf {
		
##AddRangeMetaDataToVcf

##Function : Adds range file meta data annotation headers to meta data hash.
##Returns  : ""
##Arguments: $metaDataHashRef, $vcfHeaderHashRef, $dataHashRef, $featureAnnotationCoulumnsArrayRef, $fileKey
##         : $metaDataHashRef                   => Vcf meta data {REF}
##         : $vcfHeaderHashRef                  => Vcf header meta data
##         : $dataHashRef                       => Range file hash {REF}
##         : $featureAnnotationCoulumnsArrayRef => Range columns to include {REF}
##         : $fileKey                           => Range file key used to seperate range file(s) i.e., select and range
    
    my ($argHashRef) = @_;
    
    ## Flatten argument(s)
    my $metaDataHashRef = ${$argHashRef}{'metaDataHashRef'};
    my $vcfHeaderHashRef = ${$argHashRef}{'vcfHeaderHashRef'};
    my $featureAnnotationColumnsArrayRef = ${$argHashRef}{'featureAnnotationColumnsArrayRef'};
    my $dataHashRef = ${$argHashRef}{'dataHashRef'};
    
    if (scalar(@{$featureAnnotationColumnsArrayRef}) > 0) { #RangeFile annotations
	
	for my $annotation (keys %{${$dataHashRef}{'Present'}}) {
	    
	    unless (defined(${$vcfHeaderHashRef}{'INFO'}{$annotation})) { #Unless INFO header is already present add to file 
		
		push(@{${$metaDataHashRef}{ ${$argHashRef}{'fileKey'} }{'INFO'}{$annotation}}, ${$dataHashRef}{'Present'}{$annotation}{'INFO'});  #Save specific rangeFile INFO
	    }
	}
    }
}


sub SetDefaultArg {

##SetDefaultArg
    
##Function : Set the default arguments for argHashRef using $defaultHashRef
##Returns  : ""
##Arguments: $argHashRef, $parameterValueRef, $parameterName, $associatedProgram
##         : $argHashRef     => The argument hash {REF}
##         : $defaultHashRef => The default hash {REF}

    my $argHashRef = $_[0];
    my $defaultHashRef = $_[1];
    
    foreach my $key (keys %{$defaultHashRef}) {
	
	unless (defined(${$argHashRef}{$key})) {
	    
	    ${$argHashRef}{$key} = ${$defaultHashRef}{$key};  #Set default
	}
    }
}


sub CheckTerms {

##CheckTerms
    
##Function : Check that the found terms in the vcf corresond to known terms - otherwise croak and exit.
##Returns  : ""
##Arguments: $termHashRef, $termRef, 
##         : $termHashRef => The term hash {REF}
##         : $termRef => The found term {REF}
##         : $termName => The origin of the term i.e Sift, PolyPhen, SO

    my $termHashRef = $_[0];
    my $termRef = $_[1];
    my $termName = $_[2];
    
    if ( (defined($$termRef)) && ($$termRef ne "") ) {
	
	unless (exists(${$termHashRef}{ $$termRef })) {
	    
	    warn("Could not find ".$termName." term from vcf in corresponding hash. Update hash to contain term: '".$$termRef."'\n");
	    exit 1;
	}
    }   
}

