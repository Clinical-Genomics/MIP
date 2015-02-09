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
           -of/--outputFormat format (Default: vcf)
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

my ($infile, $outputFormat, $parseVEP, $rangeFeatureFile, $selectFeatureFile, $selectFeatureMatchingColumn, $selectOutfile, $writeSoftwareTag, $padding) = ("", "vcf", 0, 0, 0, "nocmdinput", "nocmdinput", 1, 5000);
my (@metaData, @selectMetaData, @rangeFeatureAnnotationColumns, @selectFeatureAnnotationColumns); 
my (%geneAnnotation, %consequenceSeverity, %rangeData, %selectData, %snpEffCmd, %tree, %metaData);

my ($help, $version) = (0, 0);

if (defined($ARGV) && $ARGV[0]!~/^-/) { #Collect potential infile - otherwise read from STDIN
    
    $infile = $ARGV[0];
}

###User Options
GetOptions('of|outputFormat:s' => \$outputFormat,
	   'pVEP|parseVEP:s' => \$parseVEP,
	   'rf|rangeFeatures:s' => \$rangeFeatureFile,
	   'rf_ac|rangeFeatureAnnotationColumns:s'  => \@rangeFeatureAnnotationColumns, #Comma separated list
	   'sf|selectFeatures:s' => \$selectFeatureFile,
	   'sf_mc|selectFeatureMatchingColumn:n' => \$selectFeatureMatchingColumn,
	   'sf_ac|selectFeatureAnnotationColumns:s'  => \@selectFeatureAnnotationColumns, #Comma separated list
	   'sof|selectOutfile:s' => \$selectOutfile,
	   'wst|writeSoftwareTag:n' => \$writeSoftwareTag,
	   'pad|padding:n' => \$padding,
	   'h|help' => \$help,  #Display help text
	   'v|version' => \$version, #Display version number
    );

if($help) {
    
    print STDOUT "\n".$USAGE, "\n";
    exit;
}

my $vcfParserVersion = "1.1.1";

if($version) {
    
    print STDOUT "\nvcfParser.pl v".$vcfParserVersion, "\n\n";
    exit
}
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

    &ReadRangeFile($rangeFeatureFile);
}

if ($selectFeatureFile ne 0) {

    &ReadSelectFile($selectFeatureFile, $selectFeatureMatchingColumn);
}

&DefineSnpEffAnnotations();
&DefineConsequenceSeverity();

&ReadInfileVCF(\%metaData, $infile, $selectOutfile);

###
#Sub Routines
###

sub DefineSelectData {

##Defines arbitrary INFO fields based on headers in selectFile

    $selectData{'SelectFile'}{'HGNC_symbol'}{'INFO'} = q?##INFO=<ID=HGNC_symbol,Number=.,Type=String,Description="The HGNC gene symbol">?;
    $selectData{'SelectFile'}{'Ensembl_gene_id'}{'INFO'} = q?##INFO=<ID=Ensembl_gene_id,Number=.,Type=String,Description="Ensembl gene identifier">?;
    $selectData{'SelectFile'}{'OMIM_morbid'}{'INFO'} = q?##INFO=<ID=OMIM_morbid,Number=.,Type=String,Description="OMIM morbid ID associated with gene(s)">?;
    $selectData{'SelectFile'}{'Phenotypic_disease_model'}{'INFO'} = q?##INFO=<ID=Phenotypic_disease_model,Number=.,Type=String,Description="Known disease gene(s) phenotype inheritance model">?;
    $selectData{'SelectFile'}{'Clinical_db_gene_annotation'}{'INFO'} = q?##INFO=<ID=Clinical_db_gene_annotation,Number=.,Type=String,Description="Gene disease group association">?;
    $selectData{'SelectFile'}{'Reduced_penetrance'}{'INFO'} = q?##INFO=<ID=Reduced_penetrance,Number=.,Type=String,Description="Pathogenic gene which can exhibit reduced penetrance">?;
    $selectData{'SelectFile'}{'Disease_associated_transcript'}{'INFO'} = q?##INFO=<ID=Disease_associated_transcript,Number=.,Type=String,Description="Known pathogenic transcript(s) for gene">?;

}

sub DefineSnpEffAnnotations {

##Defines the snpEff annotations that can be parsed and modified
    
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

    $snpEffCmd{'SnpEff'}{'CLNSIG'}{'File'} = q?clinvar_\d+.vcf?;
    $snpEffCmd{'SnpEff'}{'CLNSIG'}{'INFO'} = q?##INFO=<ID=CLNSIG,Number=A,Type=String,Description="Variant Clinical Significance, 0 - Uncertain significance, 1 - not provided, 2 - Benign, 3 - Likely benign, 4 - Likely pathogenic, 5 - Pathogenic, 6 - drug response, 7 - histocompatibility, 255 - other">?;
    $snpEffCmd{'SnpEff'}{'CLNSIG'}{'FIX_INFO'} = q?##INFO=<ID=SnpSift_CLNSIG,Number=A,Type=String,Description="Variant Clinical Significance, 0 - Uncertain significance, 1 - not provided, 2 - Benign, 3 - Likely benign, 4 - Likely pathogenic, 5 - Pathogenic, 6 - drug response, 7 - histocompatibility, 255 - other">?;

    $snpEffCmd{'SnpEff'}{'phastCons100way_vertebrate_prediction_term'}{'File'} = q?SnpSift dbnsfp?;
    $snpEffCmd{'SnpEff'}{'phastCons100way_vertebrate_prediction_term'}{'INFO'} = q?##INFO=<ID=phastCons100way_vertebrate_prediction_term,Number=A,Type=String,Description="PhastCons conservation prediction term">?;

    $snpEffCmd{'SnpEff'}{'phyloP100way_vertebrate_prediction_term'}{'File'} = q?SnpSift dbnsfp?;
    $snpEffCmd{'SnpEff'}{'phyloP100way_vertebrate_prediction_term'}{'INFO'} = q?##INFO=<ID=phyloP100way_vertebrate_prediction_term,Number=A,Type=String,Description="PhyloP conservation prediction term">?;

    $snpEffCmd{'SnpEff'}{'GERP++_RS_prediction_term'}{'File'} = q?SnpSift dbnsfp?;
    $snpEffCmd{'SnpEff'}{'GERP++_RS_prediction_term'}{'INFO'} = q?##INFO=<ID=GERP++_RS_prediction_term,Number=A,Type=String,Description="GERP RS conservation prediction term">?;

}

sub DefineConsequenceSeverity {
##Defines the precedence of consequences

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
    $consequenceSeverity{'initiator_codon_variant'}{'Rank'} = 6;
    $consequenceSeverity{'initiator_codon_variant'}{'GeneticRegionAnnotation'} = "exonic";
    $consequenceSeverity{'inframe_insertion'}{'Rank'} = 6;
    $consequenceSeverity{'inframe_insertion'}{'GeneticRegionAnnotation'} = "exonic";
    $consequenceSeverity{'inframe_deletion'}{'Rank'} = 6;
    $consequenceSeverity{'inframe_deletion'}{'GeneticRegionAnnotation'} = "exonic";
    $consequenceSeverity{'missense_variant'}{'Rank'} = 6;
    $consequenceSeverity{'missense_variant'}{'GeneticRegionAnnotation'} = "exonic";
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

sub ReadRangeFile {
##Reads a file containg features to be annotated using range queries e.g. EnsemblGeneID

    my $infileName = $_[0];

    my @headers; #Save headers from rangeFile

    open(RRF, "<".$infileName) or die "Can't open ".$infileName.":".$!, "\n"; 

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
	    next;
	}
	if ( $_ =~/^(\S+)/ ) {	
	    
	    my @lineElements = split("\t",$_); #Loads range file line elements

##Create Interval Tree
	    if (scalar(@rangeFeatureAnnotationColumns) > 0) {#Annotate vcf with features from select file
		
		&RangeAnnotations(\@rangeFeatureAnnotationColumns, \@lineElements, \%rangeData, "RangeFile", \$infileName, \@headers, \$padding);
	    }
	}
    }
    close(RRF);
    print STDERR "Finished Reading Range file: ".$infileName,"\n";
}

sub ReadSelectFile {
##Reads a file containg features to be analysed seperately

    my $infileName = $_[0];
    my $selectFeatureColumn = $_[1];

    my @headers; #Save headers from selectFile

    open(RSF, "<".$infileName) or die "Can't open ".$infileName.":".$!, "\n"; 

    while (<RSF>) {
	
	chomp $_; #Remove newline

	if (m/^\s+$/) { # Avoid blank lines
	    
	    next;
	}
	if ($_=~/^##/) { # MetaData - Avoid
	    
	    next;
	}
	if ($_=~/^#/) { # Header/Comment
	    
	    @headers = split(/\t/, $_);
	    next;
	}
	if ( $_ =~/^(\S+)/ ) {	
	    
	    my @lineElements = split("\t",$_); # Loads select line elements
	    $lineElements[$selectFeatureColumn] =~ s/\s/_/g; # Replace whitespace with "_"
	    $selectData{$lineElements[$selectFeatureColumn]} = $lineElements[$selectFeatureColumn];

	    ## Create Interval Tree
	    if (scalar(@selectFeatureAnnotationColumns) > 0) {#Annotate vcf with features from select file
		
		&RangeAnnotations(\@selectFeatureAnnotationColumns, \@lineElements, \%selectData, "SelectFile", \$selectFeatureFile, \@headers, \$padding);
	    }
	}
    }
    close(RSF);
    print STDERR "Finished Reading Select file: ".$infileName,"\n";
}

sub ReadInfileVCF {
#Reads infile vcf format 

    my $metaDataHashRef = $_[0];
    my $infileName = $_[1];
    my $selectOutFileName = $_[2];

    my @vepFormatField;
    my %vepFormatFieldColumn;

    my @printTSVFields = ("FeatureType");
    my @featureFields;

    my %vcfHeader;

    if ($selectFeatureFile ne 0) {

	open(WOSFTSV, ">".$selectOutFileName) or die "Can't open ".$selectOutFileName.":".$!, "\n";
    }
    
    while (<>) {
	
	chomp $_;  # Remove newline
	
	if (m/^\s+$/) {	# Avoid blank lines
	    next;
	}
	if ($_=~/^##(\w+)=/) {  # MetaData

	    &ParseMetaData(\%{$metaDataHashRef}, $_);

	    if ($_=~/INFO\=\<ID\=(\w+)/) { # Collect all INFO keys
	    
		$vcfHeader{'INFO'}{$1} = $1; #Save to hash
	    }
	    if ($_=~/SnpSiftCmd\=/) { #Find SnpEff command meta line

		for my $database (keys % {$snpEffCmd{'SnpEff'}}) {

		    if ($_=~/$snpEffCmd{'SnpEff'}{$database}{'File'}/) { #SnpEff/Sift has been used to annotate input vcf
			
			unless (defined($vcfHeader{'INFO'}{$database})) { #Unless INFO header is already present add to metaDataHeader
			    
			    $snpEffCmd{'Present'}{'database'}{$database} = $database; #Save which frequency db has been used for later
			    push(@{${$metaDataHashRef}{'INFO'}{$database}}, $snpEffCmd{'SnpEff'}{$database}{'INFO'});

			    if (defined($snpEffCmd{'SnpEff'}{$database}{'FIX_INFO'})) { #If 'FIX_INFO' flag is present add to metaDataHeader

				push(@{${$metaDataHashRef}{'FIX_INFO'}{$database}}, $snpEffCmd{'SnpEff'}{$database}{'FIX_INFO'});
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

		    if ($outputFormat eq "vcf") {
			
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
		}
		next;
	    }
	    next;
	}
	if ($_=~/^#CHROM/) {

	    if (scalar(@rangeFeatureAnnotationColumns) > 0) { #RangeFile annotations
		
		for my $rangeAnnotation (keys % {$rangeData{'Present'}}) {
		    
		    unless (defined($vcfHeader{'INFO'}{$rangeAnnotation})) { #Unless INFO header is already present add to file 
			
			push(@{${$metaDataHashRef}{'Range'}{'INFO'}{$rangeAnnotation}}, $rangeData{'Present'}{$rangeAnnotation}{'INFO'});  #Save specific rangeFile INFO
		    }
		}
	    }
	    ## Select MetaData
	    if (scalar(@selectFeatureAnnotationColumns) > 0) { #SelectFile annotations
		
		for my $selectAnnotation (keys % {$selectData{'Present'}}) {
		    
		    unless (defined($vcfHeader{'INFO'}{$selectAnnotation})) { #Unless INFO header is already present add to file 

			push(@{${$metaDataHashRef}{'Select'}{'INFO'}{$selectAnnotation}}, $selectData{'Present'}{$selectAnnotation}{'INFO'});  #Save specific selectFile INFO
			
		    }
		}
	    }
	    if ($writeSoftwareTag == 1) {

		&AddProgramToMeta(\%{$metaDataHashRef});
	    }
	    if (scalar(@selectFeatureAnnotationColumns) > 0) { #SelectFile annotations

		&WriteMetaData(\%{$metaDataHashRef}, *STDOUT, *WOSFTSV);
		print STDOUT $_, "\n";  #Write header line
		print WOSFTSV $_, "\n";  #Write header line
	    }
	    else {

		&WriteMetaData(\%{$metaDataHashRef}, *STDOUT);
		print STDOUT $_, "\n";  #Write header line
	    }
	    
	    if ($parseVEP == 1) {

		if ($outputFormat eq "tsv") {
		    
		    $metaData[-1] .= "\tHGNC_transcript_info"."\tFunctional_annotation"."\tGene_annotation";
		}
		@featureFields = ("MostSevereConsequence", "GeneticRegionAnnotation");
		
		if ($vepFormatFieldColumn{'SIFT'}) {
		    
		    if ($outputFormat eq "tsv") {
			
			$metaData[-1] .="\tSIFT";
		    }
		    push(@featureFields, "Sift");
		}
		if ($vepFormatFieldColumn{'PolyPhen'}) {
		    
		    if ($outputFormat eq "tsv") {
			
			$metaData[-1] .="\tPoly_phen_hdiv";
		    }
		    push(@featureFields, "PolyPhen");
		}
	    }
	    next;
	}
	if ( $_ =~/^(\S+)/ ) {	
	    
	    my %variantData;
	    my %selectedVariantData;
	    my %consequence;
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
	    for my $database (keys % {$snpEffCmd{'Present'}{'database'}}) { #Note that the vcf should only contain 1 database entry

		if ( ($database eq "Dbsnp129LCAF") || ($database eq "DbsnpLCAF") ) {

		    if ($lineElements[7] =~/\S+_CAF=\[(.+)\]/) {
			
			my @tempMafs = sort {$a <=> $b} grep { $_ ne "." } split(",", $1); #Split on ",", remove entries containing only "." and sort remaining entries numerically

			if (scalar(@tempMafs) > 0) {

			    ## Save LCAF frequency info
			    $variantLine .= $database."=".$tempMafs[0].";";
			    $selectedVariantLine .= $database."=".$tempMafs[0].";";
			}
		    }
		}
		elsif($database eq "1000GAF") {
			
		    my @tempArray = split(/;/, $lineElements[7]);  #Split INFO field to key=value items

		    my $tempMaf = &FindAF(\@tempArray, "\\S+_AF=");

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
		    
		    my $tempMaf = &FindAF(\@tempArray, "\\S+_AF=");
		    
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
	    &TreeAnnotations("SelectFile", \@lineElements, \%selectData, \$selectedVariantLine);
	    &TreeAnnotations("RangeFile", \@lineElements, \%rangeData, \$variantLine);
	    
	    my @variantEffects = split(/;/, $lineElements[7]); #Split INFO field
	    my $CSQTranscripts;
	    
	    for (my $variantEffectCounter=0;$variantEffectCounter<scalar(@variantEffects);$variantEffectCounter++) {
		
		if ($parseVEP == 1) {
		    
		    ## Find CSQ field and extract transripts that belong to select genes
		    if ($variantEffects[$variantEffectCounter]=~/CSQ\=(\S+)/) { #Find CSQ
			
			$CSQTranscripts = $1;
			my @transcripts = split(/,/, $1); #Split into transcripts elements
			
			if ($outputFormat eq "vcf") { #Collect CSQ info and write to correct file
			    
			    for (my $fieldCounter=0;$fieldCounter<scalar(@transcripts);$fieldCounter++) { #CSQ field
				
				my @transcriptsEffects = split(/\|/, $transcripts[$fieldCounter]); #Split in "|"
				
				if ($transcriptsEffects[ $vepFormatFieldColumn{'Feature'} ] =~/N\w_\d+\.\d/) { #RefSeq
				    
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
					## Always include all Refseq transcripts in research list
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
		if ($outputFormat eq "tsv") {
		    
		    $selectedVariantLine .= "\t".$sampleIDInfo;
		    $variantLine .= "\t".$sampleIDInfo;	
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
			
			if ($transcriptsEffects[ $vepFormatFieldColumn{'Feature'} ] =~/N\w_\d+\.\d/) { #RefSeq
			    
			    my $selectedTranscriptTracker = 0; #Track if any transcripts belong to selected features
			    
			    if (defined($transcriptsEffects[ $vepFormatFieldColumn{'SYMBOL'} ])) { #Save HGNC Symbol
				
				$variantData{'Symbol'} = $transcriptsEffects[ $vepFormatFieldColumn{'SYMBOL'} ];
				
				if ($selectData{ $variantData{'Symbol'} }) { #Exists in selected Features
				    
				    $selectedTranscriptTracker = 1; #Record belongs to selected Features
				    my $alleleGeneEntry = $transcriptsEffects[ $vepFormatFieldColumn{'Allele'} ].":".$variantData{'Symbol'};
				   
				    &AddFieldToElementCounter(\$selectedTranscriptCounter, \$selectedVariantLine, ",", \$alleleGeneEntry, \$outputFormat, "HGVScp=");
				}
				## Always include all Refseq transcripts in research list
				my $alleleGeneEntry = $transcriptsEffects[ $vepFormatFieldColumn{'Allele'} ].":".$variantData{'Symbol'};
				    
				&AddFieldToElementCounter(\$transcriptsCounter, \$variantLine, ",", \$alleleGeneEntry, \$outputFormat, "HGVScp=");
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
				## Always include all Refseq transcripts in research list				    
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
				    
				    if ( defined($consequence{ $variantData{'Symbol'} }{$allele}{'Score'}) ) { #Compare to previous record
					
					if ($consequenceSeverity{$consequences[$consequencesCounter]}{'Rank'} < $consequence{ $variantData{'Symbol'} }{$allele}{'Score'}) { #Collect most severe consequence
					    
					    &AddToConsequenceHash(\$consequence{ $variantData{'Symbol'} }{$allele}{'Score'}, \$consequenceSeverity{$consequences[$consequencesCounter]}{'Rank'});
					    &AddToConsequenceHash(\$consequence{ $variantData{'Symbol'} }{$allele}{'GeneticRegionAnnotation'}, \$consequenceSeverity{$consequences[$consequencesCounter]}{'GeneticRegionAnnotation'});
					    &AddToConsequenceHash(\$consequence{ $variantData{'Symbol'} }{$allele}{'MostSevereConsequence'}, \$consequences[$consequencesCounter]);
					    
					    if (defined($vepFormatFieldColumn{'SIFT'}) ) {    
						
						&AddToConsequenceHash(\$consequence{ $variantData{'Symbol'} }{$allele}{'Sift'}, \$transcriptsEffects[ $vepFormatFieldColumn{'SIFT'} ]);
					    }
					    if (defined($vepFormatFieldColumn{'PolyPhen'}) ) {
						
						&AddToConsequenceHash(\$consequence{ $variantData{'Symbol'} }{$allele}{'PolyPhen'}, \$transcriptsEffects[ $vepFormatFieldColumn{'PolyPhen'} ]);
					    }
					}
				    }
				    else { #First pass
					
					&AddToConsequenceHash(\$consequence{ $variantData{'Symbol'} }{$allele}{'Score'}, \$consequenceSeverity{$consequences[$consequencesCounter]}{'Rank'});
					&AddToConsequenceHash(\$consequence{ $variantData{'Symbol'} }{$allele}{'GeneticRegionAnnotation'}, \$consequenceSeverity{$consequences[$consequencesCounter]}{'GeneticRegionAnnotation'});
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
			    if (defined($transcriptsEffects[ $vepFormatFieldColumn{'STRAND'} ])) { #Save strand 
				
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
			    if (defined($transcriptsEffects[ $vepFormatFieldColumn{'HGVSc'} ])) { #Save HGVS cDNA change
				
				my @cDNAChange = split(/:/, $transcriptsEffects[ $vepFormatFieldColumn{'HGVSc'} ]);
				&AddFieldToElement(\$selectedTranscriptTracker, \$selectedVariantLine, \$variantLine, ":", \$cDNAChange[1]);
			    }
			    if (defined($transcriptsEffects[ $vepFormatFieldColumn{'HGVSp'} ])) {
				
				my @pAAChange = split(/:/, $transcriptsEffects[ $vepFormatFieldColumn{'HGVSp'} ]);
				&AddFieldToElement(\$selectedTranscriptTracker, \$selectedVariantLine, \$variantLine, ":", \$pAAChange[1]);
			    }
			    if ($selectedTranscriptTracker == 1) {
				
				$selectedTranscriptCounter++;
			    }
			    ## Always include all Refseq transcripts in research list	
			    $transcriptsCounter++;
			}
		    }
		    &CollectConsequenceGenes(\%consequence, \@featureFields, \$selectedVariantLine, \$variantLine);
		}
	    }
	    if ($outputFormat eq "vcf") {

		$selectedVariantLine .= "\t".$sampleIDInfo;
		$variantLine .= "\t".$sampleIDInfo;	
	    }
	    if ($parseVEP == 1) {
		
		if ($selectedTranscriptCounter > 0) { #Write to selected file
		    
		    print WOSFTSV $selectedVariantLine, "\n";
		}
		if ($transcriptsCounter > 0) { #Write to transcript file
		    
		    print STDOUT $variantLine, "\n";
		}
		elsif ( ($selectedTranscriptCounter == 0) && ($transcriptsCounter == 0) ) {

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

sub PrintToFileHandles {
##Prints to open FILEHANDLES

    my $arrayRef = $_[0];
    my $statement = $_[1];
   
    for (my $filehandleCounter=0;$filehandleCounter<scalar(@{$arrayRef});$filehandleCounter++) {
	 
	print {${$arrayRef}[$filehandleCounter]} $statement;
    }
}

sub RemoveFileEnding {
##Removes ".fileEnding" in filename.FILENDING

    my $fileNameRef = $_[0];
    my $fileEnding = $_[1];

    my $fileNameNoEnding;
    
    if ( (defined($$fileNameRef)) && $$fileNameRef =~/(\S+)($fileEnding$|$fileEnding.gz$)/) {

	$fileNameNoEnding = $1;
    }
    return $fileNameNoEnding;
}

sub AddToConsequenceHash {
##Adds the most severe consequence or prediction to gene
    
    my $hashKeyRef = $_[0];
    my $valueRef = $_[1];
    
    if (defined($valueRef)) {
	
	$$hashKeyRef = $$valueRef;
    }
}

sub AddFieldToElementCounter {
##Adds a field to an element. Adjust addition depending on if field has been seen before
    
    my $transcriptCounterRef = $_[0];
    my $lineRef = $_[1];
    my $separator = $_[2];
    my $valueRef = $_[3];
    my $formatRef = $_[4];
    my $fieldID = $_[5];
    
    if ($$transcriptCounterRef == 0) {
	
	if ($$formatRef eq "vcf") {

	    $$lineRef .= $fieldID.$$valueRef; #First selected feature
	}
	if ($$formatRef eq "tsv") {
	    
	    $$lineRef .= "\t".$$valueRef; #First selected feature
	}
    }
    else {

	    $$lineRef .= $separator.$$valueRef; #Following features
	}
}

sub AddFieldToElement {
##Adds adds a field to an element
    
    my $selectedTranscriptTrackerRef = $_[0];
    my $selectedLineRef = $_[1];
    my $lineRef = $_[2];
    my $separator = $_[3];
    my $valueRef = $_[4];
    
    if ($$selectedTranscriptTrackerRef == 1) {
	
	$$selectedLineRef .= $separator.$$valueRef;
    }
    ## Always include all Refseq transcripts in research list
    $$lineRef .= $separator.$$valueRef;
}

sub CollectConsequenceGenes {
##Collects all consequence and predictors per gene and adds info to line to be written 
			    
    my $hashRef = $_[0];
    my $arrayRef = $_[1];
    my $selectedVariantLineRef = $_[2];
    my $variantLineRef = $_[3];
    
    my %geneCounter;
    my %selectedGeneCounter;
    my @tempFields;
    my @selectedTempFields;
    
    for (my $fieldCounter=0;$fieldCounter<scalar(@{$arrayRef});$fieldCounter++) { #Set transcript counter to "0"
	
	$selectedGeneCounter{$fieldCounter} = 0;
	$geneCounter{$fieldCounter} = 0;
    }
    for my $genes (keys %{$hashRef}) {
	
	for (my $fieldCounter=0;$fieldCounter<scalar(@{$arrayRef});$fieldCounter++) {
	    
	    if ($selectData{$genes}) { #Exists in selected Features
		
		&CollectConsequenceField(\$fieldCounter, \%selectedGeneCounter, \%{$hashRef}, \$genes, \@{$arrayRef}, \@selectedTempFields, \$outputFormat);
	    }

	    ## Always include all Refseq transcripts in research list
	    &CollectConsequenceField(\$fieldCounter, \%geneCounter, \%{$hashRef}, \$genes, \@{$arrayRef}, \@tempFields, \$outputFormat);
	}
    }
    
    &AddToLine(\@{$arrayRef}, \@selectedTempFields, \$$selectedVariantLineRef, \$outputFormat);
    &AddToLine(\@{$arrayRef}, \@tempFields, \$$variantLineRef, \$outputFormat);
}

sub CollectConsequenceField {
##Collects consequences for features in @featureFields to temporary array for adding to line once all information are collected    

    my $fieldCounterRef = $_[0];
    my $hashCounterRef = $_[1];
    my $hashRef = $_[2];
    my $geneRef = $_[3];
    my $fieldArrayRef = $_[4];
    my $selectedArrayRef = $_[5];
    my $formatRef = $_[6];

    my $alleleNumber = scalar(keys %{${$hashRef}{$$geneRef}});  #Enable tracking of multiple alleles

    if ($$hashCounterRef{$$fieldCounterRef} == 0) { #First time
	
	my $allelCounter = 1;
	for my $allele (keys %{${$hashRef}{$$geneRef}}) {  #All alleles
	    
	    if (defined($$hashRef{$$geneRef}{$allele}{$$fieldArrayRef[$$fieldCounterRef]}) && ($$hashRef{$$geneRef}{$allele}{$$fieldArrayRef[$$fieldCounterRef]} ne "")) { #If feature exists - else do nothing
		
		if ($allelCounter == 1) {
		    
		    $$selectedArrayRef[$$fieldCounterRef] .= ";".$$fieldArrayRef[$$fieldCounterRef]."=".$$geneRef.":".$allele."|".$$hashRef{$$geneRef}{$allele}{$$fieldArrayRef[$$fieldCounterRef]};
		}
		else {
		    
		    $$selectedArrayRef[$$fieldCounterRef] .= ",".$$geneRef.":".$allele."|".$$hashRef{$$geneRef}{$allele}{$$fieldArrayRef[$$fieldCounterRef]};
		}
		$$hashCounterRef{$$fieldCounterRef}++;
		$allelCounter++;
	    }
	}
    }
    else { #Subsequent passes
	
	for my $allele (keys %{${$hashRef}{$$geneRef}}) {  #All alleles
	    
	    if (defined($$hashRef{$$geneRef}{$allele}{$$fieldArrayRef[$$fieldCounterRef]}) && ($$hashRef{$$geneRef}{$allele}{$$fieldArrayRef[$$fieldCounterRef]} ne "")) { #If feature exists - else do nothing
		
		#for my $allele (keys %{${$hashRef}{$$geneRef}}) {  #All alleles
		    
		$$selectedArrayRef[$$fieldCounterRef] .= ",".$$geneRef.":".$allele."|".$$hashRef{$$geneRef}{$allele}{$$fieldArrayRef[$$fieldCounterRef]};
		#}
		#last; #All alleles processed - skip rest of loop
	    }
	}
    }    
}

sub AddToLine {
##Adds to present line
    
    my $fieldsArrayRef = $_[0];
    my $tempArrayRef = $_[1];
    my $lineRef = $_[2];
    my $formatRef = $_[3];
    
    for (my $arrayFieldCounter=0;$arrayFieldCounter<scalar(@{$fieldsArrayRef});$arrayFieldCounter++) {
	
	if (defined($$tempArrayRef[$arrayFieldCounter])) {
	    
	    $$lineRef .= $$tempArrayRef[$arrayFieldCounter];
	}
	else {

	    if ($$formatRef eq "tsv") {

		$$lineRef .= "\t-";
	    }
	}
    }
}

sub ConvertToRange {
##Converts vcf sv to corresponding range coordinates

    my $fieldRef = $_[0];
   
    my $chromosome = $$fieldRef[0];
    my $startPosition = $$fieldRef[1];
    my $referenceAllele = $$fieldRef[3];
    my $alternativeAllele = $$fieldRef[4];

    my $finalStartPosition = $startPosition; #The most "uppstream" position per variant
    my $finalStopPosition = 0; #The most "downstream" position per variant
	
    #Convert to upper case
    ($referenceAllele, $alternativeAllele) = (uc $referenceAllele, uc $alternativeAllele);
    
    if ($alternativeAllele eq ".") { #No Variant Call
	
	next;
    }
    my @alternativeAlleles = split(/,/, $$fieldRef[4]);
    
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
	##Collect largest range per variant based on all alternativeAlleles
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
##Adds arbitrary INFO fields to hash based on supplied headers unless header is already defined

    my $hashRef = $_[0];
    my $rangeFileKey = $_[1];
    my $headerRef = $_[2];
    my $positionRef = $_[3];
    my $rangeFileRef = $_[4];

    if (defined($$hashRef{$rangeFileKey}{$$headerRef})) { #Add INFO from predefined entries

	$$hashRef{'Present'}{$$headerRef}{'INFO'} = $$hashRef{$rangeFileKey}{$$headerRef}{'INFO'};
	$$hashRef{'Present'}{$$headerRef}{'ColumnOrder'} = $$positionRef; #Column position in supplied range input file
    }
    else { #Add arbitrary INFO field using input header

	$$hashRef{'Present'}{$$headerRef}{'INFO'} = q?##INFO=<ID=?.$$headerRef.q?,Number=.,Type=String,Description="String taken from ?.$$rangeFileRef.q?">?;
	$$hashRef{'Present'}{$$headerRef}{'ColumnOrder'} = $$positionRef; #Column position in supplied -sf_ac
    }
}

sub RangeAnnotations {
##Creates the interval tree(s) for range and select files. Adds corresponding INFO fields to metadata.
    
    my $rangeCoulumnsArrayRef = $_[0]; #range or select file annotation columns
    my $lineElementsArrayRef = $_[1]; #range or select line
    my $hashRef = $_[2]; #Hash to store metaData in
    my $rangeFileKey = $_[3];
    my $rangeFileRef = $_[4];
    my $headersArrayRef = $_[5]; #Headers from rangeFile
    my $paddingRef = $_[6];  #Nucleotides to pad the range
    
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
	## Add header to future INFO field
	&AddMetaDataINFO($hashRef, $rangeFileKey, \$$headersArrayRef[ $$rangeCoulumnsArrayRef[$extractColumnsCounter] ], \$extractColumnsCounter, \$$rangeFileRef);
    }
    unless(defined($tree{$rangeFileKey}{ $$lineElementsArrayRef[0] })) { #Only create once per firstKey
	
	$tree{$rangeFileKey}{ $$lineElementsArrayRef[0] } = Set::IntervalTree->new(); #Create tree
    }
    my $paddedStart = $$lineElementsArrayRef[1] - $$paddingRef;
    my $paddedStop = $$lineElementsArrayRef[2] + $$paddingRef;
    $tree{$rangeFileKey}{ $$lineElementsArrayRef[0] }->insert($features, $paddedStart, $paddedStop); #Store range and ";" sep string
}

sub TreeAnnotations {
##Checks if an interval tree exists (per chr) and collects features from input file and adds annotations to line
  
    my $rangeFileKey = $_[0];
    my $lineElementsArrayRef = $_[1]; #Infile (vcf)
    my $hashRef = $_[2];
    my $printLineRef = $_[3]; #Line to add annotations to
    
    if(defined($tree{$rangeFileKey}{ $$lineElementsArrayRef[0] }) ) { #Range annotations
	
	my $feature; #Features to be collected
	
	my ($start, $stop) = &ConvertToRange($lineElementsArrayRef); #Convert SVs to range coordinates from vcf coordinates
	
	if ($start eq $stop) { #SNV
	    
	    $feature =$tree{$rangeFileKey}{ $$lineElementsArrayRef[0] }->fetch($start, $stop+1); #Add 1 to SNV to create range input.
	}
	else {#Range input
	    
	    $feature = $tree{$rangeFileKey}{ $$lineElementsArrayRef[0] }->fetch($start, $stop);
	}
	if (defined(@{$feature})) { #Features found in tree
	    
	    my %collectedAnnotations; #Collect all features before adding to line
	    
	    for (my $featureCounter=0;$featureCounter<scalar(@{$feature});$featureCounter++) { #All features

		my @annotations = split(/;/, @{$feature}[$featureCounter]); #Split feature array ref into annotations
		
		for (my $annotationsCounter=0;$annotationsCounter<scalar(@annotations);$annotationsCounter++) { #All annotations
		    
		    if ($featureCounter == (scalar(@{$feature}-1)) ) { #Last for this feature tuple

			if ($featureCounter == 0) { #First time 
			    
			    $collectedAnnotations{$annotationsCounter} = $annotations[$annotationsCounter];
			}
			else {
			 
			    if(defined($annotations[$annotationsCounter])) {

				$collectedAnnotations{$annotationsCounter} .= ",".$annotations[$annotationsCounter];
			    }
			}
			for my $rangeAnnotation (keys % {$$hashRef{'Present'}}) { #All selected annotations
			    
			    if ($$hashRef{'Present'}{$rangeAnnotation}{'ColumnOrder'} eq $annotationsCounter) { #Correct feature
			
				if ($rangeAnnotation eq "Clinical_db_gene_annotation") {  #Special case, which is global and not gene centric

				    my $databaseAnnotationsRef = [ split(/,/, $collectedAnnotations{$annotationsCounter}) ];  #Create reference

				    ## Collect unique elements from array reference and return array reference with unique elements
				    my $uniqueRef = &UniqElements($databaseAnnotationsRef);
				    $collectedAnnotations{$annotationsCounter} = join(",", @{$uniqueRef});  #Create string with distinct databaseAnnotations
				}
				if ( ($collectedAnnotations{$annotationsCounter} ne "") && ($collectedAnnotations{$annotationsCounter}!~/^,+$/) ) {

				    $$printLineRef .= $rangeAnnotation."=".$collectedAnnotations{$annotationsCounter}.";"; #Add to corresponding line
				}
			    }
			}
		    }
		    elsif ($featureCounter == 0) { #First time but @feature size is greater than 1
			
			$collectedAnnotations{$annotationsCounter} = $annotations[$annotationsCounter];
		    }
		    else {

			$collectedAnnotations{$annotationsCounter} .= ",".$annotations[$annotationsCounter];	
		    }
		}
	    }
	}
    }
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
