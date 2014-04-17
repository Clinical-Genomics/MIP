#!/usr/bin/perl - w

use strict;
use warnings;


use Getopt::Long;
use IO::File;

use vars qw($USAGE);

BEGIN {
    $USAGE =
	qq{vep_parser.pl -i infile.vcf -o outfile.tsv
           -i/--infile Infile (vcf)
           -o/--outfile outfile (tsv)
           -sf/--selectFeatureFile (tsv)
           -sf_mc/--selectFeatureMatchingColumn
           -h/--help Display this help message    
           -v/--version Display version
        };    
}

my ($infile, $outfile, $selectFeatureFile, $selectFeatureMatchingColumn) = (0, 0, 0, "nocmdinput");
my @metaData; 
my (%geneAnnotation, %consequenceSeverity, %selectData);

my ($help, $version) = (0, 0);

###User Options
GetOptions('i|infile:s' => \$infile,
	   'o|outfile:s' => \$outfile,
	   'sf|selectFeatures:s' => \$selectFeatureFile,
	   'sf_mc|selectFeatureMatchingColumn:n' => \$selectFeatureMatchingColumn,
	   'h|help' => \$help,  #Display help text
	   'v|version' => \$version, #Display version number
    );

if($help) {
    
    print STDOUT "\n".$USAGE, "\n";
    exit;
}

if($version) {
    
    print STDOUT "\nvep_parser.pl v1.0", "\n\n";
    exit
}

if ($infile eq 0) {
    
    print STDOUT "\n".$USAGE, "\n";
    print STDERR "\n", "Need to specify  infile by using flag -i","\n\n";
    exit;
}
if ($outfile eq 0) {
    
    print STDOUT "\n".$USAGE, "\n";
    print STDERR "\n", "Need to specify  outfile by using flag -o","\n\n";
    exit;
}
if ( ($selectFeatureMatchingColumn eq "nocmdinput") && ($selectFeatureFile ne 0) ) {
    
    print STDOUT "\n".$USAGE, "\n";
    print STDERR "\n", "Need to specify which feature column to use when selecting variants by using flag -sf_mc","\n\n";
    exit;
}

###
#MAIN
###

if ($selectFeatureFile ne 0) {

    &ReadSelectFile($selectFeatureFile, $selectFeatureMatchingColumn);
}

&DefineConsequenceSeverity();

&ReadInfileVCF($infile, $outfile);

###
#Sub Routines
###

sub DefineConsequenceSeverity {
##Defines the precedence of consequences

    $consequenceSeverity{'transcript_ablation'}{'Rank'} = 1;
    $consequenceSeverity{'transcript_ablation'}{'GeneAnnotation'} = "exonic";
    $consequenceSeverity{'splice_donor_variant'}{'Rank'} = 2;
    $consequenceSeverity{'splice_donor_variant'}{'GeneAnnotation'} = "splicing";
    $consequenceSeverity{'splice_acceptor_variant'}{'Rank'} = 2;
    $consequenceSeverity{'splice_acceptor_variant'}{'GeneAnnotation'} = "splicing";
    $consequenceSeverity{'stop_gained'}{'Rank'} = 3;
    $consequenceSeverity{'stop_gained'}{'GeneAnnotation'} = "exonic";
    $consequenceSeverity{'frameshift_variant'}{'Rank'} = 4;
    $consequenceSeverity{'frameshift_variant'}{'GeneAnnotation'} = "exonic";
    $consequenceSeverity{'stop_lost'}{'Rank'} = 5;
    $consequenceSeverity{'stop_lost'}{'GeneAnnotation'} = "exonic";
    $consequenceSeverity{'initiator_codon_variant'}{'Rank'} = 6;
    $consequenceSeverity{'initiator_codon_variant'}{'GeneAnnotation'} = "exonic";
    $consequenceSeverity{'inframe_insertion'}{'Rank'} = 6;
    $consequenceSeverity{'inframe_insertion'}{'GeneAnnotation'} = "exonic";
    $consequenceSeverity{'inframe_deletion'}{'Rank'} = 6;
    $consequenceSeverity{'inframe_deletion'}{'GeneAnnotation'} = "exonic";
    $consequenceSeverity{'missense_variant'}{'Rank'} = 6;
    $consequenceSeverity{'missense_variant'}{'GeneAnnotation'} = "exonic";
    $consequenceSeverity{'transcript_amplification'}{'Rank'} = 7;
    $consequenceSeverity{'transcript_amplification'}{'GeneAnnotation'} = "exonic";
    $consequenceSeverity{'splice_region_variant'}{'Rank'} = 8;
    $consequenceSeverity{'splice_region_variant'}{'GeneAnnotation'} = "splicing";
    $consequenceSeverity{'incomplete_terminal_codon_variant'}{'Rank'} = 9;
    $consequenceSeverity{'incomplete_terminal_codon_variant'}{'GeneAnnotation'} = "exonic";
    $consequenceSeverity{'synonymous_variant'}{'Rank'} = 10;
    $consequenceSeverity{'synonymous_variant'}{'GeneAnnotation'} = "exonic";
    $consequenceSeverity{'stop_retained_variant'}{'Rank'} = 10;
    $consequenceSeverity{'stop_retained_variant'}{'GeneAnnotation'} = "exonic";
    $consequenceSeverity{'coding_sequence_variant'}{'Rank'} = 11;
    $consequenceSeverity{'coding_sequence_variant'}{'GeneAnnotation'} = "exonic";
    $consequenceSeverity{'mature_miRNA_variant'}{'Rank'} = 12;
    $consequenceSeverity{'mature_miRNA_variant'}{'GeneAnnotation'} = "ncRNA_exonic";
    $consequenceSeverity{'5_prime_UTR_variant'}{'Rank'} = 13;
    $consequenceSeverity{'5_prime_UTR_variant'}{'GeneAnnotation'} = "5UTR";
    $consequenceSeverity{'3_prime_UTR_variant'}{'Rank'} = 14;
    $consequenceSeverity{'3_prime_UTR_variant'}{'GeneAnnotation'} = "3UTR";
    $consequenceSeverity{'non_coding_exon_variant'}{'Rank'} = 15;
    $consequenceSeverity{'non_coding_exon_variant'}{'GeneAnnotation'} = "ncRNA_exonic";
    $consequenceSeverity{'nc_transcript_variant'}{'Rank'} = 15;
    $consequenceSeverity{'nc_transcript_variant'}{'GeneAnnotation'} = "ncRNA";
    $consequenceSeverity{'intron_variant'}{'Rank'} = 16;
    $consequenceSeverity{'intron_variant'}{'GeneAnnotation'} = "intronic";
    $consequenceSeverity{'NMD_transcript_variant'}{'Rank'} = 17;
    $consequenceSeverity{'NMD_transcript_variant'}{'GeneAnnotation'} = "ncRNA";
    $consequenceSeverity{'upstream_gene_variant'}{'Rank'} = 18;
    $consequenceSeverity{'upstream_gene_variant'}{'GeneAnnotation'} = "upstream";
    $consequenceSeverity{'downstream_gene_variant'}{'Rank'} = 19;
    $consequenceSeverity{'downstream_gene_variant'}{'GeneAnnotation'} = "downstream";
    $consequenceSeverity{'TFBS_ablation'}{'Rank'} = 20;
    $consequenceSeverity{'TFBS_ablation'}{'GeneAnnotation'} = "TFBS";
    $consequenceSeverity{'TFBS_amplification'}{'Rank'} = 21;
    $consequenceSeverity{'TFBS_amplification'}{'GeneAnnotation'} = "TFBS";
    $consequenceSeverity{'TF_binding_site_variant'}{'Rank'} = 22;
    $consequenceSeverity{'TF_binding_site_variant'}{'GeneAnnotation'} = "TFBS";
    $consequenceSeverity{'regulatory_region_variant'}{'Rank'} = 22;
    $consequenceSeverity{'regulatory_region_variant'}{'GeneAnnotation'} = "regulatory_region";
    $consequenceSeverity{'regulatory_region_ablation'}{'Rank'} = 23;
    $consequenceSeverity{'regulatory_region_ablation'}{'GeneAnnotation'} = "regulatory_region";
    $consequenceSeverity{'regulatory_region_amplification'}{'Rank'} = 24;
    $consequenceSeverity{'regulatory_region_amplification'}{'GeneAnnotation'} = "regulatory_region";
    $consequenceSeverity{'feature_elongation'}{'Rank'} = 25;
    $consequenceSeverity{'feature_elongation'}{'GeneAnnotation'} = "genomic_feature";
    $consequenceSeverity{'feature_truncation'}{'Rank'} = 26;
    $consequenceSeverity{'feature_truncation'}{'GeneAnnotation'} = "genomic_feature";
    $consequenceSeverity{'intergenic_variant'}{'Rank'} = 27;
    $consequenceSeverity{'intergenic_variant'}{'GeneAnnotation'} = "intergenic" 

}

sub ReadSelectFile {
##Reads a file containg features to be analysed seperately

    my $infileName = $_[0];
    my $selectFeatureColumn = $_[1];

    open(RSF, "<".$infileName) or die "Can't open ".$infileName.":".$!, "\n"; 

    while (<RSF>) {
	
	chomp $_; #Remove newline

	if (m/^\s+$/) {		# Avoid blank lines
	    
	    next;
	}
	if ($_=~/^#/) {#Comment - Avoid
	    
	    next;
	}
	if ( $_ =~/^(\S+)/ ) {	
	    
	    my @lineElements = split("\t",$_); #Loads databse elements description
	    $selectData{$lineElements[$selectFeatureColumn]} = $lineElements[$selectFeatureColumn];
	}
    }
    close(RSF);
    print STDOUT "Finished Reading Select file: ".$infileName,"\n";
}

sub ReadInfileVCF {
#Reads infile vcf format 

    my $infileName = $_[0];
    my $outfileName = $_[1];

    my ($volume,$directories,$file) = File::Spec->splitpath($outfileName);
    my $selectOutFileName = &RemoveFileEnding(\$file, ".vcf");
    $selectOutFileName = $directories.$selectOutFileName.".selected.vcf";

    my @vepFormatField;
    my %vepFormatFieldColumn;

    my @printTSVFields = ("FeatureType");
    my @featureFields;

    open(WOFTSV, ">".$outfileName) or die "Can't open ".$outfileName.":".$!, "\n";

    if ($selectFeatureFile ne 0) {

	open(WOSFTSV, ">".$selectOutFileName) or die "Can't open ".$selectOutFileName.":".$!, "\n";
    }

    open(RIFVCF, "<".$infileName) or die "Can't open ".$infileName.":".$!, "\n";    
    
    while (<RIFVCF>) {
	
	chomp $_; #Remove newline
	
	if (m/^\s+$/) {		# Avoid blank lines
	    next;
	}
	if ($_=~/^##/) {#MetaData

	    push(@metaData, $_); #Save metadata string

	    if ($_=~/INFO\=\<ID\=CSQ/) { #Find VEP INFO Field

		if ($_=~/Format:\s(\S+)"\>/) { #Locate Format within VEP INFO Field
		
		    @vepFormatField = split(/\|/, $1);  
		    
		    for (my $fieldCounter=0;$fieldCounter<scalar(@vepFormatField);$fieldCounter++) {
			
			$vepFormatFieldColumn{$vepFormatField[$fieldCounter]} = $fieldCounter; #Save the order of VEP features
		    }
		}
	    }
	    next;
	}
	if ($_=~/^#CHROM/) {

	    push(@metaData, $_); #Save string
	    $metaData[-1] .= "\tHGNC_transcript_info"."\tFunctional_annotation"."\tGene_annotation";
	    @featureFields = ("MostSevere", "GeneAnnotation");

	    if ($vepFormatFieldColumn{'SIFT'}) {

		$metaData[-1] .="\tSIFT";
		push(@featureFields, "Sift");
	    }
	    if ($vepFormatFieldColumn{'PolyPhen'}) {

		$metaData[-1] .="\tPoly_phen_hdiv";
		push(@featureFields, "PolyPhen");
	    }
	    if (@metaData) { #Print metaData if supplied
		
		for (my $metaDataCounter=0;$metaDataCounter<scalar(@metaData);$metaDataCounter++) {
		    
		    print WOFTSV $metaData[$metaDataCounter],"\n";
		    if ($selectFeatureFile ne 0) {

			print WOSFTSV $metaData[$metaDataCounter],"\n";
		    }
		}
	    }
	}
	if ( $_ =~/^(\S+)/ ) {	
	    
	    my %variantData;
	    my %selectedVariantData;
	    my %consequence;
	    my $variantLine;
	    my $foundSelectedFeatureCounter = 0;
	    my @lineElements = split("\t",$_); #Loads databse elements description
	    
	    my @variantEffects = split(/;/, $lineElements[7]); #Split INFO field
	    
	    for (my $variantEffectCounter=scalar(@variantEffects)-1;$variantEffectCounter>1;) { #Count backwards since CSQ usually is at the end
		
		if ($variantEffects[$variantEffectCounter]=~/CSQ\=/) { #Find CSQ

		    my @transcripts = split(/,/, $'); #Split into transcripts elements'

		    for (my $lineElementsCounter=0;$lineElementsCounter<scalar(@lineElements);$lineElementsCounter++) { #Add until INFO field
			
			if ($lineElementsCounter != 7) { #Do not add INFO field

			    $variantLine .= $lineElements[$lineElementsCounter]."\t";
			}
			else { 
			    
			    for (my $noCSQElements=0;$noCSQElements<scalar(@variantEffects);$noCSQElements++) { #Add all INFO field except CSQ
			
				unless ($variantEffectCounter == $noCSQElements) {
				    
				    $variantLine .=$variantEffects[$noCSQElements].";";
				}
			    }	 
			    $variantLine .= "\t";
			}		    
		    }
		    my $transcriptsCounter = 0; #Tracks the number of transcripts to enable print of ", and ;" at correct position
		    my $selectedTranscriptCounter = 0; #Tracks the number of selected transcripts to enable print of ", and ;" at correct position
		    my $selectedVariantLine = $variantLine; #Copy vcf info to enable selective print downstream

		    for (my $fieldCounter=0;$fieldCounter<scalar(@transcripts);$fieldCounter++) { #CSQ field
			
			my @transcriptsEffects = split(/\|/, $transcripts[$fieldCounter]); #Split in "|"
			
			if ($transcriptsEffects[ $vepFormatFieldColumn{'Feature'} ] =~/N\w_\d+\.\d/) { #RefSeq
			    	    
			    my $selectedTranscriptTracker = 0; #Track if any transcripts belong to selected features
			    
			    if (defined($transcriptsEffects[ $vepFormatFieldColumn{'SYMBOL'} ])) { #Save HGNC Symbol
				
				$variantData{'Symbol'} = $transcriptsEffects[ $vepFormatFieldColumn{'SYMBOL'} ];
				    
				if ($selectData{ $variantData{'Symbol'} }) { #Exists in selected Features
				    
				    $selectedTranscriptTracker = 1; #record belongs to selected Features

				    &AddFieldToElementCounter(\$selectedTranscriptCounter, \$selectedVariantLine, ",",\$variantData{'Symbol'});
				}
				else {

				    &AddFieldToElementCounter(\$transcriptsCounter, \$variantLine, ",",\$variantData{'Symbol'});
				}				
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
				else { #Not in selected genes

				    if ($transcriptsCounter == 0) { #First Gene
					
					$variantData{'FeatureType'} = $transcriptsEffects[ $vepFormatFieldColumn{'Feature_type'} ];
				    }
				    else {
					
					$variantData{'FeatureType'} .= ",".$transcriptsEffects[ $vepFormatFieldColumn{'Feature_type'} ];
				    }
				}
			    }
			    if (defined($transcriptsEffects[ $vepFormatFieldColumn{'Consequence'} ])) { #Save Consequence
				    
				&AddFieldToElement(\$selectedTranscriptTracker, \$selectedVariantLine, \$variantLine, ":", \$transcriptsEffects[ $vepFormatFieldColumn{'Consequence'} ]);
				my @consequences = split(/\&/, $transcriptsEffects[ $vepFormatFieldColumn{'Consequence'} ]); #Find "MostSevere" consequence
				
				for (my $consequencesCounter=0;$consequencesCounter<scalar(@consequences);$consequencesCounter++) {
				    
				    if ( defined($consequence{ $variantData{'Symbol'} }{'Score'}) ) { #Compare to previous record
					
					if ($consequenceSeverity{$consequences[$consequencesCounter]}{'Rank'} < $consequence{ $variantData{'Symbol'} }{'Score'}) { #Collect most severe consequence for Gene_annotation
					    
					    &AddToConsequenceHash(\$consequence{ $variantData{'Symbol'} }{'Score'}, \$consequenceSeverity{$consequences[$consequencesCounter]}{'Rank'});
					    &AddToConsequenceHash(\$consequence{ $variantData{'Symbol'} }{'GeneAnnotation'}, \$consequenceSeverity{$consequences[$consequencesCounter]}{'GeneAnnotation'});
					    &AddToConsequenceHash(\$consequence{ $variantData{'Symbol'} }{'MostSevere'}, \$consequences[$consequencesCounter]);    
					    &AddToConsequenceHash(\$consequence{ $variantData{'Symbol'} }{'Sift'}, \$transcriptsEffects[ $vepFormatFieldColumn{'SIFT'} ]);
					    &AddToConsequenceHash(\$consequence{ $variantData{'Symbol'} }{'PolyPhen'}, \$transcriptsEffects[ $vepFormatFieldColumn{'PolyPhen'} ]);
					}
				    }
				    else { #First pass
					
					&AddToConsequenceHash(\$consequence{ $variantData{'Symbol'} }{'Score'}, \$consequenceSeverity{$consequences[$consequencesCounter]}{'Rank'});
					&AddToConsequenceHash(\$consequence{ $variantData{'Symbol'} }{'GeneAnnotation'}, \$consequenceSeverity{$consequences[$consequencesCounter]}{'GeneAnnotation'});
					&AddToConsequenceHash(\$consequence{ $variantData{'Symbol'} }{'MostSevere'}, \$consequences[$consequencesCounter]);    
					&AddToConsequenceHash(\$consequence{ $variantData{'Symbol'} }{'Sift'}, \$transcriptsEffects[ $vepFormatFieldColumn{'SIFT'} ]);
					&AddToConsequenceHash(\$consequence{ $variantData{'Symbol'} }{'PolyPhen'}, \$transcriptsEffects[ $vepFormatFieldColumn{'PolyPhen'} ]);
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
			    else {
			
				$transcriptsCounter++;
			    }
			    #for (my $transcriptsEffectsCounter=0;$transcriptsEffectsCounter<scalar(@transcriptsEffects);$transcriptsEffectsCounter++) {

				#print WOFTSV $transcriptsEffects[$transcriptsEffectsCounter], "\t";
			    #}
			}
		    }
		    &CollectConsequenceGenes(\%consequence, \@featureFields, \$selectedVariantLine, \$variantLine);
			
		    if ($selectedTranscriptCounter > 0) { #Write to selected file
			
			print WOSFTSV $selectedVariantLine, "\n";
		    }
		    if ($transcriptsCounter > 0) { #Write to transcript file
			
			print WOFTSV $variantLine, "\n";
		    }
		    last; #No need to search variantEfffect longer once found
		}
		
		$variantEffectCounter = $variantEffectCounter -1;
	    }
	    
	}
    }
    close(RIFVCF);
    close(WOFTSV);
    print STDOUT "Finished Reading infile: ".$infileName,"\n";
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
    
    if ($$transcriptCounterRef == 0) {
	
	$$lineRef .= $$valueRef; #First selected feature
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
    else {
	
	$$lineRef .= $separator.$$valueRef;
    }
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
		
		&CollectConsequenceField(\$fieldCounter, \%selectedGeneCounter, \%{$hashRef}, \$genes, \@{$arrayRef}, \@selectedTempFields);
	    }
	    else { #Not selected feature
		
		&CollectConsequenceField(\$fieldCounter, \%geneCounter, \%{$hashRef}, \$genes, \@{$arrayRef}, \@tempFields);
	    }
	}
    }
    
    &AddToLine(\@{$arrayRef}, \@selectedTempFields, \$$selectedVariantLineRef);
    &AddToLine(\@{$arrayRef}, \@tempFields, \$$variantLineRef);
}

sub CollectConsequenceField {
##Collects consequences for features in @featureFields to temporary array for adding to line once all information are collected    

    my $fieldCounterRef = $_[0];
    my $hashCounterRef = $_[1];
    my $hashRef = $_[2];
    my $geneRef = $_[3];
    my $fieldArrayRef = $_[4];
    my $selectedArrayRef = $_[5];
    
    if ($$hashCounterRef{$$fieldCounterRef} == 0) { #First time
	
	if (defined($$hashRef{$$geneRef}{$$fieldArrayRef[$$fieldCounterRef]}) && ($$hashRef{$$geneRef}{$$fieldArrayRef[$$fieldCounterRef]} ne "")) { #If feature exists - else do nothing
	    
	    $$selectedArrayRef[$$fieldCounterRef] .= "\t".$$geneRef.":".$$hashRef{$$geneRef}{$$fieldArrayRef[$$fieldCounterRef]};
	    $$hashCounterRef{$$fieldCounterRef}++;
	}
    }
    else { #Subsequent passes
	
	if (defined($$hashRef{$$geneRef}{$$fieldArrayRef[$$fieldCounterRef]}) && ($$hashRef{$$geneRef}{$$fieldArrayRef[$$fieldCounterRef]} ne "")) { #If feature exists - else do nothing
	    
	    $$selectedArrayRef[$$fieldCounterRef] .= ",".$$geneRef.":".$$hashRef{$$geneRef}{$$fieldArrayRef[$$fieldCounterRef]};
	}
    }    
}

sub AddToLine {
##Adds to present line
    
    my $fieldsArrayRef = $_[0];
    my $tempArrayRef = $_[1];
    my $lineRef = $_[2];
    
    for (my $arrayFieldCounter=0;$arrayFieldCounter<scalar(@{$fieldsArrayRef});$arrayFieldCounter++) {
	
	if (defined($$tempArrayRef[$arrayFieldCounter])) {
	    
	    $$lineRef .= $$tempArrayRef[$arrayFieldCounter];
	}
	else {

	    $$lineRef .= "\t-";
	}
    }
}
