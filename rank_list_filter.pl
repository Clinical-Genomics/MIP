#!/usr/bin/perl - w

#Ranks variants according to weighted sums. Input is a tab-sep file containing all variants per family with each variant recorded on 1 line. The program can process tab-sep list(s) of interesting/dispensible genes, which should either be scored or removed from the final ranked list. A pedigree file describing the family and the family ID should be supplied. In the pedigree file a note on the most likely inheritance models should be noted so that these can be used in the scoring (presently column 11). 
#Commenting one/several subjects in the pedigree file can also be done as long as the -nos flag is set for the total nr of subjects in the master file to ensure proper columns are used.
#Copyright 2012 Henrik Stranneheim

=head1 SYNOPSIS

rank_list_filter.pl -i [infile1] -pedigree [path/file.txt] -o [outfile.txt]

=head2 COMMANDS AND OPTIONS

-i/--infile Infile(s)

-i_gidh/--Infile_Gene_Id_Header Infile geneId Header column to be used in AR_Compound analysis (Defaults to "HGNC_symbol")

-sid_aff/--sampleid_affected Sampleid(s) affected, comma sep

-sid_hea/--sampleid_healthy Sampleid(s) healthy, comma sep

-m/--mother SampleID for mother

-f/--father SampleID for father

-c/--child(ren) Child(ren), comma sep

-sid_male/--sampleid_males Sampleid(s) males, comma sep

-sid_female/--sampleid_females Sampleid(s) females, comma sep

-nos/--numberofsubjects Tracks the number of subjects to make sure that all info stays in proper column

-o/--outfile The output file (defaults to ranked_variants.txt)

-rs/--rankscore The rank score cut-off (defaults to "0")

-im_db/--Im_Db Extra weigth to genes that is in Im_db (Defaults to "0")

-im_db_file/--Im_Db_file Im_Db_CMMS file (Supply whole path)

-im_db_cc/--Im_Db_Gene_Coverage_Calculation Im_Db_CMMS file coverage calculation (Defaults to "0")

-im_db_gidc/--Im_Db_Gene_Id_Col Im_Db_CMMS file gene Id column (Mandatory if -im_db_file is supplied)

-dgf/--dispGeneFiltering Filtering of genes that should be removed from downstream processing (Defaults to "0")

-dgfl/--dispGeneList List of genes that should be removed from downstream processing (Supply whole path, Format: 1 entry per line;HGNC Symbol)

-pedigree/--pedigree_file (Supply whole path)

-familyid/--family Group id of samples to be compared (defaults to "", (Ex: 1 for IDN 1-1-1A))

-recgm/--recGeneticModels (Defaults to "0")

-annovar_dbsnp_ver/--annovar_dbsnp_version Flag for setting the version of dbsnp in annovar (defaults to snp135NonFlagged;Supporded:"132","135","snp135NonFlagged")

-annovar_1000g_ver/--annovar_1000g_version Flag for setting the version of 1000g in annovar (defaults to 1000g2012feb_all)

-tarcov/--targetcoverage Target coverage files for family members, comma sep

=head3 I/O

Input format (annovar_all_variants(tab sep) )

Output format (tab separate list of ranked variants)

=cut

use strict;
use warnings;
use Pod::Usage;
use Pod::Text;
use Getopt::Long;

use vars qw($USAGE);

BEGIN {
    $USAGE =
	qq{rank_list_filter_chr.pl -i [queryfile.txt,,n] -pedigree /pedigree.txt -o outfile.txt
	       -i/--infile Infile(s), comma sep
               -i_gidh/--Infile_Gene_Id_Header Infile geneId Header column to be used in AR_Compound analysis (Defaults to "HGNC_symbol")
               -sid_aff/--sampleid_affected Sampleid(s) affected, comma sep
               -sid_hea/--sampleid_healthy Sampleid(s) healthy, comma sep
               -m/--mother SampleID for mother
               -f/--father SampleID for father
               -c/--children SampleID(s) for Child(ren), comma sep
               -sid_male/--sampleid_males Sampleid(s) males, comma sep
               -sid_female/--sampleid_females Sampleid(s) females, comma sep
               -nos/--numberofsubjects Tracks the number of subjects to make sure that all info stays in proper column
	       -o/--outfile The output file (defaults to ranked_variants.txt)
               -rs/--rankscore The rank score cut-off (defaults to "0")
               -im_db/--Im_Db Extra weigth to genes that is in Im_db (Defaults to "0")
               -im_db_file/--Im_Db_File Im_Db_CMMS file (Supply whole path)
               -im_db_cc/--Im_Db_Gene_Coverage_Calculation Im_Db_CMMS file coverage calculation (Defaults to "0")
               -im_db_gidc/--Im_Db_Gene_Id_Col Im_Db_CMMS file gene Id column Nr (Zero-based;Mandatory if -im_db_cc is supplied)
               -dgf/--dispGeneFiltering Filtering of genes that should be removed from downstream processing (Defaults to "0")
               -dgfl/--dispGeneList List of genes that should be removed from downstream processing (Supply whole path, Format: 1 entry per line;HGNC Symbol)
               -pedigree/--pedigree_file (Supply whole path)
               -familyid/--family Group id of samples to be compared (defaults to "", (Ex: 1 for IDN 1-1-1A))
               -recgm/--recGeneticModels (Defaults to "0")
               -annovar_dbsnp_ver/--annovar_dbsnp_version Flag for setting the version of dbsnp in annovar (defaults to snp135NonFlagged;Supporded:"132","135","snp135NonFlagged")
               -annovar_1000g_ver/--annovar_1000g_version Flag for setting the version of 1000g in annovar (defaults to 1000g2012feb_all)
               -tarcov/--targetcoverage Target coverage files for family members, comma sep
	   };
}

###
#Infile flags (variables)
###
my ($i_gidh, $motherID, $fatherID,$of, $nos, $rankscore, $im_db, $im_db_file, $im_db_cc, $im_db_gidc, $dgf, $dgf_l,$pedigree,$familyid, $annovar_dbsnp_ver, $annovar_1000g_ver, $cmms_imdb, $help) = ("HGNC_symbol", 0,0,"ranked_variants.txt", 0, 0, 0, "",0,"",0,0,0,0,"snp135NonFlagged", "1000g2012feb_all",0);
my (@infn, @childID, @sid_aff, @sid_hea, @sid_male, @sid_female, @samples, @recgm,@tarcovfiles);
my @chr = ("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrMT");

GetOptions('i|infile:s'  => \@infn, #Comma separated list
	   'i_gidh|Infile_Gene_Id_Header:s'  => \$i_gidh, #Header used to find GeneName column to be used in the AR_compound filtering
	   'sid_aff|sampleid_affected:s'  => \@sid_aff, #Comma separated list of affected subjects
	   'sid_hea|samplid_healthy:s'  => \@sid_hea, #Comma separated list of healthy subjects
	   'm|mother:s'  => \$motherID, #SampleID for mother
	   'f|father:s'  => \$fatherID, #SampleID for father
	   'c|children:s'  => \@childID, #SampleID for child(ren), comma separated list
	   'sid_male|sampleid_males:s'  => \@sid_male, #Male subjects, comma separated list
	   'sid_female|sampleid_females:s'  => \@sid_female, #Female subjects, comma separated list
	   'nos|noofsubjects:n'  => \$nos, #Nr of subjects
	   'o|outfile:s'  => \$of,
	   'rsrankscore:n'  => \$rankscore, #The rank score cut-off
	   'im_db|Im_Db:n'  => \$im_db, #Enables important genes scoring
	   'im_db_file|Im_Db_file:s'  => \$im_db_file, #Db of important genes (Supply whole path)
	   'im_db_cc|Im_Db_Gene_Coverage_Calculation:n'  => \$im_db_cc, #Db of important genes coverage calculation (all features connected to overlapping genes across variant)
	   'im_db_gidc|Im_Db_Gene_Id_Col:n'  => \$im_db_gidc, #Db of important genes GeneName column nr zero-based
	   'dgf|dispGeneFiltering:n'  => \$dgf, #Enables dispensible gene filtering
	   'dgfl|dispGeneList:s'  => \$dgf_l, #List of dispensible genes (1 entry per line; HGNC Symbol)
	   'pedigree|pedigree_file:s'  => \$pedigree, #Path to pedigree file location
	   'familyid|familygroup:s' => \$familyid, #Family group ID (FDN)
	   'recgm|recGeneticModels:s'  => \@recgm, #Recommended genetic model
	   'annovar_dbsnp_ver|annovar_dbsnp_version:s' => \$annovar_dbsnp_ver, #dnSNP versions
	   'annovar_1000g_ver|annovar_1000g_version:s' => \$annovar_1000g_ver, #1000G versions
	   'h|help' => \$help,
	   'cmms_imdb|cmms_Db_pathway:n'  => \$cmms_imdb, #Enables CMMS pathways to be added to final list (specific for IEM_Db_CMMS - Secret option)
	   'tarcov|targetcoverage:s'  => \@tarcovfiles, #Target coverage files for members, comma separated list.
	   );

die $USAGE if( $help );

if (@infn == 0) {
   my $verbosity = 2;
 print"\n";
 pod2usage({-message => "Must supply an infile directory as comma separeted list.\n",
     -verbose => $verbosity
   });
}
if ( $im_db == 1 ) {
    unless ($im_db_file) { 
	my $verbosity = 2;
	print"\n";
	pod2usage({-message => "Must supply the important genes file.\n",
	    -verbose => $verbosity
		  });
    }
    if ($im_db_cc == 1) {
	unless ( $im_db_gidc ) {
	    print STDERR "\n";
	    print STDERR "Must supply  the important database GeneId column number to enable coverage calculation (zero-based) ", "\n\n";
	    die $USAGE;
	}
    }
}
if ($im_db_cc == 1) {
    unless ($im_db_file) { 
	my $verbosity = 2;
	print"\n";
	pod2usage({-message => "Must supply the important genes file.\n",
		   -verbose => $verbosity
		  });
    }
    unless ( $im_db_gidc ) {
	print STDERR "\n";
	print STDERR "Must supply  the important database GeneId column number to enable coverage calculation (zero-based) ", "\n\n";
	die $USAGE;
    }
}

if ( $dgf == 1 ) {
    unless ($dgf_l) { 
	my $verbosity = 2;
	print"\n";
	pod2usage({-message => "Must supply a dispensible gene list.\n",
	    -verbose => $verbosity
		  });
    }
}
if ( $familyid eq 0 ) {
    if ($pedigree eq 0) { #Else family ID can be estimated from  pedigree
	print STDERR "\n";
	print STDERR "Must supply a family id or a pedigree file. If not applicable supply the same familyid as the sampleid ", "\n\n";
	die $USAGE;
    }
}

###
#Program Variables
###
my ($compoundc, $pattern_overlap) = (0,0,"");
my (@allVariants, @allVariants_unique, @allVariants_sorted, @ar_comp, @weigthed_scores_seen, @weighted_scores_sorted, @infile_header);
my (%genes, %sid_aff_GT, %sid_hea_GT, %sid_aff_GT_C, %sid_hea_GT_C, %com_var_GT, %filtered, %gene_variants);
my (%allVariants, %allVariants_chr, %allVariants_chr_unique, %allVariants_chr_sorted, %allVariants_geneName_parsed, %ar_comp, %ar_comp_gene, %filtered_ar_comp_weigthed_score, %imdb,%imdb_pos, %dispgenelist,%imdb_tarcov, %weigthed_scores_seen);
my (%pedigree,%rec_genetic_models);


@infn = split(/,/,join(',',@infn)); #Enables comma separated indir(s)
@sid_aff = split(/,/,join(',',@sid_aff)); #Enables comma separated affected IDN(s)
if (@childID) {
    @childID = split(/,/,join(',',@childID)); #Enables comma separated childID(s)
}
if (@sid_male) {
@sid_male = split(/,/,join(',',@sid_male)); #Enables comma separated male IDN(s)
}
if (@sid_female) {
@sid_female = split(/,/,join(',',@sid_female)); #Enables comma separated female IDN(s)
}
if (@sid_hea) {
    @sid_hea = split(/,/,join(',',@sid_hea)); #Enables comma separated unaffected IDN(s)
}
if (@recgm) {
    @recgm = split(/,/,join(',',@recgm)); #Enables comma separated recommended genetic models(s)
    for (my $model=0;$model<scalar(@recgm);$model++) { #Saves supplied models to hash for later comparison
	$rec_genetic_models{$recgm[$model]} = $recgm[$model];
    }
}
if (@tarcovfiles) {
    @tarcovfiles = split(/,/,join(',',@tarcovfiles)); #Enables comma separated unaffected IDN(s)
}

#Collect all supplied subjects info
if (@childID) {
    @samples = @childID
}
if ($motherID){
    push (@samples, $motherID);
}
if ($fatherID){
    push (@samples, $fatherID);
}    
if ($nos == 0 && $pedigree eq 0) {
 print STDERR "\n", "Specify how many subjects that are included or you migt overwrite prescence of subjects.\n";
 print STDERR "\n", "Setting number of subjects to ". scalar(@samples)." change by using flag -nos.\n";
 $nos = scalar(@samples);
}



###
#Column Nr. in annovar_all_variants file (Decided by varcall_merge_post_annovar_master.pl). Used in filtering. If $nos is not supplied and pedigree is then nos is assigned based on the pedigree information and the sub AssignColNr is called first thing in "MAIN" to make sure that correct columns are queried in the filtering.
###
my ($chrcol,$startcol,$stopcol,$refallelecol,$altallelecol);
my ($loccol, $genecol, $hgmdcol, $syncol, $ensemblegeneidcol); 
my ($mce46waycol, $gerpelemcol, $segdupcol);
my ($thGcol, $dbsnp129col);
my ($dbsnpcol, $cg69col);
my ($avsiftcol, $pp2col, $muttascol, $gerpbasecol, $phylopcol);
my ($genmodelcol, $weigthedsumcol);


###
#Genetic Models. Used in Create Models
###
#Collect parents sampleID and disease status    
my $fatherDS = 0; #DS = Disease status, 0=Unaffected 
my $motherDS = 0;
my $childDS = 0;

my ($adom_mom,$adom_father,$adom_child, $adom_model);
my @adom_model;

my ($arecessive_mom,$arecessive_father,$arecessive_child, $arecessive_model);
my @arecessive_model;

my ($xrecessive_mom,$xrecessive_father,$xrecessive_child, $xrecessive_model);
my @xrecessive_model;

my ($denovo_dom_mom,$denovo_dom_father,$denovo_dom_child, $denovo_dom_model);
my @denovo_dom_model;

my ($denovo_rec_mom,$denovo_rec_father,$denovo_rec_child, $denovo_rec_model);
my @denovo_rec_model;

my ($denovo_x_mom,$denovo_x_father,$denovo_x_child, $denovo_x_model);
my @denovo_x_model;

my ($comp_aff_sampleid, $comp_hea_sampleid);
my (@compound_aff_model, @compound_hea_model);

###
#Main
###

if ($fatherID eq 0 && $motherID eq 0 && scalar(@childID) eq 0 && scalar(@sid_aff) eq 0) {
    if ($pedigree eq 0) { 
	print STDERR "\n";
	print STDERR "No information on pedigree\n";
	print STDERR "Must supply location of pedigree file or supply mother/father/Child(ren)/affected/(unaffected) flags", "\n\n";
	die $USAGE;
    }
    else {
	print STDOUT "Estimating information on mother/father/child(ren) affected/unaffected from pedigree file\n";
	ReadPedigreeFile($pedigree);
	AssignColNr($infn[0]); #Note: We lose option of having several infiles with different columns here!
#Create Gentic Models - MEDIUM. For all variants that are either PASS|PRES are included and it only the GT that is important. 
	CreateModels("MEDIUM");
    }
}
elsif ($fatherID || $motherID || scalar(@childID) || scalar(@sid_aff) ) { #No pedigree and user supplied family info on command line
    AssignColNr($infn[0]);		    
    CreateModelsNoPedigree("MEDIUM");
}

if ($im_db_file) {
    ReadImDbCMMS($im_db_file);
    if ($im_db_cc == 1) {    
	for (my $sampleID=0;$sampleID<scalar(@samples);$sampleID++) {
	    for (my $tar_cov_file=0;$tar_cov_file<scalar(@tarcovfiles);$tar_cov_file++) {    
		if ($tarcovfiles[$tar_cov_file]=~/$samples[$sampleID]/) {
		    ReadTargetCov($tarcovfiles[$tar_cov_file],$samples[$sampleID]);
		}
	    }
	}
    }
}

if ($dgf) {
    ReadDispGenes($dgf_l);
}

for (my $chr=0;$chr<scalar(@chr);$chr++) {
    
    if ($chr[$chr+1]) {
	for (my $inputfiles=0;$inputfiles<scalar(@infn);$inputfiles++) {
	    
#AR_compound is a special case and needs to be handled separately
	    if ( $i_gidh eq "HGNC_symbol" ) {
		ReadForARCompAnnovar($infn[$inputfiles],$chr[$chr],$chr[$chr+1]); #Read both files so that all variants are loaded
	    }
	    else {
		ReadForARComp($infn[$inputfiles],$chr[$chr],$chr[$chr+1]); #Read both files so that all variants are loaded
	    }
	}
    }
    else {

	for (my $inputfiles=0;$inputfiles<scalar(@infn);$inputfiles++) {
#AR_compound is a special case and needs to be handled separately
	    if ( $i_gidh eq "HGNC_symbol" ) {
		ReadForARCompAnnovar($infn[$inputfiles],$chr[$chr]); #Read both files so that all variants are loaded
	    }
	    else {
		ReadForARComp($infn[$inputfiles],$chr[$chr]); #Read both files so that all variants are loaded
	    }
	}
    }
    AnalyseCompound();#Filter variants that passes Ar_compound criteria and add gene to hash for checking all variants within gene
    CountGT();
    CompareGT();
    for (my $inputfiles=0;$inputfiles<scalar(@infn);$inputfiles++) {
	ReadVCF($infn[$inputfiles],$chr[$chr]); 
    }
    SortAllVariants();    
    WriteCompound($of,$chr[$chr]);
###
#Blank hashes for next chr
###
    ($compoundc, $pattern_overlap) = (0,0);
    (@allVariants, @allVariants_unique, @allVariants_sorted, @ar_comp, @weigthed_scores_seen, @weighted_scores_sorted) = ();
    (%genes, %sid_aff_GT, %sid_hea_GT, %sid_aff_GT_C, %sid_hea_GT_C, %com_var_GT, %filtered, %gene_variants) = ();
    (%allVariants, %allVariants_chr, %allVariants_chr_unique, %allVariants_chr_sorted, %allVariants_geneName_parsed, %ar_comp, %ar_comp_gene, %filtered_ar_comp_weigthed_score, %weigthed_scores_seen) = ();
}

ReadForFinalSort($of,$of); #Rewrite to same file after reading and sorting whole file


###
#Sub routines
###

sub ReadPedigreeFile {
#Reads famid_pedigree.txt file
#IDN\tSampleID\tMother\tFather\t\Child..n
#$_[0] = filename
    
    open(PEDF, "<$_[0]") or die "Can't open $_[0]:$!, \n";    
    my $temp_fatherID;
    my $temp_motherID;
    my $nrsamples;

    while (<PEDF>) {
	chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
            next;
        }
	if (m/^\#/) {		# Avoid #
            next;
        }		
	if ( ($_ =~/(\S+)/) ) {	
	    chomp($_);
	    my @temp = split("\t",$_);	    #Loads pedigree info
	    if ( $temp[0] =~ /(\d+)-(\d+|-\d+)-(\d+)(A|U)/) { #Match IDN
		$nrsamples++;
		push (@samples, $temp[0]); #Enables commenting one member and running analysis anyway if -nos is supplied
		if ($1) {
		    $familyid = $1;
		}
		if ($3 % 2 == 1) { #Male
#% modulous operator, it gives you an integer representation of the remainder of a division operation. If X / Y divides evenly (that is to say, there's no remainder (or modulus)), the result of X % Y will be zero.
		    $pedigree{$1}{$temp[0]}[0] = "M"; #Sex, M=Male
		    if ($2 eq 2) {
			$temp_fatherID= $temp[0]; #NOTE: No global variable assignation, hence not used in CreateModels. this is only for recording		    
		    }
		}
		else { #Female
		    $pedigree{$1}{$temp[0]}[0] = "F"; #Sex, F=Female
		    if ($2 eq 2) {
			$temp_motherID= $temp[0]; #NOTE: No global variable assignation, hence not used in CreateModels. this is only for recording
		    }
		}
		if ($4 eq "A") { #Affected
		    $pedigree{$1}{$temp[0]}[1] = 1; #1=Affected
		    push (@sid_aff, $temp[0]); #Add affected to enable AR_Compound filtering
		}
		else { #Unaffected
		    $pedigree{$1}{$temp[0]}[1] = 0; #0=Unaffected
		    push (@sid_hea, $temp[0]); #Add healthy to enable AR_Compound filtering
		}
		my @temp_rec_genetic_models = split("\;", $temp[10]); #Loads recommende inheritance models as suggested by MDs
		for (my $model=0;$model<scalar(@temp_rec_genetic_models);$model++) { #Saves models to hash for later comparison
		    $rec_genetic_models{$temp_rec_genetic_models[$model]} = $temp_rec_genetic_models[$model];
		}
		push(@{$pedigree{$1}{$temp[0]}},@temp[2..4,10]); #Populate hash of array for Mother/Father/Child. Full hash: hash{FDN}{IDN}[Sex,Affected,Mother/Father/Child]. Used later in CreateModels subroutine
	    }
	} 	
    }
    close(PEDF);
    print STDOUT "Read Pedigree file: $_[0]","\n";
    print STDOUT "Pedigree information\n";
    if ($nos eq 0) {
	$nos = $nrsamples;
	print STDOUT "Nr of subjects: $nos\n";
    }
    if ($familyid) {
	print STDOUT "Familyid: $familyid\n";
    }
    if ($temp_fatherID) {
	print STDOUT "FatherID: $temp_fatherID\n";
    }
    if ($temp_motherID) {
	print STDOUT "MotherID: $temp_motherID\n";
    }   
    if (@sid_aff) {
	print STDOUT "Affected Subjects:\n";
	for (my $i=0;$i<scalar(@sid_aff);$i++) {
	    print STDOUT "$sid_aff[$i],\t";
	}
	print "\n";
    }
    if (@sid_hea) {
	print STDOUT "Unaffected Subjects:\n";
	for (my $i=0;$i<scalar(@sid_hea);$i++) {
	    print STDOUT "$sid_hea[$i],\t";
	}
	print "\n";
    }
    return;
}

sub AssignColNr {
#$_[0] = First infile

    open(VCF, "<$_[0]") or die "Can't open $_[0]:$!, \n"; 
    
    while (<VCF>) {
	chomp $_;
	if ($_=~/^#/) {
	    my @temp = split("\t",$'); #' 
	    for (my $column=0;$column<scalar(@temp);$column++) {
		if ($temp[$column] eq "Chr") { $chrcol = $column; }
		if ($temp[$column] eq "Start") { $startcol = $column; }
		if ($temp[$column] eq "Stop") { $stopcol = $column; }
		if ($temp[$column] eq "Ref_allele") { $refallelecol = $column; }
		if ($temp[$column] eq "Alt_allele") { $altallelecol = $column; }
		#if ($temp[$column] eq "IDN") { } #Do nothing
		if ($temp[$column] eq $i_gidh) { $genecol = $column; } #Currently used "HGNC_symbol" (121218) in pipe (Annovar output) 
		if ($temp[$column] eq "HGMD") { $hgmdcol = $column; }
		if ($temp[$column] eq "Ensemble_GeneID") { $ensemblegeneidcol = $column; }
		if ($temp[$column] eq "Gene_annotation") { $loccol = $column; }
		if ($temp[$column] eq "Functional_annotation") { $syncol = $column; }
		if ($temp[$column] eq "phastConsElements") { $mce46waycol = $column; }
		if ($temp[$column] eq "gerp++elem") { $gerpelemcol = $column; }
		if ($temp[$column] eq "genomicSuperDups") { $segdupcol = $column; }
		if ($temp[$column] eq "1000G") { $thGcol = $column; }
		if ($temp[$column] eq "dbsnp129") { $dbsnp129col = $column; }
		if ($temp[$column] eq "dbsnp135") { $dbsnpcol = $column; }
		if ($temp[$column] eq "cg69") { $cg69col = $column;} #Currently not used, but supported
		if ($temp[$column] eq "SIFT_Whole-exome") { $avsiftcol = $column; }
		if ($temp[$column] eq "PolyPhen_version_2_HumDiv_Whole-exome") { $pp2col = $column; }
		if ($temp[$column] eq "MutationTaster_Whole-exome") { $muttascol = $column; }
		if ($temp[$column] eq "GERP++_Whole-exome") { $gerpbasecol = $column; }
		if ($temp[$column] eq "PhyloP_Whole-exome") { $phylopcol = $column; }
	    }
#New columns that are to be appended
	    $genmodelcol = scalar(@temp); #0-based
	    $weigthedsumcol =  scalar(@temp)+1; #0-based
	    last;
	}
    }
    close(VCF);
}
sub ReadImDbCMMS {
#Reads Important genes list to enable scoring of these genes. Format: tab-sep, bed-like format 1 gene per line with first four columns: chr, start, stop and GeneName. GeneName should be same as used in the annovar filtering.
#$_[0] = filename (Whole path)
    
    open(IMDB, "<$_[0]") or die "Can't open $_[0]:$!, \n";    
    
    while (<IMDB>) {
        chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
            next;
        }	
	if (m/^\#/) {		# Avoid #
            next;
        }				
	if (/(\S+)/) {  
	    my @temp = split("\t",$_);	    #Loads variant calls
	    if ($im_db_cc == 1) {
		
		if ($temp[$im_db_gidc]) { #If GeneID entry is found
		    my @geneid = split(";", $temp[$im_db_gidc]);
		    for (my $i=0;$i<scalar(@geneid);$i++ ) {
			
			$imdb_pos{$geneid[$i]} = [$temp[0],$temp[1],$temp[2]]; #Add positions
			$imdb{$geneid[$i]} = $geneid[$i]; #Add geneId
			#print $imdb{$geneid[$i]}, "\n";
		    }
		}
	    }
	}
    } 	
    close(IMDB);
    print STDOUT "Read Im_Db file: $_[0]","\n";
    return;
}

sub ReadDispGenes {
#Reads Dispensible genes list to enable removal of these genes from downstream filtering. Format: 1 entry per line with HGNC Symbol or same key as GeneName used in Annovar in the 1st column.
#$_[0] = filename (Whole path)
    
    open(DG, "<$_[0]") or die "Can't open $_[0]:$!, \n";    
    
    while (<DG>) {
        chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
            next;
        }	
	if (m/^\#/) {		# Avoid #
            next;
        }				
	if (/(\S+)/) {  
	    my @temp = split("\t",$_);	    #Loads variant calls
	    
	    if ($temp[0]) { #If entry is found
		$dispgenelist{$temp[0]} = $temp[0]; #Add entry to enable removal. Done in ReadVCF. 
	    }
	}
    } 	
    close(DG);
    print STDOUT "Read dispensible gene list: $_[0]","\n";
    return;
}

sub ReadTargetCov {
#Reads target coverage file per IDN to enable addition to Im_Db entries later in VCF. 
#$_[0] = filename (Whole path)
#$_[1] = sampleID
    
    open(TC, "<$_[0]") or die "Can't open $_[0]:$!, \n";    
    
    while (<TC>) {
        chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
            next;
        }	
	if (m/^\#/) {		# Avoid #
            next;
        }				
	if (/(\S+)/) {  
	    my @temp = split("\t",$_);	    #Loads variant calls
	    $temp[8]=~s/\,/\./;
	    my @temp_EnsembleGeneID = split(";", $temp[11]); #Can have multiple entries due to overlapping genes
	    for (my $i=0;$i<scalar(@temp_EnsembleGeneID);$i++) { #Add all entries to imdb_tarcov
		if ($temp[8]<1) { #If a fraction of bases has below 10 in coverage within feature
		    $temp[8] = sprintf("%.3f",$temp[8]); #Only print three decimal digits	
		    $imdb_tarcov{$_[1]}{$temp_EnsembleGeneID[$i]}{$temp[0]}{$temp[1]}{$temp[2]} = "$_[1]:$temp[0]:$temp[1]:$temp[2]:$temp[8]:"; #Add entry to enable annotation and pathway to add to allvariants hash later. Ensembl Gene ID symbol
		}
		else {
		    $imdb_tarcov{$_[1]}{$temp_EnsembleGeneID[$i]}{$temp[0]}{$temp[1]}{$temp[2]} = "PASS";
		}
	    }
	} 	
    }
    close(TC);
    print STDOUT "Read Target Coverage file: $_[0]","\n";
    return;
}

sub ReadForARComp {
#Reads master genes list to enable AR_compound filtering and GeneName parsing (From X;Y to separate elements ("X","Y","N" ) ), which is ised later in when matching to GeneName is required. Format: tab-sep, 1 variant perl line
#$_[0] = filename (Whole path)
#$_[1] = chr number
#$_[2] = next chr
    
    open(ARC, "<$_[0]") or die "Can't open $_[0]:$!, \n"; #Same as infile. Need to be looped to deduce all genes containing >=2 variants   
    print STDOUT "Reading for $_[1] in $_[0]\n";
    my %gene;
    my @line;

    while (<ARC>) {
        chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
            next;
        }	
	if (m/^\#/) {		# Avoid #
            next;
        }				
	if ($_[2] && ($_ =~/^(\S+)/) && ($1 eq $_[2])) { #If next chr is found return (Since all infiles are sorted this is ok)
	    print STDOUT "Finished Reading $_[1] in Infile $_[0]","\n";
	    close(ARC);
	    last;
	} 
	if ( ($_ =~/^(\S+)/) && ($1 eq $_[1]) ) {#Collect chr number and then if it matches for loop go ahead and collect
	    push(@line, $_); #Store line 
	    my $temp_line = $_;
	    my @temp = split(/\t/,$_);
	    my @element_entries = split("\;",$temp[$genecol]); #Save element entry. Parses ("N;N" entries)
	    for (my $element=0;$element<scalar(@element_entries);$element++) { #For all unique (relative to left side) genes "," i.e. splice and overlapping entries on the right side
		$gene{$element_entries[$element] }++; #Counts all overlapping entries as individual entries
		##Replace overlapping entry with individual geneName
		my @temp_line = split(/\t/,$temp_line); #Local variable
		$temp_line[$genecol] = $element_entries[$element]; #Replace spliced and overlapping entry with individual geneName
		$temp_line = join("\t",@temp_line); #Join the string again
		push(@line, $temp_line); #Add line as individual entry, do not have to remove original entry.
		#   print "Added gene: $splice_genes[$gen]\t";
	    }
	}
    } 	
    close(ARC);
    #print "Nr of lines: ", scalar(@line), "\n";
    for (my $lineCounter=0;$lineCounter<scalar(@line);$lineCounter++) { #For all lines from infile
	my @temp_gene = split(/\t/,$line[$lineCounter]);
	push ( @ {$allVariants_geneName_parsed{$temp_gene[0]}{$temp_gene[1]}{$temp_gene[4]} }, $temp_gene[$genecol]); #Save gene as key to allow multiple entries per chr,start,alt_allel
	#print $allVariants_geneName_parsed{$temp_gene[0]}{$temp_gene[1]}{$temp_gene[4]}{$temp_gene[4+$nos+2]}, "\n"
	if ($gene{$temp_gene[$genecol]} ) { #Original entries containing ";" and "," will be present but should not match (at least very rarely). Exon;splicing events that contain "(" should be correctly added since they should always occur within a gene (relevant cases) and then 1 entry will be correctly added and the line pushed and later compared.   
	    #print "$temp_gene[0]\t $temp_gene[1]\t $temp_gene[4]\t $temp_gene[$genecol]\n";
	    if ($gene{$temp_gene[$genecol]}>1) { #Gene with >=2 variants
		push(@ar_comp, $line[$lineCounter]);
		push ( @ {$ar_comp_gene{$temp_gene[0]}{$temp_gene[1]}{$temp_gene[4]} }, $temp_gene[$genecol]); #Save gene as array to allow multiple entries per chr,start,alt_allel	
		#print "$temp_gene[0]\t $temp_gene[1]\t $temp_gene[4]\t $temp_gene[$genecol]\n";
	    }
	}
    }
    #print "Nr of ar_comp: ", scalar(@ar_comp), "\n";
    print STDOUT "Identified potential AR_compound genes in list: $_[0]","\n";
    return;
}

sub ReadForARCompAnnovar {
#Reads master genes list to enable AR_compound filtering and GeneName parsing (From X;Y,N() to seperate elements ("X","Y","N" ) ), which is ised later in when matching to GeneName is required. Format: tab-sep, 1 variant perl line
#$_[0] = filename (Whole path)
#$_[1] = chr number
#$_[2] = next chr
    
    open(ARC, "<$_[0]") or die "Can't open $_[0]:$!, \n"; #Same as infile. Need to be looped to deduce all genes containing >=2 variants   
    print STDOUT "Reading for $_[1] in $_[0]\n";
    my %gene;
    my @line;

    while (<ARC>) {
        chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
            next;
        }	
	if (m/^\#/) {		# Avoid #
            next;
        }				
	if ($_[2] && ($_ =~/^(\S+)/) && ($1 eq $_[2])) { #If next chr is found return (Since all infiles are sorted this is ok)
	    print STDOUT "Finished Reading $_[1] in Infile $_[0]","\n";
	    close(ARC);
	    last;
	} 
	if ( ($_ =~/^(\S+)/) && ($1 eq $_[1]) ) {#Collect chr number and then if it matches for loop go ahead and collect
	    push(@line, $_); #Store line 
	    my $temp_line = $_;
	    my @temp = split(/\t/,$_);
	    if ($temp[$genecol]=~/\;/) { #If variant is annotated as splicing/overlapping
		my @splice_entries = split("\;",$temp[$genecol]); #Save splicing entry. Parses ("N;N" entries)
		if ( $temp[$genecol]=~/\,/ && $temp[$genecol]!~/\(/) { #If overlapping entries exists on line but not lines with "X(" entries, those are dealt with later
		    #print "$temp[$genecol]\t $splice_entries[0]\t $splice_entries[1]\n";
		    my @splice_genes_left;
		    my @splice_genes_right;
		    if ($splice_entries[0]=~/\,/) { #If splicing entry overlaps multiple genes. (Parses "X,Y..N;" entries)
			@splice_genes_left = split("\,",$splice_entries[0]); #Save each individual gene i.e. "X", "Y" to "N" on the left side of ";"
		    }
		    else {
			$splice_genes_left[0] = $splice_entries[0]; #X; entries
		    }
		    if ($splice_entries[1]=~/\,/) { #If splicing entry overlaps multiple genes. (Parses ";X,Y..N" entries)
			@splice_genes_right = split("\,",$splice_entries[1]); #Save each individual gene i.e. "X", "Y" to "N" on the right side of ";"
		    }
		    else {
			$splice_genes_right[0] = $splice_entries[1]; #;X entries
		    }
		    for (my $gen_left=0;$gen_left<scalar(@splice_genes_left);$gen_left++) { #For all genes containing ; and , i.e. splice and overlapping entries on the left side
			$gene{$splice_genes_left[$gen_left] }++; #Counts all individual splice and overlapping entries as individual entries
			##Replace spliced and overlapping entry with individual geneName
			my @temp_line = split(/\t/,$temp_line); #Local variable
			$temp_line[$genecol] = $splice_genes_left[$gen_left]; #Replace spliced and overlapping entry with individual geneName
			$temp_line = join("\t",@temp_line); #Join the string again
			push(@line, $temp_line); #Add line as individual entry, do not have to remove original entry.
		#	print "Added left side: $splice_genes_left[$gen_left]\t";
			for (my $gen_right=0;$gen_right<scalar(@splice_genes_right);$gen_right++) { #For all genes containing ; and , i.e. splice and overlapping entries on the right side
			    if ($splice_genes_right[$gen_right] && $splice_genes_left[$gen_left] eq $splice_genes_right[$gen_right]) {#remove identical entries on both sides
				$splice_genes_right[$gen_right] ="";
			    } 
			}
		    }
		    for (my $gen_right=0;$gen_right<scalar(@splice_genes_right);$gen_right++) { #For all unique (relative to left side) genes "," i.e. splice and overlapping entries on the right side
			$gene{$splice_genes_right[$gen_right] }++; #Counts all individual splice and overlapping entries as individual entries
			##Replace spliced and overlapping entry with individual geneName
			my @temp_line = split(/\t/,$temp_line); #Local variable
			$temp_line[$genecol] = $splice_genes_right[$gen_right]; #Replace spliced and overlapping entry with individual geneName
			$temp_line = join("\t",@temp_line); #Join the string again
			push(@line, $temp_line); #Add line as individual entry, do not have to remove original entry.
		#	print "Added right side: $splice_genes_right[$gen_right]\t";
		    }
		 #   print "\n";
		}
		elsif ($temp[$genecol]!~/\,/ && $temp[$genecol]!~/\(/) { #Parses X;X, X;Y
		    #print "$temp[$genecol]\t $splice_entries[0]\t $splice_entries[1]\n";
		    if ($splice_entries[0] ne $splice_entries[1]) { #Add both
			$gene{$splice_entries[0] }++; #Counts all individual splice and overlapping entries as individual entries
			##Replace spliced and overlapping entry with individual geneName
			my @temp_line = split(/\t/,$temp_line); #Local variable
			$temp_line[$genecol] = $splice_entries[0]; #Replace spliced and overlapping entry with individual geneName
			$temp_line = join("\t",@temp_line); #Join the string again
			push(@line, $temp_line); #Add line as individual entry, do not have to remove original entry.
		#	print "Added left side: $splice_entries[0]\t";
			$gene{$splice_entries[1] }++; #Counts all individual splice and overlapping entries as individual entries
			##Replace spliced and overlapping entry with individual geneName
			@temp_line = split(/\t/,$temp_line); #Local variable
			$temp_line[$genecol] = $splice_entries[1]; #Replace spliced and overlapping entry with individual geneName
			$temp_line = join("\t",@temp_line); #Join the string again
			push(@line, $temp_line); #Add line as individual entry, do not have to remove original entry.
		#	print "Added right side: $splice_entries[1]\n";
		    }
		    else {
			$gene{$splice_entries[0] }++; #Counts all individual splice and overlapping entries as individual entries
			##Replace spliced and overlapping entry with individual geneName
			my @temp_line = split(/\t/,$temp_line); #Local variable
			$temp_line[$genecol] = $splice_entries[0]; #Replace spliced and overlapping entry with individual geneName
			$temp_line = join("\t",@temp_line); #Join the string again
			push(@line, $temp_line); #Add line as individual entry, do not have to remove original entry.
		#	print "Added left side: $splice_entries[0]\n";
		    }
		}
		if ( $temp[$genecol]=~/\(/) { # All lines with "(" but still splicing ";"
		    my @splice_genes_left;
		    my @splice_genes_right;
		    my @splice_par_left;
		    my @splice_par_right;
		    #print "$temp[$genecol]\t $splice_entries[0]\t $splice_entries[1]\n";		    
		    if ($splice_entries[0]=~/\)\,/) {
			@splice_par_left = split(/\)\,/,$splice_entries[0]); #For X,X(),..,N(),..),N;. Leaves 1. [X,X(] 2. [N]  3. [X,N()] 
			for (my $par_left=0;$par_left<scalar(@splice_par_left);$par_left++) {
			    if ($splice_par_left[$par_left]=~/\(/) {
				$splice_par_left[$par_left] = $`; #Leaves 1. [X,X] 2. [N] untouched  3. [X,N] 
				#print "Spliced par left: ", $splice_par_left[$par_left], "\n";
			    }
			}
			for (my $par_left=0;$par_left<scalar(@splice_par_left);$par_left++) {
			    if ($splice_par_left[$par_left]=~/\,/) {
				@splice_genes_left = split("\,",$splice_par_left[$par_left]); #Leaves 1. [X] [X] 2. [N] 3. [X] [N]
				#print "Spliced genes left: ", $splice_genes_left[0], "\n";
			    }
			    else { #X
				$splice_genes_left[$par_left] = $splice_par_left[$par_left];
				#print "Spliced genes : ", $splice_genes[$par_left], "\n";
			    }
			}
			
		    }
		    elsif ($splice_entries[0]=~/\(/) { #For X,..,X() entries
			$splice_par_left[0] = $`; #Leaves 1. [X,..,X]
			#print "Spliced par left: ", $splice_par_left[0], "\n";
			if ($splice_par_left[0]=~/\,/) {
			    @splice_genes_left = split("\,",$splice_par_left[0]); #Leaves 1. [X] [X]
			 #   print "Spliced genes left: ", $splice_genes_left[0], "\n";
			}
			else { #;X()
			    $splice_genes_left[0] = $splice_par_left[0];
			    #print "Spliced genes left: ", $splice_genes_left[0], "\n";
			}
		    }
		    else { #Parses both X;X() and X,N;X();
			@splice_genes_left = split("\,",$splice_entries[0]);
			for (my $gen_left=0;$gen_left<scalar(@splice_genes_left);$gen_left++) {
			    #print "Spliced genes left: ", $splice_genes_left[$gen_left], "\n";
			}
		    }
		    if ($splice_entries[1]=~/\)\,/) {
			@splice_par_right = split(/\)\,/,$splice_entries[1]); #Leaves 1. X(...  2. X,N(... 
			for (my $par_right=0;$par_right<scalar(@splice_par_right);$par_right++) {
			    if ($splice_par_right[$par_right]=~/\(/) {
				$splice_par_right[$par_right] = $`; #Leaves 1. X   2. X,N
				#print "Spliced par right: ", $splice_par_right[$par_right], "\n";
			    }
			}
			for (my $par_right=0;$par_right<scalar(@splice_par_right);$par_right++) {
			    if ($splice_par_right[$par_right]=~/\,/) {
				@splice_genes_right = split("\,",$splice_par_right[$par_right]); #Leaves 2. X
				#print "Spliced genes right: ", $splice_genes_right[0], "\n";
			    }
			    else { #X
				$splice_genes_right[$par_right] = $splice_par_right[$par_right];
				#print "Spliced genes : ", $splice_genes_right[$par_right], "\n";
			    }
			}
		    }
		    elsif ($splice_entries[1]=~/\(/) { #For ;X,..,X() entries and ;X()
			$splice_par_right[0] = $`; #Leaves 1. [X,..,X]
			#print "Spliced par right: ", $splice_par_right[0], "\n";
			if ($splice_par_right[0]=~/\,/) {
			    @splice_genes_right = split("\,",$splice_par_right[0]); #Leaves 1. [X] [X]
			    for (my $gen_right=0;$gen_right<scalar(@splice_genes_right);$gen_right++) {
				#print "Spliced genes right: ", $splice_genes_right[$gen_right], "\n";
			    }
			}
			else { #;X()
			    $splice_genes_right[0] = $splice_par_right[0];
			    #print "Spliced genes right: ", $splice_genes_right[0], "\n";
			}
		    }
		    else { #Parses both X();X and X();X,N;
			@splice_genes_right = split("\,",$splice_entries[1]);
			for (my $gen_right=0;$gen_right<scalar(@splice_genes_right);$gen_right++) {
			   #print "Spliced genes right: ", $splice_genes_right[$gen_right], "\n";
			}
		    }
		    for (my $gen_left=0;$gen_left<scalar(@splice_genes_left);$gen_left++) { #For all genes containing ; and , i.e. splice and overlapping entries on the left side
			$gene{$splice_genes_left[$gen_left] }++; #Counts all individual splice and overlapping entries as individual entries
			##Replace spliced and overlapping entry with individual geneName
			my @temp_line = split(/\t/,$temp_line); #Local variable
			$temp_line[$genecol] = $splice_genes_left[$gen_left]; #Replace spliced and overlapping entry with individual geneName
			$temp_line = join("\t",@temp_line); #Join the string again
			push(@line, $temp_line); #Add line as individual entry, do not have to remove original entry.
			#print "Added left side: $splice_genes_left[$gen_left]\t";
			for (my $gen_right=0;$gen_right<scalar(@splice_genes_right);$gen_right++) { #For all genes containing ; and , i.e. splice and overlapping entries on the right side
			    if ($splice_genes_right[$gen_right] && $splice_genes_left[$gen_left] eq $splice_genes_right[$gen_right]) {#remove identical entries on both sides
				$splice_genes_right[$gen_right] ="";
			    } 
			}
		    }
		    for (my $gen_right=0;$gen_right<scalar(@splice_genes_right);$gen_right++) { #For all unique (relative to left side) genes "," i.e. splice and overlapping entries on the right side
			$gene{$splice_genes_right[$gen_right] }++; #Counts all individual splice and overlapping entries as individual entries
			##Replace spliced and overlapping entry with individual geneName
			my @temp_line = split(/\t/,$temp_line); #Local variable
			$temp_line[$genecol] = $splice_genes_right[$gen_right]; #Replace spliced and overlapping entry with individual geneName
			$temp_line = join("\t",@temp_line); #Join the string again
			push(@line, $temp_line); #Add line as individual entry, do not have to remove original entry.
			#print "Added right side: $splice_genes_right[$gen_right]\t";
		    }
		    #print "\n";
		}
	    }
	    elsif ( $temp[$genecol]=~/\(/) { # All lines with "(" but no lines with ";". They are already taken care of in previous loop
		my @splice_genes;
		my @splice_par;
		#print "$temp[$genecol]\n";		    
		if ($temp[$genecol]=~/\)\,/) {
		    @splice_par = split(/\)\,/,$temp[$genecol]); #For X,X(),..,N(),..),N. Leaves 1. [X,X(] 2. [N]  3. [X,N()] 
		    for (my $par=0;$par<scalar(@splice_par);$par++) {
			if ($splice_par[$par]=~/\(/) {
			    $splice_par[$par] = $`; #Leaves 1. [X,X] 2. [N] untouched  3. [X,N] 
			    #   print "Spliced par: ", $splice_par[$par], "\n";
			}
		    }
		    for (my $par=0;$par<scalar(@splice_par);$par++) {
			if ($splice_par[$par]=~/\,/) {
			    @splice_genes = split("\,",$splice_par[$par]); #Leaves 1. [X] [X] 2. [N] 3. [X] [N]
			    #  print "Spliced genes : ", $splice_genes[$par], "\n";
			}
			else { #X
			    $splice_genes[$par] = $splice_par[$par];
			    #  print "Spliced genes : ", $splice_genes[$par], "\n";
			}
		    }
		    
		}
		elsif ($temp[$genecol]=~/\(/) { #For X,..,X() entries
		    $splice_par[0] = $`; #Leaves 1. [X,..,X]
		    # print "Spliced par: ", $splice_par[0], "\n";
		    if ($splice_par[0]=~/\,/) {
			@splice_genes = split("\,",$splice_par[0]); #Leaves 1. [X] [X]
			#	print "Spliced genes: ", $splice_genes[0], "\n";
		    }
		    else { #X
			$splice_genes[0] = $splice_par[0];
			#	print "Spliced genes: ", $splice_genes[0], "\n";
		    }
		}
		for (my $gen=0;$gen<scalar(@splice_genes);$gen++) { #For all unique (relative to left side) genes "," i.e. splice and overlapping entries on the right side
		    $gene{$splice_genes[$gen] }++; #Counts all individual splice and overlapping entries as individual entries
		    ##Replace spliced and overlapping entry with individual geneName
		    my @temp_line = split(/\t/,$temp_line); #Local variable
		    $temp_line[$genecol] = $splice_genes[$gen]; #Replace spliced and overlapping entry with individual geneName
		    $temp_line = join("\t",@temp_line); #Join the string again
		    push(@line, $temp_line); #Add line as individual entry, do not have to remove original entry.
		 #   print "Added gene: $splice_genes[$gen]\t";
		}
		#print "\n";
	    }
	    elsif ($temp[$genecol]=~/\,/) { # All lines with "," but no lines with ";" or "(". They are already taken care of in previous loop
		my @splice_genes;
		#print "$temp[$genecol]\n";
		@splice_genes = split("\,",$temp[$genecol]);
		for (my $gen=0;$gen<scalar(@splice_genes);$gen++) { #For all unique (relative to left side) genes "," i.e. splice and overlapping entries on the right side
		    $gene{$splice_genes[$gen] }++; #Counts all individual splice and overlapping entries as individual entries
		    ##Replace spliced and overlapping entry with individual geneName
		    my @temp_line = split(/\t/,$temp_line); #Local variable
		    $temp_line[$genecol] = $splice_genes[$gen]; #Replace spliced and overlapping entry with individual geneName
		    $temp_line = join("\t",@temp_line); #Join the string again
		    push(@line, $temp_line); #Add line as individual entry, do not have to remove original entry.
		    #   print "Added gene: $splice_genes[$gen]\t";
		}
	    }
	    else { #Parses X
		$gene{$temp[$genecol]}++; #Add gene to counter. No unusal entries i.e. just HGNC symbol
		#print "Genes: ", $temp[$genecol], "\n";
		#print $temp_line, "\n";
	    }
	}
    } 	
    close(ARC);
    #print "Nr of lines: ", scalar(@line), "\n";
    for (my $l=0;$l<scalar(@line);$l++) { #For all lines from infile
	my @temp_gene = split(/\t/,$line[$l]);
	push ( @ {$allVariants_geneName_parsed{$temp_gene[0]}{$temp_gene[1]}{$temp_gene[4]} }, $temp_gene[$genecol]); #Save gene as key to allow multiple entries per chr,start,alt_allel
	#print $allVariants_geneName_parsed{$temp_gene[0]}{$temp_gene[1]}{$temp_gene[4]}{$temp_gene[4+$nos+2]}, "\n"
	if ($gene{$temp_gene[$genecol]} ) { #Original entries containing ";" and "," will be present but should not match (at least very rarely). Exon;splicing events that contain "(" should be correctly added since they should always occur within a gene (relevant cases) and then 1 entry will be correctly added and the line pushed and later compared.   
	    #print "$temp_gene[0]\t $temp_gene[1]\t $temp_gene[4]\t $temp_gene[$genecol]\n";
	    if ($gene{$temp_gene[$genecol]}>1) { #Gene with >=2 variants
		push(@ar_comp, $line[$l]);
		push ( @ {$ar_comp_gene{$temp_gene[0]}{$temp_gene[1]}{$temp_gene[4]} }, $temp_gene[$genecol]); #Save gene as array to allow multiple entries per chr,start,alt_allel	
		#print "$temp_gene[0]\t $temp_gene[1]\t $temp_gene[4]\t $temp_gene[$genecol]\n";
	    }
	}
    }
    #print "Nr of ar_comp: ", scalar(@ar_comp), "\n";
    print STDOUT "Identified potential AR_compound genes in list: $_[0]","\n";
    return;
}

sub AnalyseCompound {
#Reads potential AR_compound variants and filters each variant for nonallowed GT combinations.
#$_[0] = filename
	
    my $het_C=0; #To track if all GT are identical 
    my $hom_C=0; #To track if all GT are identical
    my $correct_GT=0; #To track how many GTs are correct per variant
    my $parent_ref_hom=0; #Track is both parents have 0/0 (i.e. AD_denovo and then removes variant from AR_compound list)
    my $parent_ref_het=0; #Track is both parents have 0/0 (i.e. AD_denovo and then removes variant from AR_compound list)
    my $skip=0;
   
    for (my $ar_line=0;$ar_line<scalar(@ar_comp);$ar_line++) { # All potential AR_compound variants
   
	    my @temp = split("\t",$ar_comp[$ar_line]);	    #Loads variant line
	    
	    for (my $sid=0;$sid<scalar( @sid_aff );$sid++)  { #All affected subjects within the same gene
		if ($ar_comp[$ar_line] =~ /($sid_aff[$sid])\:\S+\:(\d+\/\d+)/ ) { #Find affected GT (Genotype) per sampleID
		    if ($2 eq "0/0") { #No need to add variant since within a family it is unlikely that there are two different heterozygous compounds mutations causing the disease between siblings/parent, i.e. even if other affected members have a variant at that position it is still not valid. 
			$skip++; #Not relevant variant, skip it
			next;
		    }
		    if ($2 eq "1/1") {
			if ( ($temp[0] eq "chrX") || ($temp[0] eq "X") ) { #Check for chrX
			    if ( ($sid_aff[$sid] =~ /(\d+)-(\d+|-\d+)-(\d+)(A|U)/) ) {#Match sampleID
				if ($3 % 2 == 1) { #Male#% modulous operator, it gives you an integer representation of the remainder of a division operation. If X / Y divides evenly (that is to say, there's no remainder (or modulus)), the result of X % Y will be zero.
				    if ($temp[1] > 3522410 && $temp[1] <148858526 ) { # 3522410 = Upper bound for PAR1. 148858526 = Lower bound for PAR2. Hence, variants with location outside of PAR1 or PAR2 should be skiped.
					$skip++; #No heterozygous compound on chrX outside of PAR regions for affected males
				    }
				}
			    }
			}
			else {
			    $correct_GT++;
			    $hom_C++;
			}
		    }
		    if ($2 eq "0/1") {
			if ( ($temp[0] eq "chrX") || ($temp[0] eq "X") ) { #Check for chrX
			    if ( ($sid_aff[$sid] =~ /(\d+)-(\d+|-\d+)-(\d+)(A|U)/) ) {#Match sampleID
				if ($3 % 2 == 1) { #Male#% modulous operator, it gives you an integer representation of the remainder of a division operation. If X / Y divides evenly (that is to say, there's no remainder (or modulus)), the result of X % Y will be zero.
				    if ($temp[1] > 3522410 && $temp[1] <148858526 ) { # 3522410 = Upper bound for PAR1. 148858526 = Lower bound for PAR2. Hence, variants with location outside of PAR1 or PAR2 should be skiped.
					$skip++; #No heterozygous compound on chrX outside of PAR regions for affected males
				    }
				}
			    }
			}
			else {
			    $correct_GT++;
			    $het_C++;
			}
		    }
		}
	    }
	    if ($skip) { #No need to continue with variant
		$het_C=0;
		$hom_C=0;
		$correct_GT=0;
		$skip=0;
		next;
	    }
	    for (my $sid=0;$sid<scalar( @sid_hea );$sid++)  { #All healthy subjects within the same gene
		if ($ar_comp[$ar_line] =~ /($sid_hea[$sid])\:\S+\:(\d+\/\d+)/ ) { #Find healthy GT (Genotype) per sampleID
		    
		    if ($2 eq "1/1") { #No need to add variant since a healthy subject with a homozygous variant cannot be involved in a genetic compound.    
			$skip++;
		    }
		    if ($2 eq "0/0") {
			if ( ($sid_hea[$sid] =~ /(\d+)-(\d+|-\d+)-(\d+)(A|U)/) ) {#Match sampleID
			    if ($2 eq 2) { #Parent
				$parent_ref_hom++;
				if ($parent_ref_hom eq 2) { #Both parent have GT: 0/0 and cannot contribute to AR_compound
				    $skip++;
				}
			    }
			}
			$correct_GT++;
		    }
		    if ($2 eq "0/1") {
			if ( ($sid_hea[$sid] =~ /(\d+)-(\d+|-\d+)-(\d+)(A|U)/) ) {#Match sampleID
			    if ($2 eq 2) { #Parent
				$parent_ref_het++;
				if ($parent_ref_het eq 2) { #Both parent have GT: 0/1 and cannot contribute to AR_compound. If affected is 0/1 all GT:s are identical and not a true AR_compound. If affected is 1/1 then it is a AR mutation and should be recorded as such. 
				    $skip++;
				}
			    }
			}
			$correct_GT++;
			$het_C++;
		    }
		}	
	    }	
	    if ($skip) { #No need to continue with variant
		$het_C=0;
		$hom_C=0;
		$correct_GT=0;
		$parent_ref_hom=0;
		$parent_ref_het=0;
		$skip=0;
		next;
	    }
	    if (@sid_hea) {#Only valid test if a healthy subject is included
		if ( ($nos eq $het_C) || ($nos eq $hom_C) ) { #All subjects have the same GT, hence cannot contribute to genetic compound
		    $skip++;
		}
	    }
	    if ($skip) {
		#print "Not Added Subjects Identical GT: $temp[0]\t", "$temp[1]\t", "$temp[4]\t", "$temp[4+$nos+2]\t", "HetC: $het_C\t", "HomC: $hom_C", "\n";
		$het_C=0;
		$hom_C=0;
		$correct_GT=0;
		$parent_ref_hom=0;
		$parent_ref_het=0;
		$skip=0;
		next;
	    }
	    if ($dispgenelist{$temp[$genecol]}) { #Do not add dispensible genes
		$het_C=0;
		$hom_C=0;
		$correct_GT=0;
		$parent_ref_hom=0;
		$parent_ref_het=0;
		$skip=0;
		next;
	    }
	    if ( $nos eq $correct_GT ) { #All subjects have the correct GT, hence can contribute to genetic compound
		push ( @ {$genes{$temp[$genecol]} },$ar_comp[$ar_line]); #Store line in array
		$ar_comp{$temp[0]}{$temp[1]}{$temp[4]} = [@temp]; # Hash{chr}{pos}{variant}, all variants non overlapping and array [chr->unknown] All info starting from chr
		#print "$ar_comp{$temp[0]}{$temp[1]}{$temp[4]}[1]\t", $ar_comp{$temp[0]}{$temp[1]}{$temp[4]}[$genecol], "\n";
	    }
	    $het_C=0;
	    $hom_C=0;
	    $correct_GT=0;
	    $parent_ref_hom=0;
	    $parent_ref_het=0;
	    $skip=0;
    }
    print STDOUT "Filtered potential AR_compound genes per variant","\n";
    return;
}

sub CountGT {
###
#Determine genotype for every affected and healthy sampleID. All variants for all GT:s (genotype(s):0/1,1/1,0/0,NA) are stored by variant Nr within the gene.GT:s are counted for affected subjects. All variants muste be assigned a GT to ensure correct numbering for both affected and healthy subjects.
###    
    my $gender=""; #F eq female and M eq male

    for my $gene (keys %genes) { #Every gene
	
	for (my $sid=0;$sid<scalar( @sid_aff );$sid++)  { #All affected subjects within the same gene
	    
	    for (my $i=0;$i<scalar( @{$genes{$gene} } );$i++) { #All variants within the same gene for all affected sampleID
		
		if ($genes{$gene}[$i] =~ /$sid_aff[$sid]\:\S+\:(\d+\/\d+)/ ) { #Find affected GT (Genotype) per sampleID
		    my $gt = $1;
		    
		    $sid_aff_GT{$gene}{$sid_aff[$sid]}[$i] = $gt; #Save GT (Genotype) for gene, sampleID, and variant number NOTE: Not variant position
		    $sid_aff_GT_C{$gene}{$sid_aff[$sid]}{$gt}++; #Count specific GT for gene, sampleID   
		}
		else {
		    $sid_aff_GT{$gene}{$sid_aff[$sid]}[$i] = "Na"; #Save uncalled GT for gene, sampleID and variant number NOTE: Not variant position
		}
	    }
	}
	if ( @sid_hea ) { #If any healthy is supplied
	    for (my $sid=0;$sid<scalar( @sid_hea );$sid++)  { #All healthy subjects within the same gene
		
		for (my $i=0;$i<scalar( @{$genes{$gene} } );$i++) { #All variants within the same gene for all healthy samplID 
		    if ($genes{$gene}[$i] =~ /$sid_hea[$sid]\:\S+\:(\d+\/\d+)/ ) { #Find healthy GT (Genotype) per sampleID
			$sid_hea_GT{$gene}{$sid_hea[$sid]}[$i] = $1; #Save GT (Genotype) for gene, sampleID, and variant number NOTE: Not variant position
			if ($1 eq "0/0") {
			    $sid_hea_GT_C{$gene}{$sid_hea[$sid]}{$1}++;
			}
		    }
		    else {
			$sid_hea_GT{$gene}{$sid_hea[$sid]}[$i] = "Na"; #Save uncalled GT for gene, sampleID and variant number NOTE: Not variant position
		    }
		}
	    }	 	    
	}
    }
    return;
}

sub CompareGT {
###
#Compare overlap between affected and healthy per affected sampleID. If no healthy subjects are included - then all variants passed to CompareGT will be passed since no addtional information is supplied to make a pattern match. 
###
    my $hetC=0;
    my $homC=0;
    my $No_hetcomp=0; #Tracks healthy subjects GT(0/0)
    my $skip=0; #Skip adding gene to filtered
    
    for my $gene (keys %sid_aff_GT) { #Every gene found containing variants within affected subjects
	
	for my $sampleid_aff (keys %{ $sid_aff_GT{$gene} } )  { #All affected subjects sampleIDs where a variation has been found within the gene
	    for (my $i=0;$i<scalar( @{$sid_aff_GT{$gene}{$sampleid_aff} } );$i++) { #All variants within the gene per affected samplID (Same as total number of variants within gene but some might be "Na")
		#if ($gene eq "DUSP23") {
		#print STDOUT "\nAffected\n Gene $gene\t SAmpleID $sampleid_aff\t";
		#print STDOUT "Variant $sid_aff_GT{$gene}{$sampleid_aff}[$i]\n";
		#}
		if (@sid_hea) { #If any healthy is supplied - check and count non-relevant variants
		    for (my $sampleid_hea=0;$sampleid_hea<scalar( @sid_hea );$sampleid_hea++)  { #For all healthy sampleID within the gene (Same as total number of variants within gene but some might be "Na")
			if ($sid_aff_GT{$gene}{$sampleid_aff}[$i] eq $sid_hea_GT{$gene}{$sid_hea[$sampleid_hea]}[$i] ){ #Check for existence on same variant, does not matter if het or hom		    
			    $com_var_GT{$gene}{$sampleid_aff}{$sid_hea[$sampleid_hea]}++; #Store number of common variants hash{geneName}{sampleID_affected}{sampleID_healthy}
			}
		    }
		}
	    }
	    if (@sid_hea) { #If any healthy is supplied - check that there is at least 1 non-overlapping variant of all variants within affected subject compared to healthy subject. Additionally there must be a total of either two het or two hom variations in the gene for the affected subject.
		for (my $sampleid_hea=0;$sampleid_hea<scalar( @sid_hea );$sampleid_hea++)  { #For all healthy sampleID within the same gene (Same as total number of variants within gene but some might be "Na")
		    my $nr_of_aff_variants = scalar( @{$sid_aff_GT{$gene}{$sampleid_aff} } ); #Nr of variants within gene for affected sampleID (All variants)
		    my $nr_of_overlapping_hea_variants;
		    #Determine the nr of overlapping affected variants within gene between affected subject and healthy subjects
		    if ($com_var_GT{$gene}{$sampleid_aff}{$sid_hea[$sampleid_hea]}){ #If any common variant(s) or not relevant variant(s)
			$nr_of_overlapping_hea_variants = $com_var_GT{$gene}{$sampleid_aff}{$sid_hea[$sampleid_hea]};
		    }
		    else {
			$nr_of_overlapping_hea_variants = 0;
		    }
		    if ( ($nr_of_aff_variants - $nr_of_overlapping_hea_variants)>0 ) { #Require at least 1 non-overlapping variant between affected subject and healthy subjects	
			if ( ($sid_aff_GT_C{$gene}{$sampleid_aff}{"0/1"}) && $sid_aff_GT_C{$gene}{$sampleid_aff}{"0/1"} >0) { #If entry exists and there must be at least two het i.e. hom/het, Na/het is not interesting for affected		
			    $hetC = 1;		
			    if ($sid_aff_GT_C{$gene}{$sampleid_aff}{"0/1"} >1) {
				$hetC = 2;
			    }
			    #$compoundc++;
			}
			if ( ($sid_aff_GT_C{$gene}{$sampleid_aff}{"1/1"}) && $sid_aff_GT_C{$gene}{$sampleid_aff}{"1/1"} >0 ) { #If entry exists and there must be at least two hom i.e. hom/het is not interesting for affected
			    $homC = 1;
			    if ($sid_aff_GT_C{$gene}{$sampleid_aff}{"1/1"} >1) {
				$homC = 2;
			    }
			    #$compoundc++;
			}
			if ( ($hetC > 1) || ($homC > 1) ) {
			    $compoundc++;
			}
			elsif ($hetC + $homC ==2) {
			    $compoundc++;
			}
		    }
		    else { #Did not find at least 1 non-overlapping variant between affected subject and healthy subjects
			$pattern_overlap++;
		    }
		    if ( ($motherID) || ($fatherID) ) { #Only for parents
			if ( ( $sid_hea[$sampleid_hea] eq $motherID) || ( $sid_hea[$sampleid_hea] eq $fatherID) ) { #Mother or Father
			    if ( $sid_hea_GT_C{$gene}{$sid_hea[$sampleid_hea]}{"0/0"} && (scalar( @{$sid_aff_GT{$gene}{$sampleid_aff} } ) eq  $sid_hea_GT_C{$gene}{$sid_hea[$sampleid_hea]}{"0/0"}) ) { #If a parent do not have any GT ne 0/0 for the same number of variants as the affected subject then it cannot be a true ar_compound and should be skiped. 
				$No_hetcomp++;
				$skip++; #No GT for healthy and same total as affected
			    }
			}
		    }
		}
	    }
	}
	if ($skip) { #Blank counters and proceed to next variant
	    $hetC = 0; #Reset for next affected subject
	    $homC = 0; #Reset for next affected subject
	    $compoundc = 0; #Reset for next affected subject
	    $pattern_overlap = 0; #Reset for next affected subject
	    $No_hetcomp = 0; #Reset for next affected subject
	    $skip = 0; #Reset for next affected subject
	    next;
	}
	if  (@sid_hea) { #Only when healthy subjects are supplied
	    if ( ($compoundc >0 && ($pattern_overlap == 0) ) ) { #At least 1 compound i.e. 2 het or 2 hom and at least 2 variants with either 0/1 or 1/1.	
		$filtered{$gene} = $gene; #Gene passes filtering criteria
	    }
	}
	else { #No pattern check has been made since all subjects are affected i.e. all variants passed to subroutine are passed
	    $filtered{$gene} = $gene; #Gene passes filtering criteria
	}
	$hetC = 0; #Reset for next affected subject
	$homC = 0; #Reset for next affected subject
	$compoundc = 0; #Reset for next affected subject
	$pattern_overlap = 0; #Reset for next affected subject
	$No_hetcomp = 0; #Reset for next affected subject
    }
    return;
}

sub CreateModels {
#Decide genotype combinations for family depending on model, affected/healthy, and pedigree. AR_compound is handled separately
#$_[0] = MEDIUM (Adapted from wgs_var_call_wf.1.4.pl)
#MEDIUM requires either PASS or PRES and correct genotype
    
#Blank all arrays and variants to make sure that correct model is created
    @adom_model = ();
    @arecessive_model = ();
    @xrecessive_model = ();
    @denovo_dom_model = ();
    @denovo_rec_model = ();
    @denovo_x_model = ();
    (@compound_aff_model, @compound_hea_model) = ();
    
    $fatherID="";
    $motherID="";
    $comp_aff_sampleid="";
    $comp_hea_sampleid="";
    $fatherDS = 0; #DS = Disease status, 0=Unaffected 
    $motherDS = 0;
    for my $sampleid (keys % { $pedigree{$familyid} } ) {
	
	if ($pedigree{$familyid}{$sampleid}[2] == 1 ) {#Mother
	    $motherID=$sampleid;
	    if ($pedigree{$familyid}{$sampleid}[1] == 1 ) {#Affected
		$motherDS = 1;
		if ($_[0] eq "MEDIUM") { #Either PASS or PRES genotype calls for all samples
		    $adom_mom = $sampleid.q?:\S+\:0\/1?; #Het
		    $arecessive_mom = $sampleid.q?:\S+\:1\/1?; #Hom alt
		    $xrecessive_mom = $sampleid.q?:\S+\:0\/1?; #Het
		    $denovo_x_mom = $sampleid.q?:\S+\:0\/1?; #Het
		    push(@compound_aff_model,$sampleid);
		}
	    }
	    else { #healthy
		if ($_[0] eq "MEDIUM") {
		    $adom_mom = $sampleid.q?:\S+\:0\/0?; #Hom ref
		    $arecessive_mom = $sampleid.q?:\S+\:0\/1?; #Het
		    $xrecessive_mom = $sampleid.q?:\S+\:0\/(1|0)?; #Het/Hom ref
		    $denovo_dom_mom = $sampleid.q?:\S+\:0\/0?; #Hom ref
		    $denovo_rec_mom = $sampleid.q?:\S+\:0\/(1|0)?; #Het/Hom ref
		    $denovo_x_mom = $sampleid.q?:\S+\:0\/(1|0)?; #Het/Hom ref
		    push(@compound_hea_model,$sampleid);
		}
	    }
	    push(@adom_model,$adom_mom);
	    push(@arecessive_model,$arecessive_mom);
	    push(@xrecessive_model,$xrecessive_mom);
	    push(@denovo_dom_model,$denovo_dom_mom);
	    push(@denovo_rec_model,$denovo_rec_mom);
	    push(@denovo_x_model,$denovo_x_mom);
	}
	elsif ($pedigree{$familyid}{$sampleid}[3] == 1 ) {#Father
	    $fatherID=$sampleid;
	    if ($pedigree{$familyid}{$sampleid}[1] == 1 ) {#Affected
		$fatherDS = 1;
		if ($_[0] eq "MEDIUM") {
		    $adom_father = $sampleid.q?:\S+\:0\/1?; #Het
		    $arecessive_father = $sampleid.q?:\S+\:1\/1?; #Hom alt
		    $xrecessive_father = $sampleid.q?:\S+\:1\/1?; #Het (Male: GT is 1/1 for het on X)
		    push(@compound_aff_model,$sampleid);
		}
	    }
	    else { #healthy

		if ($_[0] eq "MEDIUM") {
		    $adom_father = $sampleid.q?:\S+\:0\/0?; #Hom ref
		    $arecessive_father = $sampleid.q?:\S+\:0\/1?; #Het
		    $xrecessive_father = $sampleid.q?:\S+\:0\/0?; #Hom ref
		    $denovo_dom_father = $sampleid.q?:\S+\:0\/0?; #Hom ref
		    $denovo_rec_father = $sampleid.q?:\S+\:0\/(1|0)?; #Het/Hom ref, AR_denovo-1,2
		    $denovo_x_father = $sampleid.q?:\S+\:0\/0?; #Hom reference
		    push(@compound_hea_model,$sampleid);
		}
	    }
	    push(@adom_model,$adom_father);
	    push(@arecessive_model,$arecessive_father);
	    push(@xrecessive_model,$xrecessive_father);
	    push(@denovo_dom_model,$denovo_dom_father);
	    push(@denovo_rec_model,$denovo_rec_father);
	    push(@denovo_x_model,$denovo_x_father);	    
	}
	else {#Child

	    if ($pedigree{$familyid}{$sampleid}[1] == 1 ) {#Affected

		if ($_[0] eq "MEDIUM") {
		    $arecessive_child = $sampleid.q?:\S+\:1\/1?; #Hom alt
		    $denovo_dom_child = $sampleid.q?:\S+\:0\/1?; #Het
		    $denovo_rec_child = $sampleid.q?:\S+\:1\/1?; #Hom alt
		    if ( ($fatherDS == 1) && ($motherDS == 1) ) { #Parents affected
			$adom_child = $sampleid.q?:\S+\:(0|1)\/1?; #Het/hom alt
		    }
		    else {
			$adom_child = $sampleid.q?:\S+\:0\/1?; #Het, AD-1,2 (Common)
		    }
		    push(@compound_aff_model,$sampleid);
		    if ($pedigree{$familyid}{$sampleid}[0]=~/M/i) { #Affected and Male			
			$xrecessive_child = $sampleid.q?:\S+\:1\/1?; #Het (Male: GT is 1/1 for het on X)		    
			$denovo_x_child = $sampleid.q?:\S+\:1\/1?; #Het (Male: GT is 1/1 for het on X)
		    }
		    if ($pedigree{$familyid}{$sampleid}[0]=~/F/i) { #Affected and Female
			$xrecessive_child = $sampleid.q?:\S+\:(0|1)\/1?; #Het/hom alt
			$denovo_x_child = $sampleid.q?:\S+\:(0|1)\/1?; #Het/hom alt (Hom more unlikely) )
		    }
		}
	    }
	    else { #healthy

		if ($_[0] eq "MEDIUM") {
		    $adom_child = $sampleid.q?:\S+\:0\/0?; #Hom ref
		    $arecessive_child = $sampleid.q?:\S+\:0\/(1|0)?; #Het/Hom ref
		    $denovo_dom_child = $sampleid.q?:\S+\:0\/0?; #Hom ref
		    $denovo_rec_child = $sampleid.q?:\S+\:0\/(1|0)?; #Het/Hom ref, AD_denovo-1,2
		    push(@compound_hea_model,$sampleid);
		    if ( $pedigree{$familyid}{$sampleid}[0]=~/M/i) { #Healthy and Male
			$xrecessive_child = $sampleid.q?:\S+\:0\/0?; #Hom (Male: GT is 1/1 for het on X)
			$denovo_x_child = $sampleid.q?:\S+\:0\/0?; #Hom (Male: GT is 1/1 for het on X)
		    }
		    if ($pedigree{$familyid}{$sampleid}[0]=~/F/i) { #Healthy and Female
			$xrecessive_child = $sampleid.q?:\S+\:0\/(1|0)?; #Het/Hom ref
			$denovo_x_child = $sampleid.q?:\S+\:0\/(1|0)?; #Het/Hom ref
		    }
		}
	    }
	    push(@adom_model,$adom_child);
	    push(@arecessive_model,$arecessive_child); 
	    push(@xrecessive_model,$xrecessive_child);
	    push(@denovo_dom_model,$denovo_dom_child);
	    push(@denovo_rec_model,$denovo_rec_child);
	    push(@denovo_x_model,$denovo_x_child);
	}
    }   
}

sub CreateModelsNoPedigree {
#Decide genotype combinations for family depending on model, affected/healthy, and pedigree. AR_compound is handled separately
#$_[0] = MEDIUM (Adapted from wgs_var_call_wf.1.4.pl)
#MEDIUM requires either PASS or PRES and correct genotype    
    
#Blank all arrays and variants to make sure that correct model is created
    @adom_model = ();
    @arecessive_model = ();
    @xrecessive_model = ();
    @denovo_dom_model = ();
    @denovo_rec_model = ();
    @denovo_x_model = ();
    (@compound_aff_model, @compound_hea_model) = ();
    
    $comp_aff_sampleid="";
    $comp_hea_sampleid="";        
    print STDOUT "Building Genetic Models for:\n";
    for (my $sampleid=0;$sampleid<scalar(@samples);$sampleid++ ) {
	$fatherDS = 0; #DS = Disease status, 0=Unaffected 
	$motherDS = 0; 	
	$childDS = 0;
	if ($samples[$sampleid] eq $motherID ) {#Mother
	    for (my $i=0;$i<scalar(@sid_aff);$i++) { #All affected samples	    
		if ($samples[$sampleid] eq $sid_aff[$i] ) {#Affected
		    $motherDS = 1;
		    print STDOUT "Affected Mother: $samples[$sampleid]\n";
		    if ($_[0] eq "MEDIUM") { #Either PASS or PRES genotype calls for all samples
			$adom_mom = $samples[$sampleid].q?:\S+\:0\/1?; #Het
			$arecessive_mom = $samples[$sampleid].q?:\S+\:1\/1?; #Hom alt
			$xrecessive_mom = $samples[$sampleid].q?:\S+\:0\/1?; #Het
			$denovo_x_mom = $samples[$sampleid].q?:\S+\:0\/1?; #Het
			push(@compound_aff_model,$samples[$sampleid]);
		    }
		}
	    }
	    if ($motherDS == 0) { #healthy
		print STDOUT "Unaffected Mother: $samples[$sampleid]\n";
		if ($_[0] eq "MEDIUM") {
		    $adom_mom = $samples[$sampleid].q?:\S+\:0\/0?; #Hom ref
		    $arecessive_mom = $samples[$sampleid].q?:\S+\:0\/1?; #Het
		    $xrecessive_mom = $samples[$sampleid].q?:\S+\:0\/(1|0)?; #Het/Hom ref
		    $denovo_dom_mom = $samples[$sampleid].q?:\S+\:0\/0?; #Hom ref
		    $denovo_rec_mom = $samples[$sampleid].q?:\S+\:0\/(1|0)?; #Het/Hom ref
		    $denovo_x_mom = $samples[$sampleid].q?:\S+\:0\/(1|0)?; #Het/Hom ref
		    push(@compound_hea_model,$samples[$sampleid]);
		}
	    }
	    push(@adom_model,$adom_mom);
	    push(@arecessive_model,$arecessive_mom);
	    push(@xrecessive_model,$xrecessive_mom);
	    push(@denovo_dom_model,$denovo_dom_mom);
	    push(@denovo_rec_model,$denovo_rec_mom);
	    push(@denovo_x_model,$denovo_x_mom);
	}
	elsif ($samples[$sampleid] eq $fatherID ) {#Father
	    for (my $i=0;$i<scalar(@sid_aff);$i++) { #All affected samples	    
		if ($samples[$sampleid] eq $sid_aff[$i] ) {#Affected
		    $fatherDS = 1;
		    print STDOUT "Affected Father: $samples[$sampleid]\n";
		    if ($_[0] eq "MEDIUM") {
			$adom_father = $samples[$sampleid].q?:\S+\:0\/1?; #Het
			$arecessive_father = $samples[$sampleid].q?:\S+\:1\/1?; #Hom alt
			$xrecessive_father = $samples[$sampleid].q?:\S+\:1\/1?; #Het (Male: GT is 1/1 for het on X)
			push(@compound_aff_model,$samples[$sampleid]);
		    }
		}
	    }
	    if ($fatherDS == 0) { #healthy
		print STDOUT "Unaffected Father: $samples[$sampleid]\n";
		if ($_[0] eq "MEDIUM") {
		    $adom_father = $samples[$sampleid].q?:\S+\:0\/0?; #Hom ref
		    $arecessive_father = $samples[$sampleid].q?:\S+\:0\/1?; #Het
		    $xrecessive_father = $samples[$sampleid].q?:\S+\:0\/0?; #Hom ref
		    $denovo_dom_father = $samples[$sampleid].q?:\S+\:0\/0?; #Hom ref
		    $denovo_rec_father = $samples[$sampleid].q?:\S+\:0\/(1|0)?; #Het/Hom ref, AR_denovo-1,2
		    $denovo_x_father = $samples[$sampleid].q?:\S+\:0\/0?; #Hom reference
		    push(@compound_hea_model,$samples[$sampleid]);
		}
	    }
	    push(@adom_model,$adom_father);
	    push(@arecessive_model,$arecessive_father);
	    push(@xrecessive_model,$xrecessive_father);
	    push(@denovo_dom_model,$denovo_dom_father);
	    push(@denovo_rec_model,$denovo_rec_father);
	    push(@denovo_x_model,$denovo_x_father);	    
	}
	else {#Child
	    	
	    for (my $i=0;$i<scalar(@sid_aff);$i++) { #All affected samples	    
		if ($samples[$sampleid] eq $sid_aff[$i] ) {#Affected
		    $childDS = 1;
		    if ($_[0] eq "MEDIUM") {
			$arecessive_child = $samples[$sampleid].q?:\S+\:1\/1?; #Hom alt
			$denovo_dom_child = $samples[$sampleid].q?:\S+\:0\/1?; #Het
			$denovo_rec_child = $samples[$sampleid].q?:\S+\:1\/1?; #Hom alt
			if ( ($fatherDS == 1) && ($motherDS == 1) ) { #Parents affected
			    $adom_child = $samples[$sampleid].q?:\S+\:(0|1)\/1?; #Het/hom alt
			}
			else {
			    $adom_child = $samples[$sampleid].q?:\S+\:0\/1?; #Het, AD-1,2 (Common)
			}
			push(@compound_aff_model,$samples[$sampleid]);
			if (@sid_male) {
			    for (my $males=0;$males<scalar(@sid_male);$males++) { #All male samples
				if ($samples[$sampleid] eq $sid_male[$males]) { #Affected and Male			
				    $xrecessive_child = $samples[$sampleid].q?:\S+\:1\/1?; #Het (Male: GT is 1/1 for het on X)		    
				    $denovo_x_child = $samples[$sampleid].q?:\S+\:1\/1?; #Het (Male: GT is 1/1 for het on X)
				    print STDOUT "Affected male child: $samples[$sampleid]\n";
				}
			    }
			}
			if (@sid_female) {
			    for (my $females=0;$females<scalar(@sid_female);$females++) { #All male samples
				if ($samples[$sampleid] eq $sid_female[$females]) { #Affected and Female
				    $xrecessive_child = $samples[$sampleid].q?:\S+\:(0|1)\/1?; #Het/hom alt
				    $denovo_x_child = $samples[$sampleid].q?:\S+\:(0|1)\/1?; #Het/hom alt (Hom more unlikely) )
				    print STDOUT "Affected female child: $samples[$sampleid]\n";
				}
			    }
			}
		    }
		}
	    }
	    if ($childDS == 0) { #healthy
		
		if ($_[0] eq "MEDIUM") {
		    $adom_child = $samples[$sampleid].q?:\S+\:0\/0?; #Hom ref
		    $arecessive_child = $samples[$sampleid].q?:\S+\:0\/(1|0)?; #Het/Hom ref
		    $denovo_dom_child = $samples[$sampleid].q?:\S+\:0\/0?; #Hom ref
		    $denovo_rec_child = $samples[$sampleid].q?:\S+\:0\/(1|0)?; #Het/Hom ref, AD_denovo-1,2
		    push(@compound_hea_model,$samples[$sampleid]);
		    if (@sid_male) { #Male
			for (my $males=0;$males<scalar(@sid_male);$males++) { #All male samples
			    if ($samples[$sampleid] eq $sid_male[$males]) { #Healthy and Male
				$xrecessive_child = $samples[$sampleid].q?:\S+\:0\/0?; #Hom (Male: GT is 1/1 for het on X)
				$denovo_x_child = $samples[$sampleid].q?:\S+\:0\/0?; #Hom (Male: GT is 1/1 for het on X)
				print STDOUT "Unaffected male child: $samples[$sampleid]\n";
			    }
			}
		    }
		    if (@sid_female) {
			for (my $females=0;$females<scalar(@sid_female);$females++) { #All female samples
			    if ($samples[$sampleid] eq $sid_female[$females]) {#Healthy and Female
				$xrecessive_child = $samples[$sampleid].q?:\S+\:0\/(1|0)?; #Het/Hom ref
				$denovo_x_child = $samples[$sampleid].q?:\S+\:0\/(1|0)?; #Het/Hom ref
				print STDOUT "Unaffected female child: $samples[$sampleid]\n";
			    }
			}
		    }
		}
	    }
	    push(@adom_model,$adom_child);
	    push(@arecessive_model,$arecessive_child); 
	    push(@xrecessive_model,$xrecessive_child);
	    push(@denovo_dom_model,$denovo_dom_child);
	    push(@denovo_rec_model,$denovo_rec_child);
	    push(@denovo_x_model,$denovo_x_child);
	}
    }
}

sub ReadVCF {
#Reads infile (again) and ranks all variants.
#$_[0] = filename
#$_[1] = chr number
#$_[2] = next chr
	
    my $het_C=0; #To track if all GT are identical 
    my $hom_C=0; #To track if all GT are identical
    my $correct_GT=0; #To track how many GTs are correct per variant
    my $correct_AR_model=0; #To track how many GTs are correct per variant and model
    my $thg_frequency_score=0;
    my $dbsnp129_frequency_score=0;
    my $common_frequency_score=0;
    my $pth=0; #Track passed functional prediction
    my $skip=0;
    my %genetic_model_score; #Track genetic model score
###
#Create Weighted scores
###
    my $weighted_rec_genetic_model_score =0;
    my $weighted_classic_genetic_model_score =0;
    my $weighted_genetic_model_score =0;
    my $weighted_location_score =0;
    my $weighted_frequency_score =0;
    my $weighted_funcann_score = 0;
    my $weighted_imdb_score = 0;
    my $weighted_dispgl_score = 0;
    my $weighted_pth_score = 0;
    my $weighted_gtcall_score = 0;
    my $weighted_regcons_score = 0;
    my $weighted_segdup_score = 0;
    my $weighted_basecons_score = 0;
    my $weighted_phylop_score = 0;
    my $weighted_hgmd_score = 0;
    my $weighted_tfsb_score = 0;
    my $weighted_mirna_score= 0;
    my $weighted_score = 0;
    my $weighted_score_zero_tracker = 0;
    
    open(VCF, "<$_[0]") or die "Can't open $_[0]:$!, \n";    
     
    while (<VCF>) {

	chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
            next;
        }	
	if (m/^\#/) {		#Load header and split for use later in writting to outfile
	    @infile_header = split("\t",$_);
            next;
        }			
	if ($_[2] && ($_ =~/^(\S+)/) && ($1 eq $_[2])) { #If next chr is found return (Since all infiles are sorted this is ok)
	    print STDERR "Finished Reading chr$_[1] in Infile $_[0]","\n";
	    close(VCF);
	    return;
	}
	if ( ($_ =~/^(\S+)/) && ($1 eq $_[1]) ) {#Collect chr number and then if it matches for loop go ahead and collect
	    my @temp = split("\t",$_);	    #Loads variant calls
	    $allVariants{$temp[0]}{$temp[1]}{$temp[4]} = [@temp];#saves variant line
	    
	    
###
#Check genetic models and produce genetic model score
###	    
	    
	    for (my $samples=0;$samples<scalar(@adom_model);$samples++) {
		if ( $_=~/$adom_model[$samples]/) {
		    $genetic_model_score{"AD"}++;
		}
	    }
	    for (my $samples=0;$samples<scalar(@arecessive_model);$samples++) {
		if ( $_=~/$arecessive_model[$samples]/) {
		    $genetic_model_score{"AR"}++;
		}
	    }
	    for (my $samples=0;$samples<scalar(@xrecessive_model);$samples++) {
		if ( ($_=~/$xrecessive_model[$samples]/) && ($temp[0]=~/X/) ) {#Only valid on chrX
		    $genetic_model_score{"X"}++;
		}
	    }
	    for (my $samples=0;$samples<scalar(@denovo_dom_model);$samples++) {
		if ( $_=~/$denovo_dom_model[$samples]/) {
		    $genetic_model_score{"AD_denovo"}++;
		}
	    }
	    for (my $samples=0;$samples<scalar(@denovo_rec_model);$samples++) {
		if ( $_=~/$denovo_rec_model[$samples]/) {
		    if ( ($_=~ /\d+-2-1U:\S+\:0\/1/) && ($_=~ /\d+-2-2U:\S+\:0\/1/) ) { #True AR, hence does not have to be included in the AR_denovo model
		    }
		    else {
			$genetic_model_score{"AR_denovo"}++;
		    }
		}
	    }
	    for (my $samples=0;$samples<scalar(@denovo_x_model);$samples++) {
		if ( ($_=~/$denovo_x_model[$samples]/) && ($temp[0]=~/X/) ) { #Only valid on chrX
		    $genetic_model_score{"X_denovo"}++;
		}
	    }
	    
	    if ( ($ar_comp{$temp[0]}{$temp[1]}{$temp[4]}) && $filtered{$ar_comp{$temp[0]}{$temp[1]}{$temp[4]}[$genecol]} ) {#If variant passed AR_compound check and gene passed filtering (countGT and CompareGT)
		#$genetic_model_score{"AR_compound"} = $nos;
		$genetic_model_score{"AR_compound"} = scalar(@samples);
	    }
	    for my $m ( keys( %genetic_model_score ) ) { #Determine weighted genetic model score
		if ( ($genetic_model_score{$m} eq scalar(@samples) ) && ( $rec_genetic_models{$m} ) ) { #Models is correct in all subjects and suggested by MD
	       # if ( ($genetic_model_score{$m} eq $nos) && ( $rec_genetic_models{$m} ) ) { #Models is correct in all subjects and suggested by MD
		    $allVariants{$temp[0]}{$temp[1]}{$temp[4]}[scalar(@temp)] .= $m.";"; #Add passed genetic model to variant line
		    $weighted_rec_genetic_model_score = 3;
		}
		elsif ($genetic_model_score{$m} eq scalar(@samples) ) {
		#elsif ($genetic_model_score{$m} eq $nos) {
		    $allVariants{$temp[0]}{$temp[1]}{$temp[4]}[scalar(@temp)] .= $m.";"; #Add passed genetic model to variant line
		    $weighted_classic_genetic_model_score = 1;
		}
	    }
	    if ( $weighted_rec_genetic_model_score == 3) {
		$weighted_genetic_model_score = 3; #Higher score due to MD input but max 3 points no matter how many models are passed
	    }
	    elsif ($weighted_classic_genetic_model_score == 1) {
		$weighted_genetic_model_score = 1; #Follows classical model
	    }
	    else {
		$weighted_genetic_model_score = -12; #Pretty much removes variants that does not follow a classical model
		$allVariants{$temp[0]}{$temp[1]}{$temp[4]}[scalar(@temp)] .= "Na"; #Add that none genetic model passed to variant line
	    }
###
#Location
###	 

	    if ( ($temp[$loccol]=~/exon/) || ($temp[$loccol]=~/splicing/) ) {
		$weighted_location_score = 3;
	    }
	    elsif ( ($temp[$loccol]=~/intron/) || ($temp[$loccol]=~/UTR3/) || ($temp[$loccol]=~/UTR5/) ) {
		$weighted_location_score = 1;
	    }
	    elsif ( ($temp[$loccol]=~/upstream/) || ($temp[$loccol]=~/downstream/) ) {
		$weighted_location_score = 1;
	    }
	    
###
#Frequency
###
	    if ($temp[$thGcol] eq "-") { #Not found in 1000G
		$thg_frequency_score = 3;
	    }
	    elsif ( $temp[$thGcol]=~/(\d+.\d+|\d+)/ ) {
		unless ($1) { 
		    $thg_frequency_score = 3;
		}
		if ($1) {
		    
		    if ($1<=0.005 ) { 
			$thg_frequency_score = 2;
		    } 
		    elsif ($1<=0.02 ) { # Less stringent
			$thg_frequency_score = 1; 
		    } 
		    elsif ($1>0.02 ) { # Common variant
			$common_frequency_score=1;	 
		    } 
		}
	    }
	    if ($temp[$dbsnp129col] eq "-") {
		$dbsnp129_frequency_score = 3;
	    }
	    elsif ( $temp[$dbsnp129col]=~/(\d+.\d+|\d+)/ ) { 

		unless ($1) { 
		    $dbsnp129_frequency_score = 3;
		}
		
		if ($1) {
		    
		    if ($1<=0.005 ) { 
			$dbsnp129_frequency_score = 2;
		    } 
		    elsif ($1<=0.02 ) { # Less stringent
			$dbsnp129_frequency_score = 1;
		    } 
		    elsif ($1>0.02 ) { # Common variant
			
			$common_frequency_score=1;	 
		    } 
		}
	    }
	    if ($common_frequency_score == 1) {
		$weighted_frequency_score =-12;
	    }
	    elsif ( ($thg_frequency_score + $dbsnp129_frequency_score) == 6 ){ #Very Rare
		$weighted_frequency_score =3; 
	    }
	    elsif (($thg_frequency_score + $dbsnp129_frequency_score) >=4) { #Rare
		$weighted_frequency_score =2; 
	    }
	    elsif (($thg_frequency_score + $dbsnp129_frequency_score) >=2) { 
		$weighted_frequency_score =1;
	    }
	    else {
		$weighted_frequency_score = -1;
	    }
	    #dbSNPNonFlagged
	    if ( $temp[$dbsnpcol] ne $annovar_dbsnp_ver) { #Not found in snp135NonFlagged
		$weighted_frequency_score = $weighted_frequency_score + 1;
	    }
###
#Functional Annotation
###
	    if ( ($temp[$syncol]=~/frameshift/) || ($temp[$syncol]=~/stop/) || ($temp[$syncol]=~/nonsynonymous/) ) { #Detects frameshift insertion, frameshift deletion, frameshift block substitution, stopgain, stoploss and  nonsynonymous SNV
		$weighted_funcann_score = 3;
	    }
	    elsif ( $temp[$syncol]=~/nonframeshift/) { #Detects nonframeshift insertion, nonframeshift deletion, nonframeshift block substitution
		$weighted_funcann_score = 3;
	    }
	    elsif ( $temp[$syncol]=~/unknown/) { #Detects unknown
		$weighted_funcann_score = 1;
	    }
	    else {
		$weighted_location_score = -2;
	    }
###
#IM_DB
###
	    my $gene_hit_counter =0; #Required for non-overlapping features to imdb i.e. coordinates overlap but the geneName is not interesting according to the imdb.  
	    my @temp_ensemble_Genes_Id = split(";", $temp[$ensemblegeneidcol]); #Variant can have multiple entries due to overlapping genes
 	    for (my $variant_ensemble_id_count=0;$variant_ensemble_id_count<scalar(@temp_ensemble_Genes_Id);$variant_ensemble_id_count++) { #All variants
		#print $temp_ensemble_Genes_Id[$variant_ensemble_id_count], "\n";
		if ($imdb{$temp_ensemble_Genes_Id[$variant_ensemble_id_count]} ) { #If geneName exists in imdb (only single entry)
		    $gene_hit_counter++;
		    if ($im_db) {
			$weighted_imdb_score = 3;
		    }
		    
		    my $sampleId_counter =0;
		    #print $_, "\n";     
		    for my $sampleID (keys %imdb_tarcov)  { #For all sampleIDs
			my $feature_counter =0; #Required for total sum of passed features per sampleID and Gene
			my $feature_hit_counter =0; #Required for sum of passed features per sampleID and Gene
				#print $_, "\n";
			if ($imdb_tarcov{$sampleID}{$temp_ensemble_Genes_Id[$variant_ensemble_id_count]}) { #If EnsembleGene id exists	
			    
			    if ($imdb_tarcov{$sampleID}{$temp_ensemble_Genes_Id[$variant_ensemble_id_count]}{$temp[0]}) { #Check for presence of chr in imdb_tarcov. If no poorly coverade region exist on chr then pass sampleID otherwise check overlapp
				my $chr = $temp[0]; #Collect for chr
				my $poor_coverage_counter=0; #If a feature is found within the gene range.
				
				for my $startpos (keys % { $imdb_tarcov{$sampleID}{$temp_ensemble_Genes_Id[$variant_ensemble_id_count]}{$chr} })  { #For all start    
				    
				    for my $endpos (keys % { $imdb_tarcov{$sampleID}{$temp_ensemble_Genes_Id[$variant_ensemble_id_count]}{$chr}{$startpos} })  { #For end position(s) (should only be 1)
					
					if ($startpos >= $imdb_pos{$temp_ensemble_Genes_Id[$variant_ensemble_id_count]}[1] && $startpos <= $imdb_pos{$temp_ensemble_Genes_Id[$variant_ensemble_id_count]}[2]) { #Target start is downstream of gene start and upstream of gene end i.e. within gene range therefore overlapping 
					    $feature_counter++;
					    #print "Entry\t $sampleID\t $temp[0]\t $temp[1]\t $temp[2]\t $temp_ensemble_Genes_Id[$variant_ensemble_id_count]\n";
					    #print "Feature\t $chr\t $startpos\t $endpos\t $temp_ensemble_Genes_Id[$variant_ensemble_id_count]\n";
					    #print "Gene\t $imdb_pos{$temp_ensemble_Genes_Id[$variant_ensemble_id_count]}[0]\t $imdb_pos{$temp_ensemble_Genes_Id[$variant_ensemble_id_count]}[1]\t $imdb_pos{$temp_ensemble_Genes_Id[$variant_ensemble_id_count]}[2]\n";
					    #print "$sampleID\t $chr\t $startpos\t $endpos\t $sampleId_counter\n";
					    if ($imdb_tarcov{$sampleID}{$temp_ensemble_Genes_Id[$variant_ensemble_id_count]}{$chr}{$startpos}{$endpos} eq "PASS") { #Feature overlaps gene but has adequate coverage
						$feature_hit_counter++;
					    }
					    else {
						$poor_coverage_counter++;
						#if ($cmms_imdb == 1) {
						 #   $allVariants{$temp[0]}{$temp[1]}{$temp[4]}[scalar(@temp)+3+$sampleId_counter] .= $imdb_tarcov{$sampleID}{$temp_ensemble_Genes_Id[$variant_ensemble_id_count]}{$chr}{$startpos}{$endpos};	
						#}
						#else {
						    $allVariants{$temp[0]}{$temp[1]}{$temp[4]}[scalar(@temp)+2+$sampleId_counter] .= $imdb_tarcov{$sampleID}{$temp_ensemble_Genes_Id[$variant_ensemble_id_count]}{$chr}{$startpos}{$endpos};	
						#}
						#print  $imdb_tarcov{$sampleID}{$chr}{$startpos}{$endpos}, "\n";
					    }
					    #last;
					}
					elsif ($startpos <= $imdb_pos{$temp_ensemble_Genes_Id[$variant_ensemble_id_count]}[1] && $endpos >= $imdb_pos{$temp_ensemble_Genes_Id[$variant_ensemble_id_count]}[1]) { #Target start is upstream of gene start and target end is on gene i.e. within gene range therefore overlapping
					    $feature_counter++;
					    #print "Entry found for: $sampleID\t $temp[0]\t $temp[1]\t $temp[2]\n";
					    #print "$sampleID\t $chr\t $startpos\t $endpos\n";
					    if ($imdb_tarcov{$sampleID}{$temp_ensemble_Genes_Id[$variant_ensemble_id_count]}{$chr}{$startpos}{$endpos} eq "PASS") { #Feature overlaps gene but has adequate coverage
						$feature_hit_counter++;
					    }
					    else {
						$poor_coverage_counter++;
						#if ($cmms_imdb == 1) {
						 #   $allVariants{$temp[0]}{$temp[1]}{$temp[4]}[scalar(@temp)+3+$sampleId_counter] .= $imdb_tarcov{$sampleID}{$temp_ensemble_Genes_Id[$variant_ensemble_id_count]}{$chr}{$startpos}{$endpos};
						#}
						#else {
						    $allVariants{$temp[0]}{$temp[1]}{$temp[4]}[scalar(@temp)+2+$sampleId_counter] .= $imdb_tarcov{$sampleID}{$temp_ensemble_Genes_Id[$variant_ensemble_id_count]}{$chr}{$startpos}{$endpos};
						#}
						#print  $imdb_tarcov{$sampleID}{$chr}{$startpos}{$endpos}, "\n";		
					    }
					}
				    }
				}
				if ($poor_coverage_counter == 0) { #No features with poor coverage for this sampleID
				    #print "No entry found for: $sampleID\t $temp[0]\t $temp[1]\t $temp[2]\n";
				    #if ($cmms_imdb == 1) {
				#	$allVariants{$temp[0]}{$temp[1]}{$temp[4]}[scalar(@temp)+3+$sampleId_counter] = "$sampleID:PASS";	
				 #   }
				  #  else {
					$allVariants{$temp[0]}{$temp[1]}{$temp[4]}[scalar(@temp)+2+$sampleId_counter] .= "$sampleID:PASS";
				   # }
				}
			    } #No regions with poor coverage on entire chr (rather unlikely, but checked for)
			    else {
				#if ($cmms_imdb == 1) {
				#    $allVariants{$temp[0]}{$temp[1]}{$temp[4]}[scalar(@temp)+3+$sampleId_counter] = "$sampleID:PASS";	
				#}
				#else {
				    
				    $allVariants{$temp[0]}{$temp[1]}{$temp[4]}[scalar(@temp)+2+$sampleId_counter] .= "$sampleID:PASS";
				#}
			    }
			}
			#if ($cmms_imdb == 1) {
			#    $allVariants{$temp[0]}{$temp[1]}{$temp[4]}[scalar(@temp)+3+$sampleId_counter] .= ";$feature_hit_counter/$feature_counter, ";
			#}
			#else {
			    $allVariants{$temp[0]}{$temp[1]}{$temp[4]}[scalar(@temp)+2+$sampleId_counter] .= ";$feature_hit_counter/$feature_counter, ";
			#}
			$sampleId_counter++;
		    }
		}
	    }
	    if ( $gene_hit_counter == 0 ) { #No feature that has the correct geneName according to imdb, hence not interesting
		
		#if ($cmms_imdb == 1) {
		 #   for (my $sampleID=0;$sampleID<scalar(@samples);$sampleID++) { # Add for all sampleIDs
		#	$allVariants{$temp[0]}{$temp[1]}{$temp[4]}[scalar(@temp)+2] = "#Na";
		#	$allVariants{$temp[0]}{$temp[1]}{$temp[4]}[scalar(@temp)+3+$sampleID] = "$samples[$sampleID]:No overlapping feature found";		
		 #   }
		#}
		#else {
		    for (my $sampleID=0;$sampleID<scalar(@samples);$sampleID++) { # Add for all sampleIDs
			$allVariants{$temp[0]}{$temp[1]}{$temp[4]}[scalar(@temp)+2+$sampleID] .= "$samples[$sampleID]:No overlapping feature found";		
		 #   }
		}
	    } 
	    
###
#IEM_disp_list
###
	    if ($dispgenelist{$temp[$genecol]} ) {
		$weighted_dispgl_score = -30;
	    }
###
#Functional Prediction
###	    
	    if ($temp[$avsiftcol] ne "-" || $temp[$pp2col] ne "-" || $temp[$muttascol] ne "-" ) {
		if( ($temp[$avsiftcol]=~/(\d+.\d+\S+|\d+\S+)/ ) || ($temp[$avsiftcol]=~/(\d+.\d+|\d+)/ ) ) {
		    if ($1<=0.05) {
			$pth++;
		    } 
		}
		if( ($temp[$muttascol]=~/(\d+.\d+\S+|\d+\S+)/ ) || ($temp[$muttascol]=~/(\d+.\d+|\d+)/ ) ) {
		    if ($1>=0.5) { 
			$pth++;
		    } 
		} 
		if( ($temp[$pp2col]=~/(\d+.\d+\S+|\d+\S+)/ ) || ($temp[$pp2col]=~/(\d+.\d+|\d+)/ ) ) {
		    if ($1>=0.85) { 
			$pth++;
		    } 
		}
		$weighted_pth_score = $pth; #Score correlates with how many predictors predicts the variant to be damaging.
	    } 
###
#GT CALL
###
	    if ($_=~/\:PRES\:/) {
		$weighted_gtcall_score = 1;
	    }
	    elsif ($_=~/\:PASS\:/) {
		$weighted_gtcall_score = 3;
	    }

###
#Region Conservation
###	
	    if ($temp[$mce46waycol]=~/Score/ && $temp[$gerpelemcol]=~/Name/) {
		$weighted_regcons_score = 2;
	    }
	    elsif ($temp[$mce46waycol]=~/score/ || $temp[$gerpelemcol]=~/Name/) {
		$weighted_regcons_score = 1;
	    }
###
#Segmental Duplication
###	
	    if ($temp[$segdupcol]=~/Score/) {
		$weighted_segdup_score = -2;
	    }
###
#Base conservation
###
	    if ( ($temp[$gerpbasecol]=~/(\d+.\d+\S+|\d+\S+)/) || ($temp[$gerpbasecol]=~/(\d+.\d+|\d+)/ ) ) {
		if ($1>=4) {
		    $weighted_basecons_score = 2;
		} 
		elsif ($1>=2) {
		    $weighted_basecons_score = 1;
		} 
	    }
	    if ( ($temp[$phylopcol]=~/(\d+.\d+\S+|\d+\S+)/) || ($temp[$phylopcol]=~/(\d+.\d+|\d+)/ ) ) {
		if ($1>=0.9984188612) { #Original phyloP score > 2.5. Set according to Disease gene identification strategies for exome sequencing, European Journal of Human Genetics (2012) 20, 490–497
		    $weighted_phylop_score = 2;
		} 
		elsif ($1>=0.95) { #dbNSFP: A Lightweight Database of Human Nonsynonymous SNPs and Their Functional Predictions - Boerwinkle. (DOI 10.1002/humu.21517)
		    $weighted_phylop_score = 1;
		} 
	    }
###
#HGMD
###	
	    if ($temp[$hgmdcol]=~/\d;/) {
		$weighted_hgmd_score = 1;
	    }
###
#Transcription factor binding
###	
	    #if ($_=~/tfbs\;/) {
	#	$weighted_tfsb_score = 0;
	 #   }
###
#Sno & miRNA annotations
###	
	  #  if ($_=~/mirna\;/) {
	#	$weighted_mirna_score = 0;
	 #   }
	    $weighted_score = $weighted_genetic_model_score + $weighted_location_score + $weighted_frequency_score + $weighted_funcann_score + $weighted_imdb_score + $weighted_dispgl_score + $weighted_pth_score + $weighted_gtcall_score + $weighted_regcons_score + $weighted_segdup_score + $weighted_basecons_score + $weighted_phylop_score + $weighted_hgmd_score + $weighted_tfsb_score + $weighted_mirna_score;
	    $allVariants{$temp[0]}{$temp[1]}{$temp[4]}[scalar(@temp)+1] = $weighted_score; #Add weighted score to master list
	    
	    if ( ($weighted_score == 0) || ($weigthed_scores_seen{$weighted_score}) ) {
		if ( $weighted_score_zero_tracker == 0 && $weighted_score == 0) {
		    push(@weigthed_scores_seen, $weighted_score);
		    $weighted_score_zero_tracker++;
		}
	    }
	    else {
		$weigthed_scores_seen{$weighted_score} = $weighted_score;
		push(@weigthed_scores_seen, $weighted_score);
	    }
	    if($ar_comp{$temp[0]}{$temp[1]}{$temp[4]}) { #If variant passed sub AnalyseCompound
		#print $_,"\n", $ar_comp{$temp[0]}{$temp[1]}{$temp[4]}[$genecol],"\n";
		for (my $spg=0;$spg<scalar( @{ $ar_comp_gene{$temp[0]}{$temp[1]}{$temp[4]} });$spg++) { #Check all geneNames for variant position
		
		    if ( $filtered{$ar_comp_gene{$temp[0]}{$temp[1]}{$temp[4]}[$spg]} ) {#If variant passed AR_compound_gene check and gened passed filtering. NOTE: Since only chr,start,alt allel is checked in %ar_comp - variants with complex annotations (X;y,N() etc) will not be missed. In filtered these issues have been solved and the gene names are HGNC, since the input which leads to the %filtered is %genes. 
			#print $_,"\t", $ar_comp{$temp[0]}{$temp[1]}{$temp[4]}[$genecol],"\t", $ar_comp_gene{$temp[0]}{$temp[1]}{$temp[4]}[$spg], "\t", "$weighted_score\n";
			if ( $weighted_score >= $rankscore ) { #Only want to write AR_compund variants later where at least 2 variants have decent score
			    $filtered_ar_comp_weigthed_score{$ar_comp_gene{$temp[0]}{$temp[1]}{$temp[4]}[$spg]}++; #NOTE: Counts weighted_score >= rankscore using parsed GeneNames.
			    #print $temp[$genecol], "\n";
			}
		    }
		}
	    }
	    #if ($_=~/16932490/) {
	#	print "Weighted Score: $weighted_score\t Genetic_Model_score: $weighted_genetic_model_score\t Location_score: $weighted_location_score\t Frequency_score: $weighted_frequency_score\t Funcann: $weighted_funcann_score\t Im_Db_CMMS: $weighted_imdb_score\t  IEM_dispgl: $weighted_dispgl_score\t Prediction: $weighted_pth_score\t GT Call: $weighted_gtcall_score\t Region Conservation: $weighted_regcons_score\t Segemental Dup: $weighted_segdup_score\t GerpBase: $weighted_basecons_score\t PhyloP: $weighted_phylop_score\t HGMD: $weighted_hgmd_score\t TSFB: $weighted_tfsb_score\t miRNA: $weighted_mirna_score\n";
	#	print $_,"\n";
	 #  }
	}
	%genetic_model_score = (); #Reset for next line
	$thg_frequency_score = 0; #Reset for next line
	$dbsnp129_frequency_score = 0; #Reset for next line
	$common_frequency_score = 0; #Reset for next line
	$pth = 0; #Reset for next line
	$weighted_rec_genetic_model_score = 0; #Reset for next line
	$weighted_classic_genetic_model_score = 0; #Reset for next line
	$weighted_genetic_model_score = 0; #Reset for next line
	$weighted_location_score = 0; #Reset for next line
	$weighted_frequency_score = 0; #Reset for next line
	$weighted_funcann_score = 0;  #Reset for next line
	$weighted_imdb_score = 0; #Reset for next line
	$weighted_dispgl_score = 0; #Reset for next line
	$weighted_pth_score = 0; #Reset for next line
	$weighted_gtcall_score = 0; #Reset for next line
	$weighted_regcons_score = 0; #Reset for next line
	$weighted_segdup_score = 0; #Reset for next line
	$weighted_basecons_score =0; #Reset for next line
	$weighted_phylop_score =0; #Reset for next line
	$weighted_hgmd_score = 0; #Reset for next line
	$weighted_tfsb_score = 0; #Reset for next line
	$weighted_mirna_score = 0; #Reset for next line
    } 	
    close(VCF);
    print STDOUT "Finished Reading and scoring all $_[1] variants in Infile: $_[0]","\n";
    return;
}

sub SortAllVariants {
#Creates an array of all position which are unique and in sorted ascending order
    
    for my $chr (keys %allVariants)  { #For all chr
	
	for my $pos (keys %{ $allVariants{$chr} } )  { #For all pos
	    
	    for my $variant (keys % { $allVariants{$chr}{$pos} })  { #For all variants
		push ( @{$allVariants_chr{$chr} },$pos );
	    }
	}
	my %seen = (); @{$allVariants_chr_unique{$chr} } = grep { ! $seen{$_} ++ } @{$allVariants_chr{$chr} }; #Unique entries only
	@{$allVariants_chr_sorted{$chr} } = sort { $a <=> $b } @{ $allVariants_chr_unique{$chr} }; #Sorts keys to be able to print sorted table later
    }
#Sort weighted score
    @weighted_scores_sorted = sort { $b <=> $a } @weigthed_scores_seen;
    #for (my $i=0;$i<scalar(@weighted_scores_sorted);$i++) {
	#print $weighted_scores_sorted[$i], "\n";
    #}
    #print STDOUT "Sorted all non overlapping entries per chr and position\n";
}


sub WriteCompound {
#Write filtered variants
#$_[0] = filename
#$_[1] = chr number
    
    if ($_[1] eq "chr1") {
	open (FILT, ">$_[0]") or die "Can't write to $_[0]: $!\n";
#Create header
	if ( @infile_header ) { #If the file came with a header otherwise to not write a header

	    for (my $header_elementsCounter=0;$header_elementsCounter<scalar(@infile_header);$header_elementsCounter++) { #header elements
		
		print FILT $infile_header[$header_elementsCounter], "\t"; #Recreate old header
	    }
	    print FILT "GeneModel\tRank_Score\t"; #Add new info 
	    
	    if ($im_db_cc ==1) {
		for (my $i=0;$i<$nos;$i++) {
		    print FILT "IDN:PASS;TotalFeaturesPASS||IDN:Chr:Start:Stop:Fraction_ME_Ten-Coverage_Bases:;TotalFeaturesPASS\t"; #Add new info
		}
	    }
	    if ($cmms_imdb ==1) {
		print FILT "Gene_ID(NCBI)\t"."Complete_Gene_name(NCBI)\t"."Additional_names\t"."Genome_build\t"."Isoforms\t"."Transcripts\t"."Function\t"."OMIM\t"."Expression\t"."Animal_models(knockouts)\t"."Reference_PMID\t"."Top_Prio\t";
	    }
	    print FILT "\n";
	}
    }
    else {
	open (FILT, ">>$_[0]") or die "Can't write to $_[0]: $!\n";
    }
    my $ar_comp_variant_tracker = 0;
    my $gene_not_ar_comp_tracker = 0;
    
    for (my $top_rank=0;$top_rank<scalar(@weighted_scores_sorted);$top_rank++) { #Enables top-ranked writing
	#print STDOUT $weighted_scores_sorted[$top_rank], "\n";
	for my $chr (keys %allVariants_chr_sorted) {
	    
	    for (my $i=0;$i<scalar( @{$allVariants_chr_sorted{$chr} } );$i++)  { #For all pos per chr	
		
		my $pos = $allVariants_chr_sorted{$chr}[$i]; #pos keys to hash from sorted arrray
		
		for my $variant (keys % { $allVariants{$chr}{$pos} })  { #For all variants

		    #my $weighted_score = $allVariants{$chr}{$pos}{$variant}[scalar( @{ $allVariants{$chr}{$pos}{$variant} } ) -1 ];
		    my $weighted_score = $allVariants{$chr}{$pos}{$variant}[$weigthedsumcol];

		    #print $allVariants{$chr}{$pos}{$variant}[0], "\t", $allVariants{$chr}{$pos}{$variant}[1], "\n";
		    #print "$weighted_scores_sorted[$top_rank]\t $weighted_score\n"; 
		    if ( $weighted_score >= $rankscore && ($weighted_score == $weighted_scores_sorted[$top_rank]) ) {		    
#if ( ($allVariants{$chr}{$pos}{$variant}[scalar( @{ $allVariants{$chr}{$pos}{$variant} } ) -1 ] >=10) && )  { #Print > than 10	
			
			if ( $ar_comp{$chr}{$pos}{$variant} ) {#If variant passed AR_compound check. (Can still be AR_denovo)
			    #print $ar_comp{$chr}{$pos}{$variant}[$genecol], "\n";
			    
			    for (my $spg=0;$spg<scalar( @{ $ar_comp_gene{$chr}{$pos}{$variant} });$spg++) { #Check all geneNames for variant position
				
				if ( $filtered{$ar_comp_gene{$chr}{$pos}{$variant}[$spg]} ) {#If variant passed AR_compound_gene check and gened passed filtering. NOTE: Since only chr,start,alt allel is checked in %ar_comp - variants with complex annotations (X;y,N() etc) will not be missed. In filtered these issues have been solved and the gene names are HGNC, since the input which leads to the %filtered is %genes.
				    #if ( ($filtered_ar_comp_weigthed_score{$ar_comp{$chr}{$pos}{$variant}[$genecol]}) && ($filtered_ar_comp_weigthed_score{$ar_comp{$chr}{$pos}{$variant}[$genecol]} >=2) ) { #Only want to write AR_compund variants later where at least 2 variants have decent score
			      
				    if ( ($filtered_ar_comp_weigthed_score{$ar_comp_gene{$chr}{$pos}{$variant}[$spg]}) && ($filtered_ar_comp_weigthed_score{$ar_comp_gene{$chr}{$pos}{$variant}[$spg]} >=2) ) { #Only want to write AR_compund variants later where at least 2 variants have decent score.
					if ($ar_comp_variant_tracker ==0) { #Avoid printing variant more than once if it matches more than 1 gene
					    for (my $variants=0;$variants<scalar( @{ $allVariants{$chr}{$pos}{$variant} } );$variants++)  {
						print FILT $allVariants{$chr}{$pos}{$variant}[$variants], "\t";
					    }
					    print FILT "\n";
					    $ar_comp_variant_tracker++;
					}
				    }
				}
				else {
				    $gene_not_ar_comp_tracker++;
				}
			    }
			    if ( ($ar_comp_variant_tracker == 0) && ($allVariants{$chr}{$pos}{$variant}[$genmodelcol] ne "Na" ) )  {#Not previosly printed and match to a Genetic Model other than AR_compound , hence deserves to be printed 
				unless ($allVariants{$chr}{$pos}{$variant}[$genmodelcol] eq "AR_compound;") {#Do not print if AR_compound since the >=2 over 10 condition is not satisfied
				    $allVariants{$chr}{$pos}{$variant}[$genmodelcol] =~ s/AR_compound;//;#Not a relevant AR_compound, rename the genetic model fit to other than AR_compound
				    #print $allVariants{$chr}{$pos}{$variant}[scalar( @{ $allVariants{$chr}{$pos}{$variant} } ) -2 ], "\n";
				    for (my $variants=0;$variants<scalar( @{ $allVariants{$chr}{$pos}{$variant} } );$variants++)  {
					print FILT $allVariants{$chr}{$pos}{$variant}[$variants], "\t";
				    }
				    print FILT "\n";    
				}
			    }
			    elsif ($gene_not_ar_comp_tracker eq scalar( @{ $ar_comp_gene{$chr}{$pos}{$variant} }) ) {#Passed the AR_compound check but not the filtered_gene (Can be AR_denovo)	
				for (my $variants=0;$variants<scalar( @{ $allVariants{$chr}{$pos}{$variant} } );$variants++)  {
				    print FILT $allVariants{$chr}{$pos}{$variant}[$variants], "\t";
				}
				print FILT "\n";
			    }   
			}
			else { #Not a potential AR_compound variant go ahead and print
			    for (my $variants=0;$variants<scalar( @{ $allVariants{$chr}{$pos}{$variant} } );$variants++)  {
				print FILT $allVariants{$chr}{$pos}{$variant}[$variants], "\t";
			    }
			    print FILT "\n";
			}
		    }
		    $ar_comp_variant_tracker = 0; #Reset for next variant
		    $gene_not_ar_comp_tracker = 0; #Reset for next variant
		}
	    }	
	}
    }
    close(FILT);
    print STDOUT "Wrote ranked variants to outfile: $_[0]","\n";
    return;
}

sub ReadForFinalSort {
#Reads just created gene list to enable sorting of these genes for final output.
#$_[0] = filename (Whole path).
#$_[1] = Final out filename (Whole path). 
    
    open(RFFS, "<$_[0]") or die "Can't open $_[0]:$!, \n";    
    my $header;
    my %allVariants;
    my $weighted_score_zero_tracker=0;
    while (<RFFS>) {
        chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
            next;
        }	
	if (m/^\#/) {# Avoid #
	    $header = $_;
            next;
        }				
	if (/(\S+)/) {  
	    my @temp = split("\t",$_);	    #Loads variant calls
	    
	    $allVariants{$temp[0]}{$temp[1]}{$temp[4]} = [@temp];#saves variant line
	    if ( ($temp[$weigthedsumcol] == 0) || ($weigthed_scores_seen{$temp[$weigthedsumcol]}) ) {
		if ( $weighted_score_zero_tracker == 0 && $temp[$weigthedsumcol] == 0) {
		    push(@weigthed_scores_seen, $temp[$weigthedsumcol]);
		    $weighted_score_zero_tracker++;
		}
	    }
	    else {
		$weigthed_scores_seen{$temp[$weigthedsumcol]} = $temp[$weigthedsumcol];
		push(@weigthed_scores_seen, $temp[$weigthedsumcol]); #Saves unique weigthed score
	    }
	}
    } 	
    close(RFFS);
    print STDOUT "Read file: $_[0] for sorting","\n";
    @weighted_scores_sorted = sort { $b <=> $a } @weigthed_scores_seen;
    
    open(WS, ">$_[1]") or die "Can't open $_[1]:$!, \n";
    print WS $header, "\n";
    for (my $top_rank=0;$top_rank<scalar(@weighted_scores_sorted);$top_rank++) { #Enables top-ranked writing
	for my $chr (keys %allVariants) {
	
	    for my $pos (keys % { $allVariants{$chr} })  { #For all positions

		for my $variant (keys % { $allVariants{$chr}{$pos} })  { #For all variants

		    if ($allVariants{$chr}{$pos}{$variant}[$weigthedsumcol] == $weighted_scores_sorted[$top_rank]) {
			for (my $variants=0;$variants<scalar( @{ $allVariants{$chr}{$pos}{$variant} } );$variants++)  {
			    print WS $allVariants{$chr}{$pos}{$variant}[$variants], "\t";  
			}
			print WS "\n"; 
		    }
		}
	    }
	}
    }
    close(WS);
    print STDOUT "Wrote final sorted file: $_[1]","\n";
    return;
}


