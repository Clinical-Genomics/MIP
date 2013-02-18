#!/usr/bin/perl -w

use strict;
use warnings;

#Master scripts for analysing whole genome or exome paired end reads from Illumina pipleine using samtools, GATK and custom perl scripts. The program performs variation calls on bam files and generates output in the vcf format. These files are then annotaded and filtered generating tab-separated candidate list files with meta-data having one variant per line.
#Copyright 2011 Henrik Stranneheim

=head1 SYNOPSIS
    
mip_variation.pl  -i [infile...n] -a [project ID] -s [sample ID...n] -em [e-mail] -ods [outdirscripts] -odf [outdirdata] -ids [indirscripts] -rd [referencedir]

=head2 COMMANDS AND OPTIONS

-i/--infile Infile(s) (Mandatory: Supply whole path)

-ids/--indirscript The pipeline script in dir (Mandatory: Supply whole path)

-rd/--referencesdir Reference(s) dir (Mandatory: Supply whole path)

-gref/--genomeref Flag for setting genomic reference file (defaults to "Homo_sapiens.GRCh37.57.dna.concat.fa")

-a/--projectid The project ID (Mandatory)

-s/--sampleid The sample ID(s) (Mandatory)

-em/--email

-odf/--outdirdata The data files output directory (Mandatory: Supply whole path)

-ods/--outdirscript The script files output directory (Mandatory: Supply whole path)

-familyid/--family Group id of samples to be compared (Mandatory, Ex: 1 for IDN 1-1-1A)

-pSTV_schr/--samtools_viewschr Flag running samtools view to split per chr & index (defaults to "1" (=yes))

-pGATK_REAl/--gatk_real Flag running GATK realign (defaults to "1" (=yes))

-gatk_real_knset1/--gatk_real_knownset1 GATK realign known INDEL set 1 (defaults to "1000G_phase1.indels.hg19.vcf")

-gatk_real_knset2/--gatk_real_knownset2 GATK realign known INDEL set 2 (defaults to "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf")

-gatk_path/--genatk_path  Flag for path to GATK, must be supplied for GATK (defaults to "").  

-pGATK_RECAl/--gatk_recal Flag running GATK recalibrate (defaults to "1" (=yes))

-gatk_recal_knset/--gatk_recal_knownset GATK recal known SNP set (defaults to "dbsnp_135.b37.vcf")

-pGATK_UNIGT/--gatk_unigt Flag running GATK unifiedgenotyper (defaults to "0" (=no))

-pGATK_HAPCAL/--gatk_hapcal Flag running GATK HaplotypeCaller (defaults to "1" (=yes))

-gatkugt_snp/--gatk_unigt_snp Flag running GATK unifiedgenotyper for SNPs (defaults to "1" (=yes))

-gatkugt_ind/--gatk_unigt_indel Flag running GATK unifiedgenotyper for INDELs (defaults to "1" (=yes))

-gatk_bait/--gatk_bait_il Prepared bait interval_list file for GATK_UnifiedGT. (defaults to "SureSelect_All_Exon_50mb_with_annotation_hg19_nochr.bed.pad100.interval_list")

-pGATK_VARRECAL/--gatk_varrecalibrator Flag running GATK variantrecalibrator  (defaults to "1" (=yes))

-gatk_exref_snp/--gatk_exomeref_snp Prepared exome reference file (SNVs) for GATK_Varrecal. (defaults to "all-agilent_50mb-GRCh37-SNPS_pad100_interval_list.vcf")

-gatk_exref_indel/--gatk_exomeref_indel Prepared exome reference file (INDELs) for GATK_Varrecal. (defaults to "all-agilent_50mb-GRCh37-INDELS_pad100_interval_list.vcf")

-pPIND/--pindel Flag running Pindel (defaults to "0" (=no))

-pGATK_VAREVAL_All/--gatk_vareval_all Flag running GATK varianteval for all variants  (defaults to "1" (=yes))

-pGATK_VAREVAL_Exome/--gatk_vareval_exome Flag running GATK varianteval for exonic variants  (defaults to "1" (=yes))

-pGATK_COMBVAR/--gatk_combinevariants Flag running GATK combinevariants. Use only if UnifiedGT and not HaplotypeCaller was used to call variants (SNVs & INDELs).(defaults to "0" (=no))

-pANVAR/--annovar Flag running annovar (defaults to "1" (=yes))

-anva_pa/--annovar_path  Path to annovar script directory (Supply whole path, defaults to "". NOTE:Assumes that the annovar db files are located in annovar/humandb)

-anva_gbv/--annovar_genome_build_version Annovar genome build version (defaults to "hg19")

-anva_tn/--annovar_table_names Annovar table names (comma sep)

-anva_stn/--annovar_supported_table_names Annovar MIP supported table names (defaults 0 (=no))

-anva_maf_th/--annovar_maf_threshold Flag for setting the minor allele frequency threshold in annovar (defaults to "0" (=no))

-anva_sift_th/--annovar_sift_threshold Flag for setting the avsift threshold in annovar (defaults to "0" (=no))

-pVMERGE/--variation_annotation_merge Running intersectCollect.pl to merge all annotation info to one file (defaults to "1" (=yes))

-vm_db_te/--vmerge_db_template Db template file used to create the specific family '-vm_dbf' master file (defaults to "CMMS_intersectCollect_db_master_template.txt")

-vm_dbf/--vmerge_db_file Db master file to be used in intersectCollect.pl (defaults to "FDN.intersectCollect_db_master.txt")

-pAddDP/--adddepth Flag for adding depth at nonvariant sites by mpileup and add_depth.pl (defaults to "1" (=yes))

-pRankVar/--rankvariants Flag running ranking of variants (defaults to "1" (=yes))

-rs/--rankscore The rank score cut-off (defaults to "-100")

-dgf/--dispGeneFiltering Filtering of genes that should be removed from downstream processing (Defaults to "1" (=yes))

-dgfl/--dispGeneList List of genes that should be removed from downstream processing (Supply whole path, Format: 1 entry per line;HGNC Symbol)

-all_db_file/--all_elements_Db_file All_Db file (Defaults to "mart_export_Ensembl_GeneID_key_cleaned_chr.txt")

-all_db_cc/--all_elements_Db_Gene_Coverage_Calculation All_Db file coverage calculation (Defaults to "1" (=yes))

-all_db_gidc/--all_elements_Db_Gene_Id_Col All_Db file gene Id column (zero-based, defaults to "4")

-im_db_file/--Im_Db_file Im_Db file (Defaults to "IEM_Db_CMMS_version1.2.txt")

-im_db_te/--Im_Db_template Im_Db template file used to create the specific family '-im_dbmf' master file (Defaults to "select_dbIEM_variants_db_master.txt")

-im_dbmf/--Im_db_master_file Db master file to be used when selecting variants (defaults to "FDN.intersectCollect_selectVariants_db_master.txt")

-im_dbfof/--Im_db_file_out_file The file(s) to write to when selecting variants with intersectCollect.pl. Comma sep (defaults to "$odf/$familyid/$aligner/GATK/candidates/ranking/$familyid_orphan.selectVariants, $odf/$familyid/$aligner/GATK/candidates/ranking/IEM_Db_CMMS/$familyid.selectVariants"; Supply whole path/file)

-im_db_cc/--Im_Db_Gene_Coverage_Calculation Im_Db_CMMS file coverage calculation (Defaults to "1" (=yes))

-im_db_gidc/--Im_Db_Gene_Id_Col Im_Db_CMMS file gene Id column (zero-based, defaults to "18")

-pSCheck/--samplecheck Flag running check for samples belonging to pedigree (defaults to "1" (=yes) )

-pedigree/--pedigree_file (Supply whole path, defaults to $ods/familyid/familyid_pedigree.txt)

-wgs/--whole_genome_sequencing Analysis to perform are whole genome sequencing data or not (defaults to "0" (=no))

-mc/--maximum_cores The maximum number of cores per node used in the analysis (defaults to "8")

-env_up/--environment_uppmax Sets the environment to UPPMAX. (defaults to "0" (=no))

=head3 I/O

Input format ( infiles_aligned_(merged)_sorted_pmd.bam )

Output format:

1. aligned_sorted_(merged)_chrNo.bam

2. aligned_sorted_(merged)_chrNo_real.bam(.bai)

3. aligned_sorted_(merged)_chrNo_real_recal.bam(.bai)

4. aligned_sorted_(merged)_chrNo_real_recal_resrt.bam(.bai)

5. aligned_sorted_(merged)_allchr_real_recal_resr.bam(.bai)

6. aligned_sorted_(merged)_allchr_real_recal_resr_raw_(SNVorINDEL).vcf(.idx)

7. aligned_sorted_(merged)_allchr_real_recal_resr_varrecal_(SNVorINDEL)_filt.vcf

8. aligned_sorted_(merged)_allchr_real_recal_resr_varrecal_(SNVorINDEL)_filt_annovar

9. aligned_sorted_(merged)_allchr_real_recal_resr_varrecal_(SNVorINDEL)_filt_annovar_all_variants.txt

=head4 Dependencies

Local installation of:
GATK
PicardTools
Annovar

1. varcall_merge_post_annovar_master.1.0.pl

2. compound_filter.pl

=head5 Releases


=cut

use Pod::Usage;
use Pod::Text;
use Getopt::Long;
use POSIX;
use File::Basename;
use File::Spec;

use vars qw($USAGE);

BEGIN {
    $USAGE =
	qq{mip_variation.pl -id [infile...n] -a [projectid] -s [sampleid...n] -em [e-mail] -ods [outdirscripts] -odf [outdirdata] -ids [indirscripts] -rd [referencedir]
	       -i/--infile Infile(s), comma sep (Mandatory: Supply whole path)
               -ids/--indirscript The pipeline custom script in dir (Mandatory: Supply whole path)
               -rd/--referencesdir Reference(s) dir (Mandatory: Supply whole path)
               -gref/--genomeref Flag for setting genomic reference file (defaults to "Homo_sapiens.GRCh37.57.dna.concat.fa")
	       -a/--projectid The project ID (Mandatory)
	       -s/--sampleid The sample ID,comma sep (Mandatory)
	       -em/--email e-mail
               -odf/--outdirdata The data files output directory (Mandatory: Supply whole path)
               -ods/--outdirscript The script files output directory (Mandatory: Supply whole path)
               -familyid/--family Group id of samples to be compared (Mandatory, Ex: 1 for IDN 1-1-1A )
               -pSTV_schr/--samtools_viewschr Flag running samtools view to split per chr & index (defaults to "1" (=yes))
               -pGATK_REAl/--gatk_real Flag running GATK realign (defaults to "1" (=yes))
               -gatk_real_knset1/--gatk_real_knownset1 GATK realign known INDEL set 1 (defaults to "1000G_phase1.indels.hg19.vcf")
               -gatk_real_knset2/--gatk_real_knownset2 GATK realign known INDEL set 2 (defaults to "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf")
               -gatk_path/--genatk_path  Flag for path to GATK, must be supplied for GATK (defaults to "").
               -pGATK_RECAl/--gatk_recal Flag running GATK recalibrate (defaults to "1" (=yes))
               -gatk_recal_knset/--gatk_recal_knownset GATK recal known SNP set (defaults to "dbsnp_135.b37.vcf")
               -pGATK_UNIGT/--gatk_unigt Flag running GATK unifiedgenotyper (defaults to "0" (=no))
               -gatkugt_snp/--gatk_unigt_snp Flag running GATK unifiedgenotyper for SNPs (defaults to "1" (=yes))
               -gatkugt_ind/--gatk_unigt_indel Flag running GATK unifiedgenotyper for INDELs (defaults to "1" (=yes))
               -gatk_bait/--gatk_bait_il Prepared bait interval_list file for GATK_UnifiedGT. (defaults to "SureSelect_All_Exon_50mb_with_annotation_hg19_nochr.bed.pad100.interval_list")
               -pGATK_HAPCAL/--gatk_hapcal Flag running GATK HaplotypeCaller (defaults to "1" (=yes))
               -pGATK_COMBVAR/--gatk_combinevariants Flag running GATK combinevariants. Use only if UnifiedGT and not HaplotypeCaller was used to call variants (SNVs & INDELs).(defaults to "0" (=no))
               -pGATK_VARRECAL/--gatk_varrecalibrator Flag running GATK variantrecalibrator  (defaults to "1" (=yes))
               -gatk_exref_snp/--gatk_exomeref Prepared exome reference file (SNVs) for GATK_Varrecal. (defaults to "all-agilent_50mb-GRCh37-SNPS_pad100_interval_list.vcf")
               -gatk_exref_indel/--gatk_exomeref Prepared exome reference file (INDELs) for GATK_Varrecal. (defaults to "all-agilent_50mb-GRCh37-INDELS_pad100_interval_list.vcf")
               -pPIND/--pindel Flag running Pindel (defaults to "0" (=no))
               -pGATK_VAREVAL_All/--gatk_vareval_all Flag running GATK varianteval for all variants  (defaults to "1" (=yes))
               -pGATK_VAREVAL_Exome/--gatk_vareval_exome Flag running GATK varianteval for exonic variants  (defaults to "1" (=yes))
               -pANVAR/--annovar Flag running annovar (defaults to "1" (=yes))
               -anva_pa/--annovar_path  Path to annovar script directory (Supply whole path, defaults to "". NOTE:Assumes that the annovar db files are located in annovar/humandb)
               -anva_gbv/--annovar_genome_build_version Annovar genome build version (defaults to "hg19")
               -anva_tn/--annovar_table_names Annovar table names (comma sep)
               -anva_stn/--annovar_supported_table_names Annovar MIP supported table names (defaults 0 (=no))
               -anva_maf_th/--annovar_maf_threshold Flag for setting the minor allele frequency threshold in annovar (defaults to "0" (=no))
               -anva_sift_th/--annovar_sift_threshold Flag for setting the avsift threshold in annovar (defaults to "0" (=no))
               -pVMERGE/--variation_annotation_merge Running intersectCollect.pl to merge all annotation info to one file (defaults to "1" (=yes))
               -vm_db_te/--vmerge_db_template Db template file used to create the specific family '-vm_dbf' master file (defaults to "CMMS_intersectCollect_db_master_template.txt")
               -vm_dbf/--vmerge_db_file Db master file to be used in intersectCollect.pl (defaults to "FDN.intersectCollect_db_master.txt")
               -pAddDP/--adddepth Flag for adding depth at nonvariant sites by mpileup and add_depth.pl (defaults to "1" (=yes))             
               -pRankVar/--rankvariants Flag running ranking of variants (defaults to "1" (=yes))
               -rs/--rankscore The rank score cut-off (defaults to "-100")
               -dgf/--dispGeneFiltering Filtering of genes that should be removed from downstream processing (Defaults to "1" (=yes))
               -dgfl/--dispGeneList List of genes that should be removed from downstream processing (Defaults to "IEM_dispGeneList.txt", Format: 1 entry per line;HGNC Symbol)
               -all_db_file/--all_elements_Db_file All_Db file (Defaults to "mart_export_Ensembl_GeneID_key_cleaned_chr.txt")
               -all_db_cc/--all_elements_Db_Gene_Coverage_Calculation All_Db file coverage calculation (Defaults to "1" (=yes))
               -all_db_gidc/--all_elements_Db_Gene_Id_Col All_Db file gene Id column (zero-based, defaults to "4")
               -im_db_file/--Im_Db_CMMS_file Im_Db_CMMS file (Defaults to "IEM_Db_CMMS_version1.2.txt")
               -im_db_te/--Im_Db_template Im_Db template file used to create the specific family '-im_dbmf' master file (Defaults to "select_dbIEM_variants_db_master.txt")
               -im_dbmf/--Im_db_master_file Db master file to be used when selecting variants (defaults to "FDN.intersectCollect_selectVariants_db_master.txt")
               -im_dbfof/--Im_db_file_out_file The file(s) to write to when selecting variants with intersectCollect.pl. Comma sep (defaults to 'odf/familyid/aligner/GATK/candidates/ranking/familyid_orphan.selectVariants, odf/familyid/aligner/GATK/candidates/ranking/IEM_Db_CMMS/familyid.selectVariants'; Supply whole path/file)
               -im_db_cc/--Im_Db_Gene_Coverage_Calculation Im_Db_CMMS file coverage calculation (Defaults to "1" (=yes))
               -im_db_gidc/--Im_Db_Gene_Id_Col Im_Db_CMMS file gene Id column (zero-based, defaults to "18")
               -pSCheck/--samplecheck Flag running check for samples belonging to pedigree (defaults to "1" (=yes))
               -pedigree/--pedigree_file (Supply whole path, defaults to odf/familyid/familyid_pedigree.txt)
               -wgs/--whole_genome_sequencing Analysis to perform are whole genome sequencing data or not (defaults to "0" (=no))
               -mc/--maximum_cores The maximum number of cores per node used in the analysis (defaults to "8")
               -env_up/--environment_uppmax Sets the environment to UPPMAX. (defaults to "0" (=no))
	   };
}

###
#Program parameters
###
my ($aid,$em, $ids, $rd, $annovar_path, $odf, $ods, $fnend, $annovar_genome_build_version, $annovar_supported_table_names, $annovar_maf_threshold, $annovar_sift_threshold, $genomeref, $gatk_real_knset1, $gatk_real_knset2, $gatk_recal_knset, $gatk_unigt_snp, $gatk_unigt_indel, $pedigree, $wgs, $gatk_bait, $gatk_exref_snp, $gatk_exref_indel, $vmerge_db_template, $vm_dbf, $rankscore, $dgf, $dgfl, $all_db_file, $all_db_cc, $all_db_gidc, $im_db_file, $Im_Db_template, $Im_db_master_file, $im_db_cc, $im_db_gidc, $maximum_cores,$environment_uppmax,$gatk_path, $filename, $filename2, $fnt, $fnt2, $aligner, $familyid,$help) = (0,0,0,0,0,0,0, ".sh","hg19", 0, 0, 0, "Homo_sapiens.GRCh37.57.dna.concat.fa", "1000G_phase1.indels.hg19.vcf", "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf","dbsnp_135.b37.vcf",1,1,0,0, "SureSelect_All_Exon_50mb_with_annotation_hg19_nochr.bed.pad100.interval_list", "all-agilent_50mb-GRCh37-SNPS_pad100_interval_list.vcf", "all-agilent_50mb-GRCh37-INDELS_pad100_interval_list.vcf", "CMMS_intersectCollect_db_master_template.txt","FDN.intersectCollect_db_master.txt", -100, 1, "IEM_dispGeneList.txt","mart_export_Ensembl_GeneID_key_cleaned_chr.txt",1,4,"IEM_Db_CMMS_version1.2.txt","select_dbIEM_variants_db_master.txt", "FDN.intersectCollect_selectVariants_db_master.txt", 1,18,8,0); 

###
#Arguments for project
###
my ($pSTV_schr, $pGATK_REAL, $pGATK_RECAL, $pGATK_UNIGT, $pGATK_HAPCAL, $pGATK_VARRECAL, $pGATK_COMBVAR, $pGATK_VAREVAL_All, $pGATK_VAREVAL_Exome, $pANVAR, $pVMERGE, $pAddDP, $pRankVar, $pSCheck, $pPIND) = (1,1,1,0,1,1,0,1,1,1,1,1,1,1,0); #Default arguments for running programs

###
#Staging Area
###

my (@infn,@sid, @chr, @jobID, @Im_db_file_out_file);
my (%infiles,%avcovfn, %dirname, %jobID, %pedigree);
my @script_parameters=@ARGV; #Passes over command line arguments for printing in master_logg since GetOption removes them from ARGV.


#Set default annovar table names
my @annovar_table_names = ("refgene", "mce46way", "gerp++elem", "segdup", "gwascatalog", "tfbs", "mirna", "snp137NonFlagged", "1000g2012apr_all", "hg19_esp6500si_all.txt", "avsift","ljb_pp2","ljb_mt","ljb_lrt", "ljb_gerp++","ljb_phylop");

#Set supported annovar table name filtering options
my @annovar_supported_table_names = ("refgene","knownGene","ensGene","mce46way","gerp++elem","segdup","gwascatalog","tfbs","mirna","snp137","snp135","snp132","snp131","snp130","snp129","snp137NonFlagged","snp135NonFlagged","snp132NonFlagged","snp131NonFlagged","snp130NonFlagged","1000g2012apr_all","1000g2012apr_amr","1000g2012apr_eur","1000g2012apr_asn","1000g2012apr_afr","1000g2012feb_all","hg19_esp6500si_all.txt","hg19_esp6500_all.txt","hg19_esp6500_aa.txt","hg19_esp6500_ea.txt","hg19_esp5400_all.txt","hg19_esp5400_aa.txt","hg19_esp5400_ea.txt","avsift","ljb_sift","ljb_pp2","ljb_mt","ljb_lrt","ljb_all","ljb_gerp++","ljb_phylop"); #Used to print list of supported table names

my %annovar_filtering_option = ( 
    'refgene' => "geneanno",
    'knownGene' => "geneanno",
    'ensGene' => "geneanno",
    'mce46way' => "regionanno",
    'gerp++elem' => "regionanno",
    'segdup' => "regionanno",
    'gwascatalog' => "regionanno",
    'tfbs' => "regionanno",
    'mirna' => "regionanno",
    'snp137' => "filter",
    'snp135' => "filter",
    'snp132' => "filter",
    'snp131' => "filter",
    'snp130' => "filter",
    'snp129' => "filter",
    'snp137NonFlagged' => "filter",    
    'snp135NonFlagged' => "filter",
    'snp132NonFlagged' => "filter",
    'snp131NonFlagged' => "filter",
    'snp130NonFlagged' => "filter",
    '1000g2012apr_all' => "filter",
    '1000g2012apr_amr' => "filter",
    '1000g2012apr_eur' => "filter",
    '1000g2012apr_asn' => "filter",
    '1000g2012apr_afr' => "filter",
    '1000g2012feb_all' => "filter",
    'hg19_esp6500si_all.txt' => "filter",
    'hg19_esp6500_all.txt' => "filter",
    'hg19_esp6500_aa.txt' => "filter",
    'hg19_esp6500_ea.txt' => "filter",
    'hg19_esp5400_all.txt' => "filter",
    'hg19_esp5400_aa.txt' => "filter",
    'hg19_esp5400_ea.txt' => "filter",
    'avsift' => "filter",
    'ljb_sift' => "filter",
    'ljb_pp2' => "filter",
    'ljb_mt' => "filter",
    'ljb_lrt' => "filter",
    'ljb_all' => "filter",
    'ljb_gerp++' => "filter",
    'ljb_phylop' => "filter",
    );
#Set supported annovar table name generic type
my %annovar_generic_filtering_option = ( 
    'hg19_esp6500si_all.txt' => "generic",    
    'hg19_esp6500_all.txt' => "generic",
    'hg19_esp6500_aa.txt' => "generic",
    'hg19_esp6500_ea.txt' => "generic",
    'hg19_esp5400_all.txt' => "generic",
    'hg19_esp5400_aa.txt' => "generic",
    'hg19_esp5400_ea.txt' => "generic",
    );


###
#Genetic Models
###
#Collect parents sampleID and disease status
    my $fatherID =""; #Collect sampleID for father/mother in CreateModels and to be used in compound filtering using compound_filter.pl later
    my $motherID =""; 
    my $fatherDS = 0; #DS = Disease status, 0=Unaffected 
    my $motherDS = 0;

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
#User Options
###

GetOptions('i|infile:s'  => \@infn, #Comma sepatated list
	   'ids|inscriptdir:s'  => \$ids, #Directory for custom scripts required by the pipeline
	   'rd|referencedir:s'  => \$rd, #directory containing references
	   'gref|genomeref:s'  => \$genomeref, #genomic reference file
	   'a|projectid:s'  => \$aid,
	   's|sampleid:s'  => \@sid, #Comma sepatated list, one below outdirdata
	   'em|email:s'  => \$em,
	   'odf|outdirdata:s'  => \$odf, #One above sample id
	   'ods|outdirscript:s'  => \$ods, #One above sample id
	   'familyid|familygroup:s' => \$familyid, #Family group ID (Merged to same vcf file after GATK Base Recalibration) 
	   'pSTV_schr|samtools_viewschr:n' => \$pSTV_schr, #spilt to chr.bam and index
	   'pGATK_REAL|gatk_real:n' => \$pGATK_REAL, #GATK Realign
	   'gatk_real_knset1|gatk_real_knownset1:s' => \$gatk_real_knset1, #Known INDEL set to be used in GATK ReAlign
	   'gatk_real_knset2|gatk_real_knownset2:s' => \$gatk_real_knset2, #Known INDEL set to be used in GATK ReAlign
	   'gatk_path|genatk_path:s' => \$gatk_path, #Path to GATK
	   'pGATK_RECAL|gatk_recal:n' => \$pGATK_RECAL, #GATK Recalibrate
	   'gatk_recal_knset|gatk_recal_knownset:s' => \$gatk_recal_knset, #Known SNP set to be used in GATK Recal
	   'pGATK_UNIGT|gatk_unigt:n' => \$pGATK_UNIGT, #GATK Unifiedgenotyper
	   'gatkugt_snp|gatk_unigt_snp:n' => \$gatk_unigt_snp, #GATK Unifiedgenotyper SNP mode
	   'gatkugt_ind|gatk_unigt_indel:n' => \$gatk_unigt_indel, #GATK Unifiedgenotyper INDEL mode
	   'gatk_bait|gatk_baith_il:s' => \$gatk_bait, #Padded Interval_list to GATK
	   'pGATK_HAPCAL|gatk_hapcal:n' => \$pGATK_HAPCAL, #GATK Haplotypecaller
	   'pGATK_VARRECAL|gatk_varrecalibrator:n' => \$pGATK_VARRECAL, #GATK variantrecalibrator
	   'gatk_exref_snp|gatk_exomeref_snp:s' => \$gatk_exref_snp, #File of 33 exomes to power probabalistic model GATK Varrecal (SNVs) (Recieved from Måns, 120413)
	   'gatk_exref_indel|gatk_exomeref_indel:s' => \$gatk_exref_indel, #File of 33 exomes to power probabalistic model GATK Varrecal (INDELs) (Recieved from Måns, 120413)
	   'pPIND|pindel:n' => \$pPIND, #Pindel (SV variant detection)
	   'pGATK_COMBVAR|gatk_combinevariants:n' => \$pGATK_COMBVAR, #GATK combinevariants
	   'pGATK_VAREVAL_All|gatk_vareval_all:n' => \$pGATK_VAREVAL_All, #GATK varianteval all variants
	   'pGATK_VAREVAL_Exome|gatk_vareval_exome:n' => \$pGATK_VAREVAL_Exome, #GATK varianteval only exonic variants
	   'pANVAR|annovar:n' => \$pANVAR, #Performs annovar filter gene, region and filter analysis
	   'anva_pa|annovar_path:s'  => \$annovar_path, #path to annovar script dir
	   'anva_gbv|annovar_genome_build_version:s'  => \$annovar_genome_build_version,
	   'anva_tn|annovar_table_names:s'  => \@annovar_table_names, #Comma sepatated list
	   'anva_stn|annovar_supported_table_names:n' => \$annovar_supported_table_names, #Generates a list of supported table names
	   'anva_maf_th|annovar_maf_threshold:n' => \$annovar_maf_threshold,
	   'anva_sift_th|annovar_sift_threshold:n' => \$annovar_sift_threshold,	   
	   'pVMERGE|variation_annotation_merge:n' => \$pVMERGE, #Merges annovar analysis results to one master file
	   'vm_db_te|vmerge_db_template:s' => \$vmerge_db_template, #Template file to create the specific family db master file
	   'vm_dbf|vmerge_db_file:s' => \$vm_dbf, #db master file to use when collecting external data
	   'pAddDP|adddepth:n' => \$pAddDP, #Adds depth (DP) for nonvariants to master file (annovar_all.txt)
	   'pRankVar|rankvariants:n' => \$pRankVar, #Ranking variants
	   'rsrankscore:n'  => \$rankscore, #The rank score cut-off
	   'dgf|dispGeneFiltering:n'  => \$dgf, #Enables dispensible gene filtering
	   'dgfl|dispGeneList:s'  => \$dgfl, #List of dispensible genes (1 entry per line; HGNC Symbol)
	   'all_db_file|All_elements_Db_file:s'  => \$all_db_file, #Db of all genes
	   'all_db_cc|All_elements_Db_Gene_Coverage_Calculation:n'  => \$all_db_cc, #Db of all genes for coverage calculation (all features connected to overlapping genes across variant)
	   'all_db_gidc|All_elements_Db_Gene_Id_Col:n'  => \$all_db_gidc, #Db of all genes GeneName column nr zero-based
	   'im_db_file|Im_Db_CMMS_file:s'  => \$im_db_file, #Db of important genes
	   'im_db_te|Im_Db_template:s' => \$Im_Db_template, #Template file to create the specific family selectVariants db master file
	   'im_dbmf|Im_db_master_file:s' => \$Im_db_master_file, #Specific db master file to use when collecting external dataselectingVariants 
	   'im_dbfof|Im_db_file_out_file:s' => \@Im_db_file_out_file, #The intersectCollect select variants output directorys	      
	   'im_db_cc|Im_Db_Gene_Coverage_Calculation:n'  => \$im_db_cc, #Db of important genes coverage calculation (all features connected to overlapping genes across variant)
	   'im_db_gidc|Im_Db_Gene_Id_Col:n'  => \$im_db_gidc, #Db of important genes GeneName column nr zero-based
	   'pSCheck|samplecheck:n' => \$pSCheck, #Check for samples belonging to pedigree
	   'pedigree|pedigree_file:s'  => \$pedigree, #Path to pedigree file location
	   'wgs|whole_genome_sequencing:n' => \$wgs,
	   'mc|maximum_cores:n' => \$maximum_cores, #Per node
	   'env_up|environment_uppmax:n' => \$environment_uppmax, #Sets several default paths, so that they do not have to be supplied 
	   'h|help' => \$help,
	   );

die $USAGE if( $help );

if ($annovar_supported_table_names == 1) {
    print STDOUT "\nThese annovar databases are supported by MIP:\n";
    for (my $annovar_supported_table_name_Counter=0;$annovar_supported_table_name_Counter<scalar(@annovar_supported_table_names);$annovar_supported_table_name_Counter++) {
	print STDOUT $annovar_supported_table_names[$annovar_supported_table_name_Counter], "\n";
    }
    print STDOUT "\n";
    die;
}

if (@infn == 0) {
    my $verbosity = 2;
    print"\n";
    pod2usage({-message => "Must supply an infile directory as comma separated list.\n",
	       -verbose => $verbosity
	      });
}
if ($aid eq 0) {
    
    print STDERR "\n";
    print STDERR "Must supply a project ID", "\n\n";
    die $USAGE;
}
if ( scalar(@sid) eq 0) {
    
    print STDERR "\n";
    print STDERR "Must supply a sample ID as a comma separated list", "\n\n";
    die $USAGE;
}
if ( $ids eq 0) {
    
    if ($environment_uppmax == 1) {
	print STDOUT "\n";
	$ids = "/bubo/proj/$aid/private/mip_scripts_master";
	print STDOUT "Setting the MIP scripts dir to: $ids", "\n\n";
    }
    else {
	print STDERR "\n";
	print STDERR "Must supply the MIP scripts dir", "\n\n";
	die $USAGE;
    }
}
if ( $rd eq 0) {
    
    if ($environment_uppmax == 1) {
	print STDOUT "\n";
	$rd = "/bubo/proj/$aid/private/mip_references";
	print STDOUT "Setting MIP reference dir to: $rd", "\n\n";
    }
    else {
	print STDERR "\n";
	print STDERR "Must supply a MIP reference dir", "\n\n";
	die $USAGE;
    }
}
if ( $familyid eq 0 ) {
    
    print STDERR "\n";
    print STDERR "Must supply a family id. If not applicable supply the same familyid as the sampleid ", "\n\n";
    die $USAGE;
}
if ($odf eq 0) {
    
    if ($environment_uppmax == 1) {
	print STDOUT "\n";
	if ($wgs == 1) {
	    $odf = "/bubo/proj/$aid/private/genomes";
	}
	else {
	    $odf = "/bubo/proj/$aid/private/exomes";
	}
	print STDOUT "Setting MIP output data dir to: $odf", "\n\n";
    }
    else {
	print STDERR "\n";
	print STDERR "Must supply a MIP output data dir", "\n\n";
	die $USAGE;
    }
}
if ($ods eq 0) {
    
    if ($environment_uppmax == 1) {
	print STDOUT "\n";
	if ($wgs == 1) {
	    $ods = "/bubo/proj/$aid/private/genomes_scripts";
	}
	else {
	    $ods = "/bubo/proj/$aid/private/exomes_scripts";
	}
	print STDOUT "Setting MIP output scripts dir to: $ods", "\n\n";
    }
    else {
	print STDERR "\n";
	print STDERR "Must supply a MIP output script dir", "\n\n";
	die $USAGE;
    }
}
if ($pedigree eq 0) {
    
    print STDOUT "\n";
    $pedigree = "$odf/$familyid/$familyid"."_pedigree.txt";
    print STDOUT "Assuming location of pedigree file to be: $pedigree", "\n\n";
    if (-e $pedigree) { #if file exists 
	print STDOUT "Found pedigree file at: $pedigree", "\n\n";
    }
    elsif ($pRankVar eq 1) { 
	print STDERR "Could not find pedigree file at: $pedigree \n";
	print STDERR "Must supply a pedigree file to run Ranking script", "\n\n";
	die $USAGE;
    } 
}
if ( ($pANVAR eq 1) && ($annovar_path eq 0) ) {
    
    print STDERR "\n";
    print STDERR "Must supply the path to the annovar script directory if you want to include annovar in the analysis", "\n\n";
    die $USAGE;
}


#Creates master_logg for the master script 
`mkdir -p $odf/$familyid/master_logg;`; #Creates the master_logg dir
my ($base,$script) = (`date +%Y%m%d`,`basename $0`); #Catches current date and script name
chomp($base,$script); #Remove \n;
my $master_logg_name="$odf/$familyid/master_logg/$script"."_"."$base.txt"; #concatenates master_logg filename
open (MASTERL, ">>$master_logg_name") or die "Can't write to $master_logg_name: $!\n"; #Open file run logg
#Add parameters
print MASTERL "\n$script "; #Adds script name to recontruct command line
foreach (@script_parameters) { print MASTERL "$_ " }; #Adds all passed arguments
print STDOUT "\nScript parameters and info from $script are saved in file: $master_logg_name", "\n";

@infn = split(/,/,join(',',@infn)); #Enables comma separated indir(s)
@sid = split(/,/,join(',',@sid)); #Enables comma sepatated list of sample IDs
@Im_db_file_out_file = split(/,/,join(',',@Im_db_file_out_file)); #Enables comma separated selectVariants output dir(s)
@chr = ("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT"); #Chr for filtering of bam file


InfileReFormat(); #removes .bam ending and extracts filename

if ($pRankVar eq 1) {
    if ( scalar(@Im_db_file_out_file) eq 0) {
	if ($environment_uppmax == 1) {
	    @Im_db_file_out_file = ("$odf/$familyid/$aligner/GATK/candidates/ranking/$familyid"."_orphan.selectVariants","$odf/$familyid/$aligner/GATK/candidates/ranking/IEM_Db_CMMS/$familyid".".selectVariants");
	}
	else {
	    print STDERR "\nSupply the '-Im_db_file_out_file' output file(s) if you want to run 'pRankVar'.\n";
	    die $USAGE;
	}
    }
}

#########################
###Run program part######
#########################

open (MASTERL, ">>$master_logg_name") or die "Can't write to $master_logg_name: $!\n"; #Open file run logg

if ($pSTV_schr eq 1) { #print per chromosome output, ie, from one whole genome bam file per sample, to chr bam files.

    print STDOUT "\nSamtools view split genome to chr", "\n";print MASTERL "\nSamtools view spilt genome to chr", "\n";

    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
    
	SamtoolsViewSChr($sid[$sampleid], $aligner);	
    }
}

if ($pGATK_REAL eq 1) { #Run GATK realign per chr (and sampleid within subroutine)

    print STDOUT "\nGATK RealignerTargetCreator/IndelRealigner", "\n";print MASTERL "\nGATK RealignerTargetCreator/IndelRealigner", "\n";

    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
    
	GATK_real($sid[$sampleid], $aligner);
    }
}

if ($pGATK_RECAL eq 1) { #Run GATK recalibrate per chr (and sampleid within subroutine). Samtools sort and merge is then subsequently done creating a single all_real_recal_resrt.bam file ready for variant calling

    print STDOUT "\nGATK BaseRecalibrator/PrintReads", "\n";print MASTERL "\nGATK BaseRecalibrator/PrintReads", "\n";

    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
    
	GATK_recal($sid[$sampleid], $aligner);	
    }
}

if ($pGATK_UNIGT eq 1) { #Run GATK UnifiedGenoTyper (all.bam for all samples within familyID)
  
    print STDOUT "\nGATK UnifiedGenoTyper", "\n";print MASTERL "\nGATK UnifiedGenoTyper", "\n";

    if ($gatk_unigt_snp == 1) {
	    GATK_unigt($familyid, $aligner, "SNV");    
    }
    if ($gatk_unigt_indel == 1) {
	GATK_unigt($familyid, $aligner, "INDEL");    
    }
    if ( ($gatk_unigt_snp == 0) && ($gatk_unigt_indel == 0) ) {
	print STDOUT "\n\nNOTE:You have choosen to run GATK UnifiedGenoTyper and  specified no run mode (SNV or INDEL)!\n\n";print MASTERL "\n\nYou have choosen to run GATK UnifiedGenoTyper but not specified no run mode (SNV or INDEL)!\n\n";
    }
}

if ($pGATK_HAPCAL eq 1) { #Run GATK HaplotypeCaller (all.bam for all samples within familyID)
    
    print STDOUT "\nGATK HaplotypeCaller", "\n";print MASTERL "\nGATK HaplotypeCaller", "\n";  
    if ($wgs == 0) { #Exome samples   
	GATK_hapcal($familyid, $aligner, "BOTH",0,3,8); #Argument 3 & 4 is where in @chr to start and stop processing. Arg 5 is java heap allocation (Gb).
	GATK_hapcal($familyid, $aligner, "BOTH",3,6,8);
	GATK_hapcal($familyid, $aligner, "BOTH",6,12,4);
	GATK_hapcal($familyid, $aligner, "BOTH",12,18,4);
	GATK_hapcal($familyid, $aligner, "BOTH",18,26,4);
    }
    else { #Whole genome sequencing requires more memory
	GATK_hapcal($familyid, $aligner, "BOTH",0,1,24); #Argument 3 & 4 is where in @chr to start and stop processing. Arg 5 is java heap allocation (Gb).
	GATK_hapcal($familyid, $aligner, "BOTH",1,2,24);
	GATK_hapcal($familyid, $aligner, "BOTH",2,3,24);
	GATK_hapcal($familyid, $aligner, "BOTH",3,4,24);
	GATK_hapcal($familyid, $aligner, "BOTH",4,5,24);
	GATK_hapcal($familyid, $aligner, "BOTH",5,6,24);
	GATK_hapcal($familyid, $aligner, "BOTH",6,7,24);
	GATK_hapcal($familyid, $aligner, "BOTH",7,8,24);
	GATK_hapcal($familyid, $aligner, "BOTH",8,9,24);
	GATK_hapcal($familyid, $aligner, "BOTH",9,10,24);
	GATK_hapcal($familyid, $aligner, "BOTH",10,11,24);
	GATK_hapcal($familyid, $aligner, "BOTH",11,12,24);
	GATK_hapcal($familyid, $aligner, "BOTH",12,13,24);
	GATK_hapcal($familyid, $aligner, "BOTH",13,14,24);
	GATK_hapcal($familyid, $aligner, "BOTH",14,15,24);
	GATK_hapcal($familyid, $aligner, "BOTH",15,16,24);
	GATK_hapcal($familyid, $aligner, "BOTH",16,17,24);
	GATK_hapcal($familyid, $aligner, "BOTH",17,18,24);
	GATK_hapcal($familyid, $aligner, "BOTH",18,19,24);
	GATK_hapcal($familyid, $aligner, "BOTH",19,20,24);
	GATK_hapcal($familyid, $aligner, "BOTH",20,21,24);
	GATK_hapcal($familyid, $aligner, "BOTH",21,22,24);
	GATK_hapcal($familyid, $aligner, "BOTH",22,23,24);
	GATK_hapcal($familyid, $aligner, "BOTH",23,24,24);
	GATK_hapcal($familyid, $aligner, "BOTH",24,25,24);
    }
    GATK_HapCall_ComVar($familyid, $aligner, "BOTH");
}

if ($pGATK_VARRECAL eq 1) { #Run GATK VariantRecalibrator(all_samples_raw.vcf)
  
    print STDOUT "\nGATK VariantRecalibrator/ApplyRecalibration", "\n";print MASTERL "\nGATK VariantRecalibrator/ApplyRecalibration", "\n";
    
    if ( ($gatk_unigt_snp == 1) && ($gatk_unigt_indel == 1) ) {
	GATK_varrecal($familyid, $aligner, "BOTH"); 
    }
    elsif ($gatk_unigt_snp == 1) {
	GATK_varrecal($familyid, $aligner, "SNV");    
    }
    elsif ($gatk_unigt_indel == 1) {
	GATK_varrecal($familyid, $aligner, "INDEL");    
    }
    if ( ($gatk_unigt_snp == 0) && ($gatk_unigt_indel == 0) ) {
	print STDERR "\n\nNOTE:You have choosen to run GATK VariantRecalibrator and  specified no run mode (SNV or INDEL)!\n\n";print MASTERL "\n\nYou have choosen to run GATK VariantRecalibrator but not specified no run mode (SNV or INDEL)!\n\n";
    }    
}

if ($pPIND eq 1) { #Run Pindel (all_chr.bam)
  
    print STDOUT "\nPINDEL", "\n";print MASTERL "\nPINDEL", "\n";
    
    Pindel($familyid, $aligner); #All samples within family is processed together.
}

if ($pGATK_COMBVAR == 1) { #Run GATK CombineVariants(all.bam). Should only be used if UnifiedGT was used to call variants and not HaplotypeCaller. It will overwrite the HaplotypeCaller file if such exists.
  
    print STDOUT "\nGATK CombineVariants", "\n";print MASTERL "\nGATK CombineVariants", "\n";
    
    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {
	
	GATK_combinevariants($sid[$sampleid], $aligner, "BOTH",$familyid);
    }    
}

if ($pANVAR eq 1 ) { #Run annovar
    
    print STDOUT "\nAnnovar Analysis", "\n";print MASTERL "\nAnnovar Analysis", "\n";

    if ( ($gatk_unigt_snp == 0) && ($gatk_unigt_indel == 0) ) {
	print STDERR "\n\nNOTE:You have choosen to run Annovar and specified no run mode (SNV or INDEL)!\n\n";print MASTERL "\n\nYou have choosen to run Annovar but not specified no run mode (SNV or INDEL)!\n\n";
    }
    elsif ( ($gatk_unigt_snp == 1) && ($gatk_unigt_indel == 1) ) { #Run on BOTH vcf file from HaplotypeCaller or combined vcf per sampleID from UnifiedGT (BOTH file is created from SNV and INDEL files). 
	AnnovarFilter($familyid, $aligner, "BOTH");
    }
    elsif ($gatk_unigt_snp == 1) {
	AnnovarFilter($familyid, $aligner, "SNV");    
    }
    elsif ($gatk_unigt_indel == 1) {
	AnnovarFilter($familyid, $aligner, "INDEL");    
    }
        
}

if ($pGATK_VAREVAL_All == 1) { #Run GATK VariantEval(all.bam). Moved to after ANNOVAR since the MT -> M conversion might disturb the vareval process since they use the same infile
  
    print STDOUT "\nGATK VariantEval All Variants", "\n";print MASTERL "\nGATK VariantEval All Variants", "\n";

    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {
	
	if ( ($gatk_unigt_snp == 0) && ($gatk_unigt_indel == 0) ) {
	    print STDERR "\n\nNOTE:You have choosen to run GATK VariantEval All Variants, but specified no run mode (SNV or INDEL)!\n\n";print MASTERL "\n\nYou have choosen to run GATK VariantEval, but not specified no run mode (SNV or INDEL)!\n\n";
	}
	elsif ( ($gatk_unigt_snp == 1) && ($gatk_unigt_indel == 1) ) { #Run on BOTH vcf file from HaplotypeCaller or combined vcf per sampleID from UnifiedGT (BOTH file is created from SNV and INDEL files). 
	    GATK_varianteval($sid[$sampleid], $aligner, "BOTH", 0, $familyid); #0 to keep nr of arguments (instead of EXOME)
	}
	elsif ($gatk_unigt_snp == 1) {
	    GATK_varianteval($sid[$sampleid], $aligner, "SNV", 0, $familyid);    
	}
	elsif ($gatk_unigt_indel == 1) {
	    GATK_varianteval($sid[$sampleid], $aligner, "INDEL", 0, $familyid);    
	}
    }    
}

if ($pVMERGE == 1 ) { #Run varcall_merge_post_annovar_master.pl, Merges all variants across all subjects

    print STDOUT "\nintersectCollect.pl", "\n";print MASTERL "\nintersectCollect.pl", "\n";
    
    if ( ($gatk_unigt_snp == 0) && ($gatk_unigt_indel == 0) ) {
	print STDERR "\n\nNOTE:You have choosen to run  and specified no run mode (SNV or INDEL)!\n\n";print MASTERL "\n\nYou have choosen to run Annovar but not specified no run mode (SNV or INDEL)!\n\n";
    }
    elsif ( ($gatk_unigt_snp == 1) && ($gatk_unigt_indel == 1) ) { #Run on BOTH vcf file from HaplotypeCaller or combined vcf per sampleID from UnifiedGT (BOTH file is created from SNV and INDEL files). 
	VarcallMergePostAnnovar($familyid, $aligner, "BOTH");
    }
    elsif ($gatk_unigt_snp == 1) {
	VarcallMergePostAnnovar($familyid, $aligner, "SNV");    
    }
    elsif ($gatk_unigt_indel == 1) {
	VarcallMergePostAnnovar($familyid, $aligner, "INDEL");    
    }
    
}

if ( $pGATK_VAREVAL_Exome == 1 ) { #Run GATK VariantEval(all.bam) for exome variants only (Must be called after pVMERGE and annotation)
    
    print STDOUT "\nGATK VariantEval Exome Variants", "\n";print MASTERL "\nGATK VariantEval Exome Variants", "\n";
    
    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {
	
	if ( ($gatk_unigt_snp == 0) && ($gatk_unigt_indel == 0) ) {
	    print STDERR "\n\nNOTE:You have choosen to run GATK VariantEval (Exome) and  specified no run mode (SNV or INDEL)!\n\n";print MASTERL "\n\nYou have choosen to run GATK VariantEval (Exome) but not specified no run mode (SNV or INDEL)!\n\n";
	}
	elsif ( ($gatk_unigt_snp == 1) && ($gatk_unigt_indel == 1) ) { #Run on BOTH vcf file from HaplotypeCaller or combined vcf per sampleID from UnifiedGT (BOTH file is created from SNV and INDEL files). 
	    GATK_varianteval($sid[$sampleid], $aligner, "BOTH", "EXOME", $familyid);
	}
	elsif ( ($gatk_unigt_snp == 1) && ($pGATK_COMBVAR ==0) ) {
	    GATK_varianteval($sid[$sampleid], $aligner, "SNV", "EXOME", $familyid);    
	}
	elsif ( ($gatk_unigt_indel == 1) && ($pGATK_COMBVAR ==0) ) {
	    GATK_varianteval($sid[$sampleid], $aligner, "INDEL", "EXOME", $familyid);    
	}
    }    
}

if ( $pAddDP == 1 ) { #Add depth (DP) for nonvariants to masterfile (annovar_all.txt)
    
    print STDOUT "\nadd_depth.pl for nonvariants (masterfile)", "\n";print MASTERL "\nadd_depth.pl for nonvariants (masterfile)", "\n";
    
    if ( ($gatk_unigt_snp == 0) && ($gatk_unigt_indel == 0) ) {
	print STDERR "\n\nNOTE:You have choosen to run AddDP and specified no run mode (SNV or INDEL)!\n\n";print MASTERL "\n\nYou have choosen to run Annovar but not specified no run mode (SNV or INDEL)!\n\n";
    }
    elsif ( ($gatk_unigt_snp == 1) && ($gatk_unigt_indel == 1) ) { #Run on BOTH vcf file from HaplotypeCaller or combined vcf per sampleID from UnifiedGT (BOTH file is created from SNV and INDEL files). 
	AddDP($familyid, $aligner, "BOTH");
    }
    elsif ($gatk_unigt_snp == 1) {
	AddDP($familyid, $aligner, "SNV");    
    }
    elsif ($gatk_unigt_indel == 1) {
	 AddDP($familyid, $aligner, "INDEL");    
    }
}

if ( $pRankVar == 1 ) { #Ranking of variants
    
    print STDOUT "\nRank Variants", "\n";print MASTERL "\nRank Variants", "\n";
    
    if ( ($gatk_unigt_snp == 0) && ($gatk_unigt_indel == 0) ) {
	print STDERR "\n\nNOTE:You have choosen to run Filter variants and specified no run mode (SNV or INDEL)!\n\n";print MASTERL "\n\nYou have choosen to run Annovar but not specified no run mode (SNV or INDEL)!\n\n";
    }
    elsif ( ($gatk_unigt_indel == 1) && ($gatk_unigt_snp == 1) ) {
	RankVar($familyid, $aligner, "BOTH");    
    }
    elsif ($gatk_unigt_snp == 1) {
	RankVar($familyid, $aligner, "SNV");    
    }
    elsif ($gatk_unigt_indel == 1) {
	RankVar($familyid, $aligner, "INDEL");    
    }
}

if ( $pSCheck == 1 ) { #Tests sample for correct relatives (only performed for samples with relatives defined in pedigree file) performed on sequence data. Gender is checked in the collect_info.pl script by taking the chrX/chrY ration and expecting females to have a ratio > 5.   
    
    print STDOUT "\nSample check (Gender & Relatives)", "\n";print MASTERL "\nSample check (Gender & Relatives)", "\n";
    if ( ($gatk_unigt_indel == 1) && ($gatk_unigt_snp == 1) ) {
	SampleCheck($familyid, $aligner, "BOTH");    
    }
    elsif ($gatk_unigt_snp == 1) {
	SampleCheck($familyid, $aligner, "SNV");    
    }   
}

close(MASTERL); #Close Master_logg file

######################
###Sub Routines#######
######################

sub SampleCheck { 
#Tests sample for correct relatives (only performed for samples with relatives defined in pedigree file) performed on sequence data.
#$_[0] = $familyid
#$_[1] = $aligner, to choose the correct dir depending on what aligner has been used previously
#$_[2] = SNV or BOTH

    `mkdir -p $odf/$_[0]/$_[1]/info;`; #Creates the aligner folder and info data file directory
    `mkdir -p $odf/$_[0]/$_[1]/samplecheck`; #Creates the aligner folder, Samplecheck data file directory
    `mkdir -p $ods/$_[0]/$_[1]`; #Creates the aligner folder script file directory
    $filename = "$ods/$_[0]/$_[1]/samplecheck_$_[0]_$_[2]."; 

    Checkfnexists($filename, $fnend);

#Info and Logg
    print STDOUT "Creating sbatch script Sample check (Gender & Relatives) and writing script file(s) to: ", $filename, "\n";print MASTERL "Creating sbatch script Sample check (Gender & Relatives) and writing script file(s) to: ", $filename, "\n";
    print STDOUT "Sbatch script Sample check (Gender & Relatives) data files will be written to: ", $odf,"/$_[0]/$_[1]/samplecheck", "\n";print MASTERL "Sbatch script Sample check (Gender & Relatives) data files will be written to: ", $odf,"/$_[0]/$_[1]/samplecheck", "\n";

    open (SCHECK, ">$filename") or die "Can't write to $filename: $!\n";
    
    print SCHECK "#! /bin/bash -l", "\n";
    print SCHECK "#SBATCH -A ", $aid, "\n";
    print SCHECK "#SBATCH -n 1", "\n";
    print SCHECK "#SBATCH -C thin", "\n";	
    print SCHECK "#SBATCH -t 1:00:00", "\n";

    print SCHECK "#SBATCH -J samplecheck_$_[0]_", "\n";
    print SCHECK "#SBATCH -e $odf/$_[0]/$_[1]/info/samplecheck_$_[0]_$_[2].", $fnt ,".stderr.txt", "\n";
    print SCHECK "#SBATCH -o $odf/$_[0]/$_[1]/info/samplecheck_$_[0]_$_[2].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print SCHECK "#SBATCH --mail-type=END", "\n";
	print SCHECK "#SBATCH --mail-type=FAIL", "\n";
	print SCHECK "#SBATCH --mail-user=$em", "\n\n";	
    }
    
    print SCHECK 'echo "Running on: $(hostname)"',"\n\n";
    print SCHECK "#Directories", "\n";
    print SCHECK 'inFamilyDir="', "$odf/$_[0]/$_[1]/GATK", '"', "\n";
    print SCHECK 'outFamilyDir="', "$odf/$_[0]/$_[1]/samplecheck", '"', "\n";
    
    print SCHECK "#Create Plink .ped and .map file per family using vcfTools","\n";
    print SCHECK q?vcftools --vcf ${inFamilyDir}/?.$_[0].q?_allchr_real_recal_resrt_varrecal_?.$_[2].q?_filt.vcf --plink --out ${outFamilyDir}/?.$_[0], "\n\n";
    print SCHECK "#Create vcfTools inbreeding coefficient F per family using vcfTools","\n";
    print SCHECK q?vcftools --vcf ${inFamilyDir}/?.$_[0].q?_allchr_real_recal_resrt_varrecal_?.$_[2].q?_filt.vcf --het --out ${outFamilyDir}/?.$_[0], "\n\n";
    print SCHECK "#Create Plink .mibs per family","\n"; 
    print SCHECK q?plink --noweb --ped ${outFamilyDir}/?.$_[0].q?.ped --map ${outFamilyDir}/?.$_[0].q?.map --cluster --matrix --out ${outFamilyDir}/?.$_[0], "\n\n";
    print SCHECK "#Create Plink sexcheck per family","\n"; 
    print SCHECK q?plink --noweb --ped ${outFamilyDir}/?.$_[0].q?.ped --map ${outFamilyDir}/?.$_[0].q?.map --check-sex --out ${outFamilyDir}/?.$_[0], "\n\n";
    
    print SCHECK "\n\nwait", "\n\n";    
    close(SCHECK); 
    FIDSubmitJob(0,$familyid, 2, $_[2],$filename);
    return;
}

sub RankVar { 
#Filter and Rank variants depending on mendelian inheritance, frequency and phenotype using rank_filter:chr.pl
#$_[0] = $familyid NOTE: not sampleid
#$_[1] = aligner
#$_[2] = $gatk_unigt_snps (SNV), $gatk_unigt_indels (INDEL) or BOTH
    
    #if ($_[2] eq "BOTH") {
    `mkdir -p $odf/$_[0]/$_[1]/info;`; #Creates the alignment folder and info data file directory
    `mkdir -p $ods/$_[0]/$_[1];`; #Creates the aligner script folder 
    $filename = "$ods/$_[0]/$_[1]/rank_var_$_[2]_$_[0].";
    Checkfnexists($filename, $fnend);

#Info and Logg
    print STDOUT "Creating sbatch script for Ranking Variants and writing script file(s) to: ", $filename, "\n";print MASTERL "Creating sbatch script for Ranking Variants and writing script file(s) to: ", $filename, "\n";
    
    for (my $Im_db_file_out_fileCounter=0;$Im_db_file_out_fileCounter<scalar(@Im_db_file_out_file);$Im_db_file_out_fileCounter++) {
	my ($volume,$directories,$file) = File::Spec->splitpath( $Im_db_file_out_file[$Im_db_file_out_fileCounter] );
	`mkdir -p $directories;`; #Creates the ranking/Db selection folders
	print STDOUT "Ranking Variants data files will be written to: $directories$_[0]_ranked_$_[2].txt", "\n";print MASTERL "Ranking Variants data files will be written to: $directories$_[0]_ranked_$_[2].txt", "\n";
    }    
    
    open (RV, ">$filename") or die "Can't write to $filename: $!\n";
    
    print RV "#! /bin/bash -l", "\n";
    print RV "#SBATCH -A ", $aid, "\n";
    print RV "#SBATCH -p node -n 1", "\n";
    print RV "#SBATCH -C thin", "\n";	
    print RV "#SBATCH -t 05:00:00", "\n";
    print RV "#SBATCH -J RV_$_[2]_$_[0]", "\n";
    print RV "#SBATCH -e $odf/$_[0]/$_[1]/info/rank_var_$_[2]_$_[0].", $fnt ,".stderr.txt", "\n";
    print RV "#SBATCH -o $odf/$_[0]/$_[1]/info/rank_var_$_[2]_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	print RV "#SBATCH --mail-type=END", "\n";
	print RV "#SBATCH --mail-type=FAIL", "\n";
	print RV "#SBATCH --mail-user=$em", "\n\n";
    }
    
    print RV 'echo "Running on: $(hostname)"',"\n\n";
    print RV 'referenceArchive="', "$rd", '"', "\n\n"; 
    print RV "#Directories", "\n";
    print RV 'inFamilyDir="',"$odf/$_[0]/$_[1]/GATK", '"', "\n"; # Location of _annovar_all_variants

    print RV "#Create db master file to select variants from template", "\n";
    my $nrColumns; #Total Nr of columns 
    my $nrAnnotationColumns; #The number of columns containing annotation info
    my $pNrofCol; #For perl regexp
    if (-e "$odf/$_[0]/$_[1]/GATK/$_[0]_allchr_real_recal_resrt_varrecal_$_[2]_filt_annovar_all_variants.txt") { #Check if the exists (rerun actual data to sample from) 
	$pNrofCol = q?perl -nae 'if ($_=~/^#/ ) { chomp($_); my @nr_of_columns=split("\t",$_); print scalar(@nr_of_columns);last; }' ?;
	$nrColumns = `$pNrofCol $odf/$_[0]/$_[1]/GATK/$_[0]_allchr_real_recal_resrt_varrecal_$_[2]_filt_annovar_all_variants.txt;`;
	$nrAnnotationColumns = $nrColumns - scalar(@sid);
    }
    elsif (-e "$odf/$_[0]/$vm_dbf") { #First analysis run - no actual data file exists - locate IDN columns from family specific template file (if defined)
	$pNrofCol = q?perl -nae 'if ($_=~/^outinfo/ || $_=~/^outheaders/ ) { chomp($_); my @nr_of_columns=split(",",$_); print scalar(@nr_of_columns);last; }' ?;
	$nrColumns = `$pNrofCol $odf/$_[0]/$vm_dbf;`;
	$nrAnnotationColumns = $nrColumns - scalar(@sid);
    }
    elsif (-e "$rd/$vmerge_db_template") { #No information on previous intersectCollect to create annovar_all_variants file - locate IDN columns from unspecific interSect db template file
	$pNrofCol = q?perl -nae 'if ($_=~/^outinfo/ || $_=~/^outheaders/ ) { chomp($_); my @nr_of_columns=split(",",$_); print scalar(@nr_of_columns);last; }' ?;
	$nrAnnotationColumns = `$pNrofCol $rd/$vmerge_db_template;`-1; #"-1" Since IDN is already factored in from the regexp
	$nrColumns = $nrAnnotationColumns + scalar(@sid);
    }
    else {
	print STDERR "Could not estimate location of IDN columns from variant file, nor from templates ('-vm_dbf' or '-vmerge_db_template'). Please provide this information to run 'rankvariants'.", "\n";
	die;
    }
    
    my $sampleIDcolcond = $nrColumns-1; #To write last IDN entry without "," at the end
    $Im_db_master_file =~ s/FDN/$_[0]/g; #Exchange FND for the real familyID
    #Add relative path to db_template for variant file(s) 
    my ($volume,$directories,$file) = File::Spec->splitpath( $odf );
    my @dirs = File::Spec->splitdir( $directories );
    my $regexp_odf;
    for (my $dirs_Count=1;$dirs_Count<scalar(@dirs);$dirs_Count++) {
	
	$regexp_odf .= "\\/".$dirs[$dirs_Count]; #Create escape char for / in later regexp
    }
    $regexp_odf .= $file;
    #Add relative path to db_template for reference/db files
    ($volume,$directories,$file) = File::Spec->splitpath( $rd );
    @dirs = File::Spec->splitdir( $directories );
    my $regexp_rd;	
    for (my $dirs_Count=1;$dirs_Count<scalar(@dirs);$dirs_Count++) {
	
	$regexp_rd .= "\\/".$dirs[$dirs_Count]; #Create escape char for / in later regexp
    }
    $regexp_rd .= $file;
    #Create family specific template
    print RV q?perl -nae 'if ($_=~/outinfo/i) { if ($_=~/IDN/) { my $sidstring; for (my $sampleID=?.$nrAnnotationColumns.q?;$sampleID<?.$nrColumns.q?;$sampleID++) { if ($sampleID<?.$sampleIDcolcond.q?) { $sidstring.="IDN:Filter:GT=Genotype:AD=Allelic_depths_for_the_ref_and_alt_alleles:GQ=Genotype Quality:PL=Normalized_Phred-scaled_likelihoods_for_genotypes=>0_$sampleID,"} else { $sidstring.="IDN:Filter:GT=Genotype:AD=Allelic_depths_for_the_ref_and_alt_alleles:GQ=Genotype Quality:PL=Normalized_Phred-scaled_likelihoods_for_genotypes=>0_$sampleID"} } s/IDN/$sidstring/g; print $_;} next;} if ($_=~/outcolumns/i) { if ($_=~/IDN/) { my $sidstring; for (my $sampleID=?.$nrAnnotationColumns.q?;$sampleID<?.$nrColumns.q?;$sampleID++) { if ($sampleID<?.$sampleIDcolcond.q?) { $sidstring.="0_$sampleID,"} else { $sidstring.="0_$sampleID"} } s/IDN/$sidstring/g; print $_;} next;} if ($_=~/outheaders/i) { if ($_=~/IDN/) { my $sidstring; for (my $sampleID=?.$nrAnnotationColumns.q?;$sampleID<?.$nrColumns.q?;$sampleID++) { if ($sampleID<?.$sampleIDcolcond.q?) { $sidstring.="IDN:Filter:GT=Genotype:AD=Allelic_depths_for_the_ref_and_alt_alleles:GQ=Genotype Quality:PL=Normalized_Phred-scaled_likelihoods_for_genotypes,"} else { $sidstring.="IDN:Filter:GT=Genotype:AD=Allelic_depths_for_the_ref_and_alt_alleles:GQ=Genotype Quality:PL=Normalized_Phred-scaled_likelihoods_for_genotypes"} } s/IDN/$sidstring/g; print $_;} next;} elsif ($_=~s/FDN/?.$_[0].q?/g) { if($_=~s/^ODF/?.$regexp_odf.q?/g) {} if($_=~s/ALIGNER/?.$aligner.q?/g) {} if ($_=~/IDN/) { my $sidstring; for (my $sampleID=?.$nrAnnotationColumns.q?;$sampleID<?.$nrColumns.q?;$sampleID++) { if ($sampleID<?.$sampleIDcolcond.q?) { $sidstring.="$sampleID,"} else { $sidstring.="$sampleID"} } s/IDN/$sidstring/g; print $_;} else { print $_;} } else { if($_=~s/^RD/?.$regexp_rd.q?/g) {} print $_;}' ?.$rd.q?/?.$Im_Db_template.q? > ?."$odf/$_[0]/$Im_db_master_file", "\n\n";

    
    #if ($_[2] eq "BOTH") {
    my $hapcal_both_file = "$odf/$_[0]/$_[1]/GATK/$_[0]_allchr_real_recal_resrt_varrecal_$_[2]_filt_annovar_all_variants.txt";
    if ( $pAddDP == 0 ) {
	#Add chr
	print RV "#Add chr", "\n";
	print RV q?perl -i -p -e ' if($_=~/^#/) {} else {s/^(.+)/chr$1/g }' ?, '${inFamilyDir}', "/$_[0]_allchr_real_recal_resrt_varrecal_$_[2]_filt_annovar_all_variants.txt", "\n\n";
    }
###
#Only IM_Db_genes
###
    
    if ( ($pGATK_HAPCAL eq 1) || (-e $hapcal_both_file) ) { #HaplotypeCaller has been used in present call or previously
#Create temp_file containing only IEM_Db_genes (to avoid duplicates in ranked list)
	print RV "#Create temp_file containing only IEM_Db_genes (to avoid duplicates in ranked list)", "\n";
	print RV "perl $ids/intersectCollect.pl -db $odf/$_[0]/$Im_db_master_file -s 1 -sofs ";
	for (my $Im_db_file_out_fileCounter=0;$Im_db_file_out_fileCounter<scalar(@Im_db_file_out_file);$Im_db_file_out_fileCounter++) {
	    if ($Im_db_file_out_fileCounter eq scalar(@Im_db_file_out_file)-1) {
		print RV $Im_db_file_out_file[$Im_db_file_out_fileCounter], " \n\n";
	    }
	    else {
		print RV $Im_db_file_out_file[$Im_db_file_out_fileCounter], ",";
	    }
	}
	#Ranking
	print RV "#Ranking", "\n";
	for (my $Im_db_file_out_fileCounter=1;$Im_db_file_out_fileCounter<scalar(@Im_db_file_out_file);$Im_db_file_out_fileCounter++) { #Skip orphan file and run selected files
	    print RV "perl $ids/rank_list_filter.pl -i $Im_db_file_out_file[$Im_db_file_out_fileCounter] -cmms_imdb 1 -dgf $dgf -dgfl ", '${referenceArchive}',"/$dgfl -im_db_file ", '${referenceArchive}',"/$im_db_file -im_db_cc $im_db_cc -im_db_gidc $im_db_gidc  -rs $rankscore -pedigree $pedigree -tarcov ";
	    
	    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {#For all sample ids 
		my $tempinfile = $avcovfn{$sid[$sampleid]};
		
		if ($sampleid eq scalar(@sid)-1) {
		    
		    print RV "$odf/$sid[$sampleid]/$_[1]/coverageReport/$tempinfile","_IEM_target_coverage.txt ";
		}
		else {
		    print RV "$odf/$sid[$sampleid]/$_[1]/coverageReport/$tempinfile","_IEM_target_coverage.txt,";	
		    
		}
	    }
	    ($volume,$directories,$file) = File::Spec->splitpath( $Im_db_file_out_file[$Im_db_file_out_fileCounter] ); #Create outfile
	    print RV q?-o ?.$directories.$_[0].q?_ranked_?.$_[2].q?.txt?, "\n\n";
	}
    }
    
###
#Create a Mosaic BAM file for viewing only relevant variants related to dbIEM 
###
    print RV q?cp ?.$directories.$_[0].q?_ranked_?.$_[2].q?.txt ?.$directories.$_[0].q?_ranked_?.$_[2].q?_temp.txt?, "\n\n";
    print RV q?perl -i -p -e 'unless ($_=~/^#/) { s/^chr(.+)/$1/g }' ?.$directories.$_[0].q?_ranked_?.$_[2].q?_temp.txt?, "\n\n"; #Remove chr for intersect with BAM
    
    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {#For all sampleIDs
	my $tempinfile = $avcovfn{$sid[$sampleid]};
	
	print RV q?intersectBed -wa -abam ?.$odf.q?/?.$sid[$sampleid].q?/?.$_[1].q?/?.$tempinfile.q?.bam -b ?.$directories.$_[0].q?_ranked_?.$_[2].q?_temp.txt > ?.$directories.$tempinfile.q?.bam &?, "\n\n";		
    }
    print RV "wait\n\n";
###	
#No IM_Db_genes
###
    
    #Ranking
    print RV "#Ranking", "\n";
    print RV "perl $ids/rank_list_filter.pl -i $Im_db_file_out_file[0] -dgf $dgf -dgfl ", '${referenceArchive}',"/$dgfl -im_db_file ", '${referenceArchive}',"/$all_db_file -im_db_cc $all_db_cc -im_db_gidc $all_db_gidc -rs $rankscore -pedigree $pedigree -tarcov ";
    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {#For all sampleIDs
	my $tempinfile = $avcovfn{$sid[$sampleid]};
	
	if ($sampleid eq scalar(@sid)-1) {
	    
	    print RV "$odf/$sid[$sampleid]/$_[1]/coverageReport/$tempinfile","_target_coverage.txt ";
	}
	else {
	    print RV "$odf/$sid[$sampleid]/$_[1]/coverageReport/$tempinfile","_target_coverage.txt,";	
	    
	}
    }
    ($volume,$directories,$file) = File::Spec->splitpath( $Im_db_file_out_file[0] ); #Create outfile path
    print RV q?-o ?.$directories.$_[0].q?_ranked_?.$_[2].q?.txt?, "\n\n";
    
    close(RV);   
    FIDSubmitJob(0,$familyid, 1, $_[2],$filename);
    return;
}

sub AddDP { 
#Adds depth (DP) for all nonvariants pos for all chr (and subjects) to create a master file containing all annovar information and DP for nonvariants in annovar_all.txt master file
#$_[0] = $familyid NOTE: not sampleid
#$_[1] = aligner
#$_[2] = $gatk_unigt_snps (SNV) or $gatk_unigt_indels (INDEL)

    `mkdir -p $odf/$_[0]/$_[1]/info;`; #Creates the alignment folder and info data file directory
    `mkdir -p $odf/$_[0]/$_[1]/GATK;`; #Creates the aligner and annovar folder   
    `mkdir -p $ods/$_[0]/$_[1];`; #Creates the aligner script folder    
    
    $filename = "$ods/$_[0]/$_[1]/add_depth_$_[2]_$_[0].";
    Checkfnexists($filename, $fnend);
    
#Info and Logg
    print STDOUT "Creating sbatch script add_depth.pl and writing script file(s) to: ", $filename, "\n";print MASTERL "Creating sbatch script add_depth.pl and writing script file(s) to: ", $filename, "\n";
    print STDOUT "Sbatch script add_depth.pl data files will be written to: ", $odf,"/$_[0]/$_[1]/GATK", "\n";print MASTERL "Sbatch script add_depth.pl data files will be written to: ", $odf,"/$_[0]/$_[1]/GATK", "\n";

    open (ADDDP, ">$filename") or die "Can't write to $filename: $!\n";
    
    print ADDDP "#! /bin/bash -l", "\n";
    print ADDDP "#SBATCH -A ", $aid, "\n";
    print ADDDP "#SBATCH -p node -n 1", "\n";
    print ADDDP "#SBATCH -C thin", "\n";	
    print ADDDP "#SBATCH -t 10:00:00", "\n";
    print ADDDP "#SBATCH -J ADDDP_$_[2]_$_[0]", "\n";
    print ADDDP "#SBATCH -e $odf/$_[0]/$_[1]/info/add_depth_$_[2]_$_[0].", $fnt ,".stderr.txt", "\n";
    print ADDDP "#SBATCH -o $odf/$_[0]/$_[1]/info/add_depth_$_[2]_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print ADDDP "#SBATCH --mail-type=END", "\n";
	print ADDDP "#SBATCH --mail-type=FAIL", "\n";
	print ADDDP "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print ADDDP 'echo "Running on: $(hostname)"',"\n\n";
    print ADDDP "#Directories", "\n";
    print ADDDP 'inFamilyDir="',"$odf/$_[0]/$_[1]/GATK", '"', "\n"; # All variants for all subjects have been merged here by chr
    print ADDDP 'outFamilyDir="', "$odf/$_[0]/$_[1]/GATK", '"', "\n\n";
#Find all "./." per sample ID and print chr pos to new file (mpileup -l format)
    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) { #For all sample ids, find nonvariants
	
	my $tempinfile = $avcovfn{$sid[$sampleid]};
	
	print ADDDP "#Find all - (nonvariants)", "\n";
	print ADDDP q?perl -F'\t' -nae' if ($_=~ /?.$sid[$sampleid].q?\S+\.\/\./ ) { print "$F[0] $F[1]","\n"; }' ?, '${inFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_varrecal_$_[2]_filt_annovar_all_variants.txt > ", '${outFamilyDir}', "/$sid[$sampleid]_nonvariants.txt", "\n\n";
#Remove chr (dependent on chrom name in reference)
	print ADDDP "#Remove chr", "\n";
	print ADDDP q?perl -i -p -e 's/^chr(.+)/$1/g' ?, '${inFamilyDir}', "/$sid[$sampleid]_nonvariants.txt", "\n\n";
	#Indir for sample BAM
	print ADDDP "#Samples indir (BAM)", "\n\n";
	print ADDDP 'inSampleDir_2="',"$odf/$sid[$sampleid]/$_[1]", '"', "\n";
#Find depth (Only proper pairs)
	print ADDDP "samtools mpileup -A -l ", '${inFamilyDir}', "/$sid[$sampleid]_nonvariants.txt ", '${inSampleDir_2}', "/$tempinfile", q?.bam | perl -F'\t' -nae' print $F[0],"\t", $F[1],"\t", $F[3], "\n";' > ?, '${outFamilyDir}', "/$sid[$sampleid]_mpileup_nonvariants.txt", "\n\n";
#Add chr again for annovar master file uses chr
	print ADDDP "#Add chr", "\n";
	print ADDDP q?perl -i -p -e 's/^(.+)/chr$1/g' ?, '${inFamilyDir}', "/$sid[$sampleid]_mpileup_nonvariants.txt", "\n\n";
    }
    #Add chr to master file
    print ADDDP q?perl -i -p -e ' if($_=~/^#/) {} else {s/^(.+)/chr$1/g }' ?, '${inFamilyDir}', "/$_[0]_allchr_real_recal_resrt_varrecal_$_[2]_filt_annovar_all_variants.txt", "\n\n";
    print ADDDP "perl $ids/add_depth.pl -i ", , '${inFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_varrecal_$_[2]_filt_annovar_all_variants.txt -infnv ";
    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {#For all sample ids mpileup nonvariant files
	
	if ($sampleid eq scalar(@sid)-1) {
	    
	    print ADDDP '${inFamilyDir}',"/$sid[$sampleid]_mpileup_nonvariants.txt ";
	}
	else {
	    print ADDDP '${inFamilyDir}',"/$sid[$sampleid]_mpileup_nonvariants.txt,";	
	    
	}
    }
    print ADDDP "-sid "; #SampleIDs 
    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {#For all sample ids mpileup nonvariant files
	
	if ($sampleid eq scalar(@sid)-1) {
	    
	    print ADDDP "$sid[$sampleid] ";
	}
	else {
	    print ADDDP "$sid[$sampleid],";	
	    
	}
    }
    print ADDDP "-o ", '${inFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_varrecal_$_[2]_filt_annovar_all_variants.txt", "\n\n"; #Overwrites original _annovar_all.txt file
    	    
    close(ADDDP);   
    FIDSubmitJob(0,$familyid, 1, $_[2],$filename);
    return;
}
   
sub VarcallMergePostAnnovar { 
#Merges all variants for all chr (and subjects) to create a master file containing all annovar information
#$_[0] = $familyid NOTE: not sampleid
#$_[1] = aligner
#$_[2] = $gatk_unigt_snps (SNV) or $gatk_unigt_indels (INDEL)

    `mkdir -p $odf/$_[0]/$_[1]/info;`; #Creates the alignment folder and info data file directory
    `mkdir -p $odf/$_[0]/$_[1]/GATK;`; #Creates the aligner and annovar folder   
    `mkdir -p $ods/$_[0]/$_[1];`; #Creates the aligner script folder    
    
    $filename = "$ods/$_[0]/$_[1]/intersectCollect_$_[2]_$_[0].";
    Checkfnexists($filename, $fnend);
    
#Info and Logg
    print STDOUT "Creating sbatch script intersectCollect.sh and writing script file(s) to: ", $filename, "\n";print MASTERL "Creating sbatch script intersectCollect.sh and writing script file(s) to: ", $filename, "\n";
    print STDOUT "Sbatch script intersectCollect data files will be written to: ", $odf,"/$_[0]/$_[1]/GATK", "\n";print MASTERL "Sbatch script intersectCollect data files will be written to: ", $odf,"/$_[0]/$_[1]/GATK", "\n";

    open (VMERGE, ">$filename") or die "Can't write to $filename: $!\n";
    
    print VMERGE "#! /bin/bash -l", "\n";
    print VMERGE "#SBATCH -A ", $aid, "\n";
    print VMERGE "#SBATCH -p node -n 1", "\n";
    print VMERGE "#SBATCH -C thin", "\n";	
    print VMERGE "#SBATCH -t 04:00:00", "\n";
    print VMERGE "#SBATCH -J VMERGE_$_[2]_$_[0]", "\n";
    print VMERGE "#SBATCH -e $odf/$_[0]/$_[1]/info/intersectCollect_$_[2]_$_[0].", $fnt ,".stderr.txt", "\n";
    print VMERGE "#SBATCH -o $odf/$_[0]/$_[1]/info/intersectCollect_$_[2]_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print VMERGE "#SBATCH --mail-type=END", "\n";
	print VMERGE "#SBATCH --mail-type=FAIL", "\n";
	print VMERGE "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print VMERGE 'echo "Running on: $(hostname)"',"\n\n";
    print VMERGE 'referenceArchive="', "$rd", '"', "\n\n"; 
    print VMERGE "#Directories", "\n";
    print VMERGE 'inFamilyDir="',"$odf/$_[0]/$_[1]/GATK", '"', "\n"; # All variants for all subjects have been merged here by chr
    print VMERGE 'inFamilyPrefix="',"/$_[0]", "_allchr_real_recal_resrt_varrecal_$_[2]_filt_annovar", '"', "\n";
    print VMERGE 'outFamilyDir="', "$odf/$_[0]/$_[1]/GATK", '"', "\n\n";
    
    
    print VMERGE "#Create db master file from template", "\n";
    my $sampleIDcolumns = scalar(@sid)+5; #Requires CMMS format (chr,start,stop,ref_allele,alt_allel,IDN...)
    my $sampleIDcolcond = scalar(@sid)+4;
    $vm_dbf =~ s/FDN/$_[0]/g; #Exchange FND for the real familyID
    #Add relative path to db_template for annovar files 
    my ($volume,$directories,$file) = File::Spec->splitpath( $odf );
    my @dirs = File::Spec->splitdir( $directories );
    my $regexp_odf;
    for (my $dirs_Count=1;$dirs_Count<scalar(@dirs);$dirs_Count++) {
	
	$regexp_odf .= "\\/".$dirs[$dirs_Count]; #Create escape char for / in later regexp
    }
    $regexp_odf .= $file;
    #Add relative path to db_template for reference files
    ($volume,$directories,$file) = File::Spec->splitpath( $rd );
    @dirs = File::Spec->splitdir( $directories );
    my $regexp_rd;	
    for (my $dirs_Count=1;$dirs_Count<scalar(@dirs);$dirs_Count++) {
	
	$regexp_rd .= "\\/".$dirs[$dirs_Count]; #Create escape char for / in later regexp
    }
    $regexp_rd .= $file;
    #Create family specific template
    print VMERGE q?perl -nae 'if ($_=~/outinfo/i) { if ($_=~/IDN/) { my $sidstring; for (my $sampleID=5;$sampleID<?.$sampleIDcolumns.q?;$sampleID++) { if ($sampleID<?.$sampleIDcolcond.q?) { $sidstring.="IDN:Filter:GT=Genotype:AD=Allelic_depths_for_the_ref_and_alt_alleles:GQ=Genotype Quality:PL=Normalized_Phred-scaled_likelihoods_for_genotypes=>0_$sampleID,"} else { $sidstring.="IDN:Filter:GT=Genotype:AD=Allelic_depths_for_the_ref_and_alt_alleles:GQ=Genotype Quality:PL=Normalized_Phred-scaled_likelihoods_for_genotypes=>0_$sampleID"} } s/IDN/$sidstring/g; print $_;} next;} if ($_=~/outcolumns/i) { if ($_=~/IDN/) { my $sidstring; for (my $sampleID=5;$sampleID<?.$sampleIDcolumns.q?;$sampleID++) { if ($sampleID<?.$sampleIDcolcond.q?) { $sidstring.="0_$sampleID,"} else { $sidstring.="0_$sampleID"} } s/IDN/$sidstring/g; print $_;} next;} if ($_=~/outheaders/i) { if ($_=~/IDN/) { my $sidstring; for (my $sampleID=5;$sampleID<?.$sampleIDcolumns.q?;$sampleID++) { if ($sampleID<?.$sampleIDcolcond.q?) { $sidstring.="IDN:Filter:GT=Genotype:AD=Allelic_depths_for_the_ref_and_alt_alleles:GQ=Genotype Quality:PL=Normalized_Phred-scaled_likelihoods_for_genotypes,"} else { $sidstring.="IDN:Filter:GT=Genotype:AD=Allelic_depths_for_the_ref_and_alt_alleles:GQ=Genotype Quality:PL=Normalized_Phred-scaled_likelihoods_for_genotypes"} } s/IDN/$sidstring/g; print $_;} next;} elsif ($_=~s/FDN/?.$_[0].q?/g) { if($_=~s/^ODF/?.$regexp_odf.q?/g) {} if($_=~s/ALIGNER/?.$aligner.q?/g) {} if ($_=~/IDN/) { my $sidstring; for (my $sampleID=5;$sampleID<?.$sampleIDcolumns.q?;$sampleID++) { if ($sampleID<?.$sampleIDcolcond.q?) { $sidstring.="$sampleID,"} else { $sidstring.="$sampleID"} } s/IDN/$sidstring/g; print $_;} else { print $_;} } else { if($_=~s/^RD/?.$regexp_rd.q?/g) {} print $_;}' ?.$rd.q?/?.$vmerge_db_template.q? > ?."$odf/$_[0]/$vm_dbf", "\n\n";

    print VMERGE "perl $ids/intersectCollect.pl -db $odf/$_[0]/$vm_dbf -o ", '${outFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_varrecal_$_[2]_filt_annovar_all_variants.txt", "\n\n";
    
    close(VMERGE);   
    FIDSubmitJob(0,$familyid, 1, $_[2],$filename);
    return;
}

sub AnnovarFilter { 
#Filter SNVs by gene, region and databases
#Works on the familyid_allchr_real_recal_resrt_varrecal_filt_annovar.vcf where all variants per subject is present. Prints new files for each analysis
#$_[0] = $familyid NOTE: not sampleid
#$_[1] = aligner
#$_[2] = $gatk_unigt_snps (SNV) or $gatk_unigt_indels (INDEL) or BOTH

    `mkdir -p $odf/$_[0]/$_[1]/info;`; #Creates the alignment folder and info data file directory
    `mkdir -p $odf/$_[0]/$_[1]/GATK;`; #Creates the aligner and annovar folder   
    `mkdir -p $ods/$_[0]/$_[1];`; #Creates the aligner script folder
   
    $filename = "$ods/$_[0]/$_[1]/annovar_filter_$_[2]_$_[0].";
    Checkfnexists($filename, $fnend);

#Info and Logg
    print STDOUT "Creating sbatch script Annovar and writing script file(s) to: ", $filename, "\n";print MASTERL "Creating sbatch script Annovar and writing script file(s) to: ", $filename, "\n";
    print STDOUT "Sbatch script Annovar data files will be written to: ", $odf,"/$_[0]/$_[1]/GATK", "\n";print MASTERL "Sbatch script Annovar data files will be written to: ", $odf,"/$_[0]/$_[1]/GATK", "\n";

    open (ANVARF, ">$filename") or die "Can't write to $filename: $!\n";
    
    print ANVARF "#! /bin/bash -l", "\n";
    print ANVARF "#SBATCH -A ", $aid, "\n";
    print ANVARF "#SBATCH -p node -n $maximum_cores", "\n";
    print ANVARF "#SBATCH -C thin", "\n";	
    print ANVARF "#SBATCH -t 7:00:00", "\n";
    print ANVARF "#SBATCH -J ANNOVARF_$_[2]_$_[0]", "\n";
    print ANVARF "#SBATCH -e $odf/$_[0]/$_[1]/info/annovar_filter_$_[2]_$_[0].", $fnt ,".stderr.txt", "\n";
    print ANVARF "#SBATCH -o $odf/$_[0]/$_[1]/info/annovar_filter_$_[2]_$_[0].", $fnt ,".stdout.txt", "\n";    
    unless ($em eq 0) {
	
	print ANVARF "#SBATCH --mail-type=END", "\n";
	print ANVARF "#SBATCH --mail-type=FAIL", "\n";
	print ANVARF "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print ANVARF 'echo "Running on: $(hostname)"',"\n\n";
    print ANVARF 'inRefDir="', "$annovar_path/humandb", '"', "\n\n"; #$annov path to annovar folder

###
#GATK
###
    print ANVARF "#Directories", "\n";
    print ANVARF 'inFamilyDir="',"$odf/$_[0]/$_[1]/GATK", '"', "\n"; #GATK for now
    print ANVARF 'outFamilyDir="', "$odf/$_[0]/$_[1]/GATK", '"', "\n"; #All subjects
    print ANVARF 'inFamilyPrefix="',"/$_[0]", "_allchr_real_recal_resrt_varrecal_$_[2]_filt.vcf", '"', "\n";
    print ANVARF 'outFamilyPrefix_temp="',"/$_[0]", "_allchr_real_recal_resrt_varrecal_$_[2]_filt_annovar_temp", '"', "\n";
    print ANVARF 'outFamilyPrefix="',"/$_[0]", "_allchr_real_recal_resrt_varrecal_$_[2]_filt_annovar", '"', "\n\n";

#Make sure that mitochondral genome is M and not MT
    print ANVARF "perl -i -p -e 's/^(MT)/M/g' ", '${inFamilyDir}/${inFamilyPrefix}', "\n\n";
#Prepp to annovar format from GATK vcf4	    
    print ANVARF "perl $annovar_path/convert2annovar.pl ",'${inFamilyDir}/${inFamilyPrefix}'," -format vcf4 -i > ", '${outFamilyDir}/${outFamilyPrefix_temp}', "\n\n";

    #Intersect for all samples within familyid and remake file to fit annovar format and subsequent filtering
    print ANVARF q?perl -nae ' if ($_=~/^#/) {print $_;next;} if ($_=~/;set=2/) {} else{ if($F[11] eq "PASS") {} else {$F[11] = "PRES";} print $F[0], "\t", $F[1], "\t", $F[2], "\t", $F[3], "\t", $F[4], "\t", ?;
    #print ANVARF q?perl -nae ' if ($_=~/^#/) {print $_;next;} if ($_=~/;set=2/) {} else{ if($F[6] eq "PASS") {} else {$F[6] = "PRES";} print $F[0], "\t", $F[1], "\t", $F[1], "\t", $F[3], "\t", $F[4], "\t", ?; #Old filtering, before convert2annovar.pl was used.
    my @sample_lex_order = sort @sid; #Use lexiographically sorted sample IDNs since GATK UnifiedGT assigns columns in lexigraphical order. @sid is not lexiographically sorted but taken straight from the command line. This lex sort ensures that if the user did not supply samples in lex order, there will be no sample column swaping. 
    for (my $sampleid=0;$sampleid<scalar(@sample_lex_order);$sampleid++) { #For all sample ids
	
	my $samplecolumn = 14+$sampleid; #First sample genotype starts at col 14 (start 0, perl). NOTE: Important that samples for unifiedGT has same order. Otherwise there will be a sample mix-up.
	
	if ($sampleid eq scalar(@sample_lex_order)-1) {	#Ensure correct order as long as UnifiedGT uses lex sort. 
	    print ANVARF '"',"$sample_lex_order[$sampleid]:", q?$F[11]:$F[?.$samplecolumn.q?]",?;
	}
	else {
	    print ANVARF '"',"$sample_lex_order[$sampleid]:", q?$F[11]:$F[?.$samplecolumn.q?]", "\t", ?;
	}
    }

    print ANVARF q?"\n";}' ?.'${outFamilyDir}/${outFamilyPrefix_temp}'.q? > ?, '${outFamilyDir}/${outFamilyPrefix}', "\n\n"; 
 
    print ANVARF 'inFamilyPrefix="',"/$_[0]", "_allchr_real_recal_resrt_varrecal_$_[2]_filt_annovar", '"', "\n";    
    print ANVARF 'outFamilyPrefix="',"/$_[0]", "_allchr_real_recal_resrt_varrecal_$_[2]_filt_annovar", '"', "\n\n";   	    
    
    my $core_Counter=1;
    for (my $table_names_Counter=0;$table_names_Counter<scalar(@annovar_table_names);$table_names_Counter++) { #For all specified table names
	if ($table_names_Counter eq $core_Counter*$maximum_cores) { #Using only 8 cores
	    
	    print ANVARF "wait", "\n\n";
	    $core_Counter=$core_Counter+1;
	}
	print ANVARF q?perl ?.$annovar_path.q?/annotate_variation.pl ?; 
	print ANVARF q?-?.$annovar_filtering_option{ $annovar_table_names[$table_names_Counter] }.q? ?; #Filtering option
	if ( $annovar_filtering_option{ $annovar_table_names[$table_names_Counter] } eq "geneanno" ) { #Use hgvs output style
	    print ANVARF q?-hgvs ?;
	}
	print ANVARF q?-buildver ?.$annovar_genome_build_version.q? ?;
	if ( $annovar_generic_filtering_option{ $annovar_table_names[$table_names_Counter] } ) { #Handle generic format
	    print ANVARF q?-dbtype generic -genericdbfile ?.$annovar_table_names[$table_names_Counter].q? ?;
	    print ANVARF q?--outfile ${outFamilyDir}${outFamilyPrefix}.?.$annovar_table_names[$table_names_Counter].q? ?;
	}	
	else{
	    print ANVARF q?-dbtype ?.$annovar_table_names[$table_names_Counter].q? ?;
	}
	if ( ($annovar_table_names[$table_names_Counter] =~/^snp/) || ($annovar_table_names[$table_names_Counter] =~/^1000g/) || ($annovar_table_names[$table_names_Counter] =~/_esp/) ) {#Set MAF TH
	    print ANVARF q?--maf_threshold ?.$annovar_maf_threshold.q? ?;
	}
	if ( $annovar_table_names[$table_names_Counter] =~/^avsift/ ) {#Set sift score TH
	    print ANVARF q?--sift_threshold ?.$annovar_sift_threshold.q? ?;
	}
	print ANVARF q?${inFamilyDir}${inFamilyPrefix} ?; #Infile. Outfile is named using infile prefix except for generic files 
	print ANVARF q?${inRefDir} &?, "\n\n"; #annovar/humandb
    }
    print ANVARF "wait", "\n\n";
    
    print ANVARF "rm ", '${outFamilyDir}/${outFamilyPrefix_temp}', "\n"; #Remove temp file
    close(ANVARF);
    FIDSubmitJob(0,$familyid, 1, $_[2],$filename);
    return;
}

sub GATK_combinevariants { 
#GATK CombineVariants
#$_[0]= sampleid
#$_[1]= $aligner, to choose the correct dir depending on what aligner has been used previously
#$_[2] = BOTH
#$_[3] = $familyid

    `mkdir -p $odf/$_[3]/$_[1]/info;`; #Creates the aligner folder and info data file directory
    `mkdir -p $odf/$_[3]/$_[1]/GATK`; #Creates the aligner folder, GATK data file directory
    `mkdir -p $ods/$_[3]/$_[1]`; #Creates the aligner folder script file directory
    $filename = "$ods/$_[3]/$_[1]/gatk_combinevar_$_[2]_$_[3].";   

    Checkfnexists($filename, $fnend);

#Info and Logg
    print STDOUT "Creating sbatch script GATK CombineVariants and writing script file(s) to: ", $filename, "\n";print MASTERL "Creating sbatch script CombineVariants and writing script file(s) to: ", $filename, "\n";
    print STDOUT "Sbatch script GATK CombineVariants data files will be written to: ", $odf,"/$_[0]/$_[1]/GATK", "\n";print MASTERL "Sbatch script GATK CombineVariants data files will be written to: ", $odf,"/$_[0]/$_[1]/GATK", "\n";

    open (GATK_COMBVAR, ">$filename") or die "Can't write to $filename: $!\n";
    
    print GATK_COMBVAR "#! /bin/bash -l", "\n";
    print GATK_COMBVAR "#SBATCH -A ", $aid, "\n";
    print GATK_COMBVAR "#SBATCH -p node -n 1", "\n";
    print GATK_COMBVAR "#SBATCH -C thin", "\n";	
    print GATK_COMBVAR "#SBATCH -t 2:00:00", "\n";
    print GATK_COMBVAR "#SBATCH -J GATK_CoVa_$_[2]_", $_[3], "\n";
    print GATK_COMBVAR "#SBATCH -e $odf/$_[3]/$_[1]/info/gatk_combinevar_$_[2]_$_[3].", $fnt ,".stderr.txt", "\n";
    print GATK_COMBVAR "#SBATCH -o $odf/$_[3]/$_[1]/info/gatk_combinevar_$_[2]_$_[3].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print GATK_COMBVAR "#SBATCH --mail-type=END", "\n";
	print GATK_COMBVAR "#SBATCH --mail-type=FAIL", "\n";
	print GATK_COMBVAR "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print GATK_COMBVAR 'echo "Running on: $(hostname)"',"\n\n";
    print GATK_COMBVAR "#Reference Archive", "\n";
    print GATK_COMBVAR 'referenceArchive="', "$rd", '"', "\n\n"; 
    print GATK_COMBVAR "#Samples", "\n";
    print GATK_COMBVAR 'inSampleDir="',"$odf/$_[3]/$_[1]/GATK", '"', "\n";
    print GATK_COMBVAR 'outSampleDir="', "$odf/$_[3]/$_[1]/GATK", '"', "\n\n"; 
    
    print GATK_COMBVAR "#GATK CombineVariants","\n\n";
    	   
    print GATK_COMBVAR "java -Xmx2g -jar $gatk_path/GenomeAnalysisTK.jar -l INFO -T CombineVariants -R ", '${referenceArchive}',"/$genomeref -V:$_[3]_SNV ", '${inSampleDir}', "/$_[3]", "_allchr_real_recal_resrt_varrecal_SNV_filt.vcf -V:$_[3]_INDEL ", '${inSampleDir}', "/$_[3]", "_allchr_real_recal_resrt_varrecal_INDEL_filt.vcf -o ",'${outSampleDir}', "/$_[3]", "_allchr_real_recal_resrt_varrecal_BOTH_filt.vcf", "\n\n";

    print GATK_COMBVAR "wait", "\n\n";
    close(GATK_COMBVAR);   
    FIDSubmitJob(0,$familyid, 1, $_[2],$filename);
    return;
}

sub GATK_varianteval { 
#GATK VariantEval
#$_[0]= sampleid
#$_[1]= $aligner, to choose the correct dir depending on what aligner has been used previously
#$_[2] = $gatk_unigt_snps (SNV) or $gatk_unigt_indels (INDEL), (BOTH)
#$_[3] = EXOME
#$_[4] = $familyid 

    `mkdir -p $odf/$_[0]/$_[1]/info;`; #Creates the aligner folder and info data file directory
    `mkdir -p $odf/$_[0]/$_[1]/GATK/varianteval`; #Creates the aligner folder, GATK VariantEval data file directory
    `mkdir -p $ods/$_[0]/$_[1]`; #Creates the aligner folder script file directory
    if ($_[3]) { $filename = "$ods/$_[0]/$_[1]/gatk_vareval_exome_$_[2]_$_[0].";   
    }
    else { $filename = "$ods/$_[0]/$_[1]/gatk_vareval_$_[2]_$_[0].";
    }
    Checkfnexists($filename, $fnend);

#Info and Logg
    print STDOUT "Creating sbatch script GATK VariantEval and writing script file(s) to: ", $filename, "\n";print MASTERL "Creating sbatch script GATK VariantEval and writing script file(s) to: ", $filename, "\n";
    print STDOUT "Sbatch script GATK VariantEval data files will be written to: ", $odf,"/$_[0]/$_[1]/GATK/varianteval", "\n";print MASTERL "Sbatch script GATK VariantEval data files will be written to: ", $odf,"/$_[0]/$_[1]/GATK/varianteval", "\n";

    open (GATK_VAREVAL, ">$filename") or die "Can't write to $filename: $!\n";
    
    print GATK_VAREVAL "#! /bin/bash -l", "\n";
    print GATK_VAREVAL "#SBATCH -A ", $aid, "\n";
    print GATK_VAREVAL "#SBATCH -p node -n 1", "\n";
    print GATK_VAREVAL "#SBATCH -C thin", "\n";	
    print GATK_VAREVAL "#SBATCH -t 2:00:00", "\n";
    if ($_[3]) {
	print GATK_VAREVAL "#SBATCH -J GATK_VaEv_ex_$_[2]_", $_[4], "\n";
	print GATK_VAREVAL "#SBATCH -e $odf/$_[0]/$_[1]/info/gatk_vareval_exome_$_[2]_$_[0].", $fnt ,".stderr.txt", "\n";
	print GATK_VAREVAL "#SBATCH -o $odf/$_[0]/$_[1]/info/gatk_vareval_exome_$_[2]_$_[0].", $fnt ,".stdout.txt", "\n";
    }
    else {
	print GATK_VAREVAL "#SBATCH -J GATK_VarEval_$_[2]_", $_[0], "\n";
	print GATK_VAREVAL "#SBATCH -e $odf/$_[0]/$_[1]/info/gatk_vareval_$_[2]_$_[0].", $fnt ,".stderr.txt", "\n";
	print GATK_VAREVAL "#SBATCH -o $odf/$_[0]/$_[1]/info/gatk_vareval_$_[2]_$_[0].", $fnt ,".stdout.txt", "\n";
    }
    unless ($em eq 0) {
	
	print GATK_VAREVAL "#SBATCH --mail-type=END", "\n";
	print GATK_VAREVAL "#SBATCH --mail-type=FAIL", "\n";
	print GATK_VAREVAL "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print GATK_VAREVAL 'echo "Running on: $(hostname)"',"\n\n";
    print GATK_VAREVAL "#Reference Archive", "\n";
    print GATK_VAREVAL 'referenceArchive="', "$rd", '"', "\n\n"; 
    print GATK_VAREVAL "#Samples", "\n";
    print GATK_VAREVAL 'inSampleDir="',"$odf/$_[0]/$_[1]/GATK", '"', "\n"; #SampleID
    print GATK_VAREVAL 'outSampleDir="', "$odf/$_[0]/$_[1]/GATK/varianteval", '"', "\n\n"; #SampleID
    print GATK_VAREVAL 'inSampleDir_2="',"$odf/$_[4]/$_[1]/GATK", '"', "\n"; #FamilyID

    my $tempinfile = $avcovfn{$_[0]};
   
#Needed for both all_variants and only exome variants 
    unless ($_[3]) { #Create sampleID from familyID_annovar_all_variants but only once, $_[3] is 0 for first gatk_vareval call
#GATK SelectVariants
#Select SampleID from familyID_annovar_all_variants.
	print GATK_VAREVAL "\n\n#GATK SelectVariants","\n\n";
	print GATK_VAREVAL "java -Xmx2g -jar $gatk_path/GenomeAnalysisTK.jar -l INFO -T SelectVariants -R ", '${referenceArchive}',"/$genomeref -V: ", '${inSampleDir_2}', "/$_[4]", "_allchr_real_recal_resrt_varrecal_$_[2]_filt.vcf -o ", '${inSampleDir}', "/$tempinfile", "_allchr_real_recal_resrt_varrecal_$_[2]_filt.vcf -sn $_[0]", "\n\n";
	
#Remove chr
	print GATK_VAREVAL "perl -i -p -e '", 's/^chr(.+)/$1/g',"' ", '${inSampleDir}', "/$tempinfile", "_allchr_real_recal_resrt_varrecal_$_[2]_filt.vcf", "\n\n";
    }
    
    if ($_[3]) { #Exome
###
#Prepp infile
###

#Find exonic variants (familyID_annovar_all_variants)
	my $variantType=""; #To translate BOTH to INDEL since no ..BOTH_filt_annovar_all_variants.txt is created only INDEL or SNV
	my $hapcal_both_file = "$odf/$_[4]/$_[1]/GATK/$_[4]_allchr_real_recal_resrt_varrecal_$_[2]_filt_annovar_all_variants.txt";
	if ($_[2] eq "BOTH") { 
	    if ($pGATK_HAPCAL eq 1) {
		print GATK_VAREVAL "grep exon ",'${inSampleDir_2}', "/$_[4]", "_allchr_real_recal_resrt_varrecal_$_[2]_filt_annovar_all_variants.txt > ", '${inSampleDir}', "/slask_exonic_variants.txt", "\n\n";
	    }
	    elsif (-e $hapcal_both_file) {
		print GATK_VAREVAL "grep exon ",'${inSampleDir_2}', "/$_[4]", "_allchr_real_recal_resrt_varrecal_$_[2]_filt_annovar_all_variants.txt > ", '${inSampleDir}', "/slask_exonic_variants.txt", "\n\n";
	    }
	    else {
		$variantType = "SNV"; #Collect SNVs
		print GATK_VAREVAL "grep exon ",'${inSampleDir_2}', "/$_[4]", "_allchr_real_recal_resrt_varrecal_"."$variantType"."_filt_annovar_all_variants.txt > ", '${inSampleDir}', "/slask_exonic_variants.txt", "\n\n";
		$variantType = "INDEL"; #Collect INDELs and write to the same file as above to enable concatenated analysis
		print GATK_VAREVAL "grep exon ",'${inSampleDir_2}', "/$_[4]", "_allchr_real_recal_resrt_varrecal_"."$variantType"."_filt_annovar_all_variants.txt >> ", '${inSampleDir}', "/slask_exonic_variants.txt", "\n\n"; #Note >>
	    }
	    
	}
	elsif ($_[2] eq "INDEL") { 
	    $variantType = "INDEL";
	    print GATK_VAREVAL "grep exon ",'${inSampleDir_2}', "/$_[4]", "_allchr_real_recal_resrt_varrecal_"."$variantType"."_filt_annovar_all_variants.txt > ", '${inSampleDir}', "/slask_exonic_variants.txt", "\n\n";	
	}
	else { 
	    $variantType = "SNV"; 
	    print GATK_VAREVAL "grep exon ",'${inSampleDir_2}', "/$_[4]", "_allchr_real_recal_resrt_varrecal_"."$variantType"."_filt_annovar_all_variants.txt > ", '${inSampleDir}', "/slask_exonic_variants.txt", "\n\n";
	}

#Remove chr
	print GATK_VAREVAL "perl -i -p -e '", 's/^chr(.+)/$1/g',"' ", '${inSampleDir}', "/slask_exonic_variants.txt", "\n\n";
	
#Intersect exonic variants from created sampleID vcf file (required for GATK_Vareval)
	print GATK_VAREVAL "intersectBed -a ",'${inSampleDir}', "/$tempinfile", "_allchr_real_recal_resrt_varrecal_$_[2]_filt.vcf -b ",'${inSampleDir}', "/slask_exonic_variants.txt > ",'${inSampleDir}', "/slask_exonic_variants.vcf", "\n\n";

#Supply VCF header from original sampleID vcf file (required for GATK_Vareval)
	print GATK_VAREVAL "grep ",'"',"#",'" ', '${inSampleDir}', "/$tempinfile", "_allchr_real_recal_resrt_varrecal_$_[2]_filt.vcf > ", '${inSampleDir}', "/slask_exonic_variants_header.vcf ", "\n\n";

#Combine VCF header and variants (sampleID)
	print GATK_VAREVAL "cat ", '${inSampleDir}', "/slask_exonic_variants_header.vcf ", '${inSampleDir}', "/slask_exonic_variants.vcf > ", '${inSampleDir}', "/$tempinfile", "_allchr_real_recal_resrt_varrecal_$_[2]_filt.vcf_exonic", "\n\n";

###
#VariantEval (Exonic variants)
###
	print GATK_VAREVAL "#GATK VariantEval","\n\n";
	print GATK_VAREVAL "java -Xmx2g -jar $gatk_path/GenomeAnalysisTK.jar -l INFO -T VariantEval -R ", '${referenceArchive}',"/$genomeref -D ", '${referenceArchive}',"/dbsnp_132.hg19.excluding_sites_after_129_nochr.vcf -gold ",'${referenceArchive}',"/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -o ",'${outSampleDir}', "/$tempinfile", "_allchr_real_recal_resrt_varrecal_$_[2]_filt.vcf.varianteval_exonic --eval ", '${inSampleDir}', "/$tempinfile", "_allchr_real_recal_resrt_varrecal_$_[2]_filt.vcf_exonic", "\n\n";
	
#Clean up
	print GATK_VAREVAL "rm ",'${inSampleDir}', "/slask_exonic_variants_header.vcf", "\n\n";
	print GATK_VAREVAL "rm ",'${inSampleDir}', "/slask_exonic_variants.vcf", "\n\n";
	print GATK_VAREVAL "rm ",'${inSampleDir}', "/slask_exonic_variants.txt", "\n\n";
    }
    else {
###
#VariantEval (all variants)
###
	print GATK_VAREVAL "#GATK VariantEval","\n\n";
	print GATK_VAREVAL "java -Xmx2g -jar $gatk_path/GenomeAnalysisTK.jar -l INFO -T VariantEval -R ", '${referenceArchive}',"/$genomeref -D ", '${referenceArchive}',"/dbsnp_132.hg19.excluding_sites_after_129_nochr.vcf -gold ",'${referenceArchive}',"/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -o ",'${outSampleDir}', "/$tempinfile", "_allchr_real_recal_resrt_varrecal_$_[2]_filt.vcf.varianteval --eval ", '${inSampleDir}', "/$tempinfile", "_allchr_real_recal_resrt_varrecal_$_[2]_filt.vcf", "\n\n";
    }
    close(GATK_VAREVAL);   
    FIDSubmitJob(0,$familyid, 2, $_[2],$filename); #To not add jobIDs to later jobID{chainkey}
    return;
}

sub Pindel { 
#Pindel (SV detection). Runs PINDEL on all samples within family, converts _D and _SI to vcf and then merges them using CombineVariants. These are then used as input to GATK_UnifiedGT with know allel enabled analysing only the Pindel produced sites. Last stage merges GATK_INDEL and PINDEL_GATK_INDELS into one vcf. 
#$_[0] = $familyid
#$_[1] = $aligner, to choose the correct dir depending on what aligner has been used previously

    `mkdir -p $odf/$_[0]/$_[1]/info;`; #Creates the aligner folder and info data file directory
    `mkdir -p $odf/$_[0]/$_[1]/Pindel`; #Creates the aligner folder, Pindel data file directory
    `mkdir -p $ods/$_[0]/$_[1]/Pindel`; #Creates the aligner folder script file directory
    $filename = "$ods/$_[0]/$_[1]/pindel_$_[0]."; 

    Checkfnexists($filename, $fnend);

#Info and Logg
    print STDOUT "Creating sbatch script Pindle and writing script file(s) to: ", $filename, "\n";print MASTERL "Creating sbatch script Pindel and writing script file(s) to: ", $filename, "\n";
    print STDOUT "Sbatch script Pindel data files will be written to: ", $odf,"/$_[0]/$_[1]/Pindel", "\n";print MASTERL "Sbatch script Pindel data files will be written to: ", $odf,"/$_[0]/$_[1]/Pindel", "\n";

    open (PINDEL, ">$filename") or die "Can't write to $filename: $!\n";
    
    print PINDEL "#! /bin/bash -l", "\n";
    print PINDEL "#SBATCH -A ", $aid, "\n";
    print PINDEL "#SBATCH -p node -n $maximum_cores", "\n";
    print PINDEL "#SBATCH -C thin", "\n";	
    print PINDEL "#SBATCH -t 20:00:00", "\n";

    print PINDEL "#SBATCH -J PINDEL_$_[0]_", "\n";
    print PINDEL "#SBATCH -e $odf/$_[0]/$_[1]/info/pindel_$_[0].", $fnt ,".stderr.txt", "\n";
    print PINDEL "#SBATCH -o $odf/$_[0]/$_[1]/info/pindel_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print PINDEL "#SBATCH --mail-type=END", "\n";
	print PINDEL "#SBATCH --mail-type=FAIL", "\n";
	print PINDEL "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print PINDEL 'echo "Running on: $(hostname)"',"\n\n";
    print PINDEL "#Reference Archive", "\n";
    print PINDEL 'referenceArchive="', "$rd", '"', "\n\n"; 
    print PINDEL "#Directories", "\n";
    print PINDEL 'inSampleDir="', "$odf/$_[0]/$_[1]/Pindel", '"', "\n";
    print PINDEL 'inSampleDir2="', "$odf/$_[0]/$_[1]/GATK", '"', "\n";
    print PINDEL 'outSampleDir="', "$odf/$_[0]/$_[1]/Pindel", '"', "\n";
    print PINDEL 'outSampleDir2="', "$odf/$_[0]/$_[1]/GATK", '"', "\n\n"; 

    print STDOUT "#Create Pindel bam config File per family","\n";
    my $pindel_config_file = "$odf/$_[0]/$_[1]/Pindel/pindel_bam_config_fam$_[0].txt";

    open (PINDEL_INF, ">$pindel_config_file") or die "Can't write to $pindel_config_file: $!\n";
    print STDOUT "Written in $pindel_config_file:\n";
    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) { #For all sample ids
	
	my $tempinfile = $avcovfn{$sid[$sampleid]};
	
	print PINDEL_INF "$odf/$sid[$sampleid]/$_[1]/GATK/$tempinfile","_allchr_real_recal_resrt.bam 300 $sid[$sampleid]", "\n"; #file insert size and SampleID
	print STDOUT "$odf/$sid[$sampleid]/$_[1]/GATK/$tempinfile","_allchr_real_recal_resrt.bam 300 $sid[$sampleid]", "\n";
    }    	    
    close(PINDEL_INF); #Closes Pindel config file

    print PINDEL "#Pindel\n"; #Call INDELS using PINDEL
    print PINDEL q?pindel -f ${referenceArchive}/?.$genomeref.q? -i ?.$pindel_config_file.q? -c ALL -o ${outSampleDir}/?.$_[0].q?_allchr_real_recal_resrt_pindel -T $maximum_cores?, "\n\n"; 

    print PINDEL "#Pindel2vcf\n"; #Convert PINDEL calls to vcf
    print PINDEL q?pindel2vcf -p ${inSampleDir}/?.$_[0].q?_allchr_real_recal_resrt_pindel_D -r ${referenceArchive}/?.$genomeref.q? -R NA -d NA -v ${outSampleDir}/?.$_[0].q?_allchr_real_recal_resrt_pindel_D.vcf &?,"\n\n";
    print PINDEL q?pindel2vcf -p ${inSampleDir}/?.$_[0].q?_allchr_real_recal_resrt_pindel_SI -r ${referenceArchive}/?.$genomeref.q? -R NA -d NA -v ${outSampleDir}/?.$_[0].q?_allchr_real_recal_resrt_pindel_SI.vcf &?,"\n\n";
    print PINDEL "wait", "\n"; 

    print PINDEL "\n#GATK CombineVariants","\n\n"; #Combine PINDEL deletions (D) and short insertions (SI) to 1 file
    print PINDEL q?java -Xmx2g -jar ?.$gatk_path.q?/GenomeAnalysisTK.jar -l INFO -T CombineVariants -R ${referenceArchive}/?.$genomeref.q? -V: ${inSampleDir}/?.$_[0].q?_allchr_real_recal_resrt_pindel_D.vcf -V: ${inSampleDir}/?.$_[0].q?_allchr_real_recal_resrt_pindel_SI.vcf -o ${outSampleDir}/?.$_[0].q?_allchr_real_recal_resrt_pindel_INDEL.vcf?, "\n";
    
    print PINDEL "\n#GATK UnifiedGenoTyper\n"; #Tap into GATK UnifiedGT using PINDEL split-read calls (skips GATKUnifedGT discovery phase)
    print PINDEL q?java -Xmx12g -jar ?.$gatk_path.q?/GenomeAnalysisTK.jar -l INFO -T UnifiedGenotyper -R ${referenceArchive}/?.$genomeref.q? -D ${referenceArchive}/?.$gatk_recal_knset.q? -glm INDEL -nt $maximum_cores -stand_call_conf 30.0 -stand_emit_conf 30.0 --min_base_quality_score 20 --annotation QualByDepth --annotation HaplotypeScore --annotation MappingQualityRankSumTest --annotation ReadPosRankSumTest --annotation FisherStrand --annotation RMSMappingQuality --annotation DepthOfCoverage -L ${referenceArchive}/?.$gatk_bait;
	
    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) { #For all sample ids
	
	my $tempinfile = $avcovfn{$sid[$sampleid]};
	
	print PINDEL " -I $odf/$sid[$sampleid]/$_[1]/GATK/$tempinfile","_allchr_real_recal_resrt.bam"; #Supply original .bam files per SampleID
    }    	    
    
    print PINDEL q? -o ${outSampleDir}/?.$_[0].q?_allchr_real_recal_resrt_pindel_GATK_INDEL.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -alleles ${inSampleDir}/?.$_[0].q?_allchr_real_recal_resrt_pindel_INDEL.vcf?,"\n";

    #if ( -e "$odf/$_[0]/$_[1]/GATK/$_[0]_allchr_real_recal_resrt_varrecal_INDEL_filt.vcf" ) {
#	print PINDEL "\n#GATK CombineVariants","\n\n"; #Combine all INDELs (original GATK calls and PINDEL+GATK calls) to one file 
#	print PINDEL q?java -Xmx2g -jar ?.$gatk_path.q?/GenomeAnalysisTK.jar -l INFO -T CombineVariants -R ${referenceArchive}/?.$genomeref.q? -V: ${inSampleDir}/?.$_[0].q?_allchr_real_recal_resrt_pindel_GATK_INDEL.vcf -V: ${inSampleDir2}/?.$_[0].q?_allchr_real_recal_resrt_varrecal_INDEL_filt.vcf -o ${inSampleDir}/?.$_[0].q?_allchr_real_recal_resrt_pindel_GATKx2_INDEL.vcf?, "\n";
 #   }
    print PINDEL "\n\nwait", "\n\n";    
    close(PINDEL);   
    FIDSubmitJob(0,$familyid, 2, $_[2],$filename);
    return;
}

sub GATK_varrecal { 
#GATK VariantRecalibrator/ApplyRecalibration
#$_[0]= $familyid NOTE: not sampleid
#$_[1]= $aligner, to choose the correct dir depending on what aligner has been used previously
#$_[2] = $gatk_unigt_snps (SNV) or $gatk_unigt_indels (INDEL) or BOTH

    `mkdir -p $odf/$_[0]/$_[1]/info;`; #Creates the aligner folder and info data file directory
    `mkdir -p $odf/$_[0]/$_[1]/GATK/intermediary`; #Creates the aligner folder, GATK data file directory
    `mkdir -p $odf/$_[0]/$_[1]/GATK`; #Creates the aligner folder, GATK data file directory
    `mkdir -p $ods/$_[0]/$_[1]`; #Creates the aligner folder script file directory
    $filename = "$ods/$_[0]/$_[1]/gatk_varrecal_$_[2]_$_[0].";   

    Checkfnexists($filename, $fnend);

#Info and Logg
    print STDOUT "Creating sbatch script GATK VariantRecalibrator/ApplyRecalibration and writing script file(s) to: ", $filename, "\n";print MASTERL "Creating sbatch script GATK VariantRecalibrator/ApplyRecalibration and writing script file(s) to: ", $filename, "\n";
    print STDOUT "Sbatch script GATK VariantRecalibrator/ApplyRecalibration data files will be written to: ", $odf,"/$_[0]/$_[1]/GATK", "\n";print MASTERL "Sbatch script GATK VariantRecalibrator/ApplyRecalibration data files will be written to: ", $odf,"/$_[0]/$_[1]/GATK", "\n";

    open (GATK_VARREC, ">$filename") or die "Can't write to $filename: $!\n";
    
    print GATK_VARREC "#! /bin/bash -l", "\n";
    print GATK_VARREC "#SBATCH -A ", $aid, "\n";
    print GATK_VARREC "#SBATCH -p node -n $maximum_cores", "\n";
    print GATK_VARREC "#SBATCH -C thin", "\n";	
    print GATK_VARREC "#SBATCH -t 10:00:00", "\n";
    print GATK_VARREC "#SBATCH -J GATK_VarReCal_$_[2]_", $_[0], "\n";
    print GATK_VARREC "#SBATCH -e $odf/$_[0]/$_[1]/info/gatk_varrecal_$_[2]_$_[0].", $fnt ,".stderr.txt", "\n";
    print GATK_VARREC "#SBATCH -o $odf/$_[0]/$_[1]/info/gatk_varrecal_$_[2]_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print GATK_VARREC "#SBATCH --mail-type=END", "\n";
	print GATK_VARREC "#SBATCH --mail-type=FAIL", "\n";
	print GATK_VARREC "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print GATK_VARREC 'echo "Running on: $(hostname)"',"\n\n";
    print GATK_VARREC "#Reference Archive", "\n";
    print GATK_VARREC 'referenceArchive="', "$rd", '"', "\n\n"; 
    print GATK_VARREC "#Samples", "\n";
    print GATK_VARREC 'inFamilyDir="',"$odf/$_[0]/$_[1]/GATK", '"', "\n";
    print GATK_VARREC 'inFamilyDir2="', "$odf/$_[0]/$_[1]/GATK/HapCall", '"', "\n\n";
    print GATK_VARREC 'outFamilyDir="', "$odf/$_[0]/$_[1]/GATK/intermediary", '"', "\n\n";
    print GATK_VARREC 'outFamilyDir3="',"$odf/$_[0]", '"', "\n\n"; #For .fam file 
 
    unless (-e "$odf/$_[0]/$familyid.fam") {
	print GATK_VARREC "#Generating '.fam' file for GATK HAplotypeCaller","\n\n";
	print GATK_VARREC q?perl -nae 'my %sample_info;my $mother;my $father; while (<>) { my @F = split(/\t/,$_); if ($_!~/^#/) { if( $F[0]=~/(\d+)-(\d+|-\d+)-(\d+)(A|U)/) {} if ($F[3] == 1) {$father = $F[0];} if ($F[2] == 1) {$mother = $F[0];} if($3 % 2 == 1) {push (@{ $sample_info{$1}{$F[0]} }, "1");} else {push (@{ $sample_info{$1}{$F[0]} }, "2");} if ($4 eq "A") {push (@{ $sample_info{$1}{$F[0]} }, "2");} else {push (@{ $sample_info{$1}{$F[0]} }, "1");} } } for my $familyid (keys %sample_info ) { for my $sampleid (keys %{ $sample_info{$familyid} }) {print $familyid, " ", $sampleid, " ", $father, " ", $mother," "; for (my $i=0;$i<scalar(@{ $sample_info{$familyid}{$sampleid} } );$i++) {print $sample_info{$familyid}{$sampleid}[$i], " ";}print "\n"; } } last;' ?.$pedigree.q? > ${outFamilyDir3}/?.$familyid.q?.fam?, "\n\n";
    }
    if ($_[2] eq "BOTH") {#Requires that the input sample has been produced with HaplotypeCaller   
	if ($wgs == 0) { #Exome analysis 
###
#GATK CombineVariants
###
#Needed to include reference exomes to power the building of the probabalistic model. Variants unique to these exomes will be filtered out after varrecal and applyrecal.
	    print GATK_VARREC "\n#GATK CombineVariants","\n\n";
	    print GATK_VARREC "java -Xmx4g -jar $gatk_path/GenomeAnalysisTK.jar -l INFO -T CombineVariants -R ", '${referenceArchive}',"/$genomeref -V: ", '${inFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_raw_$_[2].vcf -V: ", '${referenceArchive}',"/$gatk_exref_snp -o ", '${inFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_raw_$_[2]_comb_ref.vcf";
	    
	    print GATK_VARREC "\n\n#GATK VariantRecalibrator","\n\n";	
	    print GATK_VARREC "java -Xmx12g -jar $gatk_path/GenomeAnalysisTK.jar -l INFO -T VariantRecalibrator -R ", '${referenceArchive}',"/$genomeref -recalFile ",'${outFamilyDir}', "/$_[0]", "_allchr_varrecal_$_[2]_comb_ref.intervals -rscriptFile ", '${outFamilyDir}', "/$_[0]", "_allchr_varrecal_$_[2]_comb_ref.intervals.plots.R -tranchesFile ",'${outFamilyDir}', "/$_[0]", "_allchr_varrecal_$_[2]_comb_ref.intervals.tranches -resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 ", '${referenceArchive}',"/hapmap_3.3.b37.sites.vcf -resource:omni,VCF,known=false,training=true,truth=false,prior=12.0 ",'${referenceArchive}', "/1000G_omni2.5.b37.sites.vcf -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=8.0 ",'${referenceArchive}', "/dbsnp_135.b37.vcf -resource:mills,VCF,known=true,training=true,truth=true,prior=12.0 ", '${referenceArchive}',"/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ --mode $_[2] -nt $maximum_cores -input ", '${inFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_raw_$_[2]_comb_ref.vcf -L ", '${referenceArchive}',"/$gatk_bait ";
	    if (scalar(@sid) > 2) {
		for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) { #For all sample ids
		    if ( ($sid[$sampleid] =~ /(\d+)-(\d+|-\d+)-(\d+)(A|U)/) ) {#Match sampleID
			if ($2 eq 2) { #Parent
			    print GATK_VARREC "--pedigreeValidationType SILENT --pedigree ", '${outFamilyDir3}',"/$familyid.fam ";
			    last; #Only print once if a parent is found (required to include pedigree)
			}
		    }
		}	
	    }
	}
	else { #WGS Analysis
	    print GATK_VARREC "\n\n#GATK VariantRecalibrator","\n\n";	
	    print GATK_VARREC "java -Xmx12g -jar $gatk_path/GenomeAnalysisTK.jar -l INFO -T VariantRecalibrator -R ", '${referenceArchive}',"/$genomeref -recalFile ",'${outFamilyDir}', "/$_[0]", "_allchr_varrecal_$_[2].intervals -rscriptFile ", '${outFamilyDir}', "/$_[0]", "_allchr_varrecal_$_[2].intervals.plots.R -tranchesFile ",'${outFamilyDir}', "/$_[0]", "_allchr_varrecal_$_[2].intervals.tranches -resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 ", '${referenceArchive}',"/hapmap_3.3.b37.sites.vcf -resource:omni,VCF,known=false,training=true,truth=false,prior=12.0 ",'${referenceArchive}', "/1000G_omni2.5.b37.sites.vcf -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=8.0 ",'${referenceArchive}', "/dbsnp_135.b37.vcf -resource:mills,VCF,known=true,training=true,truth=true,prior=12.0 ", '${referenceArchive}',"/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ --mode $_[2] -nt $maximum_cores -input ", '${inFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_raw_$_[2].vcf ";
	    if (scalar(@sid) > 2) {
		for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) { #For all sample ids
		    if ( ($sid[$sampleid] =~ /(\d+)-(\d+|-\d+)-(\d+)(A|U)/) ) {#Match sampleID
			if ($2 eq 2) { #Parent
			    print GATK_VARREC "--pedigreeValidationType SILENT --pedigree ", '${outFamilyDir3}',"/$familyid.fam ";
			    last; #Only print once if a parent is found (required to include pedigree)
			}
		    }
		}		
	    }
	}
    }
    elsif ($_[2] eq "SNV") {
	
	if ($wgs == 0) { #Exome analysis
###
#GATK CombineVariants
###
#Needed to include reference exomes to power the building of the probabalistic model. Variants unique to these exomes will be filtered out after varrecal and applyrecal.
	    print GATK_VARREC "\n#GATK CombineVariants","\n\n";
	    print GATK_VARREC "java -Xmx4g -jar $gatk_path/GenomeAnalysisTK.jar -l INFO -T CombineVariants -R ", '${referenceArchive}',"/$genomeref -V: ", '${inFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_raw_$_[2].vcf -V: ", '${referenceArchive}',"/$gatk_exref_snp -o ", '${inFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_raw_$_[2]_comb_ref.vcf";
	    
	    print GATK_VARREC "\n\n#GATK VariantRecalibrator","\n\n";	
	    print GATK_VARREC "java -Xmx12g -jar $gatk_path/GenomeAnalysisTK.jar -l INFO -T VariantRecalibrator -R ", '${referenceArchive}',"/$genomeref -recalFile ",'${outFamilyDir}', "/$_[0]", "_allchr_varrecal_$_[2]_comb_ref.intervals -rscriptFile ", '${outFamilyDir}', "/$_[0]", "_allchr_varrecal_$_[2]_comb_ref.intervals.plots.R -tranchesFile ",'${outFamilyDir}', "/$_[0]", "_allchr_varrecal_$_[2]_comb_ref.intervals.tranches -resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 ", '${referenceArchive}',"/hapmap_3.3.b37.sites.vcf -resource:omni,VCF,known=false,training=true,truth=false,prior=12.0 ",'${referenceArchive}', "/1000G_omni2.5.b37.sites.vcf -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=8.0 ",'${referenceArchive}', "/dbsnp_135.b37.vcf -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ --mode SNP -nt $maximum_cores -input ", '${inFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_raw_$_[2]_comb_ref.vcf -L ", '${referenceArchive}',"/$gatk_bait ";
	    if (scalar(@sid) > 2) {
		for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) { #For all sample ids
		    if ( ($sid[$sampleid] =~ /(\d+)-(\d+|-\d+)-(\d+)(A|U)/) ) {#Match sampleID
			if ($2 eq 2) { #Parent
			    print GATK_VARREC "--pedigreeValidationType SILENT --pedigree ", '${outFamilyDir3}',"/$familyid.fam ";
			    last; #Only print once if a parent is found (required to include pedigree)
			}
		    }
		}
	    }
	}
	else { #WGS Analysis
	    print GATK_VARREC "\n\n#GATK VariantRecalibrator","\n\n";	
	    print GATK_VARREC "java -Xmx12g -jar $gatk_path/GenomeAnalysisTK.jar -l INFO -T VariantRecalibrator -R ", '${referenceArchive}',"/$genomeref -recalFile ",'${outFamilyDir}', "/$_[0]", "_allchr_varrecal_$_[2].intervals -rscriptFile ", '${outFamilyDir}', "/$_[0]", "_allchr_varrecal_$_[2].intervals.plots.R -tranchesFile ",'${outFamilyDir}', "/$_[0]", "_allchr_varrecal_$_[2].intervals.tranches -resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 ", '${referenceArchive}',"/hapmap_3.3.b37.sites.vcf -resource:omni,VCF,known=false,training=true,truth=false,prior=12.0 ",'${referenceArchive}', "/1000G_omni2.5.b37.sites.vcf -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=8.0 ",'${referenceArchive}', "/dbsnp_135.b37.vcf -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ --mode SNP -nt $maximum_cores -input ", '${inFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_raw_$_[2].vcf ";
	    if (scalar(@sid) > 2) {
		for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) { #For all sample ids
		    if ( ($sid[$sampleid] =~ /(\d+)-(\d+|-\d+)-(\d+)(A|U)/) ) {#Match sampleID
			if ($2 eq 2) { #Parent
			    print GATK_VARREC "--pedigreeValidationType SILENT --pedigree ", '${outFamilyDir3}',"/$familyid.fam ";
			    last; #Only print once if a parent is found (required to include pedigree)
			}
		    }
		}
	    }
	}
    }
    else { #INDELs
	if ($wgs == 0) { #Exome analysis - Restrict analysis to padded target file
#No need to combine variants since probabalistic model is not used	    
	    
	    print GATK_VARREC "\n\n#GATK VariantFiltration","\n\n"; #Unable to build bayesian model for exome indels due to too few indels calls, therefore hard filtering is used. 
	    print GATK_VARREC "java -Xmx2g -jar $gatk_path/GenomeAnalysisTK.jar -l INFO -T VariantFiltration -R ", '${referenceArchive}',"/$genomeref ", q?--filterExpression "QD < 2.0" --filterExpression "ReadPosRankSum < -20.0" --filterExpression "FS > 200.0" --filterName QDFilter --filterName ReadPosFilter --filterName FSFilter -V ?, '${inFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_raw_$_[2].vcf -L ", '${referenceArchive}',"/$gatk_bait -o ", '${inFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_varrecal_$_[2]_filt.vcf ";
	    if (scalar(@sid) > 2) {
		for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) { #For all sample ids
		    if ( ($sid[$sampleid] =~ /(\d+)-(\d+|-\d+)-(\d+)(A|U)/) ) {#Match sampleID
			if ($2 eq 2) { #Parent
			    print GATK_VARREC "--pedigreeValidationType SILENT --pedigree ", '${outFamilyDir3}',"/$familyid.fam ";
			    last; #Only print once if a parent is found (required to include pedigree)
			}
		    }
		}
	    }
	}
	else { #WGS
	    print GATK_VARREC "\n\n#GATK VariantRecalibrator","\n\n";
	    print GATK_VARREC "java -Xmx12g -jar $gatk_path/GenomeAnalysisTK.jar -l INFO -T VariantRecalibrator -R ", '${referenceArchive}',"/$genomeref -recalFile ",'${outFamilyDir}', "/$_[0]", "_allchr_varrecal_$_[2].intervals -rscriptFile ", '${outFamilyDir}', "/$_[0]", "_allchr_varrecal_$_[2].intervals.plots.R -tranchesFile ",'${outFamilyDir}', "/$_[0]", "_allchr_varrecal_$_[2].intervals.tranches -resource:mills,VCF,known=true,training=true,truth=true,prior=12.0 ", '${referenceArchive}',"/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -an QD -an HaplotypeScore -an ReadPosRankSum -an FS --mode $_[2] -nt $maximum_cores -input ", '${inFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_raw_$_[2].vcf ";
	    if (scalar(@sid) > 2) {
		for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) { #For all sample ids
		    if ( ($sid[$sampleid] =~ /(\d+)-(\d+|-\d+)-(\d+)(A|U)/) ) {#Match sampleID
			if ($2 eq 2) { #Parent
			    print GATK_VARREC "--pedigreeValidationType SILENT --pedigree ", '${outFamilyDir3}',"/$familyid.fam ";
			    last; #Only print once if a parent is found (required to include pedigree)
			}
		    }
		}
	    }
	}
    }

    print GATK_VARREC "\n\n#GATK ApplyRecalibration","\n\n";
    print GATK_VARREC "#Samples", "\n";
    print GATK_VARREC 'inFamilyDir="',"$odf/$_[0]/$_[1]/GATK", '"', "\n";
    print GATK_VARREC 'inFamilyDir_2="',"$odf/$_[0]/$_[1]/GATK/intermediary", '"', "\n";
    print GATK_VARREC 'outFamilyDir="', "$odf/$_[0]/$_[1]/GATK", '"', "\n\n";

    if ($_[2] eq "BOTH") {
	
	if ($wgs == 0) { #Exome analysis 
	    
	    print GATK_VARREC "java -Xmx2g -jar $gatk_path/GenomeAnalysisTK.jar -l INFO -T ApplyRecalibration -R ", '${referenceArchive}',"/$genomeref -recalFile ",'${inFamilyDir_2}', "/$_[0]", "_allchr_varrecal_$_[2]_comb_ref.intervals -tranchesFile ",'${inFamilyDir_2}', "/$_[0]", "_allchr_varrecal_$_[2]_comb_ref.intervals.tranches --ts_filter_level 99.9 -o ", '${outFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_varrecal_$_[2]_comb_ref_filt.vcf -input ", '${inFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_raw_$_[2]_comb_ref.vcf -L ", '${referenceArchive}',"/$gatk_bait ";
	    if (scalar(@sid) > 2) {
		for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) { #For all sample ids
		    if ( ($sid[$sampleid] =~ /(\d+)-(\d+|-\d+)-(\d+)(A|U)/) ) {#Match sampleID
			if ($2 eq 2) { #Parent
			    print GATK_VARREC "--pedigreeValidationType SILENT --pedigree ", '${outFamilyDir3}',"/$familyid.fam ";
			    last; #Only print once if a parent is found (required to include pedigree)
			}
		    }
		}
	    }  
	    
###
#GATK SelectVariants
###
#Removes all genotype information for exome ref and recalulates meta-data info for remaining samples in new file.
	    
	    print GATK_VARREC "\n\n#GATK SelectVariants","\n\n";
	    print GATK_VARREC "java -Xmx2g -jar $gatk_path/GenomeAnalysisTK.jar -l INFO -T SelectVariants -R ", '${referenceArchive}',"/$genomeref -V: ", '${inFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_varrecal_$_[2]_comb_ref_filt.vcf -o ", '${outFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_varrecal_$_[2]_filt.vcf ";
	    
	    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) { #For all sample ids
		
		print GATK_VARREC "-sn $sid[$sampleid] ";
	    }
	}
	else { #WGS
	    print GATK_VARREC "java -Xmx2g -jar $gatk_path/GenomeAnalysisTK.jar -l INFO -T ApplyRecalibration -R ", '${referenceArchive}',"/$genomeref -recalFile ",'${inFamilyDir_2}', "/$_[0]", "_allchr_varrecal_$_[2].intervals -tranchesFile ",'${inFamilyDir_2}', "/$_[0]", "_allchr_varrecal_$_[2].intervals.tranches --ts_filter_level 99.9 -o ", '${outFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_varrecal_$_[2]_filt.vcf -input ", '${inFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_raw_$_[2].vcf ";
	    if (scalar(@sid) > 2) {
		for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) { #For all sample ids
		    if ( ($sid[$sampleid] =~ /(\d+)-(\d+|-\d+)-(\d+)(A|U)/) ) {#Match sampleID
			if ($2 eq 2) { #Parent
			    print GATK_VARREC "--pedigreeValidationType SILENT --pedigree ", '${outFamilyDir3}',"/$familyid.fam ";
			    last; #Only print once if a parent is found (required to include pedigree)
			}
		    }
		}		
	    }
	}
    }
    elsif ($_[2] eq "SNV") {
	
	if ($wgs == 0) { #Exome analysis 
	    
	    print GATK_VARREC "java -Xmx2g -jar $gatk_path/GenomeAnalysisTK.jar -l INFO -T ApplyRecalibration -R ", '${referenceArchive}',"/$genomeref -recalFile ",'${inFamilyDir_2}', "/$_[0]", "_allchr_varrecal_$_[2]_comb_ref.intervals -tranchesFile ",'${inFamilyDir_2}', "/$_[0]", "_allchr_varrecal_$_[2]_comb_ref.intervals.tranches --ts_filter_level 99.9 -o ", '${outFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_varrecal_$_[2]_comb_ref_filt.vcf -input ", '${inFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_raw_$_[2]_comb_ref.vcf -L ", '${referenceArchive}',"/$gatk_bait ";
	    if (scalar(@sid) > 2) {
		for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) { #For all sample ids
		    if ( ($sid[$sampleid] =~ /(\d+)-(\d+|-\d+)-(\d+)(A|U)/) ) {#Match sampleID
			if ($2 eq 2) { #Parent
			    print GATK_VARREC "--pedigreeValidationType SILENT --pedigree ", '${outFamilyDir3}',"/$familyid.fam ";
			    last; #Only print once if a parent is found (required to include pedigree)
			}
		    }
		}
	    }  
	    
###
#GATK SelectVariants
###
#Removes all genotype information for exome ref and recalulates meta-data info for remaining samples in new file.
	    
	    print GATK_VARREC "\n\n#GATK SelectVariants","\n\n";
	    print GATK_VARREC "java -Xmx2g -jar $gatk_path/GenomeAnalysisTK.jar -l INFO -T SelectVariants -R ", '${referenceArchive}',"/$genomeref -V: ", '${inFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_varrecal_$_[2]_comb_ref_filt.vcf -o ", '${outFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_varrecal_$_[2]_filt.vcf ";
	    
	    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) { #For all sample ids
		
		print GATK_VARREC "-sn $sid[$sampleid] ";
	    }
	}
	else { #WGS
	    print GATK_VARREC "java -Xmx2g -jar $gatk_path/GenomeAnalysisTK.jar -l INFO -T ApplyRecalibration -R ", '${referenceArchive}',"/$genomeref -recalFile ",'${inFamilyDir_2}', "/$_[0]", "_allchr_varrecal_$_[2].intervals -tranchesFile ",'${inFamilyDir_2}', "/$_[0]", "_allchr_varrecal_$_[2].intervals.tranches --ts_filter_level 99.9 -o ", '${outFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_varrecal_$_[2]_filt.vcf -input ", '${inFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_raw_$_[2].vcf ";
	    if (scalar(@sid) > 2) {
		for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) { #For all sample ids
		    if ( ($sid[$sampleid] =~ /(\d+)-(\d+|-\d+)-(\d+)(A|U)/) ) {#Match sampleID
			if ($2 eq 2) { #Parent
			    print GATK_VARREC "--pedigreeValidationType SILENT --pedigree ", '${outFamilyDir3}',"/$familyid.fam ";
			    last; #Only print once if a parent is found (required to include pedigree)
			}
		    }
		}		
	    }
	}
    }
    else { #INDELs

	if ($wgs == 0) { #Exome analysis   

        #Filtering is already done by hard filtering for exome indels, so no need for ApplyRecalibration	    

	}
	else { #WGS
	    print GATK_VARREC "java -Xmx2g -jar $gatk_path/GenomeAnalysisTK.jar -l INFO -T ApplyRecalibration -R ", '${referenceArchive}',"/$genomeref -recalFile ",'${inFamilyDir_2}', "/$_[0]", "_allchr_varrecal_$_[2].intervals -tranchesFile ",'${inFamilyDir_2}', "/$_[0]", "_allchr_varrecal_$_[2].intervals.tranches --ts_filter_level 99.9 -o ", '${outFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_varrecal_$_[2]_filt.vcf -input ", '${inFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_raw_$_[2].vcf ";
	    if (scalar(@sid) > 2) {
		for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) { #For all sample ids
		    if ( ($sid[$sampleid] =~ /(\d+)-(\d+|-\d+)-(\d+)(A|U)/) ) {#Match sampleID
			if ($2 eq 2) { #Parent
			    print GATK_VARREC "--pedigreeValidationType SILENT --pedigree ", '${outFamilyDir3}',"/$familyid.fam ";
			    last; #Only print once if a parent is found (required to include pedigree)
			}
		    }
		}
	    }
	}
    }

    print GATK_VARREC "\n\nwait", "\n\n";
    close(GATK_VARREC);   
    FIDSubmitJob(0,$familyid, 1, $_[2],$filename);
    return;
}

sub GATK_HapCall_ComVar { 
#GATK CombineVariants
#$_[0]= familyid
#$_[1]= $aligner, to choose the correct dir depending on what aligner has been used previously
#$_[2] = BOTH

    `mkdir -p $odf/$_[0]/$_[1]/info;`; #Creates the aligner folder and info data file directory
    `mkdir -p $odf/$_[0]/$_[1]/GATK`; #Creates the aligner folder, GATK data file directory
    `mkdir -p $ods/$_[0]/$_[1]`; #Creates the aligner folder script file directory
    $filename = "$ods/$_[0]/$_[1]/gatk_hapcal_comvar_$_[2]_$_[0].";   

    Checkfnexists($filename, $fnend);

#Info and Logg
    print STDOUT "Creating sbatch script GATK HaplotypeCaller CombineVariants and writing script file(s) to: ", $filename, "\n";print MASTERL "Creating sbatch script GATK HaplotypeCaller CombineVariants and writing script file(s) to: ", $filename, "\n";
    print STDOUT "Sbatch script GATK HaplotypeCaller CombineVariants data files will be written to: ", $odf,"/$_[0]/$_[1]/GATK", "\n";print MASTERL "Sbatch script GATK HaplotypeCaller CombineVariants data files will be written to: ", $odf,"/$_[0]/$_[1]/GATK", "\n";

    open (GATK_HAPCALCOMVAR, ">$filename") or die "Can't write to $filename: $!\n";
    
    print GATK_HAPCALCOMVAR "#! /bin/bash -l", "\n";
    print GATK_HAPCALCOMVAR "#SBATCH -A ", $aid, "\n";
    print GATK_HAPCALCOMVAR "#SBATCH -p node -n 1", "\n";
    print GATK_HAPCALCOMVAR "#SBATCH -C thin", "\n";	
    print GATK_HAPCALCOMVAR "#SBATCH -t 1:00:00", "\n";
    print GATK_HAPCALCOMVAR "#SBATCH -J GATK_HapCComVar_$_[2]_", $_[0], "\n";
    print GATK_HAPCALCOMVAR "#SBATCH -e $odf/$_[0]/$_[1]/info/gatk_hapcal_comvar_$_[2]_$_[0].", $fnt ,".stderr.txt", "\n";
    print GATK_HAPCALCOMVAR "#SBATCH -o $odf/$_[0]/$_[1]/info/gatk_hapcal_comvar_$_[2]_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print GATK_HAPCALCOMVAR "#SBATCH --mail-type=END", "\n";
	print GATK_HAPCALCOMVAR "#SBATCH --mail-type=FAIL", "\n";
	print GATK_HAPCALCOMVAR "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print GATK_HAPCALCOMVAR 'echo "Running on: $(hostname)"',"\n\n";
    print GATK_HAPCALCOMVAR "#Reference Archive", "\n";
    print GATK_HAPCALCOMVAR 'referenceArchive="', "$rd", '"', "\n\n"; 
    print GATK_HAPCALCOMVAR "#Samples", "\n";
    print GATK_HAPCALCOMVAR 'inFamilyDir="',"$odf/$_[0]/$_[1]/GATK/HapCall", '"', "\n";
    print GATK_HAPCALCOMVAR 'outFamilyDir="', "$odf/$_[0]/$_[1]/GATK", '"', "\n\n"; 
    
    print GATK_HAPCALCOMVAR "#GATK CombineVariants","\n\n";
    	   
    print GATK_HAPCALCOMVAR "java -Xmx2g -jar $gatk_path/GenomeAnalysisTK.jar -l INFO -T CombineVariants -R ", '${referenceArchive}',"/$genomeref ";

    for (my $chr=0;$chr<scalar(@chr);$chr++) { #For all chr	    
	print GATK_HAPCALCOMVAR "-V ", '${inFamilyDir}', "/$_[0]", "_real_recal_resrt_raw_",$chr[$chr], "_$_[2].vcf ";   
    }
    print GATK_HAPCALCOMVAR "-o ",'${outFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_raw_$_[2].vcf", "\n\n";

    print GATK_HAPCALCOMVAR "wait", "\n\n";
    close(GATK_HAPCALCOMVAR);   
    FIDSubmitJob(0,$familyid, 1, "MAIN",$filename);    
    return;
}

sub GATK_hapcal { 
#GATK HaplotypeCaller
#$_[0] = $familyid NOTE: not sampleid
#$_[1] = $aligner, to choose the correct dir depending on what aligner has been used previously
#$_[2] = BOTH
#$_[3] = Pos in @chr to start processing
#$_[4] = Pos in @chr to stop processing
#$_[5] = Java heap allocation (Gb).

    `mkdir -p $odf/$_[0]/$_[1]/info;`; #Creates the aligner folder and info data file directory
    `mkdir -p $odf/$_[0]/$_[1]/GATK/HapCall`; #Creates the aligner folder, GATK data file directory
    `mkdir -p $ods/$_[0]/$_[1]`; #Creates the aligner folder script file directory
    my $temp_start_chr = $_[3]+1;
    my $temp_stop_chr = $_[4];
    if ($_[4] == 26) {
	$temp_stop_chr = $_[4]-1;
    } 
    $filename = "$ods/$_[0]/$_[1]/gatk_hapcal_$_[2]_$_[0]_chr$temp_start_chr"."-$temp_stop_chr."; 
    
    Checkfnexists($filename, $fnend);

#Info and Logg
    print STDOUT "Creating sbatch script GATK HaplotypeCaller and writing script file(s) to: ", $filename, "\n";print MASTERL "Creating sbatch script GATK HaplotypeCaller and writing script file(s) to: ", $filename, "\n";
    print STDOUT "Sbatch script GATK HaplotypeCaller data files will be written to: ", $odf,"/$_[0]/$_[1]/GATK", "\n";print MASTERL "Sbatch script GATK HaplotypeCaller data files will be written to: ", $odf,"/$_[0]/$_[1]/GATK", "\n";

    open (GATK_HAPCAL, ">$filename") or die "Can't write to $filename: $!\n";
    
    print GATK_HAPCAL "#! /bin/bash -l", "\n";
    print GATK_HAPCAL "#SBATCH -A ", $aid, "\n";
    print GATK_HAPCAL "#SBATCH -p node -n $maximum_cores", "\n";
    print GATK_HAPCAL "#SBATCH -C thin", "\n";	
    print GATK_HAPCAL "#SBATCH -t 50:00:00", "\n";

    print GATK_HAPCAL "#SBATCH -J GATK_HAPCAL_$_[2]_", "$_[0]_chr",$temp_start_chr,"-", "$_[4]", "\n";
    print GATK_HAPCAL "#SBATCH -e $odf/$_[0]/$_[1]/info/gatk_hapcal_$_[2]_$_[0]_chr$temp_start_chr"."-$temp_stop_chr.", $fnt ,".stderr.txt", "\n";
    print GATK_HAPCAL "#SBATCH -o $odf/$_[0]/$_[1]/info/gatk_hapcal_$_[2]_$_[0]_chr$temp_start_chr"."-$temp_stop_chr.", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print GATK_HAPCAL "#SBATCH --mail-type=END", "\n";
	print GATK_HAPCAL "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print GATK_HAPCAL 'echo "Running on: $(hostname)"',"\n\n";
    print GATK_HAPCAL "#Reference Archive", "\n";
    print GATK_HAPCAL 'referenceArchive="', "$rd", '"', "\n\n"; 
    print GATK_HAPCAL "#Samples", "\n";
    print GATK_HAPCAL 'outFamilyDir="', "$odf/$_[0]/$_[1]/GATK/HapCall", '"', "\n\n"; 
    print GATK_HAPCAL 'outFamilyDir2="',"$odf/$_[0]", '"', "\n\n"; #For .fam file
    
    if ($_[3] == 0) { #Only for the first call of subroutine GATK_hapcal.
#Generate .fam file for later use in relevant GATK walkers (HaplotypeCaller, VariantscoreRequalibration etc)
	print GATK_HAPCAL "#Generating '.fam' file for GATK HAplotypeCaller","\n\n";
	print GATK_HAPCAL q?perl -nae 'my %sample_info;my $mother;my $father; while (<>) { my @F = split(/\t/,$_); if ($_!~/^#/) { if( $F[0]=~/(\d+)-(\d+|-\d+)-(\d+)(A|U)/) {} if ($F[3] == 1) {$father = $F[0];} if ($F[2] == 1) {$mother = $F[0];} if($3 % 2 == 1) {push (@{ $sample_info{$1}{$F[0]} }, "1");} else {push (@{ $sample_info{$1}{$F[0]} }, "2");} if ($4 eq "A") {push (@{ $sample_info{$1}{$F[0]} }, "2");} else {push (@{ $sample_info{$1}{$F[0]} }, "1");} } } for my $familyid (keys %sample_info ) { for my $sampleid (keys %{ $sample_info{$familyid} }) {print $familyid, " ", $sampleid, " ", $father, " ", $mother," "; for (my $i=0;$i<scalar(@{ $sample_info{$familyid}{$sampleid} } );$i++) {print $sample_info{$familyid}{$sampleid}[$i], " ";}print "\n"; } } last;' ?.$pedigree.q? > ${outFamilyDir2}/?.$familyid.q?.fam?, "\n\n";
    }

    print GATK_HAPCAL "#GATK HaplotypeCaller","\n\n";
    if ($_[4] == 26) { #Special case to enable processing of MT as well within same node for last call, overstrecthing a bit but should be fine
	for (my $chr=$_[3];$chr<$_[4]-1;$chr++) { #Determined by chr start and stop arguments given as input	    
	    print GATK_HAPCAL "java -Xmx$_[5]g -jar $gatk_path/GenomeAnalysisTK.jar -l INFO -T HaplotypeCaller -R ", '${referenceArchive}',"/$genomeref -D ", '${referenceArchive}',"/$gatk_recal_knset -stand_call_conf 30.0 -stand_emit_conf 30.0 --annotation BaseQualityRankSumTest --annotation ChromosomeCounts --annotation DepthOfCoverage --annotation FisherStrand --annotation HaplotypeScore --annotation InbreedingCoeff --annotation MappingQualityRankSumTest --annotation MappingQualityZero --annotation QualByDepth --annotation RMSMappingQuality --annotation ReadPosRankSumTest --annotation SpanningDeletions --annotation TandemRepeatAnnotator --annotation DepthPerAlleleBySample ";
	    if (scalar(@sid) > 2) {
		for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) { #For all sample ids
		    if ( ($sid[$sampleid] =~ /(\d+)-(\d+|-\d+)-(\d+)(A|U)/) ) {#Match sampleID
			if ($2 eq 2) { #Parent
			    print GATK_HAPCAL "--pedigree ", '${outFamilyDir2}',"/$familyid.fam ";		
			    last; #Only print once if a parent is found (required to include pedigree)
			}
		    }
		}		
	    }
	    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) { #For all sample ids
		
		my $tempinfile = $avcovfn{$sid[$sampleid]};
		print GATK_HAPCAL "-I $odf/$sid[$sampleid]/$_[1]/per_chr/GATK/$tempinfile","_", "$chr[$chr]","_real_recal_resrt.bam ";
	    } 
	    if ($wgs == 0) { #Exome analysis - Restrict analysis to padded target file
		
		print GATK_HAPCAL "-L ", '${referenceArchive}',"/$gatk_bait ";
	    }
	    print GATK_HAPCAL "-o ", '${outFamilyDir}', "/$_[0]", "_real_recal_resrt_raw_",$chr[$chr],"_$_[2].vcf &", "\n\n";
	}
    }
    else {
	for (my $chr=$_[3];$chr<$_[4];$chr++) { #Determined by chr start and stop arguments given as input
	    print GATK_HAPCAL "java -Xmx$_[5]g -jar $gatk_path/GenomeAnalysisTK.jar -l INFO -T HaplotypeCaller -R ", '${referenceArchive}',"/$genomeref -D ", '${referenceArchive}',"/$gatk_recal_knset -stand_call_conf 30.0 -stand_emit_conf 30.0 --annotation BaseQualityRankSumTest --annotation ChromosomeCounts --annotation DepthOfCoverage --annotation FisherStrand --annotation HaplotypeScore --annotation InbreedingCoeff --annotation MappingQualityRankSumTest --annotation MappingQualityZero --annotation QualByDepth --annotation RMSMappingQuality --annotation ReadPosRankSumTest --annotation SpanningDeletions --annotation TandemRepeatAnnotator --annotation DepthPerAlleleBySample ";
	    if (scalar(@sid) > 2) {
		for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) { #For all sample ids
		    if ( ($sid[$sampleid] =~ /(\d+)-(\d+|-\d+)-(\d+)(A|U)/) ) {#Match sampleID
			if ($2 eq 2) { #Parent
			    print GATK_HAPCAL "--pedigree ", '${outFamilyDir2}',"/$familyid.fam ";		
			    last; #Only print once if a parent is found (required to include pedigree)
			}
		    }
		}		
	    }
	    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) { #For all sample ids
		
		my $tempinfile = $avcovfn{$sid[$sampleid]};
		print GATK_HAPCAL "-I $odf/$sid[$sampleid]/$_[1]/per_chr/GATK/$tempinfile","_", "$chr[$chr]","_real_recal_resrt.bam ";
	    } 
	    if ($wgs == 0) { #Exome analysis - Restrict analysis to padded target file
		
		print GATK_HAPCAL "-L ", '${referenceArchive}',"/$gatk_bait ";
	    }
	    print GATK_HAPCAL "-o ", '${outFamilyDir}', "/$_[0]", "_real_recal_resrt_raw_",$chr[$chr],"_$_[2].vcf &", "\n\n";
	}	
    }
    print GATK_HAPCAL "\n\nwait", "\n\n";    
    
    close(GATK_HAPCAL);  
    FIDSubmitJob(0,$familyid, 3, "MAIN",$filename); #Arg2 eq 3 for parallel execution  
    return;
}

sub GATK_unigt { 
#GATK UnifiedGenotyper
#$_[0] = $familyid NOTE: not sampleid
#$_[1] = $aligner, to choose the correct dir depending on what aligner has been used previously
#$_[2] = $gatk_unigt_snps (SNV) or $gatk_unigt_indels (INDEL)

    `mkdir -p $odf/$_[0]/$_[1]/info;`; #Creates the aligner folder and info data file directory
    `mkdir -p $odf/$_[0]/$_[1]/GATK`; #Creates the aligner folder, GATK data file directory
    `mkdir -p $ods/$_[0]/$_[1]`; #Creates the aligner folder script file directory
    $filename = "$ods/$_[0]/$_[1]/gatk_unigt_$_[2]_$_[0]."; 

    Checkfnexists($filename, $fnend);

#Info and Logg
    print STDOUT "Creating sbatch script GATK UnifiedGenotyper and writing script file(s) to: ", $filename, "\n";print MASTERL "Creating sbatch script GATK UnifiedGenotyper and writing script file(s) to: ", $filename, "\n";
    print STDOUT "Sbatch script GATK UnifiedGenotyper data files will be written to: ", $odf,"/$_[0]/$_[1]/GATK", "\n";print MASTERL "Sbatch script GATK UnifiedGenotyper data files will be written to: ", $odf,"/$_[0]/$_[1]/GATK", "\n";

    open (GATK_UNIGT, ">$filename") or die "Can't write to $filename: $!\n";
    
    print GATK_UNIGT "#! /bin/bash -l", "\n";
    print GATK_UNIGT "#SBATCH -A ", $aid, "\n";
    print GATK_UNIGT "#SBATCH -p node -n $maximum_cores", "\n";
    print GATK_UNIGT "#SBATCH -C thin", "\n";	
    print GATK_UNIGT "#SBATCH -t 40:00:00", "\n";

    print GATK_UNIGT "#SBATCH -J GATK_UNIGT_$_[2]_", $_[0], "\n";
    print GATK_UNIGT "#SBATCH -e $odf/$_[0]/$_[1]/info/gatk_unigt_$_[2]_$_[0].", $fnt ,".stderr.txt", "\n";
    print GATK_UNIGT "#SBATCH -o $odf/$_[0]/$_[1]/info/gatk_unigt_$_[2]_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print GATK_UNIGT "#SBATCH --mail-type=END", "\n";
	print GATK_UNIGT "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print GATK_UNIGT 'echo "Running on: $(hostname)"',"\n\n";
    print GATK_UNIGT "#Reference Archive", "\n";
    print GATK_UNIGT 'referenceArchive="', "$rd", '"', "\n\n"; 
    print GATK_UNIGT "#Samples", "\n";
    print GATK_UNIGT 'outFamilyDir="', "$odf/$_[0]/$_[1]/GATK", '"', "\n\n"; 
    print GATK_UNIGT "#GATK UnifiedGenotyper","\n\n";

    print GATK_UNIGT "java -Xmx12g -jar $gatk_path/GenomeAnalysisTK.jar -l INFO -T UnifiedGenotyper -R ", '${referenceArchive}',"/$genomeref -D ", '${referenceArchive}',"/$gatk_recal_knset -glm ";
    
    if ($_[2] eq "SNV") { #UnifiedGT only takes "SNP" and not "SNV" as model input
	print GATK_UNIGT "SNP ";
    }
    else { #INDEL
	print GATK_UNIGT "$_[2] ";
    }
    print GATK_UNIGT "-nt $maximum_cores -stand_call_conf 30.0 -stand_emit_conf 30.0 --min_base_quality_score 20 --annotation QualByDepth --annotation HaplotypeScore --annotation MappingQualityRankSumTest --annotation ReadPosRankSumTest --annotation FisherStrand --annotation RMSMappingQuality --annotation DepthOfCoverage ";
    if ($wgs == 0) { #Exome analysis - Restrict analysis to padded target file
	
	print GATK_UNIGT "-L ", '${referenceArchive}',"/$gatk_bait ";
    }
    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) { #For all sample ids
	
	my $tempinfile = $avcovfn{$sid[$sampleid]};
	    
	    print GATK_UNIGT "-I $odf/$sid[$sampleid]/$_[1]/GATK/$tempinfile","_allchr_real_recal_resrt.bam ";
    }    	    
    
    print GATK_UNIGT "-o ", '${outFamilyDir}', "/$_[0]", "_allchr_real_recal_resrt_raw_$_[2].vcf";

    print GATK_UNIGT "\n\nwait", "\n\n";    
    close(GATK_UNIGT);
    FIDSubmitJob(0,$familyid, 1, "MAIN",$filename);
    return;
}

sub GATK_recal { 
#GATK recalibration
#$_[0]= $sampleid
#$_[1]= $aligner, to choose the correct dir depending on what aligner has been used previously

    `mkdir -p $odf/$_[0]/$_[1]/info;`; #Creates the aligner folder and info data file directory
    `mkdir -p $odf/$_[0]/$_[1]/per_chr/GATK/intermediary`; #Creates the aligner folder, per chromosome and GATK intermediary data file directory
    `mkdir -p $odf/$_[0]/$_[1]/GATK/`; #Creates the aligner folder, GATK all chr data file directory
    `mkdir -p $ods/$_[0]/$_[1];`; #Creates the aligner folder script file directory

    $filename = "$ods/$_[0]/$_[1]/gatk_recal_$avcovfn{$_[0]}.";   

    Checkfnexists($filename, $fnend);

#Info and Logg
    print STDOUT "Creating sbatch script GATK BaseRecalibrator/PrintReads and writing script file(s) to: ", $filename, "\n";print MASTERL "Creating sbatch script GATK BaseRecalibrator/PrintReads and writing script file(s) to: ", $filename, "\n";
    print STDOUT "Sbatch script GATK BaseRecalibrator/PrintReads data files will be written to: ", $odf,"/$_[0]/$_[1]/per_chr/GATK", "\n";print MASTERL "Sbatch script GATK BaseRecalibrator/PrintReads data files will be written to: ", $odf,"/$_[0]/$_[1]/per_chr/GATK", "\n";

    open (GATK_RECAL, ">$filename") or die "Can't write to $filename: $!\n";
    
    print GATK_RECAL "#! /bin/bash -l", "\n";
    print GATK_RECAL "#SBATCH -A ", $aid, "\n";
    print GATK_RECAL "#SBATCH -p node -n $maximum_cores", "\n";
    print GATK_RECAL "#SBATCH -C thin", "\n";	
    print GATK_RECAL "#SBATCH -t 60:00:00", "\n";
    print GATK_RECAL "#SBATCH -J GATK_RECAL_", $_[0], "\n";
    print GATK_RECAL "#SBATCH -e $odf/$_[0]/$_[1]/info/gatk_recal_$_[0].", $fnt ,".stderr.txt", "\n";
    print GATK_RECAL "#SBATCH -o $odf/$_[0]/$_[1]/info/gatk_recal_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print GATK_RECAL "#SBATCH --mail-type=END", "\n";
	print GATK_RECAL "#SBATCH --mail-type=FAIL", "\n";
	print GATK_RECAL "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print GATK_RECAL 'echo "Running on: $(hostname)"',"\n\n";
    print GATK_RECAL "#Reference Archive", "\n";
    print GATK_RECAL 'referenceArchive="', "$rd", '"', "\n\n"; 
    print GATK_RECAL "#Samples", "\n";
    print GATK_RECAL 'inSampleDir="',"$odf/$_[0]/$_[1]/per_chr/GATK", '"', "\n";
    print GATK_RECAL 'outSampleDir="', "$odf/$_[0]/$_[1]/per_chr/GATK/intermediary", '"', "\n\n";   

    my $tempinfile = $avcovfn{$_[0]};
    my $core_Counter=1;

    print GATK_RECAL "#GATK CountCovariates","\n\n";
    for (my $chr=0;$chr<scalar(@chr);$chr++) { #For all chr	    
	
	if ($chr eq $core_Counter*$maximum_cores) { #Using only $maximum_cores 
	    
	    print GATK_RECAL "wait", "\n\n";
	    $core_Counter=$core_Counter+1;
	}
	print GATK_RECAL "java -Xmx3g -Djava.io.tmpdir=/proj/$aid/private/nobackup",'/$SLURM_JOB_ID/', "$chr[$chr]/ -jar $gatk_path/GenomeAnalysisTK.jar -l INFO -T BaseRecalibrator -cov ReadGroupCovariate -cov ContextCovariate -cov CycleCovariate -cov QualityScoreCovariate -cov ReadGroupCovariate -R ", '${referenceArchive}',"/$genomeref -knownSites ", '${referenceArchive}', "/$gatk_recal_knset -I ",  '${inSampleDir}', "/$tempinfile","_", "$chr[$chr]","_real.bam -o ", '${outSampleDir}', "/$tempinfile","_", "$chr[$chr]", "_real_recal.grp &", "\n\n";
    }
    $core_Counter=1; #Resetting
    print GATK_RECAL "wait", "\n\n";
    print GATK_RECAL "#GATK ","\n\n";
    print GATK_RECAL "#Samples", "\n";
    print GATK_RECAL 'inSampleDir="',"$odf/$_[0]/$_[1]/per_chr/GATK", '"', "\n";
    print GATK_RECAL 'inSampleDir_2="', "$odf/$_[0]/$_[1]/per_chr/GATK/intermediary",'"', "\n";
    print GATK_RECAL 'outSampleDir="', "$odf/$_[0]/$_[1]/per_chr/GATK", '"', "\n\n";

    for (my $chr=0;$chr<scalar(@chr);$chr++) { #For all chr	    
	
	if ($chr eq $core_Counter*$maximum_cores) { #Using only $maximum_cores cores
	    
	    print GATK_RECAL "wait", "\n\n";
	    $core_Counter=$core_Counter+1;
	}
	print GATK_RECAL "java -Xmx3g -jar $gatk_path/GenomeAnalysisTK.jar -l INFO -T PrintReads -R ", '${referenceArchive}',"/$genomeref -I ", '${inSampleDir}', "/$tempinfile","_", "$chr[$chr]","_real.bam -o ", '${outSampleDir}', "/$tempinfile","_", "$chr[$chr]", "_real_recal.bam -BQSR ", '${inSampleDir_2}', "/$tempinfile","_", "$chr[$chr]", "_real_recal.grp &", "\n\n";
	
    }
    $core_Counter=1; #Resetting for new infile
    print GATK_RECAL "wait", "\n\n";
###
#SamTools sort on real and recal files
###
for (my $chr=0;$chr<scalar(@chr);$chr++) { #For all chr	    
	
	if ($chr eq $core_Counter*$maximum_cores) { #Using only $maximum_cores cores
	    
	    print GATK_RECAL "wait", "\n\n";
	    $core_Counter=$core_Counter+1;
	}
	print GATK_RECAL "samtools sort ", '${inSampleDir}', "/$tempinfile","_", "$chr[$chr]", "_real_recal.bam ", '${outSampleDir}', "/$tempinfile","_", "$chr[$chr]", "_real_recal_resrt &","\n\n"; #samtools sort adds .bam
	
    }
    $core_Counter=1; #Resetting for new infile
    print GATK_RECAL "wait", "\n\n";
#Create @RG group    
    print GATK_RECAL "samtools view -H ", '${inSampleDir}', "/$tempinfile","_", "$chr[0]", "_real_recal_resrt.bam  > ", '${inSampleDir}', "/$tempinfile","_", "$chr[0]","_real_recal_resrt_header.sam","\n\n"; #NOTE to add @RG group header in later merge
    print GATK_RECAL 'outSampleDir="', "$odf/$_[0]/$_[1]/GATK", '"', "\n\n";
    for (my $chr=0;$chr<scalar(@chr);$chr++) { #For all chr	    
	
	if ($chr eq 0) {
	    
	    print GATK_RECAL "samtools merge -f -h ", '${inSampleDir}', "/$tempinfile","_", "$chr[$chr]","_real_recal_resrt_header.sam ";
	    print GATK_RECAL '${outSampleDir}', "/$tempinfile","_", "allchr_real_recal_resrt.bam ";
	}
	
	print GATK_RECAL '${inSampleDir}', "/$tempinfile","_", "$chr[$chr]", "_real_recal_resrt.bam ";   
    }
    print GATK_RECAL "\n\nsamtools index ", '${outSampleDir}', "/$tempinfile","_allchr_real_recal_resrt.bam", "\n\n";
    print GATK_RECAL "wait", "\n\n";

#Index recal and resrt individual files for downstream use with HaplotypeCaller
    $core_Counter=1; #Resetting for new infile

    for (my $chr=0;$chr<scalar(@chr);$chr++) { #For all chr	    
	if ($chr eq $core_Counter*$maximum_cores) { #Using only $maximum_cores cores
	    
	    print GATK_RECAL "wait", "\n\n";
	    $core_Counter=$core_Counter+1;
	}
	print GATK_RECAL "samtools index ", '${inSampleDir}', "/$tempinfile","_", "$chr[$chr]","_real_recal_resrt.bam &", "\n\n";   
    }
    $core_Counter=1; #Resetting
    print GATK_RECAL "wait", "\n\n";
    
    close(GATK_RECAL);   
    FIDSubmitJob($_[0],$familyid, 1, "MAIN",$filename);
    return;
}

sub GATK_real { 
#GATK realignment
#$_[0]= $sampleid
#$_[1]= $aligner, to choose the correct dir depending on what aligner has been used previously

    `mkdir -p $odf/$_[0]/$_[1]/info;`; #Creates the aligner folder and info data file directory
    `mkdir -p $odf/$_[0]/$_[1]/per_chr/GATK/intermediary`; #Creates the aligner folder, per chromosome and GATK intermediary data file directory
    `mkdir -p $ods/$_[0]/$_[1];`; #Creates the aligner folder script file directory
    $filename = "$ods/$_[0]/$_[1]/gatk_real_$avcovfn{$_[0]}.";   

    Checkfnexists($filename, $fnend);

#Info and Logg
    print STDOUT "Creating sbatch script GATK RealignerTargetCreator/IndelRealigner and writing script file(s) to: ", $filename, "\n";print MASTERL "Creating sbatch script GATK RealignerTargetCreator/IndelRealigner and writing script file(s) to: ", $filename, "\n";
    print STDOUT "Sbatch script GATK RealignerTargetCreator/IndelRealigner data files will be written to: ", $odf,"/$_[0]/$_[1]/per_chr/GATK", "\n";print MASTERL "Sbatch script GATK RealignerTargetCreator/IndelRealigner data files will be written to: ", $odf,"/$_[0]/$_[1]/per_chr/GATK", "\n";

    open (GATK_REAL, ">$filename") or die "Can't write to $filename: $!\n";
    
    print GATK_REAL "#! /bin/bash -l", "\n";
    print GATK_REAL "#SBATCH -A ", $aid, "\n";
    print GATK_REAL "#SBATCH -p node -n $maximum_cores", "\n";
    print GATK_REAL "#SBATCH -C thin", "\n";	
    print GATK_REAL "#SBATCH -t 40:00:00", "\n";
    print GATK_REAL "#SBATCH -J GATK_REAL_", $_[0], "\n";
    print GATK_REAL "#SBATCH -e $odf/$_[0]/$_[1]/info/gatk_real_$_[0].", $fnt ,".stderr.txt", "\n";
    print GATK_REAL "#SBATCH -o $odf/$_[0]/$_[1]/info/gatk_real_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print GATK_REAL "#SBATCH --mail-type=END", "\n";
	print GATK_REAL "#SBATCH --mail-type=FAIL", "\n";
	print GATK_REAL "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print GATK_REAL 'echo "Running on: $(hostname)"',"\n\n";
    print GATK_REAL "#Reference Archive", "\n";
    print GATK_REAL 'referenceArchive="', "$rd", '"', "\n\n"; 
    print GATK_REAL "#Samples", "\n";
    print GATK_REAL 'inSampleDir="',"$odf/$_[0]/$_[1]/per_chr", '"', "\n";
    print GATK_REAL 'outSampleDir="', "$odf/$_[0]/$_[1]/per_chr/GATK/intermediary", '"', "\n\n";   

    my $tempinfile = $avcovfn{$_[0]};
    my $core_Counter=1;

    print GATK_REAL "#GATK RealignerTargetCreator","\n\n";
    for (my $chr=0;$chr<scalar(@chr);$chr++) { #For all chr	    
	
	if ($chr eq $core_Counter*$maximum_cores) { #Using only $maximum_cores cores
	    
	    print GATK_REAL "wait", "\n\n";
	    $core_Counter=$core_Counter+1;
	}
	print GATK_REAL "java -Xmx3g -Djava.io.tmpdir=/proj/$aid/private/nobackup",'/$SLURM_JOB_ID/', "$chr[$chr]/ -jar $gatk_path/GenomeAnalysisTK.jar -l INFO -T RealignerTargetCreator -R ", '${referenceArchive}',"/$genomeref -known ", '${referenceArchive}', "/$gatk_real_knset1 -known ", '${referenceArchive}',"/$gatk_real_knset2 -I ",  '${inSampleDir}', "/$tempinfile","_", "$chr[$chr].bam -o ", '${outSampleDir}', "/$tempinfile","_", "$chr[$chr]", "_real.intervals &", "\n\n";
	
    }
    $core_Counter=1; #Resetting
    print GATK_REAL "wait", "\n\n";
    print GATK_REAL "#GATK IndelRealigner","\n\n";
    print GATK_REAL "#Samples", "\n";
    print GATK_REAL 'inSampleDir="',"$odf/$_[0]/$_[1]/per_chr", '"', "\n";
    print GATK_REAL 'inSampleDir_2="', "$odf/$_[0]/$_[1]/per_chr/GATK/intermediary", '"', "\n";
    print GATK_REAL 'outSampleDir="', "$odf/$_[0]/$_[1]/per_chr/GATK", '"', "\n\n";

    for (my $chr=0;$chr<scalar(@chr);$chr++) { #For all chr	    
	
	if ($chr eq $core_Counter*$maximum_cores) { #Using only $maximum_cores cores
	    
	    print GATK_REAL "wait", "\n\n";
	    $core_Counter=$core_Counter+1;
	}
	print GATK_REAL "java -Xmx3g -jar $gatk_path/GenomeAnalysisTK.jar -l INFO -T IndelRealigner -R ", '${referenceArchive}',"/$genomeref -known ", '${referenceArchive}', "/$gatk_real_knset1 -known ", '${referenceArchive}',"/$gatk_real_knset2 -I ", '${inSampleDir}', "/$tempinfile","_", "$chr[$chr].bam -o ", '${outSampleDir}', "/$tempinfile","_", "$chr[$chr]", "_real.bam -targetIntervals ", '${inSampleDir_2}', "/$tempinfile","_", "$chr[$chr]", "_real.intervals &", "\n\n";
	
    }
    $core_Counter=1; #Resetting for new infile
    print GATK_REAL "wait", "\n\n";
      
    close(GATK_REAL);
    FIDSubmitJob($_[0],$familyid, 1, "MAIN",$filename); 
    return;
}

sub SamtoolsViewSChr { 
#samtools view split genome.bam file to chr.bam files and index them
#$_[0]= $sampleid
#$_[1]= $aligner, to choose the correct dir depending on what aligner has been used previously

    `mkdir -p $odf/$_[0]/$_[1]/info;`; #Creates the aligner folder and info data file directory
    `mkdir -p $odf/$_[0]/$_[1]/per_chr;`; #Creates the aligner folder and per chromosome data file directory
    `mkdir -p $ods/$_[0]/$_[1];`; #Creates the aligner folder script file directory

    $filename = "$ods/$_[0]/$_[1]/samtools_view_split_to_chr_$avcovfn{$_[0]}.";
    Checkfnexists($filename, $fnend);
#Info and Logg
    print STDOUT "Creating sbatch script samtools and writing script file(s) to: ", $filename, "\n";print MASTERL "Creating sbatch script samtools and writing script file(s) to: ", $filename, "\n";
    print STDOUT "Sbatch script samtools view data files will be written to: ", $odf,"/$_[0]/$_[1]/per_chr", "\n";print MASTERL "Sbatch script samtools view data files will be written to: ", $odf,"/$_[0]/$_[1]/per_chr", "\n";

    open (STVSCHR, ">$filename") or die "Can't write to $filename: $!\n";
    
    print STVSCHR "#! /bin/bash -l", "\n";
    print STVSCHR "#SBATCH -A ", $aid, "\n";
    print STVSCHR "#SBATCH -p node -n $maximum_cores", "\n";
    print STVSCHR "#SBATCH -C thin", "\n";	
    print STVSCHR "#SBATCH -t 5:00:00", "\n"; 
    print STVSCHR "#SBATCH -J STVSCHR_", $_[0], "\n";
    print STVSCHR "#SBATCH -e $odf/$_[0]/$_[1]/info/samtools_view_split_to_chr_$avcovfn{$_[0]}.", $fnt ,".stderr.txt", "\n";
    print STVSCHR "#SBATCH -o $odf/$_[0]/$_[1]/info/samtools_view_split_to_chr_$avcovfn{$_[0]}.", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print STVSCHR "#SBATCH --mail-type=END", "\n";
	print STVSCHR "#SBATCH --mail-type=FAIL", "\n";
	print STVSCHR "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print STVSCHR 'echo "Running on: $(hostname)"',"\n\n";
    print STVSCHR "#Samples", "\n";
    
    my $core_Counter=1;
    my $tempinfile = $avcovfn{$_[0]};
    
    for (my $chr=0;$chr<scalar(@chr);$chr++) { #For all chr
	
	if ($chr eq $core_Counter*$maximum_cores) { #Using only $maximum_cores cores
	    
	    print STVSCHR "wait", "\n\n";
	    $core_Counter=$core_Counter+1;
	}
	print STVSCHR "samtools view -b -o $odf/$_[0]/$_[1]/per_chr/$tempinfile", "_", "$chr[$chr].bam $odf/$_[0]/$_[1]/$tempinfile", ".bam $chr[$chr] &", "\n\n";
    }
    
    print STVSCHR "wait", "\n\n";
    $core_Counter=1; #Reset
    for (my $chr=0;$chr<scalar(@chr);$chr++) { #For all chr
	
	if ($chr eq $core_Counter*$maximum_cores) { #Using only $maximum_cores cores
	    
	    print STVSCHR "wait", "\n\n";
	    $core_Counter=$core_Counter+1;
	}
	print STVSCHR "samtools index $odf/$_[0]/$_[1]/per_chr/$tempinfile", "_", "$chr[$chr].bam &", "\n\n";
    }
    print STVSCHR "wait", "\n\n";
    close(STVSCHR);
    FIDSubmitJob($_[0],$familyid, 0, "MAIN",$filename);
    return;
}

sub FIDSubmitJob {
#Submits all jobIDs to SLURM using SLURM dependencies. The first path is MAIN and any subsequent splits into other paths later is handled by adding relevant previous jobIDs to the new paths key in jobID{path_key} hash. The subroutine supports parallel job within each step and submission which do not leave any dependencies. Currently any path downstream of MAIN inherits the relevant previous jobIds, but it is not possible to merge to splitet paths downstream of main to each other.
#$_[0] = sampleid or 0 when family is supplied
#$_[1] = familyID
#$_[2] = Dependencies
#$_[3] = Path (MAIN, SNV, INDEL, BOTH). MAIN is before there is any split.
#$_[4] = sbatch filename to submit.

###
#Dependencies
###
#0 = Not dependent on earlier scripts
#1 = Dependent on earlier scripts (within sampleID_path or familyID_path)
#2 = Dependent on earlier scripts (within sampleID_path or familyID_path), but are self cul-de-sâcs. 
#3 = Dependent on earlier scripts and executed in parallel within step

    my $jobids=""; #Create string with all previous jobIDs
    my $previousjobcounter = 0; #Count previous submissions
    my $ret; #Return jobID
    my $sid_chainkey = "$_[0]_$_[3]"; #(sampleID_Path
    my $fid_chainkey = "$_[1]_$_[3]"; #(familyID_Path
    my $sid_parallel_chainkey = "$_[0]_parallel_$_[3]"; #(sampleID_parallel_Path
    my $fid_parallel_chainkey = "$_[1]_parallel_$_[3]"; #(familyID_parallel_Path
    my $jobID; #The jobID that is returned from submission
    
    if ($_[2] == 0) { #Initiate chain - No dependencies
	$ret = `sbatch $_[4]`; #No jobs have been run: submit
	($jobID) = ($ret =~ /Submitted batch job (\d+)/);
	push ( @{ $jobID{$sid_chainkey} }, $jobID); #Add jobID to hash{$sampleID}[]
	push ( @{ $jobID{$fid_chainkey} }, $jobID); #Add jobID to hash{$familyID}[]. Required to enable later test for all subjects sampleID_MAIN have finished before merging to family.
    }
    else { #Dependent on earlier scripts and/or parallel. JbIDs that do not leave dependencies do not get pushed to jobID hash
	
	if ($_[0]) { #BEFORE merging to familyID

	    if ( ($_[2] != 3) && ($jobID{$sid_parallel_chainkey}) ){ #If not a parallel job and a parallel job within CURRENT PATH has prev been processed. Check if previous step was parallel and adds previous parallel jobs that have previously been submitted.
		for (my $job=0;$job<scalar( @{ $jobID{$sid_parallel_chainkey} });$job++) {
		    my $seen_jobIDs_counter = 0;
		    if ($jobID{$sid_chainkey}) {#If any previous jobIds within current chain exists else go ahead and add
			for (my $current_job=0;$current_job<scalar( @{ $jobID{$sid_chainkey} });$current_job++) {
			    if ($jobID{$sid_chainkey}[$current_job] =~/$jobID{$sid_parallel_chainkey}[$job]/ ) { #Only add if not already present
				$seen_jobIDs_counter++;
			    }
			}
			if ($seen_jobIDs_counter eq 0) {#JobID was not present in CURRENT path
			    push ( @{ $jobID{$sid_chainkey} }, $jobID{$sid_parallel_chainkey}[$job]); #Add jobID to hash{$}			
			}
		    }
		    else { #Go ahead and add
			push ( @{ $jobID{$sid_chainkey} }, $jobID{$sid_parallel_chainkey}[$job]); #Add jobID to hash{$}
		    }   	
		}
	    }	   
	    if ( ($_[3] eq "MAIN") && ($jobID{$sid_chainkey})  ) { #Check for any previous jobIDs within path MAIN. Test for pevious must be done to allow initiating from broken chain
		for (my $job=0;$job<scalar( @{ $jobID{$sid_chainkey} });$job++) {	
		    if ( ($job == 0) && (scalar( @{ $jobID{$sid_chainkey} } )== 1) ) {#Only 1 previous jobID 
			$jobids .= ":$jobID{$sid_chainkey}[$job]"; #first and last jobID start with ":" and end without ":"
		    }
		    elsif ($job == 0 ) {
			$jobids .= ":$jobID{$sid_chainkey}[$job]:"; #first jobID start with :
		    }
		    elsif ($job eq ( scalar( @{ $jobID{$sid_chainkey} }) -1) ) {
			$jobids .= "$jobID{$sid_chainkey}[$job]"; #last jobID finish without :
		    }
		    else {
			$jobids .= "$jobID{$sid_chainkey}[$job]:";
		    }
		}
	    }
	    if ( $_[3] ne "MAIN" ) { #Check for any previous jobIDs within path current PATH
		my $sid_main_parallel_chainkey = "$_[0]_parallel_MAIN"; 
		if ( ($_[2] != 3) && ($jobID{$sid_main_parallel_chainkey}) ){ #If not a parallel job and a parallel job within MAIN path has prev been processed. Check if previous step was parallel and adds previous parallel jobs that have previously been submitted.
		    for (my $job=0;$job<scalar( @{ $jobID{$sid_main_parallel_chainkey} });$job++) { #All prev parallel MAIN jobIDs
			my $seen_jobIDs_counter = 0;
			if ($jobID{$sid_chainkey}) {#If any previous jobIds within current chain exists - else go ahead and add
			    for (my $current_job=0;$current_job<scalar( @{ $jobID{$sid_chainkey} });$current_job++) { #CURRENT path
				if ($jobID{$sid_chainkey}[$current_job] =~/$jobID{$sid_main_parallel_chainkey}[$job]/ ) { #Only add if not already present
				    $seen_jobIDs_counter++;
				}
			    }
			    if ($seen_jobIDs_counter eq 0) {#JobID was not present in CURRENT path
				push ( @{ $jobID{$sid_chainkey} }, $jobID{$sid_main_parallel_chainkey}[$job]); #Add jobID to hash{$}
			    }
			}
			else { #Go ahead and add
			    push ( @{ $jobID{$sid_chainkey} }, $jobID{$sid_main_parallel_chainkey}[$job]); #Add jobID to hash{$}
			}
		    }
		}
		my $sid_main_chainkey = "$_[0]_MAIN";
		if ($jobID{$sid_main_chainkey}) { #Any MAIN jobIDs necessary for broken chains, since this will be empty then
		    for (my $job=0;$job<scalar( @{ $jobID{$sid_main_chainkey} });$job++) { #Prev MAIn jobIDs
			my $seen_jobIDs_counter = 0; 
			if ($jobID{$sid_chainkey}) {#If any previous jobIds within current chain exists else go ahead and add
			    for (my $current_job=0;$current_job<scalar( @{ $jobID{$sid_chainkey} });$current_job++) { #CURRENT path
				if ($jobID{$sid_chainkey}[$current_job] =~/$jobID{$sid_main_chainkey}[$job]/) {
				    $seen_jobIDs_counter++;
				}
			    }
			    if ($seen_jobIDs_counter eq 0) {#JobID was not present in CURRENT path
				push ( @{ $jobID{$sid_chainkey} }, $jobID{$sid_main_chainkey}[$job]); #Add jobID to hash{$}
			    }
			}
			else  { #Go ahead and add
			    push ( @{ $jobID{$sid_chainkey} }, $jobID{$sid_main_chainkey}[$job]); #Add jobID to hash{$}
			}
		    }
		}
		if ($jobID{$sid_chainkey}) {
		    for (my $job=0;$job<scalar( @{ $jobID{$sid_chainkey} });$job++) {	
			if ( ($job == 0) && (scalar( @{ $jobID{$sid_chainkey} } )== 1) ) {#Only 1 previous jobID 
			    $jobids .= ":$jobID{$sid_chainkey}[$job]"; #first and last jobID start with ":" and end without ":"
			}
			elsif ($job == 0 ) {
			    $jobids .= ":$jobID{$sid_chainkey}[$job]:"; #first jobID start with :
			}
			elsif ($job eq ( scalar( @{ $jobID{$sid_chainkey} }) -1) ) {
			    $jobids .= "$jobID{$sid_chainkey}[$job]"; #last jobID finish without :
			}
			else {
			    $jobids .= "$jobID{$sid_chainkey}[$job]:";
			}
		    }
		}
	    }
	    if ($jobids) {
		$ret = `sbatch --dependency=afterok$jobids $_[4]`; #Supply with dependency of previous jobs that this one is dependent on
		($jobID) = ($ret =~ /Submitted batch job (\d+)/);
	    }
	    else {
		$ret = `sbatch $_[4]`; #No jobs have been run: submit
		($jobID) = ($ret =~ /Submitted batch job (\d+)/);
	    }
	    if ($_[2] == 1) { #Ordinary job push to array
		push ( @{ $jobID{$sid_chainkey} }, $jobID); #Add jobID to hash{$sampleID}[]
		push ( @{ $jobID{$fid_chainkey} }, $jobID); #Add jobID to hash{$familyID}[]. Required to enable later test for all subjects sampleID_MAIN have finished before merging to family.
	    }
	    if ($_[2] == 3) { #Parallel job wait to push to array until all parallel jobs are finished within step
		push ( @{ $jobID{$sid_parallel_chainkey} }, $jobID); #Add jobID to hash{$sampleID_parallel}[].
		push ( @{ $jobID{$fid_chainkey} }, $jobID); #Add jobID to hash{$familyID}[]. Required to enable later test for all subjects sampleID_MAIN have finished before merging to family.
	    }
	}
	else { #AFTER merging to familyID
	    if ( ($_[2] != 3) && ($jobID{$fid_parallel_chainkey}) ){ #If not a parallel job and a parallel job within CURRENT PATH has prev been processed. Check if previous step was parallel and adds previous parallel jobs that have previously been submitted.
		for (my $job=0;$job<scalar( @{ $jobID{$fid_parallel_chainkey} });$job++) {
		    my $seen_jobIDs_counter = 0;
		    if ($jobID{$fid_chainkey}) {#If any previous jobIds within current chain exists else go ahead and add
			for (my $current_job=0;$current_job<scalar( @{ $jobID{$fid_chainkey} });$current_job++) {
			    if ($jobID{$fid_chainkey}[$current_job] =~/$jobID{$fid_parallel_chainkey}[$job]/ ) { #Only add if not already present
				$seen_jobIDs_counter++;
			    }
			}
			if ($seen_jobIDs_counter eq 0) {#JobID was not present in CURRENT path
			    push ( @{ $jobID{$fid_chainkey} }, $jobID{$fid_parallel_chainkey}[$job]); #Add jobID to hash{$}			
			}
		    }
		    else { #Go ahead and add
			push ( @{ $jobID{$fid_chainkey} }, $jobID{$fid_parallel_chainkey}[$job]); #Add jobID to hash{$}
		    }   	
		}
	    }
	    if ( ($_[3] eq "MAIN") && ($jobID{$fid_chainkey})  ) { #Check for any previous jobIDs within path MAIN. Test for pevious must be done to allow initiating from broken chain
		for (my $job=0;$job<scalar( @{ $jobID{$fid_chainkey} });$job++) {	
		    if ( ($job == 0) && (scalar( @{ $jobID{$fid_chainkey} } )== 1) ) {#Only 1 previous jobID 
			$jobids .= ":$jobID{$fid_chainkey}[$job]"; #first and last jobID start with ":" and end without ":"
		    }
		    elsif ($job == 0 ) {
			$jobids .= ":$jobID{$fid_chainkey}[$job]:"; #first jobID start with :
		    }
		    elsif ($job eq ( scalar( @{ $jobID{$fid_chainkey} }) -1) ) {
			$jobids .= "$jobID{$fid_chainkey}[$job]"; #last jobID finish without :
		    }
		    else {
			$jobids .= "$jobID{$fid_chainkey}[$job]:";
		    }
		}
	    }
	    if ($_[3] ne "MAIN") { #Check for any previous jobIDs within MAIN path and current PATH
		my $fid_main_parallel_chainkey = "$_[1]_parallel_MAIN"; 
		if ( ($_[2] != 3) && ($jobID{$fid_main_parallel_chainkey}) ){ #If not a parallel job and a parallel job within MAIN path has prev been processed. Check if previous step was parallel and adds previous parallel jobs that have previously been submitted.
		    for (my $job=0;$job<scalar( @{ $jobID{$fid_main_parallel_chainkey} });$job++) { #All prev parallel MAIN jobIDs
			my $seen_jobIDs_counter = 0;
			if ($jobID{$fid_chainkey}) {#If any previous jobIds within current chain exists - else go ahead and add
			    for (my $current_job=0;$current_job<scalar( @{ $jobID{$fid_chainkey} });$current_job++) { #CURRENT path
				if ($jobID{$fid_chainkey}[$current_job] =~/$jobID{$fid_main_parallel_chainkey}[$job]/ ) { #Only add if not already present
				    $seen_jobIDs_counter++;
				}
			    }
			    if ($seen_jobIDs_counter eq 0) {#JobID was not present in CURRENT path
				push ( @{ $jobID{$fid_chainkey} }, $jobID{$fid_main_parallel_chainkey}[$job]); #Add jobID to hash{$}
			    }
			}
			else { #Go ahead and add
			    push ( @{ $jobID{$fid_chainkey} }, $jobID{$fid_main_parallel_chainkey}[$job]); #Add jobID to hash{$}
			}
		    }
		}
		my $fid_main_chainkey = "$_[1]_MAIN";
		if ($jobID{$fid_main_chainkey}) { #Any MAIN jobIDs necessary for broken chains, since this will be empty then
		    for (my $job=0;$job<scalar( @{ $jobID{$fid_main_chainkey} });$job++) { #Prev MAIn jobIDs
			my $seen_jobIDs_counter = 0; 
			if ($jobID{$fid_chainkey}) {#If any previous jobIds within current chain exists else go ahead and add
			    for (my $current_job=0;$current_job<scalar( @{ $jobID{$fid_chainkey} });$current_job++) { #CURRENT path
				if ($jobID{$fid_chainkey}[$current_job] =~/$jobID{$fid_main_chainkey}[$job]/) {
				    $seen_jobIDs_counter++;
				}
			    }
			    if ($seen_jobIDs_counter eq 0) {#JobID was not present in CURRENT path
				push ( @{ $jobID{$fid_chainkey} }, $jobID{$fid_main_chainkey}[$job]); #Add jobID to hash{$}
			    }
			}
			else  { #Go ahead and add
			    push ( @{ $jobID{$fid_chainkey} }, $jobID{$fid_main_chainkey}[$job]); #Add jobID to hash{$}
			}
		    }
		}
		if ($jobID{$fid_chainkey}) {
		    for (my $job=0;$job<scalar( @{ $jobID{$fid_chainkey} });$job++) {	
			if ( ($job == 0) && (scalar( @{ $jobID{$fid_chainkey} } )== 1) ) {#Only 1 previous jobID 
			    $jobids .= ":$jobID{$fid_chainkey}[$job]"; #first and last jobID start with ":" and end without ":"
			}
			elsif ($job == 0 ) {
			    $jobids .= ":$jobID{$fid_chainkey}[$job]:"; #first jobID start with :
			}
			elsif ($job eq ( scalar( @{ $jobID{$fid_chainkey} }) -1) ) {
			    $jobids .= "$jobID{$fid_chainkey}[$job]"; #last jobID finish without :
			}
			else {
			    $jobids .= "$jobID{$fid_chainkey}[$job]:";
			}
		    }
		}
	    }
	    if ($jobids) {
		$ret = `sbatch --dependency=afterok$jobids $_[4]`; #Supply with dependency of previous jobs that this one is dependent on
		($jobID) = ($ret =~ /Submitted batch job (\d+)/);
	    }
	    else {
		$ret = `sbatch $_[4]`; #No jobs have been run: submit
		($jobID) = ($ret =~ /Submitted batch job (\d+)/);
	    }
	    if ($_[2] == 1) { #Ordinary job push to array
		push ( @{ $jobID{$fid_chainkey} }, $jobID); #Add jobID to hash{$familyID}[]. Required to enable later test for all subjects sampleID_MAIN have finished before merging to family.
	    }
	    if ($_[2] == 3) { #Parallel job wait to push to array until all parallel jobs are finished within step
		push ( @{ $jobID{$fid_parallel_chainkey} }, $jobID); #Add jobID to hash{$familyID_parallel}[]. 
	    }
	}
    }
    print STDOUT "Sbatch script submitted, job id: $jobID\n";
    print STDOUT "To check status of job, please run \'jobinfo -j $jobID\'\n";
    print STDOUT "To check status of job, please run \'squeue -j $jobID\'\n";
}

sub Checkfnexists {
    
#$_[0] = complete filepath
#$_[1] = file ending
    
    my $fn;
    $fnt = 0; #Nr of sbatch with identical filenames
    for (my $i=0;$i<999;$i++) { #Number of possible files with the same name
	
	$fn = "$_[0]$i$_[1]"; #filename, filenr and fileending
	$fnt = $i; #Nr of sbatch with identical filenames, global variable
	if (-e $fn) { #if file exists 
	}
	else {
	    $i=999; #Exit loop
	}
	
    }
    $filename = $fn; #Transfer to global variable
    return;
}

sub InfileReFormat {
    
#Code needed to reformat files name removing .bam and keeping stub. 

    for (my $infile=0;$infile<scalar(@infn);$infile++) {  

	if ( basename($infn[$infile]) =~ /(\S+)\.bam$/ ) { #Parse to .bam format just filename
	   $avcovfn{$sid[$infile]} = $1; #Stores infile per sample id
	    
	}
	if ( $infn[$infile] =~ /(\S+)\.bam$/ ) { #Parse to .bam format whole path
	    
	    $infiles{$sid[$infile]} = $1;
	    $dirname{$sid[$infile]} = dirname($infn[$infile]); #Stores infile dir   
	}	
	if ( $infn[$infile] =~ /(mosaik)/ ) { #Collect what aligner was used in the first pipe to write correct path in subsequent steps
	    $aligner = $1;
	}
	if ( $infn[$infile] =~ /(bwa)/ ) { #Collect what aligner was used in the first pipe to write correct path in subsequent steps
	    $aligner = $1;
	}
    }
    return;
}

sub ReadPedigreeFile {
#Reads famid_pedigree.txt file
#IDN\tSampleID\tMother\tFather\t\Child..n
#$_[0] = filename
    
    open(PEDF, "<$_[0]") or die "Can't open $_[0]:$!, \n";    
     
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
		#print "$1-$2-$3$4", "\n"; 
		if ($3 % 2 == 1) { #Male
#% modulous operator, it gives you an integer representation of the remainder of a division operation. If X / Y divides evenly (that is to say, there's no remainder (or modulus)), the result of X % Y will be zero.
		    $pedigree{$1}{$temp[0]}[0] = "M"; #Sex, M=Male
		}
		else { #Female
		   $pedigree{$1}{$temp[0]}[0] = "F"; #Sex, F=Female
		}
		if ($4 eq "A") { #Affected
		    $pedigree{$1}{$temp[0]}[1] = 1; #1=Affected
		}
		else { #Unaffected
		    $pedigree{$1}{$temp[0]}[1] = 0; #0=Unaffected
		}
		push(@{$pedigree{$1}{$temp[0]}},@temp[2..4]); #Populate hash of array for Mother/Father/Child. Full hash: hash{FDN}{IDN}[Sex,Affected,Mother/Father/Child]
	    }
	    
	    #Old##$pedigree{$temp[0]}{$temp[1]} = [@temp[2..6]]; # Hash{fam_id}{sampleID}, all members and array [Sex,Affected,Mother,Father,Child]
	} 	
    }
    close(PEDF);
    return;
}


sub CreateModels {
#Decide genotype combinations for family depending on model, affected/healthy, and pedigree
#$_[0] = STRICT or MEDIUM or LIGHT
#STRICT requires correct PASS and correct genotype (MAF 0.005)
#MEDIUM requires either PASS or PRES and correct genotype (MAF 0.05)
#LIGTH requires either that "uncorrect" genotypes (PASS|PRES) are not called but reports all else
    
#Blank all arrays to make sure that correct model is created
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
		if ($_[0] eq "STRICT") { #Only PASS genotypes calls
		    $adom_mom = q?($_=~ /?.$sampleid.q?:PASS\:0\/1/)?; #Het, AD-2,3
		    $arecessive_mom = q?($_=~ /?.$sampleid.q?:PASS\:1\/1/)?; #Hom alt, AR-2
		    $xrecessive_mom = q?($_=~ /?.$sampleid.q?:PASS\:0\/1/)?; #Het, X-2
		    $denovo_x_mom = q?($_=~ /?.$sampleid.q?:PASS\:0\/1/)?; #Het, X_denovo-3 
		    push(@compound_aff_model,$sampleid);
		}
		if ($_[0] eq "MEDIUM") { #Either PASS or PRES genotype calls for all samples
		    $adom_mom = q?($_=~ /?.$sampleid.q?:\S+\:0\/1/)?; #Het
		    $arecessive_mom = q?($_=~ /?.$sampleid.q?:\S+\:1\/1/)?; #Hom alt
		    $xrecessive_mom = q?($_=~ /?.$sampleid.q?:\S+\:0\/1/)?; #Het
		    $denovo_x_mom = q?($_=~ /?.$sampleid.q?:\S+\:0\/1/)?; #Het
		    push(@compound_aff_model,$sampleid);
		}
		if ($_[0] eq "LIGHT") { #Only what the GenoType should not be
		    $adom_mom = q?($_!~ /?.$sampleid.q?:\S+\:(1\/1|0\/0)/)?; #Not Hom alt/ref NOTE: !~ but DP, -
		    $arecessive_mom = q?($_!~ /?.$sampleid.q?:\S+\:0\/(1|0)/)?; #Not Het/Hom ref NOTE: !~ but DP, -
		    $xrecessive_mom = q?($_!~ /?.$sampleid.q?:\S+\:(1\/1|0\/0)/)?; #Not Hom alt/ref NOTE: !~ but DP, -
		    $denovo_x_mom = q?($_!~ /?.$sampleid.q?:\S+\:(1\/1|0\/0)/)?; #Not Hom alt/ref NOTE: !~ but DP, -
		    push(@compound_aff_model,$sampleid);
		}
	    }
	    else { #healthy
		if ($_[0] eq "STRICT") {
		    $adom_mom = q?($_=~ /?.$sampleid.q?:PASS\:0\/0/)?; #Hom ref, AD-1
		    $arecessive_mom = q?($_=~ /?.$sampleid.q?:PASS\:0\/1/)?; #Het, AR-1,3
		    $xrecessive_mom = q?($_=~ /?.$sampleid.q?:PASS\:0\/(1|0)/)?; #Het/Hom ref, X-1,3,4
		    $denovo_dom_mom = q?($_=~ /?.$sampleid.q?:PASS\:0\/0/)?; #Hom ref, AD_denovo-1
		    $denovo_rec_mom = q?($_=~ /?.$sampleid.q?:PASS\:0\/(1|0)/)?; #Het/Hom ref,AR_denovo-1,2
		    $denovo_x_mom = q?($_=~ /?.$sampleid.q?:PASS\:0\/(1|0)/)?; #Het/Hom ref,X_denovo-1,2
		    push(@compound_hea_model,$sampleid);
		}
		if ($_[0] eq "MEDIUM") {
		    $adom_mom = q?($_=~ /?.$sampleid.q?:\S+\:0\/0/)?; #Hom ref
		    $arecessive_mom = q?($_=~ /?.$sampleid.q?:\S+\:0\/1/)?; #Het
		    $xrecessive_mom = q?($_=~ /?.$sampleid.q?:\S+\:0\/(1|0)/)?; #Het/Hom ref
		    $denovo_dom_mom = q?($_=~ /?.$sampleid.q?:\S+\:0\/0/)?; #Hom ref
		    $denovo_rec_mom = q?($_=~ /?.$sampleid.q?:\S+\:0\/(1|0)/)?; #Het/Hom ref
		    $denovo_x_mom = q?($_=~ /?.$sampleid.q?:\S+\:0\/(1|0)/)?; #Het/Hom ref
		    push(@compound_hea_model,$sampleid);
		}
		if ($_[0] eq "LIGHT") {
		    $adom_mom = q?($_!~ /?.$sampleid.q?:\S+\:(0|1)\/1/)?; #Not Het/Hom alt NOTE: !~
		    $arecessive_mom = q?($_!~ /?.$sampleid.q?:\S+\:(1\/1|0\/0)/)?; #Not Hom alt/ref NOTE: !~
		    $xrecessive_mom = q?($_!~ /?.$sampleid.q?:\S+\:1\/1/)?; #Not Hom alt NOTE: !~
		    $denovo_dom_mom = q?($_!~ /?.$sampleid.q?:\S+\:(0|1)\/1/)?; #Not Het/Hom alt NOTE: !~
		    $denovo_rec_mom = q?($_!~ /?.$sampleid.q?:\S+\:1\/1/)?; #Not Hom alt NOTE: !~
		    $denovo_x_mom = q?($_!~ /?.$sampleid.q?:\S+\:1\/1/)?; #Not Hom alt NOTE: !~
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
		if ($_[0] eq "STRICT") {
		    $adom_father = q?($_=~ /?.$sampleid.q?:PASS\:0\/1/)?; #Het, AD-1,3
		    $arecessive_father = q?($_=~ /?.$sampleid.q?:PASS\:1\/1/)?; #Hom alt, AR-3
		    $xrecessive_father = q?($_=~ /?.$sampleid.q?:PASS\:1\/1/)?; #Het (Male: GT is 1/1 for het on X)
		    push(@compound_aff_model,$sampleid);
		}
		if ($_[0] eq "MEDIUM") {
		    $adom_father = q?($_=~ /?.$sampleid.q?:\S+\:0\/1/)?; #Het
		    $arecessive_father = q?($_=~ /?.$sampleid.q?:\S+\:1\/1/)?; #Hom alt
		    $xrecessive_father = q?($_=~ /?.$sampleid.q?:\S+\:1\/1/)?; #Het (Male: GT is 1/1 for het on X)
		    push(@compound_aff_model,$sampleid);
		}
		if ($_[0] eq "LIGHT") {
		    $adom_father = q?($_!~ /?.$sampleid.q?:\S+\:(1\/1|0\/0)/)?; #Not Hom alt/ref NOTE: !~
		    $arecessive_father = q?($_!~ /?.$sampleid.q?:\S+\:0\/(1|0)/)?; #Not Het/Hom ref NOTE: !~
		    $xrecessive_father = q?($_!~ /?.$sampleid.q?:\S+\:0\/0/)?; #Not Hom ref NOTE: !~, (Male: GT is 1/1 for het on X) includes 0/1
		    push(@compound_aff_model,$sampleid);
		}
	    }
	    else { #healthy

		if ($_[0] eq "STRICT") {
		    $adom_father = q?($_=~ /?.$sampleid.q?:PASS\:0\/0/)?; #Hom ref, AD-2
		    $arecessive_father = q?($_=~ /?.$sampleid.q?:PASS\:0\/1/)?; #Het, AR-1,2
		    $xrecessive_father = q?($_=~ /?.$sampleid.q?:PASS\:0\/0/)?; #Hom ref
		    $denovo_dom_father = q?($_=~ /?.$sampleid.q?:PASS\:0\/0/)?; #Hom ref, AD_denovo-1
		    $denovo_rec_father = q?($_=~ /?.$sampleid.q?:PASS\:0\/(1|0)/)?; #Het/Hom ref, AR_denovo-1,2
		    $denovo_x_father = q?($_=~ /?.$sampleid.q?:PASS\:0\/0/)?; #Hom reference, X_denovo-1,2,3
		    push(@compound_hea_model,$sampleid);
		}
		if ($_[0] eq "MEDIUM") {
		    $adom_father = q?($_=~ /?.$sampleid.q?:\S+\:0\/0/)?; #Hom ref
		    $arecessive_father = q?($_=~ /?.$sampleid.q?:\S+\:0\/1/)?; #Het
		    $xrecessive_father = q?($_=~ /?.$sampleid.q?:\S+\:0\/0/)?; #Hom ref
		    $denovo_dom_father = q?($_=~ /?.$sampleid.q?:\S+\:0\/0/)?; #Hom ref
		    $denovo_rec_father = q?($_=~ /?.$sampleid.q?:\S+\:0\/(1|0)/)?; #Het/Hom ref, AR_denovo-1,2
		    $denovo_x_father = q?($_=~ /?.$sampleid.q?:\S+\:0\/0/)?; #Hom reference
		    push(@compound_hea_model,$sampleid);
		}
		if ($_[0] eq "LIGHT") {
		    $adom_father = q?($_!~ /?.$sampleid.q?:\S+\:(0|1)\/1/)?; #Not Het/Hom alt NOTE: !~
		    $arecessive_father = q?($_!~ /?.$sampleid.q?:\S+\:(1\/1|0\/0)/)?; #Not Hom alt/ref NOTE: !~
		    $xrecessive_father = q?($_!~ /?.$sampleid.q?:\S+\:1\/1/)?; #Not Hom alt NOTE: !~ DP, - includes 0/1
		    $denovo_dom_father = q?($_!~ /?.$sampleid.q?:\S+\:(0|1)\/1/)?; #Not Het/Hom alt NOTE: !~
		    $denovo_rec_father = q?($_!~ /?.$sampleid.q?:\S+\:1\/1/)?; #Not Hom alt NOTE: !~
		    $denovo_x_father = q?($_!~ /?.$sampleid.q?:\S+\:1\/1/)?; #Not Hom alt NOTE: !~ DP, - includes 0/1
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

		if ($_[0] eq "STRICT") {
		    $arecessive_child = q?($_=~ /?.$sampleid.q?:PASS\:1\/1/)?; #Hom alt, AR-1,2,3
		    $denovo_dom_child = q?($_=~ /?.$sampleid.q?:PASS\:0\/1/)?; #Het, AD_denovo-1
		    $denovo_rec_child = q?($_=~ /?.$sampleid.q?:PASS\:1\/1/)?; #Hom alt, AR_denovo-1,2
		    if ( ($fatherDS == 1) && ($motherDS == 1) ) { #Parents affected
			$adom_child = q?($_=~ /?.$sampleid.q?:PASS\:(0|1)\/1/)?; #Het/hom alt, AD-3 (Rare)
		    }
		    else {
			$adom_child = q?($_=~ /?.$sampleid.q?:PASS\:0\/1/)?; #Het, AD-1,2 (Common)
		    }
		    push(@compound_aff_model,$sampleid);
		    if ($pedigree{$familyid}{$sampleid}[0]=~/M/i) { #Affected and Male
			$xrecessive_child = q?($_=~ /?.$sampleid.q?:PASS\:1\/1/)?; #Het (Male: GT is 1/1 for het on X), X-1,2,3,4	    
			$denovo_x_child = q?($_=~ /?.$sampleid.q?:PASS\:1\/1/)?; #Het (Male: GT is 1/1 for het on X), X_denovo-1,2,3		    
		    }
		    if ($pedigree{$familyid}{$sampleid}[0]=~/F/i) { #Affected and Female
			$xrecessive_child = q?($_=~ /?.$sampleid.q?:PASS\:(0|1)\/1/)?; #Het/hom alt, X-1,2,3,4
			$denovo_x_child = q?($_=~ /?.$sampleid.q?:PASS\:(0|1)\/1/)?; #Het/hom alt (Hom more unlikely) ), X_denovo-1,2,3
		    }
		}
		if ($_[0] eq "MEDIUM") {
		    $arecessive_child = q?($_=~ /?.$sampleid.q?:\S+\:1\/1/)?; #Hom alt
		    $denovo_dom_child = q?($_=~ /?.$sampleid.q?:\S+\:0\/1/)?; #Het
		    $denovo_rec_child = q?($_=~ /?.$sampleid.q?:\S+\:1\/1/)?; #Hom alt
		    if ( ($fatherDS == 1) && ($motherDS == 1) ) { #Parents affected
			$adom_child = q?($_=~ /?.$sampleid.q?:\S+\:(0|1)\/1/)?; #Het/hom alt
		    }
		    else {
			$adom_child = q?($_=~ /?.$sampleid.q?:\S+\:0\/1/)?; #Het, AD-1,2 (Common)
		    }
		    push(@compound_aff_model,$sampleid);
		    if ($pedigree{$familyid}{$sampleid}[0]=~/M/i) { #Affected and Male			
			$xrecessive_child = q?($_=~ /?.$sampleid.q?:\S+\:1\/1/)?; #Het (Male: GT is 1/1 for het on X)		    
			$denovo_x_child = q?($_=~ /?.$sampleid.q?:\S+\:1\/1/)?; #Het (Male: GT is 1/1 for het on X)
		    }
		    if ($pedigree{$familyid}{$sampleid}[0]=~/F/i) { #Affected and Female
			$xrecessive_child = q?($_=~ /?.$sampleid.q?:\S+\:(0|1)\/1/)?; #Het/hom alt
			$denovo_x_child = q?($_=~ /?.$sampleid.q?:\S+\:(0|1)\/1/)?; #Het/hom alt (Hom more unlikely) )
		    }
		}
		if ($_[0] eq "LIGHT") {
		    $arecessive_child = q?($_!~ /?.$sampleid.q?:\S+\:0\/(1|0)/)?; #Not Het/Hom ref NOTE: !~
		    $denovo_dom_child = q?($_!~ /?.$sampleid.q?:\S+\:(1\/1|0\/0)/)?; #Not Hom alt/ref NOTE: !~
		    $denovo_rec_child = q?($_!~ /?.$sampleid.q?:\S+\:0\/(1|0)/)?; #Not Het/Hom ref NOTE: !~
		    if ( ($fatherDS == 1) && ($motherDS == 1) ) { #Parents affected
			$adom_child = q?($_!~ /?.$sampleid.q?:\S+\:0\/0/)?; #Not Hom ref
		    }
		    else {
			$adom_child = q?($_!~ /?.$sampleid.q?:\S+\:(1\/1|0\/0)/)?; #Not Hom alt/ref NOTE: !~
		    }
		    push(@compound_aff_model,$sampleid);
		    if ($pedigree{$familyid}{$sampleid}[0]=~/M/i) { #Affected and Male
			$xrecessive_child = q?($_!~ /?.$sampleid.q?:\S+\:0\/0/)?; #Not Hom ref (Male: GT is 1/1 for het on X) includes 0/1
			$denovo_x_child = q?($_!~ /?.$sampleid.q?:\S+\:0\/0/)?; #Not Hom ref (Male: GT is 1/1 for het on X) includes 0/1
		    }
		    if ($pedigree{$familyid}{$sampleid}[0]=~/F/i) { #Affected and Female
			$xrecessive_child = q?($_!~ /?.$sampleid.q?:\S+\:0\/0/)?; #Not Hom ref
			$denovo_x_child = q?($_!~ /?.$sampleid.q?:\S+\:0\/0/)?; #Not Hom ref
		    }
		}
	    }
	    else { #healthy

		if ($_[0] eq "STRICT") {
		    $adom_child = q?($_=~ /?.$sampleid.q?:PASS\:0\/0/)?; #Hom ref, AD-1,2,3
		    $arecessive_child = q?($_=~ /?.$sampleid.q?:PASS\:0\/(1|0)/)?; #Het/Hom ref, AR-1,2,3
		    $denovo_dom_child = q?($_=~ /?.$sampleid.q?:PASS\:0\/0/)?; #Hom ref, AD_denovo-1
		    $denovo_rec_child = q?($_=~ /?.$sampleid.q?:PASS\:0\/(1|0)/)?; #Het/Hom ref, AD_denovo-1,2
		    push(@compound_hea_model,$sampleid);
		    if ( $pedigree{$familyid}{$sampleid}[0]=~/M/i) { #Healthy and Male			
			$xrecessive_child = q?($_=~ /?.$sampleid.q?:PASS\:0\/0/)?; #Hom ref (Male: GT is 1/1 for het on X), X-1,2,3,4
			$denovo_x_child = q?($_=~ /?.$sampleid.q?:PASS\:0\/0/)?; #Hom ref (Male: GT is 1/1 for het on X), X_denovo-1,2,3
		    }
		    if ($pedigree{$familyid}{$sampleid}[0]=~/F/i) { #Healthy and Female
			$xrecessive_child = q?($_=~ /?.$sampleid.q?:PASS\:0\/(1|0)/)?; #Het/Hom ref, X-1,2,3,4
			$denovo_x_child = q?($_=~ /?.$sampleid.q?:PASS\:0\/(1|0)/)?; #Het/Hom ref, X_denovo-1,2,3
		    }
		}
		if ($_[0] eq "MEDIUM") {
		    $adom_child = q?($_=~ /?.$sampleid.q?:\S+\:0\/0/)?; #Hom ref
		    $arecessive_child = q?($_=~ /?.$sampleid.q?:\S+\:0\/(1|0)/)?; #Het/Hom ref
		    $denovo_dom_child = q?($_=~ /?.$sampleid.q?:\S+\:0\/0/)?; #Hom ref
		    $denovo_rec_child = q?($_=~ /?.$sampleid.q?:\S+\:0\/(1|0)/)?; #Het/Hom ref, AD_denovo-1,2
		    push(@compound_hea_model,$sampleid);
		    if ( $pedigree{$familyid}{$sampleid}[0]=~/M/i) { #Healthy and Male
			$xrecessive_child = q?($_=~ /?.$sampleid.q?:\S+\:0\/0/)?; #Hom (Male: GT is 1/1 for het on X)
			$denovo_x_child = q?($_=~ /?.$sampleid.q?:\S+\:0\/0/)?; #Hom (Male: GT is 1/1 for het on X)
		    }
		    if ($pedigree{$familyid}{$sampleid}[0]=~/F/i) { #Healthy and Female
			$xrecessive_child = q?($_=~ /?.$sampleid.q?:\S+\:0\/(1|0)/)?; #Het/Hom ref
			$denovo_x_child = q?($_=~ /?.$sampleid.q?:\S+\:0\/(1|0)/)?; #Het/Hom ref
		    }
		}
		if ($_[0] eq "LIGHT") {
		    $adom_child = q?($_!~ /?.$sampleid.q?:\S+\:(0|1)\/1/)?; #Not Het/Hom alt NOTE: !~
		    $arecessive_child = q?($_!~ /?.$sampleid.q?:\S+\:1\/1/)?; #Not Hom alt NOTE: !~
		    $denovo_dom_child = q?($_!~ /?.$sampleid.q?:\S+\:(0|1)\/1/)?; #Not Het/hom alt NOTE: !~
		    $denovo_rec_child = q?($_!~ /?.$sampleid.q?:\S+\:1\/1/)?; #Not Hom alt
		    push(@compound_hea_model,$sampleid);
		    if ( $pedigree{$familyid}{$sampleid}[0]=~/M/i) { #Healthy and Male
			$xrecessive_child = q?($_!~ /?.$sampleid.q?:\S+\:1\/1/)?; #Not Het (Male: GT is 1/1 for het on X) includes 0/1
			$denovo_x_child = q?($_!~ /?.$sampleid.q?:\S+\:1\/1/)?; #Not Het (Male: GT is 1/1 for het on X) includes 0/1
		    }
		    if ($pedigree{$familyid}{$sampleid}[0]=~/F/i) { #Healthy and Female
			$xrecessive_child = q?($_!~ /?.$sampleid.q?:\S+\:1\/1/)?; #Not Hom alt
			$denovo_x_child = q?($_!~ /?.$sampleid.q?:\S+\:1\/1/)?; #Not Hom alt
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
    
#Autosomal Dominant (AD)
    $adom_model =q?perl -nae' if ( ?;
    for (my $samples=0;$samples<scalar(@adom_model);$samples++) {
	    if ($samples == scalar(@adom_model)-1) { #Ensure that last print occurs without &&
		$adom_model .= $adom_model[$samples];
	    }
	    else {
		$adom_model .= "$adom_model[$samples] && ";
	    }
    }
    $adom_model .=q? ) { print $_;}'?;    
    
#Autosomal Recessive (AR)
    $arecessive_model =q?perl -nae' if ( ?;
    for (my $samples=0;$samples<scalar(@arecessive_model);$samples++) {
	if ($samples == scalar(@arecessive_model)-1) { #Ensure that last print occurs without &&
	    $arecessive_model .= $arecessive_model[$samples];
	}
	else {
	    $arecessive_model .= "$arecessive_model[$samples] && ";
	}
    }
    $arecessive_model .=q? ) { print $_;}'?;    
    
#X-linked (X)
    $xrecessive_model =q?perl -nae' if ( ?;
    for (my $samples=0;$samples<scalar(@xrecessive_model);$samples++) {
	if ($samples == scalar(@xrecessive_model)-1) { #Ensure that last print occurs without &&
	    $xrecessive_model .= $xrecessive_model[$samples];
	}
	else {
	    $xrecessive_model .= "$xrecessive_model[$samples] && ";
	}
    }
    $xrecessive_model .=q? ) { print $_;}'?;
    
#De novo Dominant (AD_denovo)
    $denovo_dom_model =q?perl -nae' if ( ?;
    for (my $samples=0;$samples<scalar(@denovo_dom_model);$samples++) {
	if ($denovo_dom_model[$samples]) { #If entry exist i.e subject is entered

	    if ($samples == scalar(@denovo_dom_model)-1) { #Ensure that last print occurs without &&
		$denovo_dom_model .= $denovo_dom_model[$samples];
	    }
	    else {
		$denovo_dom_model .= "$denovo_dom_model[$samples] && ";
	    }
	}
	else { #No entry for subject i.e. mother or father affected. Do not use subject in filtering
	    if ($samples == scalar(@denovo_dom_model)-1) { #Ensure that last print occurs without &&
		$denovo_dom_model .= '$_'; #If no subjects match model i.e. mother or father is affected and no other members have been sequenced. $_ = No segregation analysis (in this case)
	    }
	    else {
		$denovo_dom_model .= '$_ && ';
	    }
	}
    }
    $denovo_dom_model .=q? ) { print $_;}'?;
    
#De novo Recessive (AR_denovo)
    $denovo_rec_model =q?perl -nae' if ( ?;
    for (my $samples=0;$samples<scalar(@denovo_rec_model);$samples++) {
	if ($denovo_rec_model[$samples]) { #If entry exist i.e subject is entered
	    
	    if ($samples == scalar(@denovo_rec_model)-1) { #Ensure that last print occurs without &&
		$denovo_rec_model .= $denovo_rec_model[$samples];
	    }
	    else {
		$denovo_rec_model .= "$denovo_rec_model[$samples] && ";
	    }
	}
	else { #No entry for subject i.e. mother or father affected. Do not use subject in filtering
	    if ($samples == scalar(@denovo_rec_model)-1) { #Ensure that last print occurs without &&
		$denovo_rec_model .= '$_'; #If no subjects match model i.e. mother or father is affected and no other members have been sequenced. $_ = No segregation analysis (in this case)
	    }
	    else {
		$denovo_rec_model .= '$_ && ';
	    }
	}
    }
    $denovo_rec_model .=q? ) { print $_;}'?;
    
#De novo X (X_denovo)
    $denovo_x_model =q?perl -nae' if ( ?;
    for (my $samples=0;$samples<scalar(@denovo_x_model);$samples++) {
	if ($denovo_x_model[$samples]) { #If entry exist i.e subject is entered

	    if ($samples == scalar(@denovo_x_model)-1) { #Ensure that last print occurs without &&
		$denovo_x_model .= $denovo_x_model[$samples];
	    }
	    else {
		$denovo_x_model .= "$denovo_x_model[$samples] && ";
	    }
	}
	else { #No entry for subject i.e. mother or father affected. Do not use subject in filtering
	    if ($samples == scalar(@denovo_x_model)-1) { #Ensure that last print occurs without &&
		$denovo_x_model .= '$_'; #If no subjects match model i.e. mother or father is affected and no other members have been sequenced. $_ = No segregation analysis (in this case)
	    }
	    else {
		$denovo_x_model .= '$_ && ';
	    }
	}
    }
    $denovo_x_model .=q? ) { print $_;}'?;
    
#Genetic Compound (AR_compound)
#Affected samplesIDs
    for (my $samples=0;$samples<scalar(@compound_aff_model);$samples++) {
	if ($samples == scalar(@compound_aff_model)-1) { #Ensure that last print occurs without ,
	    $comp_aff_sampleid .= $compound_aff_model[$samples];
	}
	else {
	    $comp_aff_sampleid .= "$compound_aff_model[$samples],";
	}
    }
    
#Healthy sampleIDs
    if (@compound_hea_model) {
	for (my $samples=0;$samples<scalar(@compound_hea_model);$samples++) {
	    if ($samples == scalar(@compound_hea_model)-1) { #Ensure that last print occurs without ,
		$comp_hea_sampleid .= $compound_hea_model[$samples];
	    }
	    else {
		$comp_hea_sampleid .= "$compound_hea_model[$samples],";
	    }
	}
    }
    else {
	$comp_hea_sampleid = "";
	}
    return;
}
