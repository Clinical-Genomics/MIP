#!/usr/bin/perl -w

use strict;
use warnings;

###Master script for analysing paired end reads from the Illumina plattform in fastq(.gz) format to annotated ranked disease causing variants. The program performs QC, aligns reads using Mosaik or BWA, performs variant discovery and annotation as well as ranking the found variants according to disease potential.
 
###Copyright 2011 Henrik Stranneheim
    
use Pod::Usage;
use Pod::Text;
use Getopt::Long;
use POSIX;
use IO::File;
use YAML;

use vars qw($USAGE);

BEGIN {
    $USAGE =
	qq{
mip.pl  -ifd [inFilesDirs,.,.,.,n] -isd [inScriptDir,.,.,.,n] -rd [refdir] -p [project ID] -s [sample ID...n] -em [e-mail] -osd [outdirscripts] -odd [outDataDir] -f [familyID] -p[program]
               ####MIP
	       -ifd/--inFilesDirs Infile directory(s), comma sep (Mandatory: Supply whole path,)
               -isd/--inScriptDir The pipeline custom script in directory (Mandatory: Supply whole path)
               -rd/--referencesDir Reference(s) directory (Mandatory: Supply whole path)
	       -p/--projectID The project ID  (Mandatory)
	       -s/--sampleIDs The sample ID(s),comma sep (Mandatory)
	       -em/--email E-mail
	       -odd/--outDataDir The data files output directory (Mandatory: Supply whole path)
	       -osd/--outScriptDir The script files (.sh) output directory (Mandatory: Supply whole path)
               -f/--familyID Group id of samples to be compared (defaults to "0" (=no), (Ex: 1 for IDN 1-1-1A))
               -pedigree/--pedigreeFile (Supply whole path, defaults to "")
               -huref/--humanGenomeReference Fasta file for the human genome reference (defaults to "")
               -al/--aligner Setting which aligner was used for alignment in previous analysis (defaults to "")
               -at/--analysisType Type of analysis to perform (defaults to "exomes";Valid entries: "genomes", "exomes", "rapid")
               -mc/--maximumCores The maximum number of cores per node used in the analysis (defaults to "8")
               -c/--configFile YAML config file for script parameters (defaults to "")
               -wc/--writeConfigFile Write YAML configuration file for script parameters (defaults to "";Supply whole path)
               -si/--sampleInfoFile YAML file for sample info used in the analysis (defaults to "{outDataDir}/{familyID}/{familyID}_qc_sampleInfo.yaml")
               -dra/--dryRunAll Sets all programs to dry run mode i.e. no sbatch submission (defaults to "0" (=no))
               -julp/--javaUseLargePages Use large page memory. (-XX,hence option considered not stable and are subject to change without notice, but can be consiered when faced with Java Runtime Environment Memory issues)
               -pve/--pythonVirtualEnvironment Pyhton virtualenvironment (defaults to "")
               -h/--help Display this help message    
               -v/--version Display version of MIP            

               
               ####Programs
               -pGZ/--pGZip GZip fastq files (defaults to "1" (=yes))
	       -pFQC/--pFastQC Sequence quality analysis using FastQC (defaults to "1" (=yes))
               -pREM/--pRemovalRedundantFiles Generating sbatch script for deletion of redundant files (defaults to "1" (=yes);Note: Must be submitted manually)
               
               ##Mosaik
	       -pMoB/--pMosaikBuild  Convert reads to Mosaik format using MosaikBuild (defaults to "1" (=yes))
                -mobmfl/--mosaikBuildMedianFragLength Flag for setting the mean fragment length, mfl, (defaults to (=375) bp)
	       -pMoA/--pMosaikAlign Align reads using MosaikAlign (defaults to "1" (=yes))
                 -moaref/--mosaikAlignReference MosaikAlign reference (defaults to "{humanGenomeReference}")
                 -moaannpe/--mosaikAlignNeuralNetworkPeFile MosaikAlign Neural Network PE File (defaults to "2.1.78.pe.ann")
                 -moaannse/--mosaikAlignNeuralNetworkSeFile MosaikAlign Neural Network SE File (defaults to "2.1.78.se.ann")
                 -mojdb/--mosaikJumpDbStub MosaikJump stub (defaults to "{humanGenomeReference}")
               
               ##BWA
               -pBWA_mem/--pBwaMem Align reads using BWA Mem (defaults to "0" (=no))
                 -bwamemrdb/--bwaMemRapidDb Selection of relevant regions post alignment (Defaults to "")
               -pBWA_aln/--pBwaAln Index reads using BWA Aln (defaults to "0" (=no))
                 -bwaalnq/--bwaAlnQualityTrimming BWA Aln quality threshold for read trimming (defaults to "20")
               -pBWA_sampe/--pBwaSampe Align reads using BWA Sampe (defaults to "0" (=no))
               
               ##PicardTools
               -picardpath/--picardToolsPath Path to PicardTools. Mandatory for use of PicardTools (defaults to "")
               -picttmpd/--PicardToolsTempDirectory Temporary Directory to write to using PicardTools (defaults to "/scratch/SLURM_JOB_ID";Supply whole path)
               -pPicT_sort/--pPicardToolsSortSam Sort & index aligned reads using PicardTools SortSam & index (defaults to "1" (=yes))
               -pPicT_merge/--pPicardToolsMergeSamFiles Merge (BAM file(s) ) using PicardTools MergeSamFiles (defaults to "1" (=yes))
               -pPicT_mergerr/--pPicardToolsMergeRapidReads Merge Read batch processed (BAM file(s)) using PicardTools MergeSamFiles (Only relevant in rapid mode;defaults to "0" (=no))
                 -picT_mergeprev/--picardToolsMergeSamFilesPrevious PicardTools MergeSamFiles on merged current files and previous BAM-file(s) (Supply whole path and name, name must contain sample id, and lanes_Xn info)
               -pPicT_markdup/--pPicardToolsMarkduplicates Markduplicates using PicardTools MarkDuplicates (defaults to "1" (=yes))
               
               ##Coverage Calculations
               -pCh_B/--pChanjoBuild Chanjo build central SQLite database file (defaults to "1" (=yes))
                 -chbdb/--chanjoBuildDb  Reference database (defaults to "")
               -pCh_C/--pChanjoCalculate Chanjo coverage analysis (defaults to "1" (=yes))
                 -chccut/--chanjoCalculateCutoff Read depth cutoff (defaults to "10")
               -pCh_I/--pChanjoImport Chanjo import to collect sample info to family Db  (defaults to "1" (=yes))
               -pCC_bedgc/--pGenomeCoverageBED Genome coverage calculation using genomeCoverageBED (defaults to "1" (=yes))
               -pCC_qac/--pQaCompute Genome coverage calculation using qaCompute (defaults to "1" (=yes))
               -xcov/--xCoverage Max coverage depth when using '-pGenomeCoverageBED', '-pQaCompute' (defaults to "30")
               -pCC_picmm/--pPicardToolsCollectMultipleMetrics Metrics calculation using PicardTools collectMultipleMetrics (defaults to "1" (=yes))
               -pCCE_pichs/--pPicardToolsCalculateHSMetrics Capture calculation using PicardTools CalculateHSmetrics (defaults to "1" (=yes))
                 -extbl/--exomeTargetBedInfileLists Prepared target BED file for PicardTools CalculateHSMetrics (defaults to "". File ending should be ".infile_list") 
                 -extpbl/--exomeTargetPaddedBedInfileLists Prepared padded target BED file for PicardTools CalculateHSMetrics (defaults to "". File ending should be ".padXXX.infile_list")
               -pRCP/--pRCovPlots Plots of genome coverage using rCovPlots (defaults to "1" (=yes))
               
               ##GATK              
               -gatkpath/--genomeAnalysisToolKitPath  Path to GATK. Mandatory for use of GATK (defaults to "")
               -gatktmpd/--GATKTempDirectory Temporary Directory to write to using GATK ReAlignerTargetCreator & BaseRecalibrator (defaults to "/scratch/SLURM_JOB_ID";Supply whole path)
               -gatktpbl/--GATKTargetPaddedBedIntervalLists Target BED file interval for GATK (defaults to "". File ending should be ".padXXX.interval_list")
               -gatkdcov/--GATKDownSampleToCoverage Coverage to downsample to at any given locus (defaults to "1000")
               -pGATK_real/--pGATKRealigner Realignments of reads using GATK realign (defaults to "1" (=yes))
                 -gatkrealknset1/--GATKReAlignerINDELKnownSet1 GATK ReAlignerTargetCreator/IndelRealigner known INDEL set 1 (defaults to "1000G_phase1.indels.hg19.vcf")
                 -gatkrealknset2/--GATKReAlignerINDELKnownSet2 GATK ReAlignerTargetCreator/IndelRealigner known INDEL set 2 (defaults to "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf")
               -pGATK_baserecal/--pGATKBaseRecalibration Recalibration of bases using GATK BaseRecalibrator/PrintReads (defaults to "1" (=yes))
                 -gatkbaserecalknset/--GATKBaseReCalibrationSNPKnownSet GATK BaseReCalinbration known SNP set (defaults to "dbsnp_135.b37.vcf") 
               -pGATK_rr/--pGATKReduceReads Reduces reads in BAM file (defaults to "1" (=yes))                  
               -pGATK_hapcall/--pGATKHaploTypeCaller Variant discovery using GATK HaplotypeCaller (defaults to "1" (=yes))
                 -gatkhapcallrefbaminfile/--GATKHaploTypeCallerRefBAMInfile GATK HaplotypeCaller BAM reference infile list for joint genotyping (defaults to "")
                 -gatkhapcallsnpknset/--GATKHaploTypeCallerSNPKnownSet GATK HaplotypeCaller dbSNP set for annotating ID columns (defaults to "dbsnp_135.b37.vcf")
               -pGATK_hapcallcombine/--pGATKHaploTypeCallerCombineVariants Combine variants from HaplotypeCaller (defaults to "1" (=yes))
               -pGATK_varrecal/--pGATKVariantRecalibration Variant recalibration using GATK VariantRecalibrator/ApplyRecalibration (defaults to "1" (=yes))
                 -gatkexrefsnp/--GATKExomeReferenceSNPs Prepared exome reference file (SNVs) for GATKVariantRecalibration (defaults to "")
                 -gatkvarrecaltrhapmap/--GATKVariantReCalibrationTrainingSetHapMap GATK VariantRecalibrator HapMap training set (defaults to "hapmap_3.3.b37.sites.vcf")
                 -gatkvarrecaltrdbsnp/--GATKVariantReCalibrationTrainingSetDbSNP GATK VariantRecalibrator dbSNP training set (defaults to "dbsnp_135.b37.vcf")
                 -gatkvarrecaltrd1000Gsnp/--GATKVariantReCalibrationTrainingSet1000GSNP GATK VariantRecalibrator 1000G high confidence SNP training set (defaults to "1000G_phase1.snps.high_confidence.b37.vcf")
                 -gatkvarrecaltromni/--GATKVariantReCalibrationTrainingSet1000GOmni GATK VariantRecalibrator 1000G_omni training set (defaults to "1000G_omni2.5.b37.sites.vcf")
                 -gatkvarrecaltrdbmills/--GATKVariantReCalibrationTrainingSetMills GATK VariantRecalibrator Mills training set (defaults to "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf")
                 -gatkvarrecaltsfilterlevel/--GATKVariantReCalibrationTSFilterLevel The truth sensitivity level at which to start filtering used in GATK VariantRecalibrator (defaults to "99.9")
               -pGATK_phaseTr/--pGATKPhaseByTransmission Computes the most likely genotype and phases calls were unamibigous using GATK PhaseByTransmission (defaults to "1" (=yes))
               -pGATK_readPh/--pGATKReadBackedPhasing Performs physical phasing of SNP calls, based on sequencing reads using GATK ReadBackedPhasing (defaults to "1" (=yes))
                 -gatkreadphphaseqthr/--GATKReadBackedPhasingPhaseQualityThresh The minimum phasing quality score required to output phasing
               -pGATK_varevalall/--pGATKVariantEvalAll Variant evaluation using GATK VariantEval for all variants  (defaults to "1" (=yes))
               -pGATK_varevalexome/--pGATKVariantEvalExome Variant evaluation using GATK VariantEval for exonic variants  (defaults to "1" (=yes))
                 -gatkvarevaldbsnp/--GATKVariantEvalDbSNP DbSNP file used in GATK VariantEval (defaults to "")
                 -gatkvarevaldbgold/--GATKVariantEvalGold Gold Indel file used in GATK VariantEval (defaults to "")
               
               ##ANNOVAR
               -pANVAR/--pAnnovar Annotate variants using Annovar (defaults to "1" (=yes))
                 -anvarpath/--annovarPath  Path to Annovar script directory (Supply whole path, defaults to "". NOTE: Assumes that the annovar db files are located in annovar/humandb)
                 -anvargbv/--annovarGenomeBuildVersion Annovar genome build version (defaults to "hg19")
                 -anvartn/--annovarTableNames Annovar table names (comma sep)
                 -anvarstn/--annovarSupportedTableNames Print Annovar MIP supported table names (defaults 0 (=no))
                 -anvarmafth/--annovarMAFThreshold Sets the minor allele frequency threshold in annovar (defaults to "0")
               
               ##VMerge  
               -pMerge_anvar/--pMergeAnnotatedVariants Merge (& annotate) all annotated variants into one file using intersectCollect.pl to  (defaults to "1" (=yes))
                 -mergeanvarte/--mergeAnnotatedVariantsTemplateFile Db template file used to create the specific family '-mergeanvardbf' master file (defaults to "")
                 -mergeanvardbf/--mergeAnnotatedVariantsDbFile Db master file to be used in intersectCollect.pl (defaults to  "{outDataDir}/{familyID}/{familyID}_intersectCollect_db_master.txt";Supply whole path)

               ##Add_depth
               -pAdd_dp/--pAddDepth Adds read depth at nonvariant sites using SamTools mpileup and add_depth.pl (defaults to "1" (=yes))

               ##RankVariants
               -pRankVar/--pRankVariants Ranking of annotated variants (defaults to "1" (=yes))
                 -rs/--rankScore The rank score cut-off (defaults to "-100", .i.e. include everything
                 -imdbfile/--ImportantDbFile Important Db file (Defaults to "")
                 -imdbte/--ImportantDbTemplate Important Db template file used to create the specific family '-im_dbmf' master file (Defaults to "")
                 -imdbmf/--ImportantDbMasterFile Important Db master file to be used when selecting variants (defaults to "{outDataDir}/{familyID}/{familyID}.intersectCollect_selectVariants_db_master.txt";Supply whole path) 
                 -imdbfof/--ImportantDbFileOutFiles The file(s) to write to when selecting variants with intersectCollect.pl. Comma sep (defaults to "{outDataDir}/{familyID}/{aligner}/GATK/candidates/ranking/{familyID}_orphan.selectVariants, {outDataDir}/{familyID}/{aligner}/GATK/candidates/ranking/clinical/{familyID}.selectVariants"; Supply whole path/file)
               -pSCheck/--pSampleCheck QC for samples gender and relationship (defaults to "1" (=yes) )
               -pQCC/--pQCCollect Collect QC metrics from programs processed (defaults to "1" (=yes) )
                 -QCCsampleinfo/--QCCollectSampleInfoFile SampleInfo File containing info on what to parse from this analysis run (defaults to "{outDataDir}/{familyID}/{familyID}_qc_sampleInfo.yaml")
                 -QCCregexp/--QCCollectRegExpFile Regular expression file containing the regular expression to be used for each program (defaults to "")
	   };
}


####Script parameters

my %parameter; #Holds all parameters for MIP
my %scriptParameter; #Holds all active parameters after the value has been set

$scriptParameter{'MIP'} = 1; #Enable/activate MIP 

my @orderParameters; #To add/write parameters in the correct order

#Add timestamp for later use in mip_log and qcmetrics yaml file
my $timeStamp = (`date +%Y%m%d_%Hh%Mm`); #Catches current date and script name
chomp($timeStamp); #Remove \n;

####Set program parameters

###Project specific
##DefineParameters
##parameterName, parameterType, parameterDefault, AssociatedProgram, Check directory/file existence, parameterChain, programCheck)
##DefineParametersPath
##parameterName, parameterDefault, AssociatedProgram, Check directory/file existence, File Autovivication)

&DefineParameters("projectID", "MIP", "nodefault", "MIP");

&DefineParameters("email", "MIP", 0, "MIP");

&DefineParametersPath("familyID", "nodefault", "MIP", 0);

&DefineParameters("maximumCores", "MIP", 16, "MIP");

&DefineParametersPath("configFile", 0, "MIP", "file");

&DefineParameters("analysisType", "MIP", "exomes", "MIP");

&DefineParametersPath("outDataDir", "nodefault", "MIP", 0);

&DefineParametersPath("outScriptDir", "nodefault", "MIP", 0);

&DefineParametersPath("writeConfigFile", 0, "MIP", 0);

&DefineParametersPath("pedigreeFile", "nodefault", "MIP", "file", "noAutoBuild");

&DefineParametersPath("sampleInfoFile", "NotsetYet", "MIP", "file", "noAutoBuild");

&DefineParametersPath("inScriptDir", "nodefault", "MIP", "directory");

&DefineParametersPath("referencesDir", "nodefault", "MIP", "directory");

&DefineParameters("dryRunAll", "MIP", 0, "MIP");

my (@inFilesDirs,@sampleIDs); #Arrays for input file directorys,sampleIDs

###Programs

##GZip
&DefineParameters("pGZip", "program", 1, "MIP", "nofileEnding", "MAIN", "gzip");


##FastQC
&DefineParameters("pFastQC", "program", 1, "MIP", "nofileEnding", "RawSeqQC", "fastqc");


##RemovalRedundantFiles
&DefineParameters("pRemovalRedundantFiles", "program", 1, "MIP", "nofileEnding", "MAIN");


##Mosaik
&DefineParameters("pMosaikBuild", "program", 1, "MIP", "nofileEnding", "MAIN", "MosaikBuild");

&DefineParameters("mosaikBuildMedianFragLength", "program", 375, "pMosaikBuild");

&DefineParameters("pMosaikAlign", "program", 1, "MIP", "nofileEnding", "MAIN", "MosaikAligner");

&DefineParametersPath("mosaikAlignReference", "notSetYet", "pMosaikAlign", "file", "yesAutoBuild");

&DefineParametersPath("mosaikAlignNeuralNetworkPeFile", "2.1.78.pe.ann", "pMosaikAlign", "file", "yesAutoBuild");

&DefineParametersPath("mosaikAlignNeuralNetworkSeFile", "2.1.78.se.ann", "pMosaikAlign", "file", "yesAutoBuild");

&DefineParametersPath("mosaikJumpDbStub", "notSetYet", "pMosaikAlign", "file", "yesAutoBuild");
my @mosaikJumpDbStubFileEndings = ("_keys.jmp", "_meta.jmp", "_positions.jmp");

##BWA

&DefineParameters("pBwaMem", "program", 0, "MIP", "nofileEnding", "MAIN", "bwa");

&DefineParametersPath("bwaMemRapidDb", "nodefault", "pBwaMem", "file", "noAutoBuild");

&DefineParameters("pBwaAln", "program", 0, "MIP", "nofileEnding", "MAIN", "bwa");

&DefineParameters("bwaAlnQualityTrimming", "program", 20, "pBwaAln");

&DefineParameters("pBwaSampe", "program", 0, "MIP", "nofileEnding", "MAIN", "bwa");

&DefineParametersPath("bwaBuildReference", "notSetYet", "pBwaMem,pBwaAln,pBwaSampe", "file", "yesAutoBuild");
my @bwaBuildReferenceFileEndings = (".amb", ".ann", ".bwt", ".pac", ".sa");

##Choosen MIP Aligner

&DefineParameters("aligner", "MIP", "mosaik", "MIP");

##PicardTools

&DefineParameters("pPicardToolsSortSam", "program", 1, "MIP", "_sorted", "MAIN");

&DefineParameters("pPicardToolsMergeRapidReads", "program", 0, "MIP", "_sorted", "MAIN");#Rapid mode special case

&DefineParameters("pPicardToolsMergeSamFiles", "program", 1, "MIP", "_merged", "MAIN");

&DefineParameters("pPicardToolsMarkduplicates", "program", 1, "MIP", "_pmd", "MAIN");

&DefineParametersPath("PicardToolsTempDirectory", "/scratch/", "pBwaMem,pPicardToolsSortSam,pPicardToolsMergeSamFiles,pPicardToolsMarkduplicates", 0); #Directory created by sbatch script and '$SLURM_JOB_ID' is appended to TMP directory

my (@picardToolsMergeSamFilesPrevious); #Any previous sequencing runs

##Coverage
&DefineParameters("pChanjoBuild", "program", 1, "MIP", "nofileEnding", "CoverageReport");

&DefineParametersPath("chanjoBuildDb", "nodefault", "pChanjoBuild", "file", "noAutoBuild");

&DefineParameters("pChanjoCalculate", "program", 1, "MIP","_coverage", "CoverageReport");

&DefineParameters("chanjoCalculateCutoff", "program", 10, "pChanjoCalculate");

&DefineParameters("pChanjoImport", "program", 1, "MIP", "nofileEnding", "CoverageReport");

&DefineParameters("pGenomeCoverageBED", "program", 1, "MIP", "_genomeCoverageBed", "CoverageQC_GcovBed", "bedtools");

&DefineParameters("pQaCompute", "program", 1, "MIP", "_qaCompute", "CoverageQC_QAComp", "qaCompute");

&DefineParameters("pPicardToolsCollectMultipleMetrics", "program", 1, "MIP", "nofileEnding", "CoverageQC_PTCMM");

&DefineParameters("pPicardToolsCalculateHSMetrics", "program", 1, "MIP", "nofileEnding", "CoverageQC_PTCHSM");

&DefineParameters("xCoverage", "program", 30, "pGenomeCoverageBED,pQaCompute");

&DefineParameters("pRCovPlots", "program", 0, "MIP", "nofileEnding", "CoverageQC_RCOVP");

&DefineParametersPath("picardToolsPath", "nodefault", "pBwaMem,pPicardToolsMergeSamFiles,pPicardToolsMarkduplicates,pPicardToolsCalculateHSMetrics,pPicardToolsCollectMultipleMetrics,pGATKHaploTypeCaller,pGATKVariantRecalibration", "directory"); #pGATKHaploTypeCaller,pGATKVariantRecalibration since these jars can use merged interval_list files, which are created in MIP with picardTools

##Target definition files
my (@exomeTargetBedInfileLists, @exomeTargetPaddedBedInfileLists); #Arrays for target bed infile lists

##GATK

&DefineParameters("pGATKRealigner", "program", 1, "MIP", "_rreal", "MAIN");

&DefineParametersPath("GATKReAlignerINDELKnownSet1", "1000G_phase1.indels.b37.vcf", "pGATKRealigner", "file", "noAutoBuild");

&DefineParametersPath("GATKReAlignerINDELKnownSet2", "Mills_and_1000G_gold_standard.indels.b37.sites.vcf", "pGATKRealigner", "file", "noAutoBuild");


&DefineParameters("pGATKBaseRecalibration", "program", 1, "MIP", "_brecal", "MAIN");

&DefineParametersPath("GATKBaseReCalibrationSNPKnownSet", "dbsnp_138.b37.vcf", "pGATKBaseRecalibration", "file", "noAutoBuild");


&DefineParameters("pGATKReduceReads", "program", 1, "MIP", "_reduced", "MAIN");


&DefineParameters("pGATKHaploTypeCaller", "program", 1, "MIP", "_", "MAIN");

&DefineParametersPath("GATKHaploTypeCallerSNPKnownSet", "dbsnp_138.b37.vcf", "pGATKHaploTypeCaller", "file", "noAutoBuild");

&DefineParametersPath("GATKHaploTypeCallerRefBAMInfile", "nodefault", "pGATKHaploTypeCaller", "file", "noAutoBuild");


&DefineParameters("pGATKHaploTypeCallerCombineVariants", "program", 1, "MIP", "nofileEnding", "MAIN");

&DefineParameters("pGATKVariantRecalibration", "program", 1, "MIP", "vrecal_", "MAIN");

&DefineParametersPath("GATKExomeReferenceSNPs", "nodefault", "pGATKVariantRecalibration", "file", "noAutoBuild");

&DefineParametersPath("GATKVariantReCalibrationTrainingSetHapMap", "hapmap_3.3.b37.sites.vcf", "pGATKVariantRecalibration", "file", "noAutoBuild");

&DefineParametersPath("GATKVariantReCalibrationTrainingSetDbSNP", "dbsnp_138.b37.vcf", "pGATKVariantRecalibration", "file", "noAutoBuild");

&DefineParametersPath("GATKVariantReCalibrationTrainingSet1000GSNP", "1000G_phase1.snps.high_confidence.b37.vcf", "pGATKVariantRecalibration", "file", "noAutoBuild");

&DefineParametersPath("GATKVariantReCalibrationTrainingSet1000GOmni", "1000G_omni2.5.b37.sites.vcf", "pGATKVariantRecalibration", "file", "noAutoBuild");

&DefineParametersPath("GATKVariantReCalibrationTrainingSetMills", "Mills_and_1000G_gold_standard.indels.b37.sites.vcf", "pGATKVariantRecalibration", "file", "noAutoBuild");

&DefineParameters("GATKVariantReCalibrationTSFilterLevel", "program", 99.9, "pGATKVariantRecalibration");

 
&DefineParameters("pGATKPhaseByTransmission", "program", 1, "MIP", "phtr_", "Phasing");

&DefineParameters("pGATKReadBackedPhasing", "program", 1, "MIP", "phrb_", "Phasing");

&DefineParameters("GATKReadBackedPhasingPhaseQualityThresh", "program", 20, "pGATKReadBackedPhasing");


&DefineParameters("pGATKVariantEvalAll", "program", 1, "MIP", "nofileEnding", "AllVariantQC");

&DefineParameters("pGATKVariantEvalExome", "program", 1, "MIP", "nofileEnding", "ExomeVarintQC", "bedtools");

&DefineParametersPath("GATKVariantEvalDbSNP", "dbsnp_138.b37.excluding_sites_after_129.vcf", "pGATKVariantEvalAll,pGATKVariantEvalExome", "file", "noAutoBuild");

&DefineParametersPath("GATKVariantEvalGold", "Mills_and_1000G_gold_standard.indels.b37.sites.vcf", "pGATKVariantEvalAll,pGATKVariantEvalExome", "file", "noAutoBuild");

&DefineParametersPath("genomeAnalysisToolKitPath", "nodefault", "pGATKRealigner,pGATKBaseRecalibration,pGATKHaploTypeCaller,pGATKVariantRecalibration,pGATKPhaseByTransmission,pGATKReadBackedPhasing,pGATKVariantEvalAll,pGATKVariantEvalExome", "directory");

&DefineParametersPath("GATKTempDirectory", "/scratch/", "pGATKRealigner,pGATKBaseRecalibration,pGATKReduceReads,pGATKHaploTypeCaller,pGATKVariantRecalibration,pGATKReadBackedPhasing", 0); #Depends on -projectID input, directory created by sbatch script and '$SLURM_JOB_ID' is appended to TMP directory

&DefineParameters("GATKDownSampleToCoverage", "program", 1000, "pGATKRealigner,pGATKBaseRecalibration,pGATKHaploTypeCaller");

my (@GATKTargetPaddedBedIntervalLists); #Array for target infile lists used in GATK

&DefineParametersPath("javaUseLargePages", "no", "pGATKRealigner,pGATKBaseRecalibration,pGATKReduceReads,pGATKHaploTypeCaller");

##Annovar

&DefineParameters("pAnnovar", "program", 1, "MIP", "annovar_", "MAIN");

&DefineParametersPath("annovarPath", "nodefault", "pAnnovar", "directory"); #Note not projectID specific

&DefineParameters("annovarGenomeBuildVersion", "program", "hg19", "pAnnovar");

&DefineParameters("annovarSupportedTableNames", "program", 0, "pAnnovar");

&DefineParameters("annovarMAFThreshold", "program", 0, "pAnnovar");

my @annovarTableNames; #List of Annovar table names to be used


##VMerge

&DefineParameters("pMergeAnnotatedVariants", "program", 1, "MIP", "merged_", "MAIN");

&DefineParametersPath("mergeAnnotatedVariantsTemplateFile", "nodefault", "pMergeAnnotatedVariants", "file", "noAutoBuild");

&DefineParameters("mergeAnnotatedVariantsDbFile", "program", "notSetYet", "pMergeAnnotatedVariants"); #No file check since file is created by MIP later


##Add_depth

&DefineParameters("pAddDepth", "program", 1, "MIP", "", "MAIN");


##RankVariants

&DefineParameters("pRankVariants", "program", 1, "MIP", "nofileEnding", "MAIN");

&DefineParameters("rankScore", "program", -100, "pRankVariants");

&DefineParametersPath("ImportantDbFile", "nodefault", "pRankVariants", "file", "noAutoBuild");

&DefineParametersPath("ImportantDbTemplate", "nodefault", "pRankVariants", "file", "noAutoBuild");

&DefineParameters("ImportantDbMasterFile", "program", "notSetYet", "pRankVariants"); #No file check since file is created by MIP later

my @ImportantDbFileOutFiles; #List of db outfiles

&DefineParametersPath("pythonVirtualEnvironment", "nodefault", "pChanjoBuild,pChanjoCalculate,pChanjoImport,pRankVariants");

##SChecks
&DefineParameters("pSampleCheck", "program", 1, "MIP", "nofileEnding", "IDQC", "vcftools:plink");

##QcCollect

&DefineParameters("pQCCollect", "program", 1, "MIP", "nofileEnding", "QCMetrics");

&DefineParameters("QCCollectSampleInfoFile", "program", "notSetYet", "pQCCollect"); #No file check since file is created by MIP later

&DefineParametersPath("QCCollectRegExpFile", "qc_regexp.yaml", "pQCCollect", "file", "noAutoBuild");


##MIP

##humanGenomeReference
&DefineParametersPath("humanGenomeReference", "nodefault", "pBwaMem,pBwaAln,pBwaSampe,pGATKRealigner,pGATKBaseRecalibration,pGATKHaploTypeCaller,pGATKHaploTypeCallerCombineVariants,pGATKVariantRecalibration,pGATKVariantEvalAll,pGATKVariantEvalExome,pAnnovar,pAddDepth,pPicardToolsCalculateHSMetrics,pPicardToolsCollectMultipleMetrics", "file", "noAutoBuild");
my @humanGenomeReferenceFileEndings = (".dict", ".fasta.fai"); #Meta files

my ($humanGenomeReferenceSource, $humanGenomeReferenceVersion, $humanGenomeReferenceNameNoEnding, $humanGenomeCompressed, $fnend, $aligner, $fileName, $fileNameTracker, $version, $help) = ("nocmdinput", "nocmdinput", "nocmdinput", "nocmdinput", ".sh", "nocmdinput", "nocmdinput", 0);

my (@contigs);

my (%infile, %indirpath, %infilesLaneNoEnding, %lane, %infilesBothStrandsNoEnding, %jobID, %sampleInfo); 

####Staging/Sanity Check Area 

##Capture kits supported from pedigree file.
my %supportedCaptureKits = (
    'Nimblegen_SeqCapEZExome.V2' => "Nimblegen_SeqCapEZExome.V2.GenomeReferenceSourceVersion_targets.bed",
    'Agilent_SureSelect.V2' => "Agilent_SureSelect.V2.GenomeReferenceSourceVersion_targets.bed",
    'Agilent_SureSelect.V3' => "Agilent_SureSelect.V3.GenomeReferenceSourceVersion_targets.bed",
    'Agilent_SureSelect.V4' => "Agilent_SureSelect.V4.GenomeReferenceSourceVersion_targets.bed",
    'Agilent_SureSelect.V5' => "Agilent_SureSelect.V5.GenomeReferenceSourceVersion_targets.bed",
    'Latest' => "Agilent_SureSelect.V5.GenomeReferenceSourceVersion_targets.bed",
    );

my %referenceFileEndings = (
    'mosaikAlignReference' => ".dat",
    'mosaikJumpDbStub' => "_jdb_15",
    'bwaBuildReference' => "",
    'exomeTargetBedInfileLists' => ".infile_list",
    'exomeTargetPaddedBedInfileLists' => ".pad100.infile_list",
    'GATKTargetPaddedBedIntervalLists' => ".pad100.interval_list",
    );

##Set supported annovar table name filtering options
my @annovarSupportedTableNames = ("refGene", "knownGene", "ensGene", "mce46way", "gerp++elem", "segdup", "gwascatalog", "tfbs", "mirna", "snp137", "snp135", "snp132", "snp131", "snp130", "snp129", "snp137NonFlagged", "snp135NonFlagged", "snp132NonFlagged", "snp131NonFlagged", "snp130NonFlagged", "1000g2012apr_all", "1000g2012apr_amr", "1000g2012apr_eur", "1000g2012apr_asn", "1000g2012apr_afr", "1000g2012feb_all", "esp6500si_all", "esp6500_all", "esp6500_aa", "esp6500_ea", "esp5400_all", "esp5400_aa", "esp5400_ea","clinvar_20131105", "ljb2_sift", "ljb2_pp2hdiv", "ljb2_pp2hvar", "ljb2_mt", "ljb2_ma", "ljb2_fathmm", "ljb2_siphy", "ljb2_lrt", "ljb_all", "ljb2_gerp++", "ljb2_phylop"); #Used to print list of supported table names

my %annovarTables;

###User Options
#$parameter{''}{'value'}
GetOptions('ifd|inFilesDirs:s'  => \@inFilesDirs, #Comma separated list
	   'isd|inScriptDir:s'  => \$parameter{'inScriptDir'}{'value'}, #Directory for custom scripts required by the pipeline
	   'rd|referencesDir:s'  => \$parameter{'referencesDir'}{'value'}, #directory containing references
	   'p|projectID:s'  => \$parameter{'projectID'}{'value'},
	   's|sampleIDs:s'  => \@sampleIDs, #Comma separated list, one below outDataDir
	   'em|email:s'  => \$parameter{'email'}{'value'},
	   'odd|outDataDir:s'  => \$parameter{'outDataDir'}{'value'}, #One dir above sample id, must supply whole path i.e. /proj/...
	   'osd|outScriptDir:s'  => \$parameter{'outScriptDir'}{'value'},  #One dir above sample id, must supply whole path i.e. /proj/...
	   'f|familyID:s' => \$parameter{'familyID'}{'value'}, #Family group ID (Merged to same vcf file after GATK Base Recalibration)
	   'pedigree|pedigreeFile:s' => \$parameter{'pedigreeFile'}{'value'}, #Pedigree file
	   'huref|humanGenomeReference:s' => \$parameter{'humanGenomeReference'}{'value'}, #Human genome reference
	   'al|aligner:s' => \$parameter{'aligner'}{'value'}, #determining which aligner was used previously (if not specified)
	   'at|analysisType:s' => \$parameter{'analysisType'}{'value'}, #Type of analysis
	   'mc|maximumCores:n' => \$parameter{'maximumCores'}{'value'}, #Per node
	   'c|configFile:s' => \$parameter{'configFile'}{'value'},
	   'wc|writeConfigFile:s' => \$parameter{'writeConfigFile'}{'value'},
	   'si|sampleInfoFile:s' => \$parameter{'sampleInfoFile'}{'value'}, #Write all info on samples and run to YAML file
	   'dra|dryRunAll:n' => \$parameter{'dryRunAll'}{'value'},
	   'pve|pythonVirtualEnvironment:s' => \$parameter{'pythonVirtualEnvironment'}{'value'},
	   'julp|javaUseLargePages:s' => \$parameter{'javaUseLargePages'}{'value'},
	   'h|help' => \$help, #Display help text
	   'v|version' => \$version, #Display version number
	   'pGZ|pGZip:n' => \$parameter{'pGZip'}{'value'},
	   'pFQC|pFastQC:n' => \$parameter{'pFastQC'}{'value'},
	   'pREM|pRemovalRedundantFiles:n' => \$parameter{'pRemovalRedundantFiles'}{'value'},
	   'pMoB|pMosaikBuild:n' => \$parameter{'pMosaikBuild'}{'value'},
	   'mobmfl|mosaikBuildMedianFragLength:n' => \$parameter{'mosaikBuildMedianFragLength'}{'value'}, #for fragment length estimation and local search
	   'pMoA|pMosaikAlign:n' => \$parameter{'pMosaikAlign'}{'value'},
	   'moaref|mosaikAlignReference:s' => \$parameter{'mosaikAlignReference'}{'value'}, #MosaikAlign reference file assumes existance of jump database files in same dir
	   'moaannpe|mosaikAlignNeuralNetworkPeFile:s' => \$parameter{'mosaikAlignNeuralNetworkPeFile'}{'value'},
	   'moaannse|mosaikAlignNeuralNetworkSeFile:s' => \$parameter{'mosaikAlignNeuralNetworkSeFile'}{'value'}, 
	   'mojdb|mosaikJumpDbStub:s' => \$parameter{'mosaikJumpDbStub'}{'value'}, #Stub for MosaikJump database
	   'pBWA_mem|pBwaMem:n' => \$parameter{'pBwaMem'}{'value'},
	   'bwamemrdb|bwaMemRapidDb:s' => \$parameter{'bwaMemRapidDb'}{'value'},
	   'pBWA_aln|pBwaAln:n' => \$parameter{'pBwaAln'}{'value'},
	   'bwaalnq|bwaAlnQualityTrimming:n' => \$parameter{'bwaAlnQualityTrimming'}{'value'}, #BWA aln quality threshold for read trimming down to 35bp
	   'pBWA_sampe|pBwaSampe:n' => \$parameter{'pBwaSampe'}{'value'},
	   'pPicT_sort|pPicardToolsSortSam:n' => \$parameter{'pPicardToolsSortSam'}{'value'},
	   'pPicT_merge|pPicardToolsMergeSamFiles:n' => \$parameter{'pPicardToolsMergeSamFiles'}{'value'}, #PicardTools MergeSamFiles
	   'pPicT_mergerr|pPicardToolsMergeRapidReads:n' => \$parameter{'pPicardToolsMergeRapidReads'}{'value'}, #PicardTools MergeSamFiles - Rapid mode
	   'pictmergetmpd|PicardToolsTempDirectory:s' => \$parameter{'PicardToolsTempDirectory'}{'value'}, #PicardToolsMerge Temporary Directory
	   'pict_mergeprev|picardToolsMergeSamFilesPrevious:s' => \@picardToolsMergeSamFilesPrevious, #Comma separated list
	   'pPicT_markdup|pPicardToolsMarkduplicates:s' => \$parameter{'pPicardToolsMarkduplicates'}{'value'}, #PicardTools MarkDuplicates
	   'picardpath|picardToolsPath:s' => \$parameter{'picardToolsPath'}{'value'}, #Path to picardtools
	   'pCh_B|pChanjoBuild:n' => \$parameter{'pChanjoBuild'}{'value'},  #Build central SQLiteDatabase
	   'chbdb|chanjoBuildDb:s' => \$parameter{'chanjoBuildDb'}{'value'}, #Chanjo reference database
	   'pCh_C|pChanjoCalculate:n' => \$parameter{'pChanjoCalculate'}{'value'},  # Chanjo coverage analysis
	   'chccut|chanjoCalculateCutoff:n' => \$parameter{'chanjoCalculateCutoff'}{'value'},  # Cutoff used for completeness
	   'pCh_I|pChanjoImport:n' => \$parameter{'pChanjoImport'}{'value'},  #Build family SQLiteDatabase
	   'pCC_bedgc|pGenomeCoverageBED:n' => \$parameter{'pGenomeCoverageBED'}{'value'},
	   'pCC_qac|pQaCompute:n' => \$parameter{'pQaCompute'}{'value'},
	   'xcov|xCoverage:n' => \$parameter{'xCoverage'}{'value'}, #Sets max depth to calculate coverage
	   'pCC_picmm|pPicardToolsCollectMultipleMetrics:n' => \$parameter{'pPicardToolsCollectMultipleMetrics'}{'value'},
	   'pCCE_pichs|pPicardToolsCalculateHSMetrics:n' => \$parameter{'pPicardToolsCalculateHSMetrics'}{'value'},
	   'extbl|exomeTargetBedInfileLists:s' => \@exomeTargetBedInfileLists, #Comma separated list of target file for CalculateHsMetrics
	   'extpbl|exomeTargetPaddedBedInfileLists:s' => \@exomeTargetPaddedBedInfileLists, #Comma separated list of padded target file for CalculateHsMetrics
	   'pRCP|pRCovPlots:n' => \$parameter{'pRCovPlots'}{'value'},
	   'gatkpath|genomeAnalysisToolKitPath:s' => \$parameter{'genomeAnalysisToolKitPath'}{'value'}, #GATK whole path
	   'gatktmpd|GATKTempDirectory:s' => \$parameter{'GATKTempDirectory'}{'value'}, #GATK ReAlignerTargetCreator & BaseRecalibrator temporary directory
	   'gatktpbl|GATKTargetPaddedBedIntervalLists:s' => \@GATKTargetPaddedBedIntervalLists, #Comma separated list of padded target file set to be used in GATK
	   'gatkdcov|GATKDownSampleToCoverage:n' => \$parameter{'GATKDownSampleToCoverage'}{'value'}, #GATK downsample to coverage
	   'pGATK_real|pGATKRealigner:n' => \$parameter{'pGATKRealigner'}{'value'}, #GATK ReAlignerTargetCreator/IndelRealigner
	   'gatkrealknset1|GATKReAlignerINDELKnownSet1:s' => \$parameter{'GATKReAlignerINDELKnownSet1'}{'value'}, #Known INDEL set to be used in GATK ReAlignerTargetCreator/IndelRealigner
	   'gatkrealknset2|GATKReAlignerINDELKnownSet2:s' => \$parameter{'GATKReAlignerINDELKnownSet2'}{'value'}, #Known INDEL set to be used in GATK ReAlignerTargetCreator/IndelRealigner
	   'pGATK_baserecal|pGATKBaseRecalibration:n' => \$parameter{'pGATKBaseRecalibration'}{'value'}, #GATK BaseRecalibrator/PrintReads
	   'gatkbaserecalknset|GATKBaseReCalibrationSNPKnownSet:s' => \$parameter{'GATKBaseReCalibrationSNPKnownSet'}{'value'}, #Known SNP set to be used in GATK BaseRecalibrator/PrintReads
	   'pGATK_rr|pGATKReduceReads:n' => \$parameter{'pGATKReduceReads'}{'value'}, #GATK ReduceReads
	   'pGATK_hapcall|pGATKHaploTypeCaller:n' => \$parameter{'pGATKHaploTypeCaller'}{'value'}, #GATK Haplotypecaller
	   'gatkhapcallrefbaminfile|GATKHaploTypeCallerRefBAMInfile:s' => \$parameter{'GATKHaploTypeCallerRefBAMInfile'}{'value'}, #GATK Haplotypecaller BAM reference infiles
	   'gatkhapcallsnpknset|GATKHaploTypeCallerSNPKnownSet:s' => \$parameter{'GATKHaploTypeCallerSNPKnownSet'}{'value'}, #Known SNP set to be used in GATK HaplotypeCaller
	   'pGATK_hapcallcombine|pGATKHaploTypeCallerCombineVariants:n' => \$parameter{'pGATKHaploTypeCallerCombineVariants'}{'value'}, #Combine variants from Haplotypecaller
	   'pGATK_varrecal|pGATKVariantRecalibration:n' => \$parameter{'pGATKVariantRecalibration'}{'value'}, #GATK VariantRecalibrator/ApplyRecalibration
	   'gatkexrefsnp|GATKExomeReferenceSNPs:s' => \$parameter{'GATKExomeReferenceSNPs'}{'value'}, #File of 33 exomes to power probabalistic model GATK Varrecal (SNVs) (Recieved from Måns, 120413)
	   'gatkvarrecaltrhapmap|GATKVariantReCalibrationTrainingSetHapMap:s' => \$parameter{'GATKVariantReCalibrationTrainingSetHapMap'}{'value'}, #GATK VariantRecalibrator resource
	   'gatkvarrecaltrdbsnp|GATKVariantReCalibrationTrainingSetDbSNP:s' => \$parameter{'GATKVariantReCalibrationTrainingSetDbSNP'}{'value'}, #GATK VariantRecalibrator resource
	   'gatkvarrecaltrd1000Gsnp|GATKVariantReCalibrationTrainingSet1000GSNP:s' => \$parameter{'GATKVariantReCalibrationTrainingSet1000GSNP'}{'value'}, #GATK VariantRecalibrator resource
	   'gatkvarrecaltromni|GATKVariantReCalibrationTrainingSet1000GOmni:s' => \$parameter{'GATKVariantReCalibrationTrainingSet1000GOmni'}{'value'}, #GATK VariantRecalibrator resource
	   'gatkvarrecaltrdbmills|GATKVariantReCalibrationTrainingSetMills:s' => \$parameter{'GATKVariantReCalibrationTrainingSetMills'}{'value'}, #GATK VariantRecalibrator resource
	   'gatkvarrecaltsfilterlevel|GATKVariantReCalibrationTSFilterLevel:s' => \$parameter{'GATKVariantReCalibrationTSFilterLevel'}{'value'}, #Truth sensativity level
	   'pGATK_phaseTr|pGATKPhaseByTransmission:n' => \$parameter{'pGATKPhaseByTransmission'}{'value'}, #GATK PhaseByTransmission to produce phased genotype calls
	   'pGATK_readPh|pGATKReadBackedPhasing:n' => \$parameter{'pGATKReadBackedPhasing'}{'value'}, #GATK ReadBackedPhasing
	   'gatkreadphphaseqthr|GATKReadBackedPhasingPhaseQualityThresh' => \$parameter{'GATKReadBackedPhasingPhaseQualityThresh'}{'value'}, #quality score required to output phasing
	   'pGATK_varevalall|pGATKVariantEvalAll:n' => \$parameter{'pGATKVariantEvalAll'}{'value'}, #GATK varianteval all variants
	   'pGATK_varevalexome|pGATKVariantEvalExome:n' => \$parameter{'pGATKVariantEvalExome'}{'value'}, #GATK varianteval only exonic variants
	   'gatkvarevaldbsnp|GATKVariantEvalDbSNP:s' => \$parameter{'GATKVariantEvalDbSNP'}{'value'},
	   'gatkvarevaldbgold|GATKVariantEvalGold:s' => \$parameter{'GATKVariantReCalibrationTrainingSetMills'}{'value'},
	   'pANVAR|pAnnovar:n' => \$parameter{'pAnnovar'}{'value'}, #Performs annovar filter gene, region and filter analysis
	   'anvarpath|annovarPath:s'  => \$parameter{'annovarPath'}{'value'}, #path to annovar script dir
	   'anvargbv|annovarGenomeBuildVersion:s'  => \$parameter{'annovarGenomeBuildVersion'}{'value'},
	   'anvartn|annovarTableNames:s'  => \@annovarTableNames, #Comma sepatated list
	   'anvarstn|annovarSupportedTableNames:n' => \$parameter{'annovarSupportedTableNames'}{'value'}, #Generates a list of supported table names
	   'anvarmafth|annovarMAFThreshold:n' => \$parameter{'annovarMAFThreshold'}{'value'},
	   'pMerge_anvar|pMergeAnnotatedVariants:n' => \$parameter{'pMergeAnnotatedVariants'}{'value'}, #Merges annovar analysis results to one master file
	   'mergeanvarte|mergeAnnotatedVariantsTemplateFile:s' => \$parameter{'mergeAnnotatedVariantsTemplateFile'}{'value'}, #Template file to create the specific family db master file
	   'mergeanvardbf|mergeAnnotatedVariantsDbFile:s' => \$parameter{'mergeAnnotatedVariantsDbFile'}{'value'}, #db master file to use when collecting external data
	   'pAdd_dp|pAddDepth:n' => \$parameter{'pAddDepth'}{'value'}, #Adds depth (DP) for nonvariants to master file (annovar_merged.txt)
	   'pRankVar|pRankVariants:n' => \$parameter{'pRankVariants'}{'value'}, #Ranking variants
	   'rs|rankscore:n'  => \$parameter{'rankScore'}{'value'}, #The rank score cut-off
	   'imdbfile|ImportantDbFile:s'  => \$parameter{'ImportantDbFile'}{'value'}, #Db of important genes
	   'imdbte|ImportantDbTemplate:s' => \$parameter{'ImportantDbTemplate'}{'value'}, #Template file to create the specific family selectVariants db master file
	   'imdbmf|ImportantDbMasterFile:s' => \$parameter{'ImportantDbMasterFile'}{'value'}, #Specific db master file to use when collecting external dataselectingVariants 
	   'imdbfof|ImportantDbFileOutFiles:s' => \@ImportantDbFileOutFiles, #The intersectCollect select variants output directorys	      
	   'pSCheck|pSampleCheck:n' => \$parameter{'pSampleCheck'}{'value'}, #QC for samples gender and relationship
	   'pQCC|pQCCollect:n' => \$parameter{'pQCCollect'}{'value'}, #QCmetrics collect
	   'QCCsampleinfo|QCCollectSampleInfoFile:s' => \$parameter{'QCCollectSampleInfoFile'}{'value'}, #SampleInfo yaml file produced by MIP
	   'QCCregexp|QCCollectRegExpFile:s' => \$parameter{'QCCollectRegExpFile'}{'value'}, #Regular expression yaml file
    );

if($help) {

    print STDOUT $USAGE, "\n";
    exit;
}

if($version) {

    print STDOUT "\nMip.pl v1.5.1\n\n";
    exit;
}

if ($parameter{'configFile'}{'value'} ne "nocmdinput") { #No input from cmd

    %scriptParameter = &LoadYAML($parameter{'configFile'}{'value'}); #Load parameters from configfile

    &OverWriteConfigParamWithCMDInfo("analysisType");
    &OverWriteConfigParamWithCMDInfo("aligner");

    foreach my $orderParameterElement (@orderParameters) { #Loop through all parameters and update info   

	&UpdateYAML($orderParameterElement, $scriptParameter{'clusterConstantPath'}, $scriptParameter{'analysisConstantPath'}, $scriptParameter{'analysisType'}, $parameter{'familyID'}{'value'}, $scriptParameter{'aligner'} );
    }
}

if ($parameter{'annovarSupportedTableNames'}{'value'} eq 1) {

    print STDOUT "\nThese Annovar databases are supported by MIP:\n";

    foreach my $annovarSupportedTableName (@annovarSupportedTableNames) {

	print STDOUT $annovarSupportedTableName, "\n";
    }
    print STDOUT "\n";
    exit;
}

foreach my $orderParameterElement (@orderParameters) { #Populate scriptParameters{'parameterName'} => 'Value'
    
##3 type of variables: MIP, path or program/program_parameters each is handled in the &AddToScriptParameter subroutine.
##parameterName, parameterValue, parameterType, parameterDefault, AssociatedProgram, Check directory/file existence)    
    &AddToScriptParameter($orderParameterElement, $parameter{$orderParameterElement}{'value'}, $parameter{$orderParameterElement}{'type'}, $parameter{$orderParameterElement}{'default'}, $parameter{$orderParameterElement}{'associatedProgram'}, $parameter{$orderParameterElement}{'existsCheck'}, $parameter{$orderParameterElement}{'programNamePath'});
   
    if ($orderParameterElement eq "outDataDir") { #Set defaults depending on $scriptParameter{'outDataDir'} value that now has been set

	$parameter{'sampleInfoFile'}{'default'} = $scriptParameter{'outDataDir'}."/".$scriptParameter{'familyID'}."/".$scriptParameter{'familyID'}."_qc_sampleInfo.yaml";
	$parameter{'QCCollectSampleInfoFile'}{'default'} = $parameter{'sampleInfoFile'}{'default'};

	$parameter{'mergeAnnotatedVariantsDbFile'}{'default'} = $scriptParameter{'outDataDir'}."/".$scriptParameter{'familyID'}."/".$scriptParameter{'familyID'}."_intersectCollect_db_master.txt";
	
	$parameter{'ImportantDbMasterFile'}{'default'} = $scriptParameter{'outDataDir'}."/".$scriptParameter{'familyID'}."/".$scriptParameter{'familyID'}.".intersectCollect_selectVariants_db_master.txt";
	
    }
    if ( $orderParameterElement eq "pedigreeFile") { #Write QC for only pedigree data used in analysis                                                        
	
	if (defined($scriptParameter{'pedigreeFile'})) {
	    `mkdir -p $scriptParameter{'outDataDir'}/$scriptParameter{'familyID'};`;
	    &WriteYAML($scriptParameter{'outDataDir'}."/".$scriptParameter{'familyID'}."/qc_pedigree.yaml", \%sampleInfo);
	}
    }
    if ( $orderParameterElement eq "humanGenomeReference") { #Supply humanGenomeReference to mosaikAlignReference if required
	
	if ( (defined($scriptParameter{'humanGenomeReference'})) && (defined($humanGenomeReferenceNameNoEnding)) ) {

	    &SetAutoBuildFeature("mosaikAlignReference", \$referenceFileEndings{'mosaikAlignReference'}, \$humanGenomeReferenceNameNoEnding);
	    &SetAutoBuildFeature("mosaikJumpDbStub", \$referenceFileEndings{'mosaikJumpDbStub'}, \$humanGenomeReferenceNameNoEnding);
	    &SetAutoBuildFeature("bwaBuildReference", \$referenceFileEndings{'bwaBuildReference'}, \$humanGenomeReferenceNameNoEnding);	
	}
    }
}


##sampleIDs
&PrepareArrayParameters(\@sampleIDs, "sampleIDs", "path", "nodefault", "MIP", "");
&CheckUniqueIDNs(); #Test that smapleIDs are unique
for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) { #all sampleIDs
    
    &ScriptParameterPerSampleID(\$scriptParameter{'familyID'}, \$sampleIDs[$sampleIDCounter], "exomeTargetBedInfileLists");
    &ScriptParameterPerSampleID(\$scriptParameter{'familyID'}, \$sampleIDs[$sampleIDCounter], "exomeTargetPaddedBedInfileLists");
    &ScriptParameterPerSampleID(\$scriptParameter{'familyID'}, \$sampleIDs[$sampleIDCounter], "GATKTargetPaddedBedIntervalLists");
}

##inFileDirs
&PrepareArrayParameters(\@inFilesDirs, "inFilesDirs", "path", "notSetYet", "MIP", "directory");


##picardToolsMergeSamFilesPrevious
if ( ($scriptParameter{'pPicardToolsMergeSamFiles'} > 0) || (scalar(@picardToolsMergeSamFilesPrevious) > 0)) { #2nd term to enable write to config

    &PrepareArrayParameters(\@picardToolsMergeSamFilesPrevious, "picardToolsMergeSamFilesPrevious", "path", "nodefault", "pPicardToolsMergeSamFiles", "file");    
     
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) { #Check all samples to check, which are to be merged with previous files later

	if (scalar(@picardToolsMergeSamFilesPrevious) > 0) { #Supplied info - check for which sampleIDs  	

	    for (my $mergeFileCounter=0;$mergeFileCounter<scalar(@picardToolsMergeSamFilesPrevious);$mergeFileCounter++) {
		
		if ($picardToolsMergeSamFilesPrevious[$mergeFileCounter] =~ /$sampleIDs[$sampleIDCounter]/) { #Look for sampleID in previously generated file to be merged with current run to be able to merge correct files

		    $sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{'picardToolsMergeSamFilesPrevious'} = 1;
		}
		else {

		    $sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{'picardToolsMergeSamFilesPrevious'} = 0;
		}
	    }
	}
	else { #Not supplied - Set to 0 

	    $sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{'picardToolsMergeSamFilesPrevious'} = 0;
	}
    }
}
else { #Not supplied - Set to 0 to handle correctly in program subroutines 

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) { #Set for all sampleIDs

	$sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{'picardToolsMergeSamFilesPrevious'} = 0;
    }
}

if ($scriptParameter{'pAnnovar'} > 0) {

    &PrepareArrayParameters(\@annovarTableNames, "annovarTableNames", "path", "yes", "pAnnovar", "file"); #"yes" added to enable addition of default table names in &AddToScriptParameters  
}

if ($scriptParameter{'pRankVariants'} > 0) {
    
    &UpdateYAML("ImportantDbFileOutFiles", $scriptParameter{'clusterConstantPath'}, $scriptParameter{'analysisConstantPath'}, $scriptParameter{'analysisType'},$parameter{'familyID'}{'value'}, $scriptParameter{'aligner'} );
    &PrepareArrayParameters(\@ImportantDbFileOutFiles, "ImportantDbFileOutFiles", "program", "yes", "pRankVariants");  
}

##Set Target files
&PrepareArrayParameters(\@exomeTargetBedInfileLists, "exomeTargetBedInfileLists", "path", "notSetYet", "pPicardToolsCalculateHSMetrics", "file");

&PrepareArrayParameters(\@exomeTargetPaddedBedInfileLists, "exomeTargetPaddedBedInfileLists", "path", "notSetYet", "pPicardToolsCalculateHSMetrics", "file");
 
&PrepareArrayParameters(\@GATKTargetPaddedBedIntervalLists, "GATKTargetPaddedBedIntervalLists", "path", "notSetYet", "pGATKHaploTypeCaller,pGATKVariantRecalibration", "file");

if ($scriptParameter{'writeConfigFile'} ne 0) { #Write config file for family

    &WriteYAML($scriptParameter{'writeConfigFile'}, \%scriptParameter); #Write used settings to configfile
}

##Set chr prefix and chromosome names depending on reference used
if ($scriptParameter{'humanGenomeReference'}=~/hg\d+/) { #Refseq - prefix and M
    @contigs = ("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM"); #Chr for filtering of bam file
}
elsif ($scriptParameter{'humanGenomeReference'}=~/GRCh\d+/) { #Ensembl - no prefix and MT
    @contigs = ("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT"); #Chr for filtering of bam file
}

####Creates master_log for the master script 

my ($base, $script) = (`date +%Y%m%d`,`basename $0`); #Catches current date and script name
chomp($base,$script); #Remove \n;
`mkdir -p $scriptParameter{'outDataDir'}/$scriptParameter{'familyID'}/mip_log/$base;`; #Creates the mip_log dir
my $mipLogName = $scriptParameter{'outDataDir'}."/".$scriptParameter{'familyID'}."/mip_log/".$base."/".$script."_".$timeStamp.".txt"; #concatenates mip_log filename

open (MIPLOG, ">>".$mipLogName) or die "Can't write to ".$mipLogName.":".$!, "\n"; #Open file masterLog

##Add parameters
print MIPLOG "\n".$script." "; #Adds script name to recontruct command line

&WriteCMDMipLog();

print STDOUT "\nScript parameters and info from ".$script." are saved in file: ".$mipLogName, "\n";

####Collect infiles

for (my $inputDirectoryCounter=0;$inputDirectoryCounter<scalar(@inFilesDirs);$inputDirectoryCounter++) { #Collects inputfiles
    
    my @infiles = `cd $inFilesDirs[ $inputDirectoryCounter ];ls *.fastq*;`; #cd to input dir and collect fastq files and fastq.gz files
   
    print STDOUT "\nReads from Platform", "\n";print MIPLOG "\nReads from Platform", "\n";
    print STDOUT "\nSample ID\t".$sampleIDs[$inputDirectoryCounter],"\n";print MIPLOG "\nSample ID\t".$sampleIDs[$inputDirectoryCounter],"\n";
    print STDOUT "Inputfiles\n",@ { $infile{ $sampleIDs[$inputDirectoryCounter] }  =[@infiles] }, "\n"; #hash with sample id as key and inputfiles in dir as array 
    print MIPLOG "Inputfiles\n",@ { $infile{ $sampleIDs[$inputDirectoryCounter] }  =[@infiles] }, "\n";
    
    $indirpath{$sampleIDs[$inputDirectoryCounter]} = $inFilesDirs[ $inputDirectoryCounter ];  #Catch inputdir path
    chomp(@infiles);    #Remove newline from every entry in array
    $infile{ $sampleIDs[$inputDirectoryCounter] }  =[@infiles]; #Reload files into hash (kept above newline just for print STDOUT)
}

close(MIPLOG);

my $uncompressedFileSwitch = &InfilesReFormat(); #Required to format infiles correctly for subsequent input into aligners
    
&CreateFileEndings(); #Creates all fileendings as the samples is processed depending on the chain of modules activated

#Create .fam file to be used in variant calling analyses
my $pqFamFile = q?perl -nae 'print $F[0], "\t", $F[1], "\t", $F[2], "\t", $F[3], "\t", $F[4], "\t", $F[5], "\n";'?;
my $famFile = $scriptParameter{'outDataDir'}."/".$scriptParameter{'familyID'}."/".$scriptParameter{'familyID'}.".fam";
`$pqFamFile $scriptParameter{'pedigreeFile'} > $famFile;`;

####MAIN

open (MIPLOG, ">>".$mipLogName) or die "Can't write to ".$mipLogName.":".$!, "\n"; #Open file run log

if ( ($scriptParameter{'pGZip'} > 0) && ($uncompressedFileSwitch eq "unCompressed") ) { #GZip of fastq files

    print STDOUT "\nGZip for fastq files", "\n";print MIPLOG "\nGZip for fastq files", "\n";

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {  

	for (my $infileCounter=0;$infileCounter<scalar( @{ $infile{$sampleIDs[$sampleIDCounter]} });$infileCounter++) { #To determine which sampleID had the uncompressed files
	    
	    if ($infile{$sampleIDs[$sampleIDCounter]}[$infileCounter] =~/.fastq$/) {
	
		&GZipFastq($sampleIDs[$sampleIDCounter]);
		last; #Return to sampleID loop i.e. only call subroutine GZipFastq once per sampleID
	    }
	}
    }
}

if ($scriptParameter{'pFastQC'} > 0) { #Run FastQC
    
    print STDOUT "\nFastQC", "\n";print MIPLOG "\nFastQC", "\n";
    
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {  
	
	&FastQC($sampleIDs[$sampleIDCounter]);	
    }
}

if ($scriptParameter{'pMosaikBuild'} > 0) { #Run MosaikBuild
    
    print STDOUT "\nMosaikBuild", "\n";print MIPLOG "\nMosaikBuild", "\n";

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {  
	
	&MosaikBuild($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});	
    }
}


if ($scriptParameter{'pMosaikAlign'} > 0) { #Run MosaikAlign

    print STDOUT "\nMosaikAlign", "\n"; print MIPLOG "\nMosaikAlign", "\n";

    if ( ($parameter{'mosaikAlignReference'}{'buildFile'} eq 1) || ($parameter{'mosaikJumpDbStub'}{'buildFile'} eq 1) ) {
		
	&BuildMosaikAlignPreRequisites($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "MosaikAlign");
	
    }
    if ( ($parameter{'mosaikAlignNeuralNetworkPeFile'}{'buildFile'} eq 1) || ($parameter{'mosaikAlignNeuralNetworkSeFile'}{'buildFile'} eq 1) ){

	&MoveMosaikNN();
    }

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {  
	
	&MosaikAlign($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});	
    }
}

if ($scriptParameter{'pBwaMem'} > 0) { #Run BWA Mem
    
    print STDOUT "\nBWA Mem", "\n";print MIPLOG "\nBWA Mem", "\n";

    &CheckBuildHumanGenomePreRequisites();
    
    if ($parameter{'bwaBuildReference'}{'buildFile'} eq 1) {
	
	&BuildBwaPreRequisites($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BwaMem");
    }

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {  
	
	&BWA_Mem($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});	
	
    }    
}

if ($scriptParameter{'pPicardToolsMergeRapidReads'} > 0) { #Run PicardToolsMergeRapidReads - Relevant only in rapid mode
    
    print STDOUT "\nPicardToolsMergeRapidReads", "\n";print MIPLOG "\nPicardToolsMergeRapidReads", "\n";
    
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {  
	
        #Merge all read batch processes to 1 file again containing sorted & indexed reads matching clinical test genes
	&PicardToolsMergeRapidReads($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});
    }    
}

if ($scriptParameter{'pBwaAln'} > 0) { #Run BWA Aln
    
    print STDOUT "\nBWA Aln", "\n";print MIPLOG "\nBWA Aln", "\n";

    &CheckBuildHumanGenomePreRequisites();

    if ($parameter{'bwaBuildReference'}{'buildFile'} eq 1) {
	
	&BuildBwaPreRequisites($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BwaAln");
    }
    
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {  
	
	&BWA_Aln($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});	
    }    
}

if ($scriptParameter{'pBwaSampe'} > 0) { #Run BWA Sampe
    
    print STDOUT "\nBWA Sampe", "\n";print MIPLOG "\nBWA Sampe", "\n";
    
    &CheckBuildHumanGenomePreRequisites();

    if ($parameter{'bwaBuildReference'}{'buildFile'} eq 1) {
	
	&BuildBwaPreRequisites($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BwaSampe");
    }

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {  
	
	&BWA_Sampe($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});
    }
}

if ($scriptParameter{'pPicardToolsSortSam'} > 0) { #Run Picardtools SortSam and Index

    if ($scriptParameter{'analysisType'} ne "rapid") { #In rapid mode Sort and index is done for each batch of reads in the BWA_Mem call, since the link to infile is broken by the read batch processing. However pPicardToolsSortSam should be enabled to ensure correct fileending and merge the flow to ordinary modules.

	print STDOUT "\nPicardTools SortSam & index", "\n";

	for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {  
	    
	    &PicardToolsSortSamIndex($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});
	}
    }
}

if ($scriptParameter{'pPicardToolsMergeSamFiles'} > 0) { #Run picardtools merge

    print STDOUT "\nPicardTool MergeSamFiles", "\n";

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {  

	if ( ($sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{'picardToolsMergeSamFilesPrevious'} == 1) || (scalar( @{ $infilesLaneNoEnding{ $sampleIDs[$sampleIDCounter] } }) > 1) ) { #Sanity Check that we have something to merge with
	
	    &PicardToolsMerge($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'}, $sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{'fileEnding'});	
	}
    }
}

if ($scriptParameter{'pPicardToolsMarkduplicates'} > 0) { #PicardTools MarkDuplicates

    print STDOUT "\nPicardTools MarkDuplicates", "\n";print MIPLOG "\nPicardTools MarkDuplicates", "\n";

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {  
    
	&PicardToolsMarkDuplicates($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});	
    }
}

if ($scriptParameter{'pChanjoBuild'} > 0) {
    
    print STDOUT "\nChanjoBuild", "\n";print MIPLOG "\nChanjoBuild", "\n";
    
    &ChanjoBuild($scriptParameter{'familyID'});
}

if ($scriptParameter{'pChanjoCalculate'} > 0) {
    
    print STDOUT "\nChanjoCalculate", "\n";print MIPLOG "\nChanjoCalculate", "\n";
    
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) { #For all SampleIDs
	
	&ChanjoCalculate($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});
    }
}

if ($scriptParameter{'pChanjoImport'} > 0) {
    
    print STDOUT "\nChanjoImport", "\n";print MIPLOG "\nChanjoImport", "\n";
    
    &ChanjoImport($scriptParameter{'familyID'}, $scriptParameter{'aligner'});
}

if ($scriptParameter{'pGenomeCoverageBED'} > 0) { #Run GenomeCoverageBED
    
    print STDOUT "\nGenomeCoverageBED", "\n";print MIPLOG "\nGenomeCoverageBED", "\n";    
    
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {  

	&GenomeCoverageBED($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});
    }
}

if ($scriptParameter{'pQaCompute'} > 0) { #Run QaCompute
    
    print STDOUT "\nQaCompute", "\n";print MIPLOG "\nQaCompute", "\n";    
    
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {  

	&QaCompute($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});
    }
}

if ($scriptParameter{'pPicardToolsCollectMultipleMetrics'} > 0) { #Run PicardToolsCollectMultipleMetrics
    
    print STDOUT "\nPicardToolsCollectMultipleMetrics", "\n";print MIPLOG "\nPicardToolsCollectMultipleMetrics", "\n";    
    
    &CheckBuildHumanGenomePreRequisites();

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {  

	&PicardToolsCollectMultipleMetrics($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});
    }
}

if ($scriptParameter{'pPicardToolsCalculateHSMetrics'} > 0) { #Run PicardToolsCalculateHSMetrics
    
    print STDOUT "\nPicardToolsCalculateHSMetrics", "\n";print MIPLOG "\nPicardToolsCalculateHSMetrics", "\n";    
    
    &CheckBuildHumanGenomePreRequisites();

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {

	if ($parameter{ $scriptParameter{'familyID'} }{$sampleIDs[$sampleIDCounter]}{'exomeTargetBedInfileLists'}{'buildFile'} eq 1) {
	    
	    &BuildPTCHSMetricPreRequisites($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "PicardToolsCalculateHSMetrics");
	    last; #Will handle all build per sampleID within sbatch script
	}
	if ($parameter{ $scriptParameter{'familyID'} }{$sampleIDs[$sampleIDCounter]}{'exomeTargetPaddedBedInfileLists'}{'buildFile'} eq 1) {
	    
	    &BuildPTCHSMetricPreRequisites($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "PicardToolsCalculateHSMetrics");
	    last; #Will handle all build per sampleID within sbatch script
	}
    }

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {  

	&PicardToolsCalculateHSMetrics($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});
    }
}

if ($scriptParameter{'pRCovPlots'} > 0) { #Run Rcovplot scripts   
    print STDOUT "\nRCovPlots", "\n";print MIPLOG "\nRCovPlots", "\n";	

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {  
	
	&RCoveragePlots($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});	
    }
}

if ($scriptParameter{'pGATKRealigner'} > 0) { #Run GATK ReAlignerTargetCreator/IndelRealigner

    print STDOUT "\nGATK ReAlignerTargetCreator/IndelRealigner", "\n";print MIPLOG "\nGATK ReAlignerTargetCreator/IndelRealigner", "\n";

    &CheckBuildHumanGenomePreRequisites();

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {   
    
	&GATKReAligner($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});	
    }
}

if ($scriptParameter{'pGATKBaseRecalibration'} > 0) { #Run GATK BaseRecalibrator/PrintReads

    print STDOUT "\nGATK BaseRecalibrator/PrintReads", "\n";print MIPLOG "\nGATK BaseRecalibrator/PrintReads", "\n";

    &CheckBuildHumanGenomePreRequisites();

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {   
  
	&GATKBaseReCalibration($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});	
    }
}

if ($scriptParameter{'pGATKReduceReads'} > 0) { #Run GATK ReduceReads

    print STDOUT "\nGATK ReduceReads", "\n";print MIPLOG "\nGATK ReduceReads", "\n";

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {   
    
	&GATKReduceReads($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});	
    }
}

if ($scriptParameter{'pGATKHaploTypeCaller'} > 0) { #Run GATK HaploTypeCaller. Done per family

    print STDOUT "\nGATK HaplotypeCaller", "\n";print MIPLOG "\nGATK HaplotypeCaller", "\n";

    &CheckBuildHumanGenomePreRequisites();

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {
	
	if ($parameter{ $scriptParameter{'familyID'} }{$sampleIDs[$sampleIDCounter]}{'GATKTargetPaddedBedIntervalLists'}{'buildFile'} eq 1) {
	    
	    &BuildPTCHSMetricPreRequisites($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "GATKHaploTypeCaller");
	    last; #Will handle all build per sampleID within sbatch script
	}
    }
    for (my $chromosomeCounter=0;$chromosomeCounter<scalar(@contigs);$chromosomeCounter++) {
	    
	&GATKHaploTypeCaller($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH",$contigs[$chromosomeCounter],8); #Argument 3 is which chr to analyse and Arg 4 is java heap allocation (Gb).
    }
}

if ($scriptParameter{'pGATKHaploTypeCallerCombineVariants'} > 0) { #Run GATK HaplotypeCallerCombineVariants. Done per family

    print STDOUT "\nGATK HaplotypeCallerCombineVariants", "\n";print MIPLOG "\nGATK HaplotypeCallerCombineVariants", "\n";

    &CheckBuildHumanGenomePreRequisites();
    &GATKHaplotypeCallerCombineVariants($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH");

}

if ($scriptParameter{'pGATKVariantRecalibration'} > 0) { #Run GATK VariantRecalibrator/ApplyRecalibration. Done per family

    print STDOUT "\nGATK VariantRecalibrator/ApplyRecalibration", "\n";print MIPLOG "\nGATK VariantRecalibrator/ApplyRecalibration", "\n";
    
    &CheckBuildHumanGenomePreRequisites();

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {
	
	if ($parameter{ $scriptParameter{'familyID'} }{$sampleIDs[$sampleIDCounter]}{'GATKTargetPaddedBedIntervalLists'}{'buildFile'} eq 1) {
	    
	    &BuildPTCHSMetricPreRequisites($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "GATKVariantRecalibration");
	    last; #Will handle all build per sampleID within sbatch script
	}
    }
    &GATKVariantReCalibration($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH");
}

if ($scriptParameter{'pAnnovar'} > 0) { #Run Annovar. Done per family

    print STDOUT "\nAnnovar", "\n";print MIPLOG "\nAnnovar", "\n";

    &CheckBuildHumanGenomePreRequisites();

    for (my $tableNamesCounter=0;$tableNamesCounter<scalar(@annovarTableNames);$tableNamesCounter++) { #For all specified table names

	if ($parameter{$annovarTableNames[$tableNamesCounter]}{'buildFile'} eq 1) {

	&BuildAnnovarPreRequisites($scriptParameter{'familyID'}, "mosaik", "Annovar");
	last; #Will handle all build tables within sbatch script
	}
    }
    &Annovar($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH");
}

if ($scriptParameter{'pGATKPhaseByTransmission'} > 0) { #Run GATK PhaseByTransmission. Done per family
    
    print STDOUT "\nGATK PhaseByTransmission", "\n";print MIPLOG "\nGATK PhaseByTransmission", "\n";
    
    &CheckBuildHumanGenomePreRequisites();
    &GATKPhaseByTransmission($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH");
}

if ($scriptParameter{'pGATKReadBackedPhasing'} > 0) { #Run GATK ReadBackedPhasing. Done per family. NOTE: Needs phased calls
    
    print STDOUT "\nGATK ReadBackedPhasing", "\n";print MIPLOG "\nGATK ReadBackedPhasing", "\n";
    
    &CheckBuildHumanGenomePreRequisites();
    &GATKReadBackedPhasing($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH");
}

if ($scriptParameter{'pGATKVariantEvalAll'} > 0) { #Run GATK VariantEval for all variants. Done per sampleID

    print STDOUT "\nGATK VariantEval All", "\n";print MIPLOG "\nGATK VariantEval All", "\n";

    &CheckBuildHumanGenomePreRequisites();

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) { 
	
	&GATKVariantEvalAll($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'}, "BOTH", $scriptParameter{'familyID'});
    }
}

if ($scriptParameter{'pMergeAnnotatedVariants'} > 0) { #Run MergeAnnotationVariants using intersectCollect.pl. Done per family

    print STDOUT "\nMergeAnnotatedVariants", "\n";print MIPLOG "\nMergeAnnotatedVariants", "\n";

    &MergeAnnotatedVariants($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH");

}

if ($scriptParameter{'pGATKVariantEvalExome'} > 0) { #Run GATK VariantEval for exome variants. Done per sampleID

    print STDOUT "\nGATK VariantEval Exome", "\n";print MIPLOG "\nGATK VariantEval Exome", "\n";
    
    &CheckBuildHumanGenomePreRequisites();

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) { 
	
	&GATKVariantEvalExome($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'}, "BOTH", $scriptParameter{'familyID'});
    }
}

if ($scriptParameter{'pAddDepth'} > 0) { #Run AddDepth using add_depth.pl. Done per family

    print STDOUT "\nAddDepth", "\n";print MIPLOG "\nAddDepth", "\n";

    &AddDp($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH");
}

if ($scriptParameter{'pRankVariants'} > 0) { #Run RankVariants. Done per family

    print STDOUT "\nRankVariants", "\n";print MIPLOG "\nRankVariants", "\n";

    &RankVariants($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH");
}

if ($scriptParameter{'pSampleCheck'} > 0) { #Run SampleCheck. Done per family

    print STDOUT "\nSampleCheck", "\n";print MIPLOG "\nSampleCheck", "\n";

    &SampleCheck($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH");
}

if ($scriptParameter{'pQCCollect'} > 0) { #Run QCCollect. Done per family

    print STDOUT "\nQCCollect", "\n";print MIPLOG "\nQCCollect", "\n";

    &QCCollect($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH");
}

if ($scriptParameter{'pRemovalRedundantFiles'} > 0) { #Sbatch generation of removal of alignment files
    
    print STDOUT "\nRemoval of alignment files", "\n"; print MIPLOG "\nRemoval of alignment files", "\n"; 
	
	&RemoveRedundantFiles($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH");	
}

close(MIPLOG); #Close mip_log file

#Write QC for programs used in analysis                                                                                                                                                                                           
if ($scriptParameter{'sampleInfoFile'} ne 0) {#Write SampleInfo to yaml file

    &WriteYAML($scriptParameter{'sampleInfoFile'}, \%sampleInfo); #Write QC for sampleinfo used in analysis
}


######################
###SubRoutines#######
######################

sub RemoveRedundantFiles {
#Generates a sbatch script, which removes some alignment files.
    
    my $familyID = $_[0];
    my $aligner = $_[1];
    my $callType = $_[2]; #SNV,INDEL or BOTH

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle

    &ProgramPreRequisites($familyID, "RemovalRedundantFiles", $aligner, 0, $FILEHANDLE, 1, 1);
    
    `mkdir -p $scriptParameter{'outDataDir'}/$familyID/$aligner/info;`; #Creates the aligner and info data file directory
    `mkdir -p $scriptParameter{'outScriptDir'}/$familyID/$aligner;`; #Creates the aligner script directory

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) { 

	my $sampleID = $sampleIDs[$sampleIDCounter];
	my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner;

###Single files

	for (my $infileCounter=0;$infileCounter < scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #MosaikBuild takes both reads at once
	    
	    my $infile = $infilesLaneNoEnding{$sampleID}[$infileCounter]; 
	    
##MosaikBuild	
	    if ( ($scriptParameter{'pMosaikBuild'} > 0) || ($scriptParameter{'aligner'} eq "mosaik") ) {
		
		print $FILEHANDLE "rm ";
		print $FILEHANDLE $inSampleDirectory."/".$infile.".dat", "\n\n"; #MosaikBuild
		
	    }
##MosaikAlign
	if ( ($scriptParameter{'pMosaikAlign'} > 0) || ($scriptParameter{'aligner'} eq "mosaik") ) {
	    
	    print $FILEHANDLE "rm ";
	    print $FILEHANDLE $inSampleDirectory."/".$infile.".stat", "\n\n"; #MosaikAlign Stats
	    
	    print $FILEHANDLE "rm ";
	    print $FILEHANDLE $inSampleDirectory."/".$infile.".bam", "\n\n"; #MosaikAlign	    
	}
	    
##Remove BWA files
	    if ($scriptParameter{'pBwaAln'} > 0) {
		
		print $FILEHANDLE "rm ";
		print $FILEHANDLE $inSampleDirectory."/".$infile.".sai", "\n\n"; #BWA_Aln
	    }
	    if ($scriptParameter{'pBwaSampe'} >0) {
		
		print $FILEHANDLE "rm ";
		print $FILEHANDLE $inSampleDirectory."/".$infile.".bam", "\n\n"; #BWA_Sampe
	    }      	    
	    
##Sorted BAM
	    if ($scriptParameter{'pPicardToolsSortSam'} > 0) {
		
		my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pPicardToolsSortSam'}{'fileEnding'};
		
		print $FILEHANDLE "rm ";
		print $FILEHANDLE $inSampleDirectory."/".$infile.$outfileEnding.".ba*", "\n\n"; #Sorted BAM and bai file
	    }
	}
	
##Potentially merged files
	my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles($sampleID);        
	
	if ($PicardToolsMergeSwitch == 1) { #Files was merged previously
	    
	    if ($scriptParameter{'pPicardToolsMergeSamFiles'} > 0) {
		
		my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pPicardToolsMergeSamFiles'}{'fileEnding'};
		
		print $FILEHANDLE "rm ";
		print $FILEHANDLE $inSampleDirectory."/".$infile.$outfileEnding.".ba*", "\n\n"; #Sorted BAM and bai file
	    }	
	    if ($scriptParameter{'pPicardToolsMarkduplicates'} > 0) {
		
		my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pPicardToolsMarkduplicates'}{'fileEnding'};
		
		print $FILEHANDLE "rm ";
		print $FILEHANDLE $inSampleDirectory."/".$infile.$outfileEnding.".ba*", "\n\n"; #Dedupped BAM and bai file
	    }
	    if ($scriptParameter{'pGATKRealigner'} > 0) {
		
	    my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/GATK";   
	    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pGATKRealigner'}{'fileEnding'};
	    
	    print $FILEHANDLE "rm ";
	    print $FILEHANDLE $inSampleDirectory."/".$infile.$outfileEnding.".ba*", "\n\n"; #ReAligned BAM and bai file
	    }
	    if ($scriptParameter{'pGATKBaseRecalibration'} > 0) {
		
		my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/GATK";
		my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pGATKBaseRecalibration'}{'fileEnding'};
		
		print $FILEHANDLE "rm ";
		print $FILEHANDLE $inSampleDirectory."/".$infile.$outfileEnding.".ba*", "\n\n"; #BaseRecalibrated BAM and bai file
	    }
	}
	else {
	    
	    for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #For all infiles per lane
		
		my $infile = $infilesLaneNoEnding{$sampleID}[$infileCounter];
		
		if ($scriptParameter{'pPicardToolsMarkduplicates'} > 0) {
		    
		    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pPicardToolsMarkduplicates'}{'fileEnding'};
		    
		    print $FILEHANDLE "rm ";
		    print $FILEHANDLE $inSampleDirectory."/".$infile.$outfileEnding.".ba*", "\n\n"; #Dedupped BAM and bai file
		}
		if ($scriptParameter{'pGATKRealigner'} > 0) {
		    
		    my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/GATK";
		    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pGATKRealigner'}{'fileEnding'};
		    
		    print $FILEHANDLE "rm ";
		    print $FILEHANDLE $inSampleDirectory."/".$infile.$outfileEnding.".ba*", "\n\n"; #ReAligned BAM and bai file
		}
		if ($scriptParameter{'pGATKBaseRecalibration'} > 0) {
		    
		    my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/GATK";
		    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pGATKBaseRecalibration'}{'fileEnding'};
		    
		    print $FILEHANDLE "rm ";
		    print $FILEHANDLE $inSampleDirectory."/".$infile.$outfileEnding.".ba*", "\n\n"; #BaseRecalibrated BAM and bai file
		}
	    }
	}
    }
###Family files
    if ($scriptParameter{'pGATKHaploTypeCaller'} > 0) {
	
	my $outFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK/HaploTypeCaller";
	
	print $FILEHANDLE "rm -rf ";
	print $FILEHANDLE $outFamilyDirectory."/", "\n\n"; #Remove HaplotypeCaller individual contigfiles
	
	$outFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK"; #New outfile directory
	my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{'pGATKHaploTypeCaller'}{'fileEnding'};
	
	print $FILEHANDLE "rm ";
	print $FILEHANDLE $outFamilyDirectory."/".$familyID.$outfileEnding.$callType.".vcf*", "\n\n"; #HaplotypeCaller vcf file
    }
    if ($scriptParameter{'pAnnovar'} > 0) {
	
	my $outFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
	my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{'pAnnovar'}{'fileEnding'};
	
	print $FILEHANDLE "rm ";
	print $FILEHANDLE $outFamilyDirectory."/".$familyID.$outfileEnding.$callType."*".$scriptParameter{'annovarGenomeBuildVersion'}."_*", "\n\n"; #Annovar data files
	print $FILEHANDLE "rm ";
	print $FILEHANDLE $outFamilyDirectory."/".$familyID.$outfileEnding.$callType.".*", "\n\n"; #Annovar data files
    }  
    close($FILEHANDLE);
    return;
}

sub SampleCheck { 
###Tests sample for correct relatives (only performed for samples with relatives defined in pedigree file) performed on sequence data.

    my $familyID = $_[0]; #familyID NOTE: not sampleid 
    my $aligner = $_[1];
    my $callType = $_[2]; #SNV,INDEL or BOTH

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle

    &ProgramPreRequisites($familyID, "SampleCheck", $aligner."/samplecheck", $callType, $FILEHANDLE, 1, 1);
    
    my $inFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $outFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/samplecheck";
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{'pGATKVariantRecalibration'}{'fileEnding'};

    print $FILEHANDLE "#Create Plink .ped and .map file per family using vcfTools","\n";
    print $FILEHANDLE "vcftools ";
    print $FILEHANDLE "--vcf ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf "; #InFile
    print $FILEHANDLE "--plink "; #PLINK format
    print $FILEHANDLE "--out ".$outFamilyDirectory."/".$familyID, "\n\n"; #OutFile (.ped and .map)

    print $FILEHANDLE "#Create vcfTools inbreeding coefficient F per family using vcfTools","\n";
    print $FILEHANDLE "vcftools ";
    print $FILEHANDLE "--vcf ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf "; #InFile
    print $FILEHANDLE "--het "; #Individual inbreeding
    print $FILEHANDLE "--out ".$outFamilyDirectory."/".$familyID, "\n\n"; #Outfile

    if ( ($scriptParameter{'pSampleCheck'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
##Collect QC metadata info for later use                                                                                               
	&SampleInfoQC($scriptParameter{'familyID'}, "noSampleID", "InbreedingFactor", "NoInfile", $outFamilyDirectory, $familyID.".het", "infileDependent"); #"noSampleID is used to select correct keys for %sampleInfo"
    }

    print $FILEHANDLE "#Create Plink .mibs per family","\n"; 
    print $FILEHANDLE "plink ";
    print $FILEHANDLE "--noweb "; #No web check
    print $FILEHANDLE "--ped ".$outFamilyDirectory."/".$familyID.".ped "; #InFile
    print $FILEHANDLE "--map ".$outFamilyDirectory."/".$familyID.".map "; #InFile
    print $FILEHANDLE "--cluster "; #Perform IBS clustering
    print $FILEHANDLE "--matrix "; #Create a N x N matrix of genome-wide average IBS pairwise identities
    print $FILEHANDLE "--out ".$outFamilyDirectory."/".$familyID, "\n\n"; #OutFile

    if ( ($scriptParameter{'pSampleCheck'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
##Collect QC metadata info for later use                                                                                               
	&SampleInfoQC($scriptParameter{'familyID'}, "noSampleID", "RelationCheck", "NoInfile", $outFamilyDirectory, $familyID.".mibs", "infileDependent"); #"noSampleID is used to select correct keys for %sampleInfo"
    }

    print $FILEHANDLE "#Create Plink sexcheck per family","\n"; 
    print $FILEHANDLE "plink ";
    print $FILEHANDLE "--noweb "; #No web check
    print $FILEHANDLE "--ped ".$outFamilyDirectory."/".$familyID.".ped "; #InFile
    print $FILEHANDLE "--map ".$outFamilyDirectory."/".$familyID.".map "; #InFile
    print $FILEHANDLE "--check-sex "; #uses X chromosome data to determine sex (i.e. based on heterozygosity rates) 
    print $FILEHANDLE "--out ".$outFamilyDirectory."/".$familyID, "\n\n"; #OutFile

    if ( ($scriptParameter{'pSampleCheck'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
##Collect QC metadata info for later use                                                                                               
	&SampleInfoQC($scriptParameter{'familyID'}, "noSampleID", "SexCheck", "NoInfile", $outFamilyDirectory, $familyID.".sexcheck", "infileDependent"); #"noSampleID is used to select correct keys for %sampleInfo"
    }
    
    print $FILEHANDLE "wait", "\n\n";    
    
    close($FILEHANDLE); 
    if ( ($scriptParameter{'pSampleCheck'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	&FIDSubmitJob(0, $familyID, 2, $parameter{'pSampleCheck'}{'chain'}, $fileName, 0);
    }
    return;
}

sub QCCollect { 
###Collect qc metrics for this analysis run.

    my $familyID = $_[0]; #familyID NOTE: not sampleid 
    my $aligner = $_[1];
    my $callType = $_[2]; #SNV,INDEL or BOTH

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle

    &ProgramPreRequisites($familyID, "QCCollect", "qccollect", 0, $FILEHANDLE, 1, 1);
    
    my $infile = $scriptParameter{'outDataDir'}."/".$scriptParameter{'familyID'}."/qc_sampleinfo.yaml";
    my $inFamilyDirectory =  $scriptParameter{'outDataDir'}."/".$familyID;
    my $outFamilyDirectory =  $scriptParameter{'outDataDir'}."/".$familyID;

    print $FILEHANDLE "perl ".$scriptParameter{'inScriptDir'}."/qcCollect.pl ";
    print $FILEHANDLE "-sampleInfoFile ".$scriptParameter{'QCCollectSampleInfoFile'}." ";
    print $FILEHANDLE "-regExpFile ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'QCCollectRegExpFile'}." ";
    print $FILEHANDLE "-o ".$outFamilyDirectory."/".$familyID."_qcmetrics.yaml ", "\n\n";     
    
    close($FILEHANDLE); 
    
    if ( ($scriptParameter{'pQCCollect'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	&SampleInfoQC($scriptParameter{'familyID'}, "noSampleID", "QCCollect", "noInfile", $outFamilyDirectory, $familyID."_qcmetrics.yaml", "infileDependent"); #"noSampleID is used to select correct keys for %sampleInfo"
	&FIDSubmitJob(0, $familyID, 2, $parameter{'pQCCollect'}{'chain'}, $fileName, 0);
    }
    return;
}

sub RankVariants { 
###Filter and Rank variants depending on mendelian inheritance, frequency and phenotype
   
    my $familyID = $_[0]; #familyID NOTE: not sampleid 
    my $aligner = $_[1];
    my $callType = $_[2]; #SNV,INDEL or BOTH

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
 
    &ProgramPreRequisites($familyID, "RankVariants", $aligner, $callType, $FILEHANDLE, 1, 5);

    my $inFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{'pAddDepth'}{'fileEnding'};

    print $FILEHANDLE "#Create db master file to select variants from template", "\n";
    my $nrColumns; #Total Nr of columns 
    my $nrAnnotationColumns; #The number of columns containing annotation info
    my $pNrofCol; #For perl regexp
   
    if (-f $scriptParameter{'mergeAnnotatedVariantsDbFile'}) { #locate IDN columns from family specific template file (if defined)
	$pNrofCol = q?perl -nae 'if ($_=~/^outinfo/ || $_=~/^outheaders/ ) { chomp($_); my @nr_of_columns=split(",",$_); print scalar(@nr_of_columns);last; }' ?; #Find the number of columns
	$nrColumns = `$pNrofCol $scriptParameter{'mergeAnnotatedVariantsDbFile'};`; #perl one-liner, inFile and return nr of columns
	$nrAnnotationColumns = $nrColumns - scalar(@sampleIDs);
    }
    elsif (-f $scriptParameter{'referencesDir'}."/".$scriptParameter{'mergeAnnotatedVariantsTemplateFile'}) { #No information on previous intersectCollect to create annovar_merge file - locate IDN columns from unspecific interSect db template file
	$pNrofCol = q?perl -nae 'if ($_=~/^outinfo/ || $_=~/^outheaders/ ) { chomp($_); my @nr_of_columns=split(",",$_); print scalar(@nr_of_columns);last; }' ?;
	$nrAnnotationColumns = `$pNrofCol $scriptParameter{'referencesDir'}/$scriptParameter{'mergeAnnotatedVariantsTemplateFile'};`-1; #"-1" Since IDN is already factored in from the regexp
	$nrColumns = $nrAnnotationColumns + scalar(@sampleIDs);
    }
    elsif (-f $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK/".$familyID.$infileEnding.$callType.".txt") { #Check if the file exists (rerun actual data to sample from) 
	$pNrofCol = q?perl -nae 'if ($_=~/^#/ ) { chomp($_); my @nr_of_columns=split("\t",$_); print scalar(@nr_of_columns);last; }' ?; #Find the number of columns
	$nrColumns = `$pNrofCol $scriptParameter{'outDataDir'}/$familyID/$aligner/GATK/$familyID$infileEnding$callType.txt;`; #perl one-liner, inFile and return nr of columns
	$nrAnnotationColumns = $nrColumns - scalar(@sampleIDs);
    }
    else {
	print STDERR "Could not estimate location of SampleID columns from variant file, nor from templates ('-mergeAnnotatedVariantsDbFile' or '-mergeAnnotatedVariantsTemplateFile'). Please provide this information to run 'pRankVariants'.", "\n";
	exit;
    }
    
    my $sampleIDcolcond = $nrColumns-1; #To write last IDN entry without "," at the end
    
##Add relative path to db_template for variant file(s) 
    my ($volume,$directories,$file) = File::Spec->splitpath($scriptParameter{'outDataDir'});
    my @directories = File::Spec->splitdir($directories);#regExpOutDataFile
    my $regExpOutDataFile;
    for (my $directoryCount=1;$directoryCount<scalar(@directories);$directoryCount++) {
	
	$regExpOutDataFile .= "\\/".$directories[$directoryCount]; #Create escape char for / in later regexp
    }
    $regExpOutDataFile .= $file;
    
##Add relative path to db_template for reference/db files
    ($volume,$directories,$file) = File::Spec->splitpath($scriptParameter{'referencesDir'});
    @directories = File::Spec->splitdir( $directories );
    my $regExpReferenceDirectory;	
    for (my $directoryCount=1;$directoryCount<scalar(@directories);$directoryCount++) {
	
	$regExpReferenceDirectory .= "\\/".$directories[$directoryCount]; #Create escape char for / in later regexp
    }
    $regExpReferenceDirectory .= $file;
    
##Create family specific template
    print $FILEHANDLE q?perl -nae 'if ($_=~/outinfo/i) { if ($_=~/IDN!/) { my $sidstring; for (my $sampleID=?.$nrAnnotationColumns.q?;$sampleID<?.$nrColumns.q?;$sampleID++) { if ($sampleID<?.$sampleIDcolcond.q?) { $sidstring.="IDN_GT_Call=>0_$sampleID,"} else { $sidstring.="IDN_GT_Call=>0_$sampleID"} } s/IDN!/$sidstring/g; print $_;} next;} elsif ($_=~s/FDN\!_|FDN\!/?.$familyID.q?/g) { if($_=~s/^ODF!/?.$regExpOutDataFile.q?/g) {} if($_=~s/ALIGNER!/?.$aligner.q?/g) {} if($_=~s/FILEENDING!_/?.$infileEnding.q?/g) {} if($_=~s/CALLTYPE!/?.$callType.q?/g) {} if ($_=~/IDN!/) { my $sidstring; for (my $sampleID=?.$nrAnnotationColumns.q?;$sampleID<?.$nrColumns.q?;$sampleID++) { if ($sampleID<?.$sampleIDcolcond.q?) { $sidstring.="$sampleID,"} else { $sidstring.="$sampleID"} } s/IDN!/$sidstring/g; print $_;} else { print $_;} } else { if($_=~s/^RD!/?.$regExpReferenceDirectory.q?/g) {} print $_;}' ?;

    print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$scriptParameter{'ImportantDbTemplate'}." "; #Infile
    print $FILEHANDLE "> ".$scriptParameter{'ImportantDbMasterFile'}, "\n\n"; #OutFile
    
    my $haploTypeCallerFile = $inFamilyDirectory."/".$familyID.$infileEnding.$callType.".txt";

###Only Clinically interesting variants
    
    if ( ($scriptParameter{'pGATKHaploTypeCaller'} > 0) || (-f $haploTypeCallerFile) ) { #HaplotypeCaller has been used in present call or previously
	
	print $FILEHANDLE "#Create temp_file containing only clinically interesting variants (to avoid duplicates in ranked list)", "\n";
	print $FILEHANDLE "perl ".$scriptParameter{'inScriptDir'}."/intersectCollect.pl ";
	print $FILEHANDLE "-db ".$scriptParameter{'ImportantDbMasterFile'}." "; #A tab-sep file containing 1 db per line
	if ($humanGenomeReferenceSource eq "hg19") {

	    print $FILEHANDLE "-prechr 1 "; #Use chr prefix in rank script
	}
	print $FILEHANDLE "-sl 1 "; #Select all entries in first infile matching keys in subsequent db files
	print $FILEHANDLE "-s ";

	for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) { 

	    if ($sampleIDCounter eq scalar(@sampleIDs)-1) {
		print $FILEHANDLE $sampleIDs[$sampleIDCounter], " ";
	    }
	    else {
		print $FILEHANDLE $sampleIDs[$sampleIDCounter], ",";
	    }    
	}
	print $FILEHANDLE "-sofs "; #Selected variants and orphan db files out data directory

	for (my $ImportantDbFileOutFilesCounter=0;$ImportantDbFileOutFilesCounter<scalar(@ImportantDbFileOutFiles);$ImportantDbFileOutFilesCounter++) {

	    if ($ImportantDbFileOutFilesCounter eq scalar(@ImportantDbFileOutFiles)-1) {

		print $FILEHANDLE $ImportantDbFileOutFiles[$ImportantDbFileOutFilesCounter]." ","\n\n";
	    }
	    else {

		print $FILEHANDLE $ImportantDbFileOutFiles[$ImportantDbFileOutFilesCounter].",";
	    }
	}
	
###Ranking
	print $FILEHANDLE "#Ranking", "\n";
	
	print $FILEHANDLE "workon ".$scriptParameter{'pythonVirtualEnvironment'}, "\n\n"; #Activate python environment
	
	for (my $ImportantDbFileOutFilesCounter=1;$ImportantDbFileOutFilesCounter<scalar(@ImportantDbFileOutFiles);$ImportantDbFileOutFilesCounter++) { #Skip orphan file and run selected files
	    print $FILEHANDLE "run_mip_family_analysis.py ";
	    print $FILEHANDLE $scriptParameter{'pedigreeFile'}." "; #Pedigree file
	    print $FILEHANDLE "-tres ".$scriptParameter{'rankScore'}." "; #Rank score threshold
	    print $FILEHANDLE $ImportantDbFileOutFiles[$ImportantDbFileOutFilesCounter]." "; #InFile	    
	    ($volume,$directories,$file) = File::Spec->splitpath( $ImportantDbFileOutFiles[$ImportantDbFileOutFilesCounter] ); #Collect outfile directory
	    print $FILEHANDLE "-o ".$directories.$familyID."_ranked_".$callType.".txt", "\n\n"; #OutFile
	    print $FILEHANDLE "wait\n\n";
	    
	    if ( ($scriptParameter{'pRankVariants'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {

		$sampleInfo{$familyID}{'program'}{'RankVariants'}{'Clinical'}{'Path'} = $directories.$familyID."_ranked_".$callType.".txt";  #Save clinical candidate list path
	    }
	}
    }
   
###Research variants
    
##Ranking
    print $FILEHANDLE "#Ranking", "\n"; 
    print $FILEHANDLE "run_mip_family_analysis.py ";
    print $FILEHANDLE $scriptParameter{'pedigreeFile'}." "; #Pedigree file
    print $FILEHANDLE "-tres ".$scriptParameter{'rankScore'}." "; #Rank score threshold
    print $FILEHANDLE $ImportantDbFileOutFiles[0]." "; #InFile	    
    ($volume,$directories,$file) = File::Spec->splitpath( $ImportantDbFileOutFiles[0] ); #Collect outfile directory
    print $FILEHANDLE "-o ".$directories.$familyID."_ranked_".$callType.".txt", "\n\n"; #OutFile
    print $FILEHANDLE "wait\n\n";    
    
    if ( ($scriptParameter{'pRankVariants'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {

	$sampleInfo{$familyID}{'program'}{'RankVariants'}{'Research'}{'Path'} = $directories.$familyID."_ranked_".$callType.".txt";  #Save research candidate list path
    }
        
    for (my $ImportantDbFileOutFilesCounter=0;$ImportantDbFileOutFilesCounter<scalar(@ImportantDbFileOutFiles);$ImportantDbFileOutFilesCounter++) {
	print $FILEHANDLE "rm "; #Remove select files
	print $FILEHANDLE $ImportantDbFileOutFiles[$ImportantDbFileOutFilesCounter], "\n\n";
    }
    
    close($FILEHANDLE);   
    if ( ($scriptParameter{'pRankVariants'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	&FIDSubmitJob(0, $familyID, 1, $parameter{'pRankVariants'}{'chain'}, $fileName, 0);
    }
}

sub AddDp { 
###Adds depth (=DP) for all nonvariants pos for all chr (and subjects) to create a master file containing all annovar information and DP for nonvariants in annovar_merged.txt master file. NOTE: Overwrites current ..annovar_merged.txt file. 

    my $familyID = $_[0]; #familyID NOTE: not sampleid  
    my $aligner = $_[1];
    my $callType = $_[2]; #SNV,INDEL or BOTH 

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $nrCores = &NrofCoresPerSbatch(scalar(@sampleIDs)); #Detect the number of cores to use from number of samplesIDs

    &ProgramPreRequisites($familyID, "AddDepth", $aligner."/GATK", $callType, $FILEHANDLE, $nrCores, 10);
    
    my $inFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $outFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{'pMergeAnnotatedVariants'}{'fileEnding'};
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{'pAddDepth'}{'fileEnding'};
    my $coreCounter=1;

#Find all "./." per sample ID and print chr pos to new file (mpileup -l format)
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) { #For all sample ids, find nonvariants
	
	if ($sampleIDCounter == $coreCounter*$nrCores) { #Using only '$nrCores' cores
	    
	    print $FILEHANDLE "wait", "\n\n";
	    $coreCounter=$coreCounter+1;
	}
	print $FILEHANDLE "#Find all './.' per sampleID and print chrosome position to new file (mpileup -l format)", "\n";
	
	print $FILEHANDLE q?perl -F'\t' -nae' if ($_=~ /?.$sampleIDs[$sampleIDCounter].q?\S+\.\/\./ ) { print "$F[0] $F[1]","\n"; }' ?; #print chromosome and start for sampleID
	print $FILEHANDLE $inFamilyDirectory."/".$familyID.$infileEnding.$callType.".txt "; #InFile
	print $FILEHANDLE "> ".$outFamilyDirectory."/".$sampleIDs[$sampleIDCounter]."_nonvariants.txt &", "\n\n"; #OutFile
    }
    print $FILEHANDLE "wait", "\n\n";
    
    print $FILEHANDLE "#Samples indirectory (BAM-files)", "\n\n"; #Indirectory for sample BAM-files
    
##Find depth (Only proper pairs)
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) { #For all sample ids, find nonvariants
	
	my ($infile, $PicardToolsMergeSwitch) = CheckIfMergedFiles($sampleIDs[$sampleIDCounter]);
	my $sampleIDinfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{'pPicardToolsMarkduplicates'}{'fileEnding'};
	$coreCounter=1; #Reset
	
	if ($sampleIDCounter == $coreCounter*$nrCores) { #Using only '$nrCores' cores 
    
	    print $FILEHANDLE "wait", "\n\n";
	    $coreCounter=$coreCounter+1;
	}
	print $FILEHANDLE "samtools mpileup ";
	print $FILEHANDLE "-A "; #count anomalous read pairs
	print $FILEHANDLE "-l ".$outFamilyDirectory."/".$sampleIDs[$sampleIDCounter]."_nonvariants.txt "; #list of positions (chr pos) or regions (BED)
	
	if ($PicardToolsMergeSwitch == 1) { #Files was merged previously
	    
	    print $FILEHANDLE $scriptParameter{'outDataDir'}."/".$sampleIDs[$sampleIDCounter]."/".$aligner."/".$infile.$sampleIDinfileEnding.".bam "; #InFile (BAM-file)
	}
	else { #No previous merge - list all files at once 
	    
	    for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleIDs[$sampleIDCounter]} });$infileCounter++) { #For all infiles per lane
		
		my $infile = $infilesLaneNoEnding{$sampleIDs[$sampleIDCounter]}[$infileCounter];
		print $FILEHANDLE $scriptParameter{'outDataDir'}."/".$sampleIDs[$sampleIDCounter]."/".$aligner."/".$infile.$sampleIDinfileEnding.".bam "; #InFile (BAM-file)
	    }
	}
	print $FILEHANDLE "| "; #Pipe
	print $FILEHANDLE q?perl -F'\t' -nae' print $F[0],"\t", $F[1],"\t", $F[3], "\n";' ?; #only print chr coordinates 
	print $FILEHANDLE "> ".$outFamilyDirectory."/".$sampleIDs[$sampleIDCounter]."_mpileup_nonvariants.txt &", "\n\n"; #OutFile
	
    }
    print $FILEHANDLE "wait", "\n\n";
    
    print $FILEHANDLE "#Add depth to original file", "\n";
    print $FILEHANDLE "perl ".$scriptParameter{'inScriptDir'}."/add_depth.pl ";
    print $FILEHANDLE "-i ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType.".txt "; #InFile
    if ($humanGenomeReferenceSource eq "hg19") {
	
	print $FILEHANDLE "-prechr 1"; #Use chromosome prefix
    }
    print $FILEHANDLE "-infnv "; #No variant files from mpileup
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {#For all sample ids mpileup nonvariant files
	
	if ($sampleIDCounter eq scalar(@sampleIDs)-1) {
	    
	    print $FILEHANDLE $outFamilyDirectory."/".$sampleIDs[$sampleIDCounter]."_mpileup_nonvariants.txt ";
	}
	else {
	    
	    print $FILEHANDLE $outFamilyDirectory."/".$sampleIDs[$sampleIDCounter]."_mpileup_nonvariants.txt,";	
	}
    }
    print $FILEHANDLE "-s "; #SampleIDs 

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {#For all sample ids mpileup nonvariant files
	
	if ($sampleIDCounter eq scalar(@sampleIDs)-1) {
	    
	    print $FILEHANDLE $sampleIDs[$sampleIDCounter]." ";
	}
	else {
	    print $FILEHANDLE $sampleIDs[$sampleIDCounter].",";	
	    
	}
    }
    print $FILEHANDLE "-o ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType.".txt", "\n\n"; #Overwrites original annovar_merge.txt file
    
    close($FILEHANDLE);   
    if ( ($scriptParameter{'pAddDepth'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	&FIDSubmitJob(0, $familyID, 1, $parameter{'pAddDepth'}{'chain'}, $fileName, 0);
    }
}

sub MergeAnnotatedVariants { 
###Merges (& annotates) all variants for all sampleIDs within family to create a master file containing all annotated information
    
    my $familyID = $_[0]; #familyID NOTE: not sampleid 
    my $aligner = $_[1];
    my $callType = $_[2]; #SNV,INDEL or BOTH 

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $numberOfCores = 1; #Set the number of cores depending on exome/rapid or WGS

    if ($scriptParameter{'analysisType'} eq "genomes") { #WGS analysis
	$numberOfCores = 6; 
    }

    &ProgramPreRequisites($familyID, "MergeAnnotatedVariants", $aligner."/GATK", $callType, $FILEHANDLE, $numberOfCores, 4);

    my $inFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $outFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{'pAnnovar'}{'fileEnding'};
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{'pMergeAnnotatedVariants'}{'fileEnding'};
    
##Create db master file from template
    print $FILEHANDLE "#Create db master file from template", "\n";
    my $sampleIDColumns = scalar(@sampleIDs)+6; #Requires CMMS format (chr,start,stop,ref_allele,alt_allel,GT_Call_Filter,IDN...)
    my $sampleIDcolumnsCondition = scalar(@sampleIDs)+5;
    
##Add relative path to db_template for annovar files 
    my ($volume,$directories,$file) = File::Spec->splitpath($scriptParameter{'outDataDir'});
    my @directories = File::Spec->splitdir($directories);
    my $regExpOutDataFile;
    for (my $directoryCount=1;$directoryCount<scalar(@directories);$directoryCount++) {
       
	$regExpOutDataFile .= "\\/".$directories[$directoryCount]; #Create escape char for / in later regexp
    }
    $regExpOutDataFile .= $file;
    
##Add relative path to db_template for reference files
    ($volume,$directories,$file) = File::Spec->splitpath($scriptParameter{'referencesDir'});
    @directories = File::Spec->splitdir( $directories );
    my $regExpReferenceDirectory;	
    for (my $directoryCount=1;$directoryCount<scalar(@directories);$directoryCount++) {
	
	$regExpReferenceDirectory .= "\\/".$directories[$directoryCount]; #Create escape char for / in later regexp
    }
    $regExpReferenceDirectory .= $file;
    
##Create family specific template
    print $FILEHANDLE q?perl -nae 'if ($_=~/outinfo/i) { if ($_=~/IDN!/) { my $sidstring; for (my $sampleID=6;$sampleID<?.$sampleIDColumns.q?;$sampleID++) { if ($sampleID<?.$sampleIDcolumnsCondition.q?) { $sidstring.="IDN_GT_Call=>0_$sampleID,"} else { $sidstring.="IDN_GT_Call=>0_$sampleID"} } s/IDN!/$sidstring/g; print $_;} next;} elsif ($_=~s/FDN\!_|FDN\!/?.$familyID.q?/g) { if($_=~s/^ODF!/?.$regExpOutDataFile.q?/g) {} if($_=~s/ALIGNER!/?.$aligner.q?/g) {} if($_=~s/FILEENDING!_/?.$infileEnding.q?/g) {} if($_=~s/CALLTYPE!/?.$callType.q?/g) {} if ($_=~/IDN!/) { my $sidstring; for (my $sampleID=6;$sampleID<?.$sampleIDColumns.q?;$sampleID++) { if ($sampleID<?.$sampleIDcolumnsCondition.q?) { $sidstring.="$sampleID,"} else { $sidstring.="$sampleID"} } s/IDN!/$sidstring/g; print $_;} else { print $_;} } else { if($_=~s/^RD!/?.$regExpReferenceDirectory.q?/g) {} print $_;}' ?;
    print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$scriptParameter{'mergeAnnotatedVariantsTemplateFile'}." "; #Infile
    print $FILEHANDLE "> ".$scriptParameter{'mergeAnnotatedVariantsDbFile'}, "\n\n"; #OutFile

    print $FILEHANDLE "perl ".$scriptParameter{'inScriptDir'}."/intersectCollect.pl ";
    print $FILEHANDLE "-db ".$scriptParameter{'mergeAnnotatedVariantsDbFile'}." ";
    if ($humanGenomeReferenceSource eq "hg19") {
	print $FILEHANDLE "-prechr 1 "; #Use chromosome prefix
    }
    print $FILEHANDLE "-o ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.".txt ";
    print $FILEHANDLE "-s ";
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) { 
	if ($sampleIDCounter eq scalar(@sampleIDs)-1) {
	    print $FILEHANDLE $sampleIDs[$sampleIDCounter], "\n\n";
	}
	else {
	    print $FILEHANDLE $sampleIDs[$sampleIDCounter], ",";
	}    
    }
    close($FILEHANDLE);   
    if ( ($scriptParameter{'pMergeAnnotatedVariants'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	&FIDSubmitJob(0, $familyID, 1, $parameter{'pMergeAnnotatedVariants'}{'chain'}, $fileName, 0);
    }
    return;
}

sub GATKVariantEvalExome { 
###GATK VariantEval for exome variants

    my $sampleID = $_[0]; 
    my $aligner = $_[1];
    my $callType = $_[2]; #SNV,INDEL or BOTH 
    my $familyID = $_[3]; 

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle

    &ProgramPreRequisites($sampleID, "GATKVariantEvalExome", $aligner."/GATK/varianteval", $callType, $FILEHANDLE, 1, 2);

    my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/GATK";
    my $outSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/GATK/varianteval";
    my $inFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{'pGATKVariantRecalibration'}{'fileEnding'};
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{'pGATKVariantRecalibration'}{'fileEnding'};
    my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles($sampleID);
    
    if ($PicardToolsMergeSwitch == 1) { #Files was merged previously
    
###GATK SelectVariants
##Select SampleID from familyID vrecal vcf file
	print $FILEHANDLE "#GATK SelectVariants","\n\n";
	print $FILEHANDLE "java -Xmx2g ";
	print $FILEHANDLE "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	print $FILEHANDLE "-l INFO "; #Set the minimum level of logging
	print $FILEHANDLE "-T SelectVariants "; #Type of analysis to run
	print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file
	print $FILEHANDLE "-V: ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf "; #FamilyID inFile
	print $FILEHANDLE "-o ".$outSampleDirectory."/".$infile.$outfileEnding.$callType."_temp.vcf "; #SampleID exome outFile
	print $FILEHANDLE "-sn ".$sampleID, "\n\n"; #Include genotypes from this sample

##Prepp infile to only contain exonic variants

	my $sampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/GATK/varianteval";
	my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{'pMergeAnnotatedVariants'}{'fileEnding'};

	print $FILEHANDLE "grep exon ";
	print $FILEHANDLE $inFamilyDirectory."/".$familyID.$infileEnding.$callType.".txt "; #InFile
	print $FILEHANDLE "> ".$sampleDirectory."/".$sampleID.$infileEnding.$callType."_exonic_variants.txt", "\n\n"; #OutFile

##Intersect exonic variants from created sampleID vcf file (required for GATKVariantEval for exonic variants)
	print $FILEHANDLE "intersectBed ";
	print $FILEHANDLE "-header "; #Print the header from the A file prior to results.
	print $FILEHANDLE "-a ".$sampleDirectory."/".$infile.$outfileEnding.$callType."_temp.vcf "; #SampleID temp exome vcf inFile
	print $FILEHANDLE "-b ".$sampleDirectory."/".$sampleID.$infileEnding.$callType."_exonic_variants.txt "; #SampleID exonic variants
	print $FILEHANDLE "> ".$sampleDirectory."/".$infile.$outfileEnding.$callType."_exome.vcf", "\n\n"; #OutFile (VCF-format)

###VariantEval (exome variants)
	print $FILEHANDLE "#GATK VariantEval","\n\n";
	
	print $FILEHANDLE "java -Xmx2g ";
	print $FILEHANDLE "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	print $FILEHANDLE "-l INFO "; #Set the minimum level of logging
	print $FILEHANDLE "-T VariantEval "; #Type of analysis to run
	print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file
	print $FILEHANDLE "-D ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'GATKVariantEvalDbSNP'}." "; #dbSNP file
	print $FILEHANDLE "-gold ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'GATKVariantEvalGold'}." "; #Evaluations that count calls at sites of true variation (e.g., indel calls) will use this argument as their gold standard for comparison
	print $FILEHANDLE "--eval ".$sampleDirectory."/".$infile.$outfileEnding.$callType."_exome.vcf "; #InFile
	print $FILEHANDLE "-o ".$sampleDirectory."/".$infile.$outfileEnding.$callType."_exome.vcf.varianteval", "\n\n"; #OutFile

##Clean-up temp files
	print $FILEHANDLE "rm ";
	print $FILEHANDLE $sampleDirectory."/".$sampleID.$infileEnding.$callType."_exonic_variants.txt", "\n\n"; #SampleID exonic variants

	print $FILEHANDLE "rm ";
	print $FILEHANDLE $sampleDirectory."/".$infile.$outfileEnding.$callType."_temp.vcf", "\n\n"; #SampleID temp exome vcf inFile

	print $FILEHANDLE "rm ";
	print $FILEHANDLE $sampleDirectory."/".$infile.$outfileEnding.$callType."_temp.vcf.idx", "\n\n"; #SampleID temp exome vcf inFile

	if ( ($scriptParameter{'pGATKVariantEvalExome'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
##Collect QC metadata info for later use                                                                                 
	    &SampleInfoQC($scriptParameter{'familyID'}, $sampleID, "VariantEval_Exome", $infile, $sampleDirectory, $outfileEnding.$callType."_exome.vcf.varianteval", "infileDependent");
	}   
    }
    else { #No previous merge
###GATK SelectVariants
##Select SampleID from familyID vrecal vcf file
	
	for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #For all infiles per lane
	    
	    my $infile = $infilesLaneNoEnding{$sampleID}[$infileCounter];
	    
	    print $FILEHANDLE "#GATK SelectVariants","\n\n";
	    print $FILEHANDLE "java -Xmx2g ";
	    print $FILEHANDLE "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	    print $FILEHANDLE "-l INFO "; #Set the minimum level of logging
	    print $FILEHANDLE "-T SelectVariants "; #Type of analysis to run
	    print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file
	    print $FILEHANDLE "-V: ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf "; #FamilyID infile 
	    print $FILEHANDLE "-o ".$outSampleDirectory."/".$infile.$outfileEnding.$callType."_temp.vcf "; #SampleID outFile
	    print $FILEHANDLE "-sn ".$sampleID, "\n\n"; #Include genotypes from this sample

##Prepp infile to only contain exonic variants

	    my $sampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/GATK/varianteval";
	    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{'pMergeAnnotatedVariants'}{'fileEnding'};
	    
	    print $FILEHANDLE "grep exon ";
	    print $FILEHANDLE $inFamilyDirectory."/".$familyID.$infileEnding.$callType.".txt "; #InFile
	    print $FILEHANDLE "> ".$sampleDirectory."/".$sampleID.$infileEnding.$callType."_exonic_variants.txt", "\n\n"; #OutFile
	    
##Intersect exonic variants from created sampleID vcf file (required for GATKVariantEval for exonic variants)
	    print $FILEHANDLE "intersectBed ";
	    print $FILEHANDLE "-header "; #Print the header from the A file prior to results.
	    print $FILEHANDLE "-a ".$sampleDirectory."/".$infile.$outfileEnding.$callType."_temp.vcf "; #SampleID temp exome vcf inFile
	    print $FILEHANDLE "-b ".$sampleDirectory."/".$sampleID.$infileEnding.$callType."_exonic_variants.txt "; #SampleID exonic variants
	    print $FILEHANDLE "> ".$sampleDirectory."/".$infile.$outfileEnding.$callType."_exome.vcf", "\n\n"; #OutFile (VCF-format)
	    
###VariantEval (exome variants)
	    
	    print $FILEHANDLE "#GATK VariantEval","\n\n";
	    
	    print $FILEHANDLE "java -Xmx2g ";
	    print $FILEHANDLE "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	    print $FILEHANDLE "-l INFO "; #Set the minimum level of logging
	    print $FILEHANDLE "-T VariantEval "; #Type of analysis to run
	    print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file
	    print $FILEHANDLE "-D ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'GATKVariantEvalDbSNP'}." "; #dbSNP file
	    print $FILEHANDLE "-gold ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'GATKVariantEvalGold'}." "; #Evaluations that count calls at sites of true variation (e.g., indel calls) will use this argument as their gold standard for comparison
	    print $FILEHANDLE "--eval ".$sampleDirectory."/".$infile.$outfileEnding.$callType."_exome.vcf "; #InFile
	    print $FILEHANDLE "-o ".$sampleDirectory."/".$infile.$outfileEnding.$callType."_exome.vcf.varianteval", "\n\n"; #OutFile

	    if ( ($scriptParameter{'pGATKVariantEvalExome'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
##Collect QC metadata info for later use                                                                                    
		&SampleInfoQC($scriptParameter{'familyID'}, $sampleID, "VariantEval_Exome", $infile, $sampleDirectory, $outfileEnding.$callType."_exome.vcf.varianteval", "infileDependent");
	    }

##Clean-up temp files
	    print $FILEHANDLE "rm ";
	    print $FILEHANDLE $sampleDirectory."/".$sampleID.$infileEnding.$callType."_exonic_variants.txt", "\n\n"; #SampleID exonic variants
	    
	    print $FILEHANDLE "rm ";
	    print $FILEHANDLE $sampleDirectory."/".$infile.$outfileEnding.$callType."_temp.vcf", "\n\n"; #SampleID temp exome vcf inFile

	    print $FILEHANDLE "rm ";
	    print $FILEHANDLE $sampleDirectory."/".$infile.$outfileEnding.$callType."_temp.vcf.idx", "\n\n"; #SampleID temp exome vcf inFile
	}
    } 
    
    close($FILEHANDLE);   
    if ( ($scriptParameter{'pGATKVariantEvalExome'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	&FIDSubmitJob(0, $familyID, 2, $parameter{'pGATKVariantEvalExome'}{'chain'}, $fileName, 0); #Do not add jobIDs to later jobID{chainkey}
    }
    return;
}

sub GATKVariantEvalAll { 
###GATK VariantEval for all variants

    my $sampleID = $_[0]; 
    my $aligner = $_[1];
    my $callType = $_[2]; #SNV,INDEL or BOTH 
    my $familyID = $_[3]; 

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle

    &ProgramPreRequisites($sampleID, "GATKVariantEvalAll", $aligner."/GATK/varianteval", $callType, $FILEHANDLE, 1, 2);

    my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/GATK";
    my $outSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/GATK/varianteval";
    my $inFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{'pGATKVariantRecalibration'}{'fileEnding'};
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{'pGATKVariantRecalibration'}{'fileEnding'};
    my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles($sampleID);
    
    if ($PicardToolsMergeSwitch == 1) { #Files was merged previously
    
###GATK SelectVariants
##Select SampleID from familyID vrecal vcf file
	print $FILEHANDLE "#GATK SelectVariants","\n\n";
	print $FILEHANDLE "java -Xmx2g ";
	print $FILEHANDLE "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	print $FILEHANDLE "-l INFO "; #Set the minimum level of logging
	print $FILEHANDLE "-T SelectVariants "; #Type of analysis to run
	print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file
	print $FILEHANDLE "-V: ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf "; #FamilyID inFile
	print $FILEHANDLE "-o ".$outSampleDirectory."/".$infile.$outfileEnding.$callType.".vcf "; #SampleID outFile
	print $FILEHANDLE "-sn ".$sampleID, "\n\n"; #Include genotypes from this sample

####VariantEval (all variants)

	my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/GATK/varianteval";

	print $FILEHANDLE "#GATK VariantEval","\n\n";
	
	print $FILEHANDLE "java -Xmx2g ";
	print $FILEHANDLE "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	print $FILEHANDLE "-l INFO "; #Set the minimum level of logging
	print $FILEHANDLE "-T VariantEval "; #Type of analysis to run
	print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file
	print $FILEHANDLE "-D ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'GATKVariantEvalDbSNP'}." "; #dbSNP file
	print $FILEHANDLE "-gold ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'GATKVariantEvalGold'}." "; #Evaluations that count calls at sites of true variation (e.g., indel calls) will use this argument as their gold standard for comparison
	print $FILEHANDLE "--eval ".$inSampleDirectory."/".$infile.$infileEnding.$callType.".vcf "; #InFile
	print $FILEHANDLE "-o ".$outSampleDirectory."/".$infile.$outfileEnding.$callType.".vcf.varianteval", "\n\n"; #OutFile

	if ( ($scriptParameter{'pGATKVariantEvalAll'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
##Collect QC metadata info for later use                                                                                                
	&SampleInfoQC($scriptParameter{'familyID'}, $sampleID, "VariantEval_All", $infile, $outSampleDirectory, $outfileEnding.$callType.".vcf.varianteval", "infileDependent");
	}
    }   
    else { #No previous merge
###GATK SelectVariants
##Select SampleID from familyID vrecal vcf file
	
	for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #For all infiles per lane
	    
	    my $infile = $infilesLaneNoEnding{$sampleID}[$infileCounter];
	    
	    print $FILEHANDLE "#GATK SelectVariants","\n\n";
	    print $FILEHANDLE "java -Xmx2g ";
	    print $FILEHANDLE "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	    print $FILEHANDLE "-l INFO "; #Set the minimum level of logging
	    print $FILEHANDLE "-T SelectVariants "; #Type of analysis to run
	    print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file
	    print $FILEHANDLE "-V: ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf "; #FamilyID infile 
	    print $FILEHANDLE "-o ".$outSampleDirectory."/".$infile.$outfileEnding.$callType.".vcf "; #SampleID outFile
	    print $FILEHANDLE "-sn ".$sampleID, "\n\n"; #Include genotypes from this sample

###VariantEval (all variants)

	    my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/GATK/varianteval";

	    print $FILEHANDLE "#GATK VariantEval","\n\n";
	    
	    print $FILEHANDLE "java -Xmx2g ";
	    print $FILEHANDLE "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	    print $FILEHANDLE "-l INFO "; #Set the minimum level of logging
	    print $FILEHANDLE "-T VariantEval "; #Type of analysis to run
	    print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file
	    print $FILEHANDLE "-D ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'GATKVariantEvalDbSNP'}." "; #dbSNP file
	    print $FILEHANDLE "-gold ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'GATKVariantEvalGold'}." "; #Evaluations that count calls at sites of true variation (e.g., indel calls) will use this argument as their gold standard for comparison
	    print $FILEHANDLE "--eval ".$inSampleDirectory."/".$infile.$infileEnding.$callType.".vcf "; #InFile
	    print $FILEHANDLE "-o ".$outSampleDirectory."/".$infile.$outfileEnding.$callType.".vcf.varianteval", "\n\n"; #OutFile
	
	    if ( ($scriptParameter{'pGATKVariantEvalAll'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
##Collect QC metadata info for later use                                                                             
		&SampleInfoQC($scriptParameter{'familyID'}, $sampleID, "VariantEval_All", $infile, $outSampleDirectory, $outfileEnding.$callType.".vcf.varianteval", "infileDependent");
	    }
	} 
    }
    
    close($FILEHANDLE);   
    if ( ($scriptParameter{'pGATKVariantEvalAll'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	&FIDSubmitJob(0, $familyID, 2, $parameter{'pGATKVariantEvalAll'}{'chain'}, $fileName, 0); #Do not add jobIDs to later jobID{chainkey}
    }
    return;
}

sub Annovar { 
###Filter SNVs by gene, region and databases

    my $familyID = $_[0]; #familyID NOTE: not sampleid 
    my $aligner = $_[1];
    my $callType = $_[2]; #SNV,INDEL or BOTH 

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $nrCores = &NrofCoresPerSbatch(scalar(@annovarTableNames)); #Detect the number of cores to use from @annovarTableNames

    &ProgramPreRequisites( $familyID, "Annovar", $aligner."/GATK", $callType, $FILEHANDLE, $nrCores, 7);

    my $inFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $outFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{'pGATKVariantRecalibration'}{'fileEnding'};
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{'pAnnovar'}{'fileEnding'};
	
    print $FILEHANDLE "#Prepare infile to Annovar format from GATK vcf4", "\n";
    print $FILEHANDLE "perl ".$scriptParameter{'annovarPath'}."/convert2annovar.pl ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf ";
    print $FILEHANDLE "-format vcf4old "; #the format of the input file
    print $FILEHANDLE "-includeinfo "; #specify that the output should contain additional information in the input line
    print $FILEHANDLE "> ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType."_temp", "\n\n"; #Annovar script
    
    print $FILEHANDLE "#Intersect for all samples within familyid and remake file to fit annovar format and subsequent filtering", "\n";
    print $FILEHANDLE q?perl -nae 'my @format; my $formatInfo;chomp($_); if ($_=~/^#/) {print $_;next;} if ($_=~/;set=2/) {} else{ if($F[11] eq "PASS") {} else {$F[11] = "PRES";} @format = split(":",$F[13]); print $F[0], "\t", $F[1], "\t", $F[2], "\t", $F[3], "\t", $F[4], "\t", $F[11], "\t"; ?;
    
    my @sampleIDLexSorts = sort @sampleIDs; #Use lexiographically sorted sample IDNs since GATK HaplotypeCaller/UnifiedGT assigns columns in lexigraphical order. @sampleIDs is not lexiographically sorted if taken straight from the command line. This lex sort ensures that if the user did not supply samples in lex order, there will be no sample column swaping. 
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDLexSorts);$sampleIDCounter++) { #For all sample ids
	
	my $samplecolumn = 14+$sampleIDCounter; #First sample genotype starts at col 14 (start 0, perl). NOTE: Important that samples for HaplotypeCaller/UnifiedGT has same order. Otherwise there will be a sample mix-up.
	
	if ($sampleIDCounter eq scalar(@sampleIDLexSorts)-1) {	#Ensure correct order as long as HaplotypeCAller/UnifiedGT uses lex sort. 
	    print $FILEHANDLE q?print "?.$sampleIDLexSorts[$sampleIDCounter].q?:"; @formatInfo = split(":",$F[?.$samplecolumn.q?]); for (my $formatInfoCounter=0;$formatInfoCounter<scalar(@formatInfo);$formatInfoCounter++) { print "$format[$formatInfoCounter]=$formatInfo[$formatInfoCounter]"; if ( $formatInfoCounter<scalar(@formatInfo)-1 ) {print ":"} } print "\n"; } ?;
	}
	else {
	    print $FILEHANDLE q?print "?.$sampleIDLexSorts[$sampleIDCounter].q?:"; @formatInfo = split(":",$F[?.$samplecolumn.q?]); for (my $formatInfoCounter=0;$formatInfoCounter<scalar(@formatInfo);$formatInfoCounter++) { print "$format[$formatInfoCounter]=$formatInfo[$formatInfoCounter]"; if ( $formatInfoCounter<scalar(@formatInfo)-1 ) {print ":"} } print "\t"; ?;
	}
    }

    print $FILEHANDLE "' ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType."_temp "; #InFile from just created convert2annovar.pl outfile
    print $FILEHANDLE "> ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType, "\n\n"; #OutFile
 
    $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{'pAnnovar'}{'fileEnding'};
    my $coreCounter=1;   	    

    for (my $tableNamesCounter=0;$tableNamesCounter<scalar(@annovarTableNames);$tableNamesCounter++) { #For all specified table names
	
	if ($tableNamesCounter == $coreCounter*$nrCores) { #Using only $nrCores cores
	    
	    print $FILEHANDLE "wait", "\n\n";
	    $coreCounter=$coreCounter+1;
	}
	print $FILEHANDLE "perl ".$scriptParameter{'annovarPath'}."/annotate_variation.pl "; #Annovar script 
	print $FILEHANDLE "-".$annovarTables{$annovarTableNames[$tableNamesCounter]}{'annotation'}." "; #Annotation option	

	if ($annovarTables{$annovarTableNames[$tableNamesCounter]}{'annotation'} eq "geneanno" ) { #Use hgvs output style

	    print $FILEHANDLE "-hgvs ";
	    print $FILEHANDLE "-exonicsplicing "; #Annotate variants near intron/exonic borders
	}
	print $FILEHANDLE "-buildver ".$scriptParameter{'annovarGenomeBuildVersion'}." ";

	if($annovarTables{$annovarTableNames[$tableNamesCounter]}{'dbtype'} eq "generic") {
	    
	    print $FILEHANDLE "-dbtype generic -genericdbfile ".$annovarTables{$annovarTableNames[$tableNamesCounter]}{'file'}[0]." "; #generic db file
	    print $FILEHANDLE "--outfile ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType."_".$annovarTables{$annovarTableNames[$tableNamesCounter]}{'file'}[0]." "; #OutFile
	}	
	else{

	    print $FILEHANDLE "-dbtype ".$annovarTableNames[$tableNamesCounter]." "; #db file
	}
	if ($annovarTableNames[$tableNamesCounter] =~/^1000g/) {#Set MAF TH

	    print $FILEHANDLE "--maf_threshold ".$scriptParameter{'annovarMAFThreshold'}." ";
	}
	if ( ($annovarTableNames[$tableNamesCounter] =~/^snp/) || ($annovarTableNames[$tableNamesCounter] =~/_esp/) ) {#Set MAF TH
	    
	    print $FILEHANDLE "--score_threshold ".$scriptParameter{'annovarMAFThreshold'}." "; #score_threshold since Annovar reserved the maf_threshold for 1000G 
	}
	print $FILEHANDLE $inFamilyDirectory."/".$familyID.$infileEnding.$callType." "; #Infile. Outfile is named using infile prefix except for generic files 
	print $FILEHANDLE $scriptParameter{'annovarPath'}."/humandb &", "\n\n"; #annovar/humandb directory is assumed
    }
    print $FILEHANDLE "wait", "\n\n";
    
    print $FILEHANDLE "rm ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType."_temp", "\n"; #Remove temp file
    close($FILEHANDLE);

    if ( ($scriptParameter{'pAnnovar'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	&FIDSubmitJob(0, $familyID, 1, $parameter{'pAnnovar'}{'chain'}, $fileName, 0);
    }
    return;
}

sub GATKReadBackedPhasing {
###GATK ReadBackedPhasing performs physical phasing of SNP calls, based on sequencing reads. 
	
    my $familyID = $_[0]; #familyID NOTE: not sampleid
    my $aligner = $_[1];
    my $callType = $_[2]; #SNV,INDEL or BOTH
    
    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    &ProgramPreRequisites( $familyID, "GATKReadBackedPhasing", $aligner."/GATK", $callType, $FILEHANDLE, 1, 3);
    
    my $inFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $outFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $infileEnding;
    if ($scriptParameter{'pGATKPhaseByTransmission'} > 0) { 
	$infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{'pGATKPhaseByTransmission'}{'fileEnding'};
    }
    else {
	$infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{'pGATKVariantRecalibration'}{'fileEnding'};
    }
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{'pGATKReadBackedPhasing'}{'fileEnding'};
    
###GATK ReadBackedPhasing
    
    print $FILEHANDLE "\n#GATK ReadBackedPhasing","\n\n";
    print $FILEHANDLE "java -Xmx4g ";
    print $FILEHANDLE "-Djava.io.tmpdir=".$scriptParameter{'GATKTempDirectory'}.'$SLURM_JOB_ID'." "; #Temporary Directory
    print $FILEHANDLE "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
    print $FILEHANDLE "-l INFO "; #Set the minimum level of logging
    print $FILEHANDLE "-T ReadBackedPhasing "; #Type of analysis to run
    print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file
    print $FILEHANDLE "--phaseQualityThresh ".$scriptParameter{'GATKReadBackedPhasingPhaseQualityThresh'}." ";
    if ($scriptParameter{'pGATKPhaseByTransmission'} > 0) { 
	print $FILEHANDLE "-respectPhaseInInput "; #Already phased data - respect calls
    }
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) { #Collect infiles for all sampleIDs
	
	my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleIDs[$sampleIDCounter]."/".$aligner."/GATK";
	my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{'pGATKReduceReads'}{'fileEnding'};
	my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles($sampleIDs[$sampleIDCounter]);
	
	if ($PicardToolsMergeSwitch == 1) { #Alignment BAM-files merged previously
	    
	    print $FILEHANDLE "-I ".$inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile
	}
	else { #No previous merge of alignment BAM-files
	    
	    for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{ $sampleIDs[$sampleIDCounter] } });$infileCounter++) { #For all infiles per lane
		
		my $infile = $infilesLaneNoEnding{ $sampleIDs[$sampleIDCounter] }[$infileCounter];
		
		print $FILEHANDLE "-I ".$inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile(s)
	    }
	}
    } 
    print $FILEHANDLE "-L: ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf "; #Limit to  (family vcf)
    print $FILEHANDLE "-V: ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf "; #InFile (family vcf)
    print $FILEHANDLE "-o ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.".vcf", "\n\n"; #OutFile
 
    print $FILEHANDLE "#Remove Temp Directory\n\n";
    print $FILEHANDLE "rm ";
    print $FILEHANDLE "-rf ".$scriptParameter{'GATKTempDirectory'}.'$SLURM_JOB_ID', "\n\n"; #Remove Temp Directory
   
    close($FILEHANDLE);
    if ( ($scriptParameter{'pGATKReadBackedPhasing'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	&FIDSubmitJob(0,$familyID, 2, $parameter{'pGATKReadBackedPhasing'}{'chain'}, $fileName,0);
    }
    return;
}

sub GATKPhaseByTransmission {
###GATK PhaseByTransmission computes the most likely genotype combination and phases trios and parent/child pairs given their genotype likelihoods and a mutation prior and phases all sites were parent/child transmission can be inferred unambiguously.
    
    my $familyID = $_[0]; #familyID NOTE: not sampleid
    my $aligner = $_[1];
    my $callType = $_[2]; #SNV,INDEL or BOTH
    
    &ProgramPreRequisites( $familyID, "GATKPhaseByTransmission", $aligner."/GATK", $callType, *GATK_PHTR, 1, 3);
    
    my $FamilyFileDirectory = $scriptParameter{'outDataDir'}."/".$familyID;
    my $inFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $outFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{'pGATKVariantRecalibration'}{'fileEnding'};
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{'pGATKPhaseByTransmission'}{'fileEnding'};
    
    unless (-e $FamilyFileDirectory."/".$familyID.".fam") { #Check to see if file already exists
	print GATK_PHTR "#Generating '.fam' file for GATK PhaseByTransmission","\n\n";
	print GATK_PHTR q?perl -nae 'print $F[0], "\t", $F[1], "\t", $F[2], "\t", $F[3], "\t", $F[4], "\t", $F[5], "\n";' ?.$scriptParameter{'pedigreeFile'}." > ".$FamilyFileDirectory."/".$familyID.".fam", "\n\n";
    }
    
###GATK PhaseByTransmission
    
    print GATK_PHTR "\n#GATK PhaseByTransmission","\n\n";
    print GATK_PHTR "java -Xmx4g ";
    print GATK_PHTR "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
    print GATK_PHTR "-l INFO "; #Set the minimum level of logging
    print GATK_PHTR "-T PhaseByTransmission "; #Type of analysis to run
    print GATK_PHTR "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file
    print GATK_PHTR "-V: ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf "; #InFile (family vcf)
    &GATKPedigreeFlag(*GATK_PHTR, $FamilyFileDirectory, "SILENT"); #Passing filehandle directly to sub routine using "*". Sub routine prints "--pedigree file" for family
    print GATK_PHTR "-o ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.".vcf"; #OutFile
    
    close(GATK_PHTR);
    if ( ($scriptParameter{'pGATKPhaseByTransmission'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	&FIDSubmitJob(0, $familyID, 1, $parameter{'pGATKPhaseByTransmission'}{'chain'}, $fileName, 0);
    }
    return;
}

sub GATKVariantReCalibration { 
#GATK VariantRecalibrator/ApplyRecalibration

    my $familyID = $_[0]; #familyID NOTE: not sampleid 
    my $aligner = $_[1];
    my $callType = $_[2]; #SNV,INDEL or BOTH 
    
    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    &ProgramPreRequisites( $familyID, "GATKVariantRecalibration", $aligner."/GATK", $callType, $FILEHANDLE, $scriptParameter{'maximumCores'}, 10);

#Special case
    `mkdir -p $scriptParameter{'outDataDir'}/$familyID/$aligner/GATK/intermediary`; #Creates the aligner folder, GATK data file directory
 
    my $outFamilyFileDirectory = $scriptParameter{'outDataDir'}."/".$familyID;
    my $inFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $outFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{'pGATKHaploTypeCaller'}{'fileEnding'};
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{'pGATKVariantRecalibration'}{'fileEnding'};
    
    unless (-e $scriptParameter{'outDataDir'}."/".$familyID."/".$familyID.".fam") { #Check to see if file already exists
	print $FILEHANDLE "#Generating '.fam' file for GATK VariantRecalibrator/ApplyRecalibration","\n\n";
	print $FILEHANDLE q?perl -nae 'print $F[0], "\t", $F[1], "\t", $F[2], "\t", $F[3], "\t", $F[4], "\t", $F[5], "\n";' ?.$scriptParameter{'pedigreeFile'}." > ".$outFamilyFileDirectory."/".$familyID.".fam", "\n\n";
    }  

    my $contigIntervalListFile = &GATKTargetListFlag($FILEHANDLE);
    
    if ( ($scriptParameter{'analysisType'} eq "rapid") ) { #rapid analysis 
	
###GATK CombineVariants
	
##Needed to include reference exomes to power the building of the probabalistic model. Variants unique to these exomes will be filtered out after varrecal and applyrecal.
	print $FILEHANDLE "\n#GATK CombineVariants","\n\n";
	print $FILEHANDLE "java -Xmx4g ";
	print $FILEHANDLE "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	print $FILEHANDLE "-l INFO "; #Set the minimum level of logging
	print $FILEHANDLE "-T CombineVariants "; #Type of analysis to run
	print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file
	print $FILEHANDLE "-V: ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf "; #InFile (family vcf)
	print $FILEHANDLE "-V: ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'GATKExomeReferenceSNPs'}." "; #Infile (exome reference)
	print $FILEHANDLE "-o ".$outFamilyDirectory."/".$familyID.$infileEnding.$callType."_comb_ref.vcf", "\n\n"; #OutFile
	
    }

###GATK VariantRecalibrator
    
    my $variantRecalibratorOutFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK/intermediary";
    my @modes = ("SNP","INDEL");

    if ( ($scriptParameter{'analysisType'} eq "exomes") || ($scriptParameter{'analysisType'} eq "rapid") ) { #Exome/rapid analysis

	@modes = ("BOTH");
    }

    for (my $modeCounter=0;$modeCounter<scalar(@modes);$modeCounter++) { #SNP and INDEL will be recalibrated successively in the same file because when you specify eg SNP mode, the indels are emitted without modification, and vice-versa. Exome and Rapid will be processed using mode BOTH since there are to few INDELS to use in the recalibration model even though using 30 exome BAMS in Haplotypecaller step. 

	print $FILEHANDLE "\n\n#GATK VariantRecalibrator","\n\n";	
	print $FILEHANDLE "java -Xmx4g ";
	print $FILEHANDLE "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	print $FILEHANDLE "-l INFO "; #Set the minimum level of logging
	print $FILEHANDLE "-T VariantRecalibrator "; #Type of analysis to run
	print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file
	
	if ( ($scriptParameter{'analysisType'} eq "exomes") || ($scriptParameter{'analysisType'} eq "rapid") ) { #Exome/rapid analysis use combined reference for more power
	
	    print $FILEHANDLE "-recalFile ".$variantRecalibratorOutFamilyDirectory."/".$familyID.$infileEnding.$callType."_comb_ref.intervals "; #Recalibration outFile
	    print $FILEHANDLE "-rscriptFile ".$variantRecalibratorOutFamilyDirectory."/".$familyID.$infileEnding.$callType."_comb_ref.intervals.plots.R "; #The output rscript file generated by the VQSR to aid in visualization of the input data and learned model
	    print $FILEHANDLE "-tranchesFile ".$variantRecalibratorOutFamilyDirectory."/".$familyID.$infileEnding.$callType."_comb_ref.intervals.tranches "; #The output tranches file used by ApplyRecalibration
	print $FILEHANDLE "-L ".$contigIntervalListFile." ";#Target list file (merged or original)	
	    
	    if ($scriptParameter{'analysisType'} eq "rapid") {
		    
		print $FILEHANDLE "-input ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType."_comb_ref.vcf "; #Infile just created combined vcf
	    }
	    if ($scriptParameter{'analysisType'} eq "exomes") {
		
		    print $FILEHANDLE "-input ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf "; #Infile HaplotypeCaller combined vcf which used reference BAMs to create combined vcf file
	    }
	}
	else { #WGS
	    
	    print $FILEHANDLE "-recalFile ".$variantRecalibratorOutFamilyDirectory."/".$familyID.$infileEnding.$callType.".intervals ";
	    print $FILEHANDLE "-rscriptFile ".$variantRecalibratorOutFamilyDirectory."/".$familyID.$infileEnding.$callType.".intervals.plots.R ";
	    print $FILEHANDLE "-tranchesFile ".$variantRecalibratorOutFamilyDirectory."/".$familyID.$infileEnding.$callType.".intervals.tranches ";
	    if ($modes[$modeCounter] eq "SNP") {
	
		print $FILEHANDLE "-input ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf ";
	    }
	    if ($modes[$modeCounter] eq "INDEL") {#Use created recalibrated snp vcf as input
	
		print $FILEHANDLE "-input ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.".SNV.vcf ";
	    }
	    print $FILEHANDLE "-an DP "; #The names of the annotations which should used for calculations. NOTE: Not to be used with hybrid capture
	}
	if ( ($modes[$modeCounter] eq "SNP") || ($modes[$modeCounter] eq "BOTH") ) {
	    
	    print $FILEHANDLE "-resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'GATKVariantReCalibrationTrainingSetHapMap'}." "; #A list of sites for which to apply a prior probability of being correct but which aren't used by the algorithm
	    print $FILEHANDLE "-resource:omni,VCF,known=false,training=true,truth=false,prior=12.0 ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'GATKVariantReCalibrationTrainingSet1000GOmni'}." "; #A list of sites for which to apply a prior probability of being correct but which aren't used by the algorithm
	    print $FILEHANDLE "-resource:1000G,known=false,training=true,truth=false,prior=10.0 ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'GATKVariantReCalibrationTrainingSet1000GSNP'}." "; #A list of sites for which to apply a prior probability of being correct but which aren't used by the algorithm
	    print $FILEHANDLE "-an QD "; #The names of the annotations which should used for calculations
	}
	if ( ($modes[$modeCounter] eq "INDEL") || ($modes[$modeCounter] eq "BOTH") ) {

	    print $FILEHANDLE "-resource:mills,VCF,known=true,training=true,truth=true,prior=12.0 ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'GATKVariantReCalibrationTrainingSetMills'}." "; #A list of sites for which to apply a prior probability of being correct but which aren't used by the algorithm
	}
	print $FILEHANDLE "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'GATKVariantReCalibrationTrainingSetDbSNP'}." "; #A list of sites for which to apply a prior probability of being correct but which aren't used by the algorithm
    
	print $FILEHANDLE "-an MQRankSum "; #The names of the annotations which should used for calculations
	print $FILEHANDLE "-an ReadPosRankSum "; #The names of the annotations which should used for calculations
	print $FILEHANDLE "-an FS "; #The names of the annotations which should used for calculations
	print $FILEHANDLE "--mode ".$modes[$modeCounter]." "; #Recalibration mode to employ (SNP|INDEL|BOTH)
	print $FILEHANDLE "-nt ".$scriptParameter{'maximumCores'}." "; #How many data threads should be allocated to running this analysis    
	&GATKPedigreeFlag($FILEHANDLE, $outFamilyFileDirectory, "SILENT"); #Passing filehandle directly to sub routine using "*". Sub routine prints "--pedigree file" for family
	
###GATK ApplyRecalibration
	print $FILEHANDLE "\n\n#GATK ApplyRecalibration","\n\n";
	
	my $applyRecalibrationInFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK/intermediary";
	
	print $FILEHANDLE "java -Xmx3g ";
	print $FILEHANDLE  "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	print $FILEHANDLE "-l INFO "; #Set the minimum level of logging
	print $FILEHANDLE "-T ApplyRecalibration ";
	print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file
	
	if ( ($scriptParameter{'analysisType'} eq "exomes") || ($scriptParameter{'analysisType'} eq "rapid")) { #Exome/rapid analysis use combined reference for more power
	    
	    print $FILEHANDLE "-recalFile ".$applyRecalibrationInFamilyDirectory."/".$familyID.$infileEnding.$callType."_comb_ref.intervals "; #Recalibration outFile
	    print $FILEHANDLE "-tranchesFile ".$applyRecalibrationInFamilyDirectory."/".$familyID.$infileEnding.$callType."_comb_ref.intervals.tranches "; #The output tranches file used by ApplyRecalibration
	    print $FILEHANDLE "-L ".$contigIntervalListFile." ";#Target list file (merged or original)
	
	    if ($scriptParameter{'analysisType'} eq "rapid") {
		
		print $FILEHANDLE "-input ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType."_comb_ref.vcf "; #Infile just created combined vcf
		print $FILEHANDLE "-o ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType."_comb_ref_filtered.vcf ";		    
	    }
	    if ($scriptParameter{'analysisType'} eq "exomes") {
		
		print $FILEHANDLE "-input ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf "; #Infile HaplotypeCaller combined vcf which used reference BAMs to create combined vcf file
		print $FILEHANDLE "-o ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType."_comb_ref_filtered.vcf ";
	    }	   
	}
	else  { #WGS
	    print $FILEHANDLE "-recalFile ".$applyRecalibrationInFamilyDirectory."/".$familyID.$infileEnding.$callType.".intervals ";
	    print $FILEHANDLE "-tranchesFile ".$applyRecalibrationInFamilyDirectory."/".$familyID.$infileEnding.$callType.".intervals.tranches ";
	    if ($modes[$modeCounter] eq "SNP") {
		
		print $FILEHANDLE "-input ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf ";
		print $FILEHANDLE "-o ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.".SNV.vcf ";
	    }
	    if ($modes[$modeCounter] eq "INDEL") {#Use created recalibrated snp vcf as input
	
		print $FILEHANDLE "-input ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.".SNV.vcf ";
		print $FILEHANDLE "-o ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.".vcf ";
	    }
	}
	print $FILEHANDLE "--ts_filter_level ".$scriptParameter{'GATKVariantReCalibrationTSFilterLevel'}." ";
	&GATKPedigreeFlag($FILEHANDLE, $outFamilyFileDirectory, "SILENT"); #Passing filehandle directly to sub routine using "*". Sub routine prints "--pedigree file" for family    
	print $FILEHANDLE "--mode ".$modes[$modeCounter]." "; #Recalibration mode to employ (SNP|INDEL|BOTH)
    }
###GATK SelectVariants

##Removes all genotype information for exome ref and recalulates meta-data info for remaining samples in new file.
    if ( ($scriptParameter{'analysisType'} eq "exomes") || ($scriptParameter{'analysisType'} eq "rapid") ) { #Exome/rapid analysis
	
	print $FILEHANDLE "\n\n#GATK SelectVariants","\n\n";
	print $FILEHANDLE "java -Xmx2g ";
	print $FILEHANDLE  "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	print $FILEHANDLE "-l INFO "; #Set the minimum level of logging
	print $FILEHANDLE "-T SelectVariants "; #Type of analysis to run
	print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file	
	print $FILEHANDLE "-L ".$contigIntervalListFile." ";#Target list file (merged or original)
	print $FILEHANDLE "-env "; #Don't include loci found to be non-variant after the subsetting procedure. 
	print $FILEHANDLE "-V: ".$inFamilyDirectory."/".$familyID.$outfileEnding.$callType."_comb_ref_filtered.vcf "; #InFile
	print $FILEHANDLE "-o ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.".vcf "; #OutFile

	for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) { #For all sampleIDs
		
	    print $FILEHANDLE "-sn ".$sampleIDs[$sampleIDCounter]." "; #Include genotypes from this sample
	}
    }
    
    print $FILEHANDLE "\n\nwait", "\n\n";
    close($FILEHANDLE);   
    	
    if ( ($scriptParameter{'pGATKVariantRecalibration'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
##Collect QC metadata info for later use
	&SampleInfoQC($familyID, "noSampleID", "pedigreeCheck", "NoInfile", $outFamilyDirectory, $familyID.$outfileEnding.$callType.".vcf", "infileDependent"); #"noSampleID is used to select correct keys for %sampleInfo"
	$sampleInfo{ $scriptParameter{'familyID'} }{'MostCompleteVCF'}{'Path'} = $outFamilyDirectory."/".$familyID.$outfileEnding.$callType.".vcf";	
	&FIDSubmitJob(0, $familyID, 1, $parameter{'pGATKVariantRecalibration'}{'chain'}, $fileName, 0);
    }
    return;
}

sub GATKHaplotypeCallerCombineVariants { 
#GATK CombineVariants. Since HaplotypeCaller is presently used per contigs or batches of contigs this module will combine the vcf to 1 file. 

    my $familyID = $_[0]; #familyID NOTE: not sampleid
    my $aligner = $_[1];
    my $callType = $_[2]; #SNV,INDEL or BOTH
    
    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    &ProgramPreRequisites( $familyID, "GATKHaploTypeCallerCombineVariants", $aligner."/GATK", $callType, $FILEHANDLE, 1, 1);
 
    my $inFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK/HaploTypeCaller";
    my $outFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{'pGATKHaploTypeCaller'}{'fileEnding'};
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{'pGATKHaploTypeCaller'}{'fileEnding'};
    
    print $FILEHANDLE "#GATK CombineVariants","\n\n";
    	   
    print $FILEHANDLE "java -Xmx2g ";
    print $FILEHANDLE "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
    print $FILEHANDLE "-l INFO "; #Set the minimum level of logging
    print $FILEHANDLE "-T CombineVariants "; #Type of analysis to run
    print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file

    for (my $chromosomeCounter=0;$chromosomeCounter<scalar(@contigs);$chromosomeCounter++) { #For all chromosome	    
	
	if ( ($scriptParameter{'analysisType'} eq "exomes") || ($scriptParameter{'analysisType'} eq "rapid") ) { #Exome/rapid analysis - Restrict analysis to padded target file(s)
	    
	    unless ($contigs[$chromosomeCounter] =~/MT$|M$/i) { #Do not add MT for exome and rapid samples. NOTE should be determined by target file instead in the future
		
		print $FILEHANDLE "-V ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType."_".$contigs[$chromosomeCounter].".vcf "; #InFiles  
	    }
	}
	else {
	    print $FILEHANDLE "-V ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType."_".$contigs[$chromosomeCounter].".vcf "; #InFiles  
	}	
    }
    print $FILEHANDLE "-o ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.".vcf", "\n\n"; #OutFile

    print $FILEHANDLE "wait", "\n\n";
    close($FILEHANDLE);   
    if ( ($scriptParameter{'pGATKHaploTypeCallerCombineVariants'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {

	$sampleInfo{ $scriptParameter{'familyID'} }{'MostCompleteVCF'}{'Path'} = $outFamilyDirectory."/".$familyID.$outfileEnding.$callType.".vcf";	
	&FIDSubmitJob(0, $familyID, 1, $parameter{'pGATKHaploTypeCallerCombineVariants'}{'chain'}, $fileName, 0);    
    }
    return;
}

sub GATKHaploTypeCaller { 
#GATK HaplotypeCaller 
    
    my $familyID = $_[0]; #familyID NOTE: not sampleid
    my $aligner = $_[1];
    my $callType = $_[2]; #SNV,INDEL or BOTH
    my $chromosome = $_[3]; 
    my $javaHeapAllocation = $_[4];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle

#    &ProgramPreRequisites( $familyID, "GATKHaploTypeCaller", $aligner."/GATK/HaploTypeCaller", $callType, *GATK_HAPCAL, $scriptParameter{'maximumCores'}, 50); #Activate when Haplotypecaller is multithreaded. 
    
    `mkdir -p $scriptParameter{'outDataDir'}/$familyID/$aligner/info;`; #Creates the aligner folder and info data file directory
    `mkdir -p $scriptParameter{'outDataDir'}/$familyID/$aligner/GATK/HaploTypeCaller`; #Creates the aligner folder, GATK data file directory
    `mkdir -p $scriptParameter{'outScriptDir'}/$familyID/$aligner`; #Creates the aligner folder script file directory
    
    if ($scriptParameter{'pGATKHaploTypeCaller'} == 1) {
	$fileName = $scriptParameter{'outScriptDir'}."/".$familyID."/".$aligner."/gatk_haplotypecaller_".$familyID."_".$callType."_chr".$chromosome."."; 
    }
    elsif ($scriptParameter{'pGATKHaploTypeCaller'} == 2) { #Dry run
	$fileName = $scriptParameter{'outScriptDir'}."/".$familyID."/".$aligner."/dry_run_gatk_haplotypecaller_".$familyID."_".$callType."_chr".$chromosome."."; 
	print STDOUT "Dry Run:\n";print MIPLOG  "Dry Run:\n"; 
    }
    &Checkfnexists(\$fileName, \$fnend, \$fileNameTracker);

###Info and Log
    print STDOUT "Creating sbatch script GATK HaplotypeCaller and writing script file(s) to: ".$fileName, "\n";print MIPLOG "Creating sbatch script GATK HaplotypeCaller and writing script file(s) to: ".$fileName, "\n";
    print STDOUT "Sbatch script GATK HaplotypeCaller data files will be written to: ".$scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK/HaplotypeCaller", "\n";print MIPLOG "Sbatch script GATK HaplotypeCaller data files will be written to: ".$scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK/HaploTypeCaller", "\n";
    
    open ($FILEHANDLE, ">".$fileName) or die "Can't write to ".$fileName.":".$!, "\n";
    
    print $FILEHANDLE "#! /bin/bash -l", "\n";
    print $FILEHANDLE "#SBATCH -A ".$scriptParameter{'projectID'}, "\n";
    print $FILEHANDLE "#SBATCH -n 1 ", "\n";	
    print $FILEHANDLE "#SBATCH -t 50:00:00", "\n";
    
    print $FILEHANDLE "#SBATCH -J GATK_HAPCALL_".$familyID."_".$callType."_chr".$chromosome, "\n";
    
    if ($scriptParameter{'pGATKHaploTypeCaller'} == 1) {
	print $FILEHANDLE "#SBATCH -e ".$scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/info/gatk_haplotypecaller_".$familyID."_".$callType."_chr".$chromosome.".".$fileNameTracker.".stderr.txt", "\n";
	print $FILEHANDLE "#SBATCH -o ".$scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/info/gatk_haplotypecaller_".$familyID."_".$callType."_chr".$chromosome.".".$fileNameTracker.".stdout.txt", "\n";
    }
    elsif ($scriptParameter{'pGATKHaploTypeCaller'} == 2) { #Dry run
	print $FILEHANDLE "#SBATCH -e ".$scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/info/dry_run_gatk_haplotypecaller_".$familyID."_".$callType."_chr".$chromosome.".".$fileNameTracker.".stderr.txt", "\n";
	print $FILEHANDLE "#SBATCH -o ".$scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/info/dry_run_gatk_haplotypecaller_".$familyID."_".$callType."_chr".$chromosome.".".$fileNameTracker.".stdout.txt", "\n";
    }
    
    unless ($scriptParameter{'email'} eq 0) {
	print $FILEHANDLE "#SBATCH --mail-type=END", "\n";
	print $FILEHANDLE "#SBATCH --mail-type=FAIL", "\n";
	print $FILEHANDLE "#SBATCH --mail-user=".$scriptParameter{'email'}, "\n\n";
    }
    
    print $FILEHANDLE 'echo "Running on: $(hostname)"',"\n\n";
    
    my $outFamilyFileDirectory = $scriptParameter{'outDataDir'}."/".$familyID;
    my $outFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK/HaploTypeCaller";
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{'pGATKHaploTypeCaller'}{'fileEnding'};

    print $FILEHANDLE "mkdir -p ".$scriptParameter{'GATKTempDirectory'}.'$SLURM_JOB_ID', "\n\n"; 

    my $contigIntervalListFile = &GATKTargetListFlag($FILEHANDLE, \$chromosome);

    if ($chromosome eq 1) { #Only for the first call of subroutine GATK_hapcal.
	
#Generate .fam file for later use in relevant GATK walkers (HaploTypeCaller, VariantscoreRequalibration etc)
	print $FILEHANDLE "#Generating '.fam' file for GATK HaploTypeCaller","\n\n";
	print $FILEHANDLE q?perl -nae 'print $F[0], "\t", $F[1], "\t", $F[2], "\t", $F[3], "\t", $F[4], "\t", $F[5], "\n";' ?.$scriptParameter{'pedigreeFile'}." > ".$outFamilyFileDirectory."/".$familyID.".fam", "\n\n";
    }
    
    print $FILEHANDLE "#GATK HaplotypeCaller","\n\n";
    
    print $FILEHANDLE "java -Xmx".$javaHeapAllocation."g ";

    if ($scriptParameter{'javaUseLargePages'} ne "no") {
	
	    print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
    }
    print $FILEHANDLE "-Djava.io.tmpdir=".$scriptParameter{'GATKTempDirectory'}.'$SLURM_JOB_ID'." "; #Temporary Directory
    print $FILEHANDLE "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
    print $FILEHANDLE "-l INFO "; #Set the minimum level of logging
    print $FILEHANDLE "-T HaplotypeCaller "; #Type of analysis to run
    print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file
    print $FILEHANDLE "-D ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'GATKHaploTypeCallerSNPKnownSet'}." "; #Known SNPs to use for annotation SNPs
    print $FILEHANDLE "-stand_call_conf 30.0 "; #The minimum phred-scaled confidence threshold at which variants should be called
    print $FILEHANDLE "-stand_emit_conf 30.0 "; #The minimum phred-scaled confidence threshold at which variants should be emitted
    print $FILEHANDLE "-nct 8 "; #Number of CPU Threads per data thread
    print $FILEHANDLE "--annotation BaseQualityRankSumTest "; #annotations to apply to variant calls
    print $FILEHANDLE "--annotation ChromosomeCounts "; #annotations to apply to variant calls
    print $FILEHANDLE "--annotation Coverage "; #annotations to apply to variant calls
    print $FILEHANDLE "--annotation FisherStrand "; #annotations to apply to variant calls
    print $FILEHANDLE "--annotation InbreedingCoeff "; #annotations to apply to variant calls
    print $FILEHANDLE "--annotation MappingQualityRankSumTest "; #annotations to apply to variant calls
    print $FILEHANDLE "--annotation MappingQualityZero "; #annotations to apply to variant calls
    print $FILEHANDLE "--annotation QualByDepth "; #annotations to apply to variant calls
    print $FILEHANDLE "--annotation RMSMappingQuality "; #annotations to apply to variant calls
    print $FILEHANDLE "--annotation ReadPosRankSumTest "; #annotations to apply to variant calls
    print $FILEHANDLE "--annotation SpanningDeletions "; #annotations to apply to variant calls
    print $FILEHANDLE "--annotation TandemRepeatAnnotator " ;#annotations to apply to variant calls
    print $FILEHANDLE "--annotation DepthPerAlleleBySample "; #annotations to apply to variant calls
    print $FILEHANDLE "-dcov ".$scriptParameter{'GATKDownSampleToCoverage'}." "; #Coverage to downsample to at any given locus

    if ( ($scriptParameter{'analysisType'} eq "exomes") || ($scriptParameter{'analysisType'} eq "rapid") ) { #Exome/rapid analysis - Restrict analysis to padded target file(s)

	 print $FILEHANDLE "-L ".$contigIntervalListFile." ";#Prints the GATK -L parameter for contif specifc (multiple merged and sorted) interval lists files  
    }

    &GATKPedigreeFlag($FILEHANDLE, $outFamilyFileDirectory, "SILENT"); #Passing filehandle directly to sub routine using "*". Sub routine prints "--pedigree file" for family

    if ($scriptParameter{'analysisType'} eq "exomes") {

	print $FILEHANDLE "-I ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'GATKHaploTypeCallerRefBAMInfile'}." ";
    }

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) { #Collect infiles for all sampleIDs
	
	my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleIDs[$sampleIDCounter]."/".$aligner."/GATK";
	my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{'pGATKReduceReads'}{'fileEnding'};
	my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles($sampleIDs[$sampleIDCounter]);
	
	if ($PicardToolsMergeSwitch == 1) { #Alignment BAM-files merged previously
	    
	    print $FILEHANDLE "-I ".$inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile
	}
	else { #No previous merge of alignment BAM-files
	    
	    for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{ $sampleIDs[$sampleIDCounter] } });$infileCounter++) { #For all infiles per lane
		
		my $infile = $infilesLaneNoEnding{ $sampleIDs[$sampleIDCounter] }[$infileCounter];
		
		print $FILEHANDLE "-I ".$inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile(s)
	    } 
	}
    } 
    print $FILEHANDLE "-o ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType."_".$chromosome.".vcf", "\n\n"; #OutFile

    print $FILEHANDLE "#Remove Temp Directory\n\n";
    print $FILEHANDLE "rm ";
    print $FILEHANDLE "-rf ".$scriptParameter{'GATKTempDirectory'}.'$SLURM_JOB_ID', "\n\n"; #Remove Temp Directory
    
    close($FILEHANDLE);  
    if ( ($scriptParameter{'pGATKHaploTypeCaller'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {

	&FIDSubmitJob(0, $familyID, 3, $parameter{'pGATKHaploTypeCaller'}{'chain'}, $fileName, 0); #Arg2 eq 3 for parallel execution  
    }
    return;
}

sub GATKReduceReads { 
#GATK ReduceReads 
    
    my $sampleID = $_[0];
    my $aligner = $_[1];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle    
    &ProgramPreRequisites( $sampleID, "GATKReduceReads", $aligner."/GATK", 0, $FILEHANDLE, 1, 10);
    
    my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/GATK";
    my $outSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/GATK";
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pGATKBaseRecalibration'}{'fileEnding'};
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pGATKReduceReads'}{'fileEnding'};
    my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles($sampleID);
 
    print $FILEHANDLE "#GATK Reduce Reads","\n\n";

    if ($PicardToolsMergeSwitch == 1) { #Files was merged previously
    
	print $FILEHANDLE "java -Xmx4g ";

	if ($scriptParameter{'javaUseLargePages'} ne "no") {
	    
	    print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	}
	print $FILEHANDLE "-Djava.io.tmpdir=".$scriptParameter{'GATKTempDirectory'}.'$SLURM_JOB_ID'." "; #Temporary Directory
	print $FILEHANDLE "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	print $FILEHANDLE "-l INFO "; #Set the minimum level of logging
	print $FILEHANDLE "-T ReduceReads "; #Type of analysis to run
	print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file	    
	print $FILEHANDLE "-I ".$inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFiles  
	print $FILEHANDLE "-o ".$outSampleDirectory."/".$infile.$outfileEnding.".bam", "\n\n"; #OutFile

	if ( ($scriptParameter{'pGATKReduceReads'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {

	    $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'} = $outSampleDirectory."/".$infile.$outfileEnding.".bam";
	}
    }
    else { #no previous merge
	
	for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #For all infiles per lane
	    
	    my $infile = $infilesLaneNoEnding{$sampleID}[$infileCounter];
	    
	    print $FILEHANDLE "java -Xmx4g ";
	
	    if ($scriptParameter{'javaUseLargePages'} ne "no") {
		
		print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	    }
	    print $FILEHANDLE "-Djava.io.tmpdir=".$scriptParameter{'GATKTempDirectory'}.'$SLURM_JOB_ID'." "; #Temporary Directory
	    print $FILEHANDLE "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	    print $FILEHANDLE "-l INFO "; #Set the minimum level of logging
	    print $FILEHANDLE "-T ReduceReads "; #Type of analysis to run
	    print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file	    
	    print $FILEHANDLE "-I ".$inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFiles  
	    print $FILEHANDLE "-o ".$outSampleDirectory."/".$infile.$outfileEnding.".bam", "\n\n"; #OutFile
	    	
	    if ( ($scriptParameter{'pGATKReduceReads'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {

		$sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'} = $outSampleDirectory."/".$infile.$outfileEnding.".bam";
	    }
	}
    }
    print $FILEHANDLE "#Remove Temp Directory\n\n";
    print $FILEHANDLE "rm ";
    print $FILEHANDLE "-rf ".$scriptParameter{'GATKTempDirectory'}.'$SLURM_JOB_ID', "\n\n"; #Remove Temp Directory
    close($FILEHANDLE);   
    
    if ( ($scriptParameter{'pGATKReduceReads'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	
	&FIDSubmitJob($sampleID, $scriptParameter{'familyID'}, 1, $parameter{'pGATKReduceReads'}{'chain'}, $fileName, 0);    
    }
    return;
}

sub GATKBaseReCalibration { 
#GATK BaseRecalibrator/PrintReads to recalibrate bases before variant calling. Both BaseRecalibrator/PrintReads will be executed within the same sbatch script

    my $sampleID = $_[0];
    my $aligner = $_[1];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    &ProgramPreRequisites($sampleID, "GATKBaseRecalibration", $aligner."/GATK", 0, $FILEHANDLE, $scriptParameter{'maximumCores'}, 50);

#Special case
    `mkdir -p $scriptParameter{'outDataDir'}/$sampleID/$aligner/GATK/intermediary`; #Creates the aligner folder and GATK intermediary data file directory
    
    my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/GATK";
    my $intervalSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/GATK/intermediary";
    my $outSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/GATK";
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pGATKRealigner'}{'fileEnding'};
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pGATKBaseRecalibration'}{'fileEnding'};
    my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles($sampleID);
   
    print $FILEHANDLE "#GATK BaseRecalibrator","\n\n";
    
    if ($PicardToolsMergeSwitch == 1) { #Files was merged previously
       
	print $FILEHANDLE "java -Xmx24g ";

	if ($scriptParameter{'javaUseLargePages'} ne "no") {
	    
	    print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	}
	print $FILEHANDLE "-Djava.io.tmpdir=".$scriptParameter{'GATKTempDirectory'}.'$SLURM_JOB_ID'." "; #Temporary Directory per chr
	print $FILEHANDLE "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	print $FILEHANDLE "-l INFO "; #Set the minimum level of logging
	print $FILEHANDLE "-T BaseRecalibrator "; #Type of analysis to run
	print $FILEHANDLE "-cov ReadGroupCovariate "; #Covariates to be used in the recalibration
	print $FILEHANDLE "-cov ContextCovariate "; #Covariates to be used in the recalibration
	print $FILEHANDLE "-cov CycleCovariate "; #Covariates to be used in the recalibration
	print $FILEHANDLE "-cov QualityScoreCovariate "; #Covariates to be used in the recalibration
	print $FILEHANDLE "-cov ReadGroupCovariate "; #Covariates to be used in the recalibration
	print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file
	print $FILEHANDLE "-knownSites ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'GATKBaseReCalibrationSNPKnownSet'}." ";
	print $FILEHANDLE "-nct ".$scriptParameter{'maximumCores'}." "; #How many CPU threads should be allocated per data thread to running this analysis
	print $FILEHANDLE "-dcov ".$scriptParameter{'GATKDownSampleToCoverage'}." "; #Coverage to downsample to at any given locus
	print $FILEHANDLE "-I ".$inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile
	print $FILEHANDLE "-o ".$intervalSampleDirectory."/".$infile.$infileEnding.".grp ", "\n\n"; #Recalibration table file
	
	print $FILEHANDLE "#GATK PrintReads","\n\n";
	
	print $FILEHANDLE "java -Xmx24g ";

	if ($scriptParameter{'javaUseLargePages'} ne "no") {
	    
	    print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	}
	print $FILEHANDLE "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	print $FILEHANDLE "-l INFO "; #Set the minimum level of logging"-jar $gatk_path/GenomeAnalysisTK.
	print $FILEHANDLE "-T PrintReads "; #Type of analysis to run
	print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file
	print $FILEHANDLE "-nct ".$scriptParameter{'maximumCores'}." "; #How many CPU threads should be allocated per data thread to running this analysis	  
	print $FILEHANDLE "-dcov ".$scriptParameter{'GATKDownSampleToCoverage'}." "; #Coverage to downsample to at any given locus  
	print $FILEHANDLE "-I ".$inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile
	print $FILEHANDLE "-o ".$outSampleDirectory."/".$infile.$outfileEnding.".bam "; #OutFile
	print $FILEHANDLE "-BQSR ".$intervalSampleDirectory."/".$infile.$infileEnding.".grp ", "\n\n"; #Recalibration table file

	if ( ($scriptParameter{'pGATKBaseRecalibration'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) { 
	    
	    $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'} = $outSampleDirectory."/".$infile.$outfileEnding.".bam";
	}
    }
    else { #no previous merge
	
	for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #For all infiles per lane
	    
	    my $infile = $infilesLaneNoEnding{$sampleID}[$infileCounter];
	    
	    print $FILEHANDLE "java -Xmx24g ";
	    
	    if ($scriptParameter{'javaUseLargePages'} ne "no") {
		
		print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	    }
	    print $FILEHANDLE "-Djava.io.tmpdir=".$scriptParameter{'GATKTempDirectory'}.'$SLURM_JOB_ID'." "; #Temporary Directory per chr
	    print $FILEHANDLE "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	    print $FILEHANDLE "-l INFO "; #Set the minimum level of logging
	    print $FILEHANDLE "-T BaseRecalibrator "; #Type of analysis to run
	    print $FILEHANDLE "-cov ReadGroupCovariate "; #Covariates to be used in the recalibration
	    print $FILEHANDLE "-cov ContextCovariate "; #Covariates to be used in the recalibration
	    print $FILEHANDLE "-cov CycleCovariate "; #Covariates to be used in the recalibration
	    print $FILEHANDLE "-cov QualityScoreCovariate "; #Covariates to be used in the recalibration
	    print $FILEHANDLE "-cov ReadGroupCovariate "; #Covariates to be used in the recalibration
	    print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file
	    print $FILEHANDLE "-knownSites ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'GATKBaseReCalibrationSNPKnownSet'}." ";
	    print $FILEHANDLE "-nct ".$scriptParameter{'maximumCores'}." "; #How many CPU threads should be allocated per data thread to running this analysis
	    print $FILEHANDLE "-dcov ".$scriptParameter{'GATKDownSampleToCoverage'}." "; #Coverage to downsample to at any given locus	
	    print $FILEHANDLE "-I ".$inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile
	    print $FILEHANDLE "-o ".$intervalSampleDirectory."/".$infile.$infileEnding.".grp ", "\n\n"; #Recalibration table file

	    print $FILEHANDLE "#GATK PrintReads","\n\n";
	    
	    print $FILEHANDLE "java -Xmx24g ";

	    if ($scriptParameter{'javaUseLargePages'} ne "no") {
		
		print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	    }
	    print $FILEHANDLE "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	    print $FILEHANDLE "-l INFO "; #Set the minimum level of logging"-jar $gatk_path/GenomeAnalysisTK.
	    print $FILEHANDLE "-T PrintReads "; #Type of analysis to run
	    print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file
	    print $FILEHANDLE "-nct ".$scriptParameter{'maximumCores'}." "; #How many CPU threads should be allocated per data thread to running this analysis
	    print $FILEHANDLE "-dcov ".$scriptParameter{'GATKDownSampleToCoverage'}." "; #Coverage to downsample to at any given locus
	    print $FILEHANDLE "-I ".$inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile
	    print $FILEHANDLE "-o ".$outSampleDirectory."/".$infile.$outfileEnding.".bam "; #OutFile
	    print $FILEHANDLE "-BQSR ".$intervalSampleDirectory."/".$infile.$infileEnding.".grp ", "\n\n"; #Recalibration table file

	    if ( ($scriptParameter{'pGATKBaseRecalibration'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) { 

		$sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'} = $outSampleDirectory."/".$infile.$outfileEnding.".bam";
	    }
	}
    }

    print $FILEHANDLE "#Remove Temp Directory\n\n";
    print $FILEHANDLE "rm ";
    print $FILEHANDLE "-rf ".$scriptParameter{'GATKTempDirectory'}.'$SLURM_JOB_ID', "\n\n"; #Remove Temp Directory
    
    close($FILEHANDLE);  
    if ( ($scriptParameter{'pGATKBaseRecalibration'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) { 

	&FIDSubmitJob($sampleID, $scriptParameter{'familyID'}, 1, $parameter{'pGATKBaseRecalibration'}{'chain'}, $fileName,0);
    }
    return;
}

sub GATKReAligner { 
#GATK ReAlignerTargetCreator/IndelRealigner to rearrange reads around INDELs. Both ReAlignerTargetCreator and IndelRealigner will be executed within the same sbatch script

    my $sampleID = $_[0];
    my $aligner = $_[1];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    &ProgramPreRequisites($sampleID, "GATKRealigner", $aligner."/GATK", 0, $FILEHANDLE, $scriptParameter{'maximumCores'}, 40);

#Special case
    `mkdir -p $scriptParameter{'outDataDir'}/$sampleID/$aligner/GATK/intermediary`; #Creates the aligner folder and GATK intermediary data file directory
    
    my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner;
    my $intervalSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/GATK/intermediary";
    my $outSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/GATK";
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pPicardToolsMarkduplicates'}{'fileEnding'};
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pGATKRealigner'}{'fileEnding'};
    my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles($sampleID);

    print $FILEHANDLE "#GATK ReAlignerTargetCreator","\n\n";
    
    if ($PicardToolsMergeSwitch == 1) { #Files was merged previously
	
	print $FILEHANDLE "java -Xmx24g ";

	if ($scriptParameter{'javaUseLargePages'} ne "no") {

	    print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	}
	print $FILEHANDLE "-Djava.io.tmpdir=".$scriptParameter{'GATKTempDirectory'}.'$SLURM_JOB_ID'." "; #Temporary Directory
	print $FILEHANDLE "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	print $FILEHANDLE "-l INFO "; #Set the minimum level of logging
	print $FILEHANDLE "-T RealignerTargetCreator "; #Type of analysis to run
	print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file 
	print $FILEHANDLE "-known ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'GATKReAlignerINDELKnownSet1'}." "; #Input VCF file with known indels
	print $FILEHANDLE "-known ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'GATKReAlignerINDELKnownSet2'}." "; #Input VCF file with known indels
	print $FILEHANDLE "-nt ".$scriptParameter{'maximumCores'}." "; #How many data threads should be allocated to running this analysis.
	print $FILEHANDLE "-dcov ".$scriptParameter{'GATKDownSampleToCoverage'}." "; #Coverage to downsample to at any given locus
	print $FILEHANDLE "-I ".$inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile	    
	print $FILEHANDLE "-o ".$intervalSampleDirectory."/".$infile.$outfileEnding.".intervals ", "\n\n"; #Interval outFile
	
	print $FILEHANDLE "#GATK IndelRealigner","\n\n";
	
	print $FILEHANDLE "java -Xmx24g ";

	if ($scriptParameter{'javaUseLargePages'} ne "no") {
	    
	    print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	}
	print $FILEHANDLE "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	print $FILEHANDLE "-l INFO ";
	print $FILEHANDLE "-T IndelRealigner ";
	print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file
	print $FILEHANDLE "-known ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'GATKReAlignerINDELKnownSet1'}." "; #Input VCF file with known indels
	print $FILEHANDLE "-known ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'GATKReAlignerINDELKnownSet2'}." "; #Input VCF file with known indels	 
	print $FILEHANDLE "-dcov ".$scriptParameter{'GATKDownSampleToCoverage'}." "; #Coverage to downsample to at any given locus
	print $FILEHANDLE "-I ".$inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile	
	print $FILEHANDLE "-o ".$outSampleDirectory."/".$infile.$outfileEnding.".bam "; #OutFile
	print $FILEHANDLE "-targetIntervals ".$intervalSampleDirectory."/".$infile.$outfileEnding.".intervals ", "\n\n";

	if ( ($scriptParameter{'pGATKRealigner'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	    
	    $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'} = $outSampleDirectory."/".$infile.$outfileEnding.".bam";	    
	}	
    }
    else  { #No previous merge
	
	for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #For all infiles per lane
	    
	    my $infile = $infilesLaneNoEnding{$sampleID}[$infileCounter];
		
	    print $FILEHANDLE "java -Xmx24g ";

	    if ($scriptParameter{'javaUseLargePages'} ne "no") {
		
		print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	    }
	    print $FILEHANDLE "-Djava.io.tmpdir=".$scriptParameter{'GATKTempDirectory'}.'$SLURM_JOB_ID'." "; #Temporary Directory
	    print $FILEHANDLE "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	    print $FILEHANDLE "-l INFO "; #Set the minimum level of logging
	    print $FILEHANDLE "-T RealignerTargetCreator "; #Type of analysis to run
	    print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file 
	    print $FILEHANDLE "-known ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'GATKReAlignerINDELKnownSet1'}." "; #Input VCF file with known indels
	    print $FILEHANDLE "-known ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'GATKReAlignerINDELKnownSet2'}." "; #Input VCF file with known indels
	    print $FILEHANDLE "-nt ".$scriptParameter{'maximumCores'}." "; #How many data threads should be allocated to running this analysis.
	    print $FILEHANDLE "-dcov ".$scriptParameter{'GATKDownSampleToCoverage'}." "; #Coverage to downsample to at any given locus	 
	    print $FILEHANDLE "-I ".$inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile
	    print $FILEHANDLE "-o ".$intervalSampleDirectory."/".$infile.$outfileEnding.".intervals ", "\n\n"; #Interval outFile
	    
	    print $FILEHANDLE "#GATK IndelRealigner","\n\n";
	    
	    print $FILEHANDLE "java -Xmx24g ";

	    if ($scriptParameter{'javaUseLargePages'} ne "no") {
		
		print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	    }
	    print $FILEHANDLE "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	    print $FILEHANDLE "-l INFO ";
	    print $FILEHANDLE "-T IndelRealigner ";
	    print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file
	    print $FILEHANDLE "-known ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'GATKReAlignerINDELKnownSet1'}." "; #Input VCF file with known indels
	    print $FILEHANDLE "-known ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'GATKReAlignerINDELKnownSet2'}." "; #Input VCF file with known indels
	    print $FILEHANDLE "-dcov ".$scriptParameter{'GATKDownSampleToCoverage'}." "; #Coverage to downsample to at any given locus
	    print $FILEHANDLE "-I ".$inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile		
	    print $FILEHANDLE "-o ".$outSampleDirectory."/".$infile.$outfileEnding.".bam "; #OutFile
	    print $FILEHANDLE "-targetIntervals ".$intervalSampleDirectory."/".$infile.$outfileEnding.".intervals ", "\n\n";

	    if ( ($scriptParameter{'pGATKRealigner'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
		
		$sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'} = $outSampleDirectory."/".$infile.$outfileEnding.".bam";	    
	    }
	}
    }
    
    print $FILEHANDLE "#Remove Temp Directory\n\n";
    print $FILEHANDLE "rm ";
    print $FILEHANDLE "-rf ".$scriptParameter{'GATKTempDirectory'}.'$SLURM_JOB_ID', "\n\n"; #Remove Temp Directory

    close($FILEHANDLE);
    if ( ($scriptParameter{'pGATKRealigner'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {

	&FIDSubmitJob($sampleID, $scriptParameter{'familyID'}, 1, $parameter{'pGATKRealigner'}{'chain'}, $fileName, 0); 
    }
    return;
}

sub RCoveragePlots { 
#Generates sbatch scripts for R scripts:
#1. covplots_genome.R 
#2. covplots_exome.R
#on files generated from calculateCoverage genomeCoverageBED

    my $sampleID = $_[0]; 
    my $aligner = $_[1];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    &ProgramPreRequisites($sampleID, "RCovPlots", $aligner."/coverageReport", 0, $FILEHANDLE, 1, 1);
   
    my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/coverageReport";
    my $outSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/coverageReport";
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pPicardToolsMarkduplicates'}{'fileEnding'};
    my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles($sampleID);    
    
    if ($PicardToolsMergeSwitch == 1) { #Files was merged previously
	
	if ( defined($scriptParameter{'pGenomeCoverageBED'}) && ($scriptParameter{'pGenomeCoverageBED'} > 0) ) {
	    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pGenomeCoverageBED'}{'fileEnding'};

	    print $FILEHANDLE "Rscript ";
	    print $FILEHANDLE $scriptParameter{'inScriptDir'}."/covplots_genome.R ";
	    print $FILEHANDLE $inSampleDirectory."/".$infile.$infileEnding." "; #InFile
	    print $FILEHANDLE $infile." "; #Sample name
	    print $FILEHANDLE $scriptParameter{'xCoverage'}." "; #X-axis max scale
	    print $FILEHANDLE $outSampleDirectory, " &","\n\n"; #OutFile
	}
    }
    else { #No previous merge
	for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #For all infiles per lane
	    
	    my $infile = $infilesLaneNoEnding{$sampleID}[$infileCounter];
	    
	    if ( defined($scriptParameter{'pGenomeCoverageBED'}) && ($scriptParameter{'pGenomeCoverageBED'} > 0) ) {
		my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pGenomeCoverageBED'}{'fileEnding'};
		
		print $FILEHANDLE "Rscript ";
		print $FILEHANDLE $scriptParameter{'inScriptDir'}."/covplots_genome.R ";
		print $FILEHANDLE $inSampleDirectory."/".$infile.$infileEnding." "; #InFile
		print $FILEHANDLE $infile." "; #Sample name
		print $FILEHANDLE $scriptParameter{'xCoverage'}." "; #X-axis max scale
		print $FILEHANDLE $outSampleDirectory, " &", "\n\n"; #OutFile
	    }
	}
    }
    print $FILEHANDLE "wait", "\n\n";
    close($FILEHANDLE);
    if ( ($scriptParameter{'pRCovPlots'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	&FIDSubmitJob($sampleID, $scriptParameter{'familyID'} , 2, $parameter{'pRCovPlots'}{'chain'}, $fileName, 0);
    }
    return;
}

sub QaCompute { 
#Calculates average chromosome coverage on BAM files. 

    my $sampleID = $_[0]; 
    my $aligner = $_[1]; 

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    
    my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner;
    my $outSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/coverageReport";
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pPicardToolsMarkduplicates'}{'fileEnding'};
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pQaCompute'}{'fileEnding'};
    my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles($sampleID);
    my $coreCounter=1;

    if ($PicardToolsMergeSwitch == 1) { #Files was merged previously
	
	&ProgramPreRequisites($sampleID, "QaCompute", $aligner."/coverageReport", 0, $FILEHANDLE, 1, 4);

	print $FILEHANDLE "qaCompute ";
	print $FILEHANDLE "-m "; #Compute median coverage
	print $FILEHANDLE "-d "; #Print per-chromosome histogram
	print $FILEHANDLE "-i "; #Silent
	print $FILEHANDLE "-c ".$scriptParameter{'xCoverage'}." ";
	print $FILEHANDLE $inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile
	print $FILEHANDLE $outSampleDirectory."/".$infile.$outfileEnding." &", "\n\n"; #OutFile
	
	if ( ($scriptParameter{'pQaCompute'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
##Collect QC metadata info for later use                                                                                              
	    &SampleInfoQC($scriptParameter{'familyID'}, $sampleID, "QaCompute", $infile, $outSampleDirectory, $outfileEnding, "infileDependent");
	}
    }
    else { #No merged files
	
	my $nrCores = &NrofCoresPerSbatch(scalar( @{$lane{$sampleID}} ) ); #Detect the number of cores to from lanes	
	
	&ProgramPreRequisites($sampleID, "QaCompute", $aligner."/coverageReport", 0, $FILEHANDLE, $nrCores, 4);
	
	for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #For all files from MosaikAlign or BWA_Sampe
	    
	    if ($infileCounter == $coreCounter*$nrCores) { #Using only $nrCores
		
		print $FILEHANDLE "wait", "\n\n";
		$coreCounter=$coreCounter+1;
	    }

	    my $infile = $infilesLaneNoEnding{$sampleID}[$infileCounter];	    
	    
	    print $FILEHANDLE "qaCompute ";
	    print $FILEHANDLE "-m "; #Compute median coverage
	    print $FILEHANDLE "-d "; #Print per-chromosome histogram
	    print $FILEHANDLE "-i "; #Silent 
	    print $FILEHANDLE "-c ".$scriptParameter{'xCoverage'}." "; #Max depth to calculate coverage on
	    print $FILEHANDLE $inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile
	    print $FILEHANDLE $outSampleDirectory."/".$infile.$outfileEnding." &", "\n\n"; #OutFile
	    
	    if ( ($scriptParameter{'pQaCompute'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
##Collect QC metadata info for later use                                                                                                
		&SampleInfoQC($scriptParameter{'familyID'}, $sampleID, "QaCompute", $infile, $outSampleDirectory, $outfileEnding, "infileDependent");
	    }
	    
	}
    }
    print $FILEHANDLE "wait", "\n\n";
    
    close($FILEHANDLE);
    if ( ($scriptParameter{'pQaCompute'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	&FIDSubmitJob($sampleID, $scriptParameter{'familyID'}, 1, $parameter{'pQaCompute'}{'chain'}, $fileName, 0);
    }
    return;
}

sub GenomeCoverageBED { 
#Calculates coverage on BAM files. 

    my $sampleID = $_[0]; 
    my $aligner = $_[1]; 
    
    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle

    my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner;
    my $outSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/coverageReport";
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pPicardToolsMarkduplicates'}{'fileEnding'};
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pGenomeCoverageBED'}{'fileEnding'};
    my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles($sampleID);
    my $coreCounter=1;

    if ($PicardToolsMergeSwitch == 1) { #Files was merged previously
	
	&ProgramPreRequisites($sampleID, "GenomeCoverageBED", $aligner."/coverageReport", 0, $FILEHANDLE, 1, 4);
	
	print $FILEHANDLE "genomeCoverageBed ";
	print $FILEHANDLE "-max ".$scriptParameter{'xCoverage'}." "; #Combine all positions with a depth >= max into a single bin in the histogram.
	print $FILEHANDLE "-ibam ".$inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile
	print $FILEHANDLE "> ".$outSampleDirectory."/".$infile.$outfileEnding." ", "\n\n"; #OutFile

    }
    
    else { #No merged files
	
	my $nrCores = &NrofCoresPerSbatch(scalar( @{$lane{$sampleID}} ) ); #Detect the number of cores to from lanes	
	
	&ProgramPreRequisites($sampleID, "GenomeCoverageBED", $aligner."/coverageReport", 0, $FILEHANDLE, $nrCores, 4);
	
	for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #For all files from MosaikAlign or BWA_Sampe
	    
	    if ($infileCounter == $coreCounter*$nrCores) { #Using only $nrCores
		
		print $FILEHANDLE "wait", "\n\n";
		$coreCounter=$coreCounter+1;
	    }

	    my $infile = $infilesLaneNoEnding{$sampleID}[$infileCounter];	    
	    
	    print $FILEHANDLE "genomeCoverageBed ";
	    print $FILEHANDLE "-max ".$scriptParameter{'xCoverage'}." "; #Combine all positions with a depth >= max into a single bin in the histogram.
	    print $FILEHANDLE "-ibam ".$inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile
	    print $FILEHANDLE "> ".$outSampleDirectory."/".$infile.$outfileEnding." &", "\n\n"; #outFile
	}
    }
    print $FILEHANDLE "wait", "\n\n";
    
    close($FILEHANDLE);
    if ( ($scriptParameter{'pGenomeCoverageBED'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	&FIDSubmitJob($sampleID, $scriptParameter{'familyID'}, 1, $parameter{'pGenomeCoverageBED'}{'chain'}, $fileName, 0);
    }
    return;
}

sub PicardToolsCollectMultipleMetrics { 
#Calculates coverage and alignment metrics on BAM files. 

    my $sampleID = $_[0]; 
    my $aligner = $_[1]; 
    
    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle

    my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner;
    my $outSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/coverageReport";
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pPicardToolsMarkduplicates'}{'fileEnding'};
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pPicardToolsMarkduplicates'}{'fileEnding'};
    my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles($sampleID);
    my $coreCounter=1;

    if ($PicardToolsMergeSwitch == 1) { #Files was merged previously
	
	&ProgramPreRequisites($sampleID, "PicardToolsCollectMultipleMetrics", $aligner."/coverageReport", 0, $FILEHANDLE, 1, 4);

	print $FILEHANDLE "java -Xmx4g -jar ".$scriptParameter{'picardToolsPath'}."/CollectMultipleMetrics.jar ";
	print $FILEHANDLE "INPUT=".$inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile
	print $FILEHANDLE "OUTPUT=".$outSampleDirectory."/".$infile.$outfileEnding." "; #OutFile
	print $FILEHANDLE "R=".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." &", "\n\n"; #Reference file
	
	if ( ($scriptParameter{'pPicardToolsCollectMultipleMetrics'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
##Collect QC metadata info for later use                                                                                             
	    &SampleInfoQC($scriptParameter{'familyID'}, $sampleID, "CollectMultipleMetrics", $infile, $outSampleDirectory, $outfileEnding.".alignment_summary_metrics", "infileDependent");
	}
	
    }
    else { #No merged files
	
	my $nrCores = &NrofCoresPerSbatch(scalar( @{$lane{$sampleID}} ) ); #Detect the number of cores to from lanes	
	
	&ProgramPreRequisites($sampleID, "PicardToolsCollectMultipleMetrics", $aligner."/coverageReport", 0, $FILEHANDLE, $nrCores, 4);
	
	for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #For all files from MosaikAlign or BWA_Sampe
	    
	    if ($infileCounter == $coreCounter*$nrCores) { #Using only $nrCores
		
		print $FILEHANDLE "wait", "\n\n";
		$coreCounter=$coreCounter+1;
	    }

	    my $infile = $infilesLaneNoEnding{$sampleID}[$infileCounter];	    
	    
	    print $FILEHANDLE "java -Xmx4g -jar ".$scriptParameter{'picardToolsPath'}."/CollectMultipleMetrics.jar ";
	    print $FILEHANDLE "INPUT=".$inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile
	    print $FILEHANDLE "OUTPUT=".$outSampleDirectory."/".$infile.$outfileEnding." "; #outFile
	    print $FILEHANDLE "R=".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." &", "\n\n"; #Reference file
	    if ( ($scriptParameter{'pPicardToolsCollectMultipleMetrics'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
		##Collect QC metadata info for later use
		&SampleInfoQC($scriptParameter{'familyID'}, $sampleID, "CollectMultipleMetrics", $infile, $outSampleDirectory, $outfileEnding.".alignment_summary_metrics", "infileDependent");
	    }	    
	}
    }
    print $FILEHANDLE "wait", "\n\n";
    
    close($FILEHANDLE);
    if ( ($scriptParameter{'pPicardToolsCollectMultipleMetrics'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	&FIDSubmitJob($sampleID, $scriptParameter{'familyID'}, 1, $parameter{'pPicardToolsCollectMultipleMetrics'}{'chain'}, $fileName, 0);
    }
    return;
}

sub PicardToolsCalculateHSMetrics { 
#Calculates coverage on exonic part of BAM files. 
    
    my $sampleID = $_[0]; 
    my $aligner = $_[1]; 
    
    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
   
    my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner;
    my $outSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/coverageReport";
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pPicardToolsMarkduplicates'}{'fileEnding'};
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pPicardToolsMarkduplicates'}{'fileEnding'};
    my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles($sampleID);
    my $coreCounter=1;
    
    if ($PicardToolsMergeSwitch == 1) { #Files was merged previously
	
	&ProgramPreRequisites($sampleID, "PicardToolsCalculateHSMetrics", $aligner."/coverageReport", 0, $FILEHANDLE, 1, 4);
	
	print $FILEHANDLE "java -Xmx4g -jar ".$scriptParameter{'picardToolsPath'}."/CalculateHsMetrics.jar ";
	print $FILEHANDLE "INPUT=".$inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile
	print $FILEHANDLE "OUTPUT=".$outSampleDirectory."/".$infile.$outfileEnding."_CalculateHsMetrics "; #OutFile
	print $FILEHANDLE "REFERENCE_SEQUENCE=".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file
	print $FILEHANDLE "BAIT_INTERVALS=".$scriptParameter{'referencesDir'}."/".$scriptParameter{ $scriptParameter{'familyID'} }{$sampleID}{'exomeTargetPaddedBedInfileLists'}." "; #Capture kit padded target infile_list file
	print $FILEHANDLE "TARGET_INTERVALS=".$scriptParameter{'referencesDir'}."/".$scriptParameter{ $scriptParameter{'familyID'} }{$sampleID}{'exomeTargetBedInfileLists'}." &", "\n\n"; #Capture kit target infile_list file
	
	if ( ($scriptParameter{'pPicardToolsCalculateHSMetrics'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
##Collect QC metadata info for later use                                                                                   
	    &SampleInfoQC($scriptParameter{'familyID'}, $sampleID, "CalculateHsMetrics", $infile, $outSampleDirectory, $outfileEnding."_CalculateHsMetrics", "infileDependent");
	}
    }
    else { #No merged files
	
	my $nrCores = &NrofCoresPerSbatch(scalar( @{$lane{$sampleID}} ) ); #Detect the number of cores to from lanes	
	
	&ProgramPreRequisites($sampleID, "PicardToolsCalculateHSMetrics", $aligner."/coverageReport", 0, $FILEHANDLE, $nrCores, 4);
	
	for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #For all files from MosaikAlign or BWA_Sampe
	    
	    if ($infileCounter == $coreCounter*$nrCores) { #Using only $nrCores
		
		print $FILEHANDLE "wait", "\n\n";
		$coreCounter=$coreCounter+1;
	    }
	    
	    my $infile = $infilesLaneNoEnding{$sampleID}[$infileCounter];	    
	    
	    print $FILEHANDLE "java -Xmx4g -jar ".$scriptParameter{'picardToolsPath'}."/CalculateHsMetrics.jar ";
	    print $FILEHANDLE "INPUT=".$inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile
	    print $FILEHANDLE "OUTPUT=".$outSampleDirectory."/".$infile.$outfileEnding."_CalculateHsMetrics "; #OutFile
	    print $FILEHANDLE "REFERENCE_SEQUENCE=".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file
	    print $FILEHANDLE "BAIT_INTERVALS=".$scriptParameter{'referencesDir'}."/".$scriptParameter{ $scriptParameter{'familyID'} }{$sampleID}{'exomeTargetPaddedBedInfileLists'}." "; #Capture kit padded target infile_list file
	    print $FILEHANDLE "TARGET_INTERVALS=".$scriptParameter{'referencesDir'}."/".$scriptParameter{ $scriptParameter{'familyID'} }{$sampleID}{'exomeTargetBedInfileLists'}." &", "\n\n"; #Capture kit target infile_list file 
	    
	    if ( ($scriptParameter{'pPicardToolsCalculateHSMetrics'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
##Collect QC metadata info for later use                                                                                                 
		&SampleInfoQC($scriptParameter{'familyID'}, $sampleID, "CalculateHsMetrics", $infile, $outSampleDirectory, $outfileEnding."_CalculateHsMetrics", "infileDependent");	    
	    }
	}
    }
    print $FILEHANDLE "wait", "\n\n";
    
    close($FILEHANDLE);
    if ( ($scriptParameter{'pPicardToolsCalculateHSMetrics'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	&FIDSubmitJob($sampleID, $scriptParameter{'familyID'}, 1, $parameter{'pPicardToolsCalculateHSMetrics'}{'chain'}, $fileName, 0);
    }
    return;
}

sub ChanjoImport { 
##Loads the calculated coverage to family database 

    my $familyID = $_[0]; #familyID NOTE: not sampleid
    my $aligner = $_[1]; #Aligner

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    &ProgramPreRequisites($familyID, "ChanjoImport", "chanjoimport", 0, $FILEHANDLE, 1, 3);

    my $outFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID;

    my $coreCounter=1;

    print $FILEHANDLE "workon ".$scriptParameter{'pythonVirtualEnvironment'}, "\n\n"; #Activate python environment
    
 ##Build family database for coverage report
    print $FILEHANDLE "chanjo ";
    print $FILEHANDLE $outFamilyDirectory."/".$familyID.".sqlite "; #Central Db for family
    print $FILEHANDLE "import ";

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {   
	
	my $sampleID = $sampleIDs[$sampleIDCounter];
	my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/coverageReport";
	my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{'pChanjoCalculate'}{'fileEnding'};
	
	my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles($sampleID);	

	if ($PicardToolsMergeSwitch == 1) { #Files was merged previously
	    print $FILEHANDLE $inSampleDirectory."/".$infile.$infileEnding.".json ";      	
	}
	else {
	    
	    for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #For all files from independent of merged or not
		
		if ($infileCounter == $coreCounter*$scriptParameter{'maximumCores'}) { #Using only $scriptParameter{'maximumCores'} cores
		    
		    print $FILEHANDLE "wait", "\n\n";
		    $coreCounter=$coreCounter+1;
		}
		my $infile = $infilesLaneNoEnding{$sampleID}[$infileCounter];
		print $FILEHANDLE $inSampleDirectory."/".$infile.$infileEnding.".json ";
		
	    }
	}
    }
    print $FILEHANDLE "\n\ndeactivate ", "\n\n"; #Deactivate python environment
    close($FILEHANDLE); 
    if ( ($scriptParameter{'pChanjoImport'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	&FIDSubmitJob(0, $familyID, 5, $parameter{'pChanjoImport'}{'chain'}, $fileName, 0);
    }
    return;
}

sub ChanjoCalculate { 
#Generate coverage json outfile for each individual.

    my $sampleID = $_[0];
    my $aligner = $_[1]; 

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle

    &ProgramPreRequisites($sampleID, "ChanjoCalculate", $aligner."/coverageReport", 0, $FILEHANDLE, 1, 2);      
    
    my $outFamilyDirectory = $scriptParameter{'outDataDir'}."/".$scriptParameter{'familyID'};
    my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner;
    my $outSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/coverageReport";
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pPicardToolsMarkduplicates'}{'fileEnding'};
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pChanjoCalculate'}{'fileEnding'};

    
    my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles($sampleID);
    my $coreCounter=1;	
	
    print $FILEHANDLE "workon ".$scriptParameter{'pythonVirtualEnvironment'}, "\n\n"; #Activate python environment

    if ($PicardToolsMergeSwitch == 1) { #Files was merged previously

	print $FILEHANDLE "chanjo ";
	print $FILEHANDLE "annotate ";
	print $FILEHANDLE $outFamilyDirectory."/".$scriptParameter{'familyID'}.".sqlite "; #Central Db for family
	print $FILEHANDLE "using ";
	print $FILEHANDLE $inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile ; 
	print $FILEHANDLE "--cutoff ".$scriptParameter{'chanjoCalculateCutoff'}." "; #Read depth cutoff
	print $FILEHANDLE "--sample ".$sampleID." "; #SampleID
	print $FILEHANDLE "--splice-sites "; #Include splice sites for every exon
	print $FILEHANDLE "--group ".$scriptParameter{'familyID'}." "; #Group to annotate sample to
	print $FILEHANDLE "--force ";#Overwrite if file outFile exists
	print $FILEHANDLE "--json ".$outSampleDirectory."/".$infile.$outfileEnding.".json". "\n\n"; #OutFile	

	
	if ( ($scriptParameter{'pChanjoCalculate'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {

	    &SampleInfoQC($scriptParameter{'familyID'}, $sampleID, "ChanjoCalculate", $infile, $outSampleDirectory, $outfileEnding.".json", "infileDependent");
	}
    }
    else { #No merged files
	
	for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #For all files from independent of merged or not
	    
	    if ($infileCounter == $coreCounter*$scriptParameter{'maximumCores'}) { #Using only $scriptParameter{'maximumCores'} cores
		
		print $FILEHANDLE "wait", "\n\n";
		$coreCounter=$coreCounter+1;
	    }
	    my $infile = $infilesLaneNoEnding{$sampleID}[$infileCounter];
	    
	    print $FILEHANDLE "chanjo ";
	    print $FILEHANDLE "annotate ";
	    print $FILEHANDLE $outFamilyDirectory."/".$scriptParameter{'familyID'}.".sqlite "; #Central Db for family
	    print $FILEHANDLE "using ";
	    print $FILEHANDLE $inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile ; 
	    print $FILEHANDLE "--cutoff ".$scriptParameter{'chanjoCalculateCutoff'}." "; #Read depth cutoff
	    print $FILEHANDLE "--sample ".$sampleID." "; #SampleID
	    print $FILEHANDLE "--splice-sites "; #Include splice sites for every exon
	    print $FILEHANDLE "--group ".$scriptParameter{'familyID'}." "; #Group to annotate sample to
	    print $FILEHANDLE "--force ";#Overwrite if file outFile exists
	    print $FILEHANDLE "--json ".$outSampleDirectory."/".$infile.$outfileEnding.".json &". "\n\n"; #OutFile   

	    if ( ($scriptParameter{'pChanjoCalculate'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
		
		&SampleInfoQC($scriptParameter{'familyID'}, $sampleID, "ChanjoCalculate", $infile, $outSampleDirectory, $outfileEnding.".json", "infileDependent");
	    }
	}
	print $FILEHANDLE "wait", "\n\n";

    }
    print $FILEHANDLE "deactivate ", "\n\n"; #Deactivate python environment
    close($FILEHANDLE);
    if ( ($scriptParameter{'pChanjoCalculate'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	
	&FIDSubmitJob($sampleID, $scriptParameter{'familyID'}, 5, $parameter{'pChanjoCalculate'}{'chain'}, $fileName, 0);
    }
}

sub ChanjoBuild { 

    my $familyID = $_[0]; #familyID NOTE: not sampleid 

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle

    &ProgramPreRequisites($familyID, "ChanjoBuild", "chanjobuild", 0, $FILEHANDLE, 1, 1);

    my $outFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID;

    print $FILEHANDLE "workon ".$scriptParameter{'pythonVirtualEnvironment'}, "\n\n"; #Activate python environment
    
 ##Build new database
    print $FILEHANDLE "chanjo ";
    print $FILEHANDLE "build ";
    print $FILEHANDLE $outFamilyDirectory."/".$familyID.".sqlite ";
    print $FILEHANDLE "using ";
    print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$scriptParameter{'chanjoBuildDb'}." ";      
    print $FILEHANDLE "--force ", "\n\n";#Overwrite if file outFile exists

    print $FILEHANDLE "deactivate ", "\n\n"; #Deactivate python environment
    close($FILEHANDLE); 

    if ( ($scriptParameter{'pChanjoBuild'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {

	&SampleInfoQC($familyID, "noSampleID", "ChanjoBuild", "NoInfile", $outFamilyDirectory, $familyID.".sqlite", "infileDependent"); #"noSampleID is used to select correct keys for %sampleInfo"
	&FIDSubmitJob(0, $familyID, 5, $parameter{'pChanjoBuild'}{'chain'}, $fileName, 0);
    }
}

sub PicardToolsMarkDuplicates { 
#Mark duplicated reads using PicardTools MarkDuplicates in files generated from alignment (sorted, merged)

    my $sampleID = $_[0];
    my $aligner = $_[1]; 

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $time;
    
    if ($scriptParameter{'pPicardToolsMergeSamFiles'} eq 0) { #If No merge has been performed then time requirements goes down
	$time = 3;
    }
    else{
	$time = ceil(3*scalar( @{ $infilesBothStrandsNoEnding{$sampleID} })); #One full lane on Hiseq takes approx. 3 h to process, round up to nearest full hour.	
    }
    
    my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner;
    my $outSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner;
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pPicardToolsMergeSamFiles'}{'fileEnding'};
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pPicardToolsMarkduplicates'}{'fileEnding'};
    my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles($sampleID);
    my $coreCounter=1;

###
#PicardToolsMarkDuplicates
###
    
    if ($PicardToolsMergeSwitch == 1) { #Files was merged previously

	&ProgramPreRequisites($sampleID, "PicardToolsMarkduplicates", $aligner, 0, $FILEHANDLE, 1, $time);

	print $FILEHANDLE "java -Xmx4g ";
	print $FILEHANDLE "-jar ".$scriptParameter{'picardToolsPath'}."/MarkDuplicates.jar ";
	print $FILEHANDLE "TMP_DIR=".$scriptParameter{'PicardToolsTempDirectory'}.'$SLURM_JOB_ID'." "; #Temp Directory
	print $FILEHANDLE "ASSUME_SORTED=true ";
	print $FILEHANDLE "REMOVE_DUPLICATES=false ";
	print $FILEHANDLE "VALIDATION_STRINGENCY=STRICT ";
	print $FILEHANDLE "INPUT=".$inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile
	print $FILEHANDLE "OUTPUT=".$outSampleDirectory."/".$infile.$outfileEnding.".bam "; #OutFile
	print $FILEHANDLE "METRICS_FILE=".$outSampleDirectory."/".$infile.$outfileEnding."metric ", "\n\n"; #Metric file 
	
        #SamTools index on just created _sorted(_merged)_pmd.bam
	
	print $FILEHANDLE "samtools index ";
	print $FILEHANDLE $outSampleDirectory."/".$infile.$outfileEnding.".bam ","\n\n";
	
	if ( ($scriptParameter{'pPicardToolsMarkduplicates'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
##Collect QC metadata info for later use                       
	    &SampleInfoQC($scriptParameter{'familyID'}, $sampleID, "MarkDuplicates", $infile, $outSampleDirectory, $outfileEnding."metric", "infileDependent");
	    $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'} = $outSampleDirectory."/".$infile.$outfileEnding.".bam";
	}
    }
    else { #No merged files

	my $nrCores = &NrofCoresPerSbatch(scalar( @{$lane{$sampleID}} )); #Detect the number of cores to use from lanes
	
	&ProgramPreRequisites($sampleID, "PicardToolsMarkduplicates", $aligner, 0, $FILEHANDLE, $nrCores, $time);

	for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #For all files from independent of merged or not
	    
	    if ($infileCounter == $coreCounter*$nrCores) { #Using only '$nrCores' cores
		
		print $FILEHANDLE "wait", "\n\n";
		$coreCounter=$coreCounter+1;
	    }
	    my $infile = $infilesLaneNoEnding{$sampleID}[$infileCounter];
	    
	    print $FILEHANDLE "java -Xmx4g ";
	    print $FILEHANDLE "-jar ".$scriptParameter{'picardToolsPath'}."/MarkDuplicates.jar ";
	    print $FILEHANDLE "TMP_DIR=".$scriptParameter{'PicardToolsTempDirectory'}.'$SLURM_JOB_ID'." "; #Temp Directory
	    print $FILEHANDLE "ASSUME_SORTED=true ";
	    print $FILEHANDLE "REMOVE_DUPLICATES=false ";
	    print $FILEHANDLE "VALIDATION_STRINGENCY=STRICT ";
	    print $FILEHANDLE "INPUT=".$inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile
	    print $FILEHANDLE "OUTPUT=".$outSampleDirectory."/".$infile.$outfileEnding.".bam "; #OutFile
	    print $FILEHANDLE "METRICS_FILE=".$outSampleDirectory."/".$infile.$outfileEnding."metric &","\n\n"; #Metric file  
	    
	    if ( ($scriptParameter{'pPicardToolsMarkduplicates'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
##Collect QC metadata info for later use                                             
		&SampleInfoQC($scriptParameter{'familyID'}, $sampleID, "MarkDuplicates", $infile, $outSampleDirectory, $outfileEnding."metric", "infileDependent"); 
		$sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'} = $outSampleDirectory."/".$infile.$outfileEnding.".bam";
	    }
	}    
	
	print $FILEHANDLE "wait", "\n\n";

        #SamTools index on just created _sorted(_merged)_pmd.bam
	for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #For all files from alignment
	    
	    if ($infileCounter == $coreCounter*$nrCores) { #Using only '$nrCores' cores
	    #if ($infileCounter == $coreCounter*$scriptParameter{'maximumCores'}) { #Using only $scriptParameter{'maximumCores'} cores
		
		print $FILEHANDLE "wait", "\n\n";
		$coreCounter=$coreCounter+1;
	    }
	    
	    my $infile = $infilesLaneNoEnding{$sampleID}[$infileCounter];
	    
	    print $FILEHANDLE "samtools index ";
	    print $FILEHANDLE $outSampleDirectory."/".$infile.$outfileEnding.".bam &","\n\n"; #Just created dedupped inFile located in outSamplesDirectory
	    print $FILEHANDLE "wait", "\n\n";	    
	}
    }
    close($FILEHANDLE);
    if ( ($scriptParameter{'pPicardToolsMarkduplicates'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {

	&FIDSubmitJob($sampleID, $scriptParameter{'familyID'}, 1, $parameter{'pPicardToolsMarkduplicates'}{'chain'}, $fileName, 0);
    }
    return;
}

sub PicardToolsMerge { 
#Merges all bam files using PicardTools MergeSamFiles within each sampleid and files generated previously (option if provided with '-picardToolsMergeSamFilesPrevious'). The merged files have to be sorted before attempting to merge.
 
    my $sampleID = $_[0];
    my $aligner = $_[1];
    my $fileEnding = $_[2]; 

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    &ProgramPreRequisites($sampleID, "PicardToolsMergeSamFiles", $aligner, 0, $FILEHANDLE, 1, 20);
  
    my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner;
    my $outSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner;
    my $infileEnding;

    if ($scriptParameter{'analysisType'} ne "rapid") {
	$infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pPicardToolsSortSam'}{'fileEnding'};
    }    
    else { #Rapid mode used
	$infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pPicardToolsMergeRapidReads'}{'fileEnding'};
    }
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pPicardToolsMergeSamFiles'}{'fileEnding'};
    my $lanes = join("",@{$lane{$sampleID}}); #Extract lanes

    if (scalar( @{ $infilesLaneNoEnding{$sampleID} }) > 1) { #Check that we have something to merge and then merge current files before merging with previously merged files

	for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #For all files from 

	    my $infile = $infilesLaneNoEnding{$sampleID}[$infileCounter];

	    if ($infileCounter eq 0) {

		print $FILEHANDLE "java -Xmx4g ";
		print $FILEHANDLE "-jar ".$scriptParameter{'picardToolsPath'}."/MergeSamFiles.jar ";
		print $FILEHANDLE "TMP_DIR=".$scriptParameter{'PicardToolsTempDirectory'}.'$SLURM_JOB_ID'." "; #Temp Directory
		print $FILEHANDLE "OUTPUT=".$outSampleDirectory."/".$sampleID."_lanes_".$lanes.$outfileEnding.".bam "; #OutFile
	    }
	    
	    print $FILEHANDLE "INPUT=".$inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile
	}
	print $FILEHANDLE "\n\n";

	print $FILEHANDLE "samtools index ";
	print $FILEHANDLE $outSampleDirectory."/".$sampleID."_lanes_", @{ $lane{$sampleID} } ,$outfileEnding.".bam", "\n\n"; #InFile using just created merged outfile
	print $FILEHANDLE "wait", "\n\n";

	print $FILEHANDLE "#Remove Temp Directory\n\n";
	print $FILEHANDLE "rm ";
	print $FILEHANDLE "-rf ".$scriptParameter{'PicardToolsTempDirectory'}.'$SLURM_JOB_ID', "\n\n"; #Remove Temp Directory
	
	if ( ($scriptParameter{'pPicardToolsMergeSamFiles'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	    
	    $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'} = $outSampleDirectory."/".$sampleID."_lanes_".$lanes.$outfileEnding.".bam";
	}
    }
    if ( ($sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'picardToolsMergeSamFilesPrevious'} == 1) && (scalar( @{ $infilesLaneNoEnding{$sampleID} }) > 1) ) { #merge previously merged files with merged files generated this run
	
	for (my $mergeFileCounter=0;$mergeFileCounter<scalar(@picardToolsMergeSamFilesPrevious);$mergeFileCounter++) {
	    
	    if ($picardToolsMergeSamFilesPrevious[$mergeFileCounter] =~ /$sampleID/) { #Look for sampleID in previously generated file to be merged with current run to be able to merge correct files within sampleID
		if ($picardToolsMergeSamFilesPrevious[$mergeFileCounter] =~ /lane(\d+)|s_(\d+)/) { #Look for lanes_ or lane\d in previously generated file to be merged with current run to be able to extract previous lanes
		    
		    my $mergeLanes; if($1) {$mergeLanes = $1;} else {$mergeLanes = $2;} #Make sure to always supply lanes from previous regexp		    

		    print $FILEHANDLE "java -Xmx4g ";
		    print $FILEHANDLE "-jar ".$scriptParameter{'picardToolsPath'}."/MergeSamFiles.jar ";
		    print $FILEHANDLE "TMP_DIR=".$scriptParameter{'PicardToolsTempDirectory'}.'$SLURM_JOB_ID'." "; #Temp directory
		    print $FILEHANDLE "OUTPUT=".$outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes.$lanes.$outfileEnding.".bam "; #OutFile
		    print $FILEHANDLE "INPUT=".$inSampleDirectory."/".$sampleID."_lanes_".$lanes.$outfileEnding.".bam "; #InFile
		    print $FILEHANDLE "INPUT=".$picardToolsMergeSamFilesPrevious[$mergeFileCounter], "\n\n"; #$mergeLanes contains lane info on previous merge, $infilesLaneNoEnding{$sampleID}[0] uses @RG for very first .bam file to include read group for subsequent merges. Complete path. 
		    
		    print $FILEHANDLE "samtools index ";
		    print $FILEHANDLE $outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes.$lanes.$outfileEnding.".bam ","\n\n"; #InFile

		    print $FILEHANDLE "#Remove Temp Directory\n\n";
		    print $FILEHANDLE "rm ";
		    print $FILEHANDLE "-rf ".$scriptParameter{'PicardToolsTempDirectory'}.'$SLURM_JOB_ID', "\n\n"; #Remove Temp Directory
		
		    if ( ($scriptParameter{'pPicardToolsMergeSamFiles'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {

			$sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'} = $outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes.$lanes.$outfileEnding.".bam";
		    }
		}
	    }
	}
    }
    elsif ($sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'picardToolsMergeSamFilesPrevious'} == 1) { #merge previously merged files with single file generated this run
	
	for (my $mergeFileCounter=0;$mergeFileCounter<scalar(@picardToolsMergeSamFilesPrevious);$mergeFileCounter++) {
	    
	    if ($picardToolsMergeSamFilesPrevious[$mergeFileCounter] =~ /lane(\d+)|s_(\d+)/) { #Look for lanes_ or lane\d in previously generated file to be merged with current run to be able to extract previous lanes
		
		my $mergeLanes; if($1) {$mergeLanes = $1;} else {$mergeLanes = $2;} #Make sure to always supply lanes from previous regexp
		my $infile = $infilesLaneNoEnding{$sampleID}[0]; #Can only be 1 element in array due to previous if statement		    
		
		print $FILEHANDLE "java -Xmx4g ";
		print $FILEHANDLE "jar ".$scriptParameter{'picardToolsPath'}."/MergeSamFiles.jar ";
		print $FILEHANDLE "TMP_DIR=".$scriptParameter{'PicardToolsTempDirectory'}.'$SLURM_JOB_ID'." "; #Temp Directory
		print $FILEHANDLE "OUTPUT=".$outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes.$lanes.$outfileEnding.".bam "; #OutFile
		print $FILEHANDLE "INPUT=".$inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile
		print $FILEHANDLE "INPUT=".$picardToolsMergeSamFilesPrevious[$mergeFileCounter],"\n\n"; #$mergeLanes contains lane info on previous merge, $infilesLaneNoEnding{$sampleID}[0] uses @RG for very first .bam file to include read group for subsequent merges. Complete path. 
		
		print $FILEHANDLE "samtools index ";
		print $FILEHANDLE $outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes.$lanes.$outfileEnding.".bam", "\n\n"; #InFile
		
		print $FILEHANDLE "#Remove Temp Directory\n\n";
		print $FILEHANDLE "rm ";
		print $FILEHANDLE "-rf ".$scriptParameter{'PicardToolsTempDirectory'}.'$SLURM_JOB_ID', "\n\n"; #Remove Temp Directory

		if ( ($scriptParameter{'pPicardToolsMergeSamFiles'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
		    
		    $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'} = $outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes.$lanes.$outfileEnding.".bam";
		}
	    }
	}
    }
    close($FILEHANDLE);
    if ( ($scriptParameter{'pPicardToolsMergeSamFiles'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	&FIDSubmitJob($sampleID, $scriptParameter{'familyID'}, 1, $parameter{'pPicardToolsMergeSamFiles'}{'chain'}, $fileName, 0);
    }
    return;
}


sub PicardToolsSortSamIndex { 
#Sort and indexes bam files using PicradTools sort and index

    my $sampleID = $_[0];
    my $aligner = $_[1];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $sbatchScriptTracker=0;
    my $time=0;
    my $infileSize;

    for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #For all files

	if ($infile{$sampleID}[$infileCounter] =~/.fastq.gz$/) { #Files are already gz and presently the scalar for compression has not been investigated. Therefore no automatic time allocation can be performed.
	    if ($scriptParameter{'analysisType'} eq "genomes") {
		$time = 25;  
	    }
	    else {
		$time = 15;
	    }
	}
	else { #Files are in fastq format
	    $infileSize = -s $indirpath{$sampleID}."/".$infile{$sampleID}[$infileCounter+$sbatchScriptTracker]; # collect .fastq file size to enable estimation of time required for sort & index, +$sbatchScriptTracker for syncing multiple infiles per sampleID. Hence, filesize will be calculated on read1 (should not matter).	   
	    
	    if ($scriptParameter{'pMosaikBuild'} || $scriptParameter{'pMosaikAlign'} || ($scriptParameter{'aligner'} eq "mosaik")) {
		$time = ceil($infileSize/(1700000*60*60)); #1700000 is a constant calculated from the filesize and time needed for procesing in samtools-0.1.12-10 sort and index and 60*60 is to scale to hours.
	    }
	    if ($scriptParameter{'pBwaAln'} || $scriptParameter{'pBwaSampe'} || ($scriptParameter{'aligner'} eq "bwa")) {
		$time = ceil($infileSize/(1700000*60*60)); #1700000 is a constant calculated from the filesize and time needed for procesing in samtools-0.1.12-10 sort and index and 60*60 is to scale to hours.	    
	    }
	}
	&ProgramPreRequisites($sampleID, "PicardToolsSortSam", $aligner, 0, $FILEHANDLE, 1, $time);
    
###	
#PicardTools Sort and Index
###	
	my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner;
	my $outSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner;
	my $infile = $infilesLaneNoEnding{$sampleID}[$infileCounter];
	my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pPicardToolsSortSam'}{'fileEnding'};

	print $FILEHANDLE "#Sorting the reads\n\n";
	print $FILEHANDLE "java -Xmx4g ";
	print $FILEHANDLE "-jar ".$scriptParameter{'picardToolsPath'}."/SortSam.jar ";
	print $FILEHANDLE "TMP_DIR=".$scriptParameter{'PicardToolsTempDirectory'}.'$SLURM_JOB_ID'." "; #Temp Directory
	print $FILEHANDLE "SORT_ORDER=coordinate"." "; #Sort per contig and coordinate
	print $FILEHANDLE "CREATE_INDEX=TRUE "; #create a BAM index when writing a coordinate-sorted BAM file. 
	print $FILEHANDLE "INPUT=".$inSampleDirectory."/".$infile.".bam "; #InFile
	print $FILEHANDLE "OUTPUT=".$outSampleDirectory."/".$infile.$outfileEnding.".bam", "\n\n"; #Outfile	
	
	print $FILEHANDLE "#Remove Temp Directory\n";
	print $FILEHANDLE "rm ";
	print $FILEHANDLE "-rf ".$scriptParameter{'PicardToolsTempDirectory'}.'$SLURM_JOB_ID', "\n\n"; #Remove Temp Directory

	close($FILEHANDLE);
	if ( ($scriptParameter{'pPicardToolsSortSam'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {

	    $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'} = $outSampleDirectory."/".$infile.$outfileEnding.".bam";
	    &FIDSubmitJob($sampleID, $scriptParameter{'familyID'}, 4, $parameter{'pPicardToolsSortSam'}{'chain'}, $fileName, $sbatchScriptTracker);
	}
	$sbatchScriptTracker++; 
    }
    return;
}

sub BWA_Sampe {
#Alignments of BWA Aln index reads using BWA sampe
    
    my $sampleID = $_[0];
    my $aligner = $_[1];
    
    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $time=0;
    my $infileSize;
    my $pairedEndTracker = 0;

    for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #For all files from BWA aln but process in the same command i.e. both reads per align call

	if ($infile{$sampleID}[$infileCounter] =~/.fastq.gz$/) { #Files are already gz and presently the scalar for compression has not been investigated. Therefore no automatic time allocation can be performed.
	    if ($scriptParameter{'analysisType'} eq "genomes") {
		$time = 40;  
	    }
	    else {
		$time = 20;
	    }
	}
	else { #Files are in fastq format	
	    $infileSize = -s $indirpath{$sampleID}."/".$infile{$sampleID}[$pairedEndTracker]; # collect .fastq file size to enable estimation of time required for aligning.
	    $time = ceil(($infileSize/238)/(3000*60*60)); #238 is a scalar estimating the number of reads depending on filesize. 3500 is the number of reads/s in Bwa_sampe-0.6.1 plus samtools-0.1.12-10 view sam to bam conversion and 60*60 is to scale to hours. (4600 BWA-0.5.9)
	}
	my $sequenceRunMode = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesLaneNoEnding{ $sampleID }[$infileCounter]}{'sequenceRunType'}; #Collect paired-end or single-end sequence run mode

	&ProgramPreRequisites($sampleID, "BwaSampe", $aligner, 0, $FILEHANDLE, 1, $time);
	
	my $BWAinSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/bwa";
	my $FASTQinSampleDirectory = $indirpath{$sampleID};
	my $outSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/bwa";
	my $infile = $infile{$sampleID}[$pairedEndTracker]; #For required .fastq file

#BWA Sampe	
	print $FILEHANDLE "bwa sampe ";
	print $FILEHANDLE "-r ".'"@RG\tID:'.$infilesLaneNoEnding{$sampleID}[$infileCounter].'\tSM:'.$sampleID.'\tPL:ILLUMINA" '.$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #read group header line
	print $FILEHANDLE $BWAinSampleDirectory."/".$infilesBothStrandsNoEnding{$sampleID}[$pairedEndTracker].".sai "; #Read 1

	if ( $sequenceRunMode eq "Paired-end") {
	    $pairedEndTracker = $pairedEndTracker+1; #Increment to collect correct read 2 from %infile
	    print $FILEHANDLE $BWAinSampleDirectory."/".$infilesBothStrandsNoEnding{$sampleID}[$pairedEndTracker].".sai "; #Read 2
	}

	print $FILEHANDLE $FASTQinSampleDirectory."/".$infile." "; #Fastq read 1
	
	if ( $sequenceRunMode eq "Paired-end") { 
	    print $FILEHANDLE $FASTQinSampleDirectory."/".$infile{$sampleID}[$pairedEndTracker]." "; #Fastq read 2
	}

	print $FILEHANDLE "> ".$outSampleDirectory."/".$infilesLaneNoEnding{$sampleID}[$infileCounter].".sam", "\n\n"; #Outfile (SAM)

#Convert SAM to BAM using samTools view	
	print $FILEHANDLE "samtools view -bS ".$BWAinSampleDirectory."/".$infilesLaneNoEnding{$sampleID}[$infileCounter].".sam "; #Infile (SAM)
	print $FILEHANDLE "> ".$outSampleDirectory."/".$infilesLaneNoEnding{$sampleID}[$infileCounter].".bam", "\n\n"; #Outfile (BAM)

#Remove SAM file
	print $FILEHANDLE "Removing temporary SAM-file\n";
	print $FILEHANDLE "rm ".$BWAinSampleDirectory."/".$infilesLaneNoEnding{$sampleID}[$infileCounter].".sam";
		
	close($FILEHANDLE);
	if ( ($scriptParameter{'pBwaSampe'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {

	    $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'} = $outSampleDirectory."/".$infilesLaneNoEnding{$sampleID}[$infileCounter].".bam";
	    &FIDSubmitJob($sampleID, $scriptParameter{'familyID'}, 3, $parameter{'pBwaSampe'}{'chain'}, $fileName, $infileCounter);
	}
	$pairedEndTracker++;
    }
    return;
}

sub BWA_Aln {
#Generates BWA aln index on fastq files
    
    my $sampleID = $_[0];
    my $aligner = $_[1];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $time = ceil(2.5*scalar( @{ $infilesLaneNoEnding{$sampleID} })); #One full lane on Hiseq takes approx. 2,5 h for BWA_Aln to process, round up to nearest full hour.
    my $nrCores = 0;

    for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #For all files   

	if ($sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesLaneNoEnding{ $sampleID }[$infileCounter]}{'sequenceRunType'} eq "Paired-end") { #Second read direction if present
	    $nrCores =  $nrCores + 2; #2 processes per file
	}
	else {#Single-end
	    $nrCores = $nrCores + 1; #Only 1 file and one process
	}
    }

    $nrCores = &NrofCoresPerSbatch($nrCores ); #Make sure that the number of cores does not exceed maximum after incrementing above

    &ProgramPreRequisites($sampleID, "BwaAln", $aligner, 0, $FILEHANDLE, $nrCores, $time);

    my $inSampleDirectory =  $indirpath{$sampleID};
    my $outSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/bwa";
    my $coreCounter=1;    

    for (my $infileCounter=0;$infileCounter<scalar( @{ $infile{$sampleID} });$infileCounter++) {

	if ($infileCounter == $coreCounter*$nrCores) { #Using only nr of cores eq to lanes or maximumCores
	    
	    print $FILEHANDLE "wait", "\n\n";
	    $coreCounter=$coreCounter+1;
	}

	my $tempinfile = $infile{$sampleID}[$infileCounter];

	print $FILEHANDLE "bwa aln ";
	print $FILEHANDLE "-k 1 "; #maximum differences in the seed
	print $FILEHANDLE "-t 4 "; #number of threads
	print $FILEHANDLE "-n 3 "; #max #diff (int) or missing prob under 0.02 err rate (float)
	print $FILEHANDLE "-q ".$scriptParameter{'bwaAlnQualityTrimming'}." "; #Quality trimming
	print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference
	print $FILEHANDLE $inSampleDirectory."/".$tempinfile." "; #InFile
	print $FILEHANDLE "> ".$outSampleDirectory."/".$infilesBothStrandsNoEnding{$sampleID}[$infileCounter].".sai &", "\n\n"; #OutFile 
    }
    print $FILEHANDLE "wait", "\n\n";
    close($FILEHANDLE);
    if ( ($scriptParameter{'pBwaAln'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {   

	&FIDSubmitJob($sampleID, $scriptParameter{'familyID'}, 1, $parameter{'pBwaAln'}{'chain'}, $fileName, 0);
    }
    return;
}

sub PicardToolsMergeRapidReads { 
#Merges all batch read processes to one file using PicardTools MergeSamFiles within each sampleid. The read batch proccessed files have to be sorted before attempting to merge.
 
    my $sampleID = $_[0];
    my $aligner = $_[1];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    &ProgramPreRequisites($sampleID, "PicardToolsMergeRapidReads", $aligner, 0, $FILEHANDLE, $scriptParameter{'maximumCores'}, 20);
  
    my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner;
    my $outSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner;
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pBwaMem'}{'fileEnding'};
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pPicardToolsMergeRapidReads'}{'fileEnding'};
    my $coreCounter=1;
    my $coreTracker=0; #Required to portion out cores and files before wait and to track the MOS_BU outfiles to correct lane
    
    for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #For all files from 
	
	my $infile = $infilesLaneNoEnding{$sampleID}[$infileCounter];
	my $nrReadBatchProcesses = $sampleInfo{$scriptParameter{'familyID'}}{$sampleID}{$infilesLaneNoEnding{$sampleID}[$infileCounter]}{'pBwaMem'}{'ReadBatchProcesses'}; 

	if ($nrReadBatchProcesses > 0) { #Check that we have read batch processes to merge

	    if ($coreTracker == $coreCounter*$scriptParameter{'maximumCores'}) { #Using only $scriptParameter{'maximumCores'} nr of cores
		
		print $FILEHANDLE "wait", "\n\n";
		$coreCounter=$coreCounter+1;
	    }

	    for (my $readBatchProcessesCount=0;$readBatchProcessesCount<$nrReadBatchProcesses;$readBatchProcessesCount++) {
		
		if ($readBatchProcessesCount eq 0) {
		    
		    print $FILEHANDLE "java -Xmx4g ";
		    print $FILEHANDLE "-jar ".$scriptParameter{'picardToolsPath'}."/MergeSamFiles.jar ";
		    print $FILEHANDLE "TMP_DIR=".$scriptParameter{'PicardToolsTempDirectory'}.'$SLURM_JOB_ID'." "; #Temp Directory
		    print $FILEHANDLE "OUTPUT=".$outSampleDirectory."/".$infilesLaneNoEnding{$sampleID}[$infileCounter].$outfileEnding.".bam "; #OutFile
		}
		print $FILEHANDLE "INPUT=".$inSampleDirectory."/".$infilesLaneNoEnding{$sampleID}[$infileCounter]."_".$readBatchProcessesCount."_sorted.bam "; #InFile(s)
	    }
	    print $FILEHANDLE "CREATE_INDEX=TRUE &"; #Create a BAM index when writing a coordinate-sorted BAM file.
	    print $FILEHANDLE "\n\n";
	    $coreTracker++; #Track nr of merge calls for infiles so that wait can be printed at the correct intervals (dependent on $scriptParameter{'maximumCores'})
	}
	else { #Still needs to rename file to be included in potential merge of BAM files in next step
	    
	    print $FILEHANDLE "java -Xmx4g ";
	    print $FILEHANDLE "-jar ".$scriptParameter{'picardToolsPath'}."/MergeSamFiles.jar ";
	    print $FILEHANDLE "TMP_DIR=".$scriptParameter{'PicardToolsTempDirectory'}.'$SLURM_JOB_ID'." "; #Temp Directory
	    print $FILEHANDLE "OUTPUT=".$outSampleDirectory."/".$infilesLaneNoEnding{$sampleID}[$infileCounter].$outfileEnding.".bam "; #OutFile
	    
	    print $FILEHANDLE "INPUT=".$inSampleDirectory."/".$infilesLaneNoEnding{$sampleID}[$infileCounter]."_0_sorted_rg.bam "; #InFile
	    print $FILEHANDLE "CREATE_INDEX=TRUE &"; #Create a BAM index when writing a coordinate-sorted BAM file.
	    print $FILEHANDLE "\n\n";
	}
    }
    print $FILEHANDLE "wait", "\n\n";
    
    print $FILEHANDLE "#Remove Temp Directory\n\n";
    print $FILEHANDLE "rm ";
    print $FILEHANDLE "-rf ".$scriptParameter{'PicardToolsTempDirectory'}.'$SLURM_JOB_ID', "\n\n"; #Remove Temp Directory
    
    close($FILEHANDLE);
    if ( ($scriptParameter{'pPicardToolsMergeRapidReads'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {

	&FIDSubmitJob($sampleID, $scriptParameter{'familyID'}, 1, $parameter{'pPicardToolsMergeRapidReads'}{'chain'}, $fileName, 0); #0 since it is only 1 file that is handled in parallel.
    }
    return;
}


sub BWA_Mem {
###Alignment using BWA Mem
##Development Note: Keep number of nodes att 150000000, but increase the size of the read batch
   
    my $sampleID = $_[0];
    my $aligner = $_[1];
 
    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $infileSize;
    my $totalSbatchCounter = 0;
    my $pairedEndTracker = 0;

    for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #For all infiles but process in the same command i.e. both reads per align call
	
	my $sequenceRunMode = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesLaneNoEnding{ $sampleID }[$infileCounter]}{'sequenceRunType'}; #Collect paired-end or single-end sequence run mode
	
	if ($infile{$sampleID}[$infileCounter] =~/.fastq.gz$/) { #Files are already gz and presently the scalar for compression has not been investigated. Therefore no automatic time allocation can be performed.
	
	    if ($sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesLaneNoEnding{ $sampleID }[$infileCounter]}{'sequenceRunType'} eq "Paired-end") { #Second read direction if present
                $infileSize = -s $indirpath{$sampleID}."/".$infile{$sampleID}[$infileCounter+$infileCounter];
	    }
	    else { #Single-end
                $infileSize = -s $indirpath{$sampleID}."/".$infile{$sampleID}[$infileCounter];
	    }
        }
        else { #Files are in fastq format
	    
	    if ($sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesLaneNoEnding{ $sampleID }[$infileCounter]}{'sequenceRunType'} eq "Paired-end") { #Second read direction if present        
		$infileSize = -s $indirpath{$sampleID}."/".$infile{$sampleID}[$infileCounter+$infileCounter]; # collect .fastq file size to enable estimation of time required for aligning, +1 for syncing multiple infiles per sampleID. Hence, filesize will be calculated on read2 (should not matter).
	    }
	    else { #Single-end
                $infileSize = -s $indirpath{$sampleID}."/".$infile{$sampleID}[$infileCounter];
	    }
        }
	
	if ($scriptParameter{'analysisType'} eq "rapid") {
	    
	    my $seqLength = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesBothStrandsNoEnding{ $sampleID }[$infileCounter]}{'sequenceLength'};
	    my ($numberNodes, $ReadNrofLines) = &DetermineNrofRapidNodes($seqLength, $infileSize);
	    
	    for (my $sbatchCounter=0;$sbatchCounter<$numberNodes-1;$sbatchCounter++) { #Parallization for each file handled
		
		&ProgramPreRequisites($sampleID, "BwaMem", $aligner, 0, $FILEHANDLE, $scriptParameter{'maximumCores'}, 5);	    
		
		my $readStart = $sbatchCounter *  $ReadNrofLines; #Constant for gz files
		my $readStop = $readStart + ceil( $ReadNrofLines + 1); #Constant for gz files	
		
		my $BWAinSampleDirectory = $indirpath{$sampleID};
		my $BWAoutSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/bwa"; 
		my $infile;

		if ($sequenceRunMode eq "Paired-end") { #Second read direction if present
	
		    $infile = $infile{$sampleID}[$infileCounter+$infileCounter]; #For required .fastq file
                }
                else { #Single-end
		    
		    $infile = $infile{$sampleID}[$infileCounter]; #For required .fastq file
                }
		
#BWA Mem	
		print $FILEHANDLE "bwa mem ";
		print $FILEHANDLE "-M "; #Mark shorter split hits as secondary (for Picard compatibility). 
		print $FILEHANDLE "-t ".$scriptParameter{'maximumCores'}." "; #Number of threads 
		print $FILEHANDLE "-R ".'"@RG\tID:'.$infilesLaneNoEnding{$sampleID}[$infileCounter].'\tSM:'.$sampleID.'\tPL:ILLUMINA" '; #read group header line
		print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #reference

		print $FILEHANDLE "<( "; #Pipe to BWA Mem (Read 1)
		print $FILEHANDLE "zcat "; #decompress Read 1
		print $FILEHANDLE $BWAinSampleDirectory."/".$infile." "; #Read 1
		print $FILEHANDLE "| "; #Pipe
		print $FILEHANDLE q?perl -ne 'if ( ($.>?.$readStart.q?) && ($.<?.$readStop.q?) ) {print $_;}' ?; #Limit to sbatch script interval
		print $FILEHANDLE ") "; #End Read 1

		if ($sequenceRunMode eq "Paired-end") { #Second read direction if present
		      
		    print $FILEHANDLE "<( "; #Pipe to BWA Mem (Read 2)
		    print $FILEHANDLE "zcat "; #decompress Read 2
		    print $FILEHANDLE $BWAinSampleDirectory."/".$infile{$sampleID}[$infileCounter+$infileCounter+1]." "; #Read 2
		    print $FILEHANDLE "| "; #Pipe
		    print $FILEHANDLE q?perl -ne 'if ( ($.>?.$readStart.q?) && ($.<?.$readStop.q?) ) {print $_;}' ?; #Limit to sbatch script interval
		    print $FILEHANDLE ") "; #End Read 2
		}

		print $FILEHANDLE "| "; #Pipe SAM to BAM conversion of aligned reads
		print $FILEHANDLE "samtools view "; 
		print $FILEHANDLE "-S "; #input is SAM
		print $FILEHANDLE "-h "; #print header for the SAM output
		print $FILEHANDLE "-u "; #uncompressed BAM output
		print $FILEHANDLE "- "; #/dev/stdin
		print $FILEHANDLE "| "; #Pipe
		print $FILEHANDLE "intersectBed "; #Limit output to only clinically interesting genes
		print $FILEHANDLE "-abam stdin "; #The A input file is in BAM format.  Output will be BAM as well.
		print $FILEHANDLE "-b ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'bwaMemRapidDb'}." "; #Db file of clinically relevant variants
		print $FILEHANDLE "> ".$BWAoutSampleDirectory."/".$infilesLaneNoEnding{$sampleID}[$infileCounter]."_".$sbatchCounter.".bam", "\n\n"; #Outfile (BAM)
		
		print $FILEHANDLE "samtools sort ";
		print $FILEHANDLE $BWAoutSampleDirectory."/".$infilesLaneNoEnding{$sampleID}[$infileCounter]."_".$sbatchCounter.".bam "; #Infile
		print $FILEHANDLE $BWAoutSampleDirectory."/".$infilesLaneNoEnding{$sampleID}[$infileCounter]."_".$sbatchCounter."_sorted", "\n\n"; #OutFile

		print $FILEHANDLE "samtools index ";
		print $FILEHANDLE $BWAoutSampleDirectory."/".$infilesLaneNoEnding{$sampleID}[$infileCounter]."_".$sbatchCounter."_sorted.bam", "\n\n"; #OutFile

		close($FILEHANDLE);
		if ( ($scriptParameter{'pBwaMem'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
		    &FIDSubmitJob($sampleID, $scriptParameter{'familyID'}, 3, $parameter{'pBwaMem'}{'chain'}, $fileName, $totalSbatchCounter);
		}
		$totalSbatchCounter++;
                #Save sbatch Counter to track how many read batch processes we have engaged
		$sampleInfo{$scriptParameter{'familyID'}}{$sampleID}{$infilesLaneNoEnding{$sampleID}[$infileCounter]}{'pBwaMem'}{'ReadBatchProcesses'} = $sbatchCounter+1;#Used to be  $sbatchCounter
		$sampleInfo{$scriptParameter{'familyID'}}{$sampleID}{'pBwaMem'}{'sbatchBatchProcesses'} = $totalSbatchCounter;
	    }
	}
	else { #Not rapid mode align whole file

	    &ProgramPreRequisites($sampleID, "BwaMem", $aligner, 0, $FILEHANDLE, $scriptParameter{'maximumCores'}, 5);
	    
	    my $BWAinSampleDirectory = $indirpath{$sampleID};
	    my $BWAoutSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/bwa"; 
	    
	    my $infile = $infile{$sampleID}[$pairedEndTracker]; #For required .fastq file
	    
	    print $FILEHANDLE "bwa mem ";
	    print $FILEHANDLE "-M "; #Mark shorter split hits as secondary (for Picard compatibility). 
	    print $FILEHANDLE "-t ".$scriptParameter{'maximumCores'}." "; #Number of threads 
	    print $FILEHANDLE "-R ".'"@RG\tID:'.$infilesLaneNoEnding{$sampleID}[$infileCounter].'\tSM:'.$sampleID.'\tPL:ILLUMINA" '; #read group header line
	    print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #reference
	    print $FILEHANDLE $BWAinSampleDirectory."/".$infile{$sampleID}[$pairedEndTracker]." "; #Read 1

	    if ($sequenceRunMode eq "Paired-end") { #Second read direction if present

		$pairedEndTracker = $pairedEndTracker+1; #Increment to collect correct read 2 from %infile 
		print $FILEHANDLE $BWAinSampleDirectory."/".$infile{$sampleID}[$pairedEndTracker]." "; #Read 2
	    }
	    $pairedEndTracker++;
	    print $FILEHANDLE "| "; #Pipe SAM to BAM conversion of aligned reads
	    print $FILEHANDLE "samtools view "; 
	    print $FILEHANDLE "-S "; #input is SAM
	    print $FILEHANDLE "-h "; #print header for the SAM output
	    print $FILEHANDLE "-u "; #uncompressed BAM output
	    print $FILEHANDLE "- "; #/dev/stdin
	    print $FILEHANDLE "> ".$BWAoutSampleDirectory."/".$infilesLaneNoEnding{$sampleID}[$infileCounter].".bam", "\n\n"; #Outfile (BAM)

	    close($FILEHANDLE);
	    if ( ($scriptParameter{'pBwaMem'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {

		$sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'} = $BWAoutSampleDirectory."/".$infile.".bam";
		&FIDSubmitJob($sampleID, $scriptParameter{'familyID'}, 3, $parameter{'pBwaMem'}{'chain'}, $fileName,  $infileCounter);
	    }
	}
    }
    return;
}

sub MosaikAlign {
###Aligning reads using MosaikAlign
    
    my $sampleID = $_[0];
    my $aligner = $_[1];
    
    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $sbatchScriptTracker=0;
    my $time=0;
    my $infileSize;

    for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #For all infiles per lane
	
	if ($infile{$sampleID}[$infileCounter] =~/.fastq.gz$/) { #Files are already gz and presently the scalar for compression has not been investigated. Therefore no automatic time allocation can be performed.
	    if ($scriptParameter{'analysisType'} eq "genomes") {
		$time = 80;  
	    }
	    else {
		$time = 40;
	    }
	}
	else { #Files are in fastq format
	
	    if (-e $indirpath{$sampleID}."/".$infile{$sampleID}[$infileCounter+$sbatchScriptTracker]) {
		$infileSize = -s $indirpath{$sampleID}."/".$infile{$sampleID}[$infileCounter+$sbatchScriptTracker]; # collect .fastq file size to enable estimation of time required for aligning, +$sbatchScriptTracker for syncing multiple infiles per sampleID. Hence, filesize will be calculated on read1 (should not matter).      
		$time = ceil(($infileSize/238)/(650*60*60)); #238 is a scalar estimating the number of reads depending on filesize. 650 is the number of reads/s in MosaikAlign-2.1.52 and 60*60 is to scale to hours.
	    }	    
	} 
	#Set parameters depending on sequence length
	my $seqLength = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesBothStrandsNoEnding{ $sampleID }[$infileCounter]}{'sequenceLength'};
	my $actParameter = 35; #the alignment candidate threshold (length)
	my $bwParameter = 35; #specifies the Smith-Waterman bandwidth.

	if ($seqLength <=36) {
	    
	    $actParameter = 20;
	    $bwParameter = 13;   
	}
	if ($seqLength >36 && $seqLength <=51) {
	    
	    $actParameter = 25;
	    $bwParameter = 21;   
	}
	if ($seqLength >51 && $seqLength <=76) {
	    
	    $bwParameter = 29;   
	}

	my ($stdoutPath) = &ProgramPreRequisites($sampleID, "MosaikAlign", $aligner, 0, $FILEHANDLE, $scriptParameter{'maximumCores'}, $time);
	my ($volume,$directories,$file) = File::Spec->splitpath($stdoutPath); #Split to enable submission to &SampleInfoQC later

	print $FILEHANDLE "mkdir -p /scratch/mosaik_tmp", "\n";
	print $FILEHANDLE "export MOSAIK_TMP=/scratch/mosaik_tmp", "\n\n";

	my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/mosaik";
	my $outSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/mosaik";
	my $infile = $infilesLaneNoEnding{$sampleID}[$infileCounter];

	print $FILEHANDLE "MosaikAligner ";
	print $FILEHANDLE "-in ".$inSampleDirectory."/".$infile.".dat "; #Infile
	print $FILEHANDLE "-out ".$outSampleDirectory."/".$infile." "; #OutFile (MosaikAligner appends .bam to infile name)
	print $FILEHANDLE "-ia ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'mosaikAlignReference'}." "; #Mosaik Reference
	print $FILEHANDLE "-annse ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'mosaikAlignNeuralNetworkSeFile'}." "; #NerualNetworkSE
	print $FILEHANDLE "-hs 15 "; #hash size
	print $FILEHANDLE "-mm 4 "; #the # of mismatches allowed
	print $FILEHANDLE "-mhp 100 "; #the maximum # of positions stored per seed

	if ($sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesLaneNoEnding{ $sampleID }[$infileCounter]}{'sequenceRunType'} eq "Paired-end") { #Second read direction if present
	    print $FILEHANDLE "-annpe ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'mosaikAlignNeuralNetworkPeFile'}." "; #NerualNetworkPE
	    print $FILEHANDLE "-ls 100 "; #enable local alignment search for PE reads
	}
	print $FILEHANDLE "-act ".$actParameter." "; #the alignment candidate threshold (length)
	print $FILEHANDLE "-bw ".$bwParameter." "; #specifies the Smith-Waterman bandwidth.
	print $FILEHANDLE "-j ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'mosaikJumpDbStub'}." "; #JumpDatabase
	print $FILEHANDLE "-p ".$scriptParameter{'maximumCores'}, "\n\n"; #Nr of cores
	
	print $FILEHANDLE "rm -rf /scratch/mosaik_tmp", "\n\n"; #Cleaning up temp directory

#BAM to SAM conversion 
	print $FILEHANDLE "java -Xmx4g ";
	print $FILEHANDLE "-jar ".$scriptParameter{'picardToolsPath'}."/SamFormatConverter.jar "; #Make sure that the BAM file BIN field is correct (Mosaik v.2.2.3 does according to Picard not set the bin field correctly)
	print $FILEHANDLE "VALIDATION_STRINGENCY=SILENT "; #Disable errors print 
	print $FILEHANDLE "INPUT=".$inSampleDirectory."/".$infile.".bam "; #InFile
	print $FILEHANDLE "OUTPUT=".$outSampleDirectory."/".$infile.".sam ", "\n\n"; #OutFile
	
	#SAM to BAM conversion 
	print $FILEHANDLE "java -Xmx4g ";
	print $FILEHANDLE "-jar ".$scriptParameter{'picardToolsPath'}."/SamFormatConverter.jar "; #Make sure that the BAM file BIN field is correct (Mosaik v.2.2.3 does according to Picard not set the bin field correctly)
	print $FILEHANDLE "INPUT=".$inSampleDirectory."/".$infile.".sam "; #InFile
	print $FILEHANDLE "OUTPUT=".$outSampleDirectory."/".$infile.".bam ", "\n\n"; #OutFile

	#Remove Sam
	print $FILEHANDLE "rm ".$outSampleDirectory."/".$infile.".sam ", "\n\n"; 

	close($FILEHANDLE);
	
	if ( ($scriptParameter{'pMosaikAlign'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
##Collect QC metadata info for later use                     	
	    $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'} = $outSampleDirectory."/".$infile.".bam";	
	    &SampleInfoQC($scriptParameter{'familyID'}, $sampleID, "MosaikAligner", $infile , $directories, $file, "infoDirectory"); #Outdata
	    &FIDSubmitJob($sampleID, $scriptParameter{'familyID'}, 3, $parameter{'pMosaikAlign'}{'chain'}, $fileName, $sbatchScriptTracker);
	}
	$sbatchScriptTracker++; #Tracks nr of sbatch scripts
    }
    return;
}

sub MosaikBuild {
#Generates Mosaik hash format on reads using MosaikBuild   
    
    my $sampleID = $_[0];
    my $aligner = $_[1];
    
    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $time = ceil(2.5*scalar( @{ $infilesLaneNoEnding{$sampleID} })); #One full lane on Hiseq takes approx. 1 h for MosaikBuild to process (compressed format, uncompressed 0.5 h), round up to nearest full hour.
    
    my $nrCores = &NrofCoresPerSbatch(scalar( @{$lane{$sampleID}} )); #Detect the number of cores to use from lanes
    
    &ProgramPreRequisites($sampleID, "MosaikBuild", $aligner, 0, $FILEHANDLE, $nrCores, $time);
    
    my $inSampleDirectory = $indirpath{$sampleID};
    my $outSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/mosaik";
    my $coreCounter=1;
    
    my $stParameter = "ILLUMINA"; #Default
    my  $pairedEndTracker = 0;
   
    for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #For all files
	
	my $sequenceRunMode = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesLaneNoEnding{ $sampleID }[$infileCounter]}{'sequenceRunType'}; #Collect paired-end or single-end sequence run mode
	my $coreTracker=0; #Required to portion out cores and files before wait and to track the outfiles to correct lane
	
	if ($infileCounter == $coreCounter*$nrCores) { #Using only nr of cores eq to lanes or maximumCores
	    
	    $coreCounter=$coreCounter+1;
	    print $FILEHANDLE "wait", "\n\n";
	}
	
	print $FILEHANDLE "MosaikBuild ";
	print $FILEHANDLE "-id ".$infilesLaneNoEnding{$sampleID}[$infileCounter]." "; #Read group ID for BAM Header
	print $FILEHANDLE "-sam ".$sampleID." "; #Sample name for BAM Header
	print $FILEHANDLE "-st ".$stParameter." "; #Sequencing technology for BAM Header
	print $FILEHANDLE "-mfl ".$scriptParameter{'mosaikBuildMedianFragLength'}." "; #Median Fragment Length
	print $FILEHANDLE "-q ".$inSampleDirectory."/".$infile{$sampleID}[$pairedEndTracker]." "; #Read 1
	
	if ( $sequenceRunMode eq "Paired-end") {
	    $pairedEndTracker = $pairedEndTracker+1; #Increment to collect correct read 2 from %infile
	    print $FILEHANDLE "-q2 ".$inSampleDirectory."/".$infile{$sampleID}[$pairedEndTracker]." "; #Read 2
	} 

	$pairedEndTracker++; #Increment to correctly track both seingle-end runs and paired-end runs
	print $FILEHANDLE "-out ".$outSampleDirectory."/".$infilesLaneNoEnding{$sampleID}[$infileCounter].".dat &", "\n\n"; #OutFile
    }
    print $FILEHANDLE "wait", "\n\n";    
    close($FILEHANDLE);
    if ( ($scriptParameter{'pMosaikBuild'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) { 
	&FIDSubmitJob($sampleID, $scriptParameter{'familyID'}, 1, $parameter{'pMosaikBuild'}{'chain'}, $fileName, 0);
    }
    return;
}   

sub FastQC {
#Raw sequence quality analysis using FASTQC

    my $sampleID = $_[0];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $time = ceil(0.5*scalar( @{ $infile{$sampleID} })); #One full lane on Hiseq takes approx. 0.5 h for FASTQC to process, round up to nearest full hour.

    my $nrCores = 0;

    for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #For all files   

	if ($sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesLaneNoEnding{ $sampleID }[$infileCounter]}{'sequenceRunType'} eq "Paired-end") { #Second read direction if present
	    $nrCores =  $nrCores + 2; #2 processes per file
	}
	else {#Single-end
	    $nrCores = $nrCores + 1; #Only 1 file and one process
	}
    }

    $nrCores = &NrofCoresPerSbatch($nrCores ); #Make sure that the number of cores does not exceed maximum after incrementing above

    &ProgramPreRequisites($sampleID, "FastQC", "fastqc", 0, $FILEHANDLE , $nrCores, $time);
    
    my $inSampleDirectory = $indirpath{$sampleID};
    my $outSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/fastqc";
    my $coreCounter=1;
    
    for (my $infileCounter=0;$infileCounter<scalar( @{ $infile{$sampleID} });$infileCounter++) {
	
	if ($infileCounter == $coreCounter*$nrCores) { #Using only nr of cores eq to lanes or maximumCores
	    
	    print $FILEHANDLE "wait", "\n\n";
	    $coreCounter=$coreCounter+1;
	}

	my $infile = $infile{$sampleID}[$infileCounter];

	print $FILEHANDLE "fastqc ";
	print $FILEHANDLE $inSampleDirectory."/".$infile." "; #InFile
	print $FILEHANDLE "-o ".$outSampleDirectory. " &", "\n\n"; #OutFile

##Collect QC metadata info for active program for later use
	if ( ($scriptParameter{'pFastQC'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	    &SampleInfoQC($scriptParameter{'familyID'}, $sampleID, "FastQC", $infile, $outSampleDirectory."/".$sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesBothStrandsNoEnding{ $sampleID }[$infileCounter]}{'originalFileNameNoEnding'}."_fastqc", "fastqc_data.txt", "static");
	}
    }
    print $FILEHANDLE "wait", "\n";    
    
    close($FILEHANDLE);
    if ( ($scriptParameter{'pFastQC'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	&FIDSubmitJob($sampleID, $scriptParameter{'familyID'}, 2, $parameter{'pFastQC'}{'chain'}, $fileName, 0);
    }
    return;
}

sub GZipFastq { 
#Automatically gzips fastq files. 
    
    my $sampleID = $_[0];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $time = ceil(1.5*scalar( @{ $infile{$sampleID} })); #One full lane on Hiseq takes approx. 1.5 h for gzip to process, round up to nearest full hour.

    &ProgramPreRequisites($sampleID, "GZip", "gzip", 0, $FILEHANDLE, $scriptParameter{'maximumCores'}, $time);
   
    print $FILEHANDLE "cd ".$indirpath{$sampleID}, "\n\n";
   
    my $inSampleDirectory = $indirpath{$sampleID};
    my $coreCounter=1;
    my $uncompressedFileCounter = 0; #Used to print wait at the right times since infiles cannot be used (can be a mixture of .gz and .fast files)
    
    for (my $infileCounter=0;$infileCounter<scalar( @{ $infile{$sampleID} });$infileCounter++) {

	if ($infile{$sampleID}[$infileCounter] =~/.fastq$/) { #For files ending with .fastq required since there can be a mixture (also .fastq.gz) within the sample dir
	    if ($uncompressedFileCounter == $coreCounter*$scriptParameter{'maximumCores'}) { #Using only $scriptParameter{'maximumCores'} cores
		
		print $FILEHANDLE "wait", "\n\n";
		$coreCounter=$coreCounter+1;
	    }
	    my $infile = $infile{$sampleID}[$infileCounter];
	    print $FILEHANDLE "gzip ";
	    print $FILEHANDLE $inSampleDirectory."/".$infile," &", "\n\n"; #InFile
	    $uncompressedFileCounter++;
	    $infile{$sampleID}[$infileCounter] =~ s/.fastq/.fastq.gz/g; #Replace the .fastq ending with .fastq.gz since this will execute before fastQC screen and mosaikBuild, hence changing the original file name ending from ".fastq" to ".fastq.gz". 

	}
    }
    print $FILEHANDLE "wait", "\n\n";
    if ( ($scriptParameter{'pGZip'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	&FIDSubmitJob($sampleID, $scriptParameter{'familyID'}, 0, $parameter{'pGZip'}{'chain'}, $fileName, 0);
    }
    return;
}

sub BuildAnnovarPreRequisites {
##Creates the AnnovarPreRequisites

    my $familyID = $_[0];
    my $aligner = $_[1];
    my $program = $_[2];
    
    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    $parameter{'annovarBuildReference'}{'buildFile'} = 0; #Ensure that this subrutine is only executed once
    my $annovarTemporaryDirectory = $scriptParameter{'annovarPath'}."/humandb/Db_temporary"; #Temporary download directory
    
    &ProgramPreRequisites($familyID, $program, $aligner, 0, $FILEHANDLE, 1, 3);

    print STDOUT "\nNOTE: Will try to create required Annovar database files before executing ".$program,"\n\n";

    print $FILEHANDLE "#Make temporary download directory\n\n"; 
    print $FILEHANDLE "mkdir -p ".$annovarTemporaryDirectory."; ", "\n\n"; 

    print $FILEHANDLE "#Downloading Annovar Db files", "\n\n";

    for (my $tableNamesCounter=0;$tableNamesCounter<scalar(@annovarTableNames);$tableNamesCounter++) { #For all specified table names
	
	if ($parameter{$annovarTableNames[$tableNamesCounter]}{'buildFile'} eq 1) {
	    
	    print $FILEHANDLE "perl ".$scriptParameter{'annovarPath'}."/annotate_variation.pl "; #Annovar script 
	    print $FILEHANDLE "-buildver ".$scriptParameter{'annovarGenomeBuildVersion'}." "; #GenomeBuild version
	    print $FILEHANDLE "-downdb ".$annovarTables{$annovarTableNames[$tableNamesCounter]}{'download'}." "; #Db to download
	    
	    if (defined($annovarTables{$annovarTableNames[$tableNamesCounter]}{'ucscAlias'})) {
		
		print $FILEHANDLE "-webfrom ucsc "; #Download from ucsc
	    }
	    else {
		
		print $FILEHANDLE "-webfrom annovar "; #Download from annovar
	    }
	    print $FILEHANDLE $annovarTemporaryDirectory."/ ", "\n\n"; #annovar/humandb directory is assumed
	    
##Check file existance and move created file if lacking 
	    my $intendedFilePathRef;
	    my $temporaryFilePathRef;
	    
	    if (defined($annovarTables{ $annovarTableNames[$tableNamesCounter] }{'file'})) {
		
		for (my $filesCounter=0;$filesCounter<scalar(@{$annovarTables{ $annovarTableNames[$tableNamesCounter] }{'file'}});$filesCounter++) { #All annovarTables file(s), some tables have multiple files downloaded from the same call
		    
		    $intendedFilePathRef = \($scriptParameter{'annovarPath'}."/humandb/".$annovarTables{ $annovarTableNames[$tableNamesCounter] }{'file'}[$filesCounter]);  
		    $temporaryFilePathRef = \($annovarTemporaryDirectory."/".$annovarTables{ $annovarTableNames[$tableNamesCounter] }{'file'}[$filesCounter]);    
		    &PrintCheckExistandMoveFile($FILEHANDLE, $intendedFilePathRef, $temporaryFilePathRef);
		
		    if (defined($annovarTables{ $annovarTableNames[$tableNamesCounter] }{'indexFile'})) {
		    
			$intendedFilePathRef = \($scriptParameter{'annovarPath'}."/humandb/".$annovarTables{ $annovarTableNames[$tableNamesCounter] }{'file'}[$filesCounter].".idx");  
			$temporaryFilePathRef = \($annovarTemporaryDirectory."/".$annovarTables{ $annovarTableNames[$tableNamesCounter] }{'file'}[$filesCounter].".idx");
			&PrintCheckExistandMoveFile($FILEHANDLE, $intendedFilePathRef, $temporaryFilePathRef);	
		    }
		}		
	    }
	    elsif ((defined($annovarTables{ $annovarTableNames[$tableNamesCounter] }{'ucscAlias'}))){
	    
		$intendedFilePathRef = \($scriptParameter{'annovarPath'}."/humandb/".$scriptParameter{'annovarGenomeBuildVersion'}."_".$annovarTables{ $annovarTableNames[$tableNamesCounter] }{'ucscAlias'}.".txt");
		$temporaryFilePathRef = \($annovarTemporaryDirectory."/".$scriptParameter{'annovarGenomeBuildVersion'}."_".$annovarTables{ $annovarTableNames[$tableNamesCounter] }{'ucscAlias'}.".txt");
		&PrintCheckExistandMoveFile($FILEHANDLE, $intendedFilePathRef, $temporaryFilePathRef);
	    
		if (defined($annovarTables{ $annovarTableNames[$tableNamesCounter] }{'indexFile'})) {
		    
		    $intendedFilePathRef = \($scriptParameter{'annovarPath'}."/humandb/".$scriptParameter{'annovarGenomeBuildVersion'}."_".$annovarTables{ $annovarTableNames[$tableNamesCounter] }{'ucscAlias'}.".txt.idx");
		    $temporaryFilePathRef = \($annovarTemporaryDirectory."/".$scriptParameter{'annovarGenomeBuildVersion'}."_".$annovarTables{ $annovarTableNames[$tableNamesCounter] }{'ucscAlias'}.".txt.idx");
		    &PrintCheckExistandMoveFile($FILEHANDLE, $intendedFilePathRef, $temporaryFilePathRef);
		}
	    }
	    else {
	    
		$intendedFilePathRef = \($scriptParameter{'annovarPath'}."/humandb/".$scriptParameter{'annovarGenomeBuildVersion'}."_".$annovarTableNames[$tableNamesCounter].".txt");
		$temporaryFilePathRef = \($annovarTemporaryDirectory."/".$scriptParameter{'annovarGenomeBuildVersion'}."_".$annovarTableNames[$tableNamesCounter].".txt");    
		&PrintCheckExistandMoveFile($FILEHANDLE, $intendedFilePathRef, $temporaryFilePathRef);
	    
		if (defined($annovarTables{ $annovarTableNames[$tableNamesCounter] }{'indexFile'})) {
	
		    $intendedFilePathRef = \($scriptParameter{'annovarPath'}."/humandb/".$scriptParameter{'annovarGenomeBuildVersion'}."_".$annovarTableNames[$tableNamesCounter].".txt.idx");
		    $temporaryFilePathRef = \($annovarTemporaryDirectory."/".$scriptParameter{'annovarGenomeBuildVersion'}."_".$annovarTableNames[$tableNamesCounter].".txt.idx");    
		    &PrintCheckExistandMoveFile($FILEHANDLE, $intendedFilePathRef, $temporaryFilePathRef);	
		}				
	    }
	}
        $parameter{$annovarTableNames[$tableNamesCounter]}{'buildFile'} = 0;
    }
    
    print $FILEHANDLE "rm -rf $annovarTemporaryDirectory;", "\n\n"; #Cleaning up temp directory
    close($FILEHANDLE);
    
    if ( ($scriptParameter{"p".$program} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	&FIDSubmitJob(0, $familyID, 6, $parameter{"p".$program}{'chain'}, $fileName, 0);
    }
}

sub BuildPTCHSMetricPreRequisites {
##Creates the target "infiles_list" "padded.infile_list" and interval_list files

    my $familyID = $_[0];
    my $aligner = $_[1];
    my $program = $_[2];   
    
    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $parametersToEvaluate = 0; #The number of parameters to evaluate

    &ProgramPreRequisites($familyID, $program, $aligner, 0, $FILEHANDLE, 1, 1);

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) { #All sampleIDs

	my $sampleIDBuildSwitchInfile = $parameter{ $scriptParameter{'familyID'} }{$sampleIDs[$sampleIDCounter]}{'exomeTargetBedInfileLists'}{'buildFile'};
	my $sampleIDBuildFileInfile = $scriptParameter{ $scriptParameter{'familyID'} }{$sampleIDs[$sampleIDCounter]}{'exomeTargetBedInfileLists'};
	my $infileListNoEnding = &RemoveFileEnding(\$sampleIDBuildFileInfile, $referenceFileEndings{'exomeTargetBedInfileLists'}); #For comparison of identical filename.bed files, to avoid creating files twice

	my $sampleIDBuildSwitchPadded = $parameter{ $scriptParameter{'familyID'} }{$sampleIDs[$sampleIDCounter]}{'exomeTargetPaddedBedInfileLists'}{'buildFile'};
	my $sampleIDBuildFilePadded = $scriptParameter{ $scriptParameter{'familyID'} }{$sampleIDs[$sampleIDCounter]}{'exomeTargetPaddedBedInfileLists'};
	my $paddedInfileListNoEnding = &RemoveFileEnding(\$sampleIDBuildFilePadded , $referenceFileEndings{'exomeTargetPaddedBedInfileLists'}); #For comparison of identical filename.bed files, to avoid creating files twice

	my $sampleIDBuildSwitchPaddedInterval = $parameter{ $scriptParameter{'familyID'} }{$sampleIDs[$sampleIDCounter]}{'GATKTargetPaddedBedIntervalLists'}{'buildFile'};
	my $sampleIDBuildFilePaddedInterval = $scriptParameter{ $scriptParameter{'familyID'} }{$sampleIDs[$sampleIDCounter]}{'GATKTargetPaddedBedIntervalLists'};
	my $paddedIntervalListNoEnding = &RemoveFileEnding(\$sampleIDBuildFilePaddedInterval , $referenceFileEndings{'GATKTargetPaddedBedIntervalLists'}); #For comparison of identical filename.bed files, to avoid creating files twice

	if ( (defined($sampleIDBuildSwitchPadded)) && ($sampleIDBuildSwitchPadded == 1) ) { #If identical filename.bed and padded files need creation do not build infile_list in separate part of sbatch

	    if ($infileListNoEnding eq $paddedInfileListNoEnding) {		
	
		$sampleIDBuildSwitchInfile = 0; #Turn of separate infile_list creation
	    }	
	    $parametersToEvaluate = $sampleIDBuildSwitchPadded; #Add to parameters to evaluate (1 or 0)
	}
	if ( (defined($sampleIDBuildSwitchPaddedInterval)) && ($sampleIDBuildSwitchPaddedInterval == 1) ) { #If identical filename.bed and padded files need creation do not build .pad100.infile_list in separate part of sbatch

	    if ($paddedInfileListNoEnding eq $paddedIntervalListNoEnding) {
		
		$sampleIDBuildSwitchPadded = 0; #Turn of seperate paddded infile_list creation
	    }	
	    $parametersToEvaluate = $parametersToEvaluate + $sampleIDBuildSwitchPaddedInterval;  #Add to parameters to evaluate (1 or 0)
	} 
	if (defined($sampleIDBuildSwitchInfile)) {

	    $parametersToEvaluate = $parametersToEvaluate + $sampleIDBuildSwitchInfile;  #Add to parameters to evaluate (1 or 0)
	}
	
        #Turn of build of identical filename.bed files
	&CheckUniqueTargetFiles(\@sampleIDs, \$sampleIDCounter, \$sampleIDBuildFileInfile, "exomeTargetBedInfileLists");
	&CheckUniqueTargetFiles(\@sampleIDs, \$sampleIDCounter, \$sampleIDBuildFilePadded, "exomeTargetPaddedBedInfileLists");
	&CheckUniqueTargetFiles(\@sampleIDs, \$sampleIDCounter, \$sampleIDBuildFilePaddedInterval, "GATKTargetPaddedBedIntervalLists");
	
	for (my $parameterCounter=0;$parameterCounter<$parametersToEvaluate;$parameterCounter++) {

	    my $randomInteger = int(rand(10000)); #Generate a random integer between 0-10,000.		
	    
##Initiate general build variables used for all parameters
	    my $sampleIDBuildFile;
	    my $sampleIDBuildFileNoEnding;
	    my $sampleIDBuildFileNoEndingTemp;
	    
	    if ( (defined($sampleIDBuildSwitchInfile)) && ($sampleIDBuildSwitchInfile eq 1) ) {
		
		&SetTargetFileGeneralBuildParameter(\$sampleIDBuildFileInfile, "exomeTargetBedInfileLists", \$sampleIDBuildFile, \$sampleIDBuildFileNoEnding, \$sampleIDs[$sampleIDCounter]);
	    }
	    elsif ( (defined($sampleIDBuildSwitchPadded)) && ($sampleIDBuildSwitchPadded eq 1) ) {

		&SetTargetFileGeneralBuildParameter(\$sampleIDBuildFilePadded, "exomeTargetPaddedBedInfileLists", \$sampleIDBuildFile, \$sampleIDBuildFileNoEnding, \$sampleIDs[$sampleIDCounter]);
	    }
	    elsif ( (defined($sampleIDBuildSwitchPaddedInterval)) && ($sampleIDBuildSwitchPaddedInterval == 1) ) {
		
		&SetTargetFileGeneralBuildParameter(\$sampleIDBuildFilePaddedInterval, "GATKTargetPaddedBedIntervalLists", \$sampleIDBuildFile, \$sampleIDBuildFileNoEnding, \$sampleIDs[$sampleIDCounter]);
	    }
	    
	    if (defined($sampleIDBuildFile)) {
		
		$sampleIDBuildFileNoEndingTemp = $sampleIDBuildFileNoEnding."_".$randomInteger; #Add random integer	
		
		print STDOUT "\nNOTE: Will try to create required ".$sampleIDBuildFile." file before executing ".$program,"\n\n";
		
		print $FILEHANDLE "SampleID:".$sampleIDs[$sampleIDCounter], "\n\n";
		print $FILEHANDLE "#CreateSequenceDictionary from reference", "\n";
		print $FILEHANDLE "java -Xmx2g -jar ".$scriptParameter{'picardToolsPath'}."/CreateSequenceDictionary.jar ";
		print $FILEHANDLE "R=".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference genome
		print $FILEHANDLE "OUTPUT=".$scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict ", "\n\n"; #Output sequence dictionnary
		
		print $FILEHANDLE "#Add target file to headers from sequenceDictionary", "\n";
		print $FILEHANDLE "cat "; #Concatenate
		print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict "; #sequence dictionnary
		print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEnding." "; #Bed file
		print $FILEHANDLE "> "; #Write to
		print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict_body", "\n\n"; #Add bed body to dictionnary
		
		print $FILEHANDLE "#Remove target annotations, 'track', 'browse' and keep only 5 columns", "\n";
		print $FILEHANDLE q?perl  -nae 'if ($_=~/@/) {print $_;} elsif ($_=~/^track/) {} elsif ($_=~/^browser/) {} else {print @F[0], "\t", (@F[1] + 1), "\t", @F[2], "\t", "+", "\t", "-", "\n";}' ?;
		print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict_body "; #Infile
		print $FILEHANDLE "> "; #Write to
		print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict_body_col_5", "\n\n"; #Remove unnecessary info and reformat 
		
		print $FILEHANDLE "Create".$referenceFileEndings{'exomeTargetBedInfileLists'}, "\n";
		print $FILEHANDLE "java -Xmx2g -jar ".$scriptParameter{'picardToolsPath'}."/IntervalListTools.jar ";
		print $FILEHANDLE "INPUT=".$scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict_body_col_5 ";
		print $FILEHANDLE "OUTPUT=".$scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict_body_col_5_".$referenceFileEndings{'exomeTargetBedInfileLists'}." ", "\n\n";
		    
		print $FILEHANDLE "[ -s ".$scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEnding.$referenceFileEndings{'exomeTargetBedInfileLists'}." ] "; #Check file exists and is larger than 0
		print $FILEHANDLE "&& rm ".$scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict_body_col_5_".$referenceFileEndings{'exomeTargetBedInfileLists'}." "; #If other processes already has created file, remove temp file
		print $FILEHANDLE "|| "; #File has not been created by other processes
		print $FILEHANDLE "mv ".$scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict_body_col_5_".$referenceFileEndings{'exomeTargetBedInfileLists'}." ".$scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEnding.$referenceFileEndings{'exomeTargetBedInfileLists'}, "\n\n"; #Move file in place
	    }
	    if ( (defined($sampleIDBuildSwitchPadded) && ($sampleIDBuildSwitchPadded eq 1)) || (defined($sampleIDBuildSwitchPaddedInterval) && ($sampleIDBuildSwitchPaddedInterval eq 1)) ) {
		
		print $FILEHANDLE "Create padded interval list", "\n";
		print $FILEHANDLE "java -Xmx2g -jar ".$scriptParameter{'picardToolsPath'}."/IntervalListTools.jar ";
		print $FILEHANDLE "PADDING=100 "; #Add 100 nt on both sides of bed entry
		print $FILEHANDLE "INPUT=".$scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict_body_col_5 ";
		print $FILEHANDLE "OUTPUT=".$scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict_body_col_5".$referenceFileEndings{'exomeTargetPaddedBedInfileLists'}." ", "\n\n";
		
		print $FILEHANDLE "[ -s ".$scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEnding.$referenceFileEndings{'exomeTargetPaddedBedInfileLists'}." ] "; #Check file exists and is larger than 0
		print $FILEHANDLE "&& rm ".$scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict_body_col_5".$referenceFileEndings{'exomeTargetPaddedBedInfileLists'}." "; #If other processes already has created file, remove temp file
		print $FILEHANDLE "|| "; #File has not been created by other processes
		print $FILEHANDLE "(mv ".$scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict_body_col_5".$referenceFileEndings{'exomeTargetPaddedBedInfileLists'}." ".$scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEnding.$referenceFileEndings{'exomeTargetPaddedBedInfileLists'}."; ";
		
		if (defined($sampleIDBuildSwitchPaddedInterval) && ($sampleIDBuildSwitchPaddedInterval eq 1)) {
		    
		    ##Softlink '.interval_list' to padded .infile_list", "\n";
		    print $FILEHANDLE "ln -s "; #Softlink
		    print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEnding.$referenceFileEndings{'exomeTargetPaddedBedInfileLists'}." "; #Origin file
		    print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEnding.$referenceFileEndings{'GATKTargetPaddedBedIntervalLists'}; #interval_list file
		}
		
		print $FILEHANDLE ")", "\n\n"; #Move file in place
	    }
	    if (defined($sampleIDBuildFile)) {
		
		print $FILEHANDLE "#Remove temporary files", "\n";
		
		print $FILEHANDLE "rm ";
		print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict_body_col_5 ", "\n\n";
		
		print $FILEHANDLE "rm ";
		print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict_body ", "\n\n";
		
		print $FILEHANDLE "rm ";
		print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict ", "\n\n";
		
		if ( (defined($sampleIDBuildSwitchPadded)) && ($sampleIDBuildSwitchPadded == 0) ) {
		    
		    $sampleIDBuildSwitchPaddedInterval = 0;
		}
		if ( (defined($sampleIDBuildSwitchInfile)) && ($sampleIDBuildSwitchInfile == 0) ) {
		    
		    $sampleIDBuildSwitchPadded = 0;
		}
		$sampleIDBuildSwitchInfile = 0;
	    }
	}
    }
    close($FILEHANDLE);
    
    if ( ($scriptParameter{"p".$program} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	&FIDSubmitJob(0, $familyID, 6, "MIP", $fileName, 0); #"MIP" is required or the pPicardToolsCalulateHSMetrics jobs will start prematurely
    }
}

sub BuildBwaPreRequisites {
##Creates the BwaPreRequisites using scriptParameters{'humanGenomeReference'} as reference.

    my $familyID = $_[0];
    my $aligner = $_[1];
    my $program = $_[2];
    
    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $randomInteger = int(rand(10000)); #Generate a random integer between 0-10,000.
    
    $parameter{'bwaBuildReference'}{'buildFile'} = 0; #Ensure that this subrutine is only executed once

    &ProgramPreRequisites($familyID, $program, $aligner, 0, $FILEHANDLE, 1, 3);

    print STDOUT "\nNOTE: Will try to create required ".$scriptParameter{'bwaBuildReference'}." index files before executing ".$program,"\n\n";

    print $FILEHANDLE "#Building BWA index", "\n\n";
    print $FILEHANDLE "bwa index "; #index sequences in the FASTA format
    print $FILEHANDLE "-p ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'bwaBuildReference'}."_".$randomInteger." "; #prefix of the index
    print $FILEHANDLE "-a bwtsw "; #BWT construction algorithm
    print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'},"\n\n"; #the FASTA reference sequences file

    for (my $fileEndingsCounter=0;$fileEndingsCounter<scalar(@bwaBuildReferenceFileEndings);$fileEndingsCounter++) { #All fileEndings

	print $FILEHANDLE "[ -s ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'bwaBuildReference'}.$bwaBuildReferenceFileEndings[$fileEndingsCounter]." ] "; #Check file exists and is larger than 0
	print $FILEHANDLE "&& rm ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'bwaBuildReference'}."_".$randomInteger.$bwaBuildReferenceFileEndings[$fileEndingsCounter]." "; #If other processes already has created file, remove temp file
	print $FILEHANDLE "|| "; #File has not been created by other processes
	print $FILEHANDLE "mv ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'bwaBuildReference'}."_".$randomInteger.$bwaBuildReferenceFileEndings[$fileEndingsCounter]." ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'bwaBuildReference'}.$bwaBuildReferenceFileEndings[$fileEndingsCounter], "\n\n"; #Move file in place
    }

    close($FILEHANDLE);
    
    if ( ($scriptParameter{"p".$program} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	&FIDSubmitJob(0, $familyID, 6, $parameter{"p".$program}{'chain'}, $fileName, 0);
    }
}

sub BuildMosaikAlignPreRequisites {
##Creates the mosaikAlignPreRequisites using scriptParameters{'humanGenomeReference'} as reference.

    my $familyID = $_[0];
    my $aligner = $_[1];
    my $program = $_[2];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $randomInteger = int(rand(10000)); #Generate a random integer between 0-10,000.

    &ProgramPreRequisites($familyID, $program, $aligner, 0, $FILEHANDLE, 4, 2);
    
    if ($parameter{'mosaikAlignReference'}{'buildFile'} eq 1) {

	print STDOUT "\nNOTE: Will try to create required ".$scriptParameter{'mosaikAlignReference'}." before executing ".$program,"\n\n";

	print $FILEHANDLE "#Building MosaikAligner Reference", "\n\n";
	print $FILEHANDLE "MosaikBuild ";
	print $FILEHANDLE "-fr ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #the FASTA reference sequences file
	print $FILEHANDLE "-sn Homo_sapiens "; #Species name
	print $FILEHANDLE "-ga ".$humanGenomeReferenceSource.$humanGenomeReferenceVersion." "; #the genome assembly ID
	print $FILEHANDLE "-oa ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'mosaikAlignReference'}."_".$randomInteger, "\n\n";

	print $FILEHANDLE "[ -s ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'mosaikAlignReference'}." ] "; #Check file exists and is larger than 0
	print $FILEHANDLE "&& rm ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'mosaikAlignReference'}."_".$randomInteger." "; #If other processes already has created file, remove temp file
	print $FILEHANDLE "|| "; #File has not been created by other processes
	print $FILEHANDLE "mv ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'mosaikAlignReference'}."_".$randomInteger." ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'mosaikAlignReference'}, "\n\n"; #Move file in place

    }
    if ($parameter{'mosaikJumpDbStub'}{'buildFile'} eq 1) {

	print STDOUT "\nNOTE: Will try to create required ".$scriptParameter{'mosaikJumpDbStub'}." before executing ".$program,"\n\n";

	print $FILEHANDLE "#Building MosaikAligner JumpDatabase", "\n\n";
	
	print $FILEHANDLE "mkdir -p /scratch/mosaik_tmp", "\n";
	print $FILEHANDLE "export MOSAIK_TMP=/scratch/mosaik_tmp", "\n\n";
	
	print $FILEHANDLE "MosaikJump ";
	print $FILEHANDLE "-ia ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'mosaikAlignReference'}." "; #The input reference file  
	print $FILEHANDLE "-hs 15 "; #the hash size
	print $FILEHANDLE "-mem 24 "; #the amount memory used when sorting hashes
	print $FILEHANDLE "-out ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'mosaikJumpDbStub'}."_".$randomInteger, "\n\n"; #Mosaik JumpDbStub for the output filenames
	
##Meta
	print $FILEHANDLE "[ -s ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'mosaikJumpDbStub'}."_meta.jmp ] "; #Check file exists and is larger than 0
	print $FILEHANDLE "&& rm ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'mosaikJumpDbStub'}."_".$randomInteger."_meta.jmp "; #If other processes already has created file, remove temp file
	print $FILEHANDLE "|| "; #File has not been created by other processes
	print $FILEHANDLE "mv ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'mosaikJumpDbStub'}."_".$randomInteger."_meta.jmp ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'mosaikJumpDbStub'}."_meta.jmp ", "\n\n"; #Move file in place
	
##Keys
	print $FILEHANDLE "[ -s ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'mosaikJumpDbStub'}."_keys.jmp ] "; #Check file exists and is larger than 0
	print $FILEHANDLE "&& rm ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'mosaikJumpDbStub'}."_".$randomInteger."_keys.jmp "; #If other processes already has created file, remove temp file
	print $FILEHANDLE "|| "; #File has not been created by other processes
	print $FILEHANDLE "mv ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'mosaikJumpDbStub'}."_".$randomInteger."_keys.jmp ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'mosaikJumpDbStub'}."_keys.jmp ", "\n\n"; #Move file in place
	
##Positions
	print $FILEHANDLE "[ -s ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'mosaikJumpDbStub'}."_positions.jmp ] "; #Check file exists and is larger than 0
	print $FILEHANDLE "&& rm ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'mosaikJumpDbStub'}."_".$randomInteger."_positions.jmp "; #If other processes already has created file, remove temp file
	print $FILEHANDLE "|| "; #File has not been created by other processes
	print $FILEHANDLE "mv ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'mosaikJumpDbStub'}."_".$randomInteger."_positions.jmp ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'mosaikJumpDbStub'}."_positions.jmp ", "\n\n"; #Move file in place
	
	print $FILEHANDLE "rm -rf /scratch/mosaik_tmp", "\n\n"; #Cleaning up temp directory
    }
    close($FILEHANDLE);
    
    if ( ($scriptParameter{"p".$program} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	&FIDSubmitJob(0, $familyID, 6, $parameter{"p".$program}{'chain'}, $fileName, 0);
    }
}

sub CheckBuildHumanGenomePreRequisites {

    for (my $fileEndingsCounter=0;$fileEndingsCounter<scalar(@humanGenomeReferenceFileEndings);$fileEndingsCounter++) {
	
	if ( ($parameter{"humanGenomeReference".$humanGenomeReferenceFileEndings[$fileEndingsCounter]}{'buildFile'} eq 1) || ($humanGenomeCompressed eq "compressed") ) {
	    
	    &BuildHumanGenomePreRequisites($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "PicardToolsCalculateHSMetrics");
	    last;#Will handle all meatfiles build within sbatch script
	}
    }
    ##Collect sequence contigs from human reference
    #&CollectSeqContigs(); #Reloads if required NOTE:Preparation for future changes but not activated yet
}

sub BuildHumanGenomePreRequisites {
##Creates the humanGenomePreRequisites using scriptParameters{'humanGenomeReference'} as reference.

    my $familyID = $_[0];
    my $aligner = $_[1];
    my $program = $_[2];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $randomInteger = int(rand(10000)); #Generate a random integer between 0-10,000.

    &ProgramPreRequisites($familyID, $program, $aligner, 0, $FILEHANDLE, 1, 1);

    if ($humanGenomeCompressed eq "compressed") {

	print STDOUT "\nNOTE: Will try to decompres ".$scriptParameter{'humanGenomeReference'}." before executing ".$program,"\n\n";

	print $FILEHANDLE "gzip ";
	print $FILEHANDLE "-d "; #Decompress
	print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}, "\n\n";
	$scriptParameter{'humanGenomeReference'} =~ s/.fasta.gz/.fasta/g; #Replace the .fasta.gz ending with .fasta since this will execute before the analysis, hence changing the original file name ending from ".fastq" to ".fastq.gz".
	print STDOUT "Set humanGenomeReference to: ".$scriptParameter{'humanGenomeReference'}, "\n\n";
	$humanGenomeCompressed = "unCompressed";
    }

    

    for (my $fileEndingsCounter=0;$fileEndingsCounter<scalar(@humanGenomeReferenceFileEndings);$fileEndingsCounter++) { #All meta files    
	
	if ($parameter{"humanGenomeReference.dict"}{'buildFile'} eq 1) { #.dict file

	    print STDOUT "\nNOTE: Will try to create dict file for ".$scriptParameter{'humanGenomeReference'}." before executing ".$program,"\n\n";
	    
	    print $FILEHANDLE "#CreateSequenceDictionary from reference", "\n";
	    print $FILEHANDLE "java -Xmx2g -jar ".$scriptParameter{'picardToolsPath'}."/CreateSequenceDictionary.jar ";
	    print $FILEHANDLE "R=".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference genome
	    print $FILEHANDLE "OUTPUT=".$scriptParameter{'referencesDir'}."/".$humanGenomeReferenceNameNoEnding."_".$randomInteger.".dict ", "\n\n"; #Output sequence dictionnary
	    
	    &PrintCheckExistandMoveFile($FILEHANDLE, \($scriptParameter{'referencesDir'}."/".$humanGenomeReferenceNameNoEnding.".dict"), \($scriptParameter{'referencesDir'}."/".$humanGenomeReferenceNameNoEnding."_".$randomInteger.".dict"));
	    
	    $parameter{"humanGenomeReference.dict"}{'buildFile'} = 0; #Only create once

	}
	if ($parameter{"humanGenomeReference.fasta.fai"}{'buildFile'} eq 1) {

	    print STDOUT "\nNOTE: Will try to create .fai file for ".$scriptParameter{'humanGenomeReference'}." before executing ".$program,"\n\n";

	    print $FILEHANDLE "#Fai file from reference", "\n";
	    print $FILEHANDLE "ln -s "; #Softlink
	    print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference genome
	    print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}."_".$randomInteger, "\n\n"; #Softlink to Reference genome
	    
	    print $FILEHANDLE "samtools faidx ";#index/extract FASTA
	    print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}."_".$randomInteger, "\n\n"; #Softlink to Reference genome
	    
	    &PrintCheckExistandMoveFile($FILEHANDLE, \($scriptParameter{'referencesDir'}."/".$humanGenomeReferenceNameNoEnding.".fasta.fai"), \($scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}."_".$randomInteger.".fai"));
	
	    print $FILEHANDLE "rm "; #Remove softLink
	    print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}."_".$randomInteger, "\n\n"; #Softlink to Reference genome
	    
	    $parameter{"humanGenomeReference.fasta.fai"}{'buildFile'} = 0; #Only create once	
	}
    }
    close($FILEHANDLE);
    
    if ( ($scriptParameter{"p".$program} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	&FIDSubmitJob(0, $familyID, 6, "MIP", $fileName, 0);
    }
}


sub ReadPlinkPedigreeFile {
###Reads famid_pedigree.txt file in PLINK format
###FORMAT: FamliyID\tSampleID\tFather\tMother\tSex(1=male; 2=female; other=unknown)\tPhenotype(-9 missing, 0 missing, 1 unaffected, 2 affected)..n

    my $fileName = $_[0];      
    
    my @pedigreeFileElements = ("FamilyID", "SampleID", "Father", "Mother", "Sex", "Phenotype", );
    my $familyID;
    my $sampleID;
    
    my $userSampleIDsSwitch = &CheckUserInfoArrays(\@sampleIDs, "sampleIDs");
    my $userExomeTargetBedInfileListsSwitch = &CheckUserInfoArrays(\@exomeTargetBedInfileLists, "exomeTargetBedInfileLists"); 
    my $userExomeTargetPaddedBedInfileListSwitch = &CheckUserInfoArrays(\@exomeTargetPaddedBedInfileLists, "exomeTargetPaddedBedInfileLists");
    my $userExomeTargetPaddedBedIntervalListSwitch = &CheckUserInfoArrays(\@GATKTargetPaddedBedIntervalLists, "GATKTargetPaddedBedIntervalLists");

    open(PEDF, "<".$fileName) or die "Can't open ".$fileName.":$!, \n";    
     
    while (<PEDF>) {
	
	chomp $_;
	
	if ( ($. == 1) && ($_ =~/^\#/) ) { #Header present overwrite @pedigreeFileElements with header info
	
	    @pedigreeFileElements = split("\t", $'); #'
	    next;
	}
	if (m/^\s+$/) {	# Avoid blank lines
            next;
        }
	if (m/^\#/) {		# Avoid "#"
            next;
        }		
	if ($_ =~/(\S+)/) {	
	    
	    chomp($_);
	    my @lineInfo = split("\t",$_);	    #Loads pedigree file info
	    
##Need to parse familyID and sampleID separately since these have not been set yet
	    if ($lineInfo[0] =~/\S+/) { #familyID
		$familyID = $lineInfo[0];
	    }
	    else {
		print STDERR "File: ".$fileName." at line ".$.." cannot find FamilyID in column 1\n";
		exit;
	    }
	    if ($lineInfo[1] =~/\S+/) { #sampleID
		$sampleID = $lineInfo[1];		
		if ($userSampleIDsSwitch == 0) {
		    push(@sampleIDs, $lineInfo[1]); #Save sampleid info
		}
	    }
	    else {
		print STDERR "File: ".$fileName." at line ".$.." cannot find SampleID in column 2\n";
		exit;
	    }
	    
	    for (my $sampleElementsCounter=0;$sampleElementsCounter<scalar(@pedigreeFileElements);$sampleElementsCounter++) { #all pedigreeFileElements
		
		if ( defined($lineInfo[$sampleElementsCounter]) && ($lineInfo[$sampleElementsCounter] =~/\S+/) ) { #Check that we have an non blank entry
		    
		    my @elementInfo = split(";", $lineInfo[$sampleElementsCounter]); #Split element (if required)
		    
		    &CheckUniqueArrayElement(\@{ $sampleInfo{$familyID}{$sampleID}{$pedigreeFileElements[$sampleElementsCounter]} }, \@elementInfo); #Check if there are any new info and add it to sampleInfo if so. 
		    
		    if ($sampleInfo{$familyID}{$sampleID}{'Capture_kit'}) { #Add latest capture kit for each individual
			
			my $capture_kit = $sampleInfo{$familyID}{$sampleID}{$pedigreeFileElements[$sampleElementsCounter]}[-1]; #Use only the last capture kit since it should be the most interesting
			
			for my $supportedCaptureKit (keys %supportedCaptureKits) {
			    
			    if ($supportedCaptureKit eq $capture_kit) {
				
				if ($userExomeTargetBedInfileListsSwitch == 0) {
				    
				    $scriptParameter{$familyID}{$sampleID}{'exomeTargetBedInfileLists'} = $supportedCaptureKits{$supportedCaptureKit}.$referenceFileEndings{'exomeTargetBedInfileLists'}; #capture kit target in file_list
				}
				if ($userExomeTargetPaddedBedInfileListSwitch == 0) {
                                                                      
				    $scriptParameter{$familyID}{$sampleID}{'exomeTargetPaddedBedInfileLists'} = $supportedCaptureKits{$supportedCaptureKit}.$referenceFileEndings{'exomeTargetPaddedBedInfileLists'}; #capture kit padded target infile_list                               
				}
				if ($userExomeTargetPaddedBedIntervalListSwitch == 0) {
				                        
				    $scriptParameter{$familyID}{$sampleID}{'GATKTargetPaddedBedIntervalLists'} = $supportedCaptureKits{$supportedCaptureKit}.$referenceFileEndings{'GATKTargetPaddedBedIntervalLists'}; #capture kit padded target interval_list                          
				}
			    }
			}
		    }
		}
		else { #No entry in pedigre file element
		    
		    if ($sampleElementsCounter < 7) { #Only check mandatory elements 
			print STDERR $pedigreeFileElements[$sampleElementsCounter], "\t";
			print STDERR "File: ".$fileName." at line ".$.."\tcannot find '".$pedigreeFileElements[$sampleElementsCounter]."' entry in column ".$sampleElementsCounter, "\n";
			exit;
		    }  
		}
	    }
	}	
    }
    if ($userSampleIDsSwitch == 0) {
	@sampleIDs = sort(@sampleIDs); #Lexiographical sort to determine the correct order of ids indata
    }
    print STDOUT "Read pedigree file: ".$fileName, "\n";
    close(PEDF);
    return;
}

sub AddToJobID {
###Adds all previous jobIds per familyChainKey and chainKey to jobIDs string used to set the dependency in SLURM.

     my $familyIDChainKey = $_[0]; #familyID chain key
     my $chainKey = $_[1]; #sampleID or familyID chain key

     my $jobIDs = ""; #JobID string

     if ($jobID{$familyIDChainKey}{$chainKey}) {

     for (my $jobCounter=0;$jobCounter<scalar( @{ $jobID{$familyIDChainKey}{$chainKey} });$jobCounter++) {  #All previous jobIDs

          if ( ($jobCounter == 0) && (scalar( @{ $jobID{$familyIDChainKey}{$chainKey} }) == 1) ) {#Only 1 previous jobID 
               $jobIDs .= ":$jobID{$familyIDChainKey}{$chainKey}[$jobCounter]"; #first and last jobID start with ":" and end without ":"
          }
          elsif ($jobCounter == 0) { #First jobID
               $jobIDs .= ":$jobID{$familyIDChainKey}{$chainKey}[$jobCounter]:"; #first jobID start with :
          }
          elsif ($jobCounter eq (scalar( @{ $jobID{$familyIDChainKey}{$chainKey} }) -1) ) { #Last jobID
               $jobIDs .= "$jobID{$familyIDChainKey}{$chainKey}[$jobCounter]"; #last jobID finish without :
          }
          else { #JobIDs in the middle
               $jobIDs .= "$jobID{$familyIDChainKey}{$chainKey}[$jobCounter]:";
          }
     }
}
     return $jobIDs;
}

sub PushToJobID {
    
    my $familyIDChainKey = $_[0]; 
    my $sampleIDChainKey = $_[1];
    my $sampleID = $_[2];
    my $path = $_[3];
    my $chainKeyType = $_[4]; #Single, parallel, merged (familyID_sampleID)
    
    my $chainKey;
    
    if ($chainKeyType eq "Parallel") { #Push parallel jobs

	if ($scriptParameter{'analysisType'} eq "rapid" && $sampleInfo{$scriptParameter{'familyID'}}{$sampleID}{'pBwaMem'}{'sbatchBatchProcesses'}) {

	    for (my $sbatchCounter=0;$sbatchCounter<$sampleInfo{$scriptParameter{'familyID'}}{$sampleID}{'pBwaMem'}{'sbatchBatchProcesses'};$sbatchCounter++) {

		$chainKey = $sampleID."_parallel_".$path.$sbatchCounter; #Set key

		if ($jobID{$familyIDChainKey}{$chainKey}) { #job exists
		    
		    for (my $jobCounter=0;$jobCounter<scalar( @{ $jobID{$familyIDChainKey}{$chainKey} });$jobCounter++) { #All previous jobs i.e. jobs in this case equals to infiles in number
			
			push ( @{ $jobID{$familyIDChainKey}{$sampleIDChainKey} }, $jobID{$familyIDChainKey}{$chainKey}[$jobCounter]); #Add jobID to hash
		    }    
		}
	    }	
	    $sampleInfo{$scriptParameter{'familyID'}}{$sampleID}{'pBwaMem'}{'sbatchBatchProcesses'} = ();
	}
	else {

	    for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #All infiles
		
		$chainKey = $sampleID."_parallel_".$path.$infileCounter; #Set key
		
		if ($jobID{$familyIDChainKey}{$chainKey}) { #job exists
		    
		    for (my $jobCounter=0;$jobCounter<scalar( @{ $jobID{$familyIDChainKey}{$chainKey} });$jobCounter++) { #All previous jobs i.e. jobs in this case equals to infiles in number
			
			push ( @{ $jobID{$familyIDChainKey}{$sampleIDChainKey} }, $jobID{$familyIDChainKey}{$chainKey}[$jobCounter]); #Add jobID to hash
		    }    
		}
	    }
	}
    }
    elsif ( ($chainKeyType eq "Merged") || ($chainKeyType eq "Family_Merged")  ) { #Push merged jobs
	
	$chainKey = $familyIDChainKey."_".$sampleIDChainKey; #Set key
	
	if ($jobID{$familyIDChainKey}{$chainKey}) { #job exists
	    
	    for (my $jobCounter=0;$jobCounter<scalar( @{ $jobID{$familyIDChainKey}{$chainKey} });$jobCounter++) { #All previous jobs i.e. jobs in this case equals to infiles in number
		
		if ($chainKeyType eq "Family_Merged") { #Use $familyIDChainKey instead of $sampleIDChainKey
		    
		    push ( @{ $jobID{$familyIDChainKey}{$familyIDChainKey} }, $jobID{$familyIDChainKey}{$chainKey}[$jobCounter]); #Add jobID hash
		}
		else {
		    push ( @{ $jobID{$familyIDChainKey}{$sampleIDChainKey} }, $jobID{$familyIDChainKey}{$chainKey}[$jobCounter]); #Add jobID to hash
		}
	    }    
	}
    }
}

sub FIDSubmitJob {
###Submits all jobIDs to SLURM using SLURM dependencies. The trunk is the "MAIN path" and any subsequent splits into  branches "other paths" later is handled by adding relevant previous jobIDs to the new paths key in jobID{family_path_key} hash. The subroutine supports parallel job within each step and submission which do not leave any dependencies. Currently any path downstream of MAIN inherits the relevant previous jobIds, but it is not possible to merge to splited paths downstream of main to each other.

###Dependencies - $_[2]
    
##-1 = Not dependent on earlier scripts, and are self cul-de-sâcs
##0 = Not dependent on earlier scripts
##1 = Dependent on earlier scripts (within sampleID_path or familyID_path)
##2 = Dependent on earlier scripts (within sampleID_path or familyID_path), but are self cul-de-sâcs. 
##3 = Dependent on earlier scripts and executed in parallel within step
##4 = Dependent on earlier scripts and parallel scripts and executed in parallel within step 
##5 = Dependent on earlier scripts both family and sample and adds to both familyID and sampleId jobs
##6 = Not dependent on earlier scripts and adds to sampleId jobs, but sbatch is processed at family level i.e. affects all sampleID jobs e.g. building a reference
    
    my $sampleID = $_[0];
    my $familyID = $_[1];
    my $dependencies = $_[2]; 
    my $path = $_[3]; #Trunk or Branch
    my $sbatchFileName = $_[4]; #sbatch filename to submit.
    my $sbatchScriptTracker = $_[5]; #Track the number of parallel processes (e.g. sbatch scripts for a module)
    
    my $jobIDs=""; #Create string with all previous jobIDs
    my $jobIDsReturn; #Return jobID
    my $sampleIDChainKey = $sampleID."_".$path; #Sample chainkey
    my $familyIDChainKey = $familyID."_".$path; #Family chainkey
    my $sampleIDParallelChainKey = $sampleID."_parallel_".$path.$sbatchScriptTracker; #Sample parallel chainkey
    my $familyIDParallelChainKey = $familyID."_parallel_".$path.$sbatchScriptTracker; #Faimly parallel chainkey
    my $jobID; #The jobID that is returned from submission
    
    if ($dependencies == -1) { #Initiate chain - No dependencies, lonely program "sapling"
	
	$jobIDsReturn = `sbatch $sbatchFileName`; #No jobs have been run: submit
	($jobID) = ($jobIDsReturn =~ /Submitted batch job (\d+)/); #Just submitted jobID
    }
    if ($dependencies == 6) { #Initiate chain - No dependencies, adds to all sampleID(s)
	
	$jobIDsReturn = `sbatch $sbatchFileName`; #No jobs have been run: submit
	($jobID) = ($jobIDsReturn =~ /Submitted batch job (\d+)/); #Just submitted jobID
	
	for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {

	    my $sampleIDChainKey =  $sampleIDs[$sampleIDCounter]."_".$path;
	    push ( @{ $jobID{$familyIDChainKey}{$sampleIDChainKey} }, $jobID); #Add jobID to hash
	}
    }
    elsif ($dependencies == 0) { #Initiate chain - No dependencies, initiate Trunk (Main or other)
	
	$jobIDsReturn = `sbatch $sbatchFileName`; #No jobs have been run: submit
	($jobID) = ($jobIDsReturn =~ /Submitted batch job (\d+)/); #Just submitted jobID
	push ( @{ $jobID{$familyIDChainKey}{$sampleIDChainKey} }, $jobID); #Add jobID to hash
    }
    else { #Dependent on earlier scripts and/or parallel. JobIDs that do not leave dependencies do not get pushed to jobID hash
	
	if ($sampleID) { #Check jobs within sampleID (exception if dependencies = 5) 
	    
	    if ($dependencies == 5) { #Add familyID_sampleID jobs to current sampleID chain
		
		&PushToJobID($familyIDChainKey, $sampleIDChainKey, $sampleID, $path, "Merged");
	    }
	    if ( ($dependencies == 1) || ($dependencies == 2) ) { #not parallel jobs, but check if last job submission was parallel
		
		&PushToJobID($familyIDChainKey, $sampleIDChainKey, $sampleID, $path, "Parallel");
	    }
	    if ($path eq "MAIN") {
		
		if ( ($dependencies == 4) || ($dependencies == 3) ) { #Parallel jobs
		    
		    $jobIDs = &AddToJobID($familyIDChainKey, $sampleIDParallelChainKey); #Add to jobID string
		    
		    if ($jobID{$familyIDChainKey}{$sampleIDChainKey}) { #Check for previous single jobs - required to initiate broken chain with correct dependencies 
               
			$jobIDs .= &AddToJobID($familyIDChainKey, $sampleIDChainKey); #Add to jobID string
		    }
		    
		}
		else { #Previous job was a single job
		    
		    $jobIDs = &AddToJobID($familyIDChainKey, $sampleIDChainKey); #Add to jobID string
		}
	    }
	    if ($path ne "MAIN") { #Check for any previous jobIDs within path current PATH. Branch.
		
		if ($jobID{$familyIDChainKey}{$sampleIDChainKey}) { #second or later in branch chain
		    
		    $jobIDs = &AddToJobID($familyIDChainKey, $sampleIDChainKey);
		}
		elsif ($jobID{$familyID."_MAIN"}{$sampleID."_MAIN"}) { #Inherit from potential MAIN. Trunk
		    
		    $jobIDs = &AddToJobID($familyID."_MAIN", $sampleID."_MAIN");
		}
	    }     
	    if ($jobIDs) { #Previous jobs for chainkey exists
		$jobIDsReturn = `sbatch --dependency=afterok$jobIDs $sbatchFileName`; #Supply with dependency of previous jobs that this one is dependent on
		($jobID) = ($jobIDsReturn =~ /Submitted batch job (\d+)/); #Just submitted jobID
	    }
	    else { #No previous jobs
		$jobIDsReturn = `sbatch $sbatchFileName`; #No jobs have been run: submit
		($jobID) = ($jobIDsReturn =~ /Submitted batch job (\d+)/); #Just submitted jobID
	    }
	    if ($dependencies == 1) { #Ordinary job push to array
		
		@{ $jobID{$familyIDChainKey}{$sampleIDChainKey} } = (); #Clear latest familyID/sampleID chain submission
		
		##Clear all latest parallel jobs within chainkey
		for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) {
		    
		    my $sampleIDParallelChainKey = $sampleID."_parallel_".$path.$infileCounter; #Create key
		    
		    if ($jobID{$familyIDChainKey}{$sampleIDParallelChainKey}) { #Parallel job exists
			
			@{ $jobID{$familyIDChainKey}{$sampleIDParallelChainKey} } = (); #Clear latest familyID/sampleID chain submission
                    }
		}
		
		push ( @{ $jobID{$familyIDChainKey}{$sampleIDChainKey} }, $jobID); #Add jobID to hash{$sampleID}[]
	    }
	    if ( ($dependencies == 3) || ($dependencies == 4) ) { #Parallel job wait to push to array until all parallel jobs are finished within step
		
		push ( @{ $jobID{$familyIDChainKey}{$sampleIDParallelChainKey} }, $jobID); #Add jobID to hash
	    }
	    if ($dependencies == 5) { #Job dependent on both familyID and sampleID push to array
		
		@{ $jobID{$familyIDChainKey}{$familyIDChainKey."_".$sampleIDChainKey} } = (); #Clear latest familyID_sampleID chainkey
		@{ $jobID{$familyIDChainKey}{$sampleIDChainKey} } = (); #Clear latest sampleID chainkey
		push ( @{ $jobID{$familyIDChainKey}{$familyIDChainKey."_".$sampleIDChainKey} }, $jobID); #Add jobID to hash
	    }
	}
	else { #AFTER merging to familyID
	    
	    if ($dependencies == 5) { ##Add familyID_sampleID jobs to current familyID chain
		
		for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) { #Check jobs for sampleID          
		    
		    my $sampleIDChainKey = $sampleIDs[$sampleIDCounter]."_".$path; #Current chain
		    &PushToJobID($familyIDChainKey, $sampleIDChainKey, $sampleID, $path, "Family_Merged");
		}
	    }
	    if ( ($dependencies == 1) || ($dependencies == 2) ) { #not parallel jobs, but check if last job submission was parallel
		
		if ($jobID{$familyIDChainKey}{$familyIDParallelChainKey}) { #Parallel job exists
		    
		    for (my $jobCounter=0;$jobCounter<scalar( @{ $jobID{$familyIDChainKey}{$familyIDParallelChainKey} });$jobCounter++) {
			
			push ( @{ $jobID{$familyIDChainKey}{$familyIDChainKey} }, $jobID{$familyIDChainKey}{$familyIDParallelChainKey}[$jobCounter]); #Add jobID to hash{$} 
		    }
		}
	    }
	    if ( ($path eq "MAIN") && ($jobID{$familyIDChainKey}{$familyIDChainKey}) ) { #Check for any previous jobIDs within path MAIN. Test for previous must be done to allow initiating from broken chain. Trunk and not first in chain
		if ( ($dependencies == 4) || ($dependencies == 3) ) { #Parallel jobs
		    
		    $jobIDs = &AddToJobID($familyIDChainKey, $familyIDParallelChainKey); #Add to jobID string
		}
		else { #Previous job was a single job 
		    
		    $jobIDs = &AddToJobID($familyIDChainKey, $familyIDChainKey); #Add to jobID string
		}
	    }
	    elsif ($path eq "MAIN") { #First familyID MAIN chain 
		
          ##Add all previous jobId(s) from sampleId chainkey(s)
		for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {           
		    
		    my $sampleIDChainKey = $sampleIDs[$sampleIDCounter]."_".$path;
		    
		    if ($jobID{$familyIDChainKey}{$sampleIDChainKey}) {
			
			$jobIDs .= &AddToJobID($familyIDChainKey, $sampleIDChainKey); #Add to jobID string, while keeping previous additions
			
		    }
		    for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleIDs[$sampleIDCounter]} });$infileCounter++) {
			
			my $sampleIDParallelChainKey = $sampleIDs[$sampleIDCounter]."_parallel_".$path.$infileCounter; #Create key
			
			if ($jobID{$familyIDChainKey}{$sampleIDParallelChainKey}) { #Parallel job exists
			    
			    $jobIDs .= &AddToJobID($familyIDChainKey, $sampleIDParallelChainKey); #Add to jobID string, while keeping previous additions
			    
			}
		    }
		}
	    }
	    if ($path ne "MAIN" ) { #Check for any previous jobIDs within path current PATH. Branch
		
		if ($jobID{$familyIDChainKey}{$familyIDChainKey}) { #second or later in branch chain
		    
		    $jobIDs = &AddToJobID($familyIDChainKey, $familyIDChainKey); #Family chain
		}
		elsif ($jobID{$familyID."_MAIN"}{$familyID."_MAIN"}) { #Inherit from potential MAIN. Trunk
		    
		    $jobIDs = &AddToJobID($familyID."_MAIN", $familyID."_MAIN");
		}
		else { #First job in new path and first familyID MAIN chain 
		    
		    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {           
			
			my $familyIDChainKey = $familyID."_MAIN";
			my $sampleIDChainKey = $sampleIDs[$sampleIDCounter]."_MAIN";
			
			if ($jobID{$familyIDChainKey}{$sampleIDChainKey}) {
			    
			    $jobIDs .= &AddToJobID($familyIDChainKey, $sampleIDChainKey); 
			}
		    }
		}
	    }
	    if ($jobIDs) {
		$jobIDsReturn = `sbatch --dependency=afterok$jobIDs $sbatchFileName`; #Supply with dependency of previous jobs that this one is dependent on
		($jobID) = ($jobIDsReturn =~ /Submitted batch job (\d+)/);
	    }
	    else {
		$jobIDsReturn = `sbatch $sbatchFileName`; #No jobs have been run: submit
		($jobID) = ($jobIDsReturn =~ /Submitted batch job (\d+)/);
	    }
	    if ($dependencies == 1) { #Ordinary job push to array
		
		@{ $jobID{$familyIDChainKey}{$familyIDChainKey} } = (); #Clear latest familyID/sampleID chain submission
		@{ $jobID{$familyIDChainKey}{$familyIDParallelChainKey} } = (); #Clear latest familyID/sampleID chain submission
		push ( @{ $jobID{$familyIDChainKey}{$familyIDChainKey} }, $jobID); #Add jobID to hash{$sampleID}[]
	    }
	    if ( ($dependencies == 3) || ($dependencies == 4) ) { #Parallel job wait to push to array until all parallel jobs are finished within step
		
		push ( @{ $jobID{$familyIDChainKey}{$familyIDParallelChainKey} }, $jobID); #Add jobID to hash{$sampleID_parallel}[].
	    }    
	    if ($dependencies == 5) { #Job dependent on both familyID and sampleID push to array
		
		for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) { #Check jobs for sampleID          
		    
		    my $sampleIDChainKey = $sampleIDs[$sampleIDCounter]."_".$path; #Current chain
		    @{ $jobID{$familyIDChainKey}{$familyIDChainKey."_".$sampleIDChainKey} } = ();
		    @{ $jobID{$familyIDChainKey}{$familyIDChainKey} } = (); #Clear latest sampleID chainkey
		    push ( @{ $jobID{$familyIDChainKey}{$familyIDChainKey."_".$sampleIDChainKey} }, $jobID);   
		}
	    }
	}
    }
    if ($jobIDsReturn !~/\d+/) { #Catch errors since, propper sbatch submission should only return numbers
	print STDERR $jobIDsReturn."\n";
	print STDERR "\nMIP: Aborting run.\n\n";
	exit;
    }
    print STDOUT "Sbatch script submitted, job id: $jobID\n"; print MIPLOG "Sbatch script submitted, job id: $jobID\n";
    print STDOUT "To check status of job, please run \'squeue -j $jobID\'\n";print MIPLOG "To check status of job, please run \'squeue -j $jobID\'\n";
    print STDOUT "To cancel job, please run \'scancel $jobID\'\n";print MIPLOG "To cancel job, please run \'scancel $jobID\'\n";
    return;
}

sub NrofCoresPerSbatch {
##Set the number of cores to allocate per sbatch job
    
    my $nrCores = $_[0];
    
    if ($nrCores > $scriptParameter{'maximumCores'}) { #Set number of cores depending on how many lanes to process
	
	$nrCores = $scriptParameter{'maximumCores'}; #Set to max on cluster
    }
    return $nrCores;
}

sub InfilesReFormat {
###Reformat files for mosaik output, which have not yet been created into, correct format so that a sbatch script can be generated with the correct filenames.
    
    my $uncompressedFileCounter = 0; #Used to decide later if any inputfiles needs to be compressed before starting analysis 

    for my $sampleID (keys %infile) { #For every sampleID                                                                                       
	
        my $laneTracker=0; #Needed to be able to track when lanes are finished                                                                  
	
        for (my $infileCounter=0;$infileCounter<scalar( @ { $infile{$sampleID} });$infileCounter++) { #All inputfiles for all fastq dir and remakes format     
	    
            if ($infile{$sampleID}[$infileCounter] =~ /\/?([^\.\/]+\.[^\.]+)\.lane(\d+)_([12FfRr])\.fastq.gz/) { #Parse fastq.gz 'old' format, $2="lane", $3="Read direction"

		&AddInfileInfoOld($1, $2, $3, $sampleID, \$laneTracker, $infileCounter, "compressed"); 
            }
            elsif ($infile{$sampleID}[$infileCounter] =~ /\/?([^\.\/]+\.[^\.]+)\.lane(\d+)_([12FfRr])\.fastq/) { #Parse 'old' format                           

		&AddInfileInfoOld($1, $2, $3, $sampleID, \$laneTracker, $infileCounter, "uncompressed");
                $uncompressedFileCounter = 1; #File needs compression before starting analysis                               
            }
	    elsif ($infile{$sampleID}[$infileCounter] =~/(\d+)_(\d+)_([^_]+)_([^_]+)_index([^_]+)_(\d).fastq/) { #Parse 'new' no "index" format $1=lane, $2=date, $3=Flow-cell, $4=SampleID, $5=index,$6=direction                             
		
		my $compressedSwitch = &CheckGzipped(\$infile{$sampleID}[$infileCounter]);#Check gzipped or not
		
		if ($compressedSwitch eq "unCompressed") {
		    
		    $uncompressedFileCounter = "unCompressed"; #File needs compression before starting analysis. Note: All files are rechecked downstream and uncompressed ones are gzipped automatically           
		}
		&CheckSampleIDMatch($sampleID, $4, $infileCounter);
		&AddInfileInfo($1, $2, $3, $4, $5, $6, \$laneTracker, $infileCounter, $compressedSwitch);		
	    }
            elsif ($infile{$sampleID}[$infileCounter] =~/(\d+)_(\d+)_([^_]+)_([^_]+)_([^_]+)_(\d).fastq/) { #Parse 'new' no "index" format $1=lane, $2=date, $3=Flow-cell, $4=SampleID, $5=index,$6=direction                             
		
		my $compressedSwitch = &CheckGzipped(\$infile{$sampleID}[$infileCounter]);#Check gzipped or not

		if ($compressedSwitch eq "unCompressed") {
		   
		    $uncompressedFileCounter = "unCompressed"; #File needs compression before starting analysis. Note: All files are rechecked downstream and uncompressed ones are gzipped automatically           
		}
		&CheckSampleIDMatch($sampleID, $4, $infileCounter);
		&AddInfileInfo($1, $2, $3, $4, $5, $6, \$laneTracker, $infileCounter, $compressedSwitch);		
	    }
        }
    }
    return $uncompressedFileCounter;
}

sub CheckSampleIDMatch {
##Check that the sampleID provided and sampleID in infile name match.
    
    my $sampleID = $_[0]; #SampleID from user
    my $infileSampleID = $_[1]; #SampleID collect with regexp from infile
    my $infileCounter = $_[2];
    
    my %seen;
    $seen{$infileSampleID} = 1; #Add input as first increment
    
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {
	
	$seen{$sampleIDs[ $sampleIDCounter]}++;
    }
    unless ($seen{$infileSampleID} > 1) {

	print STDOUT "\n".$sampleID." supplied and sampleID ".$infileSampleID." found in file : ".$infile{$sampleID}[$infileCounter]." does not match. Please rename file to match sampleID: ".$sampleID."\n\n";
	exit;
    }
}

sub AddInfileInfoOld {
##Adds information derived from infile name to sampleInfo hash. Tracks the number of lanes sequenced and checks unique array elementents.  
    
    my $dateFlowCell = $_[0];
    my $lane = $_[1];
    my $direction = $_[2];
    my $sampleID = $_[3];
    my $laneTrackerRef = $_[4];
    my $infileCounter = $_[5];
    my $compressedInfo = $_[6]; #.fastq.gz or .fastq info governs zcat or cat downstream

    my $readFile;

    if ($compressedInfo eq "compressed") {
	$readFile = "zcat"; #read file in compressed format
    }
    else {
	$readFile = "cat"; #read file in uncompressed format
    }
    
    my $seqLengthRegExp = q?perl -ne 'if ($_!~/@/) {chomp($_);my $seqLength = length($_);print $seqLength;last;}' ?; #Prints sequence length and exits

    if ($direction == 1) { #Read 1
	
	push( @{$lane{$sampleID}}, $lane); #Lane
	$infilesLaneNoEnding{$sampleID}[$$laneTrackerRef]= $dateFlowCell.".".$lane; #Save old format in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry into hash for both strands and .ending is removed (.fastq).
	$sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesLaneNoEnding{ $sampleID }[$$laneTrackerRef]}{'sequenceRunType'} = "Single-end"; #Single-end until proven otherwise
	$$laneTrackerRef++;
    }
    if ($direction == 2) { #2nd read direction
	$sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesLaneNoEnding{ $sampleID }[$$laneTrackerRef-1]}{'sequenceRunType'} = "Paired-end"; #$laneTracker -1 since it gets incremented after direction eq 1. 
    }
    
    $infilesBothStrandsNoEnding{ $sampleID }[$infileCounter]= $dateFlowCell.".".$lane."_".$direction; #Save new format in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry per strand and .ending is removed (.fastq).
    $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesBothStrandsNoEnding{ $sampleID }[$infileCounter]}{'originalFileName'} = $infile{$sampleID}[$infileCounter]; #Original fileName
    $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesBothStrandsNoEnding{ $sampleID }[$infileCounter]}{'originalFileNameNoEnding'} = $dateFlowCell."lane".$lane."_".$direction; #Original fileName, but no ending
    $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesBothStrandsNoEnding{ $sampleID }[$infileCounter]}{'sampleBarcode'} = "X"; #Save barcode, but not defined
    $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesBothStrandsNoEnding{ $sampleID }[$infileCounter]}{'runBarcode'} = $lane."_".$dateFlowCell; #Save run barcode
    &CheckUniqueArrayElement(\@{ $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesBothStrandsNoEnding{ $sampleID }[$infileCounter]}{'readDirection'} }, \$direction); #Check if there are any new info and add it to sampleInfo if so. 
    $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesBothStrandsNoEnding{ $sampleID }[$infileCounter]}{'sequenceLength'} = `cd $indirpath{$sampleID};$readFile $infile{$sampleID}[$infileCounter] | $seqLengthRegExp;`; #Collect sequence length
}

sub AddInfileInfo {
##Adds information derived from infile name to sampleInfo hash. Tracks the number of lanes sequenced and checks unique array elementents.  

    my $lane = $_[0];
    my $date = $_[1];
    my $flowCell = $_[2];
    my $sampleID = $_[3];
    my $index = $_[4];
    my $direction = $_[5];
    my $laneTrackerRef = $_[6];
    my $infileCounter = $_[7];
    my $compressedInfo = $_[8]; #.fastq.gz or .fastq info governs zcat or cat downstream

    my $readFile;

    if ($compressedInfo eq "compressed") {
	$readFile = "zcat"; #read file in compressed format
    }
    else {
	$readFile = "cat"; #read file in uncompressed format
    }
    
    my $seqLengthRegExp = q?perl -ne 'if ($_!~/@/) {chomp($_);my $seqLength = length($_);print $seqLength;last;}' ?; #Prints sequence length and exits

    if ($direction == 1) { #Read 1
	
	push( @{$lane{$sampleID}}, $lane); #Lane
	$infilesLaneNoEnding{$sampleID}[$$laneTrackerRef]= $sampleID.".".$date."_".$flowCell."_".$index.".lane".$1; #Save new format (sampleID_date_flow-cell_index_lane) in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry into hash for both strands and .ending is removed (.fastq).
	$sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesLaneNoEnding{ $sampleID }[$$laneTrackerRef]}{'sequenceRunType'} = "Single-end"; #Single-end until proven otherwise
	$$laneTrackerRef++;
    }
    if ($direction == 2) { #2nd read direction
	$sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesLaneNoEnding{ $sampleID }[$$laneTrackerRef-1]}{'sequenceRunType'} = "Paired-end"; #$laneTracker -1 since it gets incremented after direction eq 1. 
    }
    
    $infilesBothStrandsNoEnding{ $sampleID }[$infileCounter]= $sampleID.".".$date."_".$flowCell."_".$index.".lane".$1."_".$direction; #Save new format in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry per strand and .ending is removed (.fastq).
    $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesBothStrandsNoEnding{ $sampleID }[$infileCounter]}{'originalFileName'} = $infile{$sampleID}[$infileCounter]; #Original fileName
    $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesBothStrandsNoEnding{ $sampleID }[$infileCounter]}{'originalFileNameNoEnding'} = $1."_".$date."_".$flowCell."_".$sampleID."_".$index."_".$direction; #Original fileName, but no ending
    $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesBothStrandsNoEnding{ $sampleID }[$infileCounter]}{'lane'} = $1; #Save sample lane                  
    $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesBothStrandsNoEnding{ $sampleID }[$infileCounter]}{'date'} = $date; #Save Sequence run date          
    $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesBothStrandsNoEnding{ $sampleID }[$infileCounter]}{'flow-cell'} = $flowCell; #Save Sequence flow-cell        
    $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesBothStrandsNoEnding{ $sampleID }[$infileCounter]}{'sampleBarcode'} = $index; #Save sample barcode
    $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesBothStrandsNoEnding{ $sampleID }[$infileCounter]}{'runBarcode'} = $date."_".$flowCell."_".$1."_".$index; #Save run barcode
    
    &CheckUniqueArrayElement(\@{  $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesBothStrandsNoEnding{ $sampleID }[$infileCounter]}{'ReadDirection'} }, \$direction); #Check if there are any new info and add it to sampleInfo if so.   
    $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesBothStrandsNoEnding{ $sampleID }[$infileCounter]}{'sequenceLength'} = `cd $indirpath{$sampleID};$readFile $infile{$sampleID}[$infileCounter] | $seqLengthRegExp;`; #Collect sequence length
}

sub Checkfnexists {
##Check if a file with with a filename consisting of $filePathRef.$fileCounter.$fileEndingRef exist. If so bumps the version number and return new filename.
    
    my $filePathRef = $_[0];
    my $fileEndingRef = $_[1];
    my $fileNameTrackerRef = $_[2];

    my $fileName;
    
    $$fileNameTrackerRef = 0; #Nr of sbatch scripts with identical filenames
  
    for (my $fileCounter=0;$fileCounter<9999;$fileCounter++) { #Number of possible files with the same name
	
	$fileName = $$filePathRef.$fileCounter.$$fileEndingRef; #filename, filenr and fileending
	$$fileNameTrackerRef = $fileCounter; #Nr of sbatch scripts with identical filenames
	if (-f $fileName) { #if file exists 
	}
	else {
	    last; #Exit loop 
	}	
    }
    $$filePathRef = $fileName; #Transfer to global variable
    return;
}

sub DefineParametersPath {
###Defines all attributes of a parameter, so that the correct value can be set and added to %scriptparameter later

    my $parameterName = $_[0]; #ParameterName
    my $parameterDefault = $_[1]; #Default setting
    my $associatedProgram = $_[2]; #The parameters program
    my $existsCheck = $_[3]; #Check if intendent file exists in reference directory
    my $buildFile = $_[4]; #Autovivication of file if it does not exists (yes or no)
    
    $parameter{$parameterName} = {
	'type' => "path",
	'value' => "nocmdinput",
	'default' => $parameterDefault,
	'associatedProgram' => $associatedProgram,
	'existsCheck' => $existsCheck,
	'buildFile' => $buildFile,
    };
    
    push(@orderParameters, $parameterName); #Add to enable later evaluation of parameters in proper order & write to master file
    
    return;
}

sub DefineParameters {
###Defines all attributes of a parameter, so that the correct value can be set and added to %scriptparameter later

    my $parameterName = $_[0]; #ParameterName
    my $parameterType = $_[1]; #MIP or program
    my $parameterDefault = $_[2]; #Default setting
    my $associatedProgram = $_[3]; #The parameters program
    my $fileEnding = $_[4]; #The filending after the module has been run
    my $parameterChain = $_[5]; #The chain to which the program belongs to
    my @programNamePath = split(":", $_[6]) if (defined($_[6])); #The path name of the program(s) for each sbatch script

    if (defined($programNamePath[0])) {
	
	$parameter{$parameterName} = {
	    'type' => $parameterType,
	    'value' => "nocmdinput",
	    'default' => $parameterDefault,
	    'associatedProgram' => $associatedProgram,
	    'fileEnding' => $fileEnding,
	    'chain' => $parameterChain,
	    'programNamePath' => \@programNamePath,
	};
    }
    else {
	
	$parameter{$parameterName} = {
	    'type' => $parameterType,
	    'value' => "nocmdinput",
	    'default' => $parameterDefault,
	    'associatedProgram' => $associatedProgram,
	    'fileEnding' => $fileEnding,
	    'chain' => $parameterChain,
	};
    }
    
    push(@orderParameters, $parameterName); #Add to enable later evaluation of parameters in proper order & write to master file
    
    return;
}

sub AddToScriptParameter {
###Checks and sets user input or default values to scriptPrameters
    
    my $parameterName = $_[0]; #ParameterName
    my $parameterValue = $_[1]; #Parameter to evaluate
    my $parameterType = $_[2]; #Path or program
    my $parameterDefault = $_[3]; #Default setting
    my @associatedPrograms = split(/,/, $_[4]); #The parameters program(s)
    my $parameterExistsCheck = $_[5]; #Check if intendent file exists in reference directory
   
##Validation
    #print "parameterName: ".$parameterName, "\n";
    #print "parameterValue: ".$parameterValue, "\n";
    #print "parameterType: ".$parameterType, "\n";
    #print "parameterDefault: ".$parameterDefault, "\n";
    #foreach my $associatedProgram (@associatedPrograms) {
#	print "associatedProgram: ".$associatedProgram, "\n";
 #   }
    
    foreach my $associatedProgram (@associatedPrograms) { #Check all programs that use parameter

	my $parameterSetSwitch = 0;

	if (defined($scriptParameter{$associatedProgram}) && ($scriptParameter{$associatedProgram} > 0) ) { #Only add active programs parameters
	    
	    $parameterSetSwitch = 1;

	    if ($parameterType eq "path") {
		
		if ($parameterValue eq "nocmdinput") { #No input from cmd
		    
		    if (defined($scriptParameter{$parameterName})) { #Input from config file
			
			if ($parameterName eq "sampleIDs") { #SampleIDs is a comma separated list 
			   
			    @sampleIDs = split(/,/, $scriptParameter{'sampleIDs'}); #Transfer to array
			}
			if ($parameterName eq "inFilesDirs") { #inFilesDirs is a comma separated list 

			    @inFilesDirs = split(/,/, $scriptParameter{'inFilesDirs'}); #Transfer to array			    
			} 
			if ($parameterName eq "picardToolsMergeSamFilesPrevious") {

			    @picardToolsMergeSamFilesPrevious = split(/,/, $scriptParameter{'picardToolsMergeSamFilesPrevious'}); #Transfer to array
			}
			if ($parameterName eq "exomeTargetBedInfileLists") { #exomeTargetBedInfileListss is a comma separated list 

			    &SetTargetandAutoBuild(\@sampleIDs, \$parameterName, \$referenceFileEndings{'exomeTargetBedInfileLists'});
			}
			if ($parameterName eq "exomeTargetPaddedBedInfileLists") { #exomeTargetPaddedBedInfileLists is a comma separated list 
			    
			    &SetTargetandAutoBuild(\@sampleIDs, \$parameterName, \$referenceFileEndings{'exomeTargetPaddedBedInfileLists'});
			}
			if ($parameterName eq "GATKTargetPaddedBedIntervalLists") { #GATKTargetPaddedBedIntervalLists is a comma separated list 
			    
			    &SetTargetandAutoBuild(\@sampleIDs, \$parameterName, \$referenceFileEndings{'GATKTargetPaddedBedIntervalLists'});
			}
			if ($parameterName eq "annovarTableNames") { #Input from config file
			    
			    @annovarTableNames = split(/,/, $scriptParameter{'annovarTableNames'}); #Transfer to array
			}
			if ($parameterName eq "humanGenomeReference") {
			    
			    ($humanGenomeReferenceVersion, $humanGenomeReferenceSource, $humanGenomeReferenceNameNoEnding) = ParseHumanGenomeReference($scriptParameter{'humanGenomeReference'});
			}
			if ($parameterName eq "pedigreeFile") {
			    
			    &ReadPlinkPedigreeFile($scriptParameter{'pedigreeFile'});
			}
		    }
		    elsif ($parameterDefault ne "nodefault") { #add default value
			
			if ($parameterName eq "inFilesDirs") {
			    
			    for (my $indirectoryCount=0;$indirectoryCount<scalar(@sampleIDs);$indirectoryCount++) {
				
				push(@inFilesDirs, $scriptParameter{'clusterConstantPath'}."/".$scriptParameter{'analysisType'}."/".$sampleIDs[$indirectoryCount]."/fastq");					
			    }
			    &EnableArrayParameter(\@inFilesDirs, \$parameterName);
			}
			elsif ($parameterName eq "exomeTargetBedInfileLists") { #Note that potential pedigree files entries will be updated with GenomeReferenceSource and version here
			    
			    &SetTargetandAutoBuild(\@sampleIDs, \$parameterName, \$referenceFileEndings{'exomeTargetBedInfileLists'});
			}
			elsif ($parameterName eq "exomeTargetPaddedBedInfileLists") { #Note that potential pedigree files entries will be updated with GenomeReferenceSource and version here
			    &SetTargetandAutoBuild(\@sampleIDs, \$parameterName, \$referenceFileEndings{'exomeTargetPaddedBedInfileLists'});
			}
			elsif ($parameterName eq "GATKTargetPaddedBedIntervalLists") { #Note that potential pedigree files entries will be updated with GenomeReferenceSource and version here 
			
			    &SetTargetandAutoBuild(\@sampleIDs, \$parameterName, \$referenceFileEndings{'GATKTargetPaddedBedIntervalLists'});
			}
			elsif ($parameterName eq "annovarTableNames") {

			    @annovarTableNames = ("refGene", "mce46way", "gerp++elem", "segdup", "tfbs", "mirna", "snp137NonFlagged", "1000g2012apr_all", "esp6500si_all", "ljb2_sift", "ljb2_pp2hdiv", "ljb2_pp2hvar", "ljb2_mt", "ljb2_lrt", "ljb2_gerp++","ljb2_phylop"); #Set default annovar table names
			    &EnableArrayParameter(\@inFilesDirs, \$parameterName);
			}
			else {
			    
			    $scriptParameter{$parameterName} = $parameterDefault; #Set default value
			}
		    }
		    else {
			
			if ($parameterName eq "picardToolsMergeSamFilesPrevious") {  

			    @picardToolsMergeSamFilesPrevious = (); #Empty to not add a 0 as a value, which will cause errors in later conditions use
			}
			elsif ( ($parameterName eq "sampleIDs") && (defined($scriptParameter{'pedigreeFile'})) ) { #Special case 
			    @sampleIDs = (); #Empty to not add a 0 as a value, which will cause errors in later conditions use
			}
			elsif ($parameterName eq "pedigreeFile") { #Special case - do nothing
			}
			elsif ($parameterName eq "mosaikAlignReference") { #Special case - do nothing, since file can be created by MIP from the humanGenomeReference if required
			}
			elsif ( ($parameterName eq "bwaMemRapidDb") && ($scriptParameter{'analysisType'} ne "rapid")) { #Do nothing since file is not required unless rapid mode is enabled
			}
			else {
			    
			    print STDERR $USAGE, "\n";
			    print STDERR "\nSupply '-".$parameterName."' if you want to run ".$associatedProgram, "\n\n";
			    exit;
			}
		    }
		}
		else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
		    
		    if ($parameterName eq "sampleIDs") {	  

			&EnableArrayParameter(\@sampleIDs, \$parameterName);  
		    }
		    elsif ($parameterName eq "inFilesDirs") {	    

			&EnableArrayParameter(\@inFilesDirs, \$parameterName);
		    }
		    elsif ($parameterName eq "picardToolsMergeSamFilesPrevious") {
			
			&EnableArrayParameter(\@picardToolsMergeSamFilesPrevious, \$parameterName);
		    }
		    elsif ($parameterName eq "exomeTargetBedInfileLists") {	    
			
			&EnableArrayParameter(\@exomeTargetBedInfileLists, \$parameterName);
			&CompareArrayElements(\@sampleIDs, \@exomeTargetBedInfileLists, "sampleIDs", $parameterName);
			&SetAutoBuildAndScriptParameterPerSample(\@sampleIDs, \@exomeTargetBedInfileLists, \$parameterName);
		    }
		    elsif ($parameterName eq "exomeTargetPaddedBedInfileLists") {	    
			
			&EnableArrayParameter(\@exomeTargetPaddedBedInfileLists, \$parameterName);
			&CompareArrayElements(\@sampleIDs, \@exomeTargetPaddedBedInfileLists, "sampleIDs", $parameterName);
			&SetAutoBuildAndScriptParameterPerSample(\@sampleIDs, \@exomeTargetPaddedBedInfileLists, \$parameterName);
		    }
		    elsif ($parameterName eq "GATKTargetPaddedBedIntervalLists") {	    
			
			&EnableArrayParameter(\@GATKTargetPaddedBedIntervalLists, \$parameterName);
			&CompareArrayElements(\@sampleIDs, \@GATKTargetPaddedBedIntervalLists, "sampleIDs", $parameterName);
			&SetAutoBuildAndScriptParameterPerSample(\@sampleIDs, \@GATKTargetPaddedBedIntervalLists, \$parameterName);
		    }
		    elsif ($parameterName eq "annovarTableNames") {
			
			@annovarTableNames = ("refGene", "mce46way", "gerp++elem", "segdup", "tfbs", "mirna", "snp137NonFlagged", "1000g2012apr_all", "esp6500si_all", "ljb2_sift", "ljb2_pp2hdiv", "ljb2_pp2hvar", "ljb2_mt", "ljb2_lrt", "ljb2_gerp++","ljb2_phylop"); #Set default annovar table names
			&EnableArrayParameter(\@inFilesDirs, \$parameterName);
		    }
		    elsif ($parameterName eq "pedigreeFile") { #Must come after arrays that can be populated from pedigree file to not overwrite user cmd input 

			$scriptParameter{$parameterName} = $parameterValue;
			&ReadPlinkPedigreeFile($scriptParameter{'pedigreeFile'});
		    }
		    else {
			
			if ($parameterName eq "humanGenomeReference") {
			    
			    ($humanGenomeReferenceVersion, $humanGenomeReferenceSource, $humanGenomeReferenceNameNoEnding) = &ParseHumanGenomeReference($parameterValue);
			}
			$scriptParameter{$parameterName} = $parameterValue;
		    }
		}
		if ( $parameterExistsCheck && ($parameterExistsCheck eq "directory") ) { #Check dir existence
		    
		    if ($parameterName eq "inFilesDirs") {
			
			@inFilesDirs = split(/,/, join(',', @inFilesDirs));
			
			for (my $indirectoryCount=0;$indirectoryCount<scalar(@inFilesDirs);$indirectoryCount++) {
			    
			    &CheckExistance(\$inFilesDirs[$indirectoryCount], \$parameterName, "d");
			}
		    }
		    else {
			
			&CheckExistance(\$scriptParameter{$parameterName}, \$parameterName, "d");

			if ($parameterName eq "genomeAnalysisToolKitPath") { #To enable addition of version to sampleInfo
			    
			    if ($scriptParameter{$parameterName}=~/GenomeAnalysisTK-([^,]+)/) {
				$sampleInfo{$scriptParameter{'familyID'}}{'program'}{"GATK"}{'Version'} = $1;
			    }
			}
			if ($parameterName eq "picardToolsPath") { #To enable addition of version to sampleInfo                                                                       

                            if ($scriptParameter{$parameterName}=~/picard-tools-([^,]+)/) {
                                $sampleInfo{$scriptParameter{'familyID'}}{'program'}{"PicardTools"}{'Version'} = $1;
                            }
                        }
		    }
		}
		elsif ( ($parameterExistsCheck) && ($parameterExistsCheck eq "file") && (defined($scriptParameter{$parameterName})) ) { #Check file existence in reference directory
		    
		    if ($parameterName eq "mosaikJumpDbStub") {
			
			&CheckFileEndingsToBeBuilt(\@mosaikJumpDbStubFileEndings, $parameterName); 
		    }
		    elsif ($parameterName eq "bwaBuildReference") {
			
			&CheckFileEndingsToBeBuilt(\@bwaBuildReferenceFileEndings, $parameterName);
		    }
		    elsif ($parameterName eq "humanGenomeReference") {
			
			&CheckExistance(\($scriptParameter{'referencesDir'}."/".$scriptParameter{$parameterName}), \$parameterName, "f");#Check reference genome
			$sampleInfo{$scriptParameter{'familyID'}}{"HumanGenomeBuild"}{'Path'} = $scriptParameter{'referencesDir'}."/".$scriptParameter{$parameterName};
			$sampleInfo{$scriptParameter{'familyID'}}{"HumanGenomeBuild"}{'Source'} = $humanGenomeReferenceSource;
			$sampleInfo{$scriptParameter{'familyID'}}{"HumanGenomeBuild"}{'Version'} = $humanGenomeReferenceVersion;

			#Enable autoBuild of metafiles 	       
			$parameter{$parameterName.".dict"}{'buildFile'} = "yesAutoBuild";
			$parameter{$parameterName.".fasta.fai"}{'buildFile'} = "yesAutoBuild";

			for (my $fileEndingsCounter=0;$fileEndingsCounter<scalar(@humanGenomeReferenceFileEndings);$fileEndingsCounter++) {
			
			    my $intendedFilePathRef = \($scriptParameter{'referencesDir'}."/".$humanGenomeReferenceNameNoEnding.$humanGenomeReferenceFileEndings[$fileEndingsCounter]);
			    
			    &CheckExistance($intendedFilePathRef, \($parameterName.$humanGenomeReferenceFileEndings[$fileEndingsCounter]), "f");
			}
		
			if ($parameter{$parameterName.".dict"}{'buildFile'} eq 0) {
			    ##Collect sequence contigs from human reference ".dict" file since it exists
			    &CollectSeqContigs(); #Preparation for future changes but not active yet
			}
		    }
		    elsif ( ($parameterName eq "exomeTargetBedInfileLists") || ($parameterName eq "exomeTargetPaddedBedInfileLists") || ($parameterName eq "GATKTargetPaddedBedIntervalLists") ) {
			
			for (my $sampleIDsCounter=0;$sampleIDsCounter<scalar(@sampleIDs);$sampleIDsCounter++) { #All sampleIDs
			    
			    &CheckExistance(\($scriptParameter{'referencesDir'}."/".$scriptParameter{ $scriptParameter{'familyID'} }{$sampleIDs[$sampleIDsCounter]}{$parameterName}), \$parameterName, "f", \$sampleIDs[$sampleIDsCounter]);			    
			
			    my $exomeTargetBedFileNoEnding = &RemoveFileEnding(\$scriptParameter{ $scriptParameter{'familyID'} }{$sampleIDs[$sampleIDsCounter]}{$parameterName} , $referenceFileEndings{$parameterName}); #Remove ".fileending" from reference filename
			    &CheckTargetExistFileBed(\$exomeTargetBedFileNoEnding, $parameterName);
			}
			undef($scriptParameter{$parameterName}); #Remove parameter to avoid unnecessary print to STDOUT and config
		    }
		    elsif ($parameterName eq "annovarTableNames") {
			
			&DefineAnnovarTables(); #Set all AnnovarTables properties
			my $intendedFilePathRef;
		
			for (my $tableNamesCounter=0;$tableNamesCounter<scalar(@annovarTableNames);$tableNamesCounter++) { #All AnnovarTables
		 
			    if (defined($annovarTables{$annovarTableNames[$tableNamesCounter]}{'file'})) {
				
				for (my $filesCounter=0;$filesCounter<scalar(@{$annovarTables{ $annovarTableNames[$tableNamesCounter] }{'file'}});$filesCounter++) { #All annovarTables file(s), some tables have multiple files downloaded from the same call
				      $intendedFilePathRef = \($scriptParameter{'annovarPath'}."/humandb/".$annovarTables{ $annovarTableNames[$tableNamesCounter] }{'file'}[$filesCounter]);
				     &CheckExistance($intendedFilePathRef, \$annovarTableNames[$tableNamesCounter], "f");
				}
			    }
			    elsif (defined($annovarTables{ $annovarTableNames[$tableNamesCounter] }{'ucscAlias'})){
				
				$intendedFilePathRef = \($scriptParameter{'annovarPath'}."/humandb/".$scriptParameter{'annovarGenomeBuildVersion'}."_".$annovarTables{ $annovarTableNames[$tableNamesCounter] }{'ucscAlias'}.".txt");
				&CheckExistance($intendedFilePathRef, \$annovarTableNames[$tableNamesCounter], "f");
			    }
			    else {
				
				$intendedFilePathRef = \($scriptParameter{'annovarPath'}."/humandb/".$scriptParameter{'annovarGenomeBuildVersion'}."_".$annovarTableNames[$tableNamesCounter].".txt");
				&CheckExistance($intendedFilePathRef, \$annovarTableNames[$tableNamesCounter], "f");
			    }
			    
			}
		    }
		    elsif ($parameterName eq "configFile") {  #Do nothing since file existence is checked by &LoadYAML
		    }
		    elsif ($parameterName eq "pedigreeFile") { #Do nothing since file existence is checked by ReadPlinkPedigreeFile
		    }
		    elsif ($parameterName eq "sampleInfoFile") {

			if (defined($scriptParameter{'sampleInfoFile'})) {

			    if (-f $scriptParameter{'sampleInfoFile'}) {

				%sampleInfo = &LoadYAML($scriptParameter{'sampleInfoFile'}); #Load parameters from previous run from sampleInfoFile	  
				
			    }
			    if (defined($scriptParameter{'pedigreeFile'}) ) {
				
				$sampleInfo{$scriptParameter{'familyID'}}{'pedigreeFile'}{'Path'} = $scriptParameter{'pedigreeFile'}; #Add pedigreeFile to sampleInfo
				$sampleInfo{$scriptParameter{'familyID'}}{'pedigreeFileAnalysis'}{'Path'} = $scriptParameter{'outDataDir'}."/".$scriptParameter{'familyID'}."/qc_pedigree.yaml"; #Add pedigreeFile info used in this analysis to SampleInfoFile
			    }
			} 
		    }
		    elsif ( ($parameterName eq "bwaMemRapidDb") && ($scriptParameter{'analysisType'} ne "rapid")) { #Do nothing since file is not required unless rapid mode is enabled
		    }
		    elsif ( ($parameterName eq "GATKHaploTypeCallerRefBAMInfile") && ($scriptParameter{'analysisType'} =~/rapid|genomes/) ) { #Do nothing since file is not required unless exome mode is enabled
		    }
		    else {
			
			&CheckExistance(\($scriptParameter{'referencesDir'}."/".$scriptParameter{$parameterName}), \$parameterName, "f");}
		    
		}
	    }
	    
	    if ($parameterType eq "MIP") {
		
		if ($parameterValue eq "nocmdinput") { #No input from cmd
		    
		    if (defined($scriptParameter{$parameterName})) { #Input from config file - do nothing
		    }
		    elsif ($parameterDefault ne "nodefault") {
		
			$scriptParameter{$parameterName} = $parameterDefault; #Set default value
		    }
		    else {
			
			if ($parameterName eq "aligner") { #set to "nocmdinput"
			
			    $scriptParameter{'aligner'} = "nocmdinput";
			}
			else {
			    
			    print STDERR $USAGE, "\n";
			    print STDERR "\nSupply '-".$parameterName."' if you want to run ".$associatedProgram, "\n\n";
			    exit;
			}
		    }
		}
		else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
		    
		    $scriptParameter{$parameterName} = $parameterValue; 
		}
	    }
	    
	    if ( $parameterType eq "program") {

		if($parameterValue eq "nocmdinput") { #No input from cmd
		    
		    if (defined($scriptParameter{$parameterName})) { #Input from config file - do nothing
			
			if ($parameterName eq "ImportantDbFileOutFiles") {
			    
			    @ImportantDbFileOutFiles = split(/,/, $scriptParameter{'ImportantDbFileOutFiles'});
			}
		    }
		    elsif ($parameterDefault ne "nodefault") {
			
			if ($parameterName eq "ImportantDbFileOutFiles") {
	
			    my $inDirectoryResearch = $scriptParameter{'outDataDir'}."/".$scriptParameter{'familyID'}."/".$scriptParameter{'aligner'}."/GATK/candidates/ranking";
			    my $inDirectoryClinical = $scriptParameter{'outDataDir'}."/".$scriptParameter{'familyID'}."/".$scriptParameter{'aligner'}."/GATK/candidates/ranking/clinical"; 
			    @ImportantDbFileOutFiles = ($inDirectoryResearch."/".$scriptParameter{'familyID'}."_orphan.selectVariants", $inDirectoryClinical."/".$scriptParameter{'familyID'}.".selectVariants");
			    $scriptParameter{'ImportantDbFileOutFiles'} = join(",", @ImportantDbFileOutFiles);
			}
			else {
			    
			    $scriptParameter{$parameterName} = $parameterDefault; #Set default value
			}
		    }
		}
		else {

		    if ($parameterName eq "ImportantDbFileOutFiles") {

			&EnableArrayParameter(\@ImportantDbFileOutFiles, \$parameterName); #Enables comma separated list of sample IDs from user supplied cmd info
		    }
		    else {

			$scriptParameter{$parameterName} = $parameterValue;
		    }
		}
		
		if (defined($parameter{$parameterName}{'programNamePath'}[0])) { #Code for checking commands in your path and executable
		
		    if ($scriptParameter{$parameterName} > 0) { #Only check path(s) for active programs
		
			for (my $programNamePathCounter=0;$programNamePathCounter<scalar(@{$parameter{$parameterName}{'programNamePath'}});$programNamePathCounter++) { #Check all binaries for sbatch program
			    
			    if ( grep { -x "$_/".$parameter{$parameterName}{'programNamePath'}[$programNamePathCounter]} split(/:/,$ENV{PATH}) ) {
				
				print STDOUT "ProgramCheck: ".$parameter{$parameterName}{'programNamePath'}[$programNamePathCounter]." installed\n"; 
			    }
			    else {
				print STDOUT "Warning: Could not detect ".$parameter{$parameterName}{'programNamePath'}[$programNamePathCounter]." in your Path\n";
				exit;
			    }
			}
		    }
		}
	    }
	    
	    if ($parameterName eq "aligner") {
		
		if ( ($scriptParameter{'pMosaikBuild'} > 0) || ($scriptParameter{'pMosaikAlign'} > 0)) { #Mosaik track
		    
		    if ( ($scriptParameter{'pBwaAln'} == 0) && ($scriptParameter{'pBwaSampe'} == 0)) {
			
			if ($scriptParameter{'aligner'} eq "bwa") {
			    $scriptParameter{'aligner'} = "mosaik";
			}
		    }
		    else {
		
			print STDERR $USAGE, "\n";
			print STDERR "You have to choose either mosaik or bwa to perform alignments or specify which aligner (-aligner 'mosaik' or 'bwa') was used if you want to only run programs after alignment.", "\n\n";
			exit;
		    }
		}
		elsif ( ($scriptParameter{'pBwaAln'} > 0) || ($scriptParameter{'pBwaSampe'} > 0)) { #BWA track
		    
		    if ( ($scriptParameter{'aligner'} eq "mosaik") || ($scriptParameter{'aligner'} =~ /bwa/i) ) {

			$scriptParameter{'aligner'} = "bwa";
		    }
		    else {
			print STDERR $USAGE, "\n";
			print STDERR "You have to choose either mosaik or bwa to perform alignments or specify which aligner (-aligner 'mosaik' or 'bwa') was used if you want to only run programs after alignment.", "\n\n";
			exit;
		    }
		}
		elsif ($scriptParameter{'aligner'} eq "nocmdinput") {
		    print STDERR $USAGE, "\n";
		    print STDERR "You have to choose either mosaik or bwa to perform alignments or specify which aligner (-aligner 'mosaik' or 'bwa') was used if you want to only run programs after alignment.", "\n\n";
		    exit;
		}
	    }
	}
	if ($parameterSetSwitch eq 1) { #No need to set parameter more than once
	    last;
	}
    }	
##Parameter set
    if (defined($scriptParameter{$parameterName})) {

	print STDOUT "Set ".$parameterName." to: ".$scriptParameter{$parameterName}, "\n";
    }
    return;
}

sub CreateFileEndings {
###Creates the fileEndings depending on which modules are used by the user to relevant chain. 
    
    my %tempFileEnding; #Used to enable seqential build-up of fileEndings between modules
    
    foreach my $orderParameterElement (@orderParameters) {
	
	if (defined($scriptParameter{$orderParameterElement})) { #Only active parameters

	    if ( ($orderParameterElement =~ /^p[A-Z]/) && ($parameter{$orderParameterElement}{'associatedProgram'}) ) { #Only process programs

		if ($parameter{$orderParameterElement}{'chain'} eq "MAIN") { #MAIN chain
		    
		    if ($parameter{$orderParameterElement}{'fileEnding'} ne "nofileEnding") { #FileEnding exist
			
###MAIN/Per sampleID
			for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {
			    
			    if ($scriptParameter{$orderParameterElement} > 0) { #Fileending should be added    
				
				if ($orderParameterElement eq "pPicardToolsMergeSamFiles") { #Special case
				    
				    if ( (@picardToolsMergeSamFilesPrevious) || (scalar( @{ $infilesLaneNoEnding{ $sampleIDs[$sampleIDCounter] } }) > 1) ) { #Sanity check that we have something to merge and hence to fileEnding should be added
					
					$sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{'pPicardToolsMergeSamFiles'}{'fileEnding'} = $tempFileEnding{$sampleIDs[$sampleIDCounter]}.$parameter{$orderParameterElement}{'fileEnding'}; #Adds from previous entry 
				    }
				    else {
					
					$sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{'pPicardToolsMergeSamFiles'}{'fileEnding'} = $tempFileEnding{$sampleIDs[$sampleIDCounter]}."";
				    }
				}
				else {
				    if (defined($tempFileEnding{$sampleIDs[$sampleIDCounter]})) {
					
					$sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{$orderParameterElement}{'fileEnding'} = $tempFileEnding{$sampleIDs[$sampleIDCounter]}.$parameter{$orderParameterElement}{'fileEnding'};
				    }
				    else  { #First module that should add filending
					$sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{$orderParameterElement}{'fileEnding'} = $parameter{$orderParameterElement}{'fileEnding'};
				    } 
				}
			    }
			    else { #Do not add new module fileEnding
				$sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{$orderParameterElement}{'fileEnding'} = $tempFileEnding{$sampleIDs[$sampleIDCounter]};
			    }
			    $tempFileEnding{$sampleIDs[$sampleIDCounter]} = $sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{$orderParameterElement}{'fileEnding'}; #To enable sequential build-up of fileending
			}
			
###MAIN/Per familyID
			if ($scriptParameter{$orderParameterElement} > 0) { #Fileending should be added 

			    if ($orderParameterElement eq "pPicardToolsMergeSamFiles") { #Special case - do nothing
			    }
			    elsif ( ($orderParameterElement eq "pPicardToolsSortSam") && ($scriptParameter{'analysisType'} eq "rapid") ) { #Special case - do nothing
			    }
			    elsif ( ($orderParameterElement eq "pPicardToolsMergeRapidReads") && ($scriptParameter{'analysisType'} ne "rapid") ) { #Special case - do nothing
			    }
			    else {
				
				if (defined($tempFileEnding{$scriptParameter{'familyID'}})) {
				    $sampleInfo{ $scriptParameter{'familyID'} }{$orderParameterElement}{'fileEnding'} = $tempFileEnding{$scriptParameter{'familyID'}}.$parameter{$orderParameterElement}{'fileEnding'};
				}
				else  { #First module that should add filending
				    $sampleInfo{ $scriptParameter{'familyID'} }{$orderParameterElement}{'fileEnding'} = $parameter{$orderParameterElement}{'fileEnding'};
				}
				$tempFileEnding{$scriptParameter{'familyID'}} = $sampleInfo{ $scriptParameter{'familyID'} }{$orderParameterElement}{'fileEnding'}; #To enable sequential build-up of fileending 
			    }		
			}
		    }
		}
		if ($parameter{$orderParameterElement}{'chain'} ne "MAIN") { #Other chain(s)
		    
		    my $chainfork = $parameter{$orderParameterElement}{'chain'}; 

		    if ($parameter{$orderParameterElement}{'fileEnding'} ne "nofileEnding") { #FileEnding exist
			
###OTHER/Per sampleID
			for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {
			    
			    if ($scriptParameter{$orderParameterElement} > 0) { #Fileending should be added    
			
				unless (defined($tempFileEnding{$chainfork}{$sampleIDs[$sampleIDCounter]})) {	
				    $tempFileEnding{$chainfork}{$sampleIDs[$sampleIDCounter]} = $tempFileEnding{$sampleIDs[$sampleIDCounter]}; #Inherit current MAIN chain. 
				}
				if (defined($tempFileEnding{$chainfork}{$sampleIDs[$sampleIDCounter]})) {
				    $sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{$orderParameterElement}{'fileEnding'} = $tempFileEnding{$chainfork}{$sampleIDs[$sampleIDCounter]}.$parameter{$orderParameterElement}{'fileEnding'};
				}
				else  { #First module that should add filending
				    $sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{$orderParameterElement}{'fileEnding'} = $parameter{$orderParameterElement}{'fileEnding'};
				} 
			    }
			    else { #Do not add new module fileEnding
				$sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{$orderParameterElement}{'fileEnding'} = $tempFileEnding{$chainfork}{$sampleIDs[$sampleIDCounter]};
			    }
			    #NOTE: No sequential build-up of fileending
			}
###Other/Per familyID

			if ($scriptParameter{$orderParameterElement} > 0) { #Fileending should be added
			    
			    unless (defined($tempFileEnding{$chainfork}{$scriptParameter{'familyID'}})) {	
				$tempFileEnding{$chainfork}{$scriptParameter{'familyID'}} =  $tempFileEnding{$scriptParameter{'familyID'}}; #Inherit current MAIN chain. 
			    }
			    if (defined($tempFileEnding{$chainfork}{$scriptParameter{'familyID'}})) {
				$sampleInfo{ $scriptParameter{'familyID'} }{$orderParameterElement}{'fileEnding'} = $tempFileEnding{$chainfork}{$scriptParameter{'familyID'}}.$parameter{$orderParameterElement}{'fileEnding'};
			    }
			    else  { #First module that should add filending
				$sampleInfo{ $scriptParameter{'familyID'} }{$orderParameterElement}{'fileEnding'} = $parameter{$orderParameterElement}{'fileEnding'};
			    }
			    $tempFileEnding{$chainfork}{$scriptParameter{'familyID'}} = $sampleInfo{ $scriptParameter{'familyID'} }{$orderParameterElement}{'fileEnding'}; #To enable sequential build-up of fileending 
			}
		    }
		}
	    }
	}
    }
}

sub ProgramPreRequisites {
###Creates program directories (info & programData & programScript), program script filenames and creates sbatch header.

    my $directoryID = $_[0]; #$samplID|$familyID
    my $programName = $_[1]; #Assigns filename to sbatch script
    my $programDirectory = $_[2]; #Builds from $directoryID/$aligner
    my $callType = $_[3]; ##SNV,INDEL or BOTH
    my $fileHandle = $_[4]; #Program filehandle
    my $nrofCores = $_[5]; #The number of cores to allocate
    my $processTime = $_[6]; #Hours 

    my $fileNamePath;
    my $dryRunFilenamePath;
    my $programDataDirectory;
    my $fileInfoPath;
    my $dryRunFileInfoPath;
###Sbatch script names and directory creation
    
    $programDataDirectory = $scriptParameter{'outDataDir'}."/".$directoryID."/".$programDirectory;
    $fileNamePath = $scriptParameter{'outScriptDir'}."/".$directoryID."/".$programDirectory."/".$programName."_".$directoryID;
    $dryRunFilenamePath = $scriptParameter{'outScriptDir'}."/".$directoryID."/".$programDirectory."/dry_run_".$programName."_".$directoryID;
    $fileInfoPath = $scriptParameter{'outDataDir'}."/".$directoryID."/".$programDirectory."/info/".$programName."_".$directoryID;
    $dryRunFileInfoPath = $scriptParameter{'outDataDir'}."/".$directoryID."/".$programDirectory."/info/dry_run_".$programName."_".$directoryID;
    
    if ($callType ne 0) {
	$fileNamePath .= "_".$callType.".";
	$dryRunFilenamePath .= "_".$callType.".";
	$fileInfoPath .= "_".$callType.".";
	$dryRunFileInfoPath .= "_".$callType.".";
    }
    else {
	$fileNamePath .= ".";
	$dryRunFilenamePath .= ".";
	$fileInfoPath .= ".";
	$dryRunFileInfoPath .= ".";
    }
    	
    `mkdir -p $scriptParameter{'outDataDir'}/$directoryID/$programDirectory/info;`; #Creates the aligner folder and info data file directory
    `mkdir -p $programDataDirectory`; #Creates the aligner folder and if supplied the program data file directory
    `mkdir -p $scriptParameter{'outScriptDir'}/$directoryID/$programDirectory`; #Creates the aligner folder script file directory

    if ( ($scriptParameter{"p".$programName} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	$fileName = $fileNamePath; 
    }
    elsif ($scriptParameter{"p".$programName} == 2) { #Dry run single program
	$fileName = $dryRunFilenamePath; 
	print STDOUT "Dry Run:\n";print MIPLOG  "Dry Run:\n";
    }
    else { #Dry run
	$fileName = $dryRunFilenamePath;
	print STDOUT "Dry Run:\n";print MIPLOG  "Dry Run:\n";
    }

    &Checkfnexists(\$fileName, \$fnend, \$fileNameTracker);

###Info and Log
    print STDOUT "Creating sbatch script for ".$programName." and writing script file(s) to: ".$fileName, "\n";print MIPLOG "Creating sbatch script for ".$programName." and writing script file(s) to: ".$fileName, "\n";

    if ($programName eq "RankVariants") { #Special case

	for (my $ImportantDbFileOutFilesCounter=0;$ImportantDbFileOutFilesCounter<scalar(@ImportantDbFileOutFiles);$ImportantDbFileOutFilesCounter++) {
	    
	    my ($volume,$directories,$file) = File::Spec->splitpath($ImportantDbFileOutFiles[$ImportantDbFileOutFilesCounter]);
	    `mkdir -p $directories;`; 
	    print STDOUT "RankVariants data files will be written to: ".$directories.$directoryID."_ranked_".$callType.".txt", "\n";print MIPLOG "RankVariants data files will be written to: ".$directories.$directoryID."_ranked_".$callType.".txt", "\n";    
	}
    }
    else {
	print STDOUT "Sbatch script ".$programName." data files will be written to: ".$programDataDirectory, "\n";print MIPLOG "Sbatch script ".$programName." data files will be written to: ".$programDataDirectory, "\n";
    }

###Sbatch header
    open ($fileHandle, ">".$fileName) or die "Can't write to ".$fileName.":".$!, "\n";
    
    print $fileHandle "#! /bin/bash -l", "\n";
    print $fileHandle "#SBATCH -A ".$scriptParameter{'projectID'}, "\n";
    print $fileHandle "#SBATCH -n ".$nrofCores, "\n";
    print $fileHandle "#SBATCH -t ".$processTime.":00:00", "\n";
    print $fileHandle "#SBATCH -J ".$programName."_".$directoryID."_".$callType, "\n";
    
    if ( ($scriptParameter{"p".$programName} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	print $fileHandle "#SBATCH -e ".$fileInfoPath.$fileNameTracker.".stderr.txt", "\n";
	print $fileHandle "#SBATCH -o ".$fileInfoPath.$fileNameTracker.".stdout.txt", "\n";
    }
    elsif ($scriptParameter{'pSampleCheck'} == 2) { #Single program dry run
	print $fileHandle "#SBATCH -e ".$dryRunFileInfoPath.$fileNameTracker.".stderr.txt", "\n";
	print $fileHandle "#SBATCH -o ".$dryRunFileInfoPath.$fileNameTracker.".stdout.txt", "\n";
    }
    else { #Dry run
	print $fileHandle "#SBATCH -e ".$dryRunFileInfoPath.$fileNameTracker.".stderr.txt", "\n";
	print $fileHandle "#SBATCH -o ".$dryRunFileInfoPath.$fileNameTracker.".stdout.txt", "\n";
    }
    
    unless ($scriptParameter{'email'} eq 0) {
	print $fileHandle "#SBATCH --mail-type=END", "\n";
	print $fileHandle "#SBATCH --mail-type=FAIL", "\n";
	print $fileHandle "#SBATCH --mail-user=".$scriptParameter{'email'}, "\n\n";	
    }
    
    print $fileHandle 'echo "Running on: $(hostname)"',"\n\n";
    return ($fileInfoPath.$fileNameTracker.".stdout.txt"); #Return stdout for QC check later
}

sub CheckIfMergedFiles {
###Check if any files for this sampleID were merged previously to set infile and PicardToolsMergeSwitch
    my $sampleID = $_[0];

    my $infile;
    my $mergeLanes; #To pick up merged lanes later 
    my $PicardToolsMergeSwitch = 0;

    if ($sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'picardToolsMergeSamFilesPrevious'} == 1) { # Files merged this round with merged file from previous round
	
	for (my $mergeFileCounter=0;$mergeFileCounter<scalar(@picardToolsMergeSamFilesPrevious);$mergeFileCounter++) {
	    
	    if ($picardToolsMergeSamFilesPrevious[$mergeFileCounter] =~ /lane(\d+)|s_(\d+)/) { #Look for lanes_ or lane\d in previously generated file to be merged with current run to be able to extract previous lanes
		
		if($1) {$mergeLanes = $1;} 
		else {$mergeLanes = $2;} #Make sure to always supply lanes from previous regexp  
		$infile = $sampleID."_lanes_".$mergeLanes;

		for (my $laneCounter=0;$laneCounter<scalar(@ { $lane{$sampleID} });$laneCounter++) {
		    $infile .= $lane{$sampleID}[$laneCounter]; #Extract lanes per sampleID
		}
		$PicardToolsMergeSwitch = 1;
	    }
	}
    }
    elsif ( ($scriptParameter{'pPicardToolsMergeSamFiles'} > 0) && (scalar( @{ $infilesLaneNoEnding{$sampleID} }) > 1) ) { #but only if there is more than one mosaikBuild/BWA_Aln file per sample ID (Sanity check)
	$infile = $sampleID."_lanes_";

	for (my $laneCounter=0;$laneCounter<scalar(@ { $lane{$sampleID} });$laneCounter++) {
	    $infile .= $lane{$sampleID}[$laneCounter]; #Extract lanes per sampleID
	}
	$PicardToolsMergeSwitch = 1;
    }
    else {
	$PicardToolsMergeSwitch = 0;
    }
    return ($infile, $PicardToolsMergeSwitch);
}

sub SampleInfoQC {
###Adds outDirectory and outFile to sampleInfo to track all files that QC metrics are to be extracted from later

    my $familyID = $_[0];
    my $sampleID = $_[1]; #SampleID or "noSampleID" for family level data
    my $programName = $_[2];
    my $infile = $_[3]; #infile or "noInFile for family level data"
    my $outDirectory = $_[4];
    my $outFileEnding = $_[5]; #Actually complete outfile for "static" & "infoDirectory" 
    my $outDataType = $_[6];

    if ($sampleID eq "noSampleID") {

	$sampleInfo{$familyID}{'program'}{$programName}{'OutDirectory'} = $outDirectory; #OutDirectory of QC File                                                            
	if ($outDataType eq "static") { #programs which add a static file in its own directory                                                                                                 

	    $sampleInfo{$familyID}{'program'}{$programName}{'OutFile'} = $outFileEnding; #Static QC outFile                                                                     
	}
	if ($outDataType eq "infoDirectory") { #QC metrics are sent to info files                                                                                                                   
	    $sampleInfo{$familyID}{'program'}{$programName}{'OutFile'} = $outFileEnding; #info stdout file                                                                      
	}
	if ($outDataType eq "infileDependent") { #Programs which Add a filending to infile                                                                                                          
	    $sampleInfo{$familyID}{'program'}{$programName}{'OutFile'} = $outFileEnding; #Infile dependent QC outFile                                                                                                                                                                                       
	}

    }
    else {
	
	$sampleInfo{$familyID}{$sampleID}{'program'}{$programName}{$infile}{'OutDirectory'} = $outDirectory; #OutDirectory of QC File                                                              

	if ($outDataType eq "static") { #programs which add a static file in its own directory 
	    
	    $sampleInfo{$familyID}{$sampleID}{'program'}{$programName}{$infile}{'OutFile'} = $outFileEnding; #Static QC outFile
	}
	if ($outDataType eq "infoDirectory") { #QC metrics are sent to info files
	    
	    $sampleInfo{$familyID}{$sampleID}{'program'}{$programName}{$infile}{'OutFile'} = $outFileEnding; #info stdout file
	}
	if ($outDataType eq "infileDependent") { #Programs which Add a filending to infile
	    
	    $sampleInfo{$familyID}{$sampleID}{'program'}{$programName}{$infile}{'OutFile'} = $infile.$outFileEnding; #Infile dependent QC outFile                                                                      
	}
    }
}

sub GATKTargetListFlag {
###Detects if there are different capture kits across sampleIDs. Creates a temporary merged interval_list for all interval_list that have been supplied and returns temporary list. Will also extract specific contigs if requested and return that list if enabled.
    
    my $FILEHANDLE = $_[0];
    my $contigRef = $_[1];

    my %GATKTargetPaddedBedIntervalListTracker;
    my @GATKTargetPaddedBedIntervalListFiles;

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) { #Collect infiles for all sampleIDs
	
	if (defined($scriptParameter{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{'GATKTargetPaddedBedIntervalLists'})) {

	    $scriptParameter{'GATKTargetPaddedBedIntervalLists'} = $scriptParameter{'referencesDir'}."/".$scriptParameter{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{'GATKTargetPaddedBedIntervalLists'};
       
	    $GATKTargetPaddedBedIntervalListTracker{ $scriptParameter{'GATKTargetPaddedBedIntervalLists'} }++;
	    
	    if ($GATKTargetPaddedBedIntervalListTracker{ $scriptParameter{'GATKTargetPaddedBedIntervalLists'} } == 1) { #Not detected previously
		
		push(@GATKTargetPaddedBedIntervalListFiles, $scriptParameter{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{'GATKTargetPaddedBedIntervalLists'});
	    }
	}
    }
    
    ##Determine file to print to module (untouched/merged and/or splited)
    my $outDirectory = $scriptParameter{'GATKTempDirectory'}.'$SLURM_JOB_ID'; #For merged and/or splitet

    if (scalar(@GATKTargetPaddedBedIntervalListFiles) > 1) { #Merge files
      
	print $FILEHANDLE "\n#Generate merged interval_list\n\n"; 
	print $FILEHANDLE "java -Xmx2g -jar ".$scriptParameter{'picardToolsPath'}."/IntervalListTools.jar ";
	print $FILEHANDLE "UNIQUE=TRUE "; #merge overlapping and adjacent intervals to create a list of unique intervals
    
	for (my $fileCounter=0;$fileCounter<scalar(@GATKTargetPaddedBedIntervalListFiles);$fileCounter++) {
	
	    print $FILEHANDLE "INPUT=".$scriptParameter{'referencesDir'}."/".$GATKTargetPaddedBedIntervalListFiles[$fileCounter]." ";
	}
	print $FILEHANDLE "OUTPUT=".$outDirectory."/merged.interval_list", "\n\n"; #Merged outfile

	if (defined($$contigRef)) {
	    
	    my $inDirectory = $scriptParameter{'GATKTempDirectory'}.'$SLURM_JOB_ID';
	    return &SplitTargetFile(*$FILEHANDLE, \$inDirectory, \$outDirectory, \("merged.interval_list"), \$$contigRef);
	}
        return $outDirectory."/merged.interval_list"; #No split
    }
    elsif (defined($$contigRef)) { #Supply original file but create splitted temp file

	return &SplitTargetFile(*$FILEHANDLE, \$scriptParameter{'referencesDir'}, \$outDirectory, \$GATKTargetPaddedBedIntervalListFiles[0], \$$contigRef); #Only 1 file for all samples
    }
    else {#No merge and no split. return original and only file
	return  $scriptParameter{'referencesDir'}."/".$GATKTargetPaddedBedIntervalListFiles[0];
    }

}

sub SplitTargetFile {
##Splits a target file into new contig specific target file

    my $FILEHANDLE = $_[0];
    my $inDirectoryRef = $_[1];
    my $outDirectoryRef = $_[2];
    my $infileRef = $_[3];
    my $contigRef = $_[4]; #The contig to extract
    
    if (defined($$contigRef)) {
	
	print $FILEHANDLE "\n#Generate contig specific interval_list\n\n"; 
	print $FILEHANDLE q?perl -nae 'if($_=~/^\@/) {print $_;} elsif($_=~/^?.$$contigRef.q?\s+/) {print $_;}' ?; #Select header and contig
	print $FILEHANDLE $$inDirectoryRef."/".$$infileRef." "; #Infile
	print $FILEHANDLE "> ".$$outDirectoryRef."/".$$contigRef."_".$$infileRef, "\n\n";#Extract genomic interval info
	
	return $$outDirectoryRef."/".$$contigRef."_".$$infileRef;
    }
}

sub GATKPedigreeFlag {
###Check if "--pedigree" and "--pedigreeValidationType" should be included in analysis
    
    my $FILEHANDLE = $_[0];
    my $outFamilyFileDirectory = $_[1];
    my $pedigreeValidationType = $_[2];

    if (scalar(@sampleIDs) > 2) {

	my $famFile = $outFamilyFileDirectory."/".$scriptParameter{'familyID'}.".fam";
	my $parentCounter;
	my $pqParentCounter = q?perl -ne 'my $parentCounter=0; while (<>) { my @line = split(/\t/, $_); unless ($_=~/^#/) { if ( ($line[2] eq 0) || ($line[3] eq 0) ) { $parentCounter++} } } print $parentCounter; last;'?;
	my $childCounter;
	my $pqChildCounter = q?perl -ne 'my $childCounter=0; while (<>) { my @line = split(/\t/, $_); unless ($_=~/^#/) { if ( ($line[2] ne 0) || ($line[3] ne 0) ) { $childCounter++} } } print $childCounter; last;'?;

	if ($FILEHANDLE eq "*main::GATK_PHTR") { #Special case - GATK PhaseByTransmission needs parent/child or trio 
	    
	    if (scalar(@sampleIDs) < 4) { #i.e.2-3 individuals in pedigree
		    
		$parentCounter = `$pqParentCounter $famFile`;
		$childCounter = `$pqChildCounter $famFile`;		    
		
		if ( ($childCounter == 1) && ($parentCounter > 0) ) { #parent/child or trio
		    print $FILEHANDLE "--pedigreeValidationType ".$pedigreeValidationType." --pedigree ".$outFamilyFileDirectory."/".$scriptParameter{'familyID'}.".fam "; #Pedigree files for samples
		}
		else {
		    $scriptParameter{'pGATKPhaseByTransmission'} = 0; #Override input since pedigree is not valid for analysis
		    print STDERR "Switched GATK PhaseByTransmission to no run mode since MIP did not detect a valid pedigree for this type of analysis. ";print MIPLOG "Switched GATK PhaseByTransmission to no run mode since MIP did not detect a valid pedigree for this type of analysis. ";
		    if ($scriptParameter{'pGATKReadBackedPhasing'} > 0) { #Broadcast
			print STDERR "MIP will still try to run GATK ReadBackedPhasing, but with the '-respectPhaseInInput' flag set to false\n";print MIPLOG "MIP will still try to run GATK ReadBackedPhasing, but with the '-respectPhaseInInput' flag set to false\n";
		    }
		    print STDERR "\n";
		}
	    }
	    else {
		$scriptParameter{'pGATKPhaseByTransmission'} = 0; #Override input since pedigree is not valid for analysis
		print STDERR "Switched GATK PhaseByTransmission to no run mode since MIP did not detect a valid pedigree for this type of analysis. ";print MIPLOG "Switched GATK PhaseByTransmission to no run mode since MIP did not detect a valid pedigree for this type of analysis. ";
		if ($scriptParameter{'pGATKReadBackedPhasing'} > 0) { #Broadcast
		    print STDERR "MIP will still try to run GATK ReadBackedPhasing, but with the '-respectPhaseInInput' flag set to false";print MIPLOG "MIP will still try to run GATK ReadBackedPhasing, but with the '-respectPhaseInInput' flag set to false\n";
		}
		print STDERR "\n";
	    }
	}
	else {

	    $parentCounter = `$pqParentCounter $famFile`;
	    
	    if ($parentCounter > 0) { #Parent
		print $FILEHANDLE "--pedigreeValidationType ".$pedigreeValidationType." --pedigree ".$outFamilyFileDirectory."/".$scriptParameter{'familyID'}.".fam "; #Pedigree files for samples		
	    }		
	}
    }
    return;
}

sub WriteCMDMipLog {
    
    open (MIPLOG, ">>".$mipLogName) or die "Can't write to ".$mipLogName.":".$!, "\n"; #Open file run log
    
    foreach my $orderParameterElement (@orderParameters) {
	
	if (defined($scriptParameter{$orderParameterElement}) ) {

	    if ( ($orderParameterElement eq "configFile") && ($scriptParameter{'configFile'} eq 0) ) { #Do not print
	    }
	    else {
		print MIPLOG "-".$orderParameterElement." ".$scriptParameter{$orderParameterElement}." ";
	    }
	}
    }
    print MIPLOG "\n\n";

    #Note FileHandle MIPLOG not closed
    return;
}

sub WriteYAML {
###Writes a YAML hash to file. 
###Note: 2nd argument should be a hash reference

    my $yamlFile = $_[0]; #Filename
    my $yamlHashRef = $_[1]; #Hash reference to write to file

    open (YAML, ">". $yamlFile) or die "can't open ".$yamlFile.":".$!, "\n";
    print YAML Dump( $yamlHashRef ), "\n";
    close(YAML);
    print STDOUT "Wrote: ".$yamlFile, "\n";
}

sub LoadYAML {
###Loads a YAML file into an arbitrary hash and returns it. Note: Currently only supports hashreferences and hashes and no mixed entries 

    my $yamlFile = $_[0];
    my %yamlHash;

    my $fileType = &DetectYamlContentType($yamlFile);

    open (YAML, "<". $yamlFile) or die "can't open ".$yamlFile.":".$!, "\n";    
        
        if ($fileType eq "reference") {
        %yamlHash = %{ YAML::LoadFile($yamlFile) }; #Load hashreference as hash
        }
        if ($fileType eq "hash") {
        %yamlHash = YAML::LoadFile($yamlFile); #File contained a hash = no workup
        }
    close(YAML);
    
    print STDOUT "Read Yaml file: ". $yamlFile, "\n";
    return %yamlHash;
}

sub DetectYamlContentType {
###Check the content of the YAML file for seperating hashreferences and hash. Return the content type.

    my $yamlFile = $_[0];
    my $fileType;

    open (YAML, "<". $yamlFile) or die "can't open ".$yamlFile.":".$!, "\n";
        
        while (<YAML>) {

            if ($. == 1 && $_=~/^---$/) { #YAML file contains a hashreference
                $fileType = "reference";
                last;
            }
            else {
                $fileType = "hash";
                last;
            }
        }
    close(YAML);
    return $fileType;
}

sub CheckUniqueArrayElement {
###Detects if there are elements in arrayQueryRef that are not present in scalarQueryRef or arrayToCheckRef. If unique adds the unique element to arrayToCheckRef

    my $arrayToCheckRef = $_[0]; #the arrayref to be queried
    my $arrayQueryRef;
    my $scalarQueryRef;

    if (ref($_[1]) eq "ARRAY") {
	$arrayQueryRef = $_[1];
	
	##For each arrayQueryRef element, loop through corresponding arrayToCheckRef element(s), add if there are none or an updated/unique entry.
	for (my $elementsInfoCounter=0;$elementsInfoCounter<scalar(@{$arrayQueryRef});$elementsInfoCounter++) { #all element(s)
	    
	    my $elementFound = 0; #Track if there element is present in arrayToCheckRef
	    
	    for (my $elementsCounter=0;$elementsCounter<scalar( @{$arrayToCheckRef});$elementsCounter++) { #all arrayToCheckRef elements
		
		if (${$arrayToCheckRef}[$elementsCounter] eq ${$arrayQueryRef}[$elementsInfoCounter]) { #Check presens
		    
		    $elementFound = 1;  #Entry is present in both arrays
		}
	    }
	    if ($elementFound == 0) { #Not seen in arrayToCheckRef
		push( @{$arrayToCheckRef}, ${$arrayQueryRef}[$elementsInfoCounter]); #Go ahead and add
	    }
	}
    }
    if (ref($_[1]) eq "SCALAR") {

	$scalarQueryRef = $_[1]; 

	my $elementFound = 0; #Track if there element is present in arrayToCheckRef
	
	for (my $elementsCounter=0;$elementsCounter<scalar( @{$arrayToCheckRef});$elementsCounter++) { #all arrayToCheckRef elements
	    
	    if (${$arrayToCheckRef}[$elementsCounter] eq $$scalarQueryRef) { #Check presens
		
		$elementFound = 1;  #Entry is present in both arrays
	    }
	}
	if ($elementFound == 0) { #Not seen in arrayToCheckRef
	    push( @{$arrayToCheckRef}, $$scalarQueryRef); #Go ahead and add
	}
    }
    return;
}

sub DetermineNrofRapidNodes {
##Determines the number of nodes to allocate depending on the sequence length, which affects the infile size.
    
    my $seqLength = $_[0]; 
    my $infileSize = $_[1];
    
    my $numberNodes = 0; #Nodes to allocate
    my $readPositionWeight = 1; #Scales the readStart and readStop position
    my $ReadNrofLines;    

    if ($seqLength > 75 && $seqLength <= 101) {
	$ReadNrofLines = 190000000; #Read batch size
	$numberNodes = floor($infileSize / (12 * $ReadNrofLines) ); #Determines the number of nodes to use, 150000000 ~ 37,5 million reads, 13 = 2 sdtdev from sample population - currently poor estimate with compression confunding calculation.	
	print STDOUT "Number of Nodes: ".$numberNodes, "\n";
    }
    if ($seqLength > 50 && $seqLength <= 75) {
	$ReadNrofLines = 190000000; #Read batch size
	$numberNodes = floor($infileSize / (9.75 * $ReadNrofLines) ); #Determines the number of nodes to use, 150000000 ~ 37,5 million reads, 13 = 2 sdtdev from sample population - currently poor estimate with compression confunding calculation.	
	print STDOUT "Number of Nodes: ".$numberNodes, "\n";
    }
    if ($seqLength >= 50 && $seqLength < 75) {
	$ReadNrofLines = 130000000; #Read batch size
	$numberNodes = floor($infileSize / (7 * $ReadNrofLines) ); #Determines the number of nodes to use, 150000000 ~ 37,5 million reads, 13 = 2 sdtdev from sample population - currently poor estimate with compression confunding calculation.
	print STDOUT "Number of Nodes: ".$numberNodes, "\n";
    }
    if ($seqLength >= 35 && $seqLength < 50) {
	$ReadNrofLines = 95000000; #Read batch size
	$numberNodes = floor($infileSize / (6 * $ReadNrofLines) ); #Determines the number of nodes to use, 150000000 ~ 37,5 million reads, 13 = 2 sdtdev from sample population - currently poor estimate with compression confunding calculation.
	print STDOUT "Number of Nodes: ".$numberNodes, "\n";
    }
    if ($numberNodes <= 1) {
	
	$numberNodes = 2; #Ensure that at least 1 readbatch is processed
    }
    return $numberNodes, $ReadNrofLines;
}

sub CheckUniqueIDNs {
##Test that the familyID and the sampleID(s) exists and are unique.

    my %seen; #Hash to test duplicate sampleIDs later

    if (scalar(@sampleIDs) == 0) {

	print STDOUT "\nPlease provide sampleID(s)\n\n";
	exit;
    }

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {

	$seen{$sampleIDs[$sampleIDCounter]}++; #Increment instance to check duplictaes later
	
	if ($scriptParameter{'familyID'} eq $sampleIDs[$sampleIDCounter]) {
	    
	    print STDOUT "\nFamilyID: ".$scriptParameter{'familyID'}." equals sampleID: ".$sampleIDs[$sampleIDCounter].". Please make sure that the familyID and sampleID(s) are unique.\n";
	    exit;
	}
	if ($seen{$sampleIDs[$sampleIDCounter]} > 1) {
	
	    print STDOUT "\nSampleID: ".$sampleIDs[$sampleIDCounter]." is not uniqe.\n\n";
	    exit;
	}
    }
    return;
}

sub UpdateYAML {
##Updates the config file to particular user/cluster for entries following specifications. Leaves other entries untouched.

    my $orderParameterElement = $_[0]; #Parameter to update
    my $clusterConstantPath = $_[1]; #Set the project specific path for this cluster
    my $analysisConstantPath = $_[2]; #Set the project specific path for this cluster
    my $analysisType =  $_[3]; #Sets the analysis run type e.g., "exomes", "genomes", "rapid"
    my $familyID = $_[4]; #Sets the familyID
    my $aligner = $_[5]; #Sets the aligner used
    
    if ($scriptParameter{$orderParameterElement}) {
	
	$scriptParameter{$orderParameterElement} =~ s/CLUSTERCONSTANTPATH!/$clusterConstantPath/gi; #Exchange CLUSTERCONSTANTPATH! for current cluster path
	$scriptParameter{$orderParameterElement} =~ s/ANALYSISCONSTANTPATH!/$analysisConstantPath/gi; #Exchange ANALYSISCONSTANTPATH! for the current analysis path
	$scriptParameter{$orderParameterElement} =~ s/ANALYSISTYPE!/$analysisType/gi; #Exchange ANALYSISTYPE! for the current analysis type
	$scriptParameter{$orderParameterElement} =~ s/FDN!/$familyID/gi; #Exchange FND! for the current familyID
	$scriptParameter{$orderParameterElement} =~ s/ALIGNER!/$aligner/gi; #Exchange ALIGNER! for the current aligner
    }
}

sub CheckAutoBuild {
##Checks if autobuild is on and returns "1" if enabled or "0" if not)

    my $parameterNameRef = $_[0];
    my $sampleIDRef = $_[1];

    if (defined($sampleIDRef)) {
	
	if ( ($parameter{ $scriptParameter{'familyID'} }{$$sampleIDRef}{$$parameterNameRef}{'buildFile'} eq "yesAutoBuild") || ($parameter{ $scriptParameter{'familyID'} }{$$sampleIDRef}{$$parameterNameRef}{'buildFile'} eq 1) ) {
	    
	    return "1"; #Flag that autobuild is needed
	}
	else {
	    return "0"; #No autobuild is needed   
	}
    }
    else {
	
	if ( ($parameter{$$parameterNameRef}{'buildFile'} eq "yesAutoBuild") || ($parameter{$$parameterNameRef}{'buildFile'} eq 1) ) { #1 for arrays
	    
	    return "1"; #Flag that autobuild is needed
	}
	else {
	    return "0"; #No autobuild is needed   
	}
    }
}

sub ParseHumanGenomeReference {
##Detect the humanGenomeReference: Source (hg19 or GRCh, Version and chromosome prefix (prefix might be removed in the future))
    
    my $humanGenomeReference = $_[0];
    
    my $humanGenomeReferenceVersion; #Version of GenomeBuild
    my $humanGenomeReferenceSource; #Ensembl or NCBI
    my $humanGenomeReferenceNameNoEnding;
    
    if ($humanGenomeReference =~/^Homo_sapiens.GRCh(\d+\.\d+|\d+)/) { #Used to change capture kit genome reference version later
	$humanGenomeReferenceVersion = $1;
	$humanGenomeReferenceSource = "GRCh"; #Ensembl
    }
    elsif ($humanGenomeReference =~/^Homo_sapiens.hg(\d+)/) { #Used to change capture kit genome reference version later
	$humanGenomeReferenceVersion = $1;
	$humanGenomeReferenceSource = "hg"; #Refseq
    }
    else {
	print STDERR "MIP cannot detect what kind of humanGenomeReference you have supplied. If you want to automatically set the capture kits used please supply the refrence on this format: [Species].[Source][Version].", "\n\n";
    }
    ($humanGenomeReferenceNameNoEnding) = &RemoveFileEnding(\$humanGenomeReference, ".fasta");
    ($humanGenomeCompressed) = &CheckGzipped(\$humanGenomeReference);
    return ($humanGenomeReferenceVersion, $humanGenomeReferenceSource, $humanGenomeReferenceNameNoEnding);    
}

sub CheckFileEndingsToBeBuilt {
##Checks files to be built by combining filename stub with fileendings. 
    
    my $fileEndingsRef = $_[0]; #Reference to the fileEndings to be added to the filename stub
    my $parameterName = $_[1]; 
    
    for (my $fileEndingsRefCounter=0;$fileEndingsRefCounter<scalar(@{$fileEndingsRef});$fileEndingsRefCounter++) { #All fileEndings
	
	&CheckExistance(\($scriptParameter{'referencesDir'}."/".$scriptParameter{$parameterName}.${$fileEndingsRef}[$fileEndingsRefCounter]), \$parameterName, "f");
    }
}

sub CheckExistance {
##Checks if a file/directory exists and if autoBuild is on or not. If file/directory does not extis and there is no autobuild, croaks and exists.

    my $itemNameRef = $_[0];
    my $parameterNameRef = $_[1];
    my $itemToCheck = $_[2];    
    my $sampleIDRef = $_[3];
    
    
    if ($itemToCheck eq "d") {

	unless (-d $$itemNameRef) { #Check existence of supplied directory
	    print STDERR $USAGE, "\n";
	    print STDERR "\nCould not find intended ".$$parameterNameRef." directory: ".$$itemNameRef, "\n\n";
	    exit;		
	}
    }
    elsif ($itemToCheck eq "f") {
	
	unless (-f $$itemNameRef) { #Check existence of supplied file in supplied reference dir
	    
	    if (defined($sampleIDRef)) { #Individual files per sampleID
		
		$parameter{ $scriptParameter{'familyID'} }{$$sampleIDRef}{$$parameterNameRef}{'buildFile'} = &CheckAutoBuild(\$$parameterNameRef, \$$sampleIDRef); #Check autoBuild or not and return value
		
		if ($parameter{ $scriptParameter{'familyID'} }{$$sampleIDRef}{$$parameterNameRef}{'buildFile'} == 0) { #No autobuild
		    
		    print STDERR $USAGE, "\n";
		    print STDERR "\nCould not find intended ".$$parameterNameRef." file: ".$$itemNameRef, "\n\n";
		    exit;		
		}
	    }
	    else {
 
		$parameter{$$parameterNameRef}{'buildFile'} = &CheckAutoBuild(\$$parameterNameRef); #Check autoBuild or not and return value
	       
		if ($parameter{$$parameterNameRef}{'buildFile'} == 0) { #No autobuild
		    
		    print STDERR $USAGE, "\n";
		    print STDERR "\nCould not find intended ".$$parameterNameRef." file: ".$$itemNameRef, "\n\n";
		    exit;		
		}
	    }
	}
	else {
	    
	    if (defined($sampleIDRef)) {
		
		$parameter{ $scriptParameter{'familyID'} }{$$sampleIDRef}{$$parameterNameRef}{'buildFile'} = 0; #File exist in this check
	    }
	    else {
		$parameter{$$parameterNameRef}{'buildFile'} =  0; #File exist in this check
	    }
	}
    }
    return;
}

 sub SetAutoBuildFeature {
##Sets parameters with autoBuild enabled to the new value dependent on $referenceName
     
     my $featureName = $_[0];
     my $referenceFileEndingRef = $_[1];
     my $referenceNameRef = $_[2];
     my $printSwitch = $_[3];
     
     if( defined($scriptParameter{$featureName}) && ($scriptParameter{$featureName} eq "notSetYet") ) {

	 $scriptParameter{$featureName} = $$referenceNameRef.$$referenceFileEndingRef;

	 if ( (defined($printSwitch)) && ($printSwitch ne "noPrint") ) {

	     print STDOUT "Set ".$featureName." to: ".$scriptParameter{$featureName}, "\n";
	 }
	 if ($featureName eq "bwaBuildReference") {

	     &CheckFileEndingsToBeBuilt(\@bwaBuildReferenceFileEndings, "bwaBuildReference");
	 }
	 elsif ($featureName eq "mosaikJumpDbStub") {

	     &CheckFileEndingsToBeBuilt(\@mosaikJumpDbStubFileEndings, "mosaikJumpDbStub");
	 }
	 else {#Complete fileName - No stubs
	    
	     &CheckExistance(\($scriptParameter{'referencesDir'}."/".$scriptParameter{$featureName}), \$featureName, "f");
         }
    }
}

sub MoveMosaikNN {
##Locate MOSAIK path and move neural network files in place if lacking

    my @paths = split(/:/,$ENV{PATH});

    for (my $pathsCounter=0;$pathsCounter<scalar(@paths);$pathsCounter++) {

	if ($paths[$pathsCounter] =~/MOSAIK/) {
	    
	   $paths[$pathsCounter] =~ s/bin\//src\/networkFile/g; #Location of NN files

	   print STDOUT "\nCould not find Mosaik Network Files in ".$scriptParameter{'referencesDir'},"\n";
	   print STDOUT "\nCopying Mosaik Network Files ".$scriptParameter{'mosaikAlignNeuralNetworkSeFile'}." and ".$scriptParameter{'mosaikAlignNeuralNetworkPeFile'}." to ".$scriptParameter{'referencesDir'}." from ".$paths[$pathsCounter], "\n\n";
	   `cp $paths[$pathsCounter]/$scriptParameter{'mosaikAlignNeuralNetworkSeFile'} $scriptParameter{'referencesDir'}/`;
	   `cp $paths[$pathsCounter]/$scriptParameter{'mosaikAlignNeuralNetworkPeFile'} $scriptParameter{'referencesDir'}/`;
	   last;
	}
    }
}

sub CheckUserInfoArrays {
##Determine if the user supplied info on array parameter
    
    my $arrayRef = $_[0];
    my $parameterName = $_[1];
    
    my $userSuppliedInfoSwitch;
    
    if (scalar(@{$arrayRef}) == 0) { ##No user supplied sample info
	
	if (defined($scriptParameter{$parameterName})) { #sampleIDs info in config file
	    
	    $userSuppliedInfoSwitch = 1; #No user supplied sample info, but present in config file do NOT overwrite using info from pedigree file
	}
	else { #No sampleIDs info in config file
	    
	    $userSuppliedInfoSwitch = 0; #No user supplied sample info, not defined $scriptParameter{'sampleIDs'} in config file, add it from pedigree file
	}
    }
    else {
	$userSuppliedInfoSwitch = 1; # User supplied sample info, do NOT overwrite using info from pedigree file	
    }
    return $userSuppliedInfoSwitch;
}

sub SetTargetFiles {
    
    
    my $familyIDRef = $_[0];
    my $sampleIDRef = $_[1];
    my $parameterNameRef = $_[2];
    my $referenceFileEndingRef = $_[3];
    
    if (defined($scriptParameter{$$familyIDRef}{$$sampleIDRef}{$$parameterNameRef})) { #Capture kit check
	
	$scriptParameter{$$familyIDRef}{$$sampleIDRef}{$$parameterNameRef} =~ s/GenomeReferenceSource/$humanGenomeReferenceSource/; #Replace with Refseq genome or Ensembl genome
	$scriptParameter{$$familyIDRef}{$$sampleIDRef}{$$parameterNameRef} =~ s/Version/$humanGenomeReferenceVersion/; #Replace with actual version 
	$sampleInfo{$$familyIDRef}{$$sampleIDRef}{$$parameterNameRef} = $scriptParameter{$$familyIDRef}{$$sampleIDRef}{$$parameterNameRef}; #Add to sampleInfo for qc print later
	
	&CheckExistance(\($scriptParameter{'referencesDir'}."/".$scriptParameter{$$familyIDRef}{$$sampleIDRef}{$$parameterNameRef}), \$$parameterNameRef, "f", \$$sampleIDRef);

	print STDOUT "Set ".$$parameterNameRef." to: ".$scriptParameter{$$familyIDRef}{$$sampleIDRef}{$$parameterNameRef}, "\n";
    }
    else {
	
	$supportedCaptureKits{'Latest'} =~ s/GenomeReferenceSource/$humanGenomeReferenceSource/; #Replace with Refseq genome or Ensembl genome
	$supportedCaptureKits{'Latest'} =~ s/Version/$humanGenomeReferenceVersion/; #Replace with actual version
	$scriptParameter{$$parameterNameRef} = "notSetYet"; #Required for autobuild
	&SetAutoBuildFeature($$parameterNameRef, \$$referenceFileEndingRef, \$supportedCaptureKits{'Latest'}, "noPrint"); #Always use the most updated capture kit when building target list
    }
}

sub PrepareArrayParameters {
##Check if user supplied cmd info and supplies arrayParameters to scriptParameters

    my $arrayRef = $_[0];
    my $parameterName = $_[1];
    my $parameterType = $_[2];
    my $parameterDefault = $_[3];
    my $associatedPrograms = $_[4]; #comma separated string
    my $parameterExistsCheck = $_[5]; #Check if intendent file exists in reference directory

    if (scalar(@{$arrayRef}) == 0) { #No input from cmd or from pedigree
	
	$parameter{$parameterName}{'value'} = "nocmdinput"; #To enable use of subroutine &AddToScriptParameter
    }
    else {
	$parameter{$parameterName}{'value'} = "SetbyUser";
	@{$arrayRef} = join(',',@{$arrayRef}); #If user supplied parameter a comma separated list
    }
    push(@orderParameters, $parameterName); #Add to enable later evaluation of parameters in proper order & write to master file
    &AddToScriptParameter($parameterName, $parameter{$parameterName}{'value'}, $parameterType, $parameterDefault, $associatedPrograms, $parameterExistsCheck);
}

sub CheckUniqueTargetFiles {
##Checks target files within parameters for identical entries
    
    my $arrayRef = $_[0]; #Array to loop in for parameter (e.g. sampleID)
    my $countRef = $_[1]; #Offset in array
    my $fileToCompareRef = $_[2]; #The file to compare against rest of array
    my $parameterName = $_[3]; #Parameter to evaluate
    
    for (my $compareCounter=($$countRef + 1);$compareCounter<scalar(@{$arrayRef});$compareCounter++) { #Compare all target files to remove autoBuild if file names are identical for each flag
	if ( (defined($$fileToCompareRef)) && (defined($scriptParameter{ $scriptParameter{'familyID'} }{${$arrayRef}[$compareCounter]}{$parameterName})) ) {

	    if ($$fileToCompareRef eq $scriptParameter{ $scriptParameter{'familyID'} }{${$arrayRef}[$compareCounter]}{$parameterName}) {
	    
		$parameter{ $scriptParameter{'familyID'} }{${$arrayRef}[$compareCounter]}{$parameterName}{'buildFile'} = 0;
	    }
	}
    }
    return;
}


sub CheckGzipped {
##Check if a file is gzipped

    my $fileNameRef = $_[0];

    my $fileCompressionStatus = "unCompressed";

    if ( (defined($$fileNameRef)) && ($$fileNameRef =~/.gz$/) ) {
	
	$fileCompressionStatus = "compressed"; 
    }
    return $fileCompressionStatus;
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

sub ScriptParameterPerSampleID {
##Enables files handled per SampleID to be processed by AddToScriptParameters

    my $familyRef = $_[0];
    my $sampleIDRef = $_[1];
    my $parameterName = $_[2];
    
    if (defined($scriptParameter{$$familyRef}{$$sampleIDRef}{$parameterName})) {
	
	$scriptParameter{$parameterName} = 1; #Define in scriptParameter so that we now that parameter is present per SampleID
    }
}

sub EnableArrayParameter {
##Adds arrayRef to scriptParameters for recreation of cmd in log and seperated input parameter string into array elements
    
    my $arrayRef = $_[0];
    my $parameterNameRef = $_[1];
    
    $scriptParameter{$$parameterNameRef} = join(',',@{$arrayRef}); #Add to enable recreation of cmd line later
    @{$arrayRef} = split(/,/,join(',', @{$arrayRef})); #Enables comma separated list of sample IDs from user supplied cmd info

    return;
}

sub SetTargetandAutoBuild {
##Set autoBuild for target files and calls SetTargetFile sub

    my $arrayRef = $_[0];
    my $parameterNameRef = $_[1];
    my $fileEndingRef = $_[2];
    
    for (my $elementsCounter=0;$elementsCounter<scalar(@{$arrayRef});$elementsCounter++) {
	
	$parameter{ $scriptParameter{'familyID'} }{${$arrayRef}[$elementsCounter]}{$$parameterNameRef}{'buildFile'} = "yesAutoBuild";
	&SetTargetFiles(\$scriptParameter{'familyID'}, \${$arrayRef}[$elementsCounter], \$$parameterNameRef, \$$fileEndingRef);
    }
return;    
}

sub CheckTargetExistFileBed {
##Check that supplied target file ends with ".bed" and exists

    my $fileRef = $_[0];
    my $parameterName = $_[1];

    if ($$fileRef !~/.bed$/) {

	print STDERR "Could not find intendended 'file ending with .bed' for target file: ".$$fileRef." in parameter '-".$parameterName."'", "\n";
	exit;
    }
    unless (-f $scriptParameter{'referencesDir'}."/".$$fileRef) {

	print STDERR "Could not find intendended '.bed' file for target file: ".$scriptParameter{'referencesDir'}."/".$$fileRef." in parameter '-".$parameterName."'", "\n\n";
	exit;
    }
    return;
}

sub CompareArrayElements {
##Compares the number of elements in two arrays 

    my $arrayRef = $_[0]; #Array to match
    my $arrayQueryRef = $_[1]; #Array to be compared
    my $parameterName = $_[2]; #Reference parameter
    my $parameterNameQuery = $_[3]; #Query parameter
    
    if (scalar(@{$arrayRef}) != scalar(@{$arrayQueryRef})) {
	
	print STDERR "\nThe supplied '-".$parameterNameQuery."' lists do not equal the number of elements in '-".$parameterName."'. Please specify a equal number of elements in both lists", "\n\n";
	exit;
    }
    return;
}

sub SetAutoBuildAndScriptParameterPerSample {
##Sets autoBuild and populates scriptParameter hash with array elements per sampleID

    my $sampleIDArrayRef = $_[0];
    my $parameterArrayRef = $_[1];
    my $parameterNameRef = $_[2];

    for (my $sampleIDsCounter=0;$sampleIDsCounter<scalar(@{$sampleIDArrayRef});$sampleIDsCounter++) { #All sampleIDs
	
	$parameter{ $scriptParameter{'familyID'} }{${$sampleIDArrayRef}[$sampleIDsCounter]}{$$parameterNameRef}{'buildFile'} = "yesAutoBuild"; #Turn on autoBuild
	$scriptParameter{ $scriptParameter{'familyID'} }{${$sampleIDArrayRef}[$sampleIDsCounter]}{$$parameterNameRef} = ${$parameterArrayRef}[$sampleIDsCounter]; #Populate hash that is used in modules
    }
    return;
}

sub SetTargetFileGeneralBuildParameter {
##Sets the general build parameters $sampleIDBuildFile and $sampleIDBuildFileNoEnding
    
    my $targetfileRef = $_[0];
    my $parameterName = $_[1];
    my $sampleIDBuildFileRef = $_[2];
    my $sampleIDBuildFileNoEndingRef = $_[3];
    my $sampleIDRef = $_[4];
    
    $$sampleIDBuildFileNoEndingRef = &RemoveFileEnding(\$$targetfileRef, $referenceFileEndings{$parameterName}); #Remove ".fileending" from reference filename
    $parameter{ $scriptParameter{'familyID'} }{$$sampleIDRef}{$parameterName}{'buildFile'} = 0; #Build once then done
    $$sampleIDBuildFileRef = $$targetfileRef;
    
}

sub PrintCheckExistandMoveFile {
##Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.

    my $FILEHANDLE = $_[0]; #Sbatch filehandle to print to
    my $intendedFilePathRef = $_[1]; #File to check
    my $temporaryFilePathRef = $_[2]; #File that has been created

    print $FILEHANDLE "[ -s ".$$intendedFilePathRef." ] "; #Check file exists and is larger than 0
    print $FILEHANDLE "&& rm ".$$temporaryFilePathRef." "; #If other processes already has created file, remove temp file
    print $FILEHANDLE "|| "; #File has not been created by other processes
    print $FILEHANDLE "mv ".$$temporaryFilePathRef." ".$$intendedFilePathRef,"\n\n"; #Move file in place
    
}

sub DefineAnnovarTables {
##Loads annovar tables parameters.
    
    my $annovarGenomeBuildVersion = $scriptParameter{'annovarGenomeBuildVersion'}; #Set the current annovar genome build

    my @annovarTablesGeneAnno = ("refGene", "knownGene", "ensGene"); #Tables using annotation option "geneanno"
    my @annovarTablesRegionAnno = ("mce46way", "gerp++elem", "segdup", "tfbs", "mirna"); #Tables using annotation option "regionanno"
    my @annovarTablesFilter = ("snp137", "snp135", "snp132", "snp131", "snp130", "snp129", "snp137NonFlagged", "snp135NonFlagged", "snp132NonFlagged", "snp131NonFlagged", "snp130NonFlagged", "1000g2012apr_all", "1000g2012apr_amr", "1000g2012apr_eur", "1000g2012apr_asn", "1000g2012apr_afr", "1000g2012feb_all", "esp6500si_all", "esp6500_all", "esp6500_aa", "esp6500_ea", "esp5400_all", "esp5400_aa", "esp5400_ea","clinvar_20131105", "ljb2_sift", "ljb2_pp2hdiv", "ljb2_pp2hvar", "ljb2_mt", "ljb2_ma", "ljb2_fathmm", "ljb2_siphy", "ljb2_lrt", "ljb_all", "ljb2_gerp++", "ljb2_phylop"); #Tables using annotation option "filter"
    my @annovarTablesUrlUcsc = ("mce46way", "segdup", "tfbs", "mirna"); #Tables using urlAlias "ucsc"
    my @annovarGenericFiltering = ("esp6500si_all", "esp6500_all", "esp6500_aa", "esp6500_ea", "esp5400_all", "esp5400_aa", "esp5400_ea","clinvar_20131105"); #Tables using generic option
    my @annovarGenericFiles = ($annovarGenomeBuildVersion."_esp6500si_all.txt", $annovarGenomeBuildVersion."_esp6500_all.txt", $annovarGenomeBuildVersion."_esp6500_aa.txt", $annovarGenomeBuildVersion."_esp6500_ea.txt", $annovarGenomeBuildVersion."_esp5400_all.txt", $annovarGenomeBuildVersion."_esp5400_aa.txt", $annovarGenomeBuildVersion."_esp5400_ea.txt", $annovarGenomeBuildVersion."_clinvar_20131105.txt"); #Generic table files
    my @annovarRefgeneFiles = ($annovarGenomeBuildVersion."_refGene.txt", $annovarGenomeBuildVersion."_refGeneMrna.fa", $annovarGenomeBuildVersion."_refLink.txt"); #Cater for multiple download
    my @annovarKnownGeneFiles = ($annovarGenomeBuildVersion."_knownGene.txt", $annovarGenomeBuildVersion."_kgXref.txt", $annovarGenomeBuildVersion."_knownGeneMrna.fa"); #Cater for multiple download
    my @annovarEnsGeneFiles = ($annovarGenomeBuildVersion."_ensGene.txt", $annovarGenomeBuildVersion."_ensGeneMrna.fa"); #Cater for multiple download

    #Set UCSC alias for download from UCSC
    $annovarTables{'mce46way'}{'ucscAlias'} = "phastConsElements46way";
    $annovarTables{'segdup'}{'ucscAlias'} = "genomicSuperDups";
    $annovarTables{'tfbs'}{'ucscAlias'} = "tfbsConsSites";
    $annovarTables{'mirna'}{'ucscAlias'} = "wgRna";

    #Set GeneAnno files
    push(@{$annovarTables{'refGene'}{'file'}}, @annovarRefgeneFiles);
    push(@{$annovarTables{'knownGene'}{'file'}}, @annovarKnownGeneFiles);
    push(@{$annovarTables{'ensGene'}{'file'}}, @annovarEnsGeneFiles); 

    for (my $tablesCounter=0;$tablesCounter<scalar(@annovarSupportedTableNames);$tablesCounter++) {
	
	&AnnovarTableParameters(\$annovarSupportedTableNames[$tablesCounter], \@annovarSupportedTableNames, "dbtype", $annovarSupportedTableNames[$tablesCounter]);
	&AnnovarTableParameters(\$annovarSupportedTableNames[$tablesCounter], \@annovarSupportedTableNames, "download", $annovarSupportedTableNames[$tablesCounter]);
	$parameter{$annovarSupportedTableNames[$tablesCounter]}{'buildFile'} = "yesAutoBuild";
    }


#Tables using different download call from dbtype call
    $annovarTables{'1000g2012apr_all'}{'download'} = "ALL.sites.2012_04";
    $annovarTables{'1000g2012feb_all'}{'download'} = "ALL.sites.2012_02";
    $annovarTables{'1000g2012apr_afr'}{'download'} = "AFR.sites.2012_04";
    $annovarTables{'1000g2012apr_amr'}{'download'} = "AMR.sites.2012_04";
    $annovarTables{'1000g2012apr_eur'}{'download'} = "EUR.sites.2012_04";
    $annovarTables{'1000g2012apr_asn'}{'download'} = "ASN.sites.2012_04";

#Set 1000G Table filename
    $annovarTables{'1000g2012apr_all'}{'file'}[0] = $annovarGenomeBuildVersion."_ALL.sites.2012_04.txt";
    $annovarTables{'1000g2012feb_all'}{'file'}[0] = $annovarGenomeBuildVersion."_ALL.sites.2012_02.txt";
    $annovarTables{'1000g2012apr_afr'}{'file'}[0] = $annovarGenomeBuildVersion."_AFR.sites.2012_04.txt";
    $annovarTables{'1000g2012apr_amr'}{'file'}[0] = $annovarGenomeBuildVersion."_AMR.sites.2012_04.txt";
    $annovarTables{'1000g2012apr_eur'}{'file'}[0] = $annovarGenomeBuildVersion."_EUR.sites.2012_04.txt";
    $annovarTables{'1000g2012apr_asn'}{'file'}[0] = $annovarGenomeBuildVersion."_ASN.sites.2012_04.txt";

    for (my $tablesCounter=0;$tablesCounter<scalar(@annovarTablesGeneAnno);$tablesCounter++) {

	&AnnovarTableParameters(\$annovarTablesGeneAnno[$tablesCounter], \@annovarTablesGeneAnno, "annotation", "geneanno");
	&AnnovarTableParameters(\$annovarTablesGeneAnno[$tablesCounter], \@annovarTablesGeneAnno, "urlAlias", "annovar");
    }
    for (my $tablesCounter=0;$tablesCounter<scalar(@annovarTablesRegionAnno);$tablesCounter++) {

	&AnnovarTableParameters(\$annovarTablesRegionAnno[$tablesCounter], \@annovarTablesRegionAnno, "annotation", "regionanno");
	&AnnovarTableParameters(\$annovarTablesRegionAnno[$tablesCounter], \@annovarTablesRegionAnno, "urlAlias", "annovar");
	&AnnovarTableParameters(\$annovarTablesRegionAnno[$tablesCounter], \@annovarTablesUrlUcsc, "urlAlias", "ucsc"); #Overwrite for ucsc tables NOTE: not all in RegionAnno
    }
    for (my $tablesCounter=0;$tablesCounter<scalar(@annovarTablesFilter);$tablesCounter++) {

	&AnnovarTableParameters(\$annovarTablesFilter[$tablesCounter], \@annovarTablesFilter, "annotation", "filter");
	&AnnovarTableParameters(\$annovarTablesFilter[$tablesCounter], \@annovarTablesFilter, "urlAlias", "annovar");
	&AnnovarTableParameters(\$annovarTablesFilter[$tablesCounter], \@annovarTablesFilter, "indexFile", ".idx"); 
    }	
    for (my $tablesCounter=0;$tablesCounter<scalar(@annovarGenericFiltering);$tablesCounter++) {
	
	&AnnovarTableParameters(\$annovarGenericFiltering[$tablesCounter], \@annovarGenericFiltering, "dbtype", "generic");
	&AnnovarTableParameters(\$annovarGenericFiltering[$tablesCounter], \@annovarGenericFiltering, "file", $annovarGenericFiles[$tablesCounter]);
    }
}

sub AnnovarTableParameters {
    
    my $tableNameRef = $_[0];
    my $arrayRef = $_[1];
    my $parameterType = $_[2];
    my $parameterValue = $_[3];
    
    for (my $tablesCounter=0;$tablesCounter<scalar(@{$arrayRef});$tablesCounter++) {
		
	if (${$arrayRef}[$tablesCounter] eq $$tableNameRef) {
	
	    if ($parameterType eq "file") { #Add as array instead, since some annovar tables have multiple files downloaded for the same call

		push(@{$annovarTables{$$tableNameRef}{$parameterType}}, $parameterValue);
	    }
	    else {
		
		$annovarTables{$$tableNameRef}{$parameterType} = $parameterValue;
	    }
	    last; #No table should be represented twice within the same array
	}
    }
}

sub CollectSeqContigs {
##Collects sequences contigs used in analysis from human genome sequence dictionnary

    my $pqSeqDict = q?perl -nae 'if($F[0]=~/^\@SQ/) { if($F[1]=~/SN\:(\S+)/) {print $1, ",";} }' ?; 
    my $SeqDictLocation = $scriptParameter{'referencesDir'}."/".$humanGenomeReferenceNameNoEnding.".dict";
    @contigs = `$pqSeqDict $SeqDictLocation `; #returns a comma seperated string of sequence contigs
    @contigs = split(/,/,join(',', @contigs));
}

sub OverWriteConfigParamWithCMDInfo {

    my $parameterName = $_[0];

    if ($parameter{$parameterName}{'value'} ne "nocmdinput") { #Overwrite config with cmd info for regExp entries
	
	$scriptParameter{$parameterName} = $parameter{$parameterName}{'value'};
    }
}

####
#Decommissioned
####
