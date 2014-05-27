#!/usr/bin/perl - w

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
mip.pl  -ifd [inFilesDirs,.,.,.,n] -isd [inScriptDir,.,.,.,n] -rd [refdir] -p [project ID] -s [sample ID,.,.,.,n] -em [e-mail] -osd [outdirscripts] -odd [outDataDir] -f [familyID] -p[program]
               ####MIP
	       -ifd/--inFilesDirs Infile directory(s), comma sep (Mandatory: Supply whole path,)
               -isd/--inScriptDir The pipeline custom script in directory (Mandatory: Supply whole path)
               -rd/--referencesDir Reference(s) directory (Mandatory: Supply whole path)
	       -p/--projectID The project ID  (Mandatory)
	       -s/--sampleIDs The sample ID(s),comma sep (Mandatory)
	       -em/--email E-mail
               -emt/--emailType E-mail type (defaults to F (=FAIL);Options: B (=BEGIN) and/or F (=FAIL) and/or E=(END))
	       -odd/--outDataDir The data files output directory (Mandatory: Supply whole path)
	       -osd/--outScriptDir The script files (.sh) output directory (Mandatory: Supply whole path)
               -f/--familyID Group id of samples to be compared (defaults to "0" (=no), (Ex: 1 for IDN 1-1-1A))
               -pedigree/--pedigreeFile (Supply whole path, defaults to "")
               -huref/--humanGenomeReference Fasta file for the human genome reference (defaults to "Homo_sapiens.GRCh37.d5.fasta;1000G decoy version 5")
               -al/--aligner Setting which aligner was used for alignment in previous analysis (defaults to "")
               -at/--analysisType Type of analysis to perform (defaults to "exomes";Valid entries: "genomes", "exomes", "rapid")
               -mc/--maximumCores The maximum number of cores per node used in the analysis (defaults to "8")
               -c/--configFile YAML config file for script parameters (defaults to "")
               -wc/--writeConfigFile Write YAML configuration file for script parameters (defaults to "";Supply whole path)
               -int/--instanceTag Tag family with instance association in sampleInfo file (defaults to "")
               -rea/--researchEthicalApproval Tag for displaying research candidates in Scout (defaults to "notApproved")
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
               -pAR/--pAnalysisRunStatus Sets the analysis run status flag to finished in sampleInfoFile (defaults to "1" (=yes))

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
                 -chbdb/--chanjoBuildDb  Reference database (defaults to "CCDS.current.txt")
               -pCh_C/--pChanjoCalculate Chanjo coverage analysis (defaults to "1" (=yes))
                 -chccut/--chanjoCalculateCutoff Read depth cutoff (defaults to "10")
               -pCh_I/--pChanjoImport Chanjo import to collect sample info to family Db  (defaults to "0" (=no))
               -pCC_bedgc/--pGenomeCoverageBED Genome coverage calculation using genomeCoverageBED (defaults to "1" (=yes))
               -xcov/--xCoverage Max coverage depth when using '-pGenomeCoverageBED' (defaults to "30")
               -pCC_picmm/--pPicardToolsCollectMultipleMetrics Metrics calculation using PicardTools collectMultipleMetrics (defaults to "1" (=yes))
               -pCCE_pichs/--pPicardToolsCalculateHSMetrics Capture calculation using PicardTools CalculateHSmetrics (defaults to "1" (=yes))
                 -extbl/--exomeTargetBedInfileLists Prepared target BED file for PicardTools CalculateHSMetrics (defaults to "". File ending should be ".infile_list") 
                 -extpbl/--exomeTargetPaddedBedInfileLists Prepared padded target BED file for PicardTools CalculateHSMetrics (defaults to "". File ending should be ".padXXX.infile_list")
               -pRCP/--pRCovPlots Plots of genome coverage using rCovPlots (defaults to "1" (=yes))
               
               ##GATK              
               -gatkpath/--genomeAnalysisToolKitPath  Path to GATK. Mandatory for use of GATK (defaults to "")
               -gatkbdv/--GATKBundleDownLoadVersion  GATK FTP bundle download version.(defaults to "2.8")
               -gatktmpd/--GATKTempDirectory Temporary Directory to write to using GATK ReAlignerTargetCreator & BaseRecalibrator (defaults to "/scratch/SLURM_JOB_ID";Supply whole path)
               -gatktpbl/--GATKTargetPaddedBedIntervalLists Target BED file interval for GATK (defaults to "". File ending should be ".padXXX.interval_list")
               -gatkdcov/--GATKDownSampleToCoverage Coverage to downsample to at any given locus (defaults to "1000")
               -pGATK_real/--pGATKRealigner Realignments of reads using GATK realign (defaults to "1" (=yes))
                 -gatkrealknset1/--GATKReAlignerINDELKnownSet1 GATK ReAlignerTargetCreator/IndelRealigner known INDEL set 1 (defaults to "1000G_phase1.indels.b37.vcf")
                 -gatkrealknset2/--GATKReAlignerINDELKnownSet2 GATK ReAlignerTargetCreator/IndelRealigner known INDEL set 2 (defaults to "Mills_and_1000G_gold_standard.indels.b37.vcf")
               -pGATK_baserecal/--pGATKBaseRecalibration Recalibration of bases using GATK BaseRecalibrator/PrintReads (defaults to "1" (=yes))
                 -gatkbaserecalknset/--GATKBaseReCalibrationSNPKnownSet GATK BaseReCalinbration known SNP set (defaults to "dbsnp_138.b37.vcf")                
               -pGATK_hapcall/--pGATKHaploTypeCaller Variant discovery using GATK HaplotypeCaller (defaults to "1" (=yes))
                 -gatkhapcallsnpknset/--GATKHaploTypeCallerSNPKnownSet GATK HaplotypeCaller dbSNP set for annotating ID columns (defaults to "dbsnp_138.b37.vcf")
               -pGATK_genotype/--pGATKGenoTypeGVCFs Merge gVCF records using GATK GenotypeGVCFs (defaults to "1" (=yes))
                 -gatkgenotyperefgvcfinfile/--GATKGenoTypeGVCFsRefGVCFInfile GATK GenoTypeGVCFs gVCF reference infile list for joint genotyping (defaults to "")
               -pGATK_varrecal/--pGATKVariantRecalibration Variant recalibration using GATK VariantRecalibrator/ApplyRecalibration (defaults to "1" (=yes))
                 -gatkexrefsnp/--GATKExomeReferenceSNPs Prepared exome reference file (SNVs) for GATKVariantRecalibration (defaults to "")
                 -gatkvarrecaltrhapmap/--GATKVariantReCalibrationTrainingSetHapMap GATK VariantRecalibrator HapMap training set (defaults to "hapmap_3.3.b37.sites.vcf")
                 -gatkvarrecaltrdbsnp/--GATKVariantReCalibrationTrainingSetDbSNP GATK VariantRecalibrator dbSNP training set (defaults to "dbsnp_138.b37.vcf")
                 -gatkvarrecaltrd1000Gsnp/--GATKVariantReCalibrationTrainingSet1000GSNP GATK VariantRecalibrator 1000G high confidence SNP training set (defaults to "1000G_phase1.snps.high_confidence.b37.vcf")
                 -gatkvarrecaltromni/--GATKVariantReCalibrationTrainingSet1000GOmni GATK VariantRecalibrator 1000G_omni training set (defaults to "1000G_omni2.5.b37.sites.vcf")
                 -gatkvarrecaltrdbmills/--GATKVariantReCalibrationTrainingSetMills GATK VariantRecalibrator Mills training set (defaults to "Mills_and_1000G_gold_standard.indels.b37.vcf")
                 -gatkvarrecaltsfilterlevel/--GATKVariantReCalibrationTSFilterLevel The truth sensitivity level at which to start filtering used in GATK VariantRecalibrator (defaults to "99.9")
               -pGATK_phaseTr/--pGATKPhaseByTransmission Computes the most likely genotype and phases calls were unamibigous using GATK PhaseByTransmission (defaults to "1" (=yes))
               -pGATK_readPh/--pGATKReadBackedPhasing Performs physical phasing of SNP calls, based on sequencing reads using GATK ReadBackedPhasing (defaults to "1" (=yes))
                 -gatkreadphphaseqthr/--GATKReadBackedPhasingPhaseQualityThresh The minimum phasing quality score required to output phasing
               -pGATK_varevalall/--pGATKVariantEvalAll Variant evaluation using GATK VariantEval for all variants  (defaults to "1" (=yes))
               -pGATK_varevalexome/--pGATKVariantEvalExome Variant evaluation using GATK VariantEval for exonic variants  (defaults to "1" (=yes))
                 -gatkvarevaldbsnp/--GATKVariantEvalDbSNP DbSNP file used in GATK VariantEval (defaults to "dbsnp_138.b37.excluding_sites_after_129.vcf")
                 -gatkvarevaldbgold/--GATKVariantEvalGold Gold Indel file used in GATK VariantEval (defaults to "Mills_and_1000G_gold_standard.indels.b37.vcf")
               
               ##ANNOTATION
               -pVEP/--pVariantEffectPredictor Annotate variants using VEP (defaults to "1" (=yes))
                 -vepdirpath/--vepDirectoryPath Path to VEP script directory (Supply whole path, defaults to "")
                 -vepdircache/vepDirectoryCache Specify the cache directory to use (Supply whole path, defaults to "") 
                 -vepfeat/--vepFeatures VEP features (comma sep)
               -pVCFPar/--pVCFParser Parse variants using vcfParser.pl (defaults to "1" (=yes))
                 -vcfparservepp/--vcfParserVEPParser Parse VEP transcript specific entries 
                 -vcfparserrf/--vcfParserRangeFeatureFile Range annotations file (tab-sep)
                 -vcfparserrfac/--vcfParserRangeFeatureAnnotationColumns Range annotations feature columns (Defaults to ""; comma sep)
                 -vcfparsersf/--vcfParserSelectFile File containging list of genes to analyse seperately (Defaults to "";tab-sep file and HGNC Symbol required)
                 -vcfparsersfmc/--vcfParserSelectFileMatchingColumn Position of HGNC Symbol column in SelectFile (Defaults to "")
                 -vcfparsersfac/--vcfParserSelectFeatureAnnotationColumns Feature columns to use in annotation (Defaults to ""; comma sep)
                -pSnpeff/--pSnpEff Variant annotation using snpEFF (defaults to "1" (=yes))
                 -snpeffpath/--snpEffPath Path to snpEff. Mandatory for use of snpEff (defaults to "")
                 -snpsiftanf/--snpSiftAnnotationFiles Annotation files to use with snpSift (comma sep)
                 -snpsiftdbnsfpf/--snpSiftDbNSFPFile DbNSFP File (defaults to "dbNSFP2.4_variant.txt.gz")
                 -snpsiftdbnsfpan/--snpSiftDbNSFPAnnotations DbNSFP annotations to use with snpSift (comma sep)
               -pANVAR/--pAnnovar Annotate variants using Annovar (defaults to "1" (=yes))
                 -anvarpath/--annovarPath  Path to Annovar script directory (Supply whole path, defaults to "". NOTE: Assumes that the annovar db files are located in annovar/humandb)
                 -anvargbv/--annovarGenomeBuildVersion Annovar genome build version (defaults to "hg19")
                 -anvartn/--annovarTableNames Annovar table names (comma sep)
                 -anvarstn/--annovarSupportedTableNames Print Annovar MIP supported table names (defaults 0 (=no))
                 -anvarmafth/--annovarMAFThreshold Sets the minor allele frequency threshold in annovar (defaults to "0")
               
               ##VMerge  
               -pMerge_anvar/--pMergeAnnotatedVariants Merge (& annotate) all annotated variants into one file using intersectCollect.pl (defaults to "1" (=yes))
                 -mergeanvarte/--mergeAnnotatedVariantsTemplateFile Db template file used to create the specific family '-mergeanvardbf' master file (defaults to "")
                 -mergeanvardbf/--mergeAnnotatedVariantsDbFile Db master file to be used in intersectCollect.pl (defaults to  "{outDataDir}/{familyID}/{familyID}_intersectCollect_db_master.txt";Supply whole path)
                 -mergeanvartese/--mergeAnnotatedVariantsTemplateSelectFile Select Db template file used to create the specific family '-mergeanvardbf' master file (defaults to "")
                 -mergeanvardbsef/--mergeAnnotatedVariantsDbSelectFile Select Db master file to be used in intersectCollect.pl (defaults to  "{outDataDir}/{familyID}/{familyID}_intersectCollect_db_master.select.txt";Supply whole path)

               ##Add_depth
               -pAdd_dp/--pAddDepth Adds read depth at nonvariant sites using SamTools mpileup and add_depth.pl (defaults to "1" (=yes))

               ##RankVariants
               -pRankVar/--pRankVariants Ranking of annotated variants (defaults to "1" (=yes))
                 -rs/--rankScore The rank score cut-off (defaults to "-100", .i.e. include everything
                 -gf/--geneFile Defines genes to use when calculating compounds
                 -cadd/--caddWGSSNVs Annotate whole genome sequencing CADD score (defaults to "0" (=no))
                 -caddf/--caddWGSSNVsFile Whole genome sequencing CADD score file (defaults to "")
                 -cadd1kg/--cadd1000Genomes 1000 Genome cadd score file (defaults to "0" (=no))
                 -cadd1kgf/--cadd1000GenomesFile 1000 Genome cadd score file (defaults to "")
                 -gene/--wholeGene Allow compound pairs in intronic regions (defaults to "1" (=yes))
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

&DefineParameters("emailType", "MIP", "F", "MIP");

&DefineParametersPath("familyID", "nodefault", "MIP", 0);

&DefineParameters("maximumCores", "MIP", 16, "MIP");

&DefineParametersPath("configFile", 0, "MIP", "file");

&DefineParameters("analysisType", "MIP", "exomes", "MIP");

&DefineParametersPath("outDataDir", "nodefault", "MIP", 0);

&DefineParametersPath("outScriptDir", "nodefault", "MIP", 0);

&DefineParametersPath("writeConfigFile", 0, "MIP", 0);

&DefineParametersPath("pedigreeFile", "nodefault", "MIP", "file", "noAutoBuild");

&DefineParametersPath("sampleInfoFile", "NotsetYet", "MIP", "file", "noAutoBuild");

&DefineParameters("instanceTag", "MIP", "Unknown", "MIP");

&DefineParameters("researchEthicalApproval", "MIP", "notApproved", "MIP");

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

##AnalysisRunStatus
&DefineParameters("pAnalysisRunStatus", "program", 1, "MIP", "", "MAIN");

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

&DefineParametersPath("chanjoBuildDb", "CCDS.current.txt", "pChanjoBuild", "file", "yesAutoDownLoad");

&DefineParameters("pChanjoCalculate", "program", 1, "MIP","_coverage", "CoverageReport");

&DefineParameters("chanjoCalculateCutoff", "program", 10, "pChanjoCalculate");

&DefineParameters("pChanjoImport", "program", 0, "MIP", "nofileEnding", "CoverageReport");

&DefineParameters("pGenomeCoverageBED", "program", 1, "MIP", "_genomeCoverageBed", "CoverageQC_GcovBed", "bedtools");

&DefineParameters("pPicardToolsCollectMultipleMetrics", "program", 1, "MIP", "nofileEnding", "CoverageQC_PTCMM");

&DefineParameters("pPicardToolsCalculateHSMetrics", "program", 1, "MIP", "nofileEnding", "CoverageQC_PTCHSM");

&DefineParameters("xCoverage", "program", 30, "pGenomeCoverageBED");

&DefineParameters("pRCovPlots", "program", 0, "MIP", "nofileEnding", "CoverageQC_RCOVP");

&DefineParametersPath("picardToolsPath", "nodefault", "pBwaMem,pPicardToolsMergeSamFiles,pPicardToolsMarkduplicates,pPicardToolsCalculateHSMetrics,pPicardToolsCollectMultipleMetrics,pGATKHaploTypeCaller,pGATKVariantRecalibration", "directory"); #pGATKHaploTypeCaller,pGATKVariantRecalibration since these jars can use merged interval_list files, which are created in MIP with picardTools

##Target definition files
my (@exomeTargetBedInfileLists, @exomeTargetPaddedBedInfileLists); #Arrays for target bed infile lists

##GATK

&DefineParameters("pGATKRealigner", "program", 1, "MIP", "_rreal", "MAIN");

&DefineParametersPath("GATKReAlignerINDELKnownSet1", "1000G_phase1.indels.b37.vcf", "pGATKRealigner", "file", "yesAutoDownLoad");

&DefineParametersPath("GATKReAlignerINDELKnownSet2", "Mills_and_1000G_gold_standard.indels.b37.vcf", "pGATKRealigner", "file", "yesAutoDownLoad");


&DefineParameters("pGATKBaseRecalibration", "program", 1, "MIP", "_brecal", "MAIN");

&DefineParametersPath("GATKBaseReCalibrationSNPKnownSet", "dbsnp_138.b37.vcf", "pGATKBaseRecalibration", "file", "yesAutoDownLoad");


&DefineParameters("pGATKHaploTypeCaller", "program", 1, "MIP", "_gvcf", "MAIN");

&DefineParametersPath("GATKHaploTypeCallerSNPKnownSet", "dbsnp_138.b37.vcf", "pGATKHaploTypeCaller", "file", "yesAutoDownLoad");


&DefineParameters("pGATKGenoTypeGVCFs", "program", 1, "MIP", "_", "MAIN");

&DefineParametersPath("GATKGenoTypeGVCFsRefGVCFInfile", "nodefault", "pGATKGenoTypeGVCFs", "file", "noAutoBuild");


&DefineParameters("pGATKVariantRecalibration", "program", 1, "MIP", "vrecal_", "MAIN");

&DefineParametersPath("GATKExomeReferenceSNPs", "nodefault", "pGATKVariantRecalibration", "file", "noAutoBuild");

&DefineParametersPath("GATKVariantReCalibrationTrainingSetHapMap", "hapmap_3.3.b37.sites.vcf", "pGATKVariantRecalibration", "file", "yesAutoDownLoad");

&DefineParametersPath("GATKVariantReCalibrationTrainingSetDbSNP", "dbsnp_138.b37.vcf", "pGATKVariantRecalibration", "file", "yesAutoDownLoad");

&DefineParametersPath("GATKVariantReCalibrationTrainingSet1000GSNP", "1000G_phase1.snps.high_confidence.b37.vcf", "pGATKVariantRecalibration", "file", "yesAutoDownLoad");

&DefineParametersPath("GATKVariantReCalibrationTrainingSet1000GOmni", "1000G_omni2.5.b37.sites.vcf", "pGATKVariantRecalibration", "file", "yesAutoDownLoad");

&DefineParametersPath("GATKVariantReCalibrationTrainingSetMills", "Mills_and_1000G_gold_standard.indels.b37.vcf", "pGATKVariantRecalibration", "file", "yesAutoDownLoad");

&DefineParameters("GATKVariantReCalibrationTSFilterLevel", "program", 99.9, "pGATKVariantRecalibration");

 
&DefineParameters("pGATKPhaseByTransmission", "program", 1, "MIP", "phtr_", "Phasing");

&DefineParameters("pGATKReadBackedPhasing", "program", 1, "MIP", "phrb_", "Phasing");

&DefineParameters("GATKReadBackedPhasingPhaseQualityThresh", "program", 20, "pGATKReadBackedPhasing");


&DefineParameters("pGATKVariantEvalAll", "program", 1, "MIP", "nofileEnding", "AllVariantQC");

&DefineParameters("pGATKVariantEvalExome", "program", 1, "MIP", "nofileEnding", "ExomeVarintQC", "bedtools");

&DefineParametersPath("GATKVariantEvalDbSNP", "dbsnp_138.b37.excluding_sites_after_129.vcf", "pGATKVariantEvalAll,pGATKVariantEvalExome", "file", "yesAutoDownLoad");

&DefineParametersPath("GATKVariantEvalGold", "Mills_and_1000G_gold_standard.indels.b37.vcf", "pGATKVariantEvalAll,pGATKVariantEvalExome", "file", "yesAutoDownLoad");

&DefineParametersPath("genomeAnalysisToolKitPath", "nodefault", "pGATKRealigner,pGATKBaseRecalibration,pGATKHaploTypeCaller,pGATKVariantRecalibration,pGATKPhaseByTransmission,pGATKReadBackedPhasing,pGATKVariantEvalAll,pGATKVariantEvalExome", "directory");

&DefineParametersPath("GATKTempDirectory", "/scratch/", "pGATKRealigner,pGATKBaseRecalibration,pGATKHaploTypeCaller,pGATKVariantRecalibration,pGATKReadBackedPhasing", 0); #Depends on -projectID input, directory created by sbatch script and '$SLURM_JOB_ID' is appended to TMP directory

&DefineParameters("GATKDownSampleToCoverage", "program", 1000, "pGATKRealigner,pGATKBaseRecalibration,pGATKHaploTypeCaller");

&DefineParameters("GATKBundleDownLoadVersion", "program", "2.8", "pBwaMem,pBwaAln,pBwaSampe,pGATKRealigner,pGATKBaseRecalibration,pGATKHaploTypeCaller,pGATKHaploTypeCallerCombineVariants,pGATKVariantRecalibration,pGATKPhaseByTransmission,pGATKReadBackedPhasing,pGATKVariantEvalAll,pGATKVariantEvalExome,pAnnovar,pAddDepth,pPicardToolsCalculateHSMetrics,pPicardToolsCollectMultipleMetrics"); #Sets the GATK FTP Bundle Download version. Needed for all programs that download the human genome reference

my (@GATKTargetPaddedBedIntervalLists); #Array for target infile lists used in GATK

##VEP

&DefineParameters("pVariantEffectPredictor", "program", 1, "MIP", "vep_", "MAIN");

&DefineParametersPath("vepDirectoryPath", "nodefault", "pVariantEffectPredictor", "directory"); #Note not projectID specific

&DefineParametersPath("vepDirectoryCache", "nodefault", "pVariantEffectPredictor", "directory");

my @vepFeatures; #List of VEP features to be used

##VCFParser

&DefineParameters("pVCFParser", "program", 1, "MIP", "parsed_", "MAIN");

&DefineParameters("vcfParserVEPParser", "program", 0, "pVCFParser");

&DefineParametersPath("vcfParserRangeFeatureFile", "noUserInfo", "pVCFParser", "file"); 

&DefineParametersPath("vcfParserSelectFile", "noUserInfo", "pVCFParser", "file"); 

&DefineParameters("vcfParserSelectFileMatchingColumn", "program", "nodefault", "pVCFParser");

my $VEPOutputFiles = 1; #To track if VEPParser was used with a vcfParserSelectFile (=2) or not (=1)
my @vcfParserRangeFeatureAnnotationColumns;
my @vcfParserSelectFeatureAnnotationColumns;


##Annovar

&DefineParameters("pAnnovar", "program", 1, "MIP", "annovar_", "Annovar");

&DefineParametersPath("annovarPath", "nodefault", "pAnnovar", "directory"); #Note not projectID specific

&DefineParameters("annovarGenomeBuildVersion", "program", "hg19", "pAnnovar");

&DefineParameters("annovarSupportedTableNames", "program", 0, "pAnnovar");

&DefineParameters("annovarMAFThreshold", "program", 0, "pAnnovar");

my @annovarTableNames; #List of Annovar table names to be used

##SnpEFF

&DefineParameters("pSnpEff", "program", 1, "MIP", "snpeff_", "MAIN");

&DefineParametersPath("snpEffPath", "nodefault", "pSnpEff", "directory");

&DefineParametersPath("snpSiftDbNSFPFile", "dbNSFP2.4_variant.txt.gz", "pSnpEff", "file");

my @snpSiftDbNSFPAnnotations;
my @snpSiftAnnotationFiles;

##Special case GATKPAth since in VEP and SnpEff module GATK CombineVariants is used to merge vcfs 
&DefineParametersPath("javaUseLargePages", "no", "pGATKRealigner,pGATKBaseRecalibration,pGATKHaploTypeCaller,pGATKHaploTypeCallerCombineVariants,pGATKVariantRecalibration,pGATKPhaseByTransmission,pGATKReadBackedPhasing,pGATKVariantEvalAll,pGATKVariantEvalExome,pVariantEffectPredictor,pSnpEff");

&DefineParametersPath("genomeAnalysisToolKitPath", "nodefault", "pGATKRealigner,pGATKBaseRecalibration,pGATKHaploTypeCaller,pGATKVariantRecalibration,pGATKPhaseByTransmission,pGATKReadBackedPhasing,pGATKVariantEvalAll,pGATKVariantEvalExome,pVariantEffectPredictor,pSnpEff", "directory");

##SChecks
&DefineParameters("pSampleCheck", "program", 1, "MIP", "nofileEnding", "IDQC", "vcftools:plink");


##VMerge

&DefineParameters("pMergeAnnotatedVariants", "program", 1, "MIP", "merged_", "Annovar");

&DefineParametersPath("mergeAnnotatedVariantsTemplateFile", "nodefault", "pMergeAnnotatedVariants", "file", "noAutoBuild");

&DefineParameters("mergeAnnotatedVariantsDbFile", "program", "notSetYet", "pMergeAnnotatedVariants"); #No file check since file is created by MIP later

&DefineParametersPath("mergeAnnotatedVariantsTemplateSelectFile", "nodefault", "pMergeAnnotatedVariants", "file", "noAutoBuild");

&DefineParameters("mergeAnnotatedVariantsDbSelectFile", "program", "notSetYet", "pMergeAnnotatedVariants"); #No file check since file is created by MIP later

##Add_depth

&DefineParameters("pAddDepth", "program", 0, "MIP", "", "Annovar");


##RankVariants

&DefineParameters("pRankVariants", "program", 1, "MIP", "ranked_", "MAIN");

&DefineParameters("rankScore", "program", -100, "pRankVariants");

&DefineParametersPath("geneFile", "hg19_refGene.txt", "pRankVariants", "file", "noAutoBuild");

&DefineParameters("caddWGSSNVs", "program", 0, "pRankVariants");

&DefineParametersPath("caddWGSSNVsFile", "nodefault", "pRankVariants", "file", "noAutoBuild");

&DefineParameters("cadd1000Genomes", "program", 0, "pRankVariants");

&DefineParametersPath("cadd1000GenomesFile", "nodefault", "pRankVariants", "file", "noAutoBuild");

&DefineParameters("wholeGene", "program", 1, "pRankVariants");

&DefineParametersPath("pythonVirtualEnvironment", "nodefault", "pChanjoBuild,pChanjoCalculate,pChanjoImport,pRankVariants");

##QcCollect

&DefineParameters("pQCCollect", "program", 1, "MIP", "nofileEnding", "MAIN");

&DefineParameters("QCCollectSampleInfoFile", "program", "notSetYet", "pQCCollect"); #No file check since file is created by MIP later

&DefineParametersPath("QCCollectRegExpFile", "qc_regexp.yaml", "pQCCollect", "file", "noAutoBuild");

##MIP

##humanGenomeReference
&DefineParametersPath("humanGenomeReference", "Homo_sapiens.GRCh37.d5.fasta", "pBwaMem,pBwaAln,pBwaSampe,pGATKRealigner,pGATKBaseRecalibration,pGATKHaploTypeCaller,pGATKHaploTypeCallerCombineVariants,pGATKVariantRecalibration,pGATKPhaseByTransmission,pGATKReadBackedPhasing,pGATKVariantEvalAll,pGATKVariantEvalExome,pAnnovar,pAddDepth,pPicardToolsCalculateHSMetrics,pPicardToolsCollectMultipleMetrics", "file", "yesAutoDownLoad");
my @humanGenomeReferenceFileEndings = (".dict", ".fasta.fai"); #Meta files

my ($humanGenomeReferenceSource, $humanGenomeReferenceVersion, $humanGenomeReferenceNameNoEnding, $humanGenomeCompressed, $fnend, $aligner, $fileName, $fileNameTracker, $version, $help) = ("nocmdinput", "nocmdinput", "nocmdinput", "nocmdinput", ".sh", "nocmdinput", "nocmdinput", 0);

my (@contigs);

my (%infile, %indirpath, %infilesLaneNoEnding, %lane, %infilesBothStrandsNoEnding, %jobID, %sampleInfo); 

####Staging/Sanity Check Area 

##Capture kits supported from pedigree file.
my %supportedCaptureKits = (
    'Nimblegen_SeqCapEZExome.V2' => "Nimblegen_SeqCapEZExome.V2.GenomeReferenceSourceVersion_targets.bed",
    'Nimblegen_SeqCapEZExome.V3' => "Nimblegen_SeqCapEZExome.V3.GenomeReferenceSourceVersion_targets.bed",
    'Agilent_SureSelect.V2' => "Agilent_SureSelect.V2.GenomeReferenceSourceVersion_targets.bed",
    'Agilent_SureSelect.V3' => "Agilent_SureSelect.V3.GenomeReferenceSourceVersion_targets.bed",
    'Agilent_SureSelect.V4' => "Agilent_SureSelect.V4.GenomeReferenceSourceVersion_targets.bed",
    'Agilent_SureSelect.V5' => "Agilent_SureSelect.V5.GenomeReferenceSourceVersion_targets.bed",
    'Latest' => "Agilent_SureSelect.V5.GenomeReferenceSourceVersion_targets.bed",
    );

my %plinkPedigree; #Holds allowed entries and positions for Plink pedigree files
my %supportedCosmidReferences; #References supported as downloads from Cosmid. Hash is populated after user options are processed

my %referenceFileEndings = (
    'mosaikAlignReference' => ".dat",
    'mosaikJumpDbStub' => "_jdb_15",
    'bwaBuildReference' => "",
    'exomeTargetBedInfileLists' => ".infile_list",
    'exomeTargetPaddedBedInfileLists' => ".pad100.infile_list",
    'GATKTargetPaddedBedIntervalLists' => ".pad100.interval_list",
    );

##Set supported annovar table name filtering options
my @annovarSupportedTableNames = ("refGene", "knownGene", "ensGene", "mce46way", "gerp++elem", "segdup", "gwascatalog", "tfbs", "mirna", "snp137", "snp135", "snp132", "snp131", "snp130", "snp129", "snp137NonFlagged", "snp135NonFlagged", "snp132NonFlagged", "snp131NonFlagged", "snp130NonFlagged", "1000g2012apr_all", "1000g2012apr_amr", "1000g2012apr_eur", "1000g2012apr_asn", "1000g2012apr_afr", "1000g2012feb_all", "esp6500si_all", "esp6500_all", "esp6500_aa", "esp6500_ea", "esp5400_all", "esp5400_aa", "esp5400_ea","clinvar_20131105", "ljb2_sift", "ljb2_pp2hdiv", "ljb2_pp2hvar", "ljb2_mt", "ljb2_ma", "ljb2_fathmm", "ljb2_siphy", "ljb2_lrt", "ljb_all", "ljb2_gerp++", "ljb2_phylop", "caddgt20", "caddgt10"); #Used to print list of supported table names

my %annovarTables;

###User Options
#$parameter{''}{'value'}
GetOptions('ifd|inFilesDirs:s'  => \@inFilesDirs, #Comma separated list
	   'isd|inScriptDir:s'  => \$parameter{'inScriptDir'}{'value'}, #Directory for custom scripts required by the pipeline
	   'rd|referencesDir:s'  => \$parameter{'referencesDir'}{'value'}, #directory containing references
	   'p|projectID:s'  => \$parameter{'projectID'}{'value'},
	   's|sampleIDs:s'  => \@sampleIDs, #Comma separated list, one below outDataDir
	   'em|email:s'  => \$parameter{'email'}{'value'}, #Email adress
	   'emt|emailType:s'  => \$parameter{'emailType'}{'value'}, #Email type 
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
	   'int|instanceTag:s' => \$parameter{'instanceTag'}{'value'},
	   'rea|researchEthicalApproval:s' => \$parameter{'researchEthicalApproval'}{'value'},
	   'dra|dryRunAll:n' => \$parameter{'dryRunAll'}{'value'},
	   'pve|pythonVirtualEnvironment:s' => \$parameter{'pythonVirtualEnvironment'}{'value'},
	   'julp|javaUseLargePages:s' => \$parameter{'javaUseLargePages'}{'value'},
	   'h|help' => \$help, #Display help text
	   'v|version' => \$version, #Display version number
	   'pGZ|pGZip:n' => \$parameter{'pGZip'}{'value'},
	   'pFQC|pFastQC:n' => \$parameter{'pFastQC'}{'value'},
	   'pREM|pRemovalRedundantFiles:n' => \$parameter{'pRemovalRedundantFiles'}{'value'},
	   'pAR|pAnalysisRunStatus:n' => \$parameter{'pAnalysisRunStatus'}{'value'}, #AnalysisRunStatus change flag in sampleInfo file if allowed to execute
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
	   'xcov|xCoverage:n' => \$parameter{'xCoverage'}{'value'}, #Sets max depth to calculate coverage
	   'pCC_picmm|pPicardToolsCollectMultipleMetrics:n' => \$parameter{'pPicardToolsCollectMultipleMetrics'}{'value'},
	   'pCCE_pichs|pPicardToolsCalculateHSMetrics:n' => \$parameter{'pPicardToolsCalculateHSMetrics'}{'value'},
	   'extbl|exomeTargetBedInfileLists:s' => \@exomeTargetBedInfileLists, #Comma separated list of target file for CalculateHsMetrics
	   'extpbl|exomeTargetPaddedBedInfileLists:s' => \@exomeTargetPaddedBedInfileLists, #Comma separated list of padded target file for CalculateHsMetrics
	   'pRCP|pRCovPlots:n' => \$parameter{'pRCovPlots'}{'value'},
	   'gatkpath|genomeAnalysisToolKitPath:s' => \$parameter{'genomeAnalysisToolKitPath'}{'value'}, #GATK whole path
	   'gatkbdv|GATKBundleDownLoadVersion:s' => \$parameter{'GATKBundleDownLoadVersion'}{'value'}, #Sets the GATK FTP Bundle Download version
	   'gatktmpd|GATKTempDirectory:s' => \$parameter{'GATKTempDirectory'}{'value'}, #GATK ReAlignerTargetCreator & BaseRecalibrator temporary directory
	   'gatktpbl|GATKTargetPaddedBedIntervalLists:s' => \@GATKTargetPaddedBedIntervalLists, #Comma separated list of padded target file set to be used in GATK
	   'gatkdcov|GATKDownSampleToCoverage:n' => \$parameter{'GATKDownSampleToCoverage'}{'value'}, #GATK downsample to coverage
	   'pGATK_real|pGATKRealigner:n' => \$parameter{'pGATKRealigner'}{'value'}, #GATK ReAlignerTargetCreator/IndelRealigner
	   'gatkrealknset1|GATKReAlignerINDELKnownSet1:s' => \$parameter{'GATKReAlignerINDELKnownSet1'}{'value'}, #Known INDEL set to be used in GATK ReAlignerTargetCreator/IndelRealigner
	   'gatkrealknset2|GATKReAlignerINDELKnownSet2:s' => \$parameter{'GATKReAlignerINDELKnownSet2'}{'value'}, #Known INDEL set to be used in GATK ReAlignerTargetCreator/IndelRealigner
	   'pGATK_baserecal|pGATKBaseRecalibration:n' => \$parameter{'pGATKBaseRecalibration'}{'value'}, #GATK BaseRecalibrator/PrintReads
	   'gatkbaserecalknset|GATKBaseReCalibrationSNPKnownSet:s' => \$parameter{'GATKBaseReCalibrationSNPKnownSet'}{'value'}, #Known SNP set to be used in GATK BaseRecalibrator/PrintReads
	   'pGATK_hapcall|pGATKHaploTypeCaller:n' => \$parameter{'pGATKHaploTypeCaller'}{'value'}, #GATK Haplotypecaller
	   'gatkhapcallsnpknset|GATKHaploTypeCallerSNPKnownSet:s' => \$parameter{'GATKHaploTypeCallerSNPKnownSet'}{'value'}, #Known SNP set to be used in GATK HaplotypeCaller
	   'pGATK_genotype|pGATKGenoTypeGVCFs:n' => \$parameter{'pGATKGenoTypeGVCFs'}{'value'}, #Merge gVCF records using GATK GenotypeGVCFs
	   'gatkgenotyperefgvcfinfile|GATKGenoTypeGVCFsRefGVCFInfile:s' => \$parameter{'GATKGenoTypeGVCFsRefGVCFInfile'}{'value'}, #GATK GenoTypeGVCFs gVCF reference infile list for joint genotyping
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
	   'gatkreadphphaseqthr|GATKReadBackedPhasingPhaseQualityThresh:n' => \$parameter{'GATKReadBackedPhasingPhaseQualityThresh'}{'value'}, #quality score required to output phasing
	   'pGATK_varevalall|pGATKVariantEvalAll:n' => \$parameter{'pGATKVariantEvalAll'}{'value'}, #GATK varianteval all variants
	   'pGATK_varevalexome|pGATKVariantEvalExome:n' => \$parameter{'pGATKVariantEvalExome'}{'value'}, #GATK varianteval only exonic variants
	   'gatkvarevaldbsnp|GATKVariantEvalDbSNP:s' => \$parameter{'GATKVariantEvalDbSNP'}{'value'},
	   'gatkvarevaldbgold|GATKVariantEvalGold:s' => \$parameter{'GATKVariantReCalibrationTrainingSetMills'}{'value'},
	   'pVEP|pVariantEffectPredictor:n' => \$parameter{'pVariantEffectPredictor'}{'value'}, #Annotation of variants using vep
	   'vepdirpath|vepDirectoryPath:s'  => \$parameter{'vepDirectoryPath'}{'value'}, #path to vep script dir
	   'vepdircache|vepDirectoryCache:s'  => \$parameter{'vepDirectoryCache'}{'value'}, #path to vep cache dir
	   'vepfeat|vepFeatures:s'  => \@vepFeatures, #Comma separated list
	   'pVCFPar|pVCFParser:n' => \$parameter{'pVCFParser'}{'value'},
	   'vepparservepp|vcfParserVEPParser:n' => \$parameter{'vcfParserVEPParser'}{'value'},
	   'vepparserrf|vcfParserRangeFeatureFile:s'  => \$parameter{'vcfParserRangeFeatureFile'}{'value'}, #path to vcfParserRangeFeatureFile
	   'vepparserrfac|vcfParserRangeFeatureAnnotationColumns:s'  => \@vcfParserRangeFeatureAnnotationColumns, #Comma separated list
	   'vepparsersf|vcfParserSelectFile:s'  => \$parameter{'vcfParserSelectFile'}{'value'}, #path to vcfParserSelectFile
	   'vepparsersfmc|vcfParserSelectFileMatchingColumn:n' => \$parameter{'vcfParserSelectFileMatchingColumn'}{'value'}, #Column of HGNC Symbol in SelectFile
	   'vepparsersfac|vcfParserSelectFeatureAnnotationColumns:s'  => \@vcfParserSelectFeatureAnnotationColumns, #Comma separated list
	   'snpeffpath|snpEffPath:s'  => \$parameter{'snpEffPath'}{'value'}, #path to snpEff directory
	   'pSnpeff|pSnpEff:n' => \$parameter{'pSnpEff'}{'value'},
	   'snpsiftanf|snpSiftAnnotationFiles:s'  => \@snpSiftAnnotationFiles, #Comma separated list
	   'snpsiftdbnsfpf|snpSiftDbNSFPFile:s'  => \$parameter{'snpSiftDbNSFPFile'}{'value'}, #DbNSFP file
	   'snpsiftdbnsfpan|snpSiftDbNSFPAnnotations:s'  => \@snpSiftDbNSFPAnnotations, #Comma separated list
	   'pANVAR|pAnnovar:n' => \$parameter{'pAnnovar'}{'value'}, #Performs annovar filter gene, region and filter analysis
	   'anvarpath|annovarPath:s'  => \$parameter{'annovarPath'}{'value'}, #path to annovar script dir
	   'anvargbv|annovarGenomeBuildVersion:s'  => \$parameter{'annovarGenomeBuildVersion'}{'value'},
	   'anvartn|annovarTableNames:s'  => \@annovarTableNames, #Comma separated list
	   'anvarstn|annovarSupportedTableNames:n' => \$parameter{'annovarSupportedTableNames'}{'value'}, #Generates a list of supported table names
	   'anvarmafth|annovarMAFThreshold:n' => \$parameter{'annovarMAFThreshold'}{'value'},
	   'pMerge_anvar|pMergeAnnotatedVariants:n' => \$parameter{'pMergeAnnotatedVariants'}{'value'}, #Merges annovar analysis results to one master file
	   'mergeanvarte|mergeAnnotatedVariantsTemplateFile:s' => \$parameter{'mergeAnnotatedVariantsTemplateFile'}{'value'}, #Template file to create the specific family db master file
	   'mergeanvardbf|mergeAnnotatedVariantsDbFile:s' => \$parameter{'mergeAnnotatedVariantsDbFile'}{'value'}, #db master file to use when collecting external data
	   'mergeanvartese|mergeAnnotatedVariantsTemplateSelectFile:s' => \$parameter{'mergeAnnotatedVariantsTemplateSelectFile'}{'value'}, #Select template file to create the specific family db master file
	   'mergeanvardbsef|mergeAnnotatedVariantsDbSelectFile:s' => \$parameter{'mergeAnnotatedVariantsDbSelectFile'}{'value'}, #Select db master file to use when collecting external data
	   'pAdd_dp|pAddDepth:n' => \$parameter{'pAddDepth'}{'value'}, #Adds depth (DP) for nonvariants to master file (annovar_merged.txt)
	   'pRankVar|pRankVariants:n' => \$parameter{'pRankVariants'}{'value'}, #Ranking variants
	   'rs|rankscore:n'  => \$parameter{'rankScore'}{'value'}, #The rank score cut-off
	   'gf|geneFile:s' => \$parameter{'geneFile'}{'value'},
	   'cadd|caddWGSSNVs:n' => \$parameter{'caddWGSSNVs'}{'value'},
	   'caddf|caddWGSSNVsFile:s' => \$parameter{'caddWGSSNVsFile'}{'value'},
	   'cadd1kg|cadd1000Genomes:n' => \$parameter{'cadd1000Genomes'}{'value'},
	   'cadd1kgf|cadd1000GenomesFile:s' => \$parameter{'cadd1000GenomesFile'}{'value'},
	   'gene|wholeGene:n'  => \$parameter{'wholeGene'}{'value'}, #Allow compound pairs in intronic regions
	   'pSCheck|pSampleCheck:n' => \$parameter{'pSampleCheck'}{'value'}, #QC for samples gender and relationship
	   'pQCC|pQCCollect:n' => \$parameter{'pQCCollect'}{'value'}, #QCmetrics collect
	   'QCCsampleinfo|QCCollectSampleInfoFile:s' => \$parameter{'QCCollectSampleInfoFile'}{'value'}, #SampleInfo yaml file produced by MIP
	   'QCCregexp|QCCollectRegExpFile:s' => \$parameter{'QCCollectRegExpFile'}{'value'}, #Regular expression yaml file

    );

if($help) {

    print STDOUT $USAGE, "\n";
    exit;
}

my $MipVersion = "v1.5.7";#Set version for log

if($version) {

    print STDOUT "\nMip.pl ".$MipVersion, "\n\n";
    exit;
}
print STDOUT "MIP Version: ".$MipVersion, "\n";

if ($parameter{'annovarSupportedTableNames'}{'value'} eq 1) {

    &PrintSupportedAnnovarTableNames();
}

if ($parameter{'configFile'}{'value'} ne "nocmdinput") { #Input from cmd

    %scriptParameter = &LoadYAML($parameter{'configFile'}{'value'}); #Load parameters from configfile

    &ReplaceConfigParamWithCMDInfo("analysisType");
    &ReplaceConfigParamWithCMDInfo("aligner");

    foreach my $orderParameterElement (@orderParameters) { #Loop through all parameters and update info   

	&UpdateYAML($orderParameterElement, $scriptParameter{'clusterConstantPath'}, $scriptParameter{'analysisConstantPath'}, $scriptParameter{'analysisType'}, $parameter{'familyID'}{'value'}, $scriptParameter{'aligner'} );
    }
}

foreach my $orderParameterElement (@orderParameters) { #Populate scriptParameters{'parameterName'} => 'Value'
    
##3 type of variables: MIP, path or program/program_parameters each is handled in the &AddToScriptParameter subroutine.
##parameterName, parameterValue, parameterType, parameterDefault, AssociatedProgram, Check directory/file existence)    
    &AddToScriptParameter($orderParameterElement, $parameter{$orderParameterElement}{'value'}, $parameter{$orderParameterElement}{'type'}, $parameter{$orderParameterElement}{'default'}, $parameter{$orderParameterElement}{'associatedProgram'}, $parameter{$orderParameterElement}{'existsCheck'}, $parameter{$orderParameterElement}{'programNamePath'});
   
    if ($orderParameterElement eq "outDataDir") { #Set defaults depending on $scriptParameter{'outDataDir'} value that now has been set

	$parameter{'sampleInfoFile'}{'default'} = $scriptParameter{'outDataDir'}."/".$scriptParameter{'familyID'}."/".$scriptParameter{'familyID'}."_qc_sampleInfo.yaml";
	$parameter{'QCCollectSampleInfoFile'}{'default'} = $parameter{'sampleInfoFile'}{'default'};

	$parameter{'mergeAnnotatedVariantsDbFile'}{'default'} = $scriptParameter{'outDataDir'}."/".$scriptParameter{'familyID'}."/".$scriptParameter{'familyID'}."_intersectCollect_db_master.txt";
	$parameter{'mergeAnnotatedVariantsDbSelectFile'}{'default'} = $scriptParameter{'outDataDir'}."/".$scriptParameter{'familyID'}."/".$scriptParameter{'familyID'}."_intersectCollect_db_master.select.txt";
	
    }
    if ($orderParameterElement eq "pedigreeFile") { #Write QC for only pedigree data used in analysis                                                        
	
	if (defined($scriptParameter{'pedigreeFile'})) {
	    `mkdir -p $scriptParameter{'outDataDir'}/$scriptParameter{'familyID'};`;
	    &WriteYAML($scriptParameter{'outDataDir'}."/".$scriptParameter{'familyID'}."/qc_pedigree.yaml", \%sampleInfo);
 
	}
    }
    if ($orderParameterElement eq "humanGenomeReference") { #Supply humanGenomeReference to mosaikAlignReference if required
	
	if ( (defined($scriptParameter{'humanGenomeReference'})) && (defined($humanGenomeReferenceNameNoEnding)) ) {

	    &SetAutoBuildFeature("mosaikAlignReference", \$referenceFileEndings{'mosaikAlignReference'}, \$humanGenomeReferenceNameNoEnding);
	    &SetAutoBuildFeature("mosaikJumpDbStub", \$referenceFileEndings{'mosaikJumpDbStub'}, \$humanGenomeReferenceNameNoEnding);
	    &SetAutoBuildFeature("bwaBuildReference", \$referenceFileEndings{'bwaBuildReference'}, \$humanGenomeReferenceNameNoEnding);	
	}
    }
    if ($orderParameterElement eq "mergeAnnotatedVariantsTemplateFile") { #Check that paths in master template exists

	if (defined($scriptParameter{'mergeAnnotatedVariantsTemplateFile'})) {

	    &CheckTemplateFilesPaths(\($scriptParameter{'referencesDir'}."/".$scriptParameter{'mergeAnnotatedVariantsTemplateFile'}), "mergeAnnotatedVariantsTemplateFile")	    
	}
    }
    if ($orderParameterElement eq "mergeAnnotatedVariantsTemplateSelectFile") { #Check that paths in master template exists
	
	if (defined($scriptParameter{'mergeAnnotatedVariantsTemplateSelectFile'})) {
	    
	    &CheckTemplateFilesPaths(\($scriptParameter{'referencesDir'}."/".$scriptParameter{'mergeAnnotatedVariantsTemplateSelectFile'}), "mergeAnnotatedVariantsTemplateSelectFile")	    
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

if ($scriptParameter{'pVariantEffectPredictor'} > 0) {

    &PrepareArrayParameters(\@vepFeatures, "vepFeatures", "path", "yes", "pVariantEffectPredictor"); #"yes" added to enable addition of default features in &AddToScriptParameters  
}
if ($scriptParameter{'pVCFParser'} > 0) {
    
    &PrepareArrayParameters(\@vcfParserRangeFeatureAnnotationColumns, "vcfParserRangeFeatureAnnotationColumns", "path", "nodefault", "pVCFParser");
    &PrepareArrayParameters(\@vcfParserSelectFeatureAnnotationColumns, "vcfParserSelectFeatureAnnotationColumns", "path", "nodefault", "pVCFParser");
}

if ($scriptParameter{'pAnnovar'} > 0) {

    &PrepareArrayParameters(\@annovarTableNames, "annovarTableNames", "path", "yes", "pAnnovar", "file"); #"yes" added to enable addition of default table names in &AddToScriptParameters  
}

if ($scriptParameter{'pSnpEff'} > 0) {

    &PrepareArrayParameters(\@snpSiftAnnotationFiles, "snpSiftAnnotationFiles", "path", "yes", "pSnpEff"); #"yes" added to enable addition of default features in &AddToScriptParameters 
    &PrepareArrayParameters(\@snpSiftDbNSFPAnnotations, "snpSiftDbNSFPAnnotations", "path", "yes", "pSnpEff"); #"yes" added to enable addition of default features in &AddToScriptParameters  
}

##Set Target files
&PrepareArrayParameters(\@exomeTargetBedInfileLists, "exomeTargetBedInfileLists", "path", "notSetYet", "pPicardToolsCalculateHSMetrics", "file");

&PrepareArrayParameters(\@exomeTargetPaddedBedInfileLists, "exomeTargetPaddedBedInfileLists", "path", "notSetYet", "pPicardToolsCalculateHSMetrics", "file");
 
&PrepareArrayParameters(\@GATKTargetPaddedBedIntervalLists, "GATKTargetPaddedBedIntervalLists", "path", "notSetYet", "pGATKHaploTypeCaller,pGATKVariantRecalibration", "file");

#Cosmid references
&DefineSupportedCosmidReferences("humanGenomeReference", "decoy", "5", \$humanGenomeReferenceVersion, "compressed");
&DefineSupportedCosmidReferences("chanjoBuildDb", "ccds", "latest", \$humanGenomeReferenceVersion, "unCompressed");
&DefineSupportedCosmidReferences("GATKReAlignerINDELKnownSet1", "indel", $scriptParameter{'GATKBundleDownLoadVersion'}."/b".$humanGenomeReferenceVersion, \$humanGenomeReferenceVersion, "unCompressed");
&DefineSupportedCosmidReferences("GATKReAlignerINDELKnownSet2", "mills", $scriptParameter{'GATKBundleDownLoadVersion'}."/b".$humanGenomeReferenceVersion, \$humanGenomeReferenceVersion, "unCompressed");
&DefineSupportedCosmidReferences("GATKBaseReCalibrationSNPKnownSet", "dbsnp", $scriptParameter{'GATKBundleDownLoadVersion'}."/b".$humanGenomeReferenceVersion, \$humanGenomeReferenceVersion, "unCompressed");
&DefineSupportedCosmidReferences("GATKHaploTypeCallerSNPKnownSet", "dbsnp", $scriptParameter{'GATKBundleDownLoadVersion'}."/b".$humanGenomeReferenceVersion, \$humanGenomeReferenceVersion, "unCompressed");
&DefineSupportedCosmidReferences("GATKVariantReCalibrationTrainingSetHapMap", "hapmap", $scriptParameter{'GATKBundleDownLoadVersion'}."/b".$humanGenomeReferenceVersion, \$humanGenomeReferenceVersion, "unCompressed");
&DefineSupportedCosmidReferences("GATKVariantReCalibrationTrainingSetMills", "mills", $scriptParameter{'GATKBundleDownLoadVersion'}."/b".$humanGenomeReferenceVersion, \$humanGenomeReferenceVersion, "unCompressed");
&DefineSupportedCosmidReferences("GATKVariantReCalibrationTrainingSet1000GOmni", "1000g_omni", $scriptParameter{'GATKBundleDownLoadVersion'}."/b".$humanGenomeReferenceVersion, \$humanGenomeReferenceVersion, "unCompressed");
&DefineSupportedCosmidReferences("GATKVariantReCalibrationTrainingSet1000GSNP", "1000g_snps", $scriptParameter{'GATKBundleDownLoadVersion'}."/b".$humanGenomeReferenceVersion, \$humanGenomeReferenceVersion, "unCompressed");
&DefineSupportedCosmidReferences("GATKVariantReCalibrationTrainingSetDbSNP", "dbsnp", $scriptParameter{'GATKBundleDownLoadVersion'}."/b".$humanGenomeReferenceVersion, \$humanGenomeReferenceVersion, "unCompressed");
&DefineSupportedCosmidReferences("GATKVariantEvalGold", "mills", $scriptParameter{'GATKBundleDownLoadVersion'}."/b".$humanGenomeReferenceVersion, \$humanGenomeReferenceVersion, "unCompressed"); 
&DefineSupportedCosmidReferences("GATKVariantEvalGold", "dbsnpex", $scriptParameter{'GATKBundleDownLoadVersion'}."/b".$humanGenomeReferenceVersion, \$humanGenomeReferenceVersion, "unCompressed");

for my $references (keys %supportedCosmidReferences) {

    &CheckCosmidInstallation(\$references);
    last; #Only need to check once per analysis run
}

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
my @printFilehandles = (*STDOUT, *MIPLOG); #Used for printing to several FILEHANDLES

##Add parameters
print MIPLOG "\n".$script." "; #Adds script name to recontruct command line

&WriteCMDMipLog();

print STDOUT "\nScript parameters and info from ".$script." are saved in file: ".$mipLogName, "\n";

##Collect infiles
&CollectInfiles();

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

    &PrintToFileHandles(\@printFilehandles, "\nGZip for fastq files\n");

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
    
    &PrintToFileHandles(\@printFilehandles, "\nFastQC\n");

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {  
	
	&FastQC($sampleIDs[$sampleIDCounter]);	
    }
}

if ($scriptParameter{'pMosaikBuild'} > 0) { #Run MosaikBuild
    
    &PrintToFileHandles(\@printFilehandles, "\nMosaikBuild\n");

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {  
	
	&MosaikBuild($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});	
    }
}


if ($scriptParameter{'pMosaikAlign'} > 0) { #Run MosaikAlign

    &PrintToFileHandles(\@printFilehandles, "\nMosaikAlign\n");

    if ( ($parameter{'humanGenomeReference'}{'buildFile'} eq 1) || ($parameter{'mosaikAlignReference'}{'buildFile'} eq 1) || ($parameter{'mosaikJumpDbStub'}{'buildFile'} eq 1) ) {
		
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
    
    &PrintToFileHandles(\@printFilehandles, "\nBWA Mem\n");
    
    if ( ($parameter{'humanGenomeReference'}{'buildFile'} eq 1) || ($parameter{'bwaBuildReference'}{'buildFile'} eq 1) ) {
	
	&BuildBwaPreRequisites($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BwaMem");
    }

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {  
	
	&BWA_Mem($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});	
	
    }    
}

if ($scriptParameter{'pPicardToolsMergeRapidReads'} > 0) { #Run PicardToolsMergeRapidReads - Relevant only in rapid mode
    
    &PrintToFileHandles(\@printFilehandles, "\nPicardToolsMergeRapidReads\n");
    
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {  
	
        #Merge all read batch processes to 1 file again containing sorted & indexed reads matching clinical test genes
	&PicardToolsMergeRapidReads($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});
    }    
}

if ($scriptParameter{'pBwaAln'} > 0) { #Run BWA Aln
    
    &PrintToFileHandles(\@printFilehandles, "\nBWA Aln\n");

    if ( ($parameter{'humanGenomeReference'}{'buildFile'} eq 1) || ($parameter{'bwaBuildReference'}{'buildFile'} eq 1) ) {
	
	&BuildBwaPreRequisites($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BwaAln");
    }
    
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {  
	
	&BWA_Aln($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});	
    }    
}

if ($scriptParameter{'pBwaSampe'} > 0) { #Run BWA Sampe
    
    &PrintToFileHandles(\@printFilehandles, "\nBWA Sampe\n");

    if ( ($parameter{'humanGenomeReference'}{'buildFile'} eq 1) || ($parameter{'bwaBuildReference'}{'buildFile'} eq 1) ) {
	
	&BuildBwaPreRequisites($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BwaSampe");
    }

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {  
	
	&BWA_Sampe($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});
    }
}

if ($scriptParameter{'pPicardToolsSortSam'} > 0) { #Run Picardtools SortSam and Index

    if ($scriptParameter{'analysisType'} ne "rapid") { #In rapid mode Sort and index is done for each batch of reads in the BWA_Mem call, since the link to infile is broken by the read batch processing. However pPicardToolsSortSam should be enabled to ensure correct fileending and merge the flow to ordinary modules.

    &PrintToFileHandles(\@printFilehandles, "\nPicardTools SortSam & index\n");
	
	for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {  
	    
	    &PicardToolsSortSamIndex($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});
	}
    }
}

if ($scriptParameter{'pPicardToolsMergeSamFiles'} > 0) { #Run picardtools merge

    &PrintToFileHandles(\@printFilehandles, "\nPicardTool MergeSamFiles\n");

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {  

	if ( ($sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{'picardToolsMergeSamFilesPrevious'} == 1) || (scalar( @{ $infilesLaneNoEnding{ $sampleIDs[$sampleIDCounter] } }) > 1) ) { #Sanity Check that we have something to merge with
	
	    &PicardToolsMerge($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'}, $sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{'fileEnding'});	
	}
    }
}

if ($scriptParameter{'pPicardToolsMarkduplicates'} > 0) { #PicardTools MarkDuplicates

    &PrintToFileHandles(\@printFilehandles, "\nPicardTools MarkDuplicates\n");

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {  
    
	&PicardToolsMarkDuplicates($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});	
    }
}

if ($scriptParameter{'pChanjoBuild'} > 0) {
    
    &PrintToFileHandles(\@printFilehandles, "\nChanjoBuild\n");

    &CheckBuildDownLoadPreRequisites("ChanjoBuild");
       
    &ChanjoBuild($scriptParameter{'familyID'});
}

if ($scriptParameter{'pChanjoCalculate'} > 0) {
    
    &PrintToFileHandles(\@printFilehandles, "\nChanjoCalculate\n");
    
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) { #For all SampleIDs
	
	&ChanjoCalculate($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});
    }
}

if ($scriptParameter{'pChanjoImport'} > 0) {
    
    &PrintToFileHandles(\@printFilehandles, "\nChanjoImport\n");
    
    &ChanjoImport($scriptParameter{'familyID'}, $scriptParameter{'aligner'});
}

if ($scriptParameter{'pGenomeCoverageBED'} > 0) { #Run GenomeCoverageBED
    
    &PrintToFileHandles(\@printFilehandles, "\nGenomeCoverageBED\n"); 
    
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {  

	&GenomeCoverageBED($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});
    }
}

if ($scriptParameter{'pPicardToolsCollectMultipleMetrics'} > 0) { #Run PicardToolsCollectMultipleMetrics
    
    &PrintToFileHandles(\@printFilehandles, "\nPicardToolsCollectMultipleMetrics\n");   
    
    &CheckBuildHumanGenomePreRequisites("PicardToolsCollectMultipleMetrics");

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {  

	&PicardToolsCollectMultipleMetrics($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});
    }
}

if ($scriptParameter{'pPicardToolsCalculateHSMetrics'} > 0) { #Run PicardToolsCalculateHSMetrics
    
    &PrintToFileHandles(\@printFilehandles, "\nPicardToolsCalculateHSMetrics\n");   
    
    &CheckBuildHumanGenomePreRequisites("PicardToolsCalculateHSMetrics");
    &CheckBuildPTCHSMetricPreRequisites("PicardToolsCalculateHSMetrics");

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {  

	&PicardToolsCalculateHSMetrics($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});
    }
}

if ($scriptParameter{'pRCovPlots'} > 0) { #Run Rcovplot scripts   

    &PrintToFileHandles(\@printFilehandles, "\nRCovPlots\n");	

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {  
	
	&RCoveragePlots($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});	
    }
}

if ($scriptParameter{'pGATKRealigner'} > 0) { #Run GATK ReAlignerTargetCreator/IndelRealigner

    &PrintToFileHandles(\@printFilehandles, "\nGATK ReAlignerTargetCreator/IndelRealigner\n");

    &CheckBuildHumanGenomePreRequisites("GATKRealigner");
    &CheckBuildDownLoadPreRequisites("GATKRealigner");

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {   
    
	&GATKReAligner($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});	
    }
}

if ($scriptParameter{'pGATKBaseRecalibration'} > 0) { #Run GATK BaseRecalibrator/PrintReads

    &PrintToFileHandles(\@printFilehandles, "\nGATK BaseRecalibrator/PrintReads\n");

    &CheckBuildHumanGenomePreRequisites("GATKBaseRecalibration");
    &CheckBuildDownLoadPreRequisites("GATKBaseRecalibration");

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {   
  
	&GATKBaseReCalibration($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});	
    }
}

if ($scriptParameter{'pGATKHaploTypeCaller'} > 0) { #Run GATK HaploTypeCaller

    &PrintToFileHandles(\@printFilehandles, "\nGATK HaplotypeCaller\n");

    &CheckBuildHumanGenomePreRequisites("GATKHaploTypeCaller");
    &CheckBuildDownLoadPreRequisites("GATKHaploTypeCaller");
    &CheckBuildPTCHSMetricPreRequisites("GATKHaploTypeCaller");

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {
	
	if ( (defined($parameter{ $scriptParameter{'familyID'} }{$sampleIDs[$sampleIDCounter]}{'GATKTargetPaddedBedIntervalLists'}{'buildFile'})) && ($parameter{ $scriptParameter{'familyID'} }{$sampleIDs[$sampleIDCounter]}{'GATKTargetPaddedBedIntervalLists'}{'buildFile'} eq 1) ){
	    
	    &BuildPTCHSMetricPreRequisites($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "GATKHaploTypeCaller");
	    last; #Will handle all build per sampleID within sbatch script
	}
    }
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {
	    
	&GATKHaploTypeCaller($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'});
    }
}

if ($scriptParameter{'pGATKGenoTypeGVCFs'} > 0) { #Run GATK GenoTypeGVCFs. Done per family

    &PrintToFileHandles(\@printFilehandles, "\nGATK GenoTypeGVCFs\n");

    &CheckBuildHumanGenomePreRequisites("GATKGenoTypeGVCFs");

    &GATKGenoTypeGVCFs($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH");
}

if ($scriptParameter{'pGATKVariantRecalibration'} > 0) { #Run GATK VariantRecalibrator/ApplyRecalibration. Done per family

    &PrintToFileHandles(\@printFilehandles, "\nGATK VariantRecalibrator/ApplyRecalibration\n");

    &CheckBuildHumanGenomePreRequisites("GATKVariantRecalibration");
    &CheckBuildDownLoadPreRequisites("GATKVariantRecalibration");
    &CheckBuildPTCHSMetricPreRequisites("GATKVariantRecalibration");

    &GATKVariantReCalibration($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH");
}

if ($scriptParameter{'pSampleCheck'} > 0) { #Run SampleCheck. Done per family

    &PrintToFileHandles(\@printFilehandles, "\nSampleCheck\n");

    &SampleCheck($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH");
}

if ($scriptParameter{'pGATKPhaseByTransmission'} > 0) { #Run GATK PhaseByTransmission. Done per family
    
    &PrintToFileHandles(\@printFilehandles, "\nGATK PhaseByTransmission\n");
    
    &CheckBuildHumanGenomePreRequisites("GATKPhaseByTransmission");
    &GATKPhaseByTransmission($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH");
}

if ($scriptParameter{'pGATKReadBackedPhasing'} > 0) { #Run GATK ReadBackedPhasing. Done per family. NOTE: Needs phased calls
    
    &PrintToFileHandles(\@printFilehandles, "\nGATK ReadBackedPhasing\n");
    
    &CheckBuildHumanGenomePreRequisites("GATKReadBackedPhasing");
    &GATKReadBackedPhasing($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH");
}

if ($scriptParameter{'pVariantEffectPredictor'} > 0) { #Run VariantEffectPredictor. Done per family

    &PrintToFileHandles(\@printFilehandles, "\nVariantEffectPredictor\n");

    &VariantEffectPredictor($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH");
}

if ($scriptParameter{'pVCFParser'} > 0) { #Run VariantEffectPredictor. Done per family

    &PrintToFileHandles(\@printFilehandles, "\nVCFParser\n");
    
    &VCFParser($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH");

    if ($scriptParameter{'vcfParserSelectFile'} ne "noUserInfo") {

	$VEPOutputFiles = 2; #Use seperate analysis for variants overlapping selected genes and orphans
    }
}

if ($scriptParameter{'pAnnovar'} > 0) { #Run Annovar. Done per family

    &PrintToFileHandles(\@printFilehandles, "\nAnnovar\n");

    &CheckBuildHumanGenomePreRequisites("Annovar");

    for (my $tableNamesCounter=0;$tableNamesCounter<scalar(@annovarTableNames);$tableNamesCounter++) { #For all specified table names

	if ($parameter{$annovarTableNames[$tableNamesCounter]}{'buildFile'} eq 1) {

	&BuildAnnovarPreRequisites($scriptParameter{'familyID'}, "mosaik", "Annovar");
	last; #Will handle all build tables within sbatch script
	}
    }
    &Annovar($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH");
}

if ($scriptParameter{'pSnpEff'} > 0) { #Run snpEff. Done per family

    &PrintToFileHandles(\@printFilehandles, "\nSnpEff\n");

    &SnpEff($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH");
}

if ($scriptParameter{'pGATKVariantEvalAll'} > 0) { #Run GATK VariantEval for all variants. Done per sampleID

    &PrintToFileHandles(\@printFilehandles, "\nGATK VariantEval All\n");

    &CheckBuildHumanGenomePreRequisites("GATKVariantEvalAll");

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) { 
	
	&GATKVariantEvalAll($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'}, "BOTH", $scriptParameter{'familyID'});
    }
}

if ($scriptParameter{'pMergeAnnotatedVariants'} > 0) { #Run MergeAnnotationVariants using intersectCollect.pl. Done per family

    &PrintToFileHandles(\@printFilehandles, "\nMergeAnnotatedVariants\n");

    &MergeAnnotatedVariants($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH");

}

if ($scriptParameter{'pGATKVariantEvalExome'} > 0) { #Run GATK VariantEval for exome variants. Done per sampleID

    &PrintToFileHandles(\@printFilehandles, "\nGATK VariantEval Exome\n");
    
    &CheckBuildHumanGenomePreRequisites("GATKVariantEvalExome");
    &CheckBuildDownLoadPreRequisites("GATKVariantEvalExome");

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) { 
	
	&GATKVariantEvalExome($sampleIDs[$sampleIDCounter], $scriptParameter{'aligner'}, "BOTH", $scriptParameter{'familyID'});
    }
}

if ($scriptParameter{'pAddDepth'} > 0) { #Run AddDepth using add_depth.pl. Done per family

    &PrintToFileHandles(\@printFilehandles, "\nAddDepth\n");

    &AddDp($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH");
}

if ($scriptParameter{'pRankVariants'} > 0) { #Run RankVariants. Done per family

    &PrintToFileHandles(\@printFilehandles, "\nRankVariants\n");

    &RankVariants($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH");
}

if ($scriptParameter{'pQCCollect'} > 0) { #Run QCCollect. Done per family

    &PrintToFileHandles(\@printFilehandles, "\nQCCollect\n");

    &QCCollect($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH");
}

if ($scriptParameter{'pRemovalRedundantFiles'} > 0) { #Sbatch generation of removal of alignment files
    
    &PrintToFileHandles(\@printFilehandles, "\nRemoval of alignment files\n");

    &RemoveRedundantFiles($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH");	
}

if ( ($scriptParameter{'pAnalysisRunStatus'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {

    $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'AnalysisRunStatus'} = "notFinished"; #Add analysis run status flag.
}

if ($scriptParameter{'pAnalysisRunStatus'} > 0) {

    &PrintToFileHandles(\@printFilehandles, "\nAnalysisRunStatus\n");

    &AnalysisRunStatus($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH");
}

close(MIPLOG); #Close mip_log file

#Write QC for programs used in analysis                                                                                                                         
if ($scriptParameter{'sampleInfoFile'} ne 0) {#Write SampleInfo to yaml file
    
    &WriteYAML($scriptParameter{'sampleInfoFile'}, \%sampleInfo); #Write QC for sampleinfo used in analysis
}


######################
####SubRoutines#######
######################

sub AnalysisRunStatus { 
###execute last in MAIN chain and sets analysis run status flag to finished.

    my $familyID = $_[0]; #familyID NOTE: not sampleid 
    my $aligner = $_[1];
    my $callType = $_[2]; #SNV,INDEL or BOTH

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle

    &ProgramPreRequisites($familyID, "AnalysisRunStatus", "analysisrunstatus", 0, $FILEHANDLE, 1, 1);

    print $FILEHANDLE q?perl -i -p -e 'if($_=~/AnalysisRunStatus\:/) { s/notFinished/finished/g }' ?.$scriptParameter{'sampleInfoFile'}.q? ?, "\n\n";  
    
    close($FILEHANDLE); 
    
    if ( ($scriptParameter{'pAnalysisRunStatus'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {

	&FIDSubmitJob(0, $familyID, 2, $parameter{'pAnalysisRunStatus'}{'chain'}, $fileName, 0);
    }
    return;
}

sub RemoveRedundantFiles {
#Generates a sbatch script, which removes redundant files.
    
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
		
		&CheckMostCompleteAndRemoveFile($FILEHANDLE, \$sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'}, \($inSampleDirectory."/".$infile.".bam"), ".bam");	    
	    }
##Remove BWA files
	    if ($scriptParameter{'pBwaAln'} > 0) {
		
		print $FILEHANDLE "rm ";
		print $FILEHANDLE $inSampleDirectory."/".$infile.".sai", "\n\n"; #BWA_Aln
	    }
	    if ($scriptParameter{'pBwaSampe'} >0) {
		
		&CheckMostCompleteAndRemoveFile($FILEHANDLE, \$sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'}, \($inSampleDirectory."/".$infile.".bam"), ".bam");
	    }  
	    if ($scriptParameter{'pBwaMem'} >0) {
	
		&CheckMostCompleteAndRemoveFile($FILEHANDLE, \$sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'}, \($inSampleDirectory."/".$infile.".bam"), ".bam");
	    }    	    
##Sorted BAM
	    if ($scriptParameter{'pPicardToolsSortSam'} > 0) {
		
		my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pPicardToolsSortSam'}{'fileEnding'};
		
		&CheckMostCompleteAndRemoveFile($FILEHANDLE, \$sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'}, \($inSampleDirectory."/".$infile.$outfileEnding.".bam"), ".bam"); #Sorted BAM and bai file
	    }
	}
	
##Potentially merged files
	my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles($sampleID);        
	
	if ($PicardToolsMergeSwitch == 1) { #Files was merged previously
	    
	    if ($scriptParameter{'pPicardToolsMergeSamFiles'} > 0) {
		
		my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pPicardToolsMergeSamFiles'}{'fileEnding'};
		
		&CheckMostCompleteAndRemoveFile($FILEHANDLE, \$sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'}, \($inSampleDirectory."/".$infile.$outfileEnding.".bam"), ".bam"); #merged BAM and bai file
	    }	
	    if ($scriptParameter{'pPicardToolsMarkduplicates'} > 0) {
		
		my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pPicardToolsMarkduplicates'}{'fileEnding'};
		
		&CheckMostCompleteAndRemoveFile($FILEHANDLE, \$sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'}, \($inSampleDirectory."/".$infile.$outfileEnding.".bam"), ".bam"); #Dedupped BAM and bai file
	    }
	    if ($scriptParameter{'pGATKRealigner'} > 0) {
		
		my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/GATK";   
		my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pGATKRealigner'}{'fileEnding'};
	   
		&CheckMostCompleteAndRemoveFile($FILEHANDLE, \$sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'}, \($inSampleDirectory."/".$infile.$outfileEnding.".bam"), ".bam"); #ReAligned BAM and bai file
	    }
	    if ($scriptParameter{'pGATKBaseRecalibration'} > 0) {
		
		my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/GATK";
		my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pGATKBaseRecalibration'}{'fileEnding'};

		&CheckMostCompleteAndRemoveFile($FILEHANDLE, \$sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'}, \($inSampleDirectory."/".$infile.$outfileEnding.".bam"), ".bam"); #BaseRecalibrated BAM and bai file
	    }
	    if ($scriptParameter{'pGATKHaploTypeCaller'} > 0) { #Always collapses all files even if there is only one
		
		my $lanes = join("",@{$lane{$sampleID}}); #Extract lanes
		my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/GATK";
		my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pGATKHaploTypeCaller'}{'fileEnding'};
		
		&CheckMostCompleteAndRemoveFile($FILEHANDLE, \$sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'MostCompleteVCF'}{'Path'}, \($inSampleDirectory."/".$infile.$outfileEnding.".vcf"), ".vcf"); #HaplotypeCaller gvcf file
	    }
	}
	else {
	    
	    for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #For all infiles per lane
		
		my $infile = $infilesLaneNoEnding{$sampleID}[$infileCounter];
		
		if ($scriptParameter{'pPicardToolsMarkduplicates'} > 0) {
		    
		    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pPicardToolsMarkduplicates'}{'fileEnding'};

		    &CheckMostCompleteAndRemoveFile($FILEHANDLE, \$sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'}, \($inSampleDirectory."/".$infile.$outfileEnding.".bam"), ".bam"); #Dedupped BAM and bai file
		}
		if ($scriptParameter{'pGATKRealigner'} > 0) {
		    
		    my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/GATK";
		    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pGATKRealigner'}{'fileEnding'};
		    
		    &CheckMostCompleteAndRemoveFile($FILEHANDLE, \$sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'}, \($inSampleDirectory."/".$infile.$outfileEnding.".bam"), ".bam"); #ReAligned BAM and bai file
		}
		if ($scriptParameter{'pGATKBaseRecalibration'} > 0) {
		    
		    my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/GATK";
		    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pGATKBaseRecalibration'}{'fileEnding'};

		    &CheckMostCompleteAndRemoveFile($FILEHANDLE, \$sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'}, \($inSampleDirectory."/".$infile.$outfileEnding.".bam"), ".bam"); #BaseRecalibrated BAM and bai file
		}
		if ($scriptParameter{'pGATKHaploTypeCaller'} > 0) { #Always collapses all files even if there is only one
		    
		    my $lanes = join("",@{$lane{$sampleID}}); #Extract lanes
		    my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/GATK";
		    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pGATKHaploTypeCaller'}{'fileEnding'};
		    
		    &CheckMostCompleteAndRemoveFile($FILEHANDLE, \$sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'MostCompleteVCF'}{'Path'}, \($inSampleDirectory."/".$sampleID."_lanes_".$lanes.$outfileEnding.".vcf"), ".vcf"); #HaplotypeCaller gvcf file
		}
	    }
	}
    }
###Family files
    if ($scriptParameter{'pGATKGenoTypeGVCFs'} > 0) {

	my $outFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK"; #New outfile directory
	my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pGATKGenoTypeGVCFs'}{'fileEnding'};
	
	&CheckMostCompleteAndRemoveFile($FILEHANDLE, \$sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'MostCompleteVCF'}{'Path'}, \($outFamilyDirectory."/".$familyID.$outfileEnding.$callType.".vcf"), ".vcf"); #GATKGenoTypeGVCFs vcf file
	
    }
    if ($scriptParameter{'pGATKVariantRecalibration'} > 0) {
	
	my $outFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK"; #New outfile directory
	my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pGATKVariantRecalibration'}{'fileEnding'};
	
	&CheckMostCompleteAndRemoveFile($FILEHANDLE, \$sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'MostCompleteVCF'}{'Path'}, \($outFamilyDirectory."/".$familyID.$outfileEnding.$callType.".vcf"), ".vcf"); #pGATKVariantRecalibration vcf file
    
    }
    if ($scriptParameter{'pAnnovar'} > 0) {
	
	my $outFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
	my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pAnnovar'}{'fileEnding'};
	
	print $FILEHANDLE "rm ";
	print $FILEHANDLE $outFamilyDirectory."/".$familyID.$outfileEnding.$callType."*".$scriptParameter{'annovarGenomeBuildVersion'}."_*", "\n\n"; #Annovar data files
	print $FILEHANDLE "rm ";
	print $FILEHANDLE $outFamilyDirectory."/".$familyID.$outfileEnding.$callType.".*", "\n\n"; #Annovar data files
    }  
    close($FILEHANDLE);
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
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pGATKVariantRecalibration'}{'fileEnding'};

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
	&FIDSubmitJob(0, $familyID, 1, $parameter{'pQCCollect'}{'chain'}, $fileName, 0);
    }
}

sub RankVariants { 
###Filter and Rank variants depending on mendelian inheritance, frequency and phenotype
   
    my $familyID = $_[0]; #familyID NOTE: not sampleid 
    my $aligner = $_[1];
    my $callType = $_[2]; #SNV,INDEL or BOTH

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
 
    &ProgramPreRequisites($familyID, "RankVariants", $aligner."/GATK/candidates/ranking", $callType, $FILEHANDLE, 1, 4);

    my $inFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $outFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK/candidates/ranking";
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pSnpEff'}{'fileEnding'};
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pRankVariants'}{'fileEnding'};
    my $analysisType = "";

###Ranking
    print $FILEHANDLE "#Ranking", "\n";
    
    print $FILEHANDLE "workon ".$scriptParameter{'pythonVirtualEnvironment'}, "\n\n"; #Activate python environment
    
    for (my $VEPOutputFilesCounter=0;$VEPOutputFilesCounter<$VEPOutputFiles;$VEPOutputFilesCounter++) {
	
	if ($VEPOutputFilesCounter == 1) {
	    
	    $analysisType = ".selected"; #SelectFile variants
	}
	    
	print $FILEHANDLE "run_genmod.py ";
	print $FILEHANDLE $scriptParameter{'pedigreeFile'}." "; #Pedigree file
	
	if ($scriptParameter{'pVariantEffectPredictor'} > 0) { #Use VEP annotations in compound models

	    print $FILEHANDLE "-vep "; 
	}
	else {
	    
	    print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$scriptParameter{'geneFile'}." "; #Gene file used for annotating AR_compounds
	}
	if ($scriptParameter{'instanceTag'} eq "CMMS") {
	    
	    print $FILEHANDLE "-family cmms "; #CMMS flag
	}
	if ($scriptParameter{'caddWGSSNVs'} == 1) {
	 
	    print $FILEHANDLE "-cadd ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'caddWGSSNVsFile'}." "; #Whole genome sequencing CADD score file
	}
	if ($scriptParameter{'cadd1000Genomes'} == 1) {
	 
	    print $FILEHANDLE "-c1kg ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'cadd1000GenomesFile'}." "; #1000G CADD score file
	}
	if ($scriptParameter{'wholeGene'} == 1) {
	 
	    print $FILEHANDLE "--whole_gene "; 
	}
	#print $FILEHANDLE "-tres ".$scriptParameter{'rankScore'}." "; #Rank score threshold
	print $FILEHANDLE $inFamilyDirectory."/".$familyID.$infileEnding.$callType.$analysisType.".vcf "; #InFile
	print $FILEHANDLE "-o ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.$analysisType.".vcf &", "\n\n"; #OutFile
	
	if ( ($scriptParameter{'pRankVariants'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	    
	    if ($VEPOutputFilesCounter == 1) {
	
		$sampleInfo{$familyID}{$familyID}{'program'}{'RankVariants'}{'Clinical'}{'Path'} = $outFamilyDirectory."/".$familyID.$outfileEnding.$callType.$analysisType.".txt";  #Save clinical candidate list path
	    }
	    else {

		$sampleInfo{$familyID}{$familyID}{'program'}{'RankVariants'}{'Research'}{'Path'} = $outFamilyDirectory."/".$familyID.$outfileEnding.$callType.$analysisType.".txt";  #Save research candidate list path
	    }
	}
    }
    print $FILEHANDLE "wait\n\n";
    print $FILEHANDLE "\n\ndeactivate ", "\n\n"; #Deactivate python environment

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
    my $nrCores = 1;

    &ProgramPreRequisites($familyID, "AddDepth", $aligner."/GATK", $callType, $FILEHANDLE, $nrCores, 10);
    
    my $inFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $outFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pMergeAnnotatedVariants'}{'fileEnding'};
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pAddDepth'}{'fileEnding'};
    my $coreCounter=1;
    my $analysisType = "";

    for (my $VEPOutputFilesCounter=0;$VEPOutputFilesCounter<$VEPOutputFiles;$VEPOutputFilesCounter++) {

	if ($VEPOutputFiles == 2) {

	    $analysisType = ".selected"; #SelectFile variants
	}
#Find all "./." per sample ID and print chr pos to new file (mpileup -l format)
	for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) { #For all sample ids, find nonvariants
	    
	    &PrintWait(\$sampleIDCounter, \$nrCores, \$coreCounter, $FILEHANDLE);
	    
	    print $FILEHANDLE "#Find all './.' per sampleID and print chrosome position to new file (mpileup -l format)", "\n";
	    
	    print $FILEHANDLE q?perl -F'\t' -nae' if ($_=~ /?.$sampleIDs[$sampleIDCounter].q?\S+\.\/\./ ) { print "$F[0] $F[1]","\n"; }' ?; #print chromosome and start for sampleID
	    print $FILEHANDLE $inFamilyDirectory."/".$familyID.$infileEnding.$callType.$analysisType.".txt "; #InFile
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
	print $FILEHANDLE "-i ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType.$analysisType.".txt "; #InFile
	
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
	print $FILEHANDLE "-o ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType.$analysisType."_adp.txt", "\n\n"; 

	print $FILEHANDLE "mv "; #Overwrites original annovar_merge.txt file
	print $FILEHANDLE $inFamilyDirectory."/".$familyID.$infileEnding.$callType.$analysisType."_adp.txt "; #Add_depth outfile
	print $FILEHANDLE $inFamilyDirectory."/".$familyID.$infileEnding.$callType.$analysisType.".txt", "\n\n"; #Original file
    }
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
    my $time = 4;
    my $numberOfCores = 1; #Set the number of cores depending on exome/rapid or WGS

    if ($scriptParameter{'analysisType'} eq "genomes") { #WGS analysis
	
	$time = 30;
	$numberOfCores = 6; 
    }

    &ProgramPreRequisites($familyID, "MergeAnnotatedVariants", $aligner."/GATK", $callType, $FILEHANDLE, $numberOfCores, $time);

    my $inFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $outFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pAnnovar'}{'fileEnding'};
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pMergeAnnotatedVariants'}{'fileEnding'};
    my @additionalVCFColumns;
    my $vcfColumns;
    my $sampleIDcolumnsCondition;
    my $analysisType = "";

    ##Determine number of VEP columns
    if ($scriptParameter{'pVariantEffectPredictor'} == 1) {

	my $AddtionalVEPOutColumnCounter = 3; #HGNC_transcript_info, Functional_annotation, Gene_annotation (constant)
	splice(@additionalVCFColumns, scalar(@additionalVCFColumns), 0, "HGNC_transcript_info", "Functional_annotation", "GeneticRegionAnnotation");

	for (my $vepFeatureCounter=0;$vepFeatureCounter<scalar(@vepFeatures);$vepFeatureCounter++) {
	 
	    if ($vepFeatures[$vepFeatureCounter] eq "sift") {
		
		$AddtionalVEPOutColumnCounter++;
		splice(@additionalVCFColumns, scalar(@additionalVCFColumns), 0, "Sift");
	    }
	    if ($vepFeatures[$vepFeatureCounter] eq "polyphen") {

		$AddtionalVEPOutColumnCounter++;
		splice(@additionalVCFColumns, scalar(@additionalVCFColumns), 0, "PolyPhen");
	    }
	}
	$vcfColumns = scalar(@sampleIDs)+6+$AddtionalVEPOutColumnCounter; #Requires CMMS format (chr,start,stop,ref_allele,alt_allel,GT_Call_Filter,VEPColumns,IDN...)
    }
    else {

	$vcfColumns = scalar(@sampleIDs)+6; #Requires CMMS format (chr,start,stop,ref_allele,alt_allel,GT_Call_Filter,IDN...)
    }    
    $sampleIDcolumnsCondition = $vcfColumns - scalar(@sampleIDs);
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

    for (my $VEPOutputFilesCounter=0;$VEPOutputFilesCounter<$VEPOutputFiles;$VEPOutputFilesCounter++) {

##Create family specific template
	print $FILEHANDLE "#Create db master file from template", "\n";
	print $FILEHANDLE q?perl -nae '?;
##Add additional columns
	my $columnHeaders = join(",", @additionalVCFColumns);
	print $FILEHANDLE q?my @additionalVCFColumns = (?.$columnHeaders;
	print $FILEHANDLE q?); my $additionalColumnCounter; ?;
	
	print $FILEHANDLE q?if ($_=~/outinfo/i) { if ($_=~/IDN!/) { my $columnString; my $sampleIDString; for (my $columnCounter=6;$columnCounter<?.$vcfColumns.q?;$columnCounter++) { if ($columnCounter<?.$sampleIDcolumnsCondition.q?) { $columnString.="$additionalVCFColumns[$additionalColumnCounter]=>0_$columnCounter,"; $additionalColumnCounter++;} elsif ($columnCounter ==?.$vcfColumns.q?-1) { $sampleIDString.="IDN_GT_Call=>0_$columnCounter";} else { $sampleIDString.="IDN_GT_Call=>0_$columnCounter,";} } $columnString .= $sampleIDString ;s/IDN!/$columnString/g; print $_; } next;} elsif ($_=~s/FDN\!_|FDN\!/?.$familyID.q?/g) { if($_=~s/^ODF!/?.$regExpOutDataFile.q?/g) {} if($_=~s/ALIGNER!/?.$aligner.q?/g) {} if($_=~s/FILEENDING!_/?.$infileEnding.q?/g) {} if($_=~s/CALLTYPE!/?.$callType.q?/g) {} if ($_=~/IDN!/) { my $columnString; for (my $columnCounter=6;$columnCounter<?.$vcfColumns.q?;$columnCounter++) { if ($columnCounter<?.$vcfColumns.q?-1) { $columnString.="$columnCounter,"} else { $columnString.="$columnCounter\t"} } s/IDN!/$columnString/g; print $_;} else { print $_;} } else { if($_=~s/^RD!/?.$regExpReferenceDirectory.q?/g) {} print $_;}' ?;

	##Include potential SelectFile variants
	if ($VEPOutputFilesCounter == 1) {
	    
	    print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$scriptParameter{'mergeAnnotatedVariantsTemplateSelectFile'}." "; #Infile
	    print $FILEHANDLE "> ".$scriptParameter{'mergeAnnotatedVariantsDbSelectFile'}, "\n\n"; #OutFile
	}
	else { #Orphan analysis
	 
	    print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$scriptParameter{'mergeAnnotatedVariantsTemplateFile'}." "; #Infile
	    print $FILEHANDLE "> ".$scriptParameter{'mergeAnnotatedVariantsDbFile'}, "\n\n"; #OutFile
	}
	print $FILEHANDLE "perl ".$scriptParameter{'inScriptDir'}."/intersectCollect.pl ";

	##Include potential SelectFile variants
	if ($VEPOutputFilesCounter == 1) {

	    $analysisType = ".selected"; #SelectFile variants
	    print $FILEHANDLE "-db ".$scriptParameter{'mergeAnnotatedVariantsDbSelectFile'}." ";
	}
	else {

	    print $FILEHANDLE "-db ".$scriptParameter{'mergeAnnotatedVariantsDbFile'}." ";
	}
	if ($humanGenomeReferenceSource eq "hg19") {
	    
	    print $FILEHANDLE "-prechr 1 "; #Use chromosome prefix
	}
	
	print $FILEHANDLE "-o ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.$analysisType.".txt ";
	print $FILEHANDLE "-s ";

	for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) { 
	    
	    if ($sampleIDCounter eq scalar(@sampleIDs)-1) {
		
		print $FILEHANDLE $sampleIDs[$sampleIDCounter], "\n\n";
	    }
	    else {
		
		print $FILEHANDLE $sampleIDs[$sampleIDCounter], ",";
	    }    
	}
    }
    close($FILEHANDLE);   
    
    if ( ($scriptParameter{'pMergeAnnotatedVariants'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {

	&FIDSubmitJob(0, $familyID, 1, $parameter{'pMergeAnnotatedVariants'}{'chain'}, $fileName, 0);
    }
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
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pGATKVariantRecalibration'}{'fileEnding'};
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pGATKVariantRecalibration'}{'fileEnding'};
    my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles($sampleID);
    
    if ($PicardToolsMergeSwitch == 1) { #Files was merged previously
    
###GATK SelectVariants
##Select SampleID from familyID vrecal vcf file
	print $FILEHANDLE "#GATK SelectVariants","\n\n";
	print $FILEHANDLE "java -Xmx2g ";
	
	if ($scriptParameter{'javaUseLargePages'} ne "no") {
	    
	    print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	}
	print $FILEHANDLE "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	print $FILEHANDLE "-l INFO "; #Set the minimum level of logging
	print $FILEHANDLE "-T SelectVariants "; #Type of analysis to run
	print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file
	print $FILEHANDLE "-V: ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf "; #FamilyID inFile
	print $FILEHANDLE "-o ".$outSampleDirectory."/".$infile.$outfileEnding.$callType."_temp.vcf "; #SampleID exome outFile
	print $FILEHANDLE "-sn ".$sampleID, "\n\n"; #Include genotypes from this sample

##Prepp infile to only contain exonic variants

	my $sampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/GATK/varianteval";
	my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pSnpEff'}{'fileEnding'};

	print $FILEHANDLE "grep exon ";
	print $FILEHANDLE $inFamilyDirectory."/".$familyID.$infileEnding.$callType.".txt "; #InFile
	print $FILEHANDLE "> ".$sampleDirectory."/".$sampleID.$infileEnding.$callType."_exonic_variants.txt", "\n\n"; #OutFile

	##Include potential SelectFile variants
	if ($VEPOutputFiles == 2) {
	    
	    my $analysisType = ".selected"; #SelectFile variants
	    print $FILEHANDLE "perl -ne ' if ( ($_=~/exonic/) || ($_=/splicing/) ) {print $_;}' ";
	    print $FILEHANDLE $inFamilyDirectory."/".$familyID.$infileEnding.$callType.$analysisType.".txt "; #InFile
	    print $FILEHANDLE "> ".$sampleDirectory."/".$sampleID.$infileEnding.$callType."_exonic_variants".$analysisType.".txt", "\n\n"; #OutFile
	    
	    #Merge orphans and selectfiles
	    print $FILEHANDLE "cat ";
	    print $FILEHANDLE $sampleDirectory."/".$sampleID.$infileEnding.$callType."_exonic_variants.txt "; #Orphan file
	    print $FILEHANDLE $sampleDirectory."/".$sampleID.$infileEnding.$callType."_exonic_variants".$analysisType.".txt "; #SelectFile variants
	    print $FILEHANDLE "> ".$sampleDirectory."/".$sampleID.$infileEnding.$callType."_exonic_variants_combined.txt", "\n\n"; #OutFile
	    
	    #Sort combined file
	    print $FILEHANDLE "sort ";
	    print $FILEHANDLE "-k1,1 -k2,2n "; #Numerically by chromosome and start position
	    print $FILEHANDLE $sampleDirectory."/".$sampleID.$infileEnding.$callType."_exonic_variants_combined.txt ";
	    print $FILEHANDLE "> ".$sampleDirectory."/".$sampleID.$infileEnding.$callType."_exonic_variants.txt", "\n\n"; #OutFile
	}

##Intersect exonic variants from created sampleID vcf file (required for GATKVariantEval for exonic variants)
	print $FILEHANDLE "intersectBed ";
	print $FILEHANDLE "-header "; #Print the header from the A file prior to results.
	print $FILEHANDLE "-a ".$sampleDirectory."/".$infile.$outfileEnding.$callType."_temp.vcf "; #SampleID temp exome vcf inFile
	print $FILEHANDLE "-b ".$sampleDirectory."/".$sampleID.$infileEnding.$callType."_exonic_variants.txt "; #SampleID exonic variants
	print $FILEHANDLE "> ".$sampleDirectory."/".$infile.$outfileEnding.$callType."_exome.vcf", "\n\n"; #OutFile (VCF-format)

###VariantEval (exome variants)
	print $FILEHANDLE "#GATK VariantEval","\n\n";
	
	print $FILEHANDLE "java -Xmx2g ";
	
	if ($scriptParameter{'javaUseLargePages'} ne "no") {
	    
	    print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	}
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

	    if ($scriptParameter{'javaUseLargePages'} ne "no") {
		
		print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	    }
	    print $FILEHANDLE "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	    print $FILEHANDLE "-l INFO "; #Set the minimum level of logging
	    print $FILEHANDLE "-T SelectVariants "; #Type of analysis to run
	    print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file
	    print $FILEHANDLE "-V: ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf "; #FamilyID infile 
	    print $FILEHANDLE "-o ".$outSampleDirectory."/".$infile.$outfileEnding.$callType."_temp.vcf "; #SampleID outFile
	    print $FILEHANDLE "-sn ".$sampleID, "\n\n"; #Include genotypes from this sample

##Prepp infile to only contain exonic variants

	    my $sampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/GATK/varianteval";
	    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pSnpEff'}{'fileEnding'};
	    
	    print $FILEHANDLE q?perl -ne ' if ( ($_=~/exonic/) || ($_=/splicing/) ) {print $_;}' ?;
	    print $FILEHANDLE $inFamilyDirectory."/".$familyID.$infileEnding.$callType.".txt "; #InFile
	    print $FILEHANDLE "> ".$sampleDirectory."/".$sampleID.$infileEnding.$callType."_exonic_variants.txt", "\n\n"; #OutFile

	    ##Include potential SelectFile variants
	    if ($VEPOutputFiles == 2) {

		my $analysisType = ".selected"; #SelectFile variants
		print $FILEHANDLE q?perl -ne ' if ( ($_=~/exonic/) || ($_=/splicing/) ) {print $_;}' ?;
		print $FILEHANDLE $inFamilyDirectory."/".$familyID.$infileEnding.$callType.$analysisType.".txt "; #InFile
		print $FILEHANDLE "> ".$sampleDirectory."/".$sampleID.$infileEnding.$callType."_exonic_variants".$analysisType.".txt", "\n\n"; #OutFile
		
		#Merge orphans and selectfiles
		print $FILEHANDLE "cat ";
		print $FILEHANDLE $sampleDirectory."/".$sampleID.$infileEnding.$callType."_exonic_variants.txt "; #Orphan file
		print $FILEHANDLE $sampleDirectory."/".$sampleID.$infileEnding.$callType."_exonic_variants".$analysisType.".txt "; #SelectFile variants
		print $FILEHANDLE "> ".$sampleDirectory."/".$sampleID.$infileEnding.$callType."_exonic_variants_combined.txt", "\n\n"; #OutFile

		#Sort combined file
		print $FILEHANDLE "sort ";
		print $FILEHANDLE "-k1,1 -k2,2n "; #Numerically by chromosome and start position
		print $FILEHANDLE $sampleDirectory."/".$sampleID.$infileEnding.$callType."_exonic_variants_combined.txt ";
		print $FILEHANDLE "> ".$sampleDirectory."/".$sampleID.$infileEnding.$callType."_exonic_variants.txt", "\n\n"; #OutFile
	    }

##Intersect exonic variants from created sampleID vcf file (required for GATKVariantEval for exonic variants)
	    print $FILEHANDLE "intersectBed ";
	    print $FILEHANDLE "-header "; #Print the header from the A file prior to results.
	    print $FILEHANDLE "-a ".$sampleDirectory."/".$infile.$outfileEnding.$callType."_temp.vcf "; #SampleID temp exome vcf inFile
	    print $FILEHANDLE "-b ".$sampleDirectory."/".$sampleID.$infileEnding.$callType."_exonic_variants.txt "; #SampleID exonic variants
	    print $FILEHANDLE "> ".$sampleDirectory."/".$infile.$outfileEnding.$callType."_exome.vcf", "\n\n"; #OutFile (VCF-format)
	    
###VariantEval (exome variants)
	    
	    print $FILEHANDLE "#GATK VariantEval","\n\n";
	    
	    print $FILEHANDLE "java -Xmx2g ";
	 
	    if ($scriptParameter{'javaUseLargePages'} ne "no") {
		
		print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	    }
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

	    if ($VEPOutputFiles == 2) { #Selected analysis has been performed
	
		print $FILEHANDLE "rm ";
		print $FILEHANDLE $sampleDirectory."/".$sampleID.$infileEnding.$callType."_exonic_variants_combined.txt", "\n\n"; #Combined Selceted and orphan file
	    }
	    
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
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pGATKVariantRecalibration'}{'fileEnding'};
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pGATKVariantRecalibration'}{'fileEnding'};
    my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles($sampleID);
    
    if ($PicardToolsMergeSwitch == 1) { #Files was merged previously
    
###GATK SelectVariants
##Select SampleID from familyID vrecal vcf file
	print $FILEHANDLE "#GATK SelectVariants","\n\n";
	print $FILEHANDLE "java -Xmx2g ";

	if ($scriptParameter{'javaUseLargePages'} ne "no") {
	    
	    print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	}
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

	if ($scriptParameter{'javaUseLargePages'} ne "no") {
	    
	    print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	}
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

	    if ($scriptParameter{'javaUseLargePages'} ne "no") {
		
		print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	    }
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

	    if ($scriptParameter{'javaUseLargePages'} ne "no") {
		
		print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	    }

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
}

sub Annovar { 
###Filter SNVs by gene, region and databases

    my $familyID = $_[0]; #familyID NOTE: not sampleid 
    my $aligner = $_[1];
    my $callType = $_[2]; #SNV,INDEL or BOTH 

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $nrCores = &NrofCoresPerSbatch(scalar(@annovarTableNames)); #Detect the number of cores to use from @annovarTableNames. 

    &ProgramPreRequisites( $familyID, "Annovar", $aligner."/GATK", $callType, $FILEHANDLE, $nrCores, 7);

    my $inFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $outFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pAnnovar'}{'fileEnding'};
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pVCFParser'}{'fileEnding'};
    my $analysisType = "";

    for (my $VEPOutputFilesCounter=0;$VEPOutputFilesCounter<$VEPOutputFiles;$VEPOutputFilesCounter++) {

	if ($VEPOutputFilesCounter == 1) {

	    $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pVCFParser'}{'fileEnding'};
	    $analysisType = ".selected"; #SelectFile variants
	}
	
	print $FILEHANDLE "#Prepare infile to Annovar format from GATK vcf4", "\n";
	print $FILEHANDLE "perl ".$scriptParameter{'annovarPath'}."/convert2annovar.pl ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType.$analysisType.".vcf ";
	print $FILEHANDLE "-format vcf4old "; #the format of the input file
	print $FILEHANDLE "-includeinfo "; #specify that the output should contain additional information in the input line
	print $FILEHANDLE "> ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.$analysisType."_temp", "\n\n"; #Annovar script
    
	print $FILEHANDLE "#Intersect for all samples within familyid and remake file to fit annovar format and subsequent filtering", "\n";
	print $FILEHANDLE q?perl -nae 'my @format; my $formatInfo;chomp($_); if ($_=~/^#/) {print $_;next;} if ($_=~/;set=2/) {} else{ if($F[11] eq "PASS") {} else {$F[11] = "PRES";} @format = split(":",$F[13]); print $F[0], "\t", $F[1], "\t", $F[2], "\t", $F[3], "\t", $F[4], "\t", $F[11], "\t"; ?;

	##Include any addtional columns from VEPParse
	my $numberofFields = 14+scalar(@sampleIDs); #First sample genotype starts at col 14 (start 0, perl).
	print $FILEHANDLE q?for (my $fieldsCounter=?.$numberofFields.q?;$fieldsCounter<scalar(@F);$fieldsCounter++) {print @F[$fieldsCounter], "\t";} ?;
    
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

	print $FILEHANDLE "' ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.$analysisType."_temp "; #InFile from just created convert2annovar.pl outfile
	print $FILEHANDLE "> ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.$analysisType, "\n\n"; #OutFile
	
	print $FILEHANDLE "#Perform sort since Annovar is not garantied to produce a numerically sorted outfile", "\n";
	print $FILEHANDLE "sort ";
	print $FILEHANDLE "-k1,1 -k2,2n "; #Numerically by chromosome and start position
	print $FILEHANDLE $outFamilyDirectory."/".$familyID.$outfileEnding.$callType.$analysisType." "; #Infile
	print $FILEHANDLE "> ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.$analysisType."_sorted ", "\n\n"; #Outfile

	print $FILEHANDLE "mv ";
	print $FILEHANDLE $outFamilyDirectory."/".$familyID.$outfileEnding.$callType.$analysisType."_sorted "; #Infile
	print $FILEHANDLE $outFamilyDirectory."/".$familyID.$outfileEnding.$callType.$analysisType, "\n\n"; #Outfile

	$infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pAnnovar'}{'fileEnding'};
	my $coreCounter=1;   	    

	for (my $tableNamesCounter=0;$tableNamesCounter<scalar(@annovarTableNames);$tableNamesCounter++) { #For all specified table names

	    &PrintWait(\$tableNamesCounter, \$nrCores, \$coreCounter, $FILEHANDLE);	
	
	    print $FILEHANDLE "perl ".$scriptParameter{'annovarPath'}."/annotate_variation.pl "; #Annovar script 
	    print $FILEHANDLE "-".$annovarTables{$annovarTableNames[$tableNamesCounter]}{'annotation'}." "; #Annotation option	

	    if ($annovarTables{$annovarTableNames[$tableNamesCounter]}{'annotation'} eq "geneanno" ) { #Use hgvs output style
		
		print $FILEHANDLE "-hgvs ";
		print $FILEHANDLE "-exonicsplicing "; #Annotate variants near intron/exonic borders
	    }
	    print $FILEHANDLE "-buildver ".$scriptParameter{'annovarGenomeBuildVersion'}." ";
	
	    if($annovarTables{$annovarTableNames[$tableNamesCounter]}{'dbtype'} eq "generic") {
	    
		print $FILEHANDLE "-dbtype generic -genericdbfile ".$annovarTables{$annovarTableNames[$tableNamesCounter]}{'file'}[0]." "; #generic db file
		print $FILEHANDLE "--outfile ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.$analysisType."_".$annovarTables{$annovarTableNames[$tableNamesCounter]}{'file'}[0]." "; #OutFile
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
	    print $FILEHANDLE $inFamilyDirectory."/".$familyID.$infileEnding.$callType.$analysisType." "; #Infile. Outfile is named using infile prefix except for generic files 
	    print $FILEHANDLE $scriptParameter{'annovarPath'}."/humandb &", "\n\n"; #annovar/humandb directory is assumed

	    if ($annovarTableNames[$tableNamesCounter] =~/ensGene|refGene/) { #Extra round to catch MT for refSeq as well
	    
		print $FILEHANDLE "grep MT "; #Only MT variants
		print $FILEHANDLE $inFamilyDirectory."/".$familyID.$infileEnding.$callType.$analysisType." "; #Infile.
		print $FILEHANDLE "> ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType.$analysisType.".GRCh37_MT"." ", "\n\n"; #outfile taht can be empty for exomes or MT for WGS

		print $FILEHANDLE "perl ".$scriptParameter{'annovarPath'}."/annotate_variation.pl "; #Annovar script 
		print $FILEHANDLE "-".$annovarTables{$annovarTableNames[$tableNamesCounter]}{'annotation'}." "; #Annotation option
		print $FILEHANDLE "-hgvs ";
		print $FILEHANDLE "-exonicsplicing "; #Annotate variants near intron/exonic borders
		print $FILEHANDLE "-buildver GRCh37_MT ";
		print $FILEHANDLE "-dbtype ensGene "; #db file. NOTE: RefSeq does not have mitochondria gene definition. So ANNOVAR use either UCSC Known Gene or Ensembl Gene.
		print $FILEHANDLE $inFamilyDirectory."/".$familyID.$infileEnding.$callType.$analysisType.".GRCh37_MT "; #Infile.
		print $FILEHANDLE "--outfile ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.$analysisType.".GRCh37_MT"." "; #OutFile prefix
		print $FILEHANDLE $scriptParameter{'annovarPath'}."/humandb &", "\n\n"; #annovar/humandb directory is assumed
		$nrCores--; #Reduce to make sure that print statement comes at correct interval since two calls are made for 1 annovar table
	    }
	}
	print $FILEHANDLE "wait", "\n\n";
    
	for (my $tableNamesCounter=0;$tableNamesCounter<scalar(@annovarTableNames);$tableNamesCounter++) { #For all specified table names
	
	    if ($annovarTableNames[$tableNamesCounter] =~/ensGene|refGene/) { #Extra round to concatenate MT to "variant function" and "exonic.variant" function
	    
		my @files = ("variant_function", "exonic_variant_function");
		
		for (my $fileCounter=0;$fileCounter<scalar(@files);$fileCounter++) { #For variant.function and exonic.function variants
		
		    print $FILEHANDLE "cat "; #Concatenate
		    print $FILEHANDLE $outFamilyDirectory."/".$familyID.$outfileEnding.$callType.$analysisType.".".$files[$fileCounter]." "; #Infile. 
		    print $FILEHANDLE $outFamilyDirectory."/".$familyID.$outfileEnding.$callType.$analysisType.".GRCh37_MT.".$files[$fileCounter]." "; # MT inFile
		    print $FILEHANDLE "> ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.$analysisType.".".$files[$fileCounter].".combined", "\n\n"; #Outfile
		    
		    print $FILEHANDLE "mv "; #replace original file with original information and MT info (if present)
		    print $FILEHANDLE $outFamilyDirectory."/".$familyID.$outfileEnding.$callType.$analysisType.".".$files[$fileCounter].".combined "; #Outfile
		    print $FILEHANDLE $outFamilyDirectory."/".$familyID.$outfileEnding.$callType.$analysisType.".".$files[$fileCounter], "\n\n"; #Original file	
		}
	    }
	}
	
	print $FILEHANDLE "rm ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.$analysisType."_temp", "\n\n"; #Remove temp file
    }
    close($FILEHANDLE);

    if ( ($scriptParameter{'pAnnovar'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {

	&FIDSubmitJob(0, $familyID, 1, $parameter{'pAnnovar'}{'chain'}, $fileName, 0);
    }
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
	$infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pGATKPhaseByTransmission'}{'fileEnding'};
    }
    else {
	$infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pGATKVariantRecalibration'}{'fileEnding'};
    }
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pGATKReadBackedPhasing'}{'fileEnding'};
    
###GATK ReadBackedPhasing
    
    print $FILEHANDLE "\n#GATK ReadBackedPhasing","\n\n";
    print $FILEHANDLE "java -Xmx4g ";

    if ($scriptParameter{'javaUseLargePages'} ne "no") {
	
	    print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
    }
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
	my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{'pGATKBaseRecalibration'}{'fileEnding'};
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
}

sub GATKPhaseByTransmission {
###GATK PhaseByTransmission computes the most likely genotype combination and phases trios and parent/child pairs given their genotype likelihoods and a mutation prior and phases all sites were parent/child transmission can be inferred unambiguously.
    
    my $familyID = $_[0]; #familyID NOTE: not sampleid
    my $aligner = $_[1];
    my $callType = $_[2]; #SNV,INDEL or BOTH
    
    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    &ProgramPreRequisites( $familyID, "GATKPhaseByTransmission", $aligner."/GATK", $callType, $FILEHANDLE, 1, 3);
    
    my $FamilyFileDirectory = $scriptParameter{'outDataDir'}."/".$familyID;
    my $inFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $outFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pGATKVariantRecalibration'}{'fileEnding'};
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pGATKPhaseByTransmission'}{'fileEnding'};
    
    unless (-f $FamilyFileDirectory."/".$familyID.".fam") { #Check to see if file already exists

	print $FILEHANDLE "#Generating '.fam' file for GATK PhaseByTransmission","\n\n";
	print $FILEHANDLE q?perl -nae 'print $F[0], "\t", $F[1], "\t", $F[2], "\t", $F[3], "\t", $F[4], "\t", $F[5], "\n";' ?.$scriptParameter{'pedigreeFile'}." > ".$FamilyFileDirectory."/".$familyID.".fam", "\n\n";
    }
    
###GATK PhaseByTransmission
    
    print $FILEHANDLE "\n#GATK PhaseByTransmission","\n\n";
    print $FILEHANDLE "java -Xmx4g ";
    
    if ($scriptParameter{'javaUseLargePages'} ne "no") {
	
	print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
    }
    print $FILEHANDLE "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
    print $FILEHANDLE "-l INFO "; #Set the minimum level of logging
    print $FILEHANDLE "-T PhaseByTransmission "; #Type of analysis to run
    print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file
    print $FILEHANDLE "-V: ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf "; #InFile (family vcf)
    &GATKPedigreeFlag($FILEHANDLE, $FamilyFileDirectory, "SILENT", "GATKPhaseByTransmission"); #Passing filehandle directly to sub routine using "*". Sub routine prints "--pedigree file" for family
    print $FILEHANDLE "-o ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.".vcf"; #OutFile
    
    close($FILEHANDLE);

    if ( ($scriptParameter{'pGATKPhaseByTransmission'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {

	&FIDSubmitJob(0, $familyID, 1, $parameter{'pGATKPhaseByTransmission'}{'chain'}, $fileName, 0);
    }
}

sub SnpEff {
###SnpEff annotates variants 
	
    my $familyID = $_[0]; #familyID NOTE: not sampleid
    my $aligner = $_[1];
    my $callType = $_[2]; #SNV,INDEL or BOTH
    
    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $nrCores = &NrofCoresPerSbatch(scalar(@snpSiftAnnotationFiles) + scalar(@contigs)); #Detect the number of cores to use from (@snpSiftAnnotationFiles and dbNSFP (=+1)
    &ProgramPreRequisites( $familyID, "SnpEff", $aligner."/GATK", $callType, $FILEHANDLE, $nrCores, 10);
    
    my $inFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $outFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pVCFParser'}{'fileEnding'};
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pSnpEff'}{'fileEnding'};
    my $analysisType = "";

    for (my $VEPOutputFilesCounter=0;$VEPOutputFilesCounter<$VEPOutputFiles;$VEPOutputFilesCounter++) {

	my $coreCounter = 1;

	if ($VEPOutputFilesCounter == 1) {
    
	    $analysisType = ".selected"; #SelectFile variants
	}
###SnpSift Annotation
	print $FILEHANDLE "\n#SnpSift Annotation","\n\n";	    

	for (my $snpSiftAnnotationFilesCounter=0;$snpSiftAnnotationFilesCounter<scalar(@snpSiftAnnotationFiles);$snpSiftAnnotationFilesCounter++) {

	    &PrintWait(\$snpSiftAnnotationFilesCounter, \$nrCores, \$coreCounter, $FILEHANDLE);
	    
	    print $FILEHANDLE "java -Xmx4g ";
	    
	    if ($scriptParameter{'javaUseLargePages'} ne "no") {
		
		print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	    }
	    print $FILEHANDLE "-jar ".$scriptParameter{'snpEffPath'}."/SnpSift.jar ";
	    print $FILEHANDLE "annotate ";
	    print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$snpSiftAnnotationFiles[$snpSiftAnnotationFilesCounter]." "; #Database
	    print $FILEHANDLE $inFamilyDirectory."/".$familyID.$infileEnding.$callType.$analysisType.".vcf "; #Infile
	    print $FILEHANDLE "| "; #Pipe
	    print $FILEHANDLE "perl ".$scriptParameter{'inScriptDir'}."/vcfParser.pl "; #Parses the vcf output
	    print $FILEHANDLE "> ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.$analysisType."_".$snpSiftAnnotationFiles[$snpSiftAnnotationFilesCounter]." &","\n\n"; #outfile	       
	}

###SnpSift dbNSFP
	if (scalar(@snpSiftDbNSFPAnnotations) > 0) { ##DbNSFP

	    my $combinedNrCores;

	    print $FILEHANDLE "\n#SnpSift dbNSFP","\n\n";

	    for (my $contigsCounter=0;$contigsCounter<scalar(@contigs);$contigsCounter++) {
		
		$combinedNrCores = scalar(@snpSiftAnnotationFiles) + $contigsCounter;
		&PrintWait(\$combinedNrCores, \$nrCores, \$coreCounter, $FILEHANDLE);
			    		
		print $FILEHANDLE "java -Xmx500m ";
		
		if ($scriptParameter{'javaUseLargePages'} ne "no") {
		    
		    print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
		}
		print $FILEHANDLE "-jar ".$scriptParameter{'snpEffPath'}."/SnpSift.jar ";
		print $FILEHANDLE "dbnsfp ";
		print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$scriptParameter{'snpSiftDbNSFPFile'}." "; #DbNSFP file
		print $FILEHANDLE "-f "; #fields to add
		
		for (my $snpSiftDbNSFPAnnotationsCounter=0;$snpSiftDbNSFPAnnotationsCounter<scalar(@snpSiftDbNSFPAnnotations);$snpSiftDbNSFPAnnotationsCounter++) {
		    
		    if ($snpSiftDbNSFPAnnotationsCounter == 0) { #First entry
			
			print $FILEHANDLE $snpSiftDbNSFPAnnotations[$snpSiftDbNSFPAnnotationsCounter]; #Database
		    }
		    elsif ($snpSiftDbNSFPAnnotationsCounter == (scalar(@snpSiftDbNSFPAnnotations)-1) ) { #Last entry
			
			print $FILEHANDLE ",".$snpSiftDbNSFPAnnotations[$snpSiftDbNSFPAnnotationsCounter]." "; #Database   
		    }
		    else {
			
			print $FILEHANDLE ",".$snpSiftDbNSFPAnnotations[$snpSiftDbNSFPAnnotationsCounter]; #Database
		    }
		}
		print $FILEHANDLE "<( "; #Pipe into SnpSift
		print $FILEHANDLE "cat "; #Read infile
		print $FILEHANDLE $inFamilyDirectory."/".$familyID.$infileEnding.$callType.$analysisType.".vcf "; #Infile
		print $FILEHANDLE "| ";
		print $FILEHANDLE "java -Xmx500m ";
		
		if ($scriptParameter{'javaUseLargePages'} ne "no") {
		    
		    print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
		}
		print $FILEHANDLE "-jar ".$scriptParameter{'snpEffPath'}."/SnpSift.jar ";
		print $FILEHANDLE "filter "; #Parallalize per contig for speed
		print $FILEHANDLE q?"( CHROM = '?.$contigs[$contigsCounter].q?' )"?;
		print $FILEHANDLE ") "; #End Pipe into SnpSift
		print $FILEHANDLE "> ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.$analysisType."_".$contigs[$contigsCounter]."_dbnsfp.vcf &","\n\n"; #outfile	       
	    }
	    print $FILEHANDLE "wait", "\n\n";
	    
	    ##Merge dbNSFP splitted vcfs 
	    print $FILEHANDLE "java -Xmx4g ";
	    
	    if ($scriptParameter{'javaUseLargePages'} ne "no") {
		
		print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	    }
	    print $FILEHANDLE "-jar ".$scriptParameter{'snpEffPath'}."/SnpSift.jar ";
	    print $FILEHANDLE "split -j "; #Joinf VCFs together

	    for (my $contigsCounter=0;$contigsCounter<scalar(@contigs);$contigsCounter++) {
	
		print $FILEHANDLE $outFamilyDirectory."/".$familyID.$outfileEnding.$callType.$analysisType."_".$contigs[$contigsCounter]."_dbnsfp.vcf "; #outfiles
	    }
	    print $FILEHANDLE "> ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.$analysisType."_dbnsfp.vcf ","\n\n"; #Merged outfile	  
	}

	##Merge vcf files to 1 
	print $FILEHANDLE "\n#GATK CombineVariants","\n\n";
	print $FILEHANDLE "java -Xmx4g ";
	
	if ($scriptParameter{'javaUseLargePages'} ne "no") {
	    
	    print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	}
	print $FILEHANDLE "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	print $FILEHANDLE "-l INFO "; #Set the minimum level of logging
	print $FILEHANDLE "-T CombineVariants "; #Type of analysis to run
	print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file
	for (my $snpSiftAnnotationFilesCounter=0;$snpSiftAnnotationFilesCounter<scalar(@snpSiftAnnotationFiles);$snpSiftAnnotationFilesCounter++) {
	    
	    print $FILEHANDLE "-V: ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.$analysisType."_".$snpSiftAnnotationFiles[$snpSiftAnnotationFilesCounter]." ";
	}
	if (scalar(@snpSiftDbNSFPAnnotations) > 0) { ##DbNSFP processed
	    
	    print $FILEHANDLE "-V: ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.$analysisType."_dbnsfp.vcf ";
	}
	print $FILEHANDLE "-o ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.$analysisType.".vcf "; #outfile
    }
    print $FILEHANDLE "\n\n";
    
##Remove Temp files
    print $FILEHANDLE "rm ";
    print $FILEHANDLE $outFamilyDirectory."/".$familyID.$outfileEnding.$callType.$analysisType."_*"."_dbnsfp.vcf","\n\n";

    close($FILEHANDLE);
    
    if ( ($scriptParameter{'pSnpEff'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	
	&FIDSubmitJob(0,$familyID, 1, $parameter{'pSnpEff'}{'chain'}, $fileName,0);
    }
}

sub VCFParser {
###VCFParser performs parsing of VariantEffectPredictor annotated variants 
	
    my $familyID = $_[0]; #familyID NOTE: not sampleid
    my $aligner = $_[1];
    my $callType = $_[2]; #SNV,INDEL or BOTH
    
    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    &ProgramPreRequisites( $familyID, "VCFParser", $aligner."/GATK", $callType, $FILEHANDLE, 1, 1);
    
    my $inFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $outFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pVariantEffectPredictor'}{'fileEnding'};
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pVCFParser'}{'fileEnding'};
    
###VCFParser

    print $FILEHANDLE "\n#VCFParser","\n\n";
    print $FILEHANDLE "perl ".$scriptParameter{'inScriptDir'}."/vcfParser.pl "; #Parses the VEP output to tab-sep format
    print $FILEHANDLE $outFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf "; #Infile
    
    if ($scriptParameter{'pVariantEffectPredictor'} > 0) {

	print $FILEHANDLE "--parseVEP ".$scriptParameter{'vcfParserVEPParser'}." "; #Parse VEP transcript specific entries
    }
    if ($scriptParameter{'vcfParserRangeFeatureFile'} ne "noUserInfo") {

	print $FILEHANDLE "-rf ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'vcfParserRangeFeatureFile'}." "; #List of genes to analyse separately	
	
	if (scalar(@vcfParserRangeFeatureAnnotationColumns) > 0) {

	    print $FILEHANDLE "-rf_ac "; #Range annotation columns
	    print $FILEHANDLE join(',', @vcfParserRangeFeatureAnnotationColumns)." ";	    
	}
    }
    if ($scriptParameter{'vcfParserSelectFile'} ne "noUserInfo") {

	print $FILEHANDLE "-sf ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'vcfParserSelectFile'}." "; #List of genes to analyse separately
	print $FILEHANDLE "-sf_mc ".$scriptParameter{'vcfParserSelectFileMatchingColumn'}." "; #Column of HGNC Symbol in SelectFile (-sf)

	if (scalar(@vcfParserSelectFeatureAnnotationColumns) > 0) {

	    print $FILEHANDLE "-sf_ac "; #Select annotation columns
	    print $FILEHANDLE join(',', @vcfParserSelectFeatureAnnotationColumns)." ";	    
	}
	print $FILEHANDLE "-sof ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.".selected.vcf ";
    }
    print $FILEHANDLE "> ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.".vcf"; #outfile
    print $FILEHANDLE "\n\n";

    close($FILEHANDLE);

    if ( ($scriptParameter{'pVCFParser'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {

	&FIDSubmitJob(0,$familyID, 1, $parameter{'pVCFParser'}{'chain'}, $fileName,0);
    }
}

sub VariantEffectPredictor {
###VariantEffectPredictor performs annotation of variants 
	
    my $familyID = $_[0]; #familyID NOTE: not sampleid
    my $aligner = $_[1];
    my $callType = $_[2]; #SNV,INDEL or BOTH
    
    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $nrCores = &NrofCoresPerSbatch(scalar(@contigs)); #Detect the number of cores to use
    &ProgramPreRequisites( $familyID, "VariantEffectPredictor", $aligner."/GATK", $callType, $FILEHANDLE, $scriptParameter{'maximumCores'}, 10);
    
    my $inFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $outFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pGATKVariantRecalibration'}{'fileEnding'};
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pVariantEffectPredictor'}{'fileEnding'};
    my $coreCounter = 1;

###VariantEffectPredictor
    
    print $FILEHANDLE "\n#VariantEffectPredictor","\n\n";

    for (my $contigsCounter=0;$contigsCounter<scalar(@contigs);$contigsCounter++) {

	&PrintWait(\$contigsCounter, \$nrCores, \$coreCounter, $FILEHANDLE);
	
	print $FILEHANDLE "perl ".$scriptParameter{'vepDirectoryPath'}."/variant_effect_predictor.pl "; #VEP script 
	print $FILEHANDLE "--dir_cache ".$scriptParameter{'vepDirectoryCache'}." "; #Specify the cache directory to use
	print $FILEHANDLE "--cache "; #Enables use of the cache.
	print $FILEHANDLE "--force_overwrite "; #force the overwrite of the existing file
	print $FILEHANDLE "--vcf "; #Writes output in VCF format.
	print $FILEHANDLE "--fork 16 "; #Enable forking, using the specified number of forks.
	print $FILEHANDLE "--buffer_size 20000 "; #Sets the internal buffer size, corresponding to the number of variations that are read in to memory simultaneously 
	print $FILEHANDLE "--refseq "; #Use the otherfeatures database to retrieve transcripts. This database contains transcript objects corresponding to RefSeq transcripts, along with CCDS and Ensembl ESTs
	print $FILEHANDLE "--chr ".$contigs[$contigsCounter]." ";

	for (my $vepFeatureCounter=0;$vepFeatureCounter<scalar(@vepFeatures);$vepFeatureCounter++) {
	    
	    print $FILEHANDLE "--".$vepFeatures[$vepFeatureCounter]." "; #Add VEP features to the output.
	    
	    if ( ($vepFeatures[$vepFeatureCounter] eq "sift") || ($vepFeatures[$vepFeatureCounter] eq "polyphen") )  { #Protein predictions
		
		print $FILEHANDLE "s "; #Add prediction score 
	    }
	}
	print $FILEHANDLE "-i ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf "; #InFile (family vcf)
	print $FILEHANDLE "-o ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType."_".$contigs[$contigsCounter].".vcf &", "\n\n"; #OutFile
    }
    print $FILEHANDLE "wait", "\n\n";

    ##Concatenate vcf files to 1 
    print $FILEHANDLE "\n#GATK CombineVariants","\n\n";
    print $FILEHANDLE "java -Xmx4g ";
    
    if ($scriptParameter{'javaUseLargePages'} ne "no") {
	
	print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
    }
    print $FILEHANDLE "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
    print $FILEHANDLE "-l INFO "; #Set the minimum level of logging
    print $FILEHANDLE "-T CombineVariants "; #Type of analysis to run
    print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file

    for (my $contigsCounter=0;$contigsCounter<scalar(@contigs);$contigsCounter++) {
	
	print $FILEHANDLE "-V: ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType."_".$contigs[$contigsCounter].".vcf "; #outfiles
    }
    print $FILEHANDLE "-o ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.".vcf", "\n\n"; #OutFile

    ##Remove Temp files
    print $FILEHANDLE "rm ";
    print $FILEHANDLE $outFamilyDirectory."/".$familyID.$outfileEnding.$callType."_*".".vcf","\n\n";
    print $FILEHANDLE $outFamilyDirectory."/".$familyID.$outfileEnding.$callType."_*".".vcf.*","\n\n"; #idx file
    
    close($FILEHANDLE);
	
    if ( ($scriptParameter{'pVariantEffectPredictor'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {

	&FIDSubmitJob(0,$familyID, 1, $parameter{'pVariantEffectPredictor'}{'chain'}, $fileName,0);
    }
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
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pGATKGenoTypeGVCFs'}{'fileEnding'};
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pGATKVariantRecalibration'}{'fileEnding'};
    
    unless (-e $scriptParameter{'outDataDir'}."/".$familyID."/".$familyID.".fam") { #Check to see if file already exists
	print $FILEHANDLE "#Generating '.fam' file for GATK VariantRecalibrator/ApplyRecalibration","\n\n";
	print $FILEHANDLE q?perl -nae 'print $F[0], "\t", $F[1], "\t", $F[2], "\t", $F[3], "\t", $F[4], "\t", $F[5], "\n";' ?.$scriptParameter{'pedigreeFile'}." > ".$outFamilyFileDirectory."/".$familyID.".fam", "\n\n";
    }  

    my $contigIntervalListFile = &GATKTargetListFlag($FILEHANDLE);

###GATK VariantRecalibrator
    
    my $variantRecalibratorOutFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK/intermediary";
    my @modes = ("SNP","INDEL");

    if ( ($scriptParameter{'analysisType'} eq "exomes") || ($scriptParameter{'analysisType'} eq "rapid") ) { #Exome/rapid analysis

	@modes = ("BOTH");
    }

    for (my $modeCounter=0;$modeCounter<scalar(@modes);$modeCounter++) { #SNP and INDEL will be recalibrated successively in the same file because when you specify eg SNP mode, the indels are emitted without modification, and vice-versa. Exome and Rapid will be processed using mode BOTH since there are to few INDELS to use in the recalibration model even though using 30 exome BAMS in Haplotypecaller step. 

	print $FILEHANDLE "\n\n#GATK VariantRecalibrator","\n\n";	
	print $FILEHANDLE "java -Xmx6g ";
	
	if ($scriptParameter{'javaUseLargePages'} ne "no") {
	    
	    print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	}
	print $FILEHANDLE "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	print $FILEHANDLE "-l INFO "; #Set the minimum level of logging
	print $FILEHANDLE "-T VariantRecalibrator "; #Type of analysis to run
	print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file
	print $FILEHANDLE "-recalFile ".$variantRecalibratorOutFamilyDirectory."/".$familyID.$infileEnding.$callType.".intervals ";
	print $FILEHANDLE "-rscriptFile ".$variantRecalibratorOutFamilyDirectory."/".$familyID.$infileEnding.$callType.".intervals.plots.R ";
	print $FILEHANDLE "-tranchesFile ".$variantRecalibratorOutFamilyDirectory."/".$familyID.$infileEnding.$callType.".intervals.tranches ";

	if ( ($scriptParameter{'analysisType'} eq "exomes") || ($scriptParameter{'analysisType'} eq "rapid") ) { #Exome/rapid analysis use combined reference for more power
	
	    print $FILEHANDLE "-L ".$contigIntervalListFile." ";#Target list file (merged or original)			
	    print $FILEHANDLE "-input ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf "; #Infile HaplotypeCaller combined vcf which used reference gVCFs to create combined vcf (30> samples gCVFs)
	}
	else { #WGS
	    
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
	&GATKPedigreeFlag($FILEHANDLE, $outFamilyFileDirectory, "SILENT", "GATKVariantRecalibration"); #Passing filehandle directly to sub routine using "*". Sub routine prints "--pedigree file" for family
	
###GATK ApplyRecalibration
	print $FILEHANDLE "\n\n#GATK ApplyRecalibration","\n\n";
	
	my $applyRecalibrationInFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK/intermediary";
	
	print $FILEHANDLE "java -Xmx6g ";

	if ($scriptParameter{'javaUseLargePages'} ne "no") {
	    
	    print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	}
	print $FILEHANDLE  "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	print $FILEHANDLE "-l INFO "; #Set the minimum level of logging
	print $FILEHANDLE "-T ApplyRecalibration ";
	print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file
	print $FILEHANDLE "-recalFile ".$applyRecalibrationInFamilyDirectory."/".$familyID.$infileEnding.$callType.".intervals ";
	print $FILEHANDLE "-tranchesFile ".$applyRecalibrationInFamilyDirectory."/".$familyID.$infileEnding.$callType.".intervals.tranches ";

	if ( ($scriptParameter{'analysisType'} eq "exomes") || ($scriptParameter{'analysisType'} eq "rapid")) { #Exome/rapid analysis use combined reference for more power
	    
	    print $FILEHANDLE "-L ".$contigIntervalListFile." ";#Target list file (merged or original)		
	    print $FILEHANDLE "-input ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf "; #Infile HaplotypeCaller combined vcf which used reference gVCFs to create combined vcf file
	    print $FILEHANDLE "-o ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType."_filtered.vcf ";
	}
	else  { #WGS
	    
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
	&GATKPedigreeFlag($FILEHANDLE, $outFamilyFileDirectory, "SILENT", "GATKVariantRecalibration"); #Passing filehandle directly to sub routine using "*". Sub routine prints "--pedigree file" for family    
	print $FILEHANDLE "--mode ".$modes[$modeCounter]." "; #Recalibration mode to employ (SNP|INDEL|BOTH)
    }
###GATK SelectVariants

##Removes all genotype information for exome ref and recalulates meta-data info for remaining samples in new file.
    if ( ($scriptParameter{'analysisType'} eq "exomes") || ($scriptParameter{'analysisType'} eq "rapid") ) { #Exome/rapid analysis
	
	print $FILEHANDLE "\n\n#GATK SelectVariants","\n\n";
	print $FILEHANDLE "java -Xmx2g ";

	if ($scriptParameter{'javaUseLargePages'} ne "no") {
	    
	    print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	}
	print $FILEHANDLE  "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	print $FILEHANDLE "-l INFO "; #Set the minimum level of logging
	print $FILEHANDLE "-T SelectVariants "; #Type of analysis to run
	print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file	
	print $FILEHANDLE "-L ".$contigIntervalListFile." ";#Target list file (merged or original)
	print $FILEHANDLE "-env "; #Don't include loci found to be non-variant after the subsetting procedure. 
	print $FILEHANDLE "-V: ".$inFamilyDirectory."/".$familyID.$outfileEnding.$callType."_filtered.vcf "; #InFile
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
	$sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'MostCompleteVCF'}{'Path'} = $outFamilyDirectory."/".$familyID.$outfileEnding.$callType.".vcf";	
	&FIDSubmitJob(0, $familyID, 1, $parameter{'pGATKVariantRecalibration'}{'chain'}, $fileName, 0);
    }
}

sub GATKGenoTypeGVCFs { 
#GATK GenoTypeGVCFs 
    
    my $familyID = $_[0]; #familyID NOTE: not sampleid
    my $aligner = $_[1];
    my $callType = $_[2]; #SNV,INDEL or BOTH

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle

    &ProgramPreRequisites( $familyID, "GATKGenoTypeGVCFs", $aligner."/GATK", $callType, $FILEHANDLE, $scriptParameter{'maximumCores'}, 10); #Activate when Haplotypecaller is multithreaded. 
    
    my $outFamilyFileDirectory = $scriptParameter{'outDataDir'}."/".$familyID;
    my $outFamilyDirectory = $scriptParameter{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'pGATKGenoTypeGVCFs'}{'fileEnding'};

    print $FILEHANDLE "mkdir -p ".$scriptParameter{'GATKTempDirectory'}.'$SLURM_JOB_ID', "\n\n"; 
    
    print $FILEHANDLE "#GATK GenoTypeGVCFs","\n\n";
    
    print $FILEHANDLE "java -Xmx4g ";

    if ($scriptParameter{'javaUseLargePages'} ne "no") {
	
	    print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
    }
    print $FILEHANDLE "-Djava.io.tmpdir=".$scriptParameter{'GATKTempDirectory'}.'$SLURM_JOB_ID'." "; #Temporary Directory
    print $FILEHANDLE "-jar ".$scriptParameter{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
    print $FILEHANDLE "-l INFO "; #Set the minimum level of logging
    print $FILEHANDLE "-T GenotypeGVCFs "; #Type of analysis to run
    print $FILEHANDLE "-R ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference file
    print $FILEHANDLE "-D ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'GATKHaploTypeCallerSNPKnownSet'}." "; #Known SNPs to use for annotation SNPs
    print $FILEHANDLE "-nt 16 "; #How many data threads should be allocated to running this analysis.

    if ($scriptParameter{'analysisType'} eq "exomes") {

	print $FILEHANDLE "-V ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'GATKGenoTypeGVCFsRefGVCFInfile'}." ";
    }

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) { #Collect infiles for all sampleIDs
	
	my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleIDs[$sampleIDCounter]."/".$aligner."/GATK";
	my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{'pGATKHaploTypeCaller'}{'fileEnding'};
	my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles($sampleIDs[$sampleIDCounter]);
	
	if ($PicardToolsMergeSwitch == 1) { #Alignment BAM-files merged previously
	    
	    print $FILEHANDLE "-V ".$inSampleDirectory."/".$infile.$infileEnding.".vcf "; #InFile
	}
	else { #No previous merge of alignment BAM-files

	    my $lanes = join("",@{$lane{ $sampleIDs[$sampleIDCounter] }}); #Extract lanes
	    print $FILEHANDLE "-V ".$inSampleDirectory."/".$sampleIDs[$sampleIDCounter]."_lanes_".$lanes.$infileEnding.".vcf "; #InFile(s)
	} 
    } 
    print $FILEHANDLE "-o ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.".vcf", "\n\n"; #OutFile

    print $FILEHANDLE "#Remove Temp Directory\n\n";
    print $FILEHANDLE "rm ";
    print $FILEHANDLE "-rf ".$scriptParameter{'GATKTempDirectory'}.'$SLURM_JOB_ID', "\n\n"; #Remove Temp Directory
    
    close($FILEHANDLE);  
    if ( ($scriptParameter{'pGATKGenoTypeGVCFs'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {

##Collect QC metadata info for later use
	$sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'MostCompleteVCF'}{'Path'} = $outFamilyDirectory."/".$familyID.$outfileEnding.$callType.".vcf";
	&FIDSubmitJob(0, $familyID, 1, $parameter{'pGATKGenoTypeGVCFs'}{'chain'}, $fileName, 0);
    }
}

sub GATKHaploTypeCaller { 
#GATK HaplotypeCaller 
    
    my $sampleID = $_[0];
    my $aligner = $_[1];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle

    &ProgramPreRequisites( $sampleID, "GATKHaploTypeCaller", $aligner."/GATK", 0, $FILEHANDLE, $scriptParameter{'maximumCores'}, 30); #Activate when Haplotypecaller is multithreaded. 
    
    my $outFamilyFileDirectory = $scriptParameter{'outDataDir'}."/".$scriptParameter{'familyID'};
    my $inSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/GATK";
    my $outSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/".$aligner."/GATK";
    my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pGATKBaseRecalibration'}{'fileEnding'};
    my $outfileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pGATKHaploTypeCaller'}{'fileEnding'};

    print $FILEHANDLE "mkdir -p ".$scriptParameter{'GATKTempDirectory'}.'$SLURM_JOB_ID', "\n\n"; 
    
    print $FILEHANDLE "#GATK HaplotypeCaller","\n\n";
    
    print $FILEHANDLE "java -Xmx4g ";

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
    print $FILEHANDLE "--emitRefConfidence GVCF "; #Mode for emitting experimental reference confidence scores. GVCF generates block summarized version of the BP_RESOLUTION data 
    print $FILEHANDLE "--variant_index_type LINEAR "; 
    print $FILEHANDLE "--variant_index_parameter 128000 ";
    print $FILEHANDLE "-pairHMM VECTOR_LOGLESS_CACHING "; #Hardware specific optmization

    if ( ($scriptParameter{'analysisType'} eq "exomes") || ($scriptParameter{'analysisType'} eq "rapid") ) { #Exome/rapid analysis

	my $contigIntervalListFile = &GATKTargetListFlag($FILEHANDLE);
	print $FILEHANDLE "-L ".$contigIntervalListFile." ";#Target list file (merged or original)
    }
    my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles($sampleID);
	
    if ($PicardToolsMergeSwitch == 1) { #Alignment BAM-files merged previously
	    
	print $FILEHANDLE "-I ".$inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile
	print $FILEHANDLE "-o ".$outSampleDirectory."/".$infile.$outfileEnding.".vcf", "\n\n"; #OutFile
    }
    else { #No previous merge of alignment BAM-files
	
	my $lanes = join("",@{$lane{$sampleID}}); #Extract lanes

	for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #For all infiles per lane
	    
	    my $infile = $infilesLaneNoEnding{$sampleID}[$infileCounter];
	    
	    print $FILEHANDLE "-I ".$inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile(s)
	} 
	print $FILEHANDLE "-o ".$outSampleDirectory."/".$sampleID."_lanes_".$lanes.$outfileEnding.".vcf", "\n\n"; #OutFile
    } 

    print $FILEHANDLE "#Remove Temp Directory\n\n";
    print $FILEHANDLE "rm ";
    print $FILEHANDLE "-rf ".$scriptParameter{'GATKTempDirectory'}.'$SLURM_JOB_ID', "\n\n"; #Remove Temp Directory
    
    close($FILEHANDLE);  
    if ( ($scriptParameter{'pGATKHaploTypeCaller'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {

	&FIDSubmitJob($sampleID, $scriptParameter{'familyID'}, 1, $parameter{'pGATKHaploTypeCaller'}{'chain'}, $fileName, 0);
    }
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

	    &PrintWait(\$infileCounter, \$nrCores, \$coreCounter, $FILEHANDLE);	    

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

	print $FILEHANDLE "java -Xmx4g ";
	if ($scriptParameter{'javaUseLargePages'} ne "no") {
	    
	    print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	}
	print $FILEHANDLE "-jar ".$scriptParameter{'picardToolsPath'}."/CollectMultipleMetrics.jar ";
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

	    &PrintWait(\$infileCounter, \$nrCores, \$coreCounter, $FILEHANDLE);	

	    my $infile = $infilesLaneNoEnding{$sampleID}[$infileCounter];	    
	    
	    print $FILEHANDLE "java -Xmx4g ";
	    if ($scriptParameter{'javaUseLargePages'} ne "no") {
		
		print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	    }
	    print $FILEHANDLE "-jar ".$scriptParameter{'picardToolsPath'}."/CollectMultipleMetrics.jar ";
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
	
	print $FILEHANDLE "java -Xmx4g ";
	if ($scriptParameter{'javaUseLargePages'} ne "no") {
	    
	    print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	}
	print $FILEHANDLE "-jar ".$scriptParameter{'picardToolsPath'}."/CalculateHsMetrics.jar ";
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

	    &PrintWait(\$infileCounter, \$nrCores, \$coreCounter, $FILEHANDLE);	    
	    
	    my $infile = $infilesLaneNoEnding{$sampleID}[$infileCounter];	    
	    
	    print $FILEHANDLE "java -Xmx4g ";
	    if ($scriptParameter{'javaUseLargePages'} ne "no") {
		
		print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	    }
	    print $FILEHANDLE "-jar ".$scriptParameter{'picardToolsPath'}."/CalculateHsMetrics.jar ";
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
	my $infileEnding = $sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'pChanjoCalculate'}{'fileEnding'};
	
	my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles($sampleID);	

	if ($PicardToolsMergeSwitch == 1) { #Files was merged previously
	   
	    print $FILEHANDLE $inSampleDirectory."/".$infile.$infileEnding.".json ";      	
	}
	else {
	    
	    for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #For all files from independent of merged or not

		&PrintWait(\$infileCounter, \$scriptParameter{'maximumCores'}, \$coreCounter, $FILEHANDLE);		
		
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
	    
	    &PrintWait(\$infileCounter, \$scriptParameter{'maximumCores'}, \$coreCounter, $FILEHANDLE);
	    
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
    else {

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
	if ($scriptParameter{'javaUseLargePages'} ne "no") {
	    
	    print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	}
	print $FILEHANDLE "-jar ".$scriptParameter{'picardToolsPath'}."/MarkDuplicates.jar ";
	print $FILEHANDLE "TMP_DIR=".$scriptParameter{'PicardToolsTempDirectory'}.'$SLURM_JOB_ID'." "; #Temp Directory
	print $FILEHANDLE "ASSUME_SORTED=true ";
	print $FILEHANDLE "CREATE_INDEX=TRUE "; #create a BAM index when writing a coordinate-sorted BAM file.
	print $FILEHANDLE "REMOVE_DUPLICATES=false ";
	print $FILEHANDLE "VALIDATION_STRINGENCY=STRICT ";
	print $FILEHANDLE "INPUT=".$inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile
	print $FILEHANDLE "OUTPUT=".$outSampleDirectory."/".$infile.$outfileEnding.".bam "; #OutFile
	print $FILEHANDLE "METRICS_FILE=".$outSampleDirectory."/".$infile.$outfileEnding."metric ", "\n\n"; #Metric file 
	
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
	    
	    &PrintWait(\$infileCounter, \$nrCores, \$coreCounter, $FILEHANDLE);
	    
	    my $infile = $infilesLaneNoEnding{$sampleID}[$infileCounter];
	    
	    print $FILEHANDLE "java -Xmx4g ";
	    if ($scriptParameter{'javaUseLargePages'} ne "no") {
		
		print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	    }
	    print $FILEHANDLE "-jar ".$scriptParameter{'picardToolsPath'}."/MarkDuplicates.jar ";
	    print $FILEHANDLE "TMP_DIR=".$scriptParameter{'PicardToolsTempDirectory'}.'$SLURM_JOB_ID'." "; #Temp Directory
	    print $FILEHANDLE "ASSUME_SORTED=true ";
	    print $FILEHANDLE "CREATE_INDEX=TRUE "; #create a BAM index when writing a coordinate-sorted BAM file.
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
    }
    close($FILEHANDLE);

    if ( ($scriptParameter{'pPicardToolsMarkduplicates'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {

	&FIDSubmitJob($sampleID, $scriptParameter{'familyID'}, 1, $parameter{'pPicardToolsMarkduplicates'}{'chain'}, $fileName, 0);
    }
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
		if ($scriptParameter{'javaUseLargePages'} ne "no") {
		    
		    print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
		}
		print $FILEHANDLE "-jar ".$scriptParameter{'picardToolsPath'}."/MergeSamFiles.jar ";
		print $FILEHANDLE "TMP_DIR=".$scriptParameter{'PicardToolsTempDirectory'}.'$SLURM_JOB_ID'." "; #Temp Directory
		print $FILEHANDLE "CREATE_INDEX=TRUE "; #create a BAM index when writing a coordinate-sorted BAM file.
		print $FILEHANDLE "OUTPUT=".$outSampleDirectory."/".$sampleID."_lanes_".$lanes.$outfileEnding.".bam "; #OutFile
	    }
	    
	    print $FILEHANDLE "INPUT=".$inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile
	}
	print $FILEHANDLE "\n\n";

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
		    if ($scriptParameter{'javaUseLargePages'} ne "no") {
			
			print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
		    }
		    print $FILEHANDLE "-jar ".$scriptParameter{'picardToolsPath'}."/MergeSamFiles.jar ";
		    print $FILEHANDLE "TMP_DIR=".$scriptParameter{'PicardToolsTempDirectory'}.'$SLURM_JOB_ID'." "; #Temp directory
		    print $FILEHANDLE "CREATE_INDEX=TRUE "; #create a BAM index when writing a coordinate-sorted BAM file.
		    print $FILEHANDLE "OUTPUT=".$outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes.$lanes.$outfileEnding.".bam "; #OutFile
		    print $FILEHANDLE "INPUT=".$inSampleDirectory."/".$sampleID."_lanes_".$lanes.$outfileEnding.".bam "; #InFile
		    print $FILEHANDLE "INPUT=".$picardToolsMergeSamFilesPrevious[$mergeFileCounter], "\n\n"; #$mergeLanes contains lane info on previous merge, $infilesLaneNoEnding{$sampleID}[0] uses @RG for very first .bam file to include read group for subsequent merges. Complete path. 

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
		if ($scriptParameter{'javaUseLargePages'} ne "no") {
		    
		    print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
		}
		print $FILEHANDLE "jar ".$scriptParameter{'picardToolsPath'}."/MergeSamFiles.jar ";
		print $FILEHANDLE "TMP_DIR=".$scriptParameter{'PicardToolsTempDirectory'}.'$SLURM_JOB_ID'." "; #Temp Directory
		print $FILEHANDLE "CREATE_INDEX=TRUE "; #create a BAM index when writing a coordinate-sorted BAM file.
		print $FILEHANDLE "OUTPUT=".$outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes.$lanes.$outfileEnding.".bam "; #OutFile
		print $FILEHANDLE "INPUT=".$inSampleDirectory."/".$infile.$infileEnding.".bam "; #InFile
		print $FILEHANDLE "INPUT=".$picardToolsMergeSamFilesPrevious[$mergeFileCounter],"\n\n"; #$mergeLanes contains lane info on previous merge, $infilesLaneNoEnding{$sampleID}[0] uses @RG for very first .bam file to include read group for subsequent merges. Complete path. 
		
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
	if ($scriptParameter{'javaUseLargePages'} ne "no") {
	    
	    print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	}
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
}

sub BWA_Aln {
#Generates BWA aln index on fastq files
    
    my $sampleID = $_[0];
    my $aligner = $_[1];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $time = ceil(2.5*scalar( @{ $infilesLaneNoEnding{$sampleID} })); #One full lane on Hiseq takes approx. 2,5 h for BWA_Aln to process, round up to nearest full hour.
    my $nrCores = 0;

    for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #For all files   

	&AdjustNrCoresToSeqMode(\$nrCores, \$sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesLaneNoEnding{ $sampleID }[$infileCounter]}{'sequenceRunType'});
    }

    $nrCores = &NrofCoresPerSbatch($nrCores ); #Make sure that the number of cores does not exceed maximum after incrementing above

    &ProgramPreRequisites($sampleID, "BwaAln", $aligner, 0, $FILEHANDLE, $nrCores, $time);

    my $inSampleDirectory =  $indirpath{$sampleID};
    my $outSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/bwa";
    my $coreCounter=1;    

    for (my $infileCounter=0;$infileCounter<scalar( @{ $infile{$sampleID} });$infileCounter++) {

	&PrintWait(\$infileCounter, \$nrCores, \$coreCounter, $FILEHANDLE);

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

	    &PrintWait(\$coreTracker, \$scriptParameter{'maximumCores'}, \$coreCounter, $FILEHANDLE);

	    for (my $readBatchProcessesCount=0;$readBatchProcessesCount<$nrReadBatchProcesses;$readBatchProcessesCount++) {
		
		if ($readBatchProcessesCount eq 0) {
		    
		    print $FILEHANDLE "java -Xmx4g ";
		    if ($scriptParameter{'javaUseLargePages'} ne "no") {
			
			print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
		    }
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
	    if ($scriptParameter{'javaUseLargePages'} ne "no") {
		
		print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	    }
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
	if ($scriptParameter{'javaUseLargePages'} ne "no") {
	    
	    print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	}
	print $FILEHANDLE "-jar ".$scriptParameter{'picardToolsPath'}."/SamFormatConverter.jar "; #Make sure that the BAM file BIN field is correct (Mosaik v.2.2.3 does according to Picard not set the bin field correctly)
	print $FILEHANDLE "VALIDATION_STRINGENCY=SILENT "; #Disable errors print 
	print $FILEHANDLE "INPUT=".$inSampleDirectory."/".$infile.".bam "; #InFile
	print $FILEHANDLE "OUTPUT=".$outSampleDirectory."/".$infile.".sam ", "\n\n"; #OutFile
	
	#SAM to BAM conversion 
	print $FILEHANDLE "java -Xmx4g ";
	if ($scriptParameter{'javaUseLargePages'} ne "no") {
	    
	    print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	}
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
	
	&PrintWait(\$infileCounter, \$nrCores, \$coreCounter, $FILEHANDLE);
	
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
}   

sub FastQC {
#Raw sequence quality analysis using FASTQC

    my $sampleID = $_[0];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $time = ceil(0.5*scalar( @{ $infile{$sampleID} })); #One full lane on Hiseq takes approx. 0.5 h for FASTQC to process, round up to nearest full hour.

    my $nrCores = 0;

    for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{$sampleID} });$infileCounter++) { #For all files   

	&AdjustNrCoresToSeqMode(\$nrCores, \$sampleInfo{ $scriptParameter{'familyID'} }{$sampleID}{'file'}{$infilesLaneNoEnding{ $sampleID }[$infileCounter]}{'sequenceRunType'});
    }

    $nrCores = &NrofCoresPerSbatch($nrCores ); #Make sure that the number of cores does not exceed maximum after incrementing above

    &ProgramPreRequisites($sampleID, "FastQC", "fastqc", 0, $FILEHANDLE , $nrCores, $time);
    
    my $inSampleDirectory = $indirpath{$sampleID};
    my $outSampleDirectory = $scriptParameter{'outDataDir'}."/".$sampleID."/fastqc";
    my $coreCounter=1;
    
    for (my $infileCounter=0;$infileCounter<scalar( @{ $infile{$sampleID} });$infileCounter++) {

	&PrintWait(\$infileCounter, \$nrCores, \$coreCounter, $FILEHANDLE);

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

    &PrintToFileHandles(\@printFilehandles, "\nNOTE: Will try to create required Annovar database files before executing ".$program."\n\n");

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

	    if ($annovarTableNames[$tableNamesCounter] =~/ensGene|refGene/) { #Special case for MT download
		
		print $FILEHANDLE "perl ".$scriptParameter{'annovarPath'}."/annotate_variation.pl "; #Annovar script 
		print $FILEHANDLE "-buildver GRCh37_MT "; #GenomeBuild version
		print $FILEHANDLE "-downdb ensGene "; #Db to download
		print $FILEHANDLE "-webfrom annovar "; #Download from annovar
		print $FILEHANDLE $annovarTemporaryDirectory."/ ", "\n\n"; #annovar/humandb directory is assumed
	    }
	    
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


sub BuildDownLoadablePreRequisites {
##Creates the downloadable resources

    my $familyID = $_[0];
    my $aligner = $_[1];
    my $program = $_[2];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    
    &ProgramPreRequisites($familyID, $program, $aligner, 0, $FILEHANDLE, 1, 4);

    print $FILEHANDLE "cd $scriptParameter{'referencesDir'}", "\n\n"; #Move to reference directory

    my $cosmidResourceDirectory = &CheckCosmidYAML();

    for my $parameterName (keys %supportedCosmidReferences) {

	if ( "p".$program eq $parameter{$parameterName}{'associatedProgram'}) {

	    if ($parameter{$parameterName}{'buildFile'} eq 1) {
	    
		&DownloadReference(\$program, $FILEHANDLE, $parameterName, \$cosmidResourceDirectory);
	    }
	}
    }
    
    close($FILEHANDLE);
    
    if ( ($scriptParameter{"p".$program} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	
	&FIDSubmitJob(0, $familyID, 6, "MIP", $fileName, 0);
    }
}

sub BuildPTCHSMetricPreRequisites {
##Creates the target "infiles_list" "padded.infile_list" and interval_list files

    my $familyID = $_[0];
    my $aligner = $_[1];
    my $program = $_[2];   
    my $FILEHANDLE = $_[3]; #Decides if a new sbatch script will be generated or handled by supplied FILEHANDLE

    my $parametersToEvaluate = 0; #The number of parameters to evaluate

    unless(defined($FILEHANDLE)) { #No supplied FILEHANDLE i.e. create new sbatch script    
	
	$FILEHANDLE = IO::Handle->new();#Create anonymous filehandle

	&ProgramPreRequisites($familyID, $program, $aligner, 0, $FILEHANDLE, 1, 1);
    }

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) { #All sampleIDs

	my $sampleIDBuildSwitchInfile = $parameter{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{'exomeTargetBedInfileLists'}{'buildFile'};
	my $sampleIDBuildFileInfile = $scriptParameter{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{'exomeTargetBedInfileLists'};
	my $infileListNoEnding = &RemoveFileEnding(\$sampleIDBuildFileInfile, $referenceFileEndings{'exomeTargetBedInfileLists'}); #For comparison of identical filename.bed files, to avoid creating files twice

	my $sampleIDBuildSwitchPadded = $parameter{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{'exomeTargetPaddedBedInfileLists'}{'buildFile'};
	my $sampleIDBuildFilePadded = $scriptParameter{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{'exomeTargetPaddedBedInfileLists'};
	my $paddedInfileListNoEnding = &RemoveFileEnding(\$sampleIDBuildFilePadded , $referenceFileEndings{'exomeTargetPaddedBedInfileLists'}); #For comparison of identical filename.bed files, to avoid creating files twice

	my $sampleIDBuildSwitchPaddedInterval = $parameter{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{'GATKTargetPaddedBedIntervalLists'}{'buildFile'};
	my $sampleIDBuildFilePaddedInterval = $scriptParameter{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{'GATKTargetPaddedBedIntervalLists'};
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
		
		&PrintToFileHandles(\@printFilehandles, "\nNOTE: Will try to create required ".$sampleIDBuildFile." file before executing ".$program."\n\n");
		
		print $FILEHANDLE "#SampleID:".$sampleIDs[$sampleIDCounter], "\n\n";
		print $FILEHANDLE "#CreateSequenceDictionary from reference", "\n";
		print $FILEHANDLE "java -Xmx2g ";
		if ($scriptParameter{'javaUseLargePages'} ne "no") {
		    
		    print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
		}
		print $FILEHANDLE "-jar ".$scriptParameter{'picardToolsPath'}."/CreateSequenceDictionary.jar ";
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
		
		print $FILEHANDLE "#Create".$referenceFileEndings{'exomeTargetBedInfileLists'}, "\n";
		print $FILEHANDLE "java -Xmx2g ";
		if ($scriptParameter{'javaUseLargePages'} ne "no") {
		    
		    print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
		}
		print $FILEHANDLE "-jar ".$scriptParameter{'picardToolsPath'}."/IntervalListTools.jar ";
		print $FILEHANDLE "INPUT=".$scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict_body_col_5 ";
		print $FILEHANDLE "OUTPUT=".$scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict_body_col_5_".$referenceFileEndings{'exomeTargetBedInfileLists'}." ", "\n\n";
		    
		my $intendedFilePathRef = \($scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEnding.$referenceFileEndings{'exomeTargetBedInfileLists'});
		my $temporaryFilePathRef = \($scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict_body_col_5_".$referenceFileEndings{'exomeTargetBedInfileLists'});    
		&PrintCheckExistandMoveFile($FILEHANDLE, $intendedFilePathRef, $temporaryFilePathRef);
	    }
	    if ( (defined($sampleIDBuildSwitchPadded) && ($sampleIDBuildSwitchPadded eq 1)) || (defined($sampleIDBuildSwitchPaddedInterval) && ($sampleIDBuildSwitchPaddedInterval eq 1)) ) {
		
		print $FILEHANDLE "#Create padded interval list", "\n";
		print $FILEHANDLE "java -Xmx2g ";
		if ($scriptParameter{'javaUseLargePages'} ne "no") {
		    
		    print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
		}
		print $FILEHANDLE "-jar ".$scriptParameter{'picardToolsPath'}."/IntervalListTools.jar ";
		print $FILEHANDLE "PADDING=100 "; #Add 100 nt on both sides of bed entry
		print $FILEHANDLE "INPUT=".$scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict_body_col_5 ";
		print $FILEHANDLE "OUTPUT=".$scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict_body_col_5".$referenceFileEndings{'exomeTargetPaddedBedInfileLists'}." ", "\n\n";
		
		my $intendedFilePathRef = \($scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEnding.$referenceFileEndings{'exomeTargetPaddedBedInfileLists'});
		my $temporaryFilePathRef = \($scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict_body_col_5".$referenceFileEndings{'exomeTargetPaddedBedInfileLists'});    
		&PrintCheckExistandMoveFile($FILEHANDLE, $intendedFilePathRef, $temporaryFilePathRef);
		
		if (defined($sampleIDBuildSwitchPaddedInterval) && ($sampleIDBuildSwitchPaddedInterval eq 1)) {
		    
		    ##Softlink '.interval_list' to padded .infile_list", "\n";
		    print $FILEHANDLE "ln -s "; #Softlink
		    print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEnding.$referenceFileEndings{'exomeTargetPaddedBedInfileLists'}." "; #Origin file
		    print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$sampleIDBuildFileNoEnding.$referenceFileEndings{'GATKTargetPaddedBedIntervalLists'}; #interval_list file
		}
		
		print $FILEHANDLE "\n\n";
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
		    &SetTargetFileGeneralBuildParameter(\$sampleIDBuildFilePaddedInterval, "GATKTargetPaddedBedIntervalLists", \$sampleIDBuildFile, \$sampleIDBuildFileNoEnding, \$sampleIDs[$sampleIDCounter]);
		}
		if ( (defined($sampleIDBuildSwitchInfile)) && ($sampleIDBuildSwitchInfile == 0) ) {
		    
		    $sampleIDBuildSwitchPadded = 0;
		    &SetTargetFileGeneralBuildParameter(\$sampleIDBuildFilePadded, "exomeTargetPaddedBedInfileLists", \$sampleIDBuildFile, \$sampleIDBuildFileNoEnding, \$sampleIDs[$sampleIDCounter]);
		}
		$sampleIDBuildSwitchInfile = 0;
		&SetTargetFileGeneralBuildParameter(\$sampleIDBuildFileInfile, "exomeTargetBedInfileLists", \$sampleIDBuildFile, \$sampleIDBuildFileNoEnding, \$sampleIDs[$sampleIDCounter]);
	    }
	}
    }
    unless($_[3]) { #Unless FILEHANDLE was supplied close it and submit 
    
	close($FILEHANDLE);
    
	if ( ($scriptParameter{"p".$program} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {
	
	    &FIDSubmitJob(0, $familyID, 6, "MIP", $fileName, 0); #"MIP" is required or the pPicardToolsCalulateHSMetrics jobs will start prematurely
	}
    }
}

sub BuildBwaPreRequisites {
##Creates the BwaPreRequisites using scriptParameters{'humanGenomeReference'} as reference.

    my $familyID = $_[0];
    my $aligner = $_[1];
    my $program = $_[2];
    
    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $randomInteger = int(rand(10000)); #Generate a random integer between 0-10,000.

    &ProgramPreRequisites($familyID, $program, $aligner, 0, $FILEHANDLE, 1, 3);

    &BuildHumanGenomePreRequisites($familyID, $aligner, $program, $FILEHANDLE, $randomInteger);

    if ($parameter{'bwaBuildReference'}{'buildFile'} eq 1) {

	&PrintToFileHandles(\@printFilehandles, "\nNOTE: Will try to create required ".$scriptParameter{'bwaBuildReference'}." index files before executing ".$program."\n\n");
	
	print $FILEHANDLE "#Building BWA index", "\n\n";
	print $FILEHANDLE "bwa index "; #index sequences in the FASTA format
	print $FILEHANDLE "-p ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'bwaBuildReference'}."_".$randomInteger." "; #prefix of the index
	print $FILEHANDLE "-a bwtsw "; #BWT construction algorithm
	print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'},"\n\n"; #the FASTA reference sequences file
	
	for (my $fileEndingsCounter=0;$fileEndingsCounter<scalar(@bwaBuildReferenceFileEndings);$fileEndingsCounter++) { #All fileEndings
	    
	    my $intendedFilePathRef = \($scriptParameter{'referencesDir'}."/".$scriptParameter{'bwaBuildReference'}.$bwaBuildReferenceFileEndings[$fileEndingsCounter]);
	    my $temporaryFilePathRef = \($scriptParameter{'referencesDir'}."/".$scriptParameter{'bwaBuildReference'}."_".$randomInteger.$bwaBuildReferenceFileEndings[$fileEndingsCounter]);    
	    &PrintCheckExistandMoveFile($FILEHANDLE, $intendedFilePathRef, $temporaryFilePathRef);
	}
	$parameter{'bwaBuildReference'}{'buildFile'} = 0; #Ensure that this subrutine is only executed once
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

    &BuildHumanGenomePreRequisites($familyID, $aligner, $program, $FILEHANDLE, $randomInteger);

    if ($parameter{'mosaikAlignReference'}{'buildFile'} eq 1) {

	&PrintToFileHandles(\@printFilehandles, "\nNOTE: Will try to create required ".$scriptParameter{'mosaikAlignReference'}." before executing ".$program."\n\n");

	print $FILEHANDLE "#Building MosaikAligner Reference", "\n\n";
	print $FILEHANDLE "MosaikBuild ";
	print $FILEHANDLE "-fr ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #the FASTA reference sequences file
	print $FILEHANDLE "-sn Homo_sapiens "; #Species name
	print $FILEHANDLE "-ga ".$humanGenomeReferenceSource.$humanGenomeReferenceVersion." "; #the genome assembly ID
	print $FILEHANDLE "-oa ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'mosaikAlignReference'}."_".$randomInteger, "\n\n";

	my $intendedFilePathRef = \($scriptParameter{'referencesDir'}."/".$scriptParameter{'mosaikAlignReference'});
	my $temporaryFilePathRef = \($scriptParameter{'referencesDir'}."/".$scriptParameter{'mosaikAlignReference'}."_".$randomInteger);    
	&PrintCheckExistandMoveFile($FILEHANDLE, $intendedFilePathRef, $temporaryFilePathRef);
    }
    if ($parameter{'mosaikJumpDbStub'}{'buildFile'} eq 1) {

	&PrintToFileHandles(\@printFilehandles, "\nNOTE: Will try to create required ".$scriptParameter{'mosaikJumpDbStub'}." before executing ".$program."\n\n");

	print $FILEHANDLE "#Building MosaikAligner JumpDatabase", "\n\n";
	
	print $FILEHANDLE "mkdir -p /scratch/mosaik_tmp", "\n";
	print $FILEHANDLE "export MOSAIK_TMP=/scratch/mosaik_tmp", "\n\n";
	
	print $FILEHANDLE "MosaikJump ";
	print $FILEHANDLE "-ia ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'mosaikAlignReference'}." "; #The input reference file  
	print $FILEHANDLE "-hs 15 "; #the hash size
	print $FILEHANDLE "-mem 24 "; #the amount memory used when sorting hashes
	print $FILEHANDLE "-out ".$scriptParameter{'referencesDir'}."/".$scriptParameter{'mosaikJumpDbStub'}."_".$randomInteger, "\n\n"; #Mosaik JumpDbStub for the output filenames

	for (my $fileEndingsCounter=0;$fileEndingsCounter<scalar(@mosaikJumpDbStubFileEndings);$fileEndingsCounter++) {

	    my $intendedFilePathRef = \($scriptParameter{'referencesDir'}."/".$scriptParameter{'mosaikJumpDbStub'}.$mosaikJumpDbStubFileEndings[$fileEndingsCounter]);
	    my $temporaryFilePathRef = \($scriptParameter{'referencesDir'}."/".$scriptParameter{'mosaikJumpDbStub'}."_".$randomInteger.$mosaikJumpDbStubFileEndings[$fileEndingsCounter]);    
	    &PrintCheckExistandMoveFile($FILEHANDLE, $intendedFilePathRef, $temporaryFilePathRef);
	}	
	
	print $FILEHANDLE "rm -rf /scratch/mosaik_tmp", "\n\n"; #Cleaning up temp directory
    }
    close($FILEHANDLE);
    
    if ( ($scriptParameter{"p".$program} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {

	&FIDSubmitJob(0, $familyID, 6, $parameter{"p".$program}{'chain'}, $fileName, 0);
    }
}

sub CheckBuildHumanGenomePreRequisites {
##Checks if the HumanGenomePreRequisites needs to be built 

    my $program = $_[0];

    for (my $fileEndingsCounter=0;$fileEndingsCounter<scalar(@humanGenomeReferenceFileEndings);$fileEndingsCounter++) {
	
	if ( ($parameter{"humanGenomeReference".$humanGenomeReferenceFileEndings[$fileEndingsCounter]}{'buildFile'} eq 1) || ($humanGenomeCompressed eq "compressed") ) {
	   
	    if ($scriptParameter{"p".$program} == 1) {
	
		&BuildHumanGenomePreRequisites($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, $program);
		last;#Will handle all meatfiles build within sbatch script
	    }
	}
    }
    ##Collect sequence contigs from human reference
    #&CollectSeqContigs(); #Reloads if required NOTE:Preparation for future changes but not activated yet
}

sub CheckBuildPTCHSMetricPreRequisites {

    my $program = $_[0];
    my $FILEHANDLE = $_[1]; #Decides if a new sbatch script will be generated or handled by supplied FILEHANDLE

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {
	
	if ($parameter{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{'exomeTargetBedInfileLists'}{'buildFile'} eq 1) {
	    
	    &BuildPTCHSMetricPreRequisites($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, $program, $FILEHANDLE);
	    last; #Will handle all build per sampleID within sbatch script
	}
	if ($parameter{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{'exomeTargetPaddedBedInfileLists'}{'buildFile'} eq 1) {
	    
	    &BuildPTCHSMetricPreRequisites($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, $program, $FILEHANDLE);
	    last; #Will handle all build per sampleID within sbatch script
	}
	if ( (defined($parameter{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{'GATKTargetPaddedBedIntervalLists'}{'buildFile'})) && ($parameter{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{'GATKTargetPaddedBedIntervalLists'}{'buildFile'} eq 1) ){
	    
	    &BuildPTCHSMetricPreRequisites($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, $program, $FILEHANDLE);
	    last; #Will handle all build per sampleID within sbatch script
	}
    }
}

sub DownloadReference {
##Downloads references using the database download manager Cosmid.     

    my $programRef = $_[0];
    my $FILEHANDLE = $_[1];
    my $parameterName = $_[2];
    my $cosmidResourceDirectoryRef = $_[3];

    if ($parameter{$parameterName}{'buildFile'} eq 1) { #Reference need to be built a.k.a downloaded
	
	&PrintToFileHandles(\@printFilehandles, "\nNOTE: Will try to download ".$scriptParameter{$parameterName}." before executing ".$$programRef."\n\n");

	print $FILEHANDLE "workon ".$scriptParameter{'pythonVirtualEnvironment'}, "\n\n"; #Activate python environment

	print $FILEHANDLE "cosmid "; #Database download manager
	print $FILEHANDLE "clone "; #clone resource
	print $FILEHANDLE $supportedCosmidReferences{$parameterName}{'cosmidName'}; #The actual reference
	unless ($supportedCosmidReferences{$parameterName}{'version'} eq "latest") { #Version to download

	    print $FILEHANDLE "#".$supportedCosmidReferences{$parameterName}{'version'},
	}
	print $FILEHANDLE "\n\n"; 

	print $FILEHANDLE "deactivate ", "\n\n"; #Deactivate python environment
	
	if ($supportedCosmidReferences{$parameterName}{'compressedSwitch'} eq "compressed") {

	    print $FILEHANDLE "gzip ";
	    print $FILEHANDLE "-d "; #Decompress
	    print $FILEHANDLE $$cosmidResourceDirectoryRef."/".$supportedCosmidReferences{$parameterName}{'cosmidName'}."/*.gz", "\n\n";
	}

	my $intendedFilePathRef = \($scriptParameter{'referencesDir'}."/".$scriptParameter{$parameterName});
	my $temporaryFilePathRef = \($$cosmidResourceDirectoryRef."/".$supportedCosmidReferences{$parameterName}{'cosmidName'}."/*");    
	&PrintCheckExistandMoveFile($FILEHANDLE, $intendedFilePathRef, $temporaryFilePathRef);
	
	#Remove temporary Cosmid resources directory
	print $FILEHANDLE "rm -rf ";
	print $FILEHANDLE $$cosmidResourceDirectoryRef."/".$supportedCosmidReferences{$parameterName}{'cosmidName'}."/;", "\n\n";
	#Remove temporary Cosmid ".cosmid.yaml" file
	print $FILEHANDLE "rm ";
	print $FILEHANDLE $$cosmidResourceDirectoryRef."/.cosmid.yaml", "\n\n";
	
	for my $supportedParameterName (keys %supportedCosmidReferences) {

	    if ($supportedCosmidReferences{$supportedParameterName}{'cosmidName'} eq $supportedCosmidReferences{$parameterName}{'cosmidName'}) { #Reset to 0 for all supportedCosmidReferences that are shared between modules 
		
		$parameter{$supportedParameterName}{'buildFile'} = 0; #Only need to download once per analysis call
	    }
	}
    }
}

sub BuildHumanGenomePreRequisites {
##Creates the humanGenomePreRequisites using scriptParameters{'humanGenomeReference'} as reference.

    my $familyID = $_[0];
    my $aligner = $_[1];
    my $program = $_[2];
    my $FILEHANDLE = $_[3]; #Decides if a new sbatch script will be generated or handled by supplied FILEHANDLE
    my $randomInteger = $_[4];

    unless(defined($FILEHANDLE)) { #No supplied FILEHANDLE i.e. create new sbatch script

	$FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
	$randomInteger = int(rand(10000)); #Generate a random integer between 0-10,000.

	&ProgramPreRequisites($familyID, $program, $aligner, 0, $FILEHANDLE, 1, 1);
    }

    print $FILEHANDLE "cd $scriptParameter{'referencesDir'}", "\n\n"; #Move to reference directory
    my $cosmidResourceDirectory = &CheckCosmidYAML();

    &DownloadReference(\$program, $FILEHANDLE, "humanGenomeReference", \$cosmidResourceDirectory);

    if ($humanGenomeCompressed eq "compressed") {

	&PrintToFileHandles(\@printFilehandles, "\nNOTE: Will try to decompres ".$scriptParameter{'humanGenomeReference'}." before executing ".$program."\n\n");

	print $FILEHANDLE "gzip ";
	print $FILEHANDLE "-d "; #Decompress
	print $FILEHANDLE $scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}, "\n\n";
	$scriptParameter{'humanGenomeReference'} =~ s/.fasta.gz/.fasta/g; #Replace the .fasta.gz ending with .fasta since this will execute before the analysis, hence changing the original file name ending from ".fastq" to ".fastq.gz".
	print STDOUT "Set humanGenomeReference to: ".$scriptParameter{'humanGenomeReference'}, "\n\n";
	$humanGenomeCompressed = "unCompressed";
    }

    &CheckBuildPTCHSMetricPreRequisites($program, $FILEHANDLE);

    for (my $fileEndingsCounter=0;$fileEndingsCounter<scalar(@humanGenomeReferenceFileEndings);$fileEndingsCounter++) { #All meta files    
	
	if ($parameter{"humanGenomeReference.dict"}{'buildFile'} eq 1) { #.dict file

	    &PrintToFileHandles(\@printFilehandles, "\nNOTE: Will try to create dict file for ".$scriptParameter{'humanGenomeReference'}." before executing ".$program."\n\n");
	    
	    print $FILEHANDLE "#CreateSequenceDictionary from reference", "\n";
	    print $FILEHANDLE "java -Xmx2g ";
	    if ($scriptParameter{'javaUseLargePages'} ne "no") {
		
		print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	    }
	    print $FILEHANDLE "-jar ".$scriptParameter{'picardToolsPath'}."/CreateSequenceDictionary.jar ";
	    print $FILEHANDLE "R=".$scriptParameter{'referencesDir'}."/".$scriptParameter{'humanGenomeReference'}." "; #Reference genome
	    print $FILEHANDLE "OUTPUT=".$scriptParameter{'referencesDir'}."/".$humanGenomeReferenceNameNoEnding."_".$randomInteger.".dict ", "\n\n"; #Output sequence dictionnary
	    
	    &PrintCheckExistandMoveFile($FILEHANDLE, \($scriptParameter{'referencesDir'}."/".$humanGenomeReferenceNameNoEnding.".dict"), \($scriptParameter{'referencesDir'}."/".$humanGenomeReferenceNameNoEnding."_".$randomInteger.".dict"));
	    
	    $parameter{"humanGenomeReference.dict"}{'buildFile'} = 0; #Only create once

	}
	if ($parameter{"humanGenomeReference.fasta.fai"}{'buildFile'} eq 1) {

	    &PrintToFileHandles(\@printFilehandles, "\nNOTE: Will try to create .fai file for ".$scriptParameter{'humanGenomeReference'}." before executing ".$program."\n\n");

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
    unless($_[3]) { #Unless FILEHANDLE was supplied close it and submit 
	
	close($FILEHANDLE);
    
	if ( ($scriptParameter{"p".$program} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {

	    &FIDSubmitJob(0, $familyID, 6, "MIP", $fileName, 0);
	}
    }
}

sub CheckCosmidInstallation {

    my $parameterNameRef = $_[0];
    
    if ($parameter{$$parameterNameRef}{'buildFile'} eq 1) {
	
	if (defined($scriptParameter{'pythonVirtualEnvironment'})) {
	
	    &PrintToFileHandles(\@printFilehandles, "Checking your Cosmid installation in preparation for download of ".$scriptParameter{$$parameterNameRef}."\n");
 
	    my $whichrReturn = `source ~/.bash_profile; workon $scriptParameter{'pythonVirtualEnvironment'};which cosmid;deactivate;`;
	    
	    if ($whichrReturn eq "") {
		print STDERR "\nMIP uses cosmid to download ".$scriptParameter{$$parameterNameRef}." and MIP could not find a cosmid installation in your python virtualenvironment".$scriptParameter{'pythonVirtualEnvironment'}." ","\n"; 
		exit;
	    }
	    else {

		&PrintToFileHandles(\@printFilehandles, "Found installation in ".$whichrReturn."\n");
	    }
	}
	else  {
	
	    print STDERR "\nCannot download".$scriptParameter{$$parameterNameRef}." without a '-pythonVirtualEnvironment'";
	    exit;
	}
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

    &DefinePlinkPedigree(); #Loads allowed entries and positons to be checked

    open(PEDF, "<".$fileName) or die "Can't open ".$fileName.":$!, \n";    
     
    while (<PEDF>) {
	
	chomp $_; #Remove newline
	
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
		    
		    my $foundElement =  &CheckEntryHashofArray(\%plinkPedigree, \$sampleElementsCounter, \$lineInfo[$sampleElementsCounter]);

		    if ($foundElement == 1) { #Invalid element found in file

			print STDERR "\nFound illegal element: '".$lineInfo[$sampleElementsCounter]."' in column '".$sampleElementsCounter."' in pedigree file: '".$fileName."' at line '".$.."'\n";
			print STDERR "\nPlease correct the entry before analysis.\n";
			print STDERR "\nMIP: Aborting run.\n\n";
			exit;
		    }
		    
		    my @elementInfo = split(";", $lineInfo[$sampleElementsCounter]); #Split element (if required)
		    
		    &CheckUniqueArrayElement(\@{ $sampleInfo{$familyID}{$sampleID}{$pedigreeFileElements[$sampleElementsCounter]} }, \@elementInfo); #Check if there are any new info and add it to sampleInfo if so. 
		   
		    if ($sampleInfo{$familyID}{$sampleID}{'Capture_kit'}) { #Add latest capture kit for each individual
			
			my $captureKit = $sampleInfo{$familyID}{$sampleID}{$pedigreeFileElements[$sampleElementsCounter]}[-1]; #Use only the last capture kit since it should be the most interesting
			
			if ($supportedCaptureKits{$captureKit}) {

			    if ($userExomeTargetBedInfileListsSwitch == 0) {
				    
				$scriptParameter{$familyID}{$sampleID}{'exomeTargetBedInfileLists'} = $supportedCaptureKits{$captureKit}.$referenceFileEndings{'exomeTargetBedInfileLists'}; #capture kit target in file_list
			    }
			    if ($userExomeTargetPaddedBedInfileListSwitch == 0) {
				
				$scriptParameter{$familyID}{$sampleID}{'exomeTargetPaddedBedInfileLists'} = $supportedCaptureKits{$captureKit}.$referenceFileEndings{'exomeTargetPaddedBedInfileLists'}; #capture kit padded target infile_list                               
			    }
			    if ($userExomeTargetPaddedBedIntervalListSwitch == 0) {
				
				$scriptParameter{$familyID}{$sampleID}{'GATKTargetPaddedBedIntervalLists'} = $supportedCaptureKits{$captureKit}.$referenceFileEndings{'GATKTargetPaddedBedIntervalLists'}; #capture kit padded target interval_list                          
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
}

sub DefinePlinkPedigree {
##Defines which entries are allowed and links them to position

    $plinkPedigree{4} = [1, 2, "other"]; #Sex
    $plinkPedigree{5} = [-9, 0, 1, 2]; #Phenotype

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
		    for (my $infileCounter=0;$infileCounter<scalar( @{ $infilesLaneNoEnding{ $sampleIDs[$sampleIDCounter] } });$infileCounter++) {
			
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
    &PrintToFileHandles(\@printFilehandles, "Sbatch script submitted, job id: $jobID\n");
    &PrintToFileHandles(\@printFilehandles, "To check status of job, please run \'squeue -j $jobID\'\n");
    &PrintToFileHandles(\@printFilehandles, "To cancel job, please run \'scancel $jobID\'\n");
}

sub NrofCoresPerSbatch {
##Set the number of cores to allocate per sbatch job
    
    my $nrCores = $_[0];
    
    if ($nrCores > $scriptParameter{'maximumCores'}) { #Set number of cores depending on how many lanes to process
	
	$nrCores = $scriptParameter{'maximumCores'}; #Set to max on cluster
    }
    return $nrCores;
}

sub CollectInfiles {
##Collects the ".fastq(.gz)" files from the supplied infiles directory. Checks if any files exist
    
    for (my $inputDirectoryCounter=0;$inputDirectoryCounter<scalar(@inFilesDirs);$inputDirectoryCounter++) { #Collects inputfiles
	
	my @infiles = `cd $inFilesDirs[ $inputDirectoryCounter ];ls *.fastq*;`; #cd to input dir and collect fastq files and fastq.gz files
	
	if (scalar(@infiles) == 0) { #No "*.fastq*" infiles
	    
	    print STDERR "Could not find any '.fastq' files in supplied infiles directory ".$inFilesDirs[ $inputDirectoryCounter ], "\n";
	    exit;
	}
	&PrintToFileHandles(\@printFilehandles, "\nReads from Platform\n");
	&PrintToFileHandles(\@printFilehandles, "\nSample ID\t".$sampleIDs[$inputDirectoryCounter]."\n");
	print STDOUT "Inputfiles:\n",@ { $infile{ $sampleIDs[$inputDirectoryCounter] }  =[@infiles] }, "\n"; #hash with sample id as key and inputfiles in dir as array 
	print MIPLOG "Inputfiles:\n",@ { $infile{ $sampleIDs[$inputDirectoryCounter] }  =[@infiles] }, "\n";
	
	$indirpath{$sampleIDs[$inputDirectoryCounter]} = $inFilesDirs[ $inputDirectoryCounter ];  #Catch inputdir path
	chomp(@infiles);    #Remove newline from every entry in array
	$infile{ $sampleIDs[$inputDirectoryCounter] }  =[@infiles]; #Reload files into hash (kept above newline just for print STDOUT)
    }
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
	    else { #No regexp match i.e. file does not follow filename convention 

		print STDERR "Could not detect MIP file name convention for file: ".$infile{$sampleID}[$infileCounter].". \n\nPlease check that the file name follows the specified convention.", "\n";
		exit;
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
			if ( ($parameterName eq "GATKTargetPaddedBedIntervalLists") && ($scriptParameter{'analysisType'} ne "genomes") ) { #GATKTargetPaddedBedIntervalLists is a comma separated list 
			    
			    &SetTargetandAutoBuild(\@sampleIDs, \$parameterName, \$referenceFileEndings{'GATKTargetPaddedBedIntervalLists'});
			}  
			if ($parameterName eq "vepFeatures") { #Input from config file
    
			    @vepFeatures = split(/,/, $scriptParameter{'vepFeatures'}); #Transfer to array
			}
			if ($parameterName eq "vcfParserRangeFeatureAnnotationColumns") { #Input from config file
    
			    @vcfParserRangeFeatureAnnotationColumns = split(/,/, $scriptParameter{'vcfParserRangeFeatureAnnotationColumns'}); #Transfer to array
			}
			if ($parameterName eq "vcfParserSelectFeatureAnnotationColumns") { #Input from config file
    
			    @vcfParserSelectFeatureAnnotationColumns = split(/,/, $scriptParameter{'vcfParserSelectFeatureAnnotationColumns'}); #Transfer to array
			}
			if ($parameterName eq "snpSiftAnnotationFiles") { #Input from config file
    
			    @snpSiftAnnotationFiles = split(/,/, $scriptParameter{'snpSiftAnnotationFiles'}); #Transfer to array
			}
			if ($parameterName eq "snpSiftDbNSFPAnnotations") { #Input from config file
    
			    @snpSiftDbNSFPAnnotations = split(/,/, $scriptParameter{'snpSiftDbNSFPAnnotations'}); #Transfer to array
			}
			if ($parameterName eq "annovarTableNames") { #Input from config file
			    
			    @annovarTableNames = split(/,/, $scriptParameter{'annovarTableNames'}); #Transfer to array
			}
			if ($parameterName eq "humanGenomeReference") {
			    
			    ($humanGenomeReferenceVersion, $humanGenomeReferenceSource, $humanGenomeReferenceNameNoEnding) = &ParseHumanGenomeReference($scriptParameter{'humanGenomeReference'});
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
			elsif ( ($parameterName eq "GATKTargetPaddedBedIntervalLists") && ($scriptParameter{'analysisType'} ne "genomes") ) { #Note that potential pedigree files entries will be updated with GenomeReferenceSource and version here 
			    
			    &SetTargetandAutoBuild(\@sampleIDs, \$parameterName, \$referenceFileEndings{'GATKTargetPaddedBedIntervalLists'});
			}
			elsif ($parameterName eq "vepFeatures") {
			    
			    @vepFeatures = ("refseq", "hgvs", "symbol", "numbers", "sift", "polyphen", "humdiv"); #Set default vep features
			    &EnableArrayParameter(\@vepFeatures, \$parameterName);
			}
			elsif ($parameterName eq "snpSiftAnnotationFiles") { #Input from config file
			    
			    @snpSiftAnnotationFiles = ("dbsnp_138.b37.excluding_sites_after_129.vcf", "dbsnp_138.b37.vcf", "1000G_phase1.indels.b37.vcf", "1000G_phase1.snps.high_confidence.b37.vcf", "ESP6500SI-V2-SSA137.updatedProteinHgvs.snps_indels.vcf"); #Set default snpSiftAnnotationFiles
			    &EnableArrayParameter(\@snpSiftAnnotationFiles, \$parameterName);
			}
			elsif ($parameterName eq "snpSiftDbNSFPAnnotations") { #Input from config file
    
			    @snpSiftDbNSFPAnnotations = ("SIFT_pred", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_pred", "LRT_pred", "MutationTaster_pred", "GERP++_NR", "GERP++_RS", "phastCons100way_vertebrate", "1000Gp1_AF", "ESP6500_AA_AF"); #Set default snpSiftDbNSFPAnnotations
			    &EnableArrayParameter(\@snpSiftDbNSFPAnnotations, \$parameterName);
			}
			elsif ($parameterName eq "annovarTableNames") {
			    
			    @annovarTableNames = ("refGene", "mce46way", "gerp++elem", "segdup", "tfbs", "mirna", "snp137NonFlagged", "1000g2012apr_all", "esp6500si_all", "ljb2_sift", "ljb2_pp2hdiv", "ljb2_pp2hvar", "ljb2_mt", "ljb2_lrt", "ljb2_gerp++","ljb2_phylop"); #Set default annovar table names
			    &EnableArrayParameter(\@annovarTableNames, \$parameterName);
			}
			else {
			    
			    if ($parameterName eq "humanGenomeReference") {
				
				($humanGenomeReferenceVersion, $humanGenomeReferenceSource, $humanGenomeReferenceNameNoEnding) = &ParseHumanGenomeReference($parameterDefault);
			    }			    
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
			elsif ( ($parameterName eq "GATKExomeReferenceSNPs") && ($scriptParameter{'analysisType'} ne "rapid")) { #Do nothing since file is not required unless rapid mode is enabled
			}
			elsif ( ($parameterName eq "bwaMemRapidDb") && ($scriptParameter{'analysisType'} ne "rapid")) { #Do nothing since file is not required unless rapid mode is enabled
			}
			elsif ( ($parameterName eq "GATKGenoTypeGVCFsRefGVCFInfile") && ($scriptParameter{'analysisType'} =~/genomes/) ) { #Do nothing since file is not required unless exome or rapid mode is enabled
			}
			elsif ( ($parameterName eq "vcfParserRangeFeatureAnnotationColumns") && ( $scriptParameter{'vcfParserRangeFeatureFile'} eq "noUserInfo") ) { #Do nothing since no SelectFile was given
			} 
			elsif ( ($parameterName eq "vcfParserSelectFeatureAnnotationColumns") && ( $scriptParameter{'vcfParserRangeFeatureFile'} eq "noUserInfo") ) { #Do nothing since no SelectFile was given
			}
			elsif ( ($parameterName eq "vcfParserSelectFileMatchingColumn") && ( $scriptParameter{'vcfParserSelectFile'} eq "noUserInfo") ) { #Do nothing since no SelectFile was given
			
			}
			elsif ( ($parameterName eq "mergeAnnotatedVariantsTemplateSelectFile") && ( $scriptParameter{'vcfParserSelectFile'} eq "noUserInfo") ) { #Do nothing since no SelectFile was given
			}
			elsif ( ($parameterName eq "geneFile") && ($scriptParameter{'pVariantEffectPredictor'} > 0) ) { #Do nothing since VEP annotations can be used			    
			
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
		    elsif ( ($parameterName eq "GATKTargetPaddedBedIntervalLists") && ($scriptParameter{'analysisType'} ne "genomes") ) {	    
			
			&EnableArrayParameter(\@GATKTargetPaddedBedIntervalLists, \$parameterName);
			&CompareArrayElements(\@sampleIDs, \@GATKTargetPaddedBedIntervalLists, "sampleIDs", $parameterName);
			&SetAutoBuildAndScriptParameterPerSample(\@sampleIDs, \@GATKTargetPaddedBedIntervalLists, \$parameterName);
		    }
		    elsif ($parameterName eq "vcfParserRangeFeatureAnnotationColumns") {
			
			&EnableArrayParameter(\@vcfParserRangeFeatureAnnotationColumns, \$parameterName);
		    }
		    elsif ($parameterName eq "vcfParserSelectFeatureAnnotationColumns") {
			
			&EnableArrayParameter(\@vcfParserSelectFeatureAnnotationColumns, \$parameterName);
		    }
		    elsif ($parameterName eq "vepFeatures") {
			
			&EnableArrayParameter(\@vepFeatures, \$parameterName);
		    }
		    elsif ($parameterName eq "snpSiftAnnotationFiles") {
			
			&EnableArrayParameter(\@snpSiftAnnotationFiles, \$parameterName);
		    }
		    elsif ($parameterName eq "snpSiftDbNSFPAnnotations") {
			
			&EnableArrayParameter(\@snpSiftDbNSFPAnnotations, \$parameterName);
		    }
		    elsif ($parameterName eq "annovarTableNames") {
			
			&EnableArrayParameter(\@annovarTableNames, \$parameterName);
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
				
				$sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'program'}{"GATK"}{'Version'} = $1;
			    }
			}
			if ($parameterName eq "picardToolsPath") { #To enable addition of version to sampleInfo                                                                       

                            if ($scriptParameter{$parameterName}=~/picard-tools-([^,]+)/) {
                                
				$sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'program'}{"PicardTools"}{'Version'} = $1;
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
			$sampleInfo{$scriptParameter{'familyID'}}{$scriptParameter{'familyID'}}{"HumanGenomeBuild"}{'Path'} = $scriptParameter{'referencesDir'}."/".$scriptParameter{$parameterName};
			$sampleInfo{$scriptParameter{'familyID'}}{$scriptParameter{'familyID'}}{"HumanGenomeBuild"}{'Source'} = $humanGenomeReferenceSource;
			$sampleInfo{$scriptParameter{'familyID'}}{$scriptParameter{'familyID'}}{"HumanGenomeBuild"}{'Version'} = $humanGenomeReferenceVersion;

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
			
			if ( ($parameterName eq "GATKTargetPaddedBedIntervalLists") && ($scriptParameter{'analysisType'} eq "genomes") ) { #No need to check since genomes does not use GATKTargetPaddedBedIntervalLists
			}
			else {
			    
			    for (my $sampleIDsCounter=0;$sampleIDsCounter<scalar(@sampleIDs);$sampleIDsCounter++) { #All sampleIDs
				
				 &CheckSupportedFileEnding(\($scriptParameter{'referencesDir'}."/".$scriptParameter{ $scriptParameter{'familyID'} }{$sampleIDs[$sampleIDsCounter]}{$parameterName}), \$referenceFileEndings{$parameterName}, \$parameterName);
				 &CheckExistance(\($scriptParameter{'referencesDir'}."/".$scriptParameter{ $scriptParameter{'familyID'} }{$sampleIDs[$sampleIDsCounter]}{$parameterName}), \$parameterName, "f", \$sampleIDs[$sampleIDsCounter]);
			    
				my $exomeTargetBedFileNoEnding = &RemoveFileEnding(\$scriptParameter{ $scriptParameter{'familyID'} }{$sampleIDs[$sampleIDsCounter]}{$parameterName} , $referenceFileEndings{$parameterName}); #Remove ".fileending" from reference filename
				 &CheckTargetExistFileBed(\$exomeTargetBedFileNoEnding, $parameterName);
			    }
			    undef($scriptParameter{$parameterName}); #Remove parameter to avoid unnecessary print to STDOUT and config
			}
		    }
		    elsif ($parameterName eq "annovarTableNames") {
			
			&DefineAnnovarTables(); #Set all AnnovarTables properties
			my $intendedFilePathRef;
			
			for (my $tableNamesCounter=0;$tableNamesCounter<scalar(@annovarTableNames);$tableNamesCounter++) { #All AnnovarTables
		 
			    if (defined($annovarTables{$annovarTableNames[$tableNamesCounter]})) { #Supported Annovar database

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
			    else { #Annovar Table not supported by MIP
				
				print STDOUT "\nNOTE: You supplied Annovar database: ".$annovarTableNames[$tableNamesCounter]." which is not supported by MIP. MIP can only process supported annovar databases\n";
				&PrintSupportedAnnovarTableNames();
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
				
				$sampleInfo{$scriptParameter{'familyID'}}{$scriptParameter{'familyID'}}{'pedigreeFile'}{'Path'} = $scriptParameter{'pedigreeFile'}; #Add pedigreeFile to sampleInfo
				$sampleInfo{$scriptParameter{'familyID'}}{$scriptParameter{'familyID'}}{'pedigreeFileAnalysis'}{'Path'} = $scriptParameter{'outDataDir'}."/".$scriptParameter{'familyID'}."/qc_pedigree.yaml"; #Add pedigreeFile info used in this analysis to SampleInfoFile
			    }
			} 
		    }
		    elsif ( ($parameterName eq "GATKExomeReferenceSNPs") && ($scriptParameter{'analysisType'} ne "rapid")) { #Do nothing since file is not required unless rapid mode is enabled
		    }
		    elsif ( ($parameterName eq "bwaMemRapidDb") && ($scriptParameter{'analysisType'} ne "rapid")) { #Do nothing since file is not required unless rapid mode is enabled
		    }
		    elsif ( ($parameterName eq "GATKGenoTypeGVCFsRefGVCFInfile") && ($scriptParameter{'analysisType'} =~/genomes/) ) { #Do nothing since file is not required unless exome mode is enabled
		    }
                    elsif ( ($parameterName eq "vcfParserRangeFeatureFile") && ( $scriptParameter{'vcfParserRangeFeatureFile'} eq "noUserInfo") ) { #Do nothing since no RangeFile was given
		    }
		    elsif ($parameterName eq "vcfParserSelectFile") {
			
			if ($scriptParameter{'vcfParserSelectFile'} eq "noUserInfo") { #Do nothing since no SelectFile was given
			    
			}
			else { #To enable addition of selectFile to sampleInfo                                                                       
			    
			    if ($scriptParameter{$parameterName}=~/v(\d+\.\d+)/) {
				
				$sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'database'}{"SelectFile"}{'Version'} = $1;
			    }
			    $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'database'}{"SelectFile"}{'File'} = $scriptParameter{$parameterName};
			}
		    }
		    elsif ( ($parameterName eq "mergeAnnotatedVariantsTemplateSelectFile") && ( $scriptParameter{'vcfParserSelectFile'} eq "noUserInfo") ) { #Do nothing since no SelectFile was given
		    }
                    elsif ( ($parameterName eq "geneFile") && ($scriptParameter{'pVariantEffectPredictor'} > 0) ) { #Do nothing since VEP annotations can be used			    
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
		if ($parameterName eq "instanceTag") {

		    $sampleInfo{$scriptParameter{'familyID'}}{$scriptParameter{'familyID'}}{$parameterName} = $scriptParameter{$parameterName};
		}
		if ($parameterName eq "researchEthicalApproval") {

		    $sampleInfo{$scriptParameter{'familyID'}}{$scriptParameter{'familyID'}}{$parameterName} = $scriptParameter{$parameterName};
		}
	    }
	    
	    if ( $parameterType eq "program") {

		if($parameterValue eq "nocmdinput") { #No input from cmd
		    
		    if (defined($scriptParameter{$parameterName})) { #Input from config file - do nothing
			
		    }
		    elsif ($parameterDefault ne "nodefault") {
			    
			$scriptParameter{$parameterName} = $parameterDefault; #Set default value
		    }
		}
		else {

		    $scriptParameter{$parameterName} = $parameterValue;
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
		    
		    if ( ($scriptParameter{'pBwaAln'} == 0) && ($scriptParameter{'pBwaSampe'} == 0) && ($scriptParameter{'pBwaMem'} == 0) ) {
			
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
		elsif ( ($scriptParameter{'pBwaAln'} > 0) || ($scriptParameter{'pBwaSampe'} > 0) || ($scriptParameter{'pBwaMem'} > 0)) { #BWA track
		    
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
					
					$sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{'pPicardToolsMergeSamFiles'}{'fileEnding'} = $tempFileEnding{ $sampleIDs[$sampleIDCounter] }.$parameter{$orderParameterElement}{'fileEnding'}; #Adds from previous entry 
				    }
				    else {
					
					$sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{'pPicardToolsMergeSamFiles'}{'fileEnding'} = $tempFileEnding{ $sampleIDs[$sampleIDCounter] }."";
				    }
				}
				else {
				    if (defined($tempFileEnding{ $sampleIDs[$sampleIDCounter] })) {
					
					$sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{$orderParameterElement}{'fileEnding'} = $tempFileEnding{ $sampleIDs[$sampleIDCounter] }.$parameter{$orderParameterElement}{'fileEnding'};
				    }
				    else  { #First module that should add filending
					$sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{$orderParameterElement}{'fileEnding'} = $parameter{$orderParameterElement}{'fileEnding'};
				    } 
				}
			    }
			    else { #Do not add new module fileEnding
				$sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{$orderParameterElement}{'fileEnding'} = $tempFileEnding{ $sampleIDs[$sampleIDCounter] };
			    }
			    $tempFileEnding{ $sampleIDs[$sampleIDCounter] } = $sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{$orderParameterElement}{'fileEnding'}; #To enable sequential build-up of fileending
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
				    $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{$orderParameterElement}{'fileEnding'} = $tempFileEnding{$scriptParameter{'familyID'}}.$parameter{$orderParameterElement}{'fileEnding'};
				}
				else  { #First module that should add filending
				    $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{$orderParameterElement}{'fileEnding'} = $parameter{$orderParameterElement}{'fileEnding'};
				}
				$tempFileEnding{$scriptParameter{'familyID'}} = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{$orderParameterElement}{'fileEnding'}; #To enable sequential build-up of fileending 
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
			
				unless (defined($tempFileEnding{$chainfork}{ $sampleIDs[$sampleIDCounter] })) {	
				    $tempFileEnding{$chainfork}{ $sampleIDs[$sampleIDCounter] } = $tempFileEnding{ $sampleIDs[$sampleIDCounter] }; #Inherit current MAIN chain. 
				}
				if (defined($tempFileEnding{$chainfork}{ $sampleIDs[$sampleIDCounter] })) {
				    $sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{$orderParameterElement}{'fileEnding'} = $tempFileEnding{$chainfork}{ $sampleIDs[$sampleIDCounter] }.$parameter{$orderParameterElement}{'fileEnding'};
				}
				else  { #First module that should add filending
				    $sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{$orderParameterElement}{'fileEnding'} = $parameter{$orderParameterElement}{'fileEnding'};
				} 
			    }
			    else { #Do not add new module fileEnding
				$sampleInfo{ $scriptParameter{'familyID'} }{ $sampleIDs[$sampleIDCounter] }{$orderParameterElement}{'fileEnding'} = $tempFileEnding{$chainfork}{ $sampleIDs[$sampleIDCounter] };
			    }
			    #NOTE: No sequential build-up of fileending
			}
###Other/Per familyID

			if ($scriptParameter{$orderParameterElement} > 0) { #Fileending should be added
			    
			    unless (defined($tempFileEnding{$chainfork}{$scriptParameter{'familyID'}})) {	
				$tempFileEnding{$chainfork}{$scriptParameter{'familyID'}} =  $tempFileEnding{$scriptParameter{'familyID'}}; #Inherit current MAIN chain. 
			    }
			    if (defined($tempFileEnding{$chainfork}{$scriptParameter{'familyID'}})) {
				$sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{$orderParameterElement}{'fileEnding'} = $tempFileEnding{$chainfork}{$scriptParameter{'familyID'}}.$parameter{$orderParameterElement}{'fileEnding'};
			    }
			    else  { #First module that should add filending
				$sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{$orderParameterElement}{'fileEnding'} = $parameter{$orderParameterElement}{'fileEnding'};
			    }
			    $tempFileEnding{$chainfork}{$scriptParameter{'familyID'}} = $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{$orderParameterElement}{'fileEnding'}; #To enable sequential build-up of fileending 
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
	&PrintToFileHandles(\@printFilehandles, "Dry Run:\n");
    }
    else { #Dry run

	$fileName = $dryRunFilenamePath;
	&PrintToFileHandles(\@printFilehandles, "Dry Run:\n");
    }

    &Checkfnexists(\$fileName, \$fnend, \$fileNameTracker);

###Info and Log
    &PrintToFileHandles(\@printFilehandles, "Creating sbatch script for ".$programName." and writing script file(s) to: ".$fileName."\n");
    &PrintToFileHandles(\@printFilehandles, "Sbatch script ".$programName." data files will be written to: ".$programDataDirectory."\n");

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
	
	if ($scriptParameter{'emailType'} =~/B/i) {

	    print $fileHandle "#SBATCH --mail-type=BEGIN", "\n";
	}
	if ($scriptParameter{'emailType'} =~/E/i) {
	 
	    print $fileHandle "#SBATCH --mail-type=END", "\n";
	}
	if ($scriptParameter{'emailType'} =~/F/i) {
	    
	    print $fileHandle "#SBATCH --mail-type=FAIL", "\n";
	}
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
		
		if($1) {
		
		    $mergeLanes = $1;
		} 
		else {
		    $mergeLanes = $2;
		} #Make sure to always supply lanes from previous regexp  
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

	$sampleInfo{$familyID}{$familyID}{'program'}{$programName}{'OutDirectory'} = $outDirectory; #OutDirectory of QC File                                                            
	if ($outDataType eq "static") { #programs which add a static file in its own directory                                                                                                 

	    $sampleInfo{$familyID}{$familyID}{'program'}{$programName}{'OutFile'} = $outFileEnding; #Static QC outFile                                                                     
	}
	if ($outDataType eq "infoDirectory") { #QC metrics are sent to info files                                                                                                                   
	    $sampleInfo{$familyID}{$familyID}{'program'}{$programName}{'OutFile'} = $outFileEnding; #info stdout file                                                                      
	}
	if ($outDataType eq "infileDependent") { #Programs which Add a filending to infile                                                                                                          
	    $sampleInfo{$familyID}{$familyID}{'program'}{$programName}{'OutFile'} = $outFileEnding; #Infile dependent QC outFile                                                                                                                                                                                       
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
	print $FILEHANDLE "java -Xmx2g ";
	if ($scriptParameter{'javaUseLargePages'} ne "no") {
	    
	    print $FILEHANDLE "-XX:-UseLargePages "; #UseLargePages for requiring large memory pages (cross-platform flag)
	}
	print $FILEHANDLE "-jar ".$scriptParameter{'picardToolsPath'}."/IntervalListTools.jar ";
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
    my $program = $_[3];
    
    my $famFile = $outFamilyFileDirectory."/".$scriptParameter{'familyID'}.".fam";
    my $parentCounter;
    my $pqParentCounter = q?perl -ne 'my $parentCounter=0; while (<>) { my @line = split(/\t/, $_); unless ($_=~/^#/) { if ( ($line[2] eq 0) || ($line[3] eq 0) ) { $parentCounter++} } } print $parentCounter; last;'?;
    my $childCounter;
    my $pqChildCounter = q?perl -ne 'my $childCounter=0; while (<>) { my @line = split(/\t/, $_); unless ($_=~/^#/) { if ( ($line[2] ne 0) || ($line[3] ne 0) ) { $childCounter++} } } print $childCounter; last;'?;
    
    $parentCounter = `$pqParentCounter $famFile`;
    $childCounter = `$pqChildCounter $famFile`;
    
    if ($program ne "GATKPhaseByTransmission") {
	
	if ($parentCounter > 0) { #Parents present
	    
	    print $FILEHANDLE "--pedigreeValidationType ".$pedigreeValidationType." --pedigree ".$outFamilyFileDirectory."/".$scriptParameter{'familyID'}.".fam "; #Pedigree files for samples		
	}
    }
    else {
	
	&CheckPedigreeMembers($FILEHANDLE, \$pedigreeValidationType, \$outFamilyFileDirectory, \$parentCounter, \$childCounter); #Special case - GATK PhaseByTransmission needs parent/child or trio 
    }
}

sub CheckPedigreeMembers {
##Detect if the pedigree file contains a valid parent/child or trio

    my $FILEHANDLE = $_[0];
    my $outFamilyFileDirectory = $_[1];
    my $pedigreeValidationType = $_[2];
    my $parentCounterRef = $_[3];
    my $childCounterRef = $_[4];
	    
    if (scalar(@sampleIDs) < 4) { #i.e.1-3 individuals in pedigree		    
		
	if ( ($childCounterRef == 1) && ($parentCounterRef > 0) ) { #parent/child or trio

	    print $FILEHANDLE "--pedigreeValidationType ".$pedigreeValidationType." --pedigree ".$outFamilyFileDirectory."/".$scriptParameter{'familyID'}.".fam "; #Pedigree files for samples
	}
	else {

	    $scriptParameter{'pGATKPhaseByTransmission'} = 0; #Override input since pedigree is not valid for analysis
	    &PrintToFileHandles(\@printFilehandles, "Switched GATK PhaseByTransmission to 'no run' mode since MIP did not detect a valid pedigree for this type of analysis.");
	    
	    if ($scriptParameter{'pGATKReadBackedPhasing'} > 0) { #Broadcast
		
		&PrintToFileHandles(\@printFilehandles, "MIP will still try to run GATK ReadBackedPhasing, but with the '-respectPhaseInInput' flag set to false\n");
	    }
	    print STDOUT "\n";
	}
    }
    else {
	
	$scriptParameter{'pGATKPhaseByTransmission'} = 0; #Override input since pedigree is not valid for analysis
	&PrintToFileHandles(\@printFilehandles, "Switched GATK PhaseByTransmission to 'no run' mode since MIP did not detect a valid pedigree for this type of analysis.");
	
	if ($scriptParameter{'pGATKReadBackedPhasing'} > 0) { #Broadcast
	    
	    &PrintToFileHandles(\@printFilehandles, "MIP will still try to run GATK ReadBackedPhasing, but with the '-respectPhaseInInput' flag set to false\n");
	}
	print STDERR "\n";
    }
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
    print MIPLOG "MIP Version: ".$MipVersion, "\n";
    #Note FileHandle MIPLOG not closed
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

    open (YAML, "<". $yamlFile) or die "can't open ".$yamlFile.":".$!, "\n";    
        
    %yamlHash = %{ YAML::LoadFile($yamlFile) }; #Load hashreference as hash
        
    close(YAML);
    
    print STDOUT "Read Yaml file: ". $yamlFile, "\n";
    return %yamlHash;
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
##Test that the familyID and the sampleID(s) exists and are unique. Check id sampleID contains "_".

    my %seen; #Hash to test duplicate sampleIDs later

    if (scalar(@sampleIDs) == 0) {

	print STDOUT "\nPlease provide sampleID(s)\n\n";
	exit;
    }

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@sampleIDs);$sampleIDCounter++) {

	$seen{ $sampleIDs[$sampleIDCounter] }++; #Increment instance to check duplictaes later
	
	if ($scriptParameter{'familyID'} eq $sampleIDs[$sampleIDCounter]) {
	    
	    print STDERR "\nFamilyID: ".$scriptParameter{'familyID'}." equals sampleID: ".$sampleIDs[$sampleIDCounter].". Please make sure that the familyID and sampleID(s) are unique.\n";
	    exit;
	}
	if ($seen{ $sampleIDs[$sampleIDCounter] } > 1) {
	
	    print STDERR "\nSampleID: ".$sampleIDs[$sampleIDCounter]." is not uniqe.\n\n";
	    exit;
	}
	if ($sampleIDs[$sampleIDCounter] =~/_/) { #SampleID contains "_", which is not allowed accrding to filename conventions

	    print STDERR "\nSampleID: ".$sampleIDs[$sampleIDCounter]." contains '_'. Please rename sampleID according to MIP's filename convention, removing the '_'.\n\n";
	    exit;	    
	}
    }
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
	
	if (defined($clusterConstantPath)) {
	 
	    $scriptParameter{$orderParameterElement} =~ s/CLUSTERCONSTANTPATH!/$clusterConstantPath/gi; #Exchange CLUSTERCONSTANTPATH! for current cluster path
	}
	if (defined($analysisConstantPath)) {
	 
	    $scriptParameter{$orderParameterElement} =~ s/ANALYSISCONSTANTPATH!/$analysisConstantPath/gi; #Exchange ANALYSISCONSTANTPATH! for the current analysis path
	}
	if (defined($analysisType)) {

	    $scriptParameter{$orderParameterElement} =~ s/ANALYSISTYPE!/$analysisType/gi; #Exchange ANALYSISTYPE! for the current analysis type
	}
	if (defined($familyID)) {

	    $scriptParameter{$orderParameterElement} =~ s/FDN!/$familyID/gi; #Exchange FND! for the current familyID
	}
	if (defined($aligner)) {

	    $scriptParameter{$orderParameterElement} =~ s/ALIGNER!/$aligner/gi; #Exchange ALIGNER! for the current aligner
	}
    }
}

sub CheckAutoBuild {
##Checks if autobuild is on and returns "1" if enabled or "0" if not

    my $parameterNameRef = $_[0];
    my $sampleIDRef = $_[1];

    if (defined($sampleIDRef)) {
	
	if ( ($parameter{ $scriptParameter{'familyID'} }{$$sampleIDRef}{$$parameterNameRef}{'buildFile'} eq "yesAutoBuild") || ($parameter{ $scriptParameter{'familyID'} }{$$sampleIDRef}{$$parameterNameRef}{'buildFile'} eq "yesAutoDownLoad") || ($parameter{ $scriptParameter{'familyID'} }{$$sampleIDRef}{$$parameterNameRef}{'buildFile'} eq 1) ) {
	    
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
	elsif ( ($parameter{$$parameterNameRef}{'buildFile'} eq "yesAutoDownLoad") || ($parameter{$$parameterNameRef}{'buildFile'} eq 1) ) {
	    
	    if ($parameter{$$parameterNameRef}{'default'} eq $scriptParameter{$$parameterNameRef}) {
		
		return "1"; #Flag that autobuild is needed
	    }
	    else {
		
		print STDERR "\nCould not find file ".$scriptParameter{'referencesDir'}."/".$scriptParameter{$$parameterNameRef}, "\n";
		print STDERR "\nMake sure that file exists or use the default for this parameter to enable automatic download via Cosmid", "\n";
		exit;
	    }
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
}

sub SetAutoBuildFeature {

##SetAutoBuildFeature
    
##Function : Sets parameters with autoBuild enabled to the new value dependent on $referenceFileNameRef
##Returns  : "" 
##Arguments: $parameterName, $referenceFileEndingRef, $referenceFileNameRef, $printSwitch
##         : $parameterName          => MIP parameter name 
##         : $referenceFileEndingRef => Reference file name ending {REF}
##         : $referenceFileNameRef   => Reference file name {REF}
##         : $printSwitch            => To print or not

     my $parameterName = $_[0];
     my $referenceFileEndingRef = $_[1];
     my $referenceFileNameRef = $_[2];
     my $printSwitch = $_[3];
     
     if( defined($scriptParameter{$parameterName}) && ($scriptParameter{$parameterName} eq "notSetYet") ) { 

	 $scriptParameter{$parameterName} = $$referenceFileNameRef.$$referenceFileEndingRef;

	 if ( (defined($printSwitch)) && ($printSwitch ne "noPrint") ) {

	     print STDOUT "Set ".$parameterName." to: ".$scriptParameter{$parameterName}, "\n";
	 }
	 if ($parameterName eq "bwaBuildReference") {

	     &CheckFileEndingsToBeBuilt(\@bwaBuildReferenceFileEndings, "bwaBuildReference");
	 }
	 elsif ($parameterName eq "mosaikJumpDbStub") {

	     &CheckFileEndingsToBeBuilt(\@mosaikJumpDbStubFileEndings, "mosaikJumpDbStub");
	 }
	 else {#Complete fileName - No stubs
	    
	     &CheckExistance(\($scriptParameter{'referencesDir'}."/".$scriptParameter{$parameterName}), \$parameterName, "f");
         }
    }
}

sub MoveMosaikNN {

##MoveMosaikNN

##Function : Locate MOSAIK path and move neural network files in place if lacking
##Returns  : "" 
##Arguments: 

    my @paths = split(/:/,$ENV{PATH}); #Find Mosaik installation path

    for (my $pathsCounter=0;$pathsCounter<scalar(@paths);$pathsCounter++) {

	if ($paths[$pathsCounter] =~/MOSAIK/) { #Select MOSAIK path
	    
	   $paths[$pathsCounter] =~ s/bin\//src\/networkFile/g; #Location of NN files

	   print STDOUT "\nCould not find Mosaik Network Files in ".$scriptParameter{'referencesDir'},"\n";
	   print STDOUT "\nCopying Mosaik Network Files ".$scriptParameter{'mosaikAlignNeuralNetworkSeFile'}." and ".$scriptParameter{'mosaikAlignNeuralNetworkPeFile'}." to ".$scriptParameter{'referencesDir'}." from ".$paths[$pathsCounter], "\n\n";
	   `cp $paths[$pathsCounter]/$scriptParameter{'mosaikAlignNeuralNetworkSeFile'} $scriptParameter{'referencesDir'}/`; #Copying files in place
	   `cp $paths[$pathsCounter]/$scriptParameter{'mosaikAlignNeuralNetworkPeFile'} $scriptParameter{'referencesDir'}/`; #Copying files in place
	   last;
	}
    }
}

sub CheckUserInfoArrays {

##CheckUserInfoArrays
    
##Function : Determine if the user supplied info on array parameter
##Returns  : "0|1" 
##Arguments: $arrayRef, $parameterName
##         : $arrayRef      => Array to loop in for parameter {REF}
##         : $parameterName => MIP parameter to evaluate
    
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
    
##SetTargetFiles
    
##Function : Sets the target files and replaces constant genomic information with the info used in present analysis. Adds target file to sampleInfo for qc print later.
##Returns  : "" 
##Arguments: $familyIDRef, $sampleIDRef, $parameterNameRef, $referenceFileEndingRef
##         : $familyIDRef            => Family ID {REF}
##         : $sampleIDRef            => Sample ID  {REF}
##         : $parameterName          => MIP parameter name
##         : $referenceFileEndingRef => File name ending {REF}
    
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

##PrepareArrayParameters
    
##Function : Check if user supplied cmd info and supplies arrayParameters to scriptParameters
##Returns  : "" 
##Arguments: $arrayRef, $parameterName, $parameterType, $parameterDefault, $associatedPrograms, $parameterExistsCheck
##         : $arrayRef             => Array to loop in for parameter {REF}
##         : $parameterName        => MIP parameter to evaluate
##         : $parameterType        => Type of MIP parameter 
##         : $parameterDefault     => The parameter default value
##         : $associatedPrograms   => Programs that use the parameter. Comma separated string
##         : $parameterExistsCheck => Check if intendent file exists in reference directory

    my $arrayRef = $_[0];
    my $parameterName = $_[1];
    my $parameterType = $_[2];
    my $parameterDefault = $_[3];
    my $associatedPrograms = $_[4]; 
    my $parameterExistsCheck = $_[5];

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

##CheckUniqueTargetFiles
    
##Function : Checks target files within parameters for identical entries
##Returns  : "" 
##Arguments: $arrayRef, $countRef, $fileToCompareRef, $parameterName
##         : $arrayRef         => Array to loop in for parameter (e.g. sampleID) {REF}
##         : $countRef         => Offset in array {REF}
##         : $fileToCompareRef => The file to compare against rest of array {REF}
##         : $parameterName    => MIP parameter to evaluate
    
    my $arrayRef = $_[0]; 
    my $countRef = $_[1]; 
    my $fileToCompareRef = $_[2];
    my $parameterName = $_[3];
    
    for (my $compareCounter=($$countRef + 1);$compareCounter<scalar(@{$arrayRef});$compareCounter++) { #Compare all target files to remove autoBuild if file names are identical for each flag
	
	if ( (defined($$fileToCompareRef)) && (defined($scriptParameter{ $scriptParameter{'familyID'} }{ ${$arrayRef}[$compareCounter] }{$parameterName})) ) {

	    if ($$fileToCompareRef eq $scriptParameter{ $scriptParameter{'familyID'} }{ ${$arrayRef}[$compareCounter] }{$parameterName}) { #Identical target files
	    
		$parameter{ $scriptParameter{'familyID'} }{ ${$arrayRef}[$compareCounter] }{$parameterName}{'buildFile'} = 0; #Turn off autoBuild
	    }
	}
    }
}

sub CheckTemplateFilesPaths {

##CheckTemplateFilesPaths
    
##Function : Checks that file paths in template files exist
##Returns  : "" 
##Arguments: $fileNameRef, $parameterName
##         : $fileNameRef   => File name {REF}
##         : $parameterName => MIP parameter name

    my $fileNameRef = $_[0];
    my $parameterName = $_[1];

    open(TF, "<".$$fileNameRef) or die "Can't open ".$$fileNameRef.":$!, \n";  

    while (<TF>) {

	chomp $_;

	if (m/^\s+$/) {	# Avoid blank lines
            next;
        }
	if (m/^\#/) {	# Avoid "#"
            next;
        }	
	if ($_ =~/(\S+)/) {	
	    
	    my $filePath = $_;

	    if ($filePath=~/^(RD!)/) { #intersectCollect file
		
		my @filePath = split('\t', $filePath);
		$filePath[0] =~ s/^RD!/$scriptParameter{'referencesDir'}/g;

		&CheckExistance(\$filePath[0], \$parameterName, "f"); #Only check paths that pointing to reference directory
	    }
	    if ($parameterName eq "GATKHaploTypeCallerRefBAMInfile") { #Only Paths should be present i.e. check all lines
	       
		&CheckExistance(\$filePath, \$parameterName, "f");
	    }
	}
    }
    close(TF);
}

sub CheckGzipped {

##CheckGzipped
    
##Function : Check if a file is gzipped.
##Returns  : "unCompressed|compressed" 
##Arguments: $fileNameRef
##         : $fileNameRef => File name {REF}

    my $fileNameRef = $_[0];

    my $fileCompressionStatus = "unCompressed";

    if ( (defined($$fileNameRef)) && ($$fileNameRef =~/.gz$/) ) {
	
	$fileCompressionStatus = "compressed"; 
    }
    return $fileCompressionStatus;
}

sub RemoveFileEnding {

##RemoveFileEnding
    
##Function : Removes ".fileEnding" in filename.FILENDING(.gz)
##Returns  : File name with supplied $fileEnding(.gz) removed
##Arguments: $fileNameRef
##         : $fileNameRef => File name {REF}
##         : $fileEnding  => File ending to be removed

    my $fileNameRef = $_[0];
    my $fileEnding = $_[1];

    my $fileNameNoEnding;
    
    if ( (defined($$fileNameRef)) && ($$fileNameRef =~/(\S+)($fileEnding$|$fileEnding.gz$)/) ) {

	$fileNameNoEnding = $1;
    }
    return $fileNameNoEnding;
}

sub ScriptParameterPerSampleID {

##ScriptParameterPerSampleID
    
##Function : Enables files handled per SampleID to be processed by AddToScriptParameters
##Returns  : ""
##Arguments: $familyIDRef, $sampleIDRef, $parameterName
##         : $familyIDRef   => Family ID {REF}
##         : $sampleIDRef   => Sample ID  {REF}
##         : $parameterName => MIP parameter name

    my $familyIDRef = $_[0];
    my $sampleIDRef = $_[1];
    my $parameterName = $_[2];
    
    if (defined($scriptParameter{$$familyIDRef}{$$sampleIDRef}{$parameterName})) {
	
	$scriptParameter{$parameterName} = 1; #Define in scriptParameter so that we now that parameter is present per SampleID
    }
}

sub EnableArrayParameter {

##EnableArrayParameter
    
##Function : Adds arrayRef to scriptParameters for recreation of cmd in log and seperated input parameter string into array elements
##Returns  : ""
##Arguments: $arrayRef, $parameterNameRef
##         : $arrayRef         => Array to loop through {REF}
##         : $parameterNameRef => MIP parameter name {REF}

    my $arrayRef = $_[0];
    my $parameterNameRef = $_[1];
    
    $scriptParameter{$$parameterNameRef} = join(',',@{$arrayRef}); #Add to enable recreation of cmd line later
    @{$arrayRef} = split(/,/,join(',', @{$arrayRef})); #Enables comma separated list of sample IDs from user supplied cmd info
}

sub SetTargetandAutoBuild {

##SetTargetandAutoBuild
    
##Function : Set autoBuild for target files and calls SetTargetFile subroutine
##Returns  : ""
##Arguments: $arrayRef, $parameterNameRef, $fileEndingRef
##         : $arrayRef         => Array to loop through {REF}
##         : $parameterNameRef => MIP parameter name {REF}
##         : $fileEndingRef    => File ending {REF}

    my $arrayRef = $_[0];
    my $parameterNameRef = $_[1];
    my $fileEndingRef = $_[2];
    
    for (my $elementsCounter=0;$elementsCounter<scalar(@{$arrayRef});$elementsCounter++) {
	
	$parameter{ $scriptParameter{'familyID'} }{${$arrayRef}[$elementsCounter]}{$$parameterNameRef}{'buildFile'} = "yesAutoBuild";
	&SetTargetFiles(\$scriptParameter{'familyID'}, \${$arrayRef}[$elementsCounter], \$$parameterNameRef, \$$fileEndingRef);
    }   
}

sub CheckTargetExistFileBed {

##CheckTargetExistFileBed
    
##Function : Check that supplied target file ends with ".bed" and exists.
##Returns  : ""
##Arguments: $fileRef, $parameterName
##         : $fileRef        => File to check for existance and file ending {REF}
##         : $parameterName  => MIP parameter name

    my $fileRef = $_[0];
    my $parameterName = $_[1];

    if (defined($$fileRef)) {

	if ($$fileRef !~/.bed$/) {
	
	print STDERR "\nCould not find intendended 'file ending with .bed' for target file: ".$$fileRef." in parameter '-".$parameterName."'", "\n";
	exit;
	}
	unless (-f $scriptParameter{'referencesDir'}."/".$$fileRef) {

	    print STDERR "\nCould not find intendended '.bed' file for target file: ".$scriptParameter{'referencesDir'}."/".$$fileRef." in parameter '-".$parameterName."'", "\n\n";
	    exit;
	}
    }
}

sub CompareArrayElements {

##CompareArrayElements
    
##Function : Compares the number of elements in two arrays and exits if the elements are not equal.
##Returns  : ""
##Arguments: $arrayRef, $arrayQueryRef, $parameterName, $parameterNameQuery
##         : $arrayRef           => Array to match {REF}
##         : $arrayQueryRef      => Array to be compared {REF}
##         : $parameterName      => MIP reference parameter
##         : $parameterNameQuery => MIP query parameter

    my $arrayRef = $_[0]; 
    my $arrayQueryRef = $_[1]; 
    my $parameterName = $_[2];
    my $parameterNameQuery = $_[3];
    
    if (scalar(@{$arrayRef}) != scalar(@{$arrayQueryRef})) {
	
	print STDERR "\nThe supplied '-".$parameterNameQuery."' lists do not equal the number of elements in '-".$parameterName."'. Please specify a equal number of elements in both lists", "\n\n";
	exit;
    }
}

sub SetAutoBuildAndScriptParameterPerSample {

##SetAutoBuildAndScriptParameterPerSample
    
##Function : Sets autoBuild and populates scriptParameter hash with array elements per sampleID.
##Returns  : ""
##Arguments: $sampleIDArrayRef, $parameterArrayRef, $parameterNameRef
##         : $sampleIDArrayRef  => SampleID array {REF}
##         : $parameterArrayRef => Parameter array {REF}
##         : $parameterNameRef  => MIP parameter name {REF}

    my $sampleIDArrayRef = $_[0];
    my $parameterArrayRef = $_[1];
    my $parameterNameRef = $_[2];

    for (my $sampleIDsCounter=0;$sampleIDsCounter<scalar(@{$sampleIDArrayRef});$sampleIDsCounter++) { #All sampleIDs
	
	$parameter{ $scriptParameter{'familyID'} }{${$sampleIDArrayRef}[$sampleIDsCounter]}{$$parameterNameRef}{'buildFile'} = "yesAutoBuild"; #Turn on autoBuild
	$scriptParameter{ $scriptParameter{'familyID'} }{ ${$sampleIDArrayRef}[$sampleIDsCounter] }{$$parameterNameRef} = ${$parameterArrayRef}[$sampleIDsCounter]; #Populate hash that is used in modules
    }
}

sub SetTargetFileGeneralBuildParameter {

##SetTargetFileGeneralBuildParameter 
    
##Function : Sets the general build parameters $$sampleIDBuildFileRef and $$sampleIDBuildFileNoEndingRef and sets buildfile key to "0".
##Returns  : ""
##Arguments: $targetfileRef, $parameterName, $sampleIDBuildFileRef, $sampleIDBuildFileNoEndingRef, $sampleIDRef
##         : $targetfileRef                => Final file {REF}
##         : $parameterName                => MIP parameter
##         : $sampleIDBuildFileRef         => File that will be created {REF}
##         : $sampleIDBuildFileNoEndingRef => File that will be created with file ending removed {REF}
##         : $sampleIDRef                  => SampleID {REF}

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

##PrintCheckExistandMoveFile
    
##Function : Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
##Returns  : ""
##Arguments: $FILEHANDLE, $intendedFilePathRef, $temporaryFilePathRef
##         : $FILEHANDLE           => FILEHANDLE to write to
##         : $intendedFilePathRef  => Path to file to check for existence {REF}
##         : $temporaryFilePathRef => File that has been created {REF}
    
    my $FILEHANDLE = $_[0]; 
    my $intendedFilePathRef = $_[1]; 
    my $temporaryFilePathRef = $_[2];
    
    print $FILEHANDLE "[ -s ".$$intendedFilePathRef." ] "; #Check file exists and is larger than 0
    print $FILEHANDLE "&& rm ".$$temporaryFilePathRef." "; #If other processes already has created file, remove temp file
    print $FILEHANDLE "|| "; #File has not been created by other processes
    print $FILEHANDLE "mv ".$$temporaryFilePathRef." ".$$intendedFilePathRef,"\n\n"; #Move file in place
    
}

sub DefineAnnovarTables {
    
##DefineAnnovarTables
    
##Function : Defines and adds annovar tables parameters to hash
##Returns  : ""
##Arguments: 

    my $annovarGenomeBuildVersion = $scriptParameter{'annovarGenomeBuildVersion'}; #Set the current annovar genome build

    ##Define annovar tables
    my @annovarTablesGeneAnno = ("refGene", "knownGene", "ensGene"); #Tables using annotation option "geneanno"
    my @annovarTablesRegionAnno = ("mce46way", "gerp++elem", "segdup", "tfbs", "mirna"); #Tables using annotation option "regionanno"
    my @annovarTablesFilter = ("snp137", "snp135", "snp132", "snp131", "snp130", "snp129", "snp137NonFlagged", "snp135NonFlagged", "snp132NonFlagged", "snp131NonFlagged", "snp130NonFlagged", "1000g2012apr_all", "1000g2012apr_amr", "1000g2012apr_eur", "1000g2012apr_asn", "1000g2012apr_afr", "1000g2012feb_all", "esp6500si_all", "esp6500_all", "esp6500_aa", "esp6500_ea", "esp5400_all", "esp5400_aa", "esp5400_ea","clinvar_20131105", "ljb2_sift", "ljb2_pp2hdiv", "ljb2_pp2hvar", "ljb2_mt", "ljb2_ma", "ljb2_fathmm", "ljb2_siphy", "ljb2_lrt", "ljb_all", "ljb2_gerp++", "ljb2_phylop", "caddgt20", "caddgt10"); #Tables using annotation option "filter"
    my @annovarTablesUrlUcsc = ("mce46way", "segdup", "tfbs", "mirna"); #Tables using urlAlias "ucsc"
    my @annovarGenericFiltering = ("esp6500si_all", "esp6500_all", "esp6500_aa", "esp6500_ea", "esp5400_all", "esp5400_aa", "esp5400_ea","clinvar_20131105"); #Tables using generic option
    my @annovarGenericFiles = ($annovarGenomeBuildVersion."_esp6500si_all.txt", $annovarGenomeBuildVersion."_esp6500_all.txt", $annovarGenomeBuildVersion."_esp6500_aa.txt", $annovarGenomeBuildVersion."_esp6500_ea.txt", $annovarGenomeBuildVersion."_esp5400_all.txt", $annovarGenomeBuildVersion."_esp5400_aa.txt", $annovarGenomeBuildVersion."_esp5400_ea.txt", $annovarGenomeBuildVersion."_clinvar_20131105.txt"); #Generic table files
    my @annovarRefgeneFiles = ($annovarGenomeBuildVersion."_refGene.txt", $annovarGenomeBuildVersion."_refGeneMrna.fa", $annovarGenomeBuildVersion."_refLink.txt", "GRCh37_MT_ensGene.txt", "GRCh37_MT_ensGeneMrna.fa"); #Cater for multiple download
    my @annovarKnownGeneFiles = ($annovarGenomeBuildVersion."_knownGene.txt", $annovarGenomeBuildVersion."_kgXref.txt", $annovarGenomeBuildVersion."_knownGeneMrna.fa"); #Cater for multiple download
    my @annovarEnsGeneFiles = ($annovarGenomeBuildVersion."_ensGene.txt", $annovarGenomeBuildVersion."_ensGeneMrna.fa", "GRCh37_MT_ensGene.txt", "GRCh37_MT_ensGeneMrna.fa"); #Cater for multiple download

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
	&AnnovarTableParameters(\$annovarTablesRegionAnno[$tablesCounter], \@annovarTablesUrlUcsc, "urlAlias", "ucsc"); #Replace for ucsc tables NOTE: not all in RegionAnno
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
   
##AnnovarTableParameters
    
##Function : Populates annovarTable hash by looping through array and adding table with identified membership and associated parameters into hash.
##           Parameters of type "file" are added as a hash of array since multiple files can be downloaded

##Returns  : ""
##Arguments: $tableNameRef, $arrayRef, $parameterType, $parameterValue
##         : $tableNameRef   => Annovar table name {REF}
##         : $arrayRef       => Array to search for membership {REF}
##         : $parameterType  => Type of table parameter
##         : $parameterValue => Parameter value
    
    my $tableNameRef = $_[0];
    my $arrayRef = $_[1];
    my $parameterType = $_[2];
    my $parameterValue = $_[3];
    
    for (my $tablesCounter=0;$tablesCounter<scalar(@{$arrayRef});$tablesCounter++) {
    
	if (${$arrayRef}[$tablesCounter] eq $$tableNameRef) { #Membership test
	    
	    if ($parameterType eq "file") { #Add as array instead, since some annovar tables have multiple files downloaded for the same call
		
		push(@{$annovarTables{$$tableNameRef}{$parameterType}}, $parameterValue); #Add as array
	    }
	    else {
		
		$annovarTables{$$tableNameRef}{$parameterType} = $parameterValue;
	    }
	    last; #No table should be represented twice within the same array
	}
    }
}

sub CollectSeqContigs {

##CollectSeqContigs
    
##Function : Collects sequences contigs used in analysis from human genome sequence dictionnary associated with $humanGenomeReference
##Returns  : ""
##Arguments: 

    my $pqSeqDict = q?perl -nae 'if($F[0]=~/^\@SQ/) { if($F[1]=~/SN\:(\S+)/) {print $1, ",";} }' ?; 
    my $SeqDictLocation = $scriptParameter{'referencesDir'}."/".$humanGenomeReferenceNameNoEnding.".dict";
    @contigs = `$pqSeqDict $SeqDictLocation `; #returns a comma seperated string of sequence contigs from dict file
    @contigs = split(/,/,join(',', @contigs));
}

sub ReplaceConfigParamWithCMDInfo {

##ReplaceConfigParamWithCMDInfo
    
##Function : Replace config parameter with cmd info for active parameter
##Returns  : Path to Cosmid Resource directory for current analysis
##Arguments: $parameterName
##         : $parameterName => MIP parameter name

    my $parameterName = $_[0];

    if ($parameter{$parameterName}{'value'} ne "nocmdinput") { #Replace config parameter with cmd info for parameter
	
	$scriptParameter{$parameterName} = $parameter{$parameterName}{'value'}; #Transfer to active parameter
    }
}

sub DefineSupportedCosmidReferences {

##DefineSupportedCosmidReferences
    
##Function : Defines the Cosmid manager hash keys and populates it from arguments 
##Returns  : ""
##Arguments: $parameterName, $cosmidResourceName, $CosmidResourceVersion, 
##         : $parameterName                  => MIP parameter name
##         : $cosmidResourceName             => Cosmid Resource name
##         : $CosmidResourceVersion          => Version of the cosmid Resource to download
##         : $humanGenomeReferenceVersionRef => The human genome build used in the analysis
##         : $compressedSwitch               => If files after download are compressed or not

    my $parameterName = $_[0];
    my $cosmidResourceName = $_[1];
    my $CosmidResourceVersion = $_[2];
    my $humanGenomeReferenceVersionRef = $_[3];
    my $compressedSwitch = $_[4];

    $supportedCosmidReferences{$parameterName} = {

	'cosmidName' => $cosmidResourceName,
	'version' => $CosmidResourceVersion,
	'humanGenomeReferenceVersion' => $$humanGenomeReferenceVersionRef,
	'compressedSwitch' => $compressedSwitch,
    };
}

sub CheckCosmidYAML {

##CheckCosmidYAML
    
##Function : Locates and sets the cosmid directory to download to
##Returns  : Path to Cosmid Resource directory for current analysis
##Arguments: 

    my %cosmidResources; #Hash to load cosmid info to

    if (-f $scriptParameter{'referencesDir'}."/cosmid.yaml") { #cosmid.yaml file exists in reference directory
	
	%cosmidResources = &LoadYAML($scriptParameter{'referencesDir'}."/cosmid.yaml"); #Load yaml file
	
	unless (defined($cosmidResources{'directory'})) { #Set Directory entry if not defined
	    
	    $cosmidResources{'directory'} = $scriptParameter{'referencesDir'}."/resources"; #Set the Cosmid default directory
	}
    }
    else { #No cosmid.yaml exist in reference directory
	
	$cosmidResources{'directory'} = $scriptParameter{'referencesDir'}."/resources"; #Set the Cosmid default directory
    }
    return $cosmidResources{'directory'};
}

sub AdjustNrCoresToSeqMode {

##AdjustNrCoresToSeqMode
    
##Function : Adjust the number of cores to be used in the analysis according to sequencing mode requirements
##Returns  : Adds to $$nrCoresRef
##Arguments: $nrCoresRef, $sequenceRunTypeRef
##         : $nrCoresRef         => The maximum number of cores to be use before printing "wait" statement {REF}
##         : $sequenceRunTypeRef => Type of sequencing [Paired-end|Single-end] {REF}

    my $nrCoresRef = $_[0];
    my $sequenceRunTypeRef = $_[1];
    
    if ($$sequenceRunTypeRef eq "Paired-end") { #Second read direction if present

	$$nrCoresRef =  $$nrCoresRef + 2; #2 processes per file
    }
    else {#Single-end

	$$nrCoresRef = $$nrCoresRef + 1; #Only 1 file and one process
    }
}

sub PrintWait {

##PrintWait
    
##Function : Calculates when to prints "wait" statement and prints "wait" to supplied FILEHANDLE when adequate. 
##Returns  : Incremented $$coreCounterRef
##Arguments: $counterRef, $nrCoresRef, $coreCounterRef
##         : $counterRef     => The number of used cores {REF}
##         : $nrCoresRef     => The maximum number of cores to be use before printing "wait" statement {REF}
##         : $coreCounterRef => Scales the number of $nrCoresRef cores used after each print "wait" statement {REF}
##         : $FILEHANDLE     => FILEHANDLE to print "wait" statment to
    
    my $counterRef = $_[0];
    my $nrCoresRef = $_[1];
    my $coreCounterRef = $_[2];
    my $FILEHANDLE = $_[3];
    
    if ($$counterRef == $$coreCounterRef * $$nrCoresRef) { #Using only nr of cores eq to lanes or maximumCores
	
	print $FILEHANDLE "wait", "\n\n";
	$$coreCounterRef=$$coreCounterRef+1; #Increase the maximum number of cores allowed to be used since "wait" was just printed 
    }    
}

sub CheckBuildDownLoadPreRequisites {

##CheckBuildDownLoadPreRequisites
    
##Function : Checks if some of the active prerequisites needs to be downloaded and calls subroutine to download them if required
##Returns  : ""
##Arguments: $programName
##         : $programName => Active program

    my $programName = $_[0];
    
    for my $parameterName (keys %supportedCosmidReferences) { #Supported cosmid references for MIP parameters
	
	if ( "p".$programName eq $parameter{$parameterName}{'associatedProgram'}) { #If the cosmid supported parameter is associated with the MIP program
	    
	    if ($parameter{$parameterName}{'buildFile'} eq 1) { #Enable autoBuild
		
		&BuildDownLoadablePreRequisites($scriptParameter{'familyID'}, $scriptParameter{'aligner'}, $programName);
		last; #Perform once
	    }
	}
    }
}

sub PrintToFileHandles {

##PrintToFileHandles
    
##Function : Prints statement to all open FILEHANDLES 
##Returns  : ""
##Arguments: $arrayRef, $statement
##         : $arrayRef   => Array holding FILEHANDLES to print to {REF}. 
##         : $statement  => Statement to be printed

    my $arrayRef = $_[0];
    my $statement = $_[1];
   
    for (my $filehandleCounter=0;$filehandleCounter<scalar(@{$arrayRef});$filehandleCounter++) {
	 
	print {${$arrayRef}[$filehandleCounter]} $statement;
    }
}

sub CheckSupportedFileEnding {

##CheckSupportedFileEnding    

##Function : Check that the supplied fileEnding is supported. Otherwise exits. 
##Returns  : ""
##Arguments: $fileNameRef, $fileEndingRef, $parameterNameRef
##         : $fileNameRef      => File name {REF}. 
##         : $fileEndingRef    => Supported file name ending {REF}
##         : $parameterNameRef => Parameter name {REF}

    my $fileNameRef = $_[0];
    my $fileEndingRef = $_[1];
    my $parameterNameRef = $_[2];

    unless ( (defined($$fileNameRef)) && ($$fileNameRef =~/(\S+)$$fileEndingRef$/) ) {
	
	print STDERR "\nThe supplied file: ".$$fileNameRef." for parameter '".$$parameterNameRef."' does not have the supported file ending '".$$fileEndingRef."'.", "\n\n";
	exit;
    }
}

sub PrintSupportedAnnovarTableNames {

##PrintSupportedAnnovarTableNames
    
##Function : Print the by MIP supported Annovar Table names to STDOUT and exists
##Returns  : ""
##Arguments: None
    
    print STDOUT "\nThese Annovar databases are supported by MIP:\n";
    
    foreach my $annovarSupportedTableName (@annovarSupportedTableNames) {
	
	print STDOUT $annovarSupportedTableName, "\n";
    }
    print STDOUT "\n";
    exit;
}

sub CheckEntryHashofArray {

##CheckEntryHashofArray
    
##Function : Test element for being part of hash of array at supplied key. 
##Returns  : Return "1" if element is not part of array
##Arguments: $hashRef, $keyRef, $elementRef
##         : $hashRef    => Hash {REF}. 
##         : $keyRef     => The key pointing to the array in the $hashRef {REF}
##         : $elementRef => Element to look for in hash of array {REF}

    my $hashRef = $_[0];
    my $keyRef = $_[1];
    my $elementRef = $_[2];

    if (defined($$hashRef{$$keyRef})) { #Information on entry present

	if ( ! ( grep /$$elementRef/, @{$$hashRef{$$keyRef}} ) ) { #If element is not part of array

	    return 1;
	}
    }
}

sub CheckMostCompleteAndRemoveFile {
    
##CheckMostCompleteAndRemoveFile
	
##Function  : Checks if the file is recorded as the "MostCompleteBAM|VCF". If false writes removal of file(s) to supplied filehandle
##Returns   : ""
##Arguments : $FILEHANDLE, $mostCompleteRef, $fileRef, $fileEnding
##          : $FILEHANDLE      => SBATCH script FILEHANDLE to print to. 
##          : $mostCompleteRef => The mostComplete file (BAM|VCF) {REF}
##          : $fileRef         => Current file {REF}
##          : $fileEnding      => File ending of $fileRef
    
    my $FILEHANDLE = $_[0];
    my $mostCompleteRef = $_[1];
    my $fileRef = $_[2];
    my $fileEnding = $_[3];
    
    if ( (defined($$mostCompleteRef)) && (defined($$fileRef)) ) { #Not to disturb first dry_run of analysis

	unless ($$mostCompleteRef eq $$fileRef) { #Do not remove mostCompleteBAM|VCF
		
	    my $fileName = &RemoveFileEnding(\$$fileRef, $fileEnding);

	    if (defined($fileName)) { #Successfully removed file ending using &RemoveFileEnding
	    
		my $end = ".*"; #Remove all files with ending with ".*"

		if ($fileEnding eq ".bam") { #For BAM files 
		
		    $end = ".ba*"; #Removes both .bam and .bai
		}
		if ($fileEnding eq ".vcf") { #For VCF files
		
		    $end = ".vcf*"; #Removes both .vcf and .vcf.idx
		}
		##Print removal of file to sbatch script 
		print $FILEHANDLE "rm ";
		print $FILEHANDLE $fileName.$end, "\n\n"; #Remove file(s)
	    }
	}
    }
}

####
#Decommissioned
####
