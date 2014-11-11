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

##Third party module(s)
use YAML;
use Log::Log4perl;

use vars qw($USAGE);

BEGIN {

    ##Check YAML dependecy
    eval { 

	require YAML; 
    };
    if($@) {

	print STDERR "NOTE: YAML not installed - Please install to run MIP.\n";
	print STDERR "NOTE: Aborting!\n";
	exit
    }
    ##Check LOG4perl dependency
    eval { 

	require Log::Log4perl; 
    };
    if($@) {

	print STDERR "NOTE: Log::Log4perl not installed - Please install to run MIP.\n";
	print STDERR "NOTE: Aborting!\n";
	exit
    }

    $USAGE =
	qq{
mip.pl  -ifd [inFilesDirs,.,.,.,n] -isd [inScriptDir,.,.,.,n] -rd [refdir] -p [project ID] -s [sample ID,.,.,.,n] -em [e-mail] -osd [outdirscripts] -odd [outDataDir] -f [familyID] -p[program]
               ####MIP
	       -ifd/--inFilesDirs Infile directory(s) (comma sep; mandatory: supply whole path,)
               -isd/--inScriptDir The pipeline custom script in directory (mandatory: supply whole path)
               -rd/--referencesDir Reference(s) directory (mandatory: supply whole path)
	       -p/--projectID The project ID  (mandatory)
	       -s/--sampleIDs The sample ID(s)(comma sep; mandatory)
	       -em/--email E-mail (defaults to "")
               -emt/--emailType E-mail type (defaults to F (=FAIL);Options: B (=BEGIN) and/or F (=FAIL) and/or E=(END))
	       -odd/--outDataDir The data files output directory (mandatory: supply whole path)
	       -osd/--outScriptDir The script files (.sh) output directory (mandatory: supply whole path)
               -f/--familyID Group id of samples to be compared (defaults to "0" (=no), (Ex: 1 for IDN 1-1-1A))
               -ped/--pedigreeFile (defaults to ""; supply whole path)
               -hgr/--humanGenomeReference Fasta file for the human genome reference (defaults to "Homo_sapiens.GRCh37.d5.fasta;1000G decoy version 5")
               -al/--aligner Setting which aligner was used for alignment in previous analysis (defaults to "")
               -at/--analysisType Type of analysis to perform (defaults to "exomes";Valid entries: "genomes", "exomes", "rapid")
               -mc/--maximumCores The maximum number of cores per node used in the analysis (defaults to "8")
               -c/--configFile YAML config file for script parameters (defaults to "")
               -wc/--writeConfigFile Write YAML configuration file for script parameters (defaults to "";supply whole path)
               -int/--instanceTag Tag family with instance association in sampleInfo file (defaults to "")
               -rea/--researchEthicalApproval Tag for displaying research candidates in Scout (defaults to "notApproved")
               -sif/--sampleInfoFile YAML file for sample info used in the analysis (defaults to "{outDataDir}/{familyID}/{familyID}_qc_sampleInfo.yaml")
               -dra/--dryRunAll Sets all programs to dry run mode i.e. no sbatch submission (defaults to "0" (=no))
               -tmd/--tempDirectory Set the temporary directory for all programs (defaults to "/scratch/SLURM_JOB_ID";supply whole path)
               -jul/--javaUseLargePages Use large page memory. (-XX,hence option considered not stable and are subject to change without notice, but can be consiered when faced with Java Runtime Environment Memory issues)
               -pve/--pythonVirtualEnvironment Pyhton virtualenvironment (defaults to "")
               -l/--logFile Mip log file (defaults to "{outDataDir}/{familyID}/mip_log/{timestamp}/{scriptname}_{timestamp}.log")
               -h/--help Display this help message    
               -v/--version Display version of MIP            

               
               ####Programs
               -pGZ/--pGZipFastq GZip fastq files (defaults to "1" (=yes))
	       -pFqC/--pFastQC Sequence quality analysis using FastQC (defaults to "1" (=yes))

               ##Mosaik
	       -pMoB/--pMosaikBuild  Convert reads to Mosaik format using MosaikBuild (defaults to "1" (=yes))
                -mobmfl/--mosaikBuildMedianFragLength Flag for setting the mean fragment length, mfl, (defaults to (=375) bp)
	       -pMoA/--pMosaikAlign Align reads using MosaikAlign (defaults to "1" (=yes))
                 -moaref/--mosaikAlignReference MosaikAlign reference (defaults to "{humanGenomeReference}")
                 -moaape/--mosaikAlignNeuralNetworkPeFile MosaikAlign Neural Network PE File (defaults to "2.1.78.pe.ann")
                 -moaase/--mosaikAlignNeuralNetworkSeFile MosaikAlign Neural Network SE File (defaults to "2.1.78.se.ann")
                 -mojdb/--mosaikJumpDbStub MosaikJump stub (defaults to "{humanGenomeReference}")
               
               ##BWA
               -pMem/--pBwaMem Align reads using BWA Mem (defaults to "0" (=no))
                 -memrdb/--bwaMemRapidDb Selection of relevant regions post alignment (defaults to "")
               -pAln/--pBwaAln Index reads using BWA Aln (defaults to "0" (=no))
                 -alnq/--bwaAlnQualityTrimming BWA Aln quality threshold for read trimming (defaults to "20")
               -pSap/--pBwaSampe Align reads using BWA Sampe (defaults to "0" (=no))
               
               ##PicardTools
               -ptp/--picardToolsPath Path to PicardTools. Mandatory for use of PicardTools (defaults to "")
               -pPtS/--pPicardToolsSortSam Sort & index aligned reads using PicardTools SortSam & index (defaults to "1" (=yes))
               -pPtM/--pPicardToolsMergeSamFiles Merge (BAM file(s) ) using PicardTools MergeSamFiles (defaults to "1" (=yes))
               -pPtMR/--pPicardToolsMergeRapidReads Merge Read batch processed (BAM file(s)) using PicardTools MergeSamFiles (Only relevant in rapid mode;defaults to "0" (=no))
                 -ptmp/--picardToolsMergeSamFilesPrevious PicardTools MergeSamFiles on merged current files and previous BAM-file(s) (supply whole path and name, name must contain sample id, and lanes_Xn info)
               -pPtMD/--pPicardToolsMarkduplicates Markduplicates using PicardTools MarkDuplicates (defaults to "1" (=yes))
               
               ##Coverage Calculations
               -pChS/--pChanjoSexCheck Predicts gender from sex chromosome coverage (defaults to "1")
               -pChB/--pChanjoBuild Chanjo build central SQLite database file (defaults to "1" (=yes))
                 -chbdb/--chanjoBuildDb  Reference database (defaults to "CCDS.current.txt")
               -pChA/--pChanjoAnnotate Chanjo coverage analysis (defaults to "1" (=yes))
                 -chacut/--chanjoAnnotateCutoff Read depth cutoff (defaults to "10")
               -pChI/--pChanjoImport Chanjo import to collect sample info to family Db  (defaults to "0" (=no))
               -pGcB/--pGenomeCoverageBED Genome coverage calculation using genomeCoverageBED (defaults to "1" (=yes))
                -gcbcov/--GenomeCoverageBEDMaxCoverage Max coverage depth when using '-pGenomeCoverageBED' (defaults to "30")
               -pPtCMM/--pPicardToolsCollectMultipleMetrics Metrics calculation using PicardTools collectMultipleMetrics (defaults to "1" (=yes))
               -pPtCHS/--pPicardToolsCalculateHSMetrics Capture calculation using PicardTools CalculateHSmetrics (defaults to "1" (=yes))
                 -ptchsetl/--exomeTargetBedInfileLists Prepared target BED file for PicardTools CalculateHSMetrics (defaults to "". File ending should be ".infile_list") 
                 -ptchsetpl/--exomeTargetPaddedBedInfileLists Prepared padded target BED file for PicardTools CalculateHSMetrics (defaults to "". File ending should be ".padXXX.infile_list")
               -pRcP/--pRCovPlots Plots of genome coverage using rCovPlots (defaults to "1" (=yes))
               
               ##GATK              
               -gtp/--genomeAnalysisToolKitPath  Path to GATK. Mandatory for use of GATK (defaults to "")
               -gbdv/--GATKBundleDownLoadVersion  GATK FTP bundle download version.(defaults to "2.8")
               -gtpl/--GATKTargetPaddedBedIntervalLists Target BED file interval for GATK (defaults to "". File ending should be ".padXXX.interval_list")
               -gdco/--GATKDownSampleToCoverage Coverage to downsample to at any given locus (defaults to "1000")
               -pGrA/--pGATKRealigner Realignments of reads using GATK realign (defaults to "1" (=yes))
                 -graks1/--GATKReAlignerINDELKnownSet1 GATK ReAlignerTargetCreator/IndelRealigner known INDEL set 1 (defaults to "1000G_phase1.indels.b37.vcf")
                 -graks2/--GATKReAlignerINDELKnownSet2 GATK ReAlignerTargetCreator/IndelRealigner known INDEL set 2 (defaults to "Mills_and_1000G_gold_standard.indels.b37.vcf")
               -pGbR/--pGATKBaseRecalibration Recalibration of bases using GATK BaseRecalibrator/PrintReads (defaults to "1" (=yes))
                 -gbrkse/--GATKBaseReCalibrationSNPKnownSet GATK BaseReCalinbration known SNP set (defaults to "dbsnp_138.b37.vcf")                
               -pGhC/--pGATKHaploTypeCaller Variant discovery using GATK HaplotypeCaller (defaults to "1" (=yes))
                 -ghckse/--GATKHaploTypeCallerSNPKnownSet GATK HaplotypeCaller dbSNP set for annotating ID columns (defaults to "dbsnp_138.b37.vcf")
               -pGgT/--pGATKGenoTypeGVCFs Merge gVCF records using GATK GenotypeGVCFs (defaults to "1" (=yes))
                 -ggtgrl/--GATKGenoTypeGVCFsRefGVCF GATK GenoTypeGVCFs gVCF reference infile list for joint genotyping (defaults to "")
               -pGvR/--pGATKVariantRecalibration Variant recalibration using GATK VariantRecalibrator/ApplyRecalibration (defaults to "1" (=yes))
                 -gvrtsh/--GATKVariantReCalibrationTrainingSetHapMap GATK VariantRecalibrator HapMap training set (defaults to "hapmap_3.3.b37.sites.vcf")
                 -gvrtss/--GATKVariantReCalibrationTrainingSetDbSNP GATK VariantRecalibrator dbSNP training set (defaults to "dbsnp_138.b37.vcf")
                 -gvrtsg/--GATKVariantReCalibrationTrainingSet1000GSNP GATK VariantRecalibrator 1000G high confidence SNP training set (defaults to "1000G_phase1.snps.high_confidence.b37.vcf")
                 -gvrtso/--GATKVariantReCalibrationTrainingSet1000GOmni GATK VariantRecalibrator 1000G_omni training set (defaults to "1000G_omni2.5.b37.sites.vcf")
                 -gvrtsm/--GATKVariantReCalibrationTrainingSetMills GATK VariantRecalibrator Mills training set (defaults to "Mills_and_1000G_gold_standard.indels.b37.vcf")
                 -gvrtsf/--GATKVariantReCalibrationTSFilterLevel The truth sensitivity level at which to start filtering used in GATK VariantRecalibrator (defaults to "99.9")
                 -gvrevf/--GATKVariantReCalibrationexcludeNonVariantsFile Produce a vcf containing non-variant loci alongside the vcf only containing non-variant loci after GATK VariantRecalibrator (defaults to "false")
               -pGpT/--pGATKPhaseByTransmission Computes the most likely genotype and phases calls were unamibigous using GATK PhaseByTransmission (defaults to "1" (=yes))
               -pGrP/--pGATKReadBackedPhasing Performs physical phasing of SNP calls, based on sequencing reads using GATK ReadBackedPhasing (defaults to "1" (=yes))
                 -grpqth/--GATKReadBackedPhasingPhaseQualityThreshold The minimum phasing quality score required to output phasing (defaults to "20")
               -pGvEA/--pGATKVariantEvalAll Variant evaluation using GATK VariantEval for all variants  (defaults to "1" (=yes))
               -pGvEE/--pGATKVariantEvalExome Variant evaluation using GATK VariantEval for exonic variants  (defaults to "1" (=yes))
                 -gveedbs/--GATKVariantEvalDbSNP DbSNP file used in GATK VariantEval (defaults to "dbsnp_138.b37.excluding_sites_after_129.vcf")
                 -gveedbg/--GATKVariantEvalGold Gold Indel file used in GATK VariantEval (defaults to "Mills_and_1000G_gold_standard.indels.b37.vcf")
               
               ##ANNOTATION
               -pVeP/--pVariantEffectPredictor Annotate variants using VEP (defaults to "1" (=yes))
                 -vepp/--vepDirectoryPath Path to VEP script directory (defaults to ""; supply whole path)
                 -vepc/vepDirectoryCache Specify the cache directory to use (supply whole path, defaults to "") 
                 -vepf/--vepFeatures VEP features (defaults to ("refseq","hgvs","symbol","numbers","sift","polyphen","humdiv"); comma sep)
               -pVcP/--pVCFParser Parse variants using vcfParser.pl (defaults to "1" (=yes))
                 -vcpvt/--vcfParserVepTranscripts Parse VEP transcript specific entries (defaults to "0" (=no))
                 -vcprff/--vcfParserRangeFeatureFile Range annotations file (defaults to ""; tab-sep)
                 -vcprfa/--vcfParserRangeFeatureAnnotationColumns Range annotations feature columns (defaults to ""; comma sep)
                 -vcpsf/--vcfParserSelectFile File containging list of genes to analyse seperately (defaults to "";tab-sep file and HGNC Symbol required)
                 -vcpsfm/--vcfParserSelectFileMatchingColumn Position of HGNC Symbol column in SelectFile (defaults to "")
                 -vcpsfa/--vcfParserSelectFeatureAnnotationColumns Feature columns to use in annotation (defaults to ""; comma sep)
               -pSnE/--pSnpEff Variant annotation using snpEFF (defaults to "1" (=yes))
                 -snep/--snpEffPath Path to snpEff. Mandatory for use of snpEff (defaults to "")
                 -snesaf2/--snpSiftAnnotationFiles Annotation files to use with snpSift (comma sep)
                 -snesdbnsfp/--snpSiftDbNSFPFile DbNSFP File (defaults to "dbNSFP2.6.txt.gz")
                 -snesdbnsfpa/--snpSiftDbNSFPAnnotations DbNSFP annotations to use with snpSift (defaults to ("SIFT_pred","Polyphen2_HDIV_pred","Polyphen2_HVAR_pred","LRT_pred","MutationTaster_pred","GERP++_NR","GERP++_RS","phastCons100way_vertebrate","1000Gp1_AF","ESP6500_AA_AF"); comma sep)
               -pAnV/--pAnnovar Annotate variants using Annovar (defaults to "1" (=yes))
                 -anvp/--annovarPath  Path to Annovar script directory (supply whole path, defaults to "". NOTE: Assumes that the annovar db files are located in annovar/humandb)
                 -anvgbv/--annovarGenomeBuildVersion Annovar genome build version (defaults to "hg19")
                 -anvtn/--annovarTableNames Annovar table names (defaults to ("refGene","mce46way","gerp++elem","segdup","tfbs","mirna","snp137NonFlagged","1000g2012apr_all","esp6500si_all","ljb2_sift","ljb2_pp2hdiv","ljb2_pp2hvar","ljb2_mt","ljb2_lrt","ljb2_gerp++","ljb2_phylop"); comma sep)
                 -anvstn/--annovarSupportedTableNames Print Annovar MIP supported table names (defaults 0 (=no))
                 -anvarmafth/--annovarMAFThreshold Sets the minor allele frequency threshold in annovar (defaults to "0")

               ##RankVariants
               -pRaV/--pRankVariants Ranking of annotated variants (defaults to "1" (=yes))
                 -ravgf/--geneFile Defines genes to use when calculating compounds (defaults to "hg19_refGene.txt")
                 -ravcs/--caddWGSSNVs Annotate whole genome sequencing CADD score (defaults to "0" (=no))
                 -ravcsf/--caddWGSSNVsFile Whole genome sequencing CADD score file (defaults to "whole_genome_SNVs.tsv.gz")
                 -ravc1kg/--cadd1000Genomes 1000 Genome cadd score file (defaults to "0" (=no))
                 -ravc1kgf/--cadd1000GenomesFile 1000 Genome cadd score file (defaults to "1000G.tsv.gz")
                 -ravwg/--wholeGene Allow compound pairs in intronic regions (defaults to "1" (=yes))
                 -ravrm/--rankModelFile Rank model config file (defaults to "")
               
               ##Utility
               -pScK/--pSampleCheck QC for samples gender and relationship (defaults to "1" (=yes) )
               -pQcC/--pQCCollect Collect QC metrics from programs processed (defaults to "1" (=yes) )
                 -qccsi/--QCCollectSampleInfoFile SampleInfo File containing info on what to parse from this analysis run (defaults to "{outDataDir}/{familyID}/{familyID}_qc_sampleInfo.yaml")
                 -qccref/--QCCollectRegExpFile Regular expression file containing the regular expression to be used for each program (defaults to "")
               -pReM/--pRemoveRedundantFiles Generating sbatch script for deletion of redundant files (defaults to "1" (=yes);Note: Must be submitted manually)
               -pArS/--pAnalysisRunStatus Sets the analysis run status flag to finished in sampleInfoFile (defaults to "1" (=yes))
	   };
}


####Script parameters

my %parameter;  #Holds all parameters for MIP
my %scriptParameter;  #Holds all active parameters after the value has been set

$scriptParameter{'MIP'} = 1;  #Enable/activate MIP

my $logger;  #Will hold the logger object for the MIP log
my @orderParameters;  #To add/write parameters in the correct order
my @broadcasts;  #Holds all set parameters info after AddToScriptParameter

##Add dateTimestamp for later use in log and qcmetrics yaml file
my $dateTimeStamp = (`date +%Y%m%d_%Hh%Mm`);  #Catches current date, time and script name
chomp($dateTimeStamp);  #Remove \n
my ($base, $script) = (`date +%Y%m%d`,`basename $0`);  #Catches current date and script name
chomp($base, $script);  #Remove \n;

####Set program parameters

###Project specific
##DefineParameters
##parameterName, parameterType, parameterDefault, AssociatedProgram, Check directory/file existence, parameterChain, programCheck)
##DefineParametersPath
##parameterName, parameterDefault, AssociatedProgram, Check directory/file existence, File Autovivication)

&DefineParametersPath(\%parameter, \@orderParameters, "familyID", "nodefault", "MIP", 0);

&DefineParametersPath(\%parameter, \@orderParameters, "outDataDir", "nodefault", "MIP", 0);

&DefineParametersPath(\%parameter, \@orderParameters, "logFile", "NotsetYet", "MIP", "file", "noAutoBuild");

&DefineParameters(\%parameter, \@orderParameters, "projectID", "MIP", "nodefault", "MIP");

&DefineParameters(\%parameter, \@orderParameters, "email", "MIP", 0, "MIP");

&DefineParameters(\%parameter, \@orderParameters, "emailType", "MIP", "F", "MIP");

&DefineParameters(\%parameter, \@orderParameters, "maximumCores", "MIP", 16, "MIP");

&DefineParametersPath(\%parameter, \@orderParameters, "configFile", 0, "MIP", "file");

&DefineParameters(\%parameter, \@orderParameters, "analysisType", "MIP", "exomes", "MIP");

&DefineParametersPath(\%parameter, \@orderParameters, "outScriptDir", "nodefault", "MIP", 0);

&DefineParametersPath(\%parameter, \@orderParameters, "writeConfigFile", 0, "MIP", 0);

&DefineParametersPath(\%parameter, \@orderParameters, "pedigreeFile", "nodefault", "MIP", "file", "noAutoBuild");

&DefineParametersPath(\%parameter, \@orderParameters, "sampleInfoFile", "NotsetYet", "MIP", "file", "noAutoBuild");

&DefineParameters(\%parameter, \@orderParameters, "instanceTag", "MIP", "Unknown", "MIP");

&DefineParameters(\%parameter, \@orderParameters, "researchEthicalApproval", "MIP", "notApproved", "MIP");

&DefineParametersPath(\%parameter, \@orderParameters, "inScriptDir", "nodefault", "MIP", "directory");

&DefineParametersPath(\%parameter, \@orderParameters, "referencesDir", "nodefault", "MIP", "directory");

&DefineParameters(\%parameter, \@orderParameters, "dryRunAll", "MIP", 0, "MIP");

&DefineParametersPath(\%parameter, \@orderParameters, "tempDirectory", "/scratch/".'$SLURM_JOB_ID', "MIP", 0);

###Programs

##GZip
&DefineParameters(\%parameter, \@orderParameters, "pGZipFastq", "program", 1, "MIP", "nofileEnding", "MAIN", "gzip");


##FastQC
&DefineParameters(\%parameter, \@orderParameters, "pFastQC", "program", 1, "MIP", "nofileEnding", "RawSeqQC", "fastqc");


##Mosaik
&DefineParameters(\%parameter, \@orderParameters, "pMosaikBuild", "program", 1, "MIP", "nofileEnding", "MAIN", "MosaikBuild");

&DefineParameters(\%parameter, \@orderParameters, "mosaikBuildMedianFragLength", "program", 375, "pMosaikBuild");

&DefineParameters(\%parameter, \@orderParameters, "pMosaikAlign", "program", 1, "MIP", "nofileEnding", "MAIN", "MosaikAligner");

&DefineParametersPath(\%parameter, \@orderParameters, "mosaikAlignReference", "notSetYet", "pMosaikAlign", "file", "yesAutoBuild");

&DefineParametersPath(\%parameter, \@orderParameters, "mosaikAlignNeuralNetworkPeFile", "2.1.78.pe.ann", "pMosaikAlign", "file", "yesAutoBuild");

&DefineParametersPath(\%parameter, \@orderParameters, "mosaikAlignNeuralNetworkSeFile", "2.1.78.se.ann", "pMosaikAlign", "file", "yesAutoBuild");

&DefineParametersPath(\%parameter, \@orderParameters, "mosaikJumpDbStub", "notSetYet", "pMosaikAlign", "file", "yesAutoBuild");

my @mosaikJumpDbStubFileEndings = ("_keys.jmp", "_meta.jmp", "_positions.jmp");


##BWA
&DefineParameters(\%parameter, \@orderParameters, "pBwaMem", "program", 0, "MIP", "nofileEnding", "MAIN", "bwa");

&DefineParametersPath(\%parameter, \@orderParameters, "bwaMemRapidDb", "nodefault", "pBwaMem", "file", "noAutoBuild");

&DefineParameters(\%parameter, \@orderParameters, "pBwaAln", "program", 0, "MIP", "nofileEnding", "MAIN", "bwa");

&DefineParameters(\%parameter, \@orderParameters, "bwaAlnQualityTrimming", "program", 20, "pBwaAln");

&DefineParameters(\%parameter, \@orderParameters, "pBwaSampe", "program", 0, "MIP", "nofileEnding", "MAIN", "bwa");

&DefineParametersPath(\%parameter, \@orderParameters, "bwaBuildReference", "notSetYet", "pBwaMem,pBwaAln,pBwaSampe", "file", "yesAutoBuild");

my @bwaBuildReferenceFileEndings = (".amb", ".ann", ".bwt", ".pac", ".sa");


##Choosen MIP Aligner
&DefineParameters(\%parameter, \@orderParameters, "aligner", "MIP", "mosaik", "MIP");

##PicardTools
&DefineParameters(\%parameter, \@orderParameters, "pPicardToolsSortSam", "program", 1, "MIP", "_sorted", "MAIN");

&DefineParameters(\%parameter, \@orderParameters, "pPicardToolsMergeRapidReads", "program", 0, "MIP", "_sorted", "MAIN");  #Rapid mode special case

&DefineParameters(\%parameter, \@orderParameters, "pPicardToolsMergeSamFiles", "program", 1, "MIP", "_merged", "MAIN");

&DefineParameters(\%parameter, \@orderParameters, "pPicardToolsMarkduplicates", "program", 1, "MIP", "_pmd", "MAIN");


##Coverage
&DefineParameters(\%parameter, \@orderParameters, "pChanjoSexCheck", "program", 1, "MIP",".sexcheck", "CoverageReport_Gender");

&DefineParameters(\%parameter, \@orderParameters, "pChanjoBuild", "program", 1, "MIP", "nofileEnding", "CoverageReport");

&DefineParametersPath(\%parameter, \@orderParameters, "chanjoBuildDb", "CCDS.current.txt", "pChanjoBuild", "file", "yesAutoDownLoad");

&DefineParameters(\%parameter, \@orderParameters, "pChanjoAnnotate", "program", 1, "MIP","_coverage", "CoverageReport");

&DefineParameters(\%parameter, \@orderParameters, "chanjoAnnotateCutoff", "program", 10, "pChanjoAnnotate");

&DefineParameters(\%parameter, \@orderParameters, "pChanjoImport", "program", 0, "MIP", "nofileEnding", "CoverageReport");

&DefineParameters(\%parameter, \@orderParameters, "pGenomeCoverageBED", "program", 1, "MIP", "_genomeCoverageBed", "CoverageQC_GcovBed", "bedtools");

&DefineParameters(\%parameter, \@orderParameters, "pPicardToolsCollectMultipleMetrics", "program", 1, "MIP", "nofileEnding", "CoverageQC_PTCMM");

&DefineParameters(\%parameter, \@orderParameters, "pPicardToolsCalculateHSMetrics", "program", 1, "MIP", "nofileEnding", "CoverageQC_PTCHSM");

&DefineParameters(\%parameter, \@orderParameters, "GenomeCoverageBEDMaxCoverage", "program", 30, "pGenomeCoverageBED");

&DefineParameters(\%parameter, \@orderParameters, "pRCovPlots", "program", 0, "MIP", "nofileEnding", "CoverageQC_RCOVP");

&DefineParametersPath(\%parameter, \@orderParameters, "picardToolsPath", "nodefault", "pBwaMem,pPicardToolsMergeSamFiles,pPicardToolsMarkduplicates,pPicardToolsCalculateHSMetrics,pPicardToolsCollectMultipleMetrics,pGATKHaploTypeCaller,pGATKVariantRecalibration", "directory");  #pGATKHaploTypeCaller,pGATKVariantRecalibration since these jars can use merged interval_list files, which are created in MIP with picardTools

##Target definition files
my (@exomeTargetBedInfileLists, @exomeTargetPaddedBedInfileLists);  #Arrays for target bed infile lists


##GATK
&DefineParameters(\%parameter, \@orderParameters, "pGATKRealigner", "program", 1, "MIP", "_rreal", "MAIN");

&DefineParametersPath(\%parameter, \@orderParameters, "GATKReAlignerINDELKnownSet1", "1000G_phase1.indels.b37.vcf", "pGATKRealigner", "file", "yesAutoDownLoad");

&DefineParametersPath(\%parameter, \@orderParameters, "GATKReAlignerINDELKnownSet2", "Mills_and_1000G_gold_standard.indels.b37.vcf", "pGATKRealigner", "file", "yesAutoDownLoad");


&DefineParameters(\%parameter, \@orderParameters, "pGATKBaseRecalibration", "program", 1, "MIP", "_brecal", "MAIN");

&DefineParametersPath(\%parameter, \@orderParameters, "GATKBaseReCalibrationSNPKnownSet", "dbsnp_138.b37.vcf", "pGATKBaseRecalibration", "file", "yesAutoDownLoad");


&DefineParameters(\%parameter, \@orderParameters, "pGATKHaploTypeCaller", "program", 1, "MIP", "_gvcf", "MAIN");

&DefineParametersPath(\%parameter, \@orderParameters, "GATKHaploTypeCallerSNPKnownSet", "dbsnp_138.b37.vcf", "pGATKHaploTypeCaller", "file", "yesAutoDownLoad");


&DefineParameters(\%parameter, \@orderParameters, "pGATKGenoTypeGVCFs", "program", 1, "MIP", "_", "MAIN");

&DefineParametersPath(\%parameter, \@orderParameters, "GATKGenoTypeGVCFsRefGVCF", "nodefault", "pGATKGenoTypeGVCFs", "file", "noAutoBuild");


&DefineParameters(\%parameter, \@orderParameters, "pGATKVariantRecalibration", "program", 1, "MIP", "vrecal_", "MAIN");

&DefineParametersPath(\%parameter, \@orderParameters, "GATKVariantReCalibrationTrainingSetHapMap", "hapmap_3.3.b37.sites.vcf", "pGATKVariantRecalibration", "file", "yesAutoDownLoad");

&DefineParametersPath(\%parameter, \@orderParameters, "GATKVariantReCalibrationTrainingSetDbSNP", "dbsnp_138.b37.vcf", "pGATKVariantRecalibration", "file", "yesAutoDownLoad");

&DefineParametersPath(\%parameter, \@orderParameters, "GATKVariantReCalibrationTrainingSet1000GSNP", "1000G_phase1.snps.high_confidence.b37.vcf", "pGATKVariantRecalibration", "file", "yesAutoDownLoad");

&DefineParametersPath(\%parameter, \@orderParameters, "GATKVariantReCalibrationTrainingSet1000GOmni", "1000G_omni2.5.b37.sites.vcf", "pGATKVariantRecalibration", "file", "yesAutoDownLoad");

&DefineParametersPath(\%parameter, \@orderParameters, "GATKVariantReCalibrationTrainingSetMills", "Mills_and_1000G_gold_standard.indels.b37.vcf", "pGATKVariantRecalibration", "file", "yesAutoDownLoad");

&DefineParameters(\%parameter, \@orderParameters, "GATKVariantReCalibrationTSFilterLevel", "program", 99.9, "pGATKVariantRecalibration");

&DefineParameters(\%parameter, \@orderParameters, "GATKVariantReCalibrationexcludeNonVariantsFile", "program", "false", "pGATKVariantRecalibration");

 
&DefineParameters(\%parameter, \@orderParameters, "pGATKPhaseByTransmission", "program", 1, "MIP", "phtr_", "Phasing");

&DefineParameters(\%parameter, \@orderParameters, "pGATKReadBackedPhasing", "program", 1, "MIP", "phrb_", "Phasing");

&DefineParameters(\%parameter, \@orderParameters, "GATKReadBackedPhasingPhaseQualityThreshold", "program", 20, "pGATKReadBackedPhasing");


&DefineParameters(\%parameter, \@orderParameters, "pGATKVariantEvalAll", "program", 1, "MIP", "nofileEnding", "AllVariantQC");

&DefineParameters(\%parameter, \@orderParameters, "pGATKVariantEvalExome", "program", 1, "MIP", "nofileEnding", "ExomeVariantQC", "bedtools");

&DefineParametersPath(\%parameter, \@orderParameters, "GATKVariantEvalDbSNP", "dbsnp_138.b37.excluding_sites_after_129.vcf", "pGATKVariantEvalAll,pGATKVariantEvalExome", "file", "yesAutoDownLoad");

&DefineParametersPath(\%parameter, \@orderParameters, "GATKVariantEvalGold", "Mills_and_1000G_gold_standard.indels.b37.vcf", "pGATKVariantEvalAll,pGATKVariantEvalExome", "file", "yesAutoDownLoad");

&DefineParameters(\%parameter, \@orderParameters, "GATKDownSampleToCoverage", "program", 1000, "pGATKRealigner,pGATKBaseRecalibration,pGATKHaploTypeCaller");

&DefineParameters(\%parameter, \@orderParameters, "GATKBundleDownLoadVersion", "program", "2.8", "pBwaMem,pBwaAln,pBwaSampe,pGATKRealigner,pGATKBaseRecalibration,pGATKHaploTypeCaller,pGATKGenoTypeGVCFs,pGATKVariantRecalibration,pGATKPhaseByTransmission,pGATKReadBackedPhasing,pGATKVariantEvalAll,pGATKVariantEvalExome,pAnnovar,pAddDepth,pPicardToolsCalculateHSMetrics,pPicardToolsCollectMultipleMetrics");  #Sets the GATK FTP Bundle Download version. Needed for all programs that download the human genome reference

my (@GATKTargetPaddedBedIntervalLists);  #Array for target infile lists used in GATK


##VEP
&DefineParameters(\%parameter, \@orderParameters, "pVariantEffectPredictor", "program", 1, "MIP", "vep_", "MAIN");

&DefineParametersPath(\%parameter, \@orderParameters, "vepDirectoryPath", "nodefault", "pVariantEffectPredictor", "directory");  #Note not projectID specific

&DefineParametersPath(\%parameter, \@orderParameters, "vepDirectoryCache", "nodefault", "pVariantEffectPredictor", "directory");


##VCFParser
&DefineParameters(\%parameter, \@orderParameters, "pVCFParser", "program", 1, "MIP", "parsed_", "MAIN");

&DefineParameters(\%parameter, \@orderParameters, "vcfParserVepTranscripts", "program", 0, "pVCFParser");

&DefineParametersPath(\%parameter, \@orderParameters, "vcfParserRangeFeatureFile", "noUserInfo", "pVCFParser", "file"); 

&DefineParametersPath(\%parameter, \@orderParameters, "vcfParserSelectFile", "noUserInfo", "pVCFParser", "file"); 

&DefineParameters(\%parameter, \@orderParameters, "vcfParserSelectFileMatchingColumn", "program", "nodefault", "pVCFParser");

my $VEPOutputFiles = 1;  #To track if VEPParser was used with a vcfParserSelectFile (=2) or not (=1)


##SnpEFF
&DefineParameters(\%parameter, \@orderParameters, "pSnpEff", "program", 1, "MIP", "snpeff_", "MAIN");

&DefineParametersPath(\%parameter, \@orderParameters, "snpEffPath", "nodefault", "pSnpEff", "directory");

&DefineParametersPath(\%parameter, \@orderParameters, "snpSiftDbNSFPFile", "dbNSFP2.6.txt.gz", "pSnpEff", "file");


##Annovar
&DefineParameters(\%parameter, \@orderParameters, "pAnnovar", "program", 1, "MIP", "annovar_", "MAIN");

&DefineParametersPath(\%parameter, \@orderParameters, "annovarPath", "nodefault", "pAnnovar", "directory");  #Note not projectID specific

&DefineParameters(\%parameter, \@orderParameters, "annovarGenomeBuildVersion", "program", "hg19", "pAnnovar");

&DefineParameters(\%parameter, \@orderParameters, "annovarSupportedTableNames", "program", 0, "pAnnovar");

&DefineParameters(\%parameter, \@orderParameters, "annovarMAFThreshold", "program", 0, "pAnnovar");


##Special case GATKPath since in VEP, SnpEff and Annovar modules use GATK CombineVariants to merge vcfs 
&DefineParametersPath(\%parameter, \@orderParameters, "javaUseLargePages", "no", "pGATKRealigner,pGATKBaseRecalibration,pGATKHaploTypeCaller,pGATKGenoTypeGVCFs,pGATKVariantRecalibration,pGATKPhaseByTransmission,pGATKReadBackedPhasing,pGATKVariantEvalAll,pGATKVariantEvalExome,pVariantEffectPredictor,pSnpEff,pAnnovar");

&DefineParametersPath(\%parameter, \@orderParameters, "genomeAnalysisToolKitPath", "nodefault", "pGATKRealigner,pGATKBaseRecalibration,pGATKHaploTypeCaller,pGATKGenoTypeGVCFs,pGATKVariantRecalibration,pGATKPhaseByTransmission,pGATKReadBackedPhasing,pGATKVariantEvalAll,pGATKVariantEvalExome,pVariantEffectPredictor,pSnpEff", "directory");


##SChecks
&DefineParameters(\%parameter, \@orderParameters, "pSampleCheck", "program", 1, "MIP", "nofileEnding", "IDQC", "vcftools:plink");


##RankVariants
&DefineParameters(\%parameter, \@orderParameters, "pRankVariants", "program", 1, "MIP", "ranked_", "MAIN");

&DefineParametersPath(\%parameter, \@orderParameters, "geneFile", "hg19_refGene.txt", "pRankVariants", "file", "noAutoBuild");

&DefineParameters(\%parameter, \@orderParameters, "caddWGSSNVs", "program", 0, "pRankVariants");

&DefineParametersPath(\%parameter, \@orderParameters, "caddWGSSNVsFile", "whole_genome_SNVs.tsv.gz", "pRankVariants", "file", "noAutoBuild");

&DefineParameters(\%parameter, \@orderParameters, "cadd1000Genomes", "program", 0, "pRankVariants");

&DefineParametersPath(\%parameter, \@orderParameters, "cadd1000GenomesFile", "1000G.tsv.gz", "pRankVariants", "file", "noAutoBuild");

&DefineParameters(\%parameter, \@orderParameters, "wholeGene", "program", 1, "pRankVariants");

&DefineParametersPath(\%parameter, \@orderParameters, "rankModelFile", "noUserInfo", "pRankVariants", "file", "noAutoBuild");

&DefineParametersPath(\%parameter, \@orderParameters, "pythonVirtualEnvironment", "nodefault", "pChanjoBuild,pChanjoAnnotate,pChanjoImport,pRankVariants");


##QcCollect
&DefineParameters(\%parameter, \@orderParameters, "pQCCollect", "program", 1, "MIP", "nofileEnding", "MAIN");

&DefineParameters(\%parameter, \@orderParameters, "QCCollectSampleInfoFile", "program", "notSetYet", "pQCCollect");  #No file check since file is created by MIP later

&DefineParametersPath(\%parameter, \@orderParameters, "QCCollectRegExpFile", "qc_regexp.yaml", "pQCCollect", "file", "noAutoBuild");


##RemoveRedundantFiles
&DefineParameters(\%parameter, \@orderParameters, "pRemoveRedundantFiles", "program", 1, "MIP", "nofileEnding", "MAIN");


##AnalysisRunStatus
&DefineParameters(\%parameter, \@orderParameters, "pAnalysisRunStatus", "program", 1, "MIP", "", "MAIN");


##MIP

##humanGenomeReference
&DefineParametersPath(\%parameter, \@orderParameters, "humanGenomeReference", "Homo_sapiens.GRCh37.d5.fasta", "pBwaMem,pBwaAln,pBwaSampe,pGATKRealigner,pGATKBaseRecalibration,pGATKHaploTypeCaller,pGATKGenoTypeGVCFs,pGATKVariantRecalibration,pGATKPhaseByTransmission,pGATKReadBackedPhasing,pGATKVariantEvalAll,pGATKVariantEvalExome,pAnnovar,pPicardToolsCalculateHSMetrics,pPicardToolsCollectMultipleMetrics", "file", "yesAutoDownLoad");

my @humanGenomeReferenceFileEndings = (".dict", ".fasta.fai");  #Meta files

my ($aligner, $version, $help) = ("nocmdinput");

my (@contigs);  #Holds all contigs, not just chromosomes

my (%infile, %inDirPath, %infilesLaneNoEnding, %lane, %infilesBothStrandsNoEnding, %jobID, %sampleInfo, %fileInfo); 


####Staging/Sanity Check Area 

##Capture kit aliases supported from pedigree file.
my %supportedCaptureKit = (
    'Nimblegen_SeqCapEZExome.V2' => "Nimblegen_SeqCapEZExome.V2.GenomeReferenceSourceVersion_targets.bed",
    'Nimblegen_SeqCapEZExome.V3' => "Nimblegen_SeqCapEZExome.V3.GenomeReferenceSourceVersion_targets.bed",
    'Agilent_SureSelect.V2' => "Agilent_SureSelect.V2.GenomeReferenceSourceVersion_targets.bed",
    'Agilent_SureSelect.V3' => "Agilent_SureSelect.V3.GenomeReferenceSourceVersion_targets.bed",
    'Agilent_SureSelect.V4' => "Agilent_SureSelect.V4.GenomeReferenceSourceVersion_targets.bed",
    'Agilent_SureSelect.V5' => "Agilent_SureSelect.V5.GenomeReferenceSourceVersion_targets.bed",
    'Agilent_SureSelectCRE.V1' => "Agilent_SureSelectCRE.V1.GenomeReferenceSourceVersion_targets.bed",
    'Latest' => "Agilent_SureSelect.V5.GenomeReferenceSourceVersion_targets.bed",
    );

my %supportedCosmidReference;  #References supported as downloads from Cosmid. Hash is populated after user options are processed

my %referenceFileEndings = (
    'mosaikAlignReference' => ".dat",
    'mosaikJumpDbStub' => "_jdb_15",
    'bwaBuildReference' => "",
    'exomeTargetBedInfileLists' => ".infile_list",
    'exomeTargetPaddedBedInfileLists' => ".pad100.infile_list",
    'GATKTargetPaddedBedIntervalLists' => ".pad100.interval_list",
    );

##Set supported annovar table name filtering options
my @annovarSupportedTableNames = ("refGene", "knownGene", "ensGene", "mce46way", "gerp++elem", "segdup", "gwascatalog", "tfbs", "mirna", "snp137", "snp135", "snp132", "snp131", "snp130", "snp129", "snp137NonFlagged", "snp135NonFlagged", "snp132NonFlagged", "snp131NonFlagged", "snp130NonFlagged", "1000g2012apr_all", "1000g2012apr_amr", "1000g2012apr_eur", "1000g2012apr_asn", "1000g2012apr_afr", "1000g2012feb_all", "esp6500si_all", "esp6500_all", "esp6500_aa", "esp6500_ea", "esp5400_all", "esp5400_aa", "esp5400_ea","clinvar_20131105", "ljb2_sift", "ljb2_pp2hdiv", "ljb2_pp2hvar", "ljb2_mt", "ljb2_ma", "ljb2_fathmm", "ljb2_siphy", "ljb2_lrt", "ljb_all", "ljb2_gerp++", "ljb2_phylop", "caddgt20", "caddgt10");  #Used to print list of supported table names

my %annovarTable;  #Holds annovar tables and features

###User Options
GetOptions('ifd|inFilesDirs:s'  => \@{$parameter{'inFilesDirs'}{'value'}},  #Comma separated list
	   'isd|inScriptDir:s'  => \$parameter{'inScriptDir'}{'value'},  #Directory for custom scripts required by the pipeline
	   'rd|referencesDir:s'  => \$parameter{'referencesDir'}{'value'},  #directory containing references
	   'p|projectID:s'  => \$parameter{'projectID'}{'value'},
	   's|sampleIDs:s'  => \@{$parameter{'sampleIDs'}{'value'}},  #Comma separated list, one below outDataDir
	   'em|email:s'  => \$parameter{'email'}{'value'},  #Email adress
	   'emt|emailType:s'  => \$parameter{'emailType'}{'value'},  #Email type 
	   'odd|outDataDir:s'  => \$parameter{'outDataDir'}{'value'},  #One dir above sample id, must supply whole path i.e. /proj/...
	   'osd|outScriptDir:s'  => \$parameter{'outScriptDir'}{'value'},   #One dir above sample id, must supply whole path i.e. /proj/...
	   'f|familyID:s' => \$parameter{'familyID'}{'value'},  #Family group ID (Merged to same vcf file after GATK Base Recalibration)
	   'ped|pedigreeFile:s' => \$parameter{'pedigreeFile'}{'value'},  #Pedigree file
	   'hgr|humanGenomeReference:s' => \$parameter{'humanGenomeReference'}{'value'},  #Human genome reference
	   'al|aligner:s' => \$parameter{'aligner'}{'value'},  #determining which aligner was used previously (if not specified)
	   'at|analysisType:s' => \$parameter{'analysisType'}{'value'},  #Type of analysis
	   'mc|maximumCores:n' => \$parameter{'maximumCores'}{'value'},  #Per node
	   'c|configFile:s' => \$parameter{'configFile'}{'value'},
	   'wc|writeConfigFile:s' => \$parameter{'writeConfigFile'}{'value'},
	   'sif|sampleInfoFile:s' => \$parameter{'sampleInfoFile'}{'value'},  #Write all info on samples and run to YAML file
	   'int|instanceTag:s' => \$parameter{'instanceTag'}{'value'},
	   'rea|researchEthicalApproval:s' => \$parameter{'researchEthicalApproval'}{'value'},
	   'dra|dryRunAll:n' => \$parameter{'dryRunAll'}{'value'},
	   'tmd|tempDirectory:s' => \$parameter{'tempDirectory'}{'value'},
	   'pve|pythonVirtualEnvironment:s' => \$parameter{'pythonVirtualEnvironment'}{'value'},
	   'jul|javaUseLargePages:s' => \$parameter{'javaUseLargePages'}{'value'},
	   'l|logFile:s' => \$parameter{'logFile'}{'value'},
	   'h|help' => \$help,  #Display help text
	   'v|version' => \$version,  #Display version number
	   'pGZ|pGZipFastq:n' => \$parameter{'pGZipFastq'}{'value'},
	   'pFqC|pFastQC:n' => \$parameter{'pFastQC'}{'value'},
	   'pMoB|pMosaikBuild:n' => \$parameter{'pMosaikBuild'}{'value'},
	   'mobmfl|mosaikBuildMedianFragLength:n' => \$parameter{'mosaikBuildMedianFragLength'}{'value'},  #for fragment length estimation and local search
	   'pMoA|pMosaikAlign:n' => \$parameter{'pMosaikAlign'}{'value'},
	   'moaref|mosaikAlignReference:s' => \$parameter{'mosaikAlignReference'}{'value'},  #MosaikAlign reference file assumes existance of jump database files in same dir
	   'moaape|mosaikAlignNeuralNetworkPeFile:s' => \$parameter{'mosaikAlignNeuralNetworkPeFile'}{'value'},
	   'moaase|mosaikAlignNeuralNetworkSeFile:s' => \$parameter{'mosaikAlignNeuralNetworkSeFile'}{'value'}, 
	   'mojdb|mosaikJumpDbStub:s' => \$parameter{'mosaikJumpDbStub'}{'value'},  #Stub for MosaikJump database
	   'pMem|pBwaMem:n' => \$parameter{'pBwaMem'}{'value'},
	   'memrdb|bwaMemRapidDb:s' => \$parameter{'bwaMemRapidDb'}{'value'},
	   'pAln|pBwaAln:n' => \$parameter{'pBwaAln'}{'value'},
	   'alnq|bwaAlnQualityTrimming:n' => \$parameter{'bwaAlnQualityTrimming'}{'value'},  #BWA aln quality threshold for read trimming down to 35bp
	   'pSap|pBwaSampe:n' => \$parameter{'pBwaSampe'}{'value'},
	   'pPtS|pPicardToolsSortSam:n' => \$parameter{'pPicardToolsSortSam'}{'value'},
	   'pPtM|pPicardToolsMergeSamFiles:n' => \$parameter{'pPicardToolsMergeSamFiles'}{'value'},  #PicardTools mergeSamFiles
	   'pPtMR|pPicardToolsMergeRapidReads:n' => \$parameter{'pPicardToolsMergeRapidReads'}{'value'},  #PicardTools mergeSamFiles - rapid mode
	   'ptmp|picardToolsMergeSamFilesPrevious:s' => \@{$parameter{'picardToolsMergeSamFilesPrevious'}{'value'}},  #Comma separated list
	   'pPtMD|pPicardToolsMarkduplicates:s' => \$parameter{'pPicardToolsMarkduplicates'}{'value'},  #PicardTools MarkDuplicates
	   'ptp|picardToolsPath:s' => \$parameter{'picardToolsPath'}{'value'},  #Path to picardtools
	   'pChS|pChanjoSexCheck:n' => \$parameter{'pChanjoSexCheck'}{'value'},   #Chanjo coverage analysis on sex chromosomes
	   'pChB|pChanjoBuild:n' => \$parameter{'pChanjoBuild'}{'value'},   #Build central SQLiteDatabase
	   'chbdb|chanjoBuildDb:s' => \$parameter{'chanjoBuildDb'}{'value'},  #Chanjo reference database
	   'pChA|pChanjoAnnotate:n' => \$parameter{'pChanjoAnnotate'}{'value'},   #Chanjo coverage analysis
	   'chacut|chanjoAnnotateCutoff:n' => \$parameter{'chanjoAnnotateCutoff'}{'value'},   # Cutoff used for completeness
	   'pChI|pChanjoImport:n' => \$parameter{'pChanjoImport'}{'value'},   #Build family SQLiteDatabase
	   'pGcB|pGenomeCoverageBED:n' => \$parameter{'pGenomeCoverageBED'}{'value'},
	   'xcov|GenomeCoverageBEDMaxCoverage:n' => \$parameter{'GenomeCoverageBEDMaxCoverage'}{'value'},  #Sets max depth to calculate coverage
	   'pPtCMM|pPicardToolsCollectMultipleMetrics:n' => \$parameter{'pPicardToolsCollectMultipleMetrics'}{'value'},
	   'pPtCHS|pPicardToolsCalculateHSMetrics:n' => \$parameter{'pPicardToolsCalculateHSMetrics'}{'value'},
	   'ptchsetl|exomeTargetBedInfileLists:s' => \@exomeTargetBedInfileLists,  #Comma separated list of target file for CalculateHsMetrics
	   'ptchsetpl|exomeTargetPaddedBedInfileLists:s' => \@exomeTargetPaddedBedInfileLists,  #Comma separated list of padded target file for CalculateHsMetrics
	   'pRcP|pRCovPlots:n' => \$parameter{'pRCovPlots'}{'value'},
	   'gtp|genomeAnalysisToolKitPath:s' => \$parameter{'genomeAnalysisToolKitPath'}{'value'},  #GATK whole path
	   'gbdv|GATKBundleDownLoadVersion:s' => \$parameter{'GATKBundleDownLoadVersion'}{'value'},  #Sets the GATK FTP Bundle Download version
	   'gtpl|GATKTargetPaddedBedIntervalLists:s' => \@GATKTargetPaddedBedIntervalLists,  #Comma separated list of padded target file set to be used in GATK
	   'gdco|GATKDownSampleToCoverage:n' => \$parameter{'GATKDownSampleToCoverage'}{'value'},  #GATK downsample to coverage
	   'pGrA|pGATKRealigner:n' => \$parameter{'pGATKRealigner'}{'value'},  #GATK ReAlignerTargetCreator/IndelRealigner
	   'graks1|GATKReAlignerINDELKnownSet1:s' => \$parameter{'GATKReAlignerINDELKnownSet1'}{'value'},  #Known INDEL set to be used in GATK ReAlignerTargetCreator/IndelRealigner
	   'graks2|GATKReAlignerINDELKnownSet2:s' => \$parameter{'GATKReAlignerINDELKnownSet2'}{'value'},  #Known INDEL set to be used in GATK ReAlignerTargetCreator/IndelRealigner
	   'pGbR|pGATKBaseRecalibration:n' => \$parameter{'pGATKBaseRecalibration'}{'value'},  #GATK BaseRecalibrator/PrintReads
	   'gbrkse|GATKBaseReCalibrationSNPKnownSet:s' => \$parameter{'GATKBaseReCalibrationSNPKnownSet'}{'value'},  #Known SNP set to be used in GATK BaseRecalibrator/PrintReads
	   'pGhC|pGATKHaploTypeCaller:n' => \$parameter{'pGATKHaploTypeCaller'}{'value'},  #GATK Haplotypecaller
	   'ghckse|GATKHaploTypeCallerSNPKnownSet:s' => \$parameter{'GATKHaploTypeCallerSNPKnownSet'}{'value'},  #Known SNP set to be used in GATK HaplotypeCaller
	   'pGgT|pGATKGenoTypeGVCFs:n' => \$parameter{'pGATKGenoTypeGVCFs'}{'value'},  #Merge gVCF records using GATK GenotypeGVCFs
	   'ggtgrl|GATKGenoTypeGVCFsRefGVCF:s' => \$parameter{'GATKGenoTypeGVCFsRefGVCF'}{'value'},  #GATK GenoTypeGVCFs gVCF reference infile list for joint genotyping
	   'pGvR|pGATKVariantRecalibration:n' => \$parameter{'pGATKVariantRecalibration'}{'value'},  #GATK VariantRecalibrator/ApplyRecalibration
	   'gvrtsh|GATKVariantReCalibrationTrainingSetHapMap:s' => \$parameter{'GATKVariantReCalibrationTrainingSetHapMap'}{'value'},  #GATK VariantRecalibrator resource
	   'gvrtss|GATKVariantReCalibrationTrainingSetDbSNP:s' => \$parameter{'GATKVariantReCalibrationTrainingSetDbSNP'}{'value'},  #GATK VariantRecalibrator resource
	   'gvrtsg|GATKVariantReCalibrationTrainingSet1000GSNP:s' => \$parameter{'GATKVariantReCalibrationTrainingSet1000GSNP'}{'value'},  #GATK VariantRecalibrator resource
	   'gvrtso|GATKVariantReCalibrationTrainingSet1000GOmni:s' => \$parameter{'GATKVariantReCalibrationTrainingSet1000GOmni'}{'value'},  #GATK VariantRecalibrator resource
	   'gvrtsm|GATKVariantReCalibrationTrainingSetMills:s' => \$parameter{'GATKVariantReCalibrationTrainingSetMills'}{'value'},  #GATK VariantRecalibrator resource
	   'gvrtsf|GATKVariantReCalibrationTSFilterLevel:s' => \$parameter{'GATKVariantReCalibrationTSFilterLevel'}{'value'},  #Truth sensativity level
	   'gvrevf|GATKVariantReCalibrationexcludeNonVariantsFile:s' => \$parameter{'GATKVariantReCalibrationexcludeNonVariantsFile'}{'value'},  #Produce a vcf containing non-variant loci alongside the vcf only containing non-variant loci after GATK VariantRecalibrator (defaults to "false")
	   'pGpT|pGATKPhaseByTransmission:n' => \$parameter{'pGATKPhaseByTransmission'}{'value'},  #GATK PhaseByTransmission to produce phased genotype calls
	   'pGrP|pGATKReadBackedPhasing:n' => \$parameter{'pGATKReadBackedPhasing'}{'value'},  #GATK ReadBackedPhasing
	   'grpqth|GATKReadBackedPhasingPhaseQualityThreshold:n' => \$parameter{'GATKReadBackedPhasingPhaseQualityThreshold'}{'value'},  #quality score required to output phasing
	   'pGvEA|pGATKVariantEvalAll:n' => \$parameter{'pGATKVariantEvalAll'}{'value'},  #GATK varianteval all variants
	   'pGvEE|pGATKVariantEvalExome:n' => \$parameter{'pGATKVariantEvalExome'}{'value'},  #GATK varianteval only exonic variants
	   'gveedbs|GATKVariantEvalDbSNP:s' => \$parameter{'GATKVariantEvalDbSNP'}{'value'},
	   'gveedbg|GATKVariantEvalGold:s' => \$parameter{'GATKVariantReCalibrationTrainingSetMills'}{'value'},
	   'pVeP|pVariantEffectPredictor:n' => \$parameter{'pVariantEffectPredictor'}{'value'},  #Annotation of variants using vep
	   'vepp|vepDirectoryPath:s'  => \$parameter{'vepDirectoryPath'}{'value'},  #path to vep script dir
	   'vepc|vepDirectoryCache:s'  => \$parameter{'vepDirectoryCache'}{'value'},  #path to vep cache dir
	   'vepf|vepFeatures:s'  => \@{$parameter{'vepFeatures'}{'value'}},  #Comma separated list
	   'pVcP|pVCFParser:n' => \$parameter{'pVCFParser'}{'value'},
	   'vcpvt|vcfParserVepTranscripts:n' => \$parameter{'vcfParserVepTranscripts'}{'value'},
	   'vcprff|vcfParserRangeFeatureFile:s'  => \$parameter{'vcfParserRangeFeatureFile'}{'value'},  #path to vcfParserRangeFeatureFile
	   'vcprfa|vcfParserRangeFeatureAnnotationColumns:s'  => \@{$parameter{'vcfParserRangeFeatureAnnotationColumns'}{'value'}},  #Comma separated list
	   'vcpsf|vcfParserSelectFile:s'  => \$parameter{'vcfParserSelectFile'}{'value'},  #path to vcfParserSelectFile
	   'vcpsfm|vcfParserSelectFileMatchingColumn:n' => \$parameter{'vcfParserSelectFileMatchingColumn'}{'value'},  #Column of HGNC Symbol in SelectFile
	   'vcpsfa|vcfParserSelectFeatureAnnotationColumns:s'  => \@{$parameter{'vcfParserSelectFeatureAnnotationColumns'}{'value'}},  #Comma separated list
	   'snep|snpEffPath:s'  => \$parameter{'snpEffPath'}{'value'},  #path to snpEff directory
	   'pSnE|pSnpEff:n' => \$parameter{'pSnpEff'}{'value'},
	   'snesaf2|snpSiftAnnotationFiles:s'  => \$parameter{'snpSiftAnnotationFiles'}{'value'},  #Comma separated list hash entry
	   'snesdbnsfp|snpSiftDbNSFPFile:s'  => \$parameter{'snpSiftDbNSFPFile'}{'value'},  #DbNSFP file
	   'snesdbnsfpa|snpSiftDbNSFPAnnotations:s'  => \@{$parameter{'snpSiftDbNSFPAnnotations'}{'value'}},  #Comma separated list
	   'pAnV|pAnnovar:n' => \$parameter{'pAnnovar'}{'value'},  #Performs annovar filter gene, region and filter analysis
	   'anvp|annovarPath:s'  => \$parameter{'annovarPath'}{'value'},  #path to annovar script dir
	   'anvgbv|annovarGenomeBuildVersion:s'  => \$parameter{'annovarGenomeBuildVersion'}{'value'},
	   'anvtn|annovarTableNames:s'  => \@{$parameter{'annovarTableNames'}{'value'}},  #Comma separated list
	   'anvstn|annovarSupportedTableNames:n' => \$parameter{'annovarSupportedTableNames'}{'value'},  #Generates a list of supported table names
	   'anvarmafth|annovarMAFThreshold:n' => \$parameter{'annovarMAFThreshold'}{'value'},
	   'pRaV|pRankVariants:n' => \$parameter{'pRankVariants'}{'value'},  #Ranking variants
	   'ravgf|geneFile:s' => \$parameter{'geneFile'}{'value'},
	   'ravcs|caddWGSSNVs:n' => \$parameter{'caddWGSSNVs'}{'value'},
	   'ravcsf|caddWGSSNVsFile:s' => \$parameter{'caddWGSSNVsFile'}{'value'},
	   'ravc1kg|cadd1000Genomes:n' => \$parameter{'cadd1000Genomes'}{'value'},
	   'ravc1kgf|cadd1000GenomesFile:s' => \$parameter{'cadd1000GenomesFile'}{'value'},
	   'ravwg|wholeGene:n'  => \$parameter{'wholeGene'}{'value'},  #Allow compound pairs in intronic regions
	   'ravrm|rankModelFile:s' => \$parameter{'rankModelFile'}{'value'},  #The rank modell config.ini path
	   'pScK|pSampleCheck:n' => \$parameter{'pSampleCheck'}{'value'},  #QC for samples gender and relationship
	   'pQcC|pQCCollect:n' => \$parameter{'pQCCollect'}{'value'},  #QCmetrics collect
	   'qccsi|QCCollectSampleInfoFile:s' => \$parameter{'QCCollectSampleInfoFile'}{'value'},  #SampleInfo yaml file produced by MIP
	   'qccref|QCCollectRegExpFile:s' => \$parameter{'QCCollectRegExpFile'}{'value'},  #Regular expression yaml file
	   'pReM|pRemoveRedundantFiles:n' => \$parameter{'pRemoveRedundantFiles'}{'value'},
	   'pArS|pAnalysisRunStatus:n' => \$parameter{'pAnalysisRunStatus'}{'value'},  #AnalysisRunStatus change flag in sampleInfo file if allowed to execute
    );

if($help) {

    print STDOUT $USAGE, "\n";
    exit;
}

my $mipVersion = "v2.0.0";  #Set version

if($version) {

    print STDOUT "\nMip.pl ".$mipVersion, "\n\n";
    exit;
}

if ($parameter{'annovarSupportedTableNames'}{'value'} eq 1) {

    &PrintSupportedAnnovarTableNames(\%scriptParameter, \@annovarSupportedTableNames);
}

if ($parameter{'configFile'}{'value'} ne "nocmdinput") {  #Input from cmd

    %scriptParameter = &LoadYAML(\%scriptParameter, $parameter{'configFile'}{'value'});  #Load parameters from configfile

    &ReplaceConfigParamWithCMDInfo(\%parameter, \%scriptParameter, "analysisType");
    &ReplaceConfigParamWithCMDInfo(\%parameter, \%scriptParameter, "aligner");

    foreach my $orderParameterElement (@orderParameters) {  #Loop through all parameters and update info   

	&UpdateYAML(\%scriptParameter, \$orderParameterElement, \$parameter{'familyID'}{'value'});
    }
}

###Populate scriptParameters{'parameterName'} => 'Value'
foreach my $orderParameterElement (@orderParameters) {
    
##3 type of variables: MIP, path or program/program_parameters each is handled in the &AddToScriptParameter subroutine.
    &AddToScriptParameter(\%parameter, \%scriptParameter, \%sampleInfo, \%fileInfo, \%referenceFileEndings, \@broadcasts,
			   {'parameterName' => $orderParameterElement,
			    'parameterValue' => $parameter{$orderParameterElement}{'value'},
			    'parameterType' => $parameter{$orderParameterElement}{'type'},
			    'parameterDefault' => $parameter{$orderParameterElement}{'default'},
			    'associatedPrograms' => $parameter{$orderParameterElement}{'associatedProgram'},
			    'parameterExistsCheck' => $parameter{$orderParameterElement}{'existsCheck'},
			    'programNamePath' => $parameter{$orderParameterElement}{'programNamePath'}
			   });

    ##Special case for parameters that are dependent on other parameters values
    if ($orderParameterElement eq "outDataDir") {  #Set defaults depending on $scriptParameter{'outDataDir'} value that now has been set

	$parameter{'sampleInfoFile'}{'default'} = $scriptParameter{'outDataDir'}."/".$scriptParameter{'familyID'}."/".$scriptParameter{'familyID'}."_qc_sampleInfo.yaml";
	$parameter{'logFile'}{'default'} = &DeafultLog4perlFile(\%scriptParameter, \$parameter{'logFile'}{'value'}, \$script, \$base, \$dateTimeStamp);

	$parameter{'QCCollectSampleInfoFile'}{'default'} = $parameter{'sampleInfoFile'}{'default'};
    }
    if ($orderParameterElement eq "logFile") {

	###Creates log for the master script
	my $conf = &CreateLog4perlCongfig(\$scriptParameter{'logFile'});
	Log::Log4perl->init(\$conf);
	$logger = Log::Log4perl->get_logger("MIPLogger");
    }
    if ($orderParameterElement eq "pedigreeFile") {  #Write QC for only pedigree data used in analysis                                                        
	
	if (defined($scriptParameter{'pedigreeFile'})) {

	    `mkdir -p $scriptParameter{'outDataDir'}/$scriptParameter{'familyID'};`;  #Create family directory
	    my $yamlFile = $scriptParameter{'outDataDir'}."/".$scriptParameter{'familyID'}."/qc_pedigree.yaml";
	    &WriteYAML(\%sampleInfo, \$yamlFile);
	    &RemovePedigreeElements(\%sampleInfo);  #NOTE: Removes all elements at hash third level except 'Capture_kit'
	}
    }
    if ($orderParameterElement eq "humanGenomeReference") {  #Supply humanGenomeReference to mosaikAlignReference if required
	
	if ( (defined($scriptParameter{'humanGenomeReference'})) && (defined($fileInfo{'humanGenomeReferenceNameNoEnding'})) ) {

	    &SetAutoBuildFeature(\%scriptParameter, "mosaikAlignReference", \$referenceFileEndings{'mosaikAlignReference'}, \$fileInfo{'humanGenomeReferenceNameNoEnding'});
	    &SetAutoBuildFeature(\%scriptParameter, "mosaikJumpDbStub", \$referenceFileEndings{'mosaikJumpDbStub'}, \$fileInfo{'humanGenomeReferenceNameNoEnding'});
	    &SetAutoBuildFeature(\%scriptParameter, "bwaBuildReference", \$referenceFileEndings{'bwaBuildReference'}, \$scriptParameter{'humanGenomeReference'});	
	}
    }
} 


##sampleIDs
&DefineArrayParameters(\%parameter, \%scriptParameter, \@{$parameter{'sampleIDs'}{'value'}}, \@orderParameters, \@broadcasts, "sampleIDs", "path", "nodefault", "MIP", "");

&CheckUniqueIDNs(\%scriptParameter, \@{$scriptParameter{'sampleIDs'}});  #Test that sampleIDs are unique

for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{$scriptParameter{'sampleIDs'}});$sampleIDCounter++) {  #all sampleIDs
    
    &ScriptParameterPerSampleID(\%scriptParameter, \$scriptParameter{'familyID'}, \$scriptParameter{'sampleIDs'}[$sampleIDCounter], "exomeTargetBedInfileLists");
    &ScriptParameterPerSampleID(\%scriptParameter, \$scriptParameter{'familyID'}, \$scriptParameter{'sampleIDs'}[$sampleIDCounter], "exomeTargetPaddedBedInfileLists");
    &ScriptParameterPerSampleID(\%scriptParameter, \$scriptParameter{'familyID'}, \$scriptParameter{'sampleIDs'}[$sampleIDCounter], "GATKTargetPaddedBedIntervalLists");
}

##inFileDirs
&DefineArrayParameters(\%parameter, \%scriptParameter, \@{$parameter{'inFilesDirs'}{'value'}}, \@orderParameters, \@broadcasts, "inFilesDirs", "path", "notSetYet", "MIP", "directory");

##Compares the number of elements in two arrays and exits if the elements are not equal
&CompareArrayElements(\@{$scriptParameter{'sampleIDs'}}, \@{$scriptParameter{'inFilesDirs'}}, "sampleIDs", "inFileDirs");

##picardToolsMergeSamFilesPrevious
if ($scriptParameter{'pPicardToolsMergeSamFiles'} > 0) {
    
    if( (scalar(@{$parameter{'picardToolsMergeSamFilesPrevious'}{'value'}}) > 0) ) {

	&DefineArrayParameters(\%parameter, \%scriptParameter, \@{$parameter{'picardToolsMergeSamFilesPrevious'}{'value'}}, \@orderParameters, \@broadcasts, "picardToolsMergeSamFilesPrevious", "path", "nodefault", "pPicardToolsMergeSamFiles", "file");    
	
	&CheckMergePicardToolsMergeSamFilesPrevious(\%scriptParameter, \%sampleInfo);
    }
    else {  #Not supplied - Set to 0 to handle correctly in program subroutines 
	
	for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{$scriptParameter{'sampleIDs'}});$sampleIDCounter++) {  #Set for all sampleIDs
	    
	    $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'sampleIDs'}[$sampleIDCounter] }{'picardToolsMergeSamFilesPrevious'} = 0;
	}
    }
}


##pVariantEffectPredictor
if ($scriptParameter{'pVariantEffectPredictor'} > 0) {

    &DefineArrayParameters(\%parameter, \%scriptParameter, \@{$parameter{'vepFeatures'}{'value'}}, \@orderParameters, \@broadcasts, "vepFeatures", "path", "yes", "pVariantEffectPredictor", "");
}


##pVCFParser
if ($scriptParameter{'pVCFParser'} > 0) {
    
    &DefineArrayParameters(\%parameter, \%scriptParameter, \@{$parameter{'vcfParserRangeFeatureAnnotationColumns'}{'value'}}, \@orderParameters, \@broadcasts, "vcfParserRangeFeatureAnnotationColumns", "path", "nodefault", "pVCFParser", "");
    &DefineArrayParameters(\%parameter, \%scriptParameter, \@{$parameter{'vcfParserSelectFeatureAnnotationColumns'}{'value'}}, \@orderParameters, \@broadcasts, "vcfParserSelectFeatureAnnotationColumns", "path", "nodefault", "pVCFParser", "");
}


##pSnpEff
if ($scriptParameter{'pSnpEff'} > 0) {

    &DefineHashParameters(\%parameter, \%scriptParameter, \@orderParameters, \@broadcasts, "snpSiftAnnotationFiles", "path", "yes", "pSnpEff", "file");
    &DefineArrayParameters(\%parameter, \%scriptParameter, \@{$parameter{'snpSiftDbNSFPAnnotations'}{'value'}}, \@orderParameters, \@broadcasts, "snpSiftDbNSFPAnnotations", "path", "yes", "pSnpEff", "");  #"yes" added to enable addition of default features in &AddToScriptParameters  
}


##pAnnovar
if ($scriptParameter{'pAnnovar'} > 0) {

    %annovarTable = &DefineAnnovarTables(\%parameter, \$scriptParameter{'annovarGenomeBuildVersion'}); #Set all AnnovarTables properties
    &DefineArrayParameters(\%parameter, \%scriptParameter, \@{$parameter{'annovarTableNames'}{'value'}}, \@orderParameters, \@broadcasts, "annovarTableNames", "path", "yes", "pAnnovar", "file");  
}


##Set Target files
&PrepareArrayParameters(\%parameter, \@exomeTargetBedInfileLists, \@orderParameters, \@broadcasts, "exomeTargetBedInfileLists", "path", "notSetYet", "pPicardToolsCalculateHSMetrics", "file");

&PrepareArrayParameters(\%parameter, \@exomeTargetPaddedBedInfileLists, \@orderParameters, \@broadcasts, "exomeTargetPaddedBedInfileLists", "path", "notSetYet", "pPicardToolsCalculateHSMetrics", "file");
 
&PrepareArrayParameters(\%parameter, \@GATKTargetPaddedBedIntervalLists, \@orderParameters, \@broadcasts, "GATKTargetPaddedBedIntervalLists", "path", "notSetYet", "pGATKHaploTypeCaller,pGATKVariantRecalibration", "file");

##Broadcast set parameters info
foreach my $parameterInfo (@broadcasts) {

    $logger->info($parameterInfo, "\n");
}

##Cosmid references
&DefineSupportedCosmidReferences(\%supportedCosmidReference, "humanGenomeReference", "decoy", "5", \$fileInfo{'humanGenomeReferenceVersion'}, "compressed");
&DefineSupportedCosmidReferences(\%supportedCosmidReference, "chanjoBuildDb", "ccds", "latest", \$fileInfo{'humanGenomeReferenceVersion'}, "unCompressed");
&DefineSupportedCosmidReferences(\%supportedCosmidReference, "GATKReAlignerINDELKnownSet1", "indels", $scriptParameter{'GATKBundleDownLoadVersion'}."/b".$fileInfo{'humanGenomeReferenceVersion'}, \$fileInfo{'humanGenomeReferenceVersion'}, "unCompressed");
&DefineSupportedCosmidReferences(\%supportedCosmidReference, "GATKReAlignerINDELKnownSet2", "mills", $scriptParameter{'GATKBundleDownLoadVersion'}."/b".$fileInfo{'humanGenomeReferenceVersion'}, \$fileInfo{'humanGenomeReferenceVersion'}, "unCompressed");
&DefineSupportedCosmidReferences(\%supportedCosmidReference, "GATKBaseReCalibrationSNPKnownSet", "dbsnp", $scriptParameter{'GATKBundleDownLoadVersion'}."/b".$fileInfo{'humanGenomeReferenceVersion'}, \$fileInfo{'humanGenomeReferenceVersion'}, "unCompressed");
&DefineSupportedCosmidReferences(\%supportedCosmidReference, "GATKHaploTypeCallerSNPKnownSet", "dbsnp", $scriptParameter{'GATKBundleDownLoadVersion'}."/b".$fileInfo{'humanGenomeReferenceVersion'}, \$fileInfo{'humanGenomeReferenceVersion'}, "unCompressed");
&DefineSupportedCosmidReferences(\%supportedCosmidReference, "GATKVariantReCalibrationTrainingSetHapMap", "hapmap", $scriptParameter{'GATKBundleDownLoadVersion'}."/b".$fileInfo{'humanGenomeReferenceVersion'}, \$fileInfo{'humanGenomeReferenceVersion'}, "unCompressed");
&DefineSupportedCosmidReferences(\%supportedCosmidReference, "GATKVariantReCalibrationTrainingSetMills", "mills", $scriptParameter{'GATKBundleDownLoadVersion'}."/b".$fileInfo{'humanGenomeReferenceVersion'}, \$fileInfo{'humanGenomeReferenceVersion'}, "unCompressed");
&DefineSupportedCosmidReferences(\%supportedCosmidReference, "GATKVariantReCalibrationTrainingSet1000GOmni", "1000g_omni", $scriptParameter{'GATKBundleDownLoadVersion'}."/b".$fileInfo{'humanGenomeReferenceVersion'}, \$fileInfo{'humanGenomeReferenceVersion'}, "unCompressed");
&DefineSupportedCosmidReferences(\%supportedCosmidReference, "GATKVariantReCalibrationTrainingSet1000GSNP", "1000g_snps", $scriptParameter{'GATKBundleDownLoadVersion'}."/b".$fileInfo{'humanGenomeReferenceVersion'}, \$fileInfo{'humanGenomeReferenceVersion'}, "unCompressed");
&DefineSupportedCosmidReferences(\%supportedCosmidReference, "GATKVariantReCalibrationTrainingSetDbSNP", "dbsnp", $scriptParameter{'GATKBundleDownLoadVersion'}."/b".$fileInfo{'humanGenomeReferenceVersion'}, \$fileInfo{'humanGenomeReferenceVersion'}, "unCompressed");
&DefineSupportedCosmidReferences(\%supportedCosmidReference, "GATKVariantEvalGold", "mills", $scriptParameter{'GATKBundleDownLoadVersion'}."/b".$fileInfo{'humanGenomeReferenceVersion'}, \$fileInfo{'humanGenomeReferenceVersion'}, "unCompressed"); 
&DefineSupportedCosmidReferences(\%supportedCosmidReference, "GATKVariantEvalDbSNP", "dbsnpex", $scriptParameter{'GATKBundleDownLoadVersion'}."/b".$fileInfo{'humanGenomeReferenceVersion'}, \$fileInfo{'humanGenomeReferenceVersion'}, "unCompressed");

##Flag -> array parameters to enable multiple download via Cosmid using the same flag 
&DefineSupportedCosmidReferences(\%supportedCosmidReference, "dbsnp_138.b37.vcf", "dbsnp", $scriptParameter{'GATKBundleDownLoadVersion'}."/b".$fileInfo{'humanGenomeReferenceVersion'}, \$fileInfo{'humanGenomeReferenceVersion'}, "unCompressed");
&DefineSupportedCosmidReferences(\%supportedCosmidReference, "dbsnp_138.b37.excluding_sites_after_129.vcf", "dbsnpex", $scriptParameter{'GATKBundleDownLoadVersion'}."/b".$fileInfo{'humanGenomeReferenceVersion'}, \$fileInfo{'humanGenomeReferenceVersion'}, "unCompressed");
&DefineSupportedCosmidReferences(\%supportedCosmidReference, "1000G_phase1.indels.b37.vcf", "1000g_omni", $scriptParameter{'GATKBundleDownLoadVersion'}."/b".$fileInfo{'humanGenomeReferenceVersion'}, \$fileInfo{'humanGenomeReferenceVersion'}, "unCompressed");
&DefineSupportedCosmidReferences(\%supportedCosmidReference, "1000G_phase1.snps.high_confidence.b37.vcf", "1000g_snps", $scriptParameter{'GATKBundleDownLoadVersion'}."/b".$fileInfo{'humanGenomeReferenceVersion'}, \$fileInfo{'humanGenomeReferenceVersion'}, "unCompressed");


for my $references (keys %supportedCosmidReference) {

    &CheckCosmidInstallation(\%parameter, \%scriptParameter, \$references);
    last;  #Only need to check once per analysis run
}

if ($scriptParameter{'writeConfigFile'} ne 0) {  #Write config file for family

    &WriteYAML(\%scriptParameter, \$scriptParameter{'writeConfigFile'});  #Write used settings to configfile
}

##Set chr prefix and chromosome names depending on reference used
if ($scriptParameter{'humanGenomeReference'}=~/hg\d+/) {  #Refseq - prefix and M

    @contigs = ("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM");  #Chr for filtering of bam file
}
elsif ($scriptParameter{'humanGenomeReference'}=~/GRCh\d+/) {  #Ensembl - no prefix and MT

    @contigs = ("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT");  #Chr for filtering of bam file
}

## Write CMD to MIP log file
&WriteCMDMipLog(\%parameter, \%scriptParameter, \@orderParameters, \$script, \$scriptParameter{'logFile'}, \$mipVersion);

## Collects the ".fastq(.gz)" files from the supplied infiles directory. Checks if any of the files exist
&CollectInfiles(\%scriptParameter, \%inDirPath, \%infile);

## Reformat files for MIP output, which have not yet been created into, correct format so that a sbatch script can be generated with the correct filenames
my $uncompressedFileSwitch = &InfilesReFormat(\%infile);  #Required to format infiles correctly for subsequent input into aligners
    
&CreateFileEndings(\%parameter, \%scriptParameter, \%sampleInfo, \%infilesLaneNoEnding, \@orderParameters);  #Creates all fileendings as the samples is processed depending on the chain of modules activated

## Create .fam file to be used in variant calling analyses
&CreateFamFile(\%scriptParameter);

####MAIN

if ( ($scriptParameter{'pGZipFastq'} > 0) && ($uncompressedFileSwitch eq "unCompressed") ) {  #GZip of fastq files

    $logger->info("\n[GZip for fastq files]\n");

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{$scriptParameter{'sampleIDs'}});$sampleIDCounter++) {  

	for (my $infileCounter=0;$infileCounter<scalar( @{ $infile{$scriptParameter{'sampleIDs'}[$sampleIDCounter]} });$infileCounter++) {  #To determine which sampleID had the uncompressed files
	    
	    if ($infile{$scriptParameter{'sampleIDs'}[$sampleIDCounter]}[$infileCounter] =~/.fastq$/) {
	
		## Automatically gzips fastq files
		&GZipFastq(\%parameter, \%scriptParameter, \%infile, \%inDirPath, $scriptParameter{'sampleIDs'}[$sampleIDCounter]);
		last;  #Return to sampleID loop i.e. only call subroutine GZipFastq once per sampleID
	    }
	}
    }
}

if ($scriptParameter{'pFastQC'} > 0) {  #Run FastQC
    
    $logger->info("[FastQC]\n");

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{$scriptParameter{'sampleIDs'}});$sampleIDCounter++) {  
	
	&FastQC(\%parameter, \%scriptParameter, \%sampleInfo, \%infile, \%inDirPath, \%infilesLaneNoEnding, \%infilesBothStrandsNoEnding, $scriptParameter{'sampleIDs'}[$sampleIDCounter], "FastQC");	
    }
}

if ($scriptParameter{'pMosaikBuild'} > 0) {  #Run MosaikBuild
    
    $logger->info("[MosaikBuild]\n");

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{$scriptParameter{'sampleIDs'}});$sampleIDCounter++) {  
	
	&MosaikBuild(\%parameter, \%scriptParameter, \%sampleInfo, \%infile, \%inDirPath, \%infilesLaneNoEnding, \%lane, $scriptParameter{'sampleIDs'}[$sampleIDCounter], $scriptParameter{'aligner'}, "MosaikBuild");	
    }
}


if ($scriptParameter{'pMosaikAlign'} > 0) {  #Run MosaikAlign

    $logger->info("[MosaikAlign]\n");

    if ($scriptParameter{'dryRunAll'} != 1) {

	if ( ($parameter{'humanGenomeReference'}{'buildFile'} eq 1) || ($parameter{'mosaikAlignReference'}{'buildFile'} eq 1) || ($parameter{'mosaikJumpDbStub'}{'buildFile'} eq 1) ) {
	    
	    &BuildMosaikAlignPreRequisites(\%parameter, \%scriptParameter, \@mosaikJumpDbStubFileEndings, \$fileInfo{'humanGenomeReferenceSource'}, \$fileInfo{'humanGenomeReferenceVersion'}, $scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "MosaikAlign");
	    
	}
	if ( ($parameter{'mosaikAlignNeuralNetworkPeFile'}{'buildFile'} eq 1) || ($parameter{'mosaikAlignNeuralNetworkSeFile'}{'buildFile'} eq 1) ){
	    
	    &MoveMosaikNN(\%scriptParameter);
	}
    }
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{$scriptParameter{'sampleIDs'}});$sampleIDCounter++) {  
	
	&MosaikAlign(\%parameter, \%scriptParameter, \%sampleInfo, \%infile, \%inDirPath, \%infilesLaneNoEnding, \%infilesBothStrandsNoEnding, $scriptParameter{'sampleIDs'}[$sampleIDCounter], $scriptParameter{'aligner'}, "MosaikAlign");	
    }
}

if ($scriptParameter{'pBwaMem'} > 0) {  #Run BWA Mem
    
    $logger->info("[BWA Mem]\n");
    
    if ($scriptParameter{'dryRunAll'} != 1) {

	if ( ($parameter{'humanGenomeReference'}{'buildFile'} eq 1) || ($parameter{'bwaBuildReference'}{'buildFile'} eq 1) ) {
	    
	    &BuildBwaPreRequisites(\%parameter, \%scriptParameter, \@bwaBuildReferenceFileEndings, $scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BwaMem");
	}
    }
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{$scriptParameter{'sampleIDs'}});$sampleIDCounter++) {  
	
	&BWA_Mem(\%parameter, \%scriptParameter, \%sampleInfo, \%infile, \%inDirPath, \%infilesLaneNoEnding, \%infilesBothStrandsNoEnding, $scriptParameter{'sampleIDs'}[$sampleIDCounter], $scriptParameter{'aligner'}, "BwaMem");	
	
    }    
}

if ($scriptParameter{'pPicardToolsMergeRapidReads'} > 0) {  #Run PicardToolsMergeRapidReads - Relevant only in rapid mode
    
    $logger->info("[PicardToolsMergeRapidReads]\n");
    
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{$scriptParameter{'sampleIDs'}});$sampleIDCounter++) {  
	
        #Merge all read batch processes to 1 file again containing sorted & indexed reads matching clinical test genes
	&PicardToolsMergeRapidReads(\%parameter, \%scriptParameter, \%sampleInfo, \%infilesLaneNoEnding, \%lane, $scriptParameter{'sampleIDs'}[$sampleIDCounter], $scriptParameter{'aligner'}, "PicardToolsMergeRapidReads");
    }    
}

if ($scriptParameter{'pBwaAln'} > 0) {  #Run BWA Aln
    
    $logger->info("[BWA Aln]\n");

    if ($scriptParameter{'dryRunAll'} != 1) {

	if ( ($parameter{'humanGenomeReference'}{'buildFile'} eq 1) || ($parameter{'bwaBuildReference'}{'buildFile'} eq 1) ) {
	    
	    &BuildBwaPreRequisites(\%parameter, \%scriptParameter, \@bwaBuildReferenceFileEndings, $scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BwaAln");
	}
    }
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{$scriptParameter{'sampleIDs'}});$sampleIDCounter++) {  
	
	&BWA_Aln(\%parameter, \%scriptParameter, \%sampleInfo, \%infile, \%inDirPath, \%infilesLaneNoEnding, \%infilesBothStrandsNoEnding, $scriptParameter{'sampleIDs'}[$sampleIDCounter], $scriptParameter{'aligner'}, "BwaAln");	
    }    
}

if ($scriptParameter{'pBwaSampe'} > 0) {  #Run BWA Sampe
    
    $logger->info("[BWA Sampe]\n");

    if ($scriptParameter{'dryRunAll'} != 1) {

	if ( ($parameter{'humanGenomeReference'}{'buildFile'} eq 1) || ($parameter{'bwaBuildReference'}{'buildFile'} eq 1) ) {
	    
	    &BuildBwaPreRequisites(\%parameter, \%scriptParameter, \@bwaBuildReferenceFileEndings, $scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BwaSampe");
	}
    }
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{$scriptParameter{'sampleIDs'}});$sampleIDCounter++) {  
	
	&BWA_Sampe(\%parameter, \%scriptParameter, \%sampleInfo, \%infile, \%inDirPath, \%infilesLaneNoEnding, \%infilesBothStrandsNoEnding, $scriptParameter{'sampleIDs'}[$sampleIDCounter], $scriptParameter{'aligner'}, "BwaSampe");
    }
}

if ($scriptParameter{'pPicardToolsSortSam'} > 0) {  #Run Picardtools SortSam and Index

    if ($scriptParameter{'analysisType'} ne "rapid") {  #In rapid mode Sort and index is done for each batch of reads in the BWA_Mem call, since the link to infile is broken by the read batch processing. However pPicardToolsSortSam should be enabled to ensure correct fileending and merge the flow to ordinary modules.

    $logger->info("[PicardTools SortSam & index]\n");
	
	for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{$scriptParameter{'sampleIDs'}});$sampleIDCounter++) {  
	    
	    &PicardToolsSortSamIndex(\%parameter, \%scriptParameter, \%sampleInfo, \%infile, \%inDirPath, \%infilesLaneNoEnding, $scriptParameter{'sampleIDs'}[$sampleIDCounter], $scriptParameter{'aligner'}, "PicardToolsSortSam");
	}
    }
}

if ($scriptParameter{'pPicardToolsMergeSamFiles'} > 0) {  #Run picardtools merge

    $logger->info("[PicardTool MergeSamFiles]\n");

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{$scriptParameter{'sampleIDs'}});$sampleIDCounter++) {  

	if ( ($sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'sampleIDs'}[$sampleIDCounter] }{'picardToolsMergeSamFilesPrevious'} == 1) || (scalar( @{ $infilesLaneNoEnding{ $scriptParameter{'sampleIDs'}[$sampleIDCounter] } }) > 1) ) {  #Sanity Check that we have something to merge with
	
	    &PicardToolsMerge(\%parameter, \%scriptParameter, \%sampleInfo, \%infilesLaneNoEnding, \%lane, $scriptParameter{'sampleIDs'}[$sampleIDCounter], $scriptParameter{'aligner'}, $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'sampleIDs'}[$sampleIDCounter] }{'fileEnding'}, "PicardToolsMergeSamFiles");	
	}
    }
}

if ($scriptParameter{'pPicardToolsMarkduplicates'} > 0) {  #PicardTools MarkDuplicates

    $logger->info("[PicardTools MarkDuplicates]\n");

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{$scriptParameter{'sampleIDs'}});$sampleIDCounter++) {  
    
	&PicardToolsMarkDuplicates(\%parameter, \%scriptParameter, \%sampleInfo, \%infilesLaneNoEnding, \%infilesBothStrandsNoEnding, \%lane, $scriptParameter{'sampleIDs'}[$sampleIDCounter], $scriptParameter{'aligner'}, "PicardToolsMarkduplicates");	
    }
}

if ($scriptParameter{'pChanjoSexCheck'} > 0) {
    
    $logger->info("[ChanjoSexCheck]\n");
    
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{$scriptParameter{'sampleIDs'}});$sampleIDCounter++) {  #For all SampleIDs
	
	&ChanjoSexCheck(\%parameter, \%scriptParameter, \%sampleInfo, \%infilesLaneNoEnding, \%lane, $scriptParameter{'sampleIDs'}[$sampleIDCounter], $scriptParameter{'aligner'}, "ChanjoSexCheck");
    }
}

if ($scriptParameter{'pChanjoBuild'} > 0) {
    
    $logger->info("[ChanjoBuild]\n");

    &CheckBuildDownLoadPreRequisites(\%parameter, \%scriptParameter, \%supportedCosmidReference, "ChanjoBuild");
       
    &ChanjoBuild(\%parameter, \%scriptParameter, \%sampleInfo, $scriptParameter{'familyID'}, "ChanjoBuild");
}

if ($scriptParameter{'pChanjoAnnotate'} > 0) {
    
    $logger->info("[ChanjoAnnotate]\n");
    
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{$scriptParameter{'sampleIDs'}});$sampleIDCounter++) {  #For all SampleIDs
	
	&ChanjoAnnotate(\%parameter, \%scriptParameter, \%sampleInfo, \%infilesLaneNoEnding, \%lane, $scriptParameter{'sampleIDs'}[$sampleIDCounter], $scriptParameter{'aligner'}, "ChanjoAnnotate");
    }
}

if ($scriptParameter{'pChanjoImport'} > 0) {
    
    $logger->info("[ChanjoImport]\n");
    
    &ChanjoImport(\%parameter, \%scriptParameter, \%sampleInfo, $scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "ChanjoImport");
}

if ($scriptParameter{'pGenomeCoverageBED'} > 0) {  #Run GenomeCoverageBED
    
    $logger->info("[GenomeCoverageBED]\n"); 
    
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{$scriptParameter{'sampleIDs'}});$sampleIDCounter++) {  

	&GenomeCoverageBED(\%parameter, \%scriptParameter, \%sampleInfo, \%infilesLaneNoEnding, \%lane, $scriptParameter{'sampleIDs'}[$sampleIDCounter], $scriptParameter{'aligner'}, "GenomeCoverageBED");
    }
}

if ($scriptParameter{'pPicardToolsCollectMultipleMetrics'} > 0) {  #Run PicardToolsCollectMultipleMetrics
    
    $logger->info("[PicardToolsCollectMultipleMetrics]\n");   
    
    &CheckBuildHumanGenomePreRequisites(\%parameter, \%scriptParameter, \%fileInfo, \@humanGenomeReferenceFileEndings, "PicardToolsCollectMultipleMetrics");

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{$scriptParameter{'sampleIDs'}});$sampleIDCounter++) {  

	&PicardToolsCollectMultipleMetrics(\%parameter, \%scriptParameter, \%sampleInfo, \%infilesLaneNoEnding, \%lane, $scriptParameter{'sampleIDs'}[$sampleIDCounter], $scriptParameter{'aligner'}, "PicardToolsCollectMultipleMetrics");
    }
}

if ($scriptParameter{'pPicardToolsCalculateHSMetrics'} > 0) {  #Run PicardToolsCalculateHSMetrics
    
    $logger->info("[PicardToolsCalculateHSMetrics]\n");   
    
    &CheckBuildHumanGenomePreRequisites(\%parameter, \%scriptParameter, \%fileInfo, \@humanGenomeReferenceFileEndings, "PicardToolsCalculateHSMetrics");
    if ($scriptParameter{'dryRunAll'} != 1) {

	&CheckBuildPTCHSMetricPreRequisites(\%parameter, \%scriptParameter, "PicardToolsCalculateHSMetrics");
    }
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{$scriptParameter{'sampleIDs'}});$sampleIDCounter++) {  

	&PicardToolsCalculateHSMetrics(\%parameter, \%scriptParameter, \%sampleInfo, \%infilesLaneNoEnding, \%lane, $scriptParameter{'sampleIDs'}[$sampleIDCounter], $scriptParameter{'aligner'}, "PicardToolsCalculateHSMetrics");
    }
}

if ($scriptParameter{'pRCovPlots'} > 0) {  #Run Rcovplot scripts   

    $logger->info("[RCovPlots]\n");	

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{$scriptParameter{'sampleIDs'}});$sampleIDCounter++) {  
	
	&RCoveragePlots(\%parameter, \%scriptParameter, \%sampleInfo, $scriptParameter{'sampleIDs'}[$sampleIDCounter], $scriptParameter{'aligner'}, "RCovPlots");	
    }
}

if ($scriptParameter{'pGATKRealigner'} > 0) {  #Run GATK ReAlignerTargetCreator/IndelRealigner

    $logger->info("[GATK ReAlignerTargetCreator/IndelRealigner]\n");

    &CheckBuildHumanGenomePreRequisites(\%parameter, \%scriptParameter, \%fileInfo, \@humanGenomeReferenceFileEndings, "GATKRealigner");
    &CheckBuildDownLoadPreRequisites(\%parameter, \%scriptParameter, \%supportedCosmidReference, "GATKRealigner");

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{$scriptParameter{'sampleIDs'}});$sampleIDCounter++) {   
    
	&GATKReAligner(\%parameter, \%scriptParameter, \%sampleInfo, \%infilesLaneNoEnding, $scriptParameter{'sampleIDs'}[$sampleIDCounter], $scriptParameter{'aligner'}, "GATKRealigner");	
    }
}

if ($scriptParameter{'pGATKBaseRecalibration'} > 0) {  #Run GATK BaseRecalibrator/PrintReads

    $logger->info("[GATK BaseRecalibrator/PrintReads]\n");

    &CheckBuildHumanGenomePreRequisites(\%parameter, \%scriptParameter, \%fileInfo, \@humanGenomeReferenceFileEndings, "GATKBaseRecalibration");
    &CheckBuildDownLoadPreRequisites(\%parameter, \%scriptParameter, \%supportedCosmidReference, "GATKBaseRecalibration");

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{$scriptParameter{'sampleIDs'}});$sampleIDCounter++) {   
  
	&GATKBaseReCalibration(\%parameter, \%scriptParameter, \%sampleInfo, \%infilesLaneNoEnding, $scriptParameter{'sampleIDs'}[$sampleIDCounter], $scriptParameter{'aligner'}, "GATKBaseRecalibration");	
    }
}

if ($scriptParameter{'pGATKHaploTypeCaller'} > 0) {  #Run GATK HaploTypeCaller

    $logger->info("[GATK HaplotypeCaller]\n");

    &CheckBuildHumanGenomePreRequisites(\%parameter, \%scriptParameter, \%fileInfo, \@humanGenomeReferenceFileEndings, "GATKHaploTypeCaller");
    &CheckBuildDownLoadPreRequisites(\%parameter, \%scriptParameter, \%supportedCosmidReference, "GATKHaploTypeCaller");
   
    if ($scriptParameter{'dryRunAll'} != 1) {

	&CheckBuildPTCHSMetricPreRequisites(\%parameter, \%scriptParameter, "GATKHaploTypeCaller");
    }
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{$scriptParameter{'sampleIDs'}});$sampleIDCounter++) {
	
	if ( (defined($parameter{ $scriptParameter{'familyID'} }{$scriptParameter{'sampleIDs'}[$sampleIDCounter]}{'GATKTargetPaddedBedIntervalLists'}{'buildFile'})) && ($parameter{ $scriptParameter{'familyID'} }{$scriptParameter{'sampleIDs'}[$sampleIDCounter]}{'GATKTargetPaddedBedIntervalLists'}{'buildFile'} eq 1) ){
	    
	    if ($scriptParameter{'dryRunAll'} != 1) {
		
		&BuildPTCHSMetricPreRequisites(\%parameter, \%scriptParameter, \%referenceFileEndings, $scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "GATKHaploTypeCaller");
		last;  #Will handle all build per sampleID within sbatch script
	    }
	}
    }
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{$scriptParameter{'sampleIDs'}});$sampleIDCounter++) {
	    
	&GATKHaploTypeCaller(\%parameter, \%scriptParameter, \%sampleInfo, \%infilesLaneNoEnding, \%lane, $scriptParameter{'sampleIDs'}[$sampleIDCounter], $scriptParameter{'aligner'}, "GATKHaploTypeCaller");
    }
}

if ($scriptParameter{'pGATKGenoTypeGVCFs'} > 0) {  #Run GATK GenoTypeGVCFs. Done per family

    $logger->info("[GATK GenoTypeGVCFs]\n");

    &CheckBuildHumanGenomePreRequisites(\%parameter, \%scriptParameter, \%fileInfo, \@humanGenomeReferenceFileEndings, "GATKGenoTypeGVCFs");

    &GATKGenoTypeGVCFs(\%parameter, \%scriptParameter, \%sampleInfo, \%lane, $scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH", "GATKGenoTypeGVCFs");
}

if ($scriptParameter{'pGATKVariantRecalibration'} > 0) {  #Run GATK VariantRecalibrator/ApplyRecalibration. Done per family

    $logger->info("[GATK VariantRecalibrator/ApplyRecalibration]\n");

    &CheckBuildHumanGenomePreRequisites(\%parameter, \%scriptParameter, \%fileInfo, \@humanGenomeReferenceFileEndings, "GATKVariantRecalibration");
    &CheckBuildDownLoadPreRequisites(\%parameter, \%scriptParameter, \%supportedCosmidReference, "GATKVariantRecalibration");
    if ($scriptParameter{'dryRunAll'} != 1) {

	&CheckBuildPTCHSMetricPreRequisites(\%parameter, \%scriptParameter, "GATKVariantRecalibration");
    }
    &GATKVariantReCalibration(\%parameter, \%scriptParameter, \%sampleInfo, $scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH", "GATKVariantRecalibration");
}

if ($scriptParameter{'pSampleCheck'} > 0) {  #Run SampleCheck. Done per family

    $logger->info("[SampleCheck]\n");

    &SampleCheck(\%parameter, \%scriptParameter, \%sampleInfo, $scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH", "SampleCheck");
}

if ($scriptParameter{'pGATKPhaseByTransmission'} > 0) {  #Run GATK PhaseByTransmission. Done per family
    
    $logger->info("[GATK PhaseByTransmission]\n");
    
    &CheckBuildHumanGenomePreRequisites(\%parameter, \%scriptParameter, \%fileInfo, \@humanGenomeReferenceFileEndings, "GATKPhaseByTransmission");
    &GATKPhaseByTransmission(\%parameter, \%scriptParameter, \%sampleInfo, $scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH", "GATKPhaseByTransmission");
}

if ($scriptParameter{'pGATKReadBackedPhasing'} > 0) {  #Run GATK ReadBackedPhasing. Done per family. NOTE: Needs phased calls
    
    $logger->info("[GATK ReadBackedPhasing]\n");
    
    &CheckBuildHumanGenomePreRequisites(\%parameter, \%scriptParameter, \%fileInfo, \@humanGenomeReferenceFileEndings, "GATKReadBackedPhasing");
    &GATKReadBackedPhasing(\%parameter, \%scriptParameter, \%sampleInfo, \%infilesLaneNoEnding, $scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH", "GATKReadBackedPhasing");
}

if ($scriptParameter{'pVariantEffectPredictor'} > 0) {  #Run VariantEffectPredictor. Done per family

    $logger->info("[VariantEffectPredictor]\n");

    &VariantEffectPredictor(\%parameter, \%scriptParameter, \%sampleInfo, $scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH", "VariantEffectPredictor");
}

if ($scriptParameter{'pVCFParser'} > 0) {  #Run VariantEffectPredictor. Done per family

    $logger->info("[VCFParser]\n");
    
    &VCFParser(\%parameter, \%scriptParameter, \%sampleInfo, $scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH", "VCFParser");

    if ($scriptParameter{'vcfParserSelectFile'} ne "noUserInfo") {

	$VEPOutputFiles = 2;  #Use seperate analysis for variants overlapping selected genes and orphans
    }
}

if ($scriptParameter{'pSnpEff'} > 0) {  #Run snpEff. Done per family

    $logger->info("[SnpEff]\n");

    &CheckBuildDownLoadPreRequisites(\%parameter, \%scriptParameter, \%supportedCosmidReference, "SnpEff");
    &SnpEff(\%parameter, \%scriptParameter, \%sampleInfo, \%fileInfo, $scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH", "SnpEff");
}

if ($scriptParameter{'pAnnovar'} > 0) {  #Run Annovar. Done per family

    $logger->info("[Annovar]\n");

    &CheckBuildHumanGenomePreRequisites(\%parameter, \%scriptParameter, \%fileInfo, \@humanGenomeReferenceFileEndings, "Annovar");

    for (my $tableNamesCounter=0;$tableNamesCounter<scalar(@{$scriptParameter{'annovarTableNames'}});$tableNamesCounter++) {  #For all specified table names

	if ($parameter{ $scriptParameter{'annovarTableNames'}[$tableNamesCounter] }{'buildFile'} eq 1) {

	&BuildAnnovarPreRequisites(\%parameter, \%scriptParameter, \%annovarTable, $scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "Annovar");
	last;  #Will handle all build tables within sbatch script
	}
    }
    &Annovar(\%parameter, \%scriptParameter, \%sampleInfo, \%annovarTable, $scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH", "Annovar");
}

if ($scriptParameter{'pGATKVariantEvalAll'} > 0) {  #Run GATK VariantEval for all variants. Done per sampleID

    $logger->info("[GATK VariantEval All]\n");

    &CheckBuildHumanGenomePreRequisites(\%parameter, \%scriptParameter, \%fileInfo, \@humanGenomeReferenceFileEndings, "GATKVariantEvalAll");
    &CheckBuildDownLoadPreRequisites(\%parameter, \%scriptParameter, \%supportedCosmidReference, "GATKVariantEvalAll");

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{$scriptParameter{'sampleIDs'}});$sampleIDCounter++) { 
	
	&GATKVariantEvalAll(\%parameter, \%scriptParameter, \%sampleInfo, \%infilesLaneNoEnding, $scriptParameter{'sampleIDs'}[$sampleIDCounter], $scriptParameter{'aligner'}, "BOTH", $scriptParameter{'familyID'}, "GATKVariantEvalAll");
    }
}

if ($scriptParameter{'pGATKVariantEvalExome'} > 0) {  #Run GATK VariantEval for exome variants. Done per sampleID

    $logger->info("[GATK VariantEval Exome]\n");
    
    &CheckBuildHumanGenomePreRequisites(\%parameter, \%scriptParameter, \%fileInfo, \@humanGenomeReferenceFileEndings, "GATKVariantEvalExome");
    &CheckBuildDownLoadPreRequisites(\%parameter, \%scriptParameter, \%supportedCosmidReference, "GATKVariantEvalExome");

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{$scriptParameter{'sampleIDs'}});$sampleIDCounter++) { 
	
	&GATKVariantEvalExome(\%parameter, \%scriptParameter, \%sampleInfo, \%infilesLaneNoEnding, $scriptParameter{'sampleIDs'}[$sampleIDCounter], $scriptParameter{'aligner'}, "BOTH", $scriptParameter{'familyID'}, "GATKVariantEvalExome");
    }
}

if ($scriptParameter{'pRankVariants'} > 0) {  #Run RankVariants. Done per family

    $logger->info("[RankVariants]\n");

    &RankVariants(\%parameter, \%scriptParameter, \%sampleInfo, $scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH", "RankVariants");
}

if ($scriptParameter{'pQCCollect'} > 0) {  #Run QCCollect. Done per family

    $logger->info("[QCCollect]\n");

    &QCCollect(\%parameter, \%scriptParameter, $scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH", "QCCollect");
}

if ($scriptParameter{'pRemoveRedundantFiles'} > 0) {  #Sbatch generation of removal of alignment files
    
    $logger->info("[Removal of alignment files]\n");

    &RemoveRedundantFiles(\%parameter, \%scriptParameter, \%sampleInfo, \%infilesLaneNoEnding, \%lane, $scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH", "RemoveRedundantFiles");	
}

if ( ($scriptParameter{'pAnalysisRunStatus'} == 1) && ($scriptParameter{'dryRunAll'} == 0) ) {

    $sampleInfo{ $scriptParameter{'familyID'} }{ $scriptParameter{'familyID'} }{'AnalysisRunStatus'} = "notFinished";  #Add analysis run status flag.
}

if ($scriptParameter{'pAnalysisRunStatus'} > 0) {

    $logger->info("[AnalysisRunStatus]\n");

    &AnalysisRunStatus(\%parameter, \%scriptParameter, \%sampleInfo, $scriptParameter{'familyID'}, $scriptParameter{'aligner'}, "BOTH", "AnalysisRunStatus");
}

#Write QC for programs used in analysis                                                                                                                         
if ($scriptParameter{'sampleInfoFile'} ne 0) {#Write SampleInfo to yaml file
    
    &WriteYAML(\%sampleInfo, \$scriptParameter{'sampleInfoFile'});  #Write QC for sampleinfo used in analysis
}


######################
####SubRoutines#######
######################

sub AnalysisRunStatus { 

##AnalysisRunStatus
    
##Function : Execute last in MAIN chain and sets analysis run status flag to finished.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $familyID, $aligner, $callType, $programName
##         : $parameterHashRef       => The parameter hash {REF}
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef      => Info on samples and family hash {REF}
##         : $familyID               => The familyID
##         : $aligner                => The aligner used in the analysis
##         : $callType               => The variant call type
##         : $programName            => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $familyID = $_[3];
    my $aligner = $_[4];
    my $callType = $_[5];
    my $programName = $_[6];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					   {'directoryID' => $familyID,
					    'programName' => $programName,
					    'programDirectory' => lc($programName),
					   });
    
    print $FILEHANDLE q?status="0"?, "\n";  #Set status flagg so that perl notFinished remains in sampleInfoFile

    ###Test all file that are supposed to exists as they are present in the sampleInfo file
    my @pathsArrayRef;

    ## Collects all programs file path(s) created by MIP located in %sampleInfo
    &CollectPathEntries(\%{$sampleInfoHashRef}, \@pathsArrayRef);
    ## Collects all programs outfile path(s) created by MIP as OutDirectory->value and outFile->value located in %sampleInfo.
    &CollectOutDataPathsEntries(\%{$sampleInfoHashRef}, \@pathsArrayRef);
    print $FILEHANDLE q?files=(?;  #Create bash array
    foreach my $path (@pathsArrayRef) {

	if (defined($path)) {  #First analysis and dry run will otherwise cause try to print uninitialized values
	    
	    print $FILEHANDLE q?"?.$path.q?" ?;  #Add to array
	}
    }
    print $FILEHANDLE ")", "\n";  #Close bash array
    print $FILEHANDLE q?for file in ${files[@]}?, "\n";  #loop over files
    print $FILEHANDLE "do ", "\n";  #for each element in array do
    print $FILEHANDLE "\t".q?if [ -s $file ]; then?, "\n";  #file exists and is larger than zero
    print $FILEHANDLE "\t\t".q?echo "Found file $file"?, "\n";  #Echo
    print $FILEHANDLE "\t".q?else?, "\n";
    print $FILEHANDLE "\t\t".q?echo "Could not find $file" >&2?, "\n";  #Redirect to STDERR
    print $FILEHANDLE "\t\t".q?status="1"?, "\n";  #Set status flagg so that perl notFinished remains in sampleInfoFile
    print $FILEHANDLE "\t".q?fi?, "\n";
    print $FILEHANDLE q?done ?, "\n\n";

    ## Test VariantEffectPredictor fork status. If VariantEffectPredictor is unable to fork it will prematurely end the analysis and we will lose variants. 
    if (defined(${$sampleInfoHashRef}{$familyID}{$familyID}{'program'}{"VariantEffectPredictor"}{'OutFile'})) {

	my $variantEffectPredictorFile = ${$sampleInfoHashRef}{$familyID}{$familyID}{'program'}{"VariantEffectPredictor"}{'OutDirectory'}."/".${$sampleInfoHashRef}{$familyID}{$familyID}{'program'}{"VariantEffectPredictor"}{'OutFile'};

	print $FILEHANDLE q?if grep -q "WARNING Unable to fork" ?;  #not output the matched text only return the exit status code
	print $FILEHANDLE $variantEffectPredictorFile.q?; then?, "\n";  #Infile
	print $FILEHANDLE "\t".q?status="1"?, "\n";  #Found pattern
	print $FILEHANDLE "\t".q?echo "VariantEffectorPredictor fork status=FAILED for file: ?.$variantEffectPredictorFile.q?" >&2?, "\n";  #Echo
	print $FILEHANDLE q?else?, "\n";  #Infile is clean
	print $FILEHANDLE "\t".q?echo "VariantEffectorPredictor fork status=PASSED for file: ?.$variantEffectPredictorFile.q?" >&2?, "\n";  #Echo
	print $FILEHANDLE q?fi?, "\n\n";
    }	
    print $FILEHANDLE q?if [ $status -ne 1 ]; then?, "\n";  #eval status flag
    print $FILEHANDLE "\t".q?perl -i -p -e 'if($_=~/AnalysisRunStatus\:/) { s/notFinished/finished/g }' ?.${$scriptParameterHashRef}{'sampleInfoFile'}.q? ?, "\n\n";  
    print $FILEHANDLE q?fi?, "\n";

    close($FILEHANDLE); 
    
    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {
	
	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'dependencies' => 2, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
    return;
}

sub RemoveRedundantFiles {

##RemoveRedundantFiles
    
##Function : Generates a sbatch script, which removes redundant files.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $infilesLaneNoEndingHashRef, $laneHashRef, $familyID, $aligner, $callType, $programName
##         : $parameterHashRef           => The parameter hash {REF}
##         : $scriptParameterHashRef     => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef          => Info on samples and family hash {REF}
##         : $infilesLaneNoEndingHashRef => The infile(s) without the ".ending" {REF}
##         : $laneHashRef                => The lane info hash {REF}
##         : $familyID                   => The familyID
##         : $aligner                    => The aligner used in the analysis
##         : $callType                   => The variant call type
##         : $programName                => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $infilesLaneNoEndingHashRef = $_[3];
    my $laneHashRef = $_[4];
    my $familyID = $_[5];
    my $aligner = $_[6];
    my $callType = $_[7];
    my $programName = $_[8];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					   {'directoryID' => $familyID,
					    'programName' => $programName,
					    'programDirectory' => $aligner,
					   });
    
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{${$scriptParameterHashRef}{'sampleIDs'}});$sampleIDCounter++) { 
	
	my $sampleID = ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter];
	my $inSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner;
	
###Single files
	
	for (my $infileCounter=0;$infileCounter < scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} });$infileCounter++) {  #MosaikBuild takes both reads at once
	    
	    my $infile = ${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter]; 
	    
##MosaikBuild	
	    if ( (${$scriptParameterHashRef}{'pMosaikBuild'} > 0) || (${$scriptParameterHashRef}{'aligner'} eq "mosaik") ) {
		
		print $FILEHANDLE "rm ";
		print $FILEHANDLE $inSampleDirectory."/".$infile.".dat", "\n\n";  #MosaikBuild
		
	    }
##MosaikAlign
	    if ( (${$scriptParameterHashRef}{'pMosaikAlign'} > 0) || (${$scriptParameterHashRef}{'aligner'} eq "mosaik") ) {
		
		print $FILEHANDLE "rm ";
		print $FILEHANDLE $inSampleDirectory."/".$infile.".stat", "\n\n";  #MosaikAlign Stats
		
		&CheckMostCompleteAndRemoveFile($FILEHANDLE, \${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'}, \($inSampleDirectory."/".$infile.".bam"), ".bam");	    
	    }
##Remove BWA files
	    if (${$scriptParameterHashRef}{'pBwaAln'} > 0) {
		
		print $FILEHANDLE "rm ";
		print $FILEHANDLE $inSampleDirectory."/".$infile.".sai", "\n\n";  #BWA_Aln
	    }
	    if (${$scriptParameterHashRef}{'pBwaSampe'} >0) {
		
		&CheckMostCompleteAndRemoveFile($FILEHANDLE, \${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'}, \($inSampleDirectory."/".$infile.".bam"), ".bam");
	    }  
	    if (${$scriptParameterHashRef}{'pBwaMem'} >0) {
		
		&CheckMostCompleteAndRemoveFile($FILEHANDLE, \${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'}, \($inSampleDirectory."/".$infile.".bam"), ".bam");
	    }    	    
##Sorted BAM
	    if (${$scriptParameterHashRef}{'pPicardToolsSortSam'} > 0) {
		
		my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pPicardToolsSortSam'}{'fileEnding'};
		
		&CheckMostCompleteAndRemoveFile($FILEHANDLE, \${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'}, \($inSampleDirectory."/".$infile.$outfileEnding.".bam"), ".bam");  #Sorted BAM and bai file
	    }
	}
	
	## Check if any files for this sampleID were merged previously to set infile and PicardToolsMergeSwitch to enable correct handling of number of infiles to process
	my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles(\%{$scriptParameterHashRef}, \%sampleInfo, \%lane, \%infilesLaneNoEnding, $sampleID);        
	
	if ($PicardToolsMergeSwitch == 1) {  #Files were merged previously
	    
	    if (${$scriptParameterHashRef}{'pPicardToolsMergeSamFiles'} > 0) {
		
		my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pPicardToolsMergeSamFiles'}{'fileEnding'};
		
		&CheckMostCompleteAndRemoveFile($FILEHANDLE, \${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'}, \($inSampleDirectory."/".$infile.$outfileEnding.".bam"), ".bam");  #merged BAM and bai file
	    }	
	    if (${$scriptParameterHashRef}{'pPicardToolsMarkduplicates'} > 0) {
		
		my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pPicardToolsMarkduplicates'}{'fileEnding'};
		
		&CheckMostCompleteAndRemoveFile($FILEHANDLE, \${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'}, \($inSampleDirectory."/".$infile.$outfileEnding.".bam"), ".bam");  #Dedupped BAM and bai file
	    }
	    if (${$scriptParameterHashRef}{'pGATKRealigner'} > 0) {
		
		my $inSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner."/GATK";   
		my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pGATKRealigner'}{'fileEnding'};
	   
		&CheckMostCompleteAndRemoveFile($FILEHANDLE, \${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'}, \($inSampleDirectory."/".$infile.$outfileEnding.".bam"), ".bam");  #ReAligned BAM and bai file
	    }
	    if (${$scriptParameterHashRef}{'pGATKBaseRecalibration'} > 0) {
		
		my $inSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner."/GATK";
		my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pGATKBaseRecalibration'}{'fileEnding'};

		&CheckMostCompleteAndRemoveFile($FILEHANDLE, \${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'}, \($inSampleDirectory."/".$infile.$outfileEnding.".bam"), ".bam");  #BaseRecalibrated BAM and bai file
	    }
	    if (${$scriptParameterHashRef}{'pGATKHaploTypeCaller'} > 0) {  #Always collapses all files even if there is only one
		
		my $lanes = join("",@{${$laneHashRef}{$sampleID}});  #Extract lanes
		my $inSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner."/GATK";
		my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pGATKHaploTypeCaller'}{'fileEnding'};
		
		&CheckMostCompleteAndRemoveFile($FILEHANDLE, \${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'MostCompleteVCF'}{'Path'}, \($inSampleDirectory."/".$infile.$outfileEnding.".vcf"), ".vcf");  #HaplotypeCaller gvcf file
	    }
	}
	else {
	    
	    for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} });$infileCounter++) {  #For all infiles per lane
		
		my $infile = ${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter];
		
		if (${$scriptParameterHashRef}{'pPicardToolsMarkduplicates'} > 0) {
		    
		    my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pPicardToolsMarkduplicates'}{'fileEnding'};

		    &CheckMostCompleteAndRemoveFile($FILEHANDLE, \${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'}, \($inSampleDirectory."/".$infile.$outfileEnding.".bam"), ".bam");  #Dedupped BAM and bai file
		}
		if (${$scriptParameterHashRef}{'pGATKRealigner'} > 0) {
		    
		    my $inSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner."/GATK";
		    my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pGATKRealigner'}{'fileEnding'};
		    
		    &CheckMostCompleteAndRemoveFile($FILEHANDLE, \${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'}, \($inSampleDirectory."/".$infile.$outfileEnding.".bam"), ".bam");  #ReAligned BAM and bai file
		}
		if (${$scriptParameterHashRef}{'pGATKBaseRecalibration'} > 0) {
		    
		    my $inSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner."/GATK";
		    my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pGATKBaseRecalibration'}{'fileEnding'};

		    &CheckMostCompleteAndRemoveFile($FILEHANDLE, \${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'}, \($inSampleDirectory."/".$infile.$outfileEnding.".bam"), ".bam");  #BaseRecalibrated BAM and bai file
		}
		if (${$scriptParameterHashRef}{'pGATKHaploTypeCaller'} > 0) {  #Always collapses all files even if there is only one
		    
		    my $lanes = join("",@{${$laneHashRef}{$sampleID}});  #Extract lanes
		    my $inSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner."/GATK";
		    my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pGATKHaploTypeCaller'}{'fileEnding'};
		    
		    &CheckMostCompleteAndRemoveFile($FILEHANDLE, \${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'MostCompleteVCF'}{'Path'}, \($inSampleDirectory."/".$sampleID."_lanes_".$lanes.$outfileEnding.".vcf"), ".vcf");  #HaplotypeCaller gvcf file
		}
	    }
	}
    }
###Family files
    if (${$scriptParameterHashRef}{'pGATKGenoTypeGVCFs'} > 0) {

	my $outFamilyDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID."/".$aligner."/GATK";  #New outfile directory
	my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'pGATKGenoTypeGVCFs'}{'fileEnding'};
	
	&CheckMostCompleteAndRemoveFile($FILEHANDLE, \${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'MostCompleteVCF'}{'Path'}, \($outFamilyDirectory."/".$familyID.$outfileEnding.$callType.".vcf"), ".vcf");  #GATKGenoTypeGVCFs vcf file
	
    }
    if (${$scriptParameterHashRef}{'pGATKVariantRecalibration'} > 0) {
	
	my $outFamilyDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID."/".$aligner."/GATK";  #New outfile directory
	my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'pGATKVariantRecalibration'}{'fileEnding'};
	
	&CheckMostCompleteAndRemoveFile($FILEHANDLE, \${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'MostCompleteVCF'}{'Path'}, \($outFamilyDirectory."/".$familyID.$outfileEnding.$callType.".vcf"), ".vcf");  #pGATKVariantRecalibration vcf file
    
    }
    if (${$scriptParameterHashRef}{'pAnnovar'} > 0) {
	
	my $outFamilyDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID."/".$aligner."/GATK";
	my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'pAnnovar'}{'fileEnding'};
	
	print $FILEHANDLE "rm ";
	print $FILEHANDLE $outFamilyDirectory."/".$familyID.$outfileEnding.$callType."*".${$scriptParameterHashRef}{'annovarGenomeBuildVersion'}."_*", "\n\n";  #Annovar data files
	print $FILEHANDLE "rm ";
	print $FILEHANDLE $outFamilyDirectory."/".$familyID.$outfileEnding.$callType.".*", "\n\n";  #Annovar data files
    }  
    close($FILEHANDLE);
}

sub SampleCheck { 

##SampleCheck
    
##Function : Tests sample for correct relatives (only performed for samples with relatives defined in pedigree file) performed on sequence data.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $familyID, $aligner, $callType, $programName
##         : $parameterHashRef       => The parameter hash {REF}
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef      => Info on samples and family hash {REF}
##         : $familyID               => The familyID
##         : $aligner                => The aligner used in the analysis
##         : $callType               => The variant call type
##         : $programName            => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $familyID = $_[3]; 
    my $aligner = $_[4];
    my $callType = $_[5];
    my $programName = $_[6];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					   {'directoryID' => $familyID,
					    'programName' => $programName,
					    'programDirectory' => lc($aligner."/".$programName),
					    'callType' => $callType,
					   });
    
    my $inFamilyDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID."/".$aligner."/gatk";
    my $outFamilyDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID."/".$aligner."/".lc($programName);
    my $infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'pGATKVariantRecalibration'}{'fileEnding'};

    print $FILEHANDLE "#Create Plink .ped and .map file per family using vcfTools","\n";
    print $FILEHANDLE "vcftools ";
    print $FILEHANDLE "--vcf ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf ";  #InFile
    print $FILEHANDLE "--plink ";  #PLINK format
    print $FILEHANDLE "--out ".$outFamilyDirectory."/".$familyID, "\n\n";  #OutFile (.ped and .map)

    print $FILEHANDLE "#Create vcfTools inbreeding coefficient F per family using vcfTools","\n";
    print $FILEHANDLE "vcftools ";
    print $FILEHANDLE "--vcf ".$inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf ";  #InFile
    print $FILEHANDLE "--het ";  #Individual inbreeding
    print $FILEHANDLE "--out ".$outFamilyDirectory."/".$familyID, "\n\n";  #Outfile

    if ( (${$scriptParameterHashRef}{'pSampleCheck'} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	## Collect QC metadata info for later use                                                                                               
	&SampleInfoQC(\%sampleInfo,
		       {'familyID' => ${$scriptParameterHashRef}{'familyID'},
			'programName' => "InbreedingFactor",
			'outDirectory' => $outFamilyDirectory,
			'outFileEnding' => $familyID.".het",
			'outDataType' => "infileDependent"
		       });
    }

    print $FILEHANDLE "#Create Plink .mibs per family","\n"; 
    print $FILEHANDLE "plink ";
    print $FILEHANDLE "--noweb ";  #No web check
    print $FILEHANDLE "--ped ".$outFamilyDirectory."/".$familyID.".ped ";  #InFile
    print $FILEHANDLE "--map ".$outFamilyDirectory."/".$familyID.".map ";  #InFile
    print $FILEHANDLE "--cluster ";  #Perform IBS clustering
    print $FILEHANDLE "--matrix ";  #Create a N x N matrix of genome-wide average IBS pairwise identities
    print $FILEHANDLE "--out ".$outFamilyDirectory."/".$familyID, "\n\n";  #OutFile

    if ( (${$scriptParameterHashRef}{'pSampleCheck'} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	## Collect QC metadata info for later use                                                                                               
	&SampleInfoQC(\%sampleInfo,
		       {'familyID' => ${$scriptParameterHashRef}{'familyID'},
			'programName' => "RelationCheck",
			'outDirectory' => $outFamilyDirectory,
			'outFileEnding' => $familyID.".mibs",
			'outDataType' => "infileDependent"
		       });
    }

    print $FILEHANDLE "#Create Plink sexcheck per family","\n"; 
    print $FILEHANDLE "plink ";
    print $FILEHANDLE "--noweb ";  #No web check
    print $FILEHANDLE "--ped ".$outFamilyDirectory."/".$familyID.".ped ";  #InFile
    print $FILEHANDLE "--map ".$outFamilyDirectory."/".$familyID.".map ";  #InFile
    print $FILEHANDLE "--check-sex ";  #uses X chromosome data to determine sex (i.e. based on heterozygosity rates) 
    print $FILEHANDLE "--out ".$outFamilyDirectory."/".$familyID, "\n\n";  #OutFile

    if ( (${$scriptParameterHashRef}{'pSampleCheck'} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	## Collect QC metadata info for later use                                                                                               
	&SampleInfoQC(\%sampleInfo,
		       {'familyID' => ${$scriptParameterHashRef}{'familyID'},
			'programName' => "SexCheck",
			'outDirectory' => $outFamilyDirectory,
			'outFileEnding' => $familyID.".sexcheck",
			'outDataType' => "infileDependent"
		       });
    }
    
    print $FILEHANDLE "wait", "\n\n";    
    
    close($FILEHANDLE); 

    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {
	
	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'dependencies' => 2, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}

sub QCCollect { 

##QCCollect
    
##Function : Collect qc metrics for this analysis run.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $familyID, $aligner, $callType, $programName
##         : $parameterHashRef       => The parameter hash {REF}
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $familyID               => The familyID
##         : $aligner                => The aligner used in the analysis
##         : $callType               => The variant call type
##         : $programName            => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $familyID = $_[2]; 
    my $aligner = $_[3];
    my $callType = $_[4];
    my $programName = $_[5];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					   {'directoryID' => $familyID,
					    'programName' => $programName,
					    'programDirectory' => lc($programName),
					   });
    
    my $infile = ${$scriptParameterHashRef}{'outDataDir'}."/".${$scriptParameterHashRef}{'familyID'}."/qc_sampleinfo.yaml";
    my $inFamilyDirectory =  ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID;
    my $outFamilyDirectory =  ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID;

    print $FILEHANDLE "perl ".${$scriptParameterHashRef}{'inScriptDir'}."/qcCollect.pl ";
    print $FILEHANDLE "-sampleInfoFile ".${$scriptParameterHashRef}{'QCCollectSampleInfoFile'}." ";
    print $FILEHANDLE "-regExpFile ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'QCCollectRegExpFile'}." ";
    print $FILEHANDLE "-o ".$outFamilyDirectory."/".$familyID."_qcmetrics.yaml ", "\n\n";     
    
    close($FILEHANDLE); 
    
    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	## Collect QC metadata info for later use
	&SampleInfoQC(\%sampleInfo,
		       {'familyID' => ${$scriptParameterHashRef}{'familyID'},
			'programName' => "QCCollect",
			'outDirectory' => $outFamilyDirectory,
			'outFileEnding' => $familyID."_qcmetrics.yaml",
			'outDataType' => "infileDependent"
		       });

	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'dependencies' => 1, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}

sub RankVariants { 

##RankVariants
    
##Function : Filter and Rank variants depending on mendelian inheritance, frequency and phenotype.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $familyID, $aligner, $callType, $programName
##         : $parameterHashRef       => The parameter hash {REF}
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef      => Info on samples and family hash {REF}
##         : $familyID               => The familyID
##         : $aligner                => The aligner used in the analysis
##         : $callType               => The variant call type
##         : $programName            => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $familyID = $_[3]; 
    my $aligner = $_[4];
    my $callType = $_[5];
    my $programName = $_[6];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					   {'directoryID' => $familyID,
					    'programName' => $programName,
					    'programDirectory' => lc($aligner."/gatk/".$programName),
					    'processTime' => 4,
					    'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
					   });

    ## Assign directories
    my $inFamilyDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID."/".$aligner."/gatk";
    my $outFamilyDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID."/".$aligner."/gatk/";

    my $infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'pAnnovar'}{'fileEnding'};
    my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{"p".$programName}{'fileEnding'};
    my $analysisType = "";

    ## Gene models and ranking  
    print $FILEHANDLE "workon ".${$scriptParameterHashRef}{'pythonVirtualEnvironment'}, "\n\n";  #Activate python environment
    
    for (my $VEPOutputFilesCounter=0;$VEPOutputFilesCounter<$VEPOutputFiles;$VEPOutputFilesCounter++) {
	
	if ($VEPOutputFilesCounter == 1) {
	    
	    $analysisType = ".selected";  #SelectFile variants
	}

	## Copy file(s) to temporary directory
	print $FILEHANDLE "## Copy file(s) to temporary directory\n";
	&MigrateFileToTemp($FILEHANDLE, 
			    {'path' => $inFamilyDirectory."/".$familyID.$infileEnding.$callType.$analysisType.".vcf*",
			     'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
			    });
	print $FILEHANDLE "wait", "\n\n";
	
	## Calculate Gene Models
	print $FILEHANDLE "## Calculate Gene Models", "\n";    
	print $FILEHANDLE "genmod annotate ";
	print $FILEHANDLE ${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$infileEnding.$callType.$analysisType.".vcf ";  #InFile
	print $FILEHANDLE "--family_file ".${$scriptParameterHashRef}{'pedigreeFile'}." ";  #Pedigree file
	
	if (${$scriptParameterHashRef}{'pVariantEffectPredictor'} > 0) {  #Use VEP annotations in compound models

	    print $FILEHANDLE "--vep "; 
	}
	else {
	    
	    print $FILEHANDLE ${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'geneFile'}." ";  #Gene file used for annotating AR_compounds
	}
	if (${$scriptParameterHashRef}{'instanceTag'} eq "CMMS") {
	    
	    print $FILEHANDLE "--family_type cmms ";  #CMMS flag
	}
	if (${$scriptParameterHashRef}{'caddWGSSNVs'} == 1) {
	 
	    print $FILEHANDLE "--cadd_file ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'caddWGSSNVsFile'}." ";  #Whole genome sequencing CADD score file
	}
	if (${$scriptParameterHashRef}{'cadd1000Genomes'} == 1) {
	 
	    print $FILEHANDLE "--cadd_1000g ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'cadd1000GenomesFile'}." ";  #1000G CADD score file
	}
	if (${$scriptParameterHashRef}{'wholeGene'} == 1) {
	 
	    print $FILEHANDLE "--whole_gene "; 
	}

	## Ranking
	print $FILEHANDLE "| ";  #Pipe
	print $FILEHANDLE "score_mip_variants ";
	print $FILEHANDLE "- ";  #Expect infile stream
	print $FILEHANDLE ${$scriptParameterHashRef}{'pedigreeFile'}." ";  #Pedigree file

	if (${$scriptParameterHashRef}{'rankModelFile'} ne "noUserInfo") {

	    print $FILEHANDLE "--plugin_file ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'rankModelFile'}." ";  #Rank model config.ini file 
	}
	print $FILEHANDLE "-o ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType.$analysisType.".vcf ", "\n\n";  #Outfile

	## Copies file from temporary directory.
	print $FILEHANDLE "## Copy file from temporary directory\n";
	&MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType.$analysisType.".vcf", $outFamilyDirectory."/", $FILEHANDLE);
	print $FILEHANDLE "wait", "\n\n";

	if ( (${$scriptParameterHashRef}{'pRankVariants'} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {
	    
	    if ($VEPOutputFilesCounter == 1) {
	
		${$sampleInfoHashRef}{$familyID}{$familyID}{'program'}{'RankVariants'}{'Clinical'}{'Path'} = $outFamilyDirectory."/".$familyID.$outfileEnding.$callType.$analysisType.".vcf";   #Save clinical candidate list path
	    }
	    else {

		${$sampleInfoHashRef}{$familyID}{$familyID}{'program'}{'RankVariants'}{'Research'}{'Path'} = $outFamilyDirectory."/".$familyID.$outfileEnding.$callType.$analysisType.".vcf";   #Save research candidate list path
	    }
	}
    }
    print $FILEHANDLE "deactivate ", "\n\n";  #Deactivate python environment

    &RemoveDirectory(\${$scriptParameterHashRef}{'tempDirectory'}, $FILEHANDLE);

    close($FILEHANDLE);   

    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {
	
	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'dependencies' => 1, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}


sub GATKVariantEvalExome { 

##GATKVariantEvalExome
    
##Function : GATK VariantEval for exome variants.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $infilesLaneNoEndingHashRef, $sampleID, $aligner, $callType, $familyID, $programName
##         : $parameterHashRef           => The parameter hash {REF}
##         : $scriptParameterHashRef     => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef          => Info on samples and family hash {REF}
##         : $infilesLaneNoEndingHashRef => The infile(s) without the ".ending" {REF}
##         : $sampleID                   => The sampleID
##         : $aligner                    => The aligner used in the analysis
##         : $callType                   => The variant call type
##         : $familyID                   => The familyID
##         : $programName                => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $infilesLaneNoEndingHashRef =$_[3];
    my $sampleID = $_[4]; 
    my $aligner = $_[5];
    my $callType = $_[6];
    my $familyID = $_[7];
    my $programName = $_[8];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					   {'directoryID' => $sampleID,
					    'programName' => $programName,
					    'programDirectory' => lc($aligner."/gatk/varianteval"),
					    'callType' => $callType,
					    'processTime' => 2,
					    'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
					   });

    ## Assign directories
    my $outSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner."/gatk/varianteval";
    my $inFamilyDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID."/".$aligner."/gatk";

    my $infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'pGATKVariantRecalibration'}{'fileEnding'};
    my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'pGATKVariantRecalibration'}{'fileEnding'};

    ## Check if any files for this sampleID were merged previously to set infile and PicardToolsMergeSwitch to enable correct handling of number of infiles to process
    my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles(\%{$scriptParameterHashRef}, \%sampleInfo, \%lane, \%infilesLaneNoEnding, $sampleID);

    ## Copy file(s) to temporary directory
    print $FILEHANDLE "## Copy file(s) to temporary directory\n";
    &MigrateFileToTemp($FILEHANDLE, 
			{'path' => $inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf*",
			 'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
			});

    print $FILEHANDLE "wait", "\n\n";

    if ($PicardToolsMergeSwitch == 1) {  #Files were merged previously

	## Select SampleID from familyID vcf file

	## GATK SelectVariants
	print $FILEHANDLE "## GATK SelectVariants","\n";

	## Writes java core commands to filehandle.
	&JavaCore($FILEHANDLE,
		  {'memoryAllocation' => "Xmx2g",
		   'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
		   'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		  });

	print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	print $FILEHANDLE "-l INFO ";  #Set the minimum level of logging
	print $FILEHANDLE "-T SelectVariants ";  #Type of analysis to run
	print $FILEHANDLE "-R ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference file
	print $FILEHANDLE "-V: ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$infileEnding.$callType.".vcf ";  #FamilyID inFile
	print $FILEHANDLE "-o ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.$callType."_temp.vcf ";  #SampleID exome outFile
	print $FILEHANDLE "-sn ".$sampleID, "\n\n";  #Include genotypes from this sample
 
	my $infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'pSnpEff'}{'fileEnding'};

	## Copy file(s) to temporary directory
	print $FILEHANDLE "## Copy file(s) to temporary directory\n";
	&MigrateFileToTemp($FILEHANDLE, 
		       {'path' => $inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf*",
			'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		       });
	print $FILEHANDLE "wait", "\n\n";

	## Extract exonic variants
	print $FILEHANDLE "## Extract exonic variants\n";
	print $FILEHANDLE q?perl -ne ' if ( ($_=~/exonic/) || ($_=/splicing/) ) {print $_;}' ?;
	print $FILEHANDLE ${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$infileEnding.$callType.".vcf ";  #InFile
	print $FILEHANDLE "> ".${$scriptParameterHashRef}{'tempDirectory'}."/".$sampleID.$infileEnding.$callType."_exonic_variants.vcf", "\n\n";  #OutFile

	## Include potential SelectFile variants
	if ($VEPOutputFiles == 2) {
	    
	    my $analysisType = ".selected";  #SelectFile variants

	    ## Copy file(s) to temporary directory
	    print $FILEHANDLE "## Copy file(s) to temporary directory\n";
	    &MigrateFileToTemp($FILEHANDLE, 
				{'path' => $inFamilyDirectory."/".$familyID.$infileEnding.$callType.$analysisType.".vcf*",
				 'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
				});
	    print $FILEHANDLE "wait", "\n\n";

	    ## Extract exonic variants
	    print $FILEHANDLE "## Extract exonic variants\n";
	    print $FILEHANDLE q?perl -ne ' if ( ($_=~/exonic/) || ($_=~/splicing/) ) {print $_;}' ?;
	    print $FILEHANDLE ${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$infileEnding.$callType.$analysisType.".vcf ";  #InFile
	    print $FILEHANDLE "> ".${$scriptParameterHashRef}{'tempDirectory'}."/".$sampleID.$infileEnding.$callType."_exonic_variants".$analysisType.".vcf", "\n\n";  #OutFile
	    
	    ## Merge orphans and selectfiles
	    print $FILEHANDLE "## Merge orphans and selectfile(s)\n";
	    print $FILEHANDLE "cat ";
	    print $FILEHANDLE ${$scriptParameterHashRef}{'tempDirectory'}."/".$sampleID.$infileEnding.$callType."_exonic_variants.vcf ";  #Orphan file
	    print $FILEHANDLE ${$scriptParameterHashRef}{'tempDirectory'}."/".$sampleID.$infileEnding.$callType."_exonic_variants".$analysisType.".vcf ";  #SelectFile variants
	    print $FILEHANDLE "> ".${$scriptParameterHashRef}{'tempDirectory'}."/".$sampleID.$infileEnding.$callType."_exonic_variants_combined.vcf", "\n\n";  #OutFile
	    
	    ## Sort combined file
	    print $FILEHANDLE "## Sort combined file\n";
	    print $FILEHANDLE "sort ";
	    print $FILEHANDLE "-k1,1 -k2,2n ";  #Numerically by chromosome and start position
	    print $FILEHANDLE ${$scriptParameterHashRef}{'tempDirectory'}."/".$sampleID.$infileEnding.$callType."_exonic_variants_combined.vcf ";
	    print $FILEHANDLE "> ".${$scriptParameterHashRef}{'tempDirectory'}."/".$sampleID.$infileEnding.$callType."_exonic_variants.vcf", "\n\n";  #OutFile
	}

	print $FILEHANDLE q?perl -ne ' if ($_=~/^#/) {print $_;}' ?;
	print $FILEHANDLE ${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$infileEnding.$callType.".vcf ";  #InFile
	print $FILEHANDLE "| ";  #Pipe
	print $FILEHANDLE "cat ";
	print $FILEHANDLE "- ";
	print $FILEHANDLE ${$scriptParameterHashRef}{'tempDirectory'}."/".$sampleID.$infileEnding.$callType."_exonic_variants.vcf ";  #inFile
	print $FILEHANDLE "> ".${$scriptParameterHashRef}{'tempDirectory'}."/".$sampleID.$infileEnding.$callType."_exonic_variants_head.vcf", "\n\n";  #OutFile

	## Intersect exonic variants from created sampleID vcf file (required for GATKVariantEval for exonic variants)
	print $FILEHANDLE "## Intersect exonic variants from created sampleID vcf file\n";
	print $FILEHANDLE "intersectBed ";
	print $FILEHANDLE "-header ";  #Print the header from the A file prior to results.
	print $FILEHANDLE "-a ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.$callType."_temp.vcf ";  #SampleID temp exome vcf inFile
	print $FILEHANDLE "-b ".${$scriptParameterHashRef}{'tempDirectory'}."/".$sampleID.$infileEnding.$callType."_exonic_variants_head.vcf ";  #SampleID exonic variants
	print $FILEHANDLE "> ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.$callType."_exome.vcf", "\n\n";  #OutFile (VCF-format)

	## VariantEval
	print $FILEHANDLE "## GATK VariantEval","\n";
	
	## Writes java core commands to filehandle.
	&JavaCore($FILEHANDLE,
		  {'memoryAllocation' => "Xmx2g",
		   'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
		   'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		  });
	
	print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	print $FILEHANDLE "-l INFO ";  #Set the minimum level of logging
	print $FILEHANDLE "-T VariantEval ";  #Type of analysis to run
	print $FILEHANDLE "-R ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference file
	print $FILEHANDLE "-D ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'GATKVariantEvalDbSNP'}." ";  #dbSNP file
	print $FILEHANDLE "-gold ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'GATKVariantEvalGold'}." ";  #Evaluations that count calls at sites of true variation (e.g., indel calls) will use this argument as their gold standard for comparison
	print $FILEHANDLE "--eval ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.$callType."_exome.vcf ";  #InFile
	print $FILEHANDLE "-o ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.$callType."_exome.vcf.varianteval", "\n\n";  #OutFile

	## Copies file from temporary directory.
	print $FILEHANDLE "## Copy file from temporary directory\n";
	&MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.$callType."_exome.vcf.varianteval", $outSampleDirectory."/", $FILEHANDLE);
	print $FILEHANDLE "wait", "\n\n";

	if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	    ## Collect QC metadata info for later use                                                                                 
	    &SampleInfoQC(\%sampleInfo,
			   {'familyID' => ${$scriptParameterHashRef}{'familyID'},
			    'sampleID' => $sampleID,
			    'programName' => "VariantEval_Exome",
			    'infile' => $infile,
			    'outDirectory' => $outSampleDirectory,
			    'outFileEnding' => $outfileEnding.$callType."_exome.vcf.varianteval",
			    'outDataType' => "infileDependent"
			   });
	}   
    }
    else {  #No previous merge

	## Select SampleID from familyID vrecal vcf file
	
	for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} });$infileCounter++) {  #For all infiles per lane
	    
	    my $infile = ${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter];
	    
	    ## GATK SelectVariants
	    print $FILEHANDLE "## GATK SelectVariants","\n";

	    ## Writes java core commands to filehandle.
	    &JavaCore($FILEHANDLE,
		      {'memoryAllocation' => "Xmx2g",
		       'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
		       'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		      });
	    
	    print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	    print $FILEHANDLE "-l INFO ";  #Set the minimum level of logging
	    print $FILEHANDLE "-T SelectVariants ";  #Type of analysis to run
	    print $FILEHANDLE "-R ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference file
	    print $FILEHANDLE "-V: ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$infileEnding.$callType.".vcf ";  #FamilyID infile 
	    print $FILEHANDLE "-o ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.$callType."_temp.vcf ";  #SampleID outFile
	    print $FILEHANDLE "-sn ".$sampleID, "\n\n";  #Include genotypes from this sample

	    my $infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'pSnpEff'}{'fileEnding'};

	    ## Copy file(s) to temporary directory
	    print $FILEHANDLE "## Copy file(s) to temporary directory\n";
	    &MigrateFileToTemp($FILEHANDLE, 
		       {'path' => $inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf*",
			'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		       });
	    print $FILEHANDLE "wait", "\n\n";

	    ## Extract exonic variants
	    print $FILEHANDLE "## Extract exonic variants\n";
	    print $FILEHANDLE q?perl -ne ' if ( ($_=~/exonic/) || ($_=/splicing/) ) {print $_;}' ?;
	    print $FILEHANDLE ${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$infileEnding.$callType.".vcf ";  #InFile
	    print $FILEHANDLE "> ".${$scriptParameterHashRef}{'tempDirectory'}."/".$sampleID.$infileEnding.$callType."_exonic_variants.vcf", "\n\n";  #OutFile

	    ## Include potential SelectFile variants
	    if ($VEPOutputFiles == 2) {

		my $analysisType = ".selected";  #SelectFile variants

		## Copy file(s) to temporary directory
		print $FILEHANDLE "## Copy file(s) to temporary directory\n";
		&MigrateFileToTemp($FILEHANDLE, 
				    {'path' => $inFamilyDirectory."/".$familyID.$infileEnding.$callType.$analysisType.".vcf*",
				     'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
				    });
		print $FILEHANDLE "wait", "\n\n";

		## Extract exonic variants
		print $FILEHANDLE "## Extract exonic variants\n";
		print $FILEHANDLE q?perl -ne ' if ( ($_=~/exonic/) || ($_=/splicing/) ) {print $_;}' ?;
		print $FILEHANDLE ${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$infileEnding.$callType.$analysisType.".vcf ";  #InFile
		print $FILEHANDLE "> ".${$scriptParameterHashRef}{'tempDirectory'}."/".$sampleID.$infileEnding.$callType."_exonic_variants".$analysisType.".vcf", "\n\n";  #OutFile
		
		## Merge orphans and selectfile(s)
		print $FILEHANDLE "## Merge orphans and selectfile(s)\n";
		print $FILEHANDLE "cat ";
		print $FILEHANDLE ${$scriptParameterHashRef}{'tempDirectory'}."/".$sampleID.$infileEnding.$callType."_exonic_variants.vcf ";  #Orphan file
		print $FILEHANDLE ${$scriptParameterHashRef}{'tempDirectory'}."/".$sampleID.$infileEnding.$callType."_exonic_variants".$analysisType.".vcf ";  #SelectFile variants
		print $FILEHANDLE "> ".${$scriptParameterHashRef}{'tempDirectory'}."/".$sampleID.$infileEnding.$callType."_exonic_variants_combined.vcf", "\n\n";  #OutFile

		## Sort combined file
		print $FILEHANDLE "## Sort combined file\n";
		print $FILEHANDLE "sort ";
		print $FILEHANDLE "-k1,1 -k2,2n ";  #Numerically by chromosome and start position
		print $FILEHANDLE ${$scriptParameterHashRef}{'tempDirectory'}."/".$sampleID.$infileEnding.$callType."_exonic_variants_combined.vcf ";
		print $FILEHANDLE "> ".${$scriptParameterHashRef}{'tempDirectory'}."/".$sampleID.$infileEnding.$callType."_exonic_variants.vcf", "\n\n";  #OutFile
	    }

	    print $FILEHANDLE q?perl -ne ' if ($_=~/^#/) {print $_;}' ?;
	    print $FILEHANDLE ${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$infileEnding.$callType.".vcf ";  #InFile
	    print $FILEHANDLE "| ";  #Pipe
	    print $FILEHANDLE "cat ";
	    print $FILEHANDLE "- ";
	    print $FILEHANDLE ${$scriptParameterHashRef}{'tempDirectory'}."/".$sampleID.$infileEnding.$callType."_exonic_variants.vcf ";  #inFile
	    print $FILEHANDLE "> ".${$scriptParameterHashRef}{'tempDirectory'}."/".$sampleID.$infileEnding.$callType."_exonic_variants_head.vcf", "\n\n";  #OutFile
	    
	    ## Intersect exonic variants from created sampleID vcf file (required for GATKVariantEval for exonic variants)
	    print $FILEHANDLE "## Intersect exonic variants from created sampleID vcf file\n";
	    print $FILEHANDLE "intersectBed ";
	    print $FILEHANDLE "-header ";  #Print the header from the A file prior to results.
	    print $FILEHANDLE "-a ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.$callType."_temp.vcf ";  #SampleID temp exome vcf inFile
	    print $FILEHANDLE "-b ".${$scriptParameterHashRef}{'tempDirectory'}."/".$sampleID.$infileEnding.$callType."_exonic_variants_head.vcf ";  #SampleID exonic variants
	    print $FILEHANDLE "> ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.$callType."_exome.vcf", "\n\n";  #OutFile (VCF-format)
	    
	    ## VariantEval
	    print $FILEHANDLE "## GATK VariantEval","\n";
	    
	    ## Writes java core commands to filehandle.
	    &JavaCore($FILEHANDLE,
		      {'memoryAllocation' => "Xmx2g",
		       'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
		       'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		      });
	    
	    print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	    print $FILEHANDLE "-l INFO ";  #Set the minimum level of logging
	    print $FILEHANDLE "-T VariantEval ";  #Type of analysis to run
	    print $FILEHANDLE "-R ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference file
	    print $FILEHANDLE "-D ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'GATKVariantEvalDbSNP'}." ";  #dbSNP file
	    print $FILEHANDLE "-gold ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'GATKVariantEvalGold'}." ";  #Evaluations that count calls at sites of true variation (e.g., indel calls) will use this argument as their gold standard for comparison
	    print $FILEHANDLE "--eval ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.$callType."_exome.vcf ";  #InFile
	    print $FILEHANDLE "-o ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.$callType."_exome.vcf.varianteval", "\n\n";  #OutFile

	     ## Copies file from temporary directory.
	    print $FILEHANDLE "## Copy file from temporary directory\n";
	    &MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.$callType."_exome.vcf.varianteval", $outSampleDirectory."/", $FILEHANDLE);
	    print $FILEHANDLE "wait", "\n\n";

	    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

		## Collect QC metadata info for later use                                                                                    
		&SampleInfoQC(\%sampleInfo,
			       {'familyID' => ${$scriptParameterHashRef}{'familyID'},
				'sampleID' => $sampleID,
				'programName' => "VariantEval_Exome",
				'infile' => $infile,
				'outDirectory' => $outSampleDirectory,
				'outFileEnding' => $outfileEnding.$callType."_exome.vcf.varianteval",
				'outDataType' => "infileDependent"
			       });
	    }
	}
    } 

    &RemoveDirectory(\${$scriptParameterHashRef}{'tempDirectory'}, $FILEHANDLE);

    close($FILEHANDLE);   
 
    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {
	
	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'dependencies' => 2, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}


sub GATKVariantEvalAll { 

##GATKVariantEvalAll
    
##Function : GATK VariantEval for all variants.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $infilesLaneNoEndingHashRef, $sampleID, $aligner, $callType, $familyID, $programName
##         : $parameterHashRef           => The parameter hash {REF}
##         : $scriptParameterHashRef     => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef          => Info on samples and family hash {REF}
##         : $infilesLaneNoEndingHashRef => The infile(s) without the ".ending" {REF}
##         : $sampleID                   => The sampleID
##         : $aligner                    => The aligner used in the analysis
##         : $callType                   => The variant call type
##         : $familyID                   => The familyID
##         : $programName                => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $infilesLaneNoEndingHashRef = $_[3];
    my $sampleID = $_[4]; 
    my $aligner = $_[5];
    my $callType = $_[6]; 
    my $familyID = $_[7];
    my $programName =$_[8];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header    
    my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					    {'directoryID' => $sampleID,
					     'programName' => $programName,
					     'programDirectory' => lc($aligner."/gatk/varianteval"),
					     'callType' => $callType,
					     'processTime' => 2,
					     'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
					    });

    ## Assign directories
    my $inSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner."/gatk";
    my $outSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner."/gatk/varianteval";
    my $inFamilyDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID."/".$aligner."/gatk";

    my $infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'pGATKVariantRecalibration'}{'fileEnding'};
    my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'pGATKVariantRecalibration'}{'fileEnding'};

    ## Check if any files for this sampleID were merged previously to set infile and PicardToolsMergeSwitch to enable correct handling of number of infiles to process
    my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles(\%{$scriptParameterHashRef}, \%sampleInfo, \%lane, \%infilesLaneNoEnding, $sampleID);
    
    ## Copy file(s) to temporary directory
    print $FILEHANDLE "## Copy file(s) to temporary directory\n"; 
    &MigrateFileToTemp($FILEHANDLE, 
			{'path' => $inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf*",
			 'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
			});
    print $FILEHANDLE "wait", "\n\n";

    if ($PicardToolsMergeSwitch == 1) {  #Files were merged previously

	## Select SampleID from familyID vcf file

	## GATK SelectVariants
	print $FILEHANDLE "## GATK SelectVariants","\n";

	## Writes java core commands to filehandle.
	&JavaCore($FILEHANDLE,
		  {'memoryAllocation' => "Xmx2g",
		   'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
		   'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		  });
	
	print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	print $FILEHANDLE "-l INFO ";  #Set the minimum level of logging
	print $FILEHANDLE "-T SelectVariants ";  #Type of analysis to run
	print $FILEHANDLE "-R ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference file
	print $FILEHANDLE "-V: ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$infileEnding.$callType.".vcf ";  #FamilyID inFile
	print $FILEHANDLE "-o ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.$callType.".vcf ";  #SampleID outFile
	print $FILEHANDLE "-sn ".$sampleID, "\n\n";  #Include genotypes from this sample

	## GATK VariantEval
	print $FILEHANDLE "## GATK VariantEval","\n";
	
	## Writes java core commands to filehandle.
	&JavaCore($FILEHANDLE,
		  {'memoryAllocation' => "Xmx2g",
		   'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
		   'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		  });
	
	print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	print $FILEHANDLE "-l INFO ";  #Set the minimum level of logging
	print $FILEHANDLE "-T VariantEval ";  #Type of analysis to run
	print $FILEHANDLE "-R ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference file
	print $FILEHANDLE "-D ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'GATKVariantEvalDbSNP'}." ";  #dbSNP file
	print $FILEHANDLE "-gold ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'GATKVariantEvalGold'}." ";  #Evaluations that count calls at sites of true variation (e.g., indel calls) will use this argument as their gold standard for comparison
	print $FILEHANDLE "--eval ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$infileEnding.$callType.".vcf ";  #InFile
	print $FILEHANDLE "-o ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.$callType.".vcf.varianteval", "\n\n";  #OutFile

	## Copies file from temporary directory.
	print $FILEHANDLE "## Copy file from temporary directory\n";
	&MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.$callType.".vcf.varianteval", $outSampleDirectory."/", $FILEHANDLE);
	print $FILEHANDLE "wait", "\n\n";

	if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	    ## Collect QC metadata info for later use                                                                                                
	    &SampleInfoQC(\%sampleInfo,
			   {'familyID' => ${$scriptParameterHashRef}{'familyID'},
			    'sampleID' => $sampleID,
			    'programName' => "VariantEval_All",
			    'infile' => $infile,
			    'outDirectory' => $outSampleDirectory,
			    'outFileEnding' => $outfileEnding.$callType.".vcf.varianteval",
			    'outDataType' => "infileDependent"
			   });
	}
    }   
    else {  #No previous merge

	## Select SampleID from familyID vrecal vcf file
	for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} });$infileCounter++) {  #For all infiles per lane
	    
	    my $infile = ${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter];

	    ## GATK SelectVariants
	    print $FILEHANDLE "## GATK SelectVariants","\n";

	    ## Writes java core commands to filehandle.
	    &JavaCore($FILEHANDLE,
		      {'memoryAllocation' => "Xmx2g",
		       'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
		       'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		      });

	    print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	    print $FILEHANDLE "-l INFO ";  #Set the minimum level of logging
	    print $FILEHANDLE "-T SelectVariants ";  #Type of analysis to run
	    print $FILEHANDLE "-R ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference file
	    print $FILEHANDLE "-V: ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$infileEnding.$callType.".vcf ";  #FamilyID infile 
	    print $FILEHANDLE "-o ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.$callType.".vcf ";  #SampleID outFile
	    print $FILEHANDLE "-sn ".$sampleID, "\n\n";  #Include genotypes from this sample

	    ## GATK VariantEval
	    print $FILEHANDLE "## GATK VariantEval","\n";
	    
	    ## Writes java core commands to filehandle.
	    &JavaCore($FILEHANDLE,
		      {'memoryAllocation' => "Xmx2g",
		       'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
		       'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		      });
	    
	    print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	    print $FILEHANDLE "-l INFO ";  #Set the minimum level of logging
	    print $FILEHANDLE "-T VariantEval ";  #Type of analysis to run
	    print $FILEHANDLE "-R ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference file
	    print $FILEHANDLE "-D ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'GATKVariantEvalDbSNP'}." ";  #dbSNP file
	    print $FILEHANDLE "-gold ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'GATKVariantEvalGold'}." ";  #Evaluations that count calls at sites of true variation (e.g., indel calls) will use this argument as their gold standard for comparison
	    print $FILEHANDLE "--eval ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$infileEnding.$callType.".vcf ";  #InFile
	    print $FILEHANDLE "-o ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.$callType.".vcf.varianteval", "\n\n";  #OutFile

	    ## Copies file from temporary directory.
	    print $FILEHANDLE "## Copy file from temporary directory\n";
	    &MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.$callType.".vcf.varianteval", $outSampleDirectory."/", $FILEHANDLE);
	    print $FILEHANDLE "wait", "\n\n";

	    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

		## Collect QC metadata info for later use                                                                             
		&SampleInfoQC(\%sampleInfo,
			   {'familyID' => ${$scriptParameterHashRef}{'familyID'},
			    'sampleID' => $sampleID,
			    'programName' => "VariantEval_All",
			    'infile' => $infile,
			    'outDirectory' => $outSampleDirectory,
			    'outFileEnding' => $outfileEnding.$callType.".vcf.varianteval",
			    'outDataType' => "infileDependent"
			   });
	    }
	} 
    }

    &RemoveDirectory(\${$scriptParameterHashRef}{'tempDirectory'}, $FILEHANDLE);

    close($FILEHANDLE);   

    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'dependencies' => 2, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}


sub Annovar { 

##Annovar
    
##Function : Annotate and filter SNVs by gene, region and databases.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $annovarTableHashRef, $familyID, $aligner, $callType, $programName
##         : $parameterHashRef       => The parameter hash {REF}
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $annovarTableHashRef    => annovarTableHashRef {REF}
##         : $sampleInfoHashRef      => Info on samples and family hash {REF}
##         : $familyID               => The familyID
##         : $aligner                => The aligner used in the analysis
##         : $callType               => The variant call type
##         : $programName            => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $annovarTableHashRef = $_[3];
    my $familyID = $_[4];
    my $aligner = $_[5];
    my $callType = $_[6];
    my $programName = $_[7];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle

    ## Set the number of cores to allocate per sbatch job.
    my $nrCores = &NrofCoresPerSbatch(\%{$scriptParameterHashRef}, scalar(@{${$scriptParameterHashRef}{'annovarTableNames'}}));  #Detect the number of cores to use from @annovarTableNames. 
    
    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header    
    my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					   {'directoryID' => $familyID,
					    'programName' => $programName,
					    'programDirectory' => lc($aligner."/gatk"),
					    'callType' => $callType,
					    'nrofCores' => $nrCores,
					    'processTime' => 7,
					    'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
					   });

    ## Assign directories
    my $inFamilyDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID."/".$aligner."/gatk";
    my $outFamilyDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID."/".$aligner."/gatk";

    my $infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'pSnpEff'}{'fileEnding'};
    my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{"p".$programName}{'fileEnding'};
    my $analysisType = "";

    for (my $VEPOutputFilesCounter=0;$VEPOutputFilesCounter<$VEPOutputFiles;$VEPOutputFilesCounter++) {

	if ($VEPOutputFilesCounter == 1) {

	    $analysisType = ".selected";  #SelectFile variants
	}
	my $coreCounter=1;   	    	

	## Copy file(s) to temporary directory
	print $FILEHANDLE "## Copy file(s) to temporary directory\n"; 
	&MigrateFileToTemp($FILEHANDLE, 
			    {'path' => $inFamilyDirectory."/".$familyID.$infileEnding.$callType.$analysisType.".vcf*",
			     'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
			    });
	print $FILEHANDLE "wait", "\n\n";

	## Annovar
	print $FILEHANDLE "## Annovar","\n";
	print $FILEHANDLE "perl ".${$scriptParameterHashRef}{'annovarPath'}."/table_annovar.pl ";  #Annovar script 
	print $FILEHANDLE ${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$infileEnding.$callType.$analysisType.".vcf ";  #Infile
	print $FILEHANDLE ${$scriptParameterHashRef}{'annovarPath'}."/humandb ";  #annovar/humandb directory is assumed
	print $FILEHANDLE "-buildver ".${$scriptParameterHashRef}{'annovarGenomeBuildVersion'}." ";  #Genome build version
	print $FILEHANDLE "-vcfinput ";  #Input format
	print $FILEHANDLE "--remove ";  #Remove all temporary files
	print $FILEHANDLE "-out ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType.$analysisType.".vcf ";  #Outfile prefix
	print $FILEHANDLE "-protocol ";  #Comma-delimited string specifying database protocol

	print $FILEHANDLE join(',', @{${$scriptParameterHashRef}{'annovarTableNames'}})." ";  #Databases to use

	print $FILEHANDLE "-operation ";  #Comma-delimited string specifying type of operation	
	for (my $tableNamesCounter=0;$tableNamesCounter<scalar(@{${$scriptParameterHashRef}{'annovarTableNames'}});$tableNamesCounter++) {  #For all specified table names

	    if (${$annovarTableHashRef}{${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter]}{'annotation'} eq "geneanno") {
	
		print $FILEHANDLE "g";
	    }
	    elsif (${$annovarTableHashRef}{${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter]}{'annotation'} eq "regionanno") {
		
		print $FILEHANDLE "r";
	    }
	    elsif (${$annovarTableHashRef}{${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter]}{'annotation'} eq "filter") {
		
		print $FILEHANDLE "f";
	    }
	    unless ($tableNamesCounter == scalar(@{${$scriptParameterHashRef}{'annovarTableNames'}}) - 1) {
		
		print $FILEHANDLE ",";
	    }
	}
	print $FILEHANDLE " ";
	
	print $FILEHANDLE "-argument ";
	for (my $tableNamesCounter=0;$tableNamesCounter<scalar(@{${$scriptParameterHashRef}{'annovarTableNames'}});$tableNamesCounter++) {  #For all specified table names
	    
	    if (${$annovarTableHashRef}{${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter]}{'annotation'} eq "geneanno" ) {  #Use hgvs output style
		
		print $FILEHANDLE "'--hgvs ";  #Use hgvs annotation
		print $FILEHANDLE "--exonicsplicing'";  #Annotate variants near intron/exonic borders
	    }
	    if (${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter] =~/^1000g/) {#Set MAF TH

		print $FILEHANDLE "'--maf_threshold ".${$scriptParameterHashRef}{'annovarMAFThreshold'}."'";
	    }
	    if ( (${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter] =~/^snp/) || (${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter] =~/_esp/) ) {#Set MAF TH

		print $FILEHANDLE "'--score_threshold ".${$scriptParameterHashRef}{'annovarMAFThreshold'}."'"; #score_threshold since Annovar reserved the maf_threshold for 1000G
	    }
	    unless ($tableNamesCounter == scalar(@{${$scriptParameterHashRef}{'annovarTableNames'}}) - 1) {

		print $FILEHANDLE ",";
	    }
	}
	print $FILEHANDLE "\n\n";

	for (my $tableNamesCounter=0;$tableNamesCounter<scalar(@{${$scriptParameterHashRef}{'annovarTableNames'}});$tableNamesCounter++) {  #For all specified table names

	    if (${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter] =~/ensGene|refGene/) {  #Extra round to catch MT for refSeq as well

		## Annovar MT annotation
		print $FILEHANDLE "## Annovar MT annotation","\n";
		print $FILEHANDLE "perl ".${$scriptParameterHashRef}{'annovarPath'}."/table_annovar.pl ";  #Annovar script 
		print $FILEHANDLE ${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$infileEnding.$callType.$analysisType.".vcf ";  #Infile
		print $FILEHANDLE ${$scriptParameterHashRef}{'annovarPath'}."/humandb ";  #annovar/humandb directory is assumed
		print $FILEHANDLE "-vcfinput ";  #Input format
		print $FILEHANDLE "--remove ";  #Remove all temporary files
		print $FILEHANDLE "-buildver GRCh37_MT ";
		print $FILEHANDLE "-protocol ensGene ";  #Comma-delimited string specifying database protocol. NOTE: RefSeq does not have mitochondria gene definition. So ANNOVAR use either UCSC Known Gene or Ensembl Gene.
		print $FILEHANDLE "-operation g ";  #Comma-delimited string specifying type of operation	
		print $FILEHANDLE "-argument -exonicsplicing ";  #Annotate variants near intron/exonic borders 
		print $FILEHANDLE "-out ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType.$analysisType,"\n\n";  #OutFile prefix
	    }
	}
	print $FILEHANDLE "wait", "\n\n";

	##Merge vcf files to 1 
	my @mitochondria = ("vcf.".${$scriptParameterHashRef}{'annovarGenomeBuildVersion'}."_multianno", "GRCh37_MT_multianno");

	
	&CombineVariants(\%{$scriptParameterHashRef}, $FILEHANDLE, \@mitochondria, ${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType.$analysisType.".", ".vcf", ${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType.$analysisType.".vcf");

	## Copies file from temporary directory.
	print $FILEHANDLE "## Copy file from temporary directory\n";
	&MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType.$analysisType.".vcf", $outFamilyDirectory."/", $FILEHANDLE);
	&MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType.$analysisType.".vcf.idx", $outFamilyDirectory."/", $FILEHANDLE);
	print $FILEHANDLE "wait", "\n\n";
    }

    &RemoveDirectory(\${$scriptParameterHashRef}{'tempDirectory'}, $FILEHANDLE);

    close($FILEHANDLE);

    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'dependencies' => 1, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}


sub GATKReadBackedPhasing {
 
##GATKReadBackedPhasing
    
##Function : GATK ReadBackedPhasing performs physical phasing of SNP calls, based on sequencing reads.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $infilesLaneNoEndingHashRef, $familyID, $aligner, $callType, $programName
##         : $parameterHashRef           => The parameter hash {REF}
##         : $scriptParameterHashRef     => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef          => Info on samples and family hash {REF}
##         : $infilesLaneNoEndingHashRef => The infile(s) without the ".ending" {REF}
##         : $familyID                   => The familyID
##         : $aligner                    => The aligner used in the analysis
##         : $callType                   => The variant call type
##         : $programName                => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $infilesLaneNoEndingHashRef = $_[3];
    my $familyID = $_[4];
    my $aligner = $_[5];
    my $callType = $_[6];
    my $programName = $_[7];
    
    my $FILEHANDLE = IO::Handle->new();  #Create anonymous filehandle
    my $nrCores = 1;

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header    
    my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					   {'directoryID' => $familyID,
					    'programName' => $programName,
					    'programDirectory' => lc($aligner."/gatk"),
					    'callType' => $callType,
					    'nrofCores' => $nrCores,
					    'processTime' => 3,
					    'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
					   });

    ## Assign directories
    my $inFamilyDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID."/".$aligner."/gatk";
    my $outFamilyDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID."/".$aligner."/gatk";

    my $infileEnding;

    ## Choose infile ending depending on GATK PhasebyTransmission swith
    if (${$scriptParameterHashRef}{'pGATKPhaseByTransmission'} > 0) { 

	$infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'pGATKPhaseByTransmission'}{'fileEnding'};
    }
    else {

	$infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'pGATKVariantRecalibration'}{'fileEnding'};
    }
    my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{"p".$programName}{'fileEnding'};

    ## Copy VCF file(s) to temporary directory
    print $FILEHANDLE "## Copy file(s) to temporary directory\n"; 
    &MigrateFileToTemp($FILEHANDLE, 
			{'path' => $inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf*",
			 'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
			});
    print $FILEHANDLE "wait", "\n\n";

    ## Copy BAM file(s) to temporary directory
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{${$scriptParameterHashRef}{'sampleIDs'}});$sampleIDCounter++) {  #Collect infiles for all sampleIDs

	my $inSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter]."/".$aligner."/gatk";
	my $infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }{'pGATKBaseRecalibration'}{'fileEnding'};

	## Check if any files for this sampleID were merged previously to set infile and PicardToolsMergeSwitch to enable correct handling of number of infiles to process
	my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles(\%{$scriptParameterHashRef}, \%sampleInfo, \%lane, \%infilesLaneNoEnding, ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter]);

	if ($PicardToolsMergeSwitch == 1) {  #Alignment BAM-files merged previously
	    
	    ## Copy file(s) to temporary directory
	    print $FILEHANDLE "## Copy file(s) to temporary directory\n";
	    &MigrateFileToTemp($FILEHANDLE, 
				{'path' => $inSampleDirectory."/".$infile.$infileEnding.".b*",
				 'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
				});
	    print $FILEHANDLE "wait", "\n\n";
	}
	else {  #No previous merge of alignment BAM-files
	    
	    ## Copies files from source to temporary folder. Loop over files specified by $arrayRef and collects files from $extractArrayRef.
	    &MigrateFilesToTemp(\%{$scriptParameterHashRef}, \@{ ${$infilesLaneNoEndingHashRef}{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] } }, \@{ ${$infilesLaneNoEndingHashRef}{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] } }, $FILEHANDLE,
				{'inSampleDirectory' => $inSampleDirectory,
				 'nrCores' => $nrCores,
				 'fileEnding' => $infileEnding.".b*"
				});
	}
    }

    ## GATK ReadBackedPhasing
    print $FILEHANDLE "## GATK ReadBackedPhasing","\n";

    ## Writes java core commands to filehandle.
    &JavaCore($FILEHANDLE,
	      {'memoryAllocation' => "Xmx4g",
	       'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
	       'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
	      });
    
    print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
    print $FILEHANDLE "-l INFO ";  #Set the minimum level of logging
    print $FILEHANDLE "-T ReadBackedPhasing ";  #Type of analysis to run
    print $FILEHANDLE "-R ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference file
    print $FILEHANDLE "--phaseQualityThresh ".${$scriptParameterHashRef}{'GATKReadBackedPhasingPhaseQualityThreshold'}." ";

    if (${$scriptParameterHashRef}{'pGATKPhaseByTransmission'} > 0) { 

	print $FILEHANDLE "-respectPhaseInInput ";  #Already phased data - respect calls
    }

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{${$scriptParameterHashRef}{'sampleIDs'}});$sampleIDCounter++) {  #Collect infiles for all sampleIDs
	
	my $infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }{'pGATKBaseRecalibration'}{'fileEnding'};

	## Check if any files for this sampleID were merged previously to set infile and PicardToolsMergeSwitch to enable correct handling of number of infiles to process
	my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles(\%{$scriptParameterHashRef}, \%sampleInfo, \%lane, \%infilesLaneNoEnding, ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter]);
	
	if ($PicardToolsMergeSwitch == 1) {  #Alignment BAM-files merged previously
	    
	    print $FILEHANDLE "-I ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$infileEnding.".bam ";  #InFile
	}
	else {  #No previous merge of alignment BAM-files
	    
	    for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infilesLaneNoEndingHashRef}{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] } });$infileCounter++) {  #For all infiles per lane
		
		my $infile = ${$infilesLaneNoEndingHashRef}{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }[$infileCounter];
		
		print $FILEHANDLE "-I ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$infileEnding.".bam ";  #InFile(s)
	    }
	}
    } 
    print $FILEHANDLE "-L: ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$infileEnding.$callType.".vcf ";  #Limit to  (family vcf)
    print $FILEHANDLE "-V: ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$infileEnding.$callType.".vcf ";  #InFile (family vcf)
    print $FILEHANDLE "-o ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType.".vcf", "\n\n";  #OutFile

    ## Copies file from temporary directory.
    print $FILEHANDLE "## Copy file from temporary directory\n";
    &MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType.".vcf*", $outFamilyDirectory."/", $FILEHANDLE);
    print $FILEHANDLE "wait", "\n\n";

    &RemoveDirectory(\${$scriptParameterHashRef}{'tempDirectory'}, $FILEHANDLE);
   
    close($FILEHANDLE);

    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'dependencies' => 2, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}


sub GATKPhaseByTransmission {

##GATKPhaseByTransmission
    
##Function : GATK PhaseByTransmission computes the most likely genotype combination and phases trios and parent/child pairs given their genotype likelihoods and a mutation prior and phases all sites were parent/child transmission can be inferred unambiguously.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $familyID, $aligner, $callType, $programName
##         : $parameterHashRef       => The parameter hash {REF}
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef      => Info on samples and family hash {REF}
##         : $familyID               => The familyID
##         : $aligner                => The aligner used in the analysis
##         : $callType               => The variant call type
##         : $programName            => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $familyID = $_[3];
    my $aligner = $_[4];
    my $callType = $_[5];
    my $programName = $_[6];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header    
    my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					   {'directoryID' => $familyID,
					    'programName' => $programName,
					    'programDirectory' => lc($aligner."/gatk"),
					    'callType' => $callType,
					    'processTime' => 3,
					    'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
					   });
    
    ## Assign directories
    my $familyFileDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID;
    my $inFamilyDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID."/".$aligner."/gatk";
    my $outFamilyDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID."/".$aligner."/gatk";

    my $infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'pGATKVariantRecalibration'}{'fileEnding'};
    my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{"p".$programName}{'fileEnding'};
    
    unless (-f $familyFileDirectory."/".$familyID.".fam") {  #Check to see if file already exists

	print $FILEHANDLE "#Generating '.fam' file for GATK PhaseByTransmission","\n\n";
	print $FILEHANDLE q?perl -nae 'print $F[0], "\t", $F[1], "\t", $F[2], "\t", $F[3], "\t", $F[4], "\t", $F[5], "\n";' ?.${$scriptParameterHashRef}{'pedigreeFile'}." > ".$familyFileDirectory."/".$familyID.".fam", "\n\n";
    }

    ## Copy file(s) to temporary directory
    print $FILEHANDLE "## Copy file(s) to temporary directory\n"; 
    &MigrateFileToTemp($FILEHANDLE, 
			{'path' => $inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf*",
			 'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
			});
    print $FILEHANDLE "wait", "\n\n";

    ## GATK PhaseByTransmission    
    print $FILEHANDLE "## GATK PhaseByTransmission","\n";

    ## Writes java core commands to filehandle.
    &JavaCore($FILEHANDLE,
	      {'memoryAllocation' => "Xmx4g",
	       'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
	       'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
	      });
    
    print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
    print $FILEHANDLE "-l INFO ";  #Set the minimum level of logging
    print $FILEHANDLE "-T PhaseByTransmission ";  #Type of analysis to run
    print $FILEHANDLE "-R ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference file
    print $FILEHANDLE "-V: ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$infileEnding.$callType.".vcf ";  #InFile (family vcf)

    ## Check if "--pedigree" and "--pedigreeValidationType" should be included in analysis
    &GATKPedigreeFlag(\%{$scriptParameterHashRef}, $FILEHANDLE, $familyFileDirectory, "SILENT", "GATKPhaseByTransmission");

    print $FILEHANDLE "-o ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType.".vcf", "\n\n";  #OutFile

    ## Copies file from temporary directory.
    print $FILEHANDLE "## Copy file from temporary directory\n";
    &MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType.".vcf*", $outFamilyDirectory."/", $FILEHANDLE);
    print $FILEHANDLE "wait", "\n\n";

    &RemoveDirectory(\${$scriptParameterHashRef}{'tempDirectory'}, $FILEHANDLE);

    close($FILEHANDLE);

    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'dependencies' => 1, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}


sub SnpEff {
 
##SnpEff
    
##Function : SnpEff annotates variants from different sources.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $fileInfoHashRef, $familyID, $aligner, $callType, $programName
##         : $parameterHashRef       => The parameter hash {REF}
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef      => Info on samples and family hash {REF}
##         : $fileInfoHashRef        => The fileInfo hash {REF}
##         : $familyID               => The familyID
##         : $aligner                => The aligner used in the analysis
##         : $callType               => The variant call type
##         : $programName            => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $fileInfoHashRef = $_[3];
    my $familyID = $_[4]; 
    my $aligner = $_[5];
    my $callType = $_[6];
    my $programName = $_[7];
    
    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle

    ## Set the number of cores to allocate per sbatch job.
    my $nrCores = &NrofCoresPerSbatch(\%{$scriptParameterHashRef}, scalar(keys %{${$fileInfoHashRef}{'pSnpEff'}{'snpSiftAnnotationFiles'}}) + scalar(@contigs));  #Detect the number of cores to use from (snpSiftAnnotationFiles and dbNSFP (=+1)
    
    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header    
    my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					   {'directoryID' => $familyID,
					    'programName' => $programName,
					    'programDirectory' => lc($aligner."/gatk"),
					    'callType' => $callType,
					    'nrofCores' => $nrCores,
					    'processTime' => 10,
					    'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
					   });

    ## Assign directories
    my $inFamilyDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID."/".$aligner."/gatk";
    my $outFamilyDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID."/".$aligner."/gatk";

    my $infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'pVCFParser'}{'fileEnding'};
    my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{"p".$programName}{'fileEnding'};
    my $analysisType = "";

    for (my $VEPOutputFilesCounter=0;$VEPOutputFilesCounter<$VEPOutputFiles;$VEPOutputFilesCounter++) {

	my $coreCounter = 1;

	if ($VEPOutputFilesCounter == 1) {
    
	    $analysisType = ".selected";  #SelectFile variants
	}

	## Copy file(s) to temporary directory
	print $FILEHANDLE "## Copy file(s) to temporary directory\n"; 
	&MigrateFileToTemp($FILEHANDLE, 
			    {'path' => $inFamilyDirectory."/".$familyID.$infileEnding.$callType.$analysisType.".vcf*",
			     'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
			    });
	print $FILEHANDLE "wait", "\n\n";

	## SnpSift Annotation
	print $FILEHANDLE "## SnpSift Annotation","\n";
	my $annotationFileCounter = 0;

	for my $annotationFile (keys %{${$fileInfoHashRef}{'pSnpEff'}{'snpSiftAnnotationFiles'}}) {

	    my $annotationFileNoEnding = &RemoveFileEnding(\$annotationFile, ".gz");

	    &PrintWait(\$annotationFileCounter, \$nrCores, \$coreCounter, $FILEHANDLE);
	    
	    ## Writes java core commands to filehandle.
	    &JavaCore($FILEHANDLE,
		      {'memoryAllocation' => "Xmx4g",
		       'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
		       'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		      });
	    
	    print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'snpEffPath'}."/SnpSift.jar ";
	    print $FILEHANDLE "annotate ";
	    if (defined(${$fileInfoHashRef}{'pSnpEff'}{'snpSiftAnnotationFiles'}{$annotationFile})) {

		print $FILEHANDLE "-name SnpSift_ ";  #Prepend 'str' to all annotated INFO fields 
		print $FILEHANDLE "-info ".join(',', @{${$fileInfoHashRef}{'pSnpEff'}{'snpSiftAnnotationFiles'}{$annotationFile}})." ";  #Database
	    }
	    print $FILEHANDLE ${$scriptParameterHashRef}{'referencesDir'}."/".$annotationFile." ";  #Database
	    print $FILEHANDLE ${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$infileEnding.$callType.$analysisType.".vcf ";  #Infile
	    print $FILEHANDLE "| ";  #Pipe
	    print $FILEHANDLE "perl ".${$scriptParameterHashRef}{'inScriptDir'}."/vcfParser.pl ";  #Parses the vcf output
	    print $FILEHANDLE "> ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType.$analysisType."_".$annotationFileNoEnding." &","\n\n";  #outfile
	    $annotationFileCounter++;  #Increment counter
	}

	## SnpSift dbNSFP
	if (scalar(@{${$scriptParameterHashRef}{'snpSiftDbNSFPAnnotations'}}) > 0) {

	    my $combinedNrCores;

	    print $FILEHANDLE "## SnpSift dbNSFP","\n";

	    for (my $contigsCounter=0;$contigsCounter<scalar(@contigs);$contigsCounter++) {
		
		$combinedNrCores = $annotationFileCounter + $contigsCounter;

		&PrintWait(\$combinedNrCores, \$nrCores, \$coreCounter, $FILEHANDLE);

		## Writes java core commands to filehandle.
		&JavaCore($FILEHANDLE,
			  {'memoryAllocation' => "Xmx500m",
			   'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
			   'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
			  });
		
		print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'snpEffPath'}."/SnpSift.jar ";
		print $FILEHANDLE "dbnsfp ";
		print $FILEHANDLE "-db ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'snpSiftDbNSFPFile'}." ";  #DbNSFP file
		print $FILEHANDLE "-f ";  #fields to add
		print $FILEHANDLE join(',', @{${$scriptParameterHashRef}{'snpSiftDbNSFPAnnotations'}})." ";  #Databases
		print $FILEHANDLE "<( ";  #Pipe into SnpSift
		print $FILEHANDLE "cat ";  #Read infile
		print $FILEHANDLE ${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$infileEnding.$callType.$analysisType.".vcf ";  #Infile
		print $FILEHANDLE "| ";  #Pipe

		## Writes java core commands to filehandle.
		&JavaCore($FILEHANDLE,
			  {'memoryAllocation' => "Xmx500m",
			   'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
			   'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
			  });
		
		print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'snpEffPath'}."/SnpSift.jar ";
		print $FILEHANDLE "filter ";  #Parallalize per contig for speed
		print $FILEHANDLE q?"( CHROM = '?.$contigs[$contigsCounter].q?' )"?;
		print $FILEHANDLE ") ";  #End Pipe into SnpSift
		print $FILEHANDLE "> ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType.$analysisType."_".$contigs[$contigsCounter]."_dbnsfp.vcf &","\n\n";  #outfile	       
	    }
	    print $FILEHANDLE "wait", "\n\n";
	    
	    ## Merge dbNSFP splitted vcfs 
	    print $FILEHANDLE "## Merge dbNSFP splitted vcfs","\n";

	    &ConcatenateVCFs(\%{$scriptParameterHashRef}, $FILEHANDLE, \@contigs, ${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType.$analysisType."_", "_dbnsfp.vcf", ${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType.$analysisType."_dbnsfp.vcf", "reOrderHeader");	  
	}

	## GATK CombineVariants
	print $FILEHANDLE "## GATK CombineVariants","\n";  #We will loose the AD:PL info but there is currently no way around this usng bcftools, snpSift split -j
	
	## Writes java core commands to filehandle.
	&JavaCore($FILEHANDLE,
		  {'memoryAllocation' => "Xmx4g",
		   'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
		   'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		  });

	print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	print $FILEHANDLE "-l INFO ";  #Set the minimum level of logging
	print $FILEHANDLE "-T CombineVariants ";  #Type of analysis to run
	print $FILEHANDLE "-R ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference file
	print $FILEHANDLE "-genotypeMergeOptions UNSORTED ";  #Take the genotypes in any order. Should be fine since the same order and number of samples exists in all files
	for my $annotationFile (keys %{${$fileInfoHashRef}{'pSnpEff'}{'snpSiftAnnotationFiles'}}) {

	    my $annotationFileNoEnding = &RemoveFileEnding(\$annotationFile, ".gz");
	    print $FILEHANDLE "-V: ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType.$analysisType."_".$annotationFileNoEnding." ";
	}
	if (scalar(@{${$scriptParameterHashRef}{'snpSiftDbNSFPAnnotations'}}) > 0) { #DbNSFP processed
	    
	    print $FILEHANDLE "-V: ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType.$analysisType."_dbnsfp.vcf ";
	}
	print $FILEHANDLE "-o ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType.$analysisType.".vcf ", "\n\n";  #OutFile

	## Copies file from temporary directory.
	print $FILEHANDLE "## Copy file from temporary directory\n";
	&MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType.$analysisType.".vcf*", $outFamilyDirectory."/", $FILEHANDLE);
	print $FILEHANDLE "wait", "\n\n";
    }

    &RemoveDirectory(\${$scriptParameterHashRef}{'tempDirectory'}, $FILEHANDLE);

    close($FILEHANDLE);
    
    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {
	
	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'dependencies' => 1, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}


sub VCFParser {
 
##VCFParser
    
##Function : VCFParser performs parsing of VariantEffectPredictor annotated variants
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $sampleInfoHashRef$familyID, $aligner, $callType, $programName
##         : $parameterHashRef       => The parameter hash {REF}
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef      => Info on samples and family hash {REF}
##         : $familyID               => The familyID
##         : $aligner                => The aligner used in the analysis
##         : $callType               => The variant call type
##         : $programName            => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $familyID = $_[3]; 
    my $aligner = $_[4];
    my $callType = $_[5];
    my $programName = $_[6];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header    
    my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					   {'directoryID' => $familyID,
					    'programName' => $programName,
					    'programDirectory' => lc($aligner."/gatk"),
					    'callType' => $callType,
					    'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
					   });
    
    ## Assign directories
    my $inFamilyDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID."/".$aligner."/gatk";
    my $outFamilyDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID."/".$aligner."/gatk";

    my $infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'pVariantEffectPredictor'}{'fileEnding'};
    my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{"p".$programName}{'fileEnding'};

    ## VCFParser
    print $FILEHANDLE "## VCFParser","\n\n";
    print $FILEHANDLE "perl ".${$scriptParameterHashRef}{'inScriptDir'}."/vcfParser.pl ";  #Parses the VEP output to tab-sep format
    print $FILEHANDLE $outFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf ";  #Infile
    
    if (${$scriptParameterHashRef}{'pVariantEffectPredictor'} > 0) {

	print $FILEHANDLE "--parseVEP ".${$scriptParameterHashRef}{'vcfParserVepTranscripts'}." ";  #Parse VEP transcript specific entries
    }
    if (${$scriptParameterHashRef}{'vcfParserRangeFeatureFile'} ne "noUserInfo") {

	print $FILEHANDLE "-rf ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'vcfParserRangeFeatureFile'}." ";  #List of genes to analyse separately	
	
	if (scalar(@{${$scriptParameterHashRef}{'vcfParserRangeFeatureAnnotationColumns'}}) > 0) {

	    print $FILEHANDLE "-rf_ac ";  #Range annotation columns
	    print $FILEHANDLE join(',', @{${$scriptParameterHashRef}{'vcfParserRangeFeatureAnnotationColumns'}})." ";	    
	}
    }
    if (${$scriptParameterHashRef}{'vcfParserSelectFile'} ne "noUserInfo") {

	print $FILEHANDLE "-sf ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'vcfParserSelectFile'}." ";  #List of genes to analyse separately
	print $FILEHANDLE "-sf_mc ".${$scriptParameterHashRef}{'vcfParserSelectFileMatchingColumn'}." ";  #Column of HGNC Symbol in SelectFile (-sf)

	if (scalar(@{${$scriptParameterHashRef}{'vcfParserSelectFeatureAnnotationColumns'}}) > 0) {

	    print $FILEHANDLE "-sf_ac ";  #Select annotation columns
	    print $FILEHANDLE join(',', @{${$scriptParameterHashRef}{'vcfParserSelectFeatureAnnotationColumns'}})." ";	    
	}
	print $FILEHANDLE "-sof ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.".selected.vcf ";
    }
    print $FILEHANDLE "> ".$outFamilyDirectory."/".$familyID.$outfileEnding.$callType.".vcf";  #outfile
    print $FILEHANDLE "\n\n";

    close($FILEHANDLE);

    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'dependencies' => 1, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}


sub VariantEffectPredictor {
 
##VariantEffectPredictor
    
##Function : VariantEffectPredictor performs annotation of variants.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $familyID, $aligner, $callType, $programName
##         : $parameterHashRef       => The parameter hash {REF}
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef      => Info on samples and family hash {REF}
##         : $familyID               => The familyID
##         : $aligner                => The aligner used in the analysis
##         : $callType               => The variant call type
##         : $programName            => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $familyID = $_[3]; 
    my $aligner = $_[4];
    my $callType = $_[5];
    my $programName = $_[6];

    my $FILEHANDLE = IO::Handle->new();  #Create anonymous filehandle
    my $nrForkes = 4;

    ## Set the number of cores to allocate per sbatch job.
    my $nrCores = &NrofCoresPerSbatch(\%{$scriptParameterHashRef}, scalar(@contigs));  #Detect the number of cores to use
    $nrCores = floor($nrCores / $nrForkes);  #Adjust for the number of forks 

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($fileName, $stdoutPath, $stderrPath) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
								     {'directoryID' => $familyID,
								      'programName' => $programName,
								      'programDirectory' => lc($aligner."/gatk"),
								      'callType' => $callType,
								      'nrofCores' => ${$scriptParameterHashRef}{'maximumCores'},
								      'processTime' => 10,
								      'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
								     });
    my ($volume, $directories, $stderrFile) = File::Spec->splitpath($stderrPath);  #Split to enable submission to &SampleInfoQC later
    
    ## Assign directories
    my $inFamilyDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID."/".$aligner."/gatk";
    my $outFamilyDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID."/".$aligner."/gatk";

    my $infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'pGATKVariantRecalibration'}{'fileEnding'};
    my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{"p".$programName}{'fileEnding'};
    my $coreCounter = 1;

    ## Copy file(s) to temporary directory
    print $FILEHANDLE "## Copy file(s) to temporary directory\n"; 
    &MigrateFileToTemp($FILEHANDLE, 
			{'path' => $inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf*",
			 'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
			});
    print $FILEHANDLE "wait", "\n\n";

    ## VariantEffectPredictor
    print $FILEHANDLE "## VariantEffectPredictor","\n";

    for (my $contigsCounter=0;$contigsCounter<scalar(@contigs);$contigsCounter++) {

	&PrintWait(\$contigsCounter, \$nrCores, \$coreCounter, $FILEHANDLE);
	
	print $FILEHANDLE "perl ".${$scriptParameterHashRef}{'vepDirectoryPath'}."/variant_effect_predictor.pl ";  #VEP script 
	print $FILEHANDLE "--dir_cache ".${$scriptParameterHashRef}{'vepDirectoryCache'}." ";  #Specify the cache directory to use
	print $FILEHANDLE "--cache ";  #Enables use of the cache.
	print $FILEHANDLE "--force_overwrite ";  #force the overwrite of the existing file
	print $FILEHANDLE "--vcf ";  #Writes output in VCF format.
	print $FILEHANDLE "--fork ".$nrForkes." ";  #Enable forking, using the specified number of forks.
	print $FILEHANDLE "--buffer_size 20000 ";  #Sets the internal buffer size, corresponding to the number of variations that are read in to memory simultaneously 
	print $FILEHANDLE "--offline ";  #Use installed assembly version
	print $FILEHANDLE "--fasta ".${$scriptParameterHashRef}{'vepDirectoryCache'}." ";  #Use local fasta reference file
	print $FILEHANDLE "--chr ".$contigs[$contigsCounter]." ";

	for (my $vepFeatureCounter=0;$vepFeatureCounter<scalar(@{${$scriptParameterHashRef}{'vepFeatures'}});$vepFeatureCounter++) {
	    
	    print $FILEHANDLE "--".${$scriptParameterHashRef}{'vepFeatures'}[$vepFeatureCounter]." ";  #Add VEP features to the output.
	    
	    if ( (${$scriptParameterHashRef}{'vepFeatures'}[$vepFeatureCounter] eq "sift") || (${$scriptParameterHashRef}{'vepFeatures'}[$vepFeatureCounter] eq "polyphen") )  {  #Protein predictions
		
		print $FILEHANDLE "p ";  #Add prediction term 
	    }
	}
	print $FILEHANDLE "-i ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$infileEnding.$callType.".vcf ";  #InFile (family vcf)
	print $FILEHANDLE "-o ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType."_".$contigs[$contigsCounter].".vcf &", "\n\n";  #OutFile
    }
    print $FILEHANDLE "wait", "\n\n";

    ## Combine vcf files to 1
    &CombineVariants(\%{$scriptParameterHashRef}, $FILEHANDLE, \@contigs, ${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType."_", ".vcf", ${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType.".vcf");

    ## Copies file from temporary directory.
    print $FILEHANDLE "## Copy file from temporary directory\n";
    &MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType.".vcf*", $outFamilyDirectory."/", $FILEHANDLE);
    &MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType."_*.vcf_s**", $outFamilyDirectory."/", $FILEHANDLE);
    print $FILEHANDLE "wait", "\n\n";

    &RemoveDirectory(\${$scriptParameterHashRef}{'tempDirectory'}, $FILEHANDLE);

    close($FILEHANDLE);
	
    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	## Collect QC metadata info for later use                     	
	&SampleInfoQC(\%sampleInfo,
		       {'familyID' => ${$scriptParameterHashRef}{'familyID'},
			'programName' => $programName,
			'outDirectory' => $directories,
			'outFileEnding' => $stderrFile,
			'outDataType' => "infoDirectory"
		       });

	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'dependencies' => 1, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}


sub GATKVariantReCalibration { 

##GATKVariantReCalibration
    
##Function : GATK VariantRecalibrator/ApplyRecalibration.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $familyID, $aligner, $callType, $programName
##         : $parameterHashRef       => The parameter hash {REF}
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef      => Info on samples and family hash {REF}
##         : $familyID               => The familyID
##         : $aligner                => The aligner used in the analysis
##         : $callType               => The variant call type
##         : $programName            => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $familyID = $_[3]; 
    my $aligner = $_[4];
    my $callType = $_[5];
    my $programName = $_[6];

    my $FILEHANDLE = IO::Handle->new();  #Create anonymous filehandle
    my $nrCores = ${$scriptParameterHashRef}{'maximumCores'};

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					   {'directoryID' => $familyID,
					    'programName' => $programName,
					    'programDirectory' => lc($aligner."/gatk"),
					    'callType' => $callType,
					    'nrofCores' => $nrCores,
					    'processTime' => 10,
					    'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}."/gatk/intermediary"
					   });

    ## Assign directories
    my $outFamilyFileDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID;  #For ".fam" file
    my $inFamilyDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID."/".$aligner."/gatk";
    my $intermediarySampleDirectory = ${$scriptParameterHashRef}{'tempDirectory'}."/gatk/intermediary";
    my $outFamilyDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID."/".$aligner."/gatk";

    my $infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'pGATKGenoTypeGVCFs'}{'fileEnding'};
    my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{"p".$programName}{'fileEnding'};
    ## GATK ".fam" file check
    unless (-e ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID."/".$familyID.".fam") {  #Check to see if file already exists

	print $FILEHANDLE "#Generating '.fam' file for GATK VariantRecalibrator/ApplyRecalibration","\n\n";
	print $FILEHANDLE q?perl -nae 'print $F[0], "\t", $F[1], "\t", $F[2], "\t", $F[3], "\t", $F[4], "\t", $F[5], "\n";' ?.${$scriptParameterHashRef}{'pedigreeFile'}." > ".$outFamilyFileDirectory."/".$familyID.".fam", "\n\n";
    }

    ## Detects if there are different capture kits across sampleIDs. Creates a temporary merged interval_list for all interval_list that have been supplied and returns temporary list. Will also extract specific contigs if requested and return that list if enabled.
    my $contigIntervalListFile = &GATKTargetListFlag(\%scriptParameter, $FILEHANDLE);
    
    ## Copy file(s) to temporary directory
    print $FILEHANDLE "## Copy file(s) to temporary directory\n"; 
    &MigrateFileToTemp($FILEHANDLE, 
			{'path' => $inFamilyDirectory."/".$familyID.$infileEnding.$callType.".vcf*",
			 'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
			});
    print $FILEHANDLE "wait", "\n\n";
    
    ## GATK VariantRecalibrator
    my @modes = ("SNP","INDEL");

    if ( (${$scriptParameterHashRef}{'analysisType'} eq "exomes") || (${$scriptParameterHashRef}{'analysisType'} eq "rapid") ) {  #Exome/rapid analysis

	@modes = ("BOTH");
    }

    for (my $modeCounter=0;$modeCounter<scalar(@modes);$modeCounter++) {  #SNP and INDEL will be recalibrated successively in the same file because when you specify eg SNP mode, the indels are emitted without modification, and vice-versa. Exome and Rapid will be processed using mode BOTH since there are to few INDELS to use in the recalibration model even though using 30 exome BAMS in Haplotypecaller step. 

	print $FILEHANDLE "## GATK VariantRecalibrator","\n";

	## Writes java core commands to filehandle.
	&JavaCore($FILEHANDLE,
		  {'memoryAllocation' => "Xmx6g",
		   'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
		   'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		  });

	print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	print $FILEHANDLE "-l INFO ";  #Set the minimum level of logging
	print $FILEHANDLE "-T VariantRecalibrator ";  #Type of analysis to run
	print $FILEHANDLE "-R ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference file
	print $FILEHANDLE "-recalFile ".$intermediarySampleDirectory."/".$familyID.$infileEnding.$callType.".intervals ";
	print $FILEHANDLE "-rscriptFile ".$intermediarySampleDirectory."/".$familyID.$infileEnding.$callType.".intervals.plots.R ";
	print $FILEHANDLE "-tranchesFile ".$intermediarySampleDirectory."/".$familyID.$infileEnding.$callType.".intervals.tranches ";

	if ( (${$scriptParameterHashRef}{'analysisType'} eq "exomes") || (${$scriptParameterHashRef}{'analysisType'} eq "rapid") ) {  #Exome/rapid analysis use combined reference for more power
	
	    print $FILEHANDLE "-L ".$contigIntervalListFile." ";  #Target list file (merged or original)			
	    print $FILEHANDLE "-input ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$infileEnding.$callType.".vcf ";  #Infile HaplotypeCaller combined vcf which used reference gVCFs to create combined vcf (30> samples gCVFs)
	}
	else {  #WGS
	    
	    if ($modes[$modeCounter] eq "SNP") {
	
		print $FILEHANDLE "-input ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$infileEnding.$callType.".vcf ";
	    }
	    if ($modes[$modeCounter] eq "INDEL") {#Use created recalibrated snp vcf as input
	
		print $FILEHANDLE "-input ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType.".SNV.vcf ";
	    }
	    print $FILEHANDLE "-an DP ";  #The names of the annotations which should used for calculations. NOTE: Not to be used with hybrid capture
	}
	if ( ($modes[$modeCounter] eq "SNP") || ($modes[$modeCounter] eq "BOTH") ) {
	    
	    print $FILEHANDLE "-resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'GATKVariantReCalibrationTrainingSetHapMap'}." ";  #A list of sites for which to apply a prior probability of being correct but which aren't used by the algorithm
	    print $FILEHANDLE "-resource:omni,VCF,known=false,training=true,truth=false,prior=12.0 ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'GATKVariantReCalibrationTrainingSet1000GOmni'}." ";  #A list of sites for which to apply a prior probability of being correct but which aren't used by the algorithm
	    print $FILEHANDLE "-resource:1000G,known=false,training=true,truth=false,prior=10.0 ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'GATKVariantReCalibrationTrainingSet1000GSNP'}." ";  #A list of sites for which to apply a prior probability of being correct but which aren't used by the algorithm
	    print $FILEHANDLE "-an QD ";  #The names of the annotations which should used for calculations
	}
	if ( ($modes[$modeCounter] eq "INDEL") || ($modes[$modeCounter] eq "BOTH") ) {
	    
	    print $FILEHANDLE "-resource:mills,VCF,known=true,training=true,truth=true,prior=12.0 ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'GATKVariantReCalibrationTrainingSetMills'}." ";  #A list of sites for which to apply a prior probability of being correct but which aren't used by the algorithm
	}
	print $FILEHANDLE "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'GATKVariantReCalibrationTrainingSetDbSNP'}." ";  #A list of sites for which to apply a prior probability of being correct but which aren't used by the algorithm
    
	print $FILEHANDLE "-an MQRankSum ";  #The names of the annotations which should used for calculations
	print $FILEHANDLE "-an ReadPosRankSum ";  #The names of the annotations which should used for calculations
	print $FILEHANDLE "-an FS ";  #The names of the annotations which should used for calculations
	print $FILEHANDLE "--mode ".$modes[$modeCounter]." ";  #Recalibration mode to employ (SNP|INDEL|BOTH)
	print $FILEHANDLE "-nt ".${$scriptParameterHashRef}{'maximumCores'}." ";  #How many data threads should be allocated to running this analysis    

	## Check if "--pedigree" and "--pedigreeValidationType" should be included in analysis
	&GATKPedigreeFlag(\%{$scriptParameterHashRef}, $FILEHANDLE, $outFamilyFileDirectory, "SILENT", "GATKVariantRecalibration");  #Passing filehandle directly to sub routine using "*". Sub routine prints "--pedigree file" for family
	
	## GATK ApplyRecalibration
	print $FILEHANDLE "\n\n## GATK ApplyRecalibration","\n";

	## Writes java core commands to filehandle.
	&JavaCore($FILEHANDLE,
		  {'memoryAllocation' => "Xmx6g",
		   'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
		   'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		  });

	print $FILEHANDLE  "-jar ".${$scriptParameterHashRef}{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	print $FILEHANDLE "-l INFO ";  #Set the minimum level of logging
	print $FILEHANDLE "-T ApplyRecalibration ";
	print $FILEHANDLE "-R ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference file
	print $FILEHANDLE "-recalFile ".$intermediarySampleDirectory."/".$familyID.$infileEnding.$callType.".intervals ";
	print $FILEHANDLE "-tranchesFile ".$intermediarySampleDirectory."/".$familyID.$infileEnding.$callType.".intervals.tranches ";

	if ( (${$scriptParameterHashRef}{'analysisType'} eq "exomes") || (${$scriptParameterHashRef}{'analysisType'} eq "rapid")) {  #Exome/rapid analysis use combined reference for more power
	    
	    print $FILEHANDLE "-L ".$contigIntervalListFile." ";#Target list file (merged or original)		
	    print $FILEHANDLE "-input ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$infileEnding.$callType.".vcf ";  #Infile HaplotypeCaller combined vcf which used reference gVCFs to create combined vcf file
	    print $FILEHANDLE "-o ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType."_filtered.vcf ";
	}
	else  {  #WGS
	    
	    if ($modes[$modeCounter] eq "SNP") {
		
		print $FILEHANDLE "-input ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$infileEnding.$callType.".vcf ";
		print $FILEHANDLE "-o ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType.".SNV.vcf ";
	    }
	    if ($modes[$modeCounter] eq "INDEL") {#Use created recalibrated snp vcf as input
	
		print $FILEHANDLE "-input ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType.".SNV.vcf ";
		print $FILEHANDLE "-o ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType.".vcf ";
	    }
	}
	print $FILEHANDLE "--ts_filter_level ".${$scriptParameterHashRef}{'GATKVariantReCalibrationTSFilterLevel'}." ";

	## Check if "--pedigree" and "--pedigreeValidationType" should be included in analysis
	&GATKPedigreeFlag(\%{$scriptParameterHashRef}, $FILEHANDLE, $outFamilyFileDirectory, "SILENT", "GATKVariantRecalibration");  #Passing filehandle directly to sub routine using "*". Sub routine prints "--pedigree file" for family    
	print $FILEHANDLE "--mode ".$modes[$modeCounter]." ";  #Recalibration mode to employ (SNP|INDEL|BOTH)
    }

    ## GATK SelectVariants

    ## Removes all genotype information for exome ref and recalulates meta-data info for remaining samples in new file.
    if ( (${$scriptParameterHashRef}{'analysisType'} eq "exomes") || (${$scriptParameterHashRef}{'analysisType'} eq "rapid") ) {  #Exome/rapid analysis
	
	print $FILEHANDLE "\n\n## GATK SelectVariants","\n";

	## Writes java core commands to filehandle.
	&JavaCore($FILEHANDLE,
		   {'memoryAllocation' => "Xmx2g",
		    'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
		    'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		   });

	print $FILEHANDLE  "-jar ".${$scriptParameterHashRef}{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	print $FILEHANDLE "-l INFO ";  #Set the minimum level of logging
	print $FILEHANDLE "-T SelectVariants ";  #Type of analysis to run
	print $FILEHANDLE "-R ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference file	
	print $FILEHANDLE "-L ".$contigIntervalListFile." ";#Target list file (merged or original)
	print $FILEHANDLE "-env ";  #Don't include loci found to be non-variant after the subsetting procedure. 
	print $FILEHANDLE "-V: ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType."_filtered.vcf ";  #InFile

	for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{${$scriptParameterHashRef}{'sampleIDs'}});$sampleIDCounter++) {  #For all sampleIDs
		
	    print $FILEHANDLE "-sn ".${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter]." ";  #Include genotypes from this sample
	}
	print $FILEHANDLE "-o ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType.".vcf ";  #OutFile
	print $FILEHANDLE " &"; 
	
	## Produces another vcf file containing non-variant loci (useful for example in MAF comparisons), but is not used downstream in MIP
	if (${$scriptParameterHashRef}{'GATKVariantReCalibrationexcludeNonVariantsFile'} ne "false") {

	    print $FILEHANDLE "\n\n#GATK SelectVariants","\n\n";

	    ## Writes java core commands to filehandle.
	    &JavaCore($FILEHANDLE,
		      {'memoryAllocation' => "Xmx2g",
		       'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
		       'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		      });
	    
	    print $FILEHANDLE  "-jar ".${$scriptParameterHashRef}{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	    print $FILEHANDLE "-l INFO ";  #Set the minimum level of logging
	    print $FILEHANDLE "-T SelectVariants ";  #Type of analysis to run
	    print $FILEHANDLE "-R ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference file	
	    print $FILEHANDLE "-L ".$contigIntervalListFile." ";#Target list file (merged or original)
	    print $FILEHANDLE "-V: ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType."_filtered.vcf ";  #InFile
	    print $FILEHANDLE "-o ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType."_incnonvariantloci.vcf ";  #OutFile
	    
	    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{${$scriptParameterHashRef}{'sampleIDs'}});$sampleIDCounter++) {  #For all sampleIDs
		
		print $FILEHANDLE "-sn ".${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter]." ";  #Include genotypes from this sample
	    }
	    print $FILEHANDLE "\n\nwait\n\n";

	    ## Copies file from temporary directory.
	    print $FILEHANDLE "## Copy file from temporary directory\n"; 
	    &MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType."_incnonvariantloci.vcf*", $outFamilyDirectory."/", $FILEHANDLE);
	}
	print $FILEHANDLE "\n\nwait\n\n";
    }

    ## Copies file from temporary directory.
    print $FILEHANDLE "## Copy file from temporary directory\n";
    &MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType.".vcf*", $outFamilyDirectory."/", $FILEHANDLE);
    &MigrateFileFromTemp($intermediarySampleDirectory."/".$familyID.$infileEnding.$callType.".intervals.tranches.pdf", $outFamilyDirectory."/", $FILEHANDLE);
    print $FILEHANDLE "wait", "\n\n";

    &RemoveDirectory(\${$scriptParameterHashRef}{'tempDirectory'}, $FILEHANDLE);

    close($FILEHANDLE);   
    	
    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	## Collect QC metadata info for later use
	&SampleInfoQC(\%sampleInfo,
		       {'familyID' => ${$scriptParameterHashRef}{'familyID'},
			'programName' => "pedigreeCheck",
			'outDirectory' => $outFamilyDirectory,
			'outFileEnding' => $familyID.$outfileEnding.$callType.".vcf",
			'outDataType' => "infileDependent"
		       });
	${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'MostCompleteVCF'}{'Path'} = $outFamilyDirectory."/".$familyID.$outfileEnding.$callType.".vcf";	
	
	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'dependencies' => 1, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}


sub GATKGenoTypeGVCFs { 

##GATKGenoTypeGVCFs
    
##Function : GATK GenoTypeGVCFs. 
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $laneHashRef, $familyID, $aligner, $callType, $programName
##         : $parameterHashRef       => The parameter hash {REF}
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef      => Info on samples and family hash {REF}
##         : $laneHashRef            => The lane info hash {REF}
##         : $familyID               => The familyID
##         : $aligner                => The aligner used in the analysis
##         : $callType               => The variant call type
##         : $programName            => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $laneHashRef = $_[3];
    my $familyID = $_[4]; 
    my $aligner = $_[5];
    my $callType = $_[6];
    my $programName = $_[7];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $nrCores = ${$scriptParameterHashRef}{'maximumCores'};
    
    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					   {'directoryID' => $familyID,
					    'programName' => $programName,
					    'programDirectory' => lc($aligner."/gatk"),
					    'callType' => $callType,
					    'nrofCores' => $nrCores,
					    'processTime' => 10,
					    'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
					   });

    ## Assign directories
    my $outFamilyDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID."/".$aligner."/gatk";
    my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{"p".$programName}{'fileEnding'}; 

    ## Collect infiles for all sampleIDs to enable migration to temporary directory
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{${$scriptParameterHashRef}{'sampleIDs'}});$sampleIDCounter++) {
	
	## Check if any files for this sampleID were merged previously to set infile and PicardToolsMergeSwitch to enable correct handling of number of infiles to process
	my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles(\%{$scriptParameterHashRef}, \%sampleInfo, \%lane, \%infilesLaneNoEnding, ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter]);
	
	my $inSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter]."/".$aligner."/gatk";
	my $infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }{'pGATKHaploTypeCaller'}{'fileEnding'};
	
	if ($PicardToolsMergeSwitch == 1) {  #Files were merged previously
	    
	    ## Copy file(s) to temporary directory
	    print $FILEHANDLE "## Copy file(s) to temporary directory\n"; 
	    &MigrateFileToTemp($FILEHANDLE, 
				{'path' => $inSampleDirectory."/".$infile.$infileEnding.".vcf*",
				 'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
				});
	    print $FILEHANDLE "wait", "\n\n";
	}
	else {
	    
	    my $lanes = join("",@{${$laneHashRef}{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }});  #Extract lanes
	    
	    ## Copy file(s) to temporary directory
	    print $FILEHANDLE "## Copy file(s) to temporary directory\n"; 
	    &MigrateFileToTemp($FILEHANDLE, 
				{'path' => $inSampleDirectory."/".${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter]."_lanes_".$lanes.$infileEnding.".vcf*",
				 'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
				});
	    print $FILEHANDLE "wait", "\n\n";
	}
    }

    ## GATK GenoTypeGVCFs
    print $FILEHANDLE "## GATK GenoTypeGVCFs","\n";
    
    ## Writes java core commands to filehandle.
    &JavaCore($FILEHANDLE,
	      {'memoryAllocation' => "Xmx4g",
	       'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
	       'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
	      });
    
    print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
    print $FILEHANDLE "-l INFO ";  #Set the minimum level of logging
    print $FILEHANDLE "-T GenotypeGVCFs ";  #Type of analysis to run
    print $FILEHANDLE "-R ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference file
    print $FILEHANDLE "-D ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'GATKHaploTypeCallerSNPKnownSet'}." ";  #Known SNPs to use for annotation SNPs
    print $FILEHANDLE "-nt 16 ";  #How many data threads should be allocated to running this analysis.

    if ( (${$scriptParameterHashRef}{'analysisType'} eq "exomes") || (${$scriptParameterHashRef}{'analysisType'} eq "rapid") ) {

	print $FILEHANDLE "-V ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'GATKGenoTypeGVCFsRefGVCF'}." ";
    }

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{${$scriptParameterHashRef}{'sampleIDs'}});$sampleIDCounter++) {  #Collect infiles for all sampleIDs

	## Check if any files for this sampleID were merged previously to set infile and PicardToolsMergeSwitch to enable correct handling of number of infiles to process
	my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles(\%{$scriptParameterHashRef}, \%sampleInfo, \%lane, \%infilesLaneNoEnding, ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter]);

	my $infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }{'pGATKHaploTypeCaller'}{'fileEnding'};
	
	if ($PicardToolsMergeSwitch == 1) {  #Alignment BAM-files merged previously
	    
	    print $FILEHANDLE "-V ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$infileEnding.".vcf ";  #InFile
	}
	else {  #No previous merge of alignment BAM-files

	    my $lanes = join("",@{${$laneHashRef}{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }});  #Extract lanes
	    print $FILEHANDLE "-V ".${$scriptParameterHashRef}{'tempDirectory'}."/".${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter]."_lanes_".$lanes.$infileEnding.".vcf ";  #InFile(s)
	} 
    } 
    print $FILEHANDLE "-o ".${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType.".vcf", "\n\n";  #OutFile

    ## Copies file from temporary directory.
    print $FILEHANDLE "## Copy file from temporary directory\n"; 
    &MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$familyID.$outfileEnding.$callType.".vcf*", $outFamilyDirectory."/", $FILEHANDLE);
    print $FILEHANDLE "wait", "\n\n";

    &RemoveDirectory(\${$scriptParameterHashRef}{'tempDirectory'}, $FILEHANDLE);
    
    close($FILEHANDLE);  

    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	## Collect QC metadata info for later use
	${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'MostCompleteVCF'}{'Path'} = $outFamilyDirectory."/".$familyID.$outfileEnding.$callType.".vcf";

	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'dependencies' => 1, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}


sub GATKHaploTypeCaller { 

##GATKHaploTypeCaller
    
##Function : GATKHaploTypeCaller. 
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $infilesLaneNoEndingHashRef, $laneHashRef, $sampleID, $aligner, $programName
##         : $parameterHashRef           => The parameter hash {REF}
##         : $scriptParameterHashRef     => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef          => Info on samples and family hash {REF}
##         : $infilesLaneNoEndingHashRef => The infile(s) without the ".ending" {REF}
##         : $laneHashRef                => The lane info hash {REF}
##         : $sampleID                   => The sampleID
##         : $aligner                    => The aligner used in the analysis
##         : $programName                => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $infilesLaneNoEndingHashRef = $_[3];
    my $laneHashRef = $_[4];
    my $sampleID = $_[5];
    my $aligner = $_[6];
    my $programName = $_[7];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $nrCores = ${$scriptParameterHashRef}{'maximumCores'};

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					   {'directoryID' => $sampleID,
					    'programName' => $programName,
					    'programDirectory' => lc($aligner."/gatk"),
					    'nrofCores' => $nrCores,
					    'processTime' => 30,
					    'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
					   });
    
    ## Assign directories
    my $inSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner."/gatk";
    my $outSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner."/gatk";

    my $infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pGATKBaseRecalibration'}{'fileEnding'};
    my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{"p".$programName}{'fileEnding'};

    ## Check if any files for this sampleID were merged previously to set infile and PicardToolsMergeSwitch to enable correct handling of number of infiles to process
    my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles(\%{$scriptParameterHashRef}, \%sampleInfo, \%lane, \%infilesLaneNoEnding, $sampleID);

    ## Detects if there are different capture kits across sampleIDs. Creates a temporary merged interval_list for all interval_list that have been supplied and returns temporary list. Will also extract specific contigs if requested and return that list if enabled.
    my $contigIntervalListFile = &GATKTargetListFlag(\%{$scriptParameterHashRef}, $FILEHANDLE);
    
    if ($PicardToolsMergeSwitch == 1) {  #Files were merged previously

	## Copy file(s) to temporary directory
	print $FILEHANDLE "## Copy file(s) to temporary directory\n"; 
	&MigrateFileToTemp($FILEHANDLE, 
			    {'path' => $inSampleDirectory."/".$infile.$infileEnding.".b*",
			     'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
			    });
	print $FILEHANDLE "wait", "\n\n";
    }
    else {

	## Copies files from source to temporary folder. Loop over files specified by $arrayRef and collects files from $extractArrayRef.
	&MigrateFilesToTemp(\%{$scriptParameterHashRef}, \@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, \@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, $FILEHANDLE,
			    {'inSampleDirectory' => $inSampleDirectory,
			     'nrCores' => $nrCores,
			     'fileEnding' => $infileEnding.".b*"
			    });
    }

    ## GATK HaplotypeCaller
    print $FILEHANDLE "## GATK HaplotypeCaller","\n";

    ## Writes java core commands to filehandle.
    &JavaCore($FILEHANDLE,
	      {'memoryAllocation' => "Xmx12g",
	       'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
	       'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
	      });
    
    print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
    print $FILEHANDLE "-l INFO ";  #Set the minimum level of logging
    print $FILEHANDLE "-T HaplotypeCaller ";  #Type of analysis to run
    print $FILEHANDLE "-R ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference file
    print $FILEHANDLE "-D ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'GATKHaploTypeCallerSNPKnownSet'}." ";  #Known SNPs to use for annotation SNPs
    print $FILEHANDLE "-stand_call_conf 30.0 ";  #The minimum phred-scaled confidence threshold at which variants should be called
    print $FILEHANDLE "-stand_emit_conf 30.0 ";  #The minimum phred-scaled confidence threshold at which variants should be emitted
    print $FILEHANDLE "-nct 8 ";  #Number of CPU Threads per data thread
    ##annotations to apply to variant calls
    print $FILEHANDLE "--annotation BaseQualityRankSumTest ";  
    print $FILEHANDLE "--annotation ChromosomeCounts ";
    print $FILEHANDLE "--annotation Coverage ";
    print $FILEHANDLE "--annotation FisherStrand ";
    print $FILEHANDLE "--annotation InbreedingCoeff ";
    print $FILEHANDLE "--annotation MappingQualityRankSumTest ";
    print $FILEHANDLE "--annotation MappingQualityZero ";
    print $FILEHANDLE "--annotation QualByDepth ";
    print $FILEHANDLE "--annotation RMSMappingQuality ";
    print $FILEHANDLE "--annotation ReadPosRankSumTest ";
    print $FILEHANDLE "--annotation SpanningDeletions ";
    print $FILEHANDLE "--annotation TandemRepeatAnnotator " ;
    print $FILEHANDLE "--annotation DepthPerAlleleBySample ";
    print $FILEHANDLE "-dcov ".${$scriptParameterHashRef}{'GATKDownSampleToCoverage'}." ";  #Coverage to downsample to at any given locus
    print $FILEHANDLE "--emitRefConfidence GVCF ";  #Mode for emitting experimental reference confidence scores. GVCF generates block summarized version of the BP_RESOLUTION data 
    print $FILEHANDLE "--variant_index_type LINEAR "; 
    print $FILEHANDLE "--variant_index_parameter 128000 ";
    print $FILEHANDLE "-pairHMM VECTOR_LOGLESS_CACHING ";  #Hardware specific optmization

    if ( (${$scriptParameterHashRef}{'analysisType'} eq "exomes") || (${$scriptParameterHashRef}{'analysisType'} eq "rapid") ) {  #Exome/rapid analysis

	print $FILEHANDLE "-L ".$contigIntervalListFile." ";  #Target list file (merged or original)
    }	
    if ($PicardToolsMergeSwitch == 1) {  #Alignment BAM-files merged previously
	    
	print $FILEHANDLE "-I ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$infileEnding.".bam ";  #InFile
	print $FILEHANDLE "-o ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.".vcf", "\n\n";  #OutFile

	## Copies file from temporary directory.
	print $FILEHANDLE "## Copy file from temporary directory\n"; 
	&MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.".vcf*", $outSampleDirectory."/", $FILEHANDLE);
	print $FILEHANDLE "wait", "\n\n";
    }
    else {  #No previous merge of alignment BAM-files
	
	my $lanes = join("",@{${$laneHashRef}{$sampleID}});  #Extract lanes

	for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} });$infileCounter++) {  #For all infiles per lane
	    
	    my $infile = ${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter];
	    
	    print $FILEHANDLE "-I ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$infileEnding.".bam ";  #InFile(s)
	} 
	print $FILEHANDLE "-o ".${$scriptParameterHashRef}{'tempDirectory'}."/".$sampleID."_lanes_".$lanes.$outfileEnding.".vcf", "\n\n";  #OutFile

	## Copies file from temporary directory.
	print $FILEHANDLE "## Copy file from temporary directory\n"; 
	&MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$sampleID."_lanes_".$lanes.$outfileEnding.".vcf*", $outSampleDirectory."/", $FILEHANDLE);
	print $FILEHANDLE "wait", "\n\n";
    } 

    &RemoveDirectory(\${$scriptParameterHashRef}{'tempDirectory'}, $FILEHANDLE);

    close($FILEHANDLE);  

    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'sampleID' => $sampleID,
			'dependencies' => 1, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}

sub GATKBaseReCalibration { 

##GATKBaseReCalibration
    
##Function : GATK BaseRecalibrator/PrintReads to recalibrate bases before variant calling. Both BaseRecalibrator/PrintReads will be executed within the same sbatch script.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $infilesLaneNoEndingHashRef, $sampleID, $aligner, $programName
##         : $parameterHashRef           => The parameter hash {REF}
##         : $scriptParameterHashRef     => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef          => Info on samples and family hash {REF}
##         : $infilesLaneNoEndingHashRef => The infile(s) without the ".ending" {REF}
##         : $sampleID                   => The sampleID
##         : $aligner                    => The aligner used in the analysis
##         : $programName                => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $infilesLaneNoEndingHashRef = $_[3];
    my $sampleID = $_[4];
    my $aligner = $_[5];
    my $programName = $_[6];

    my $FILEHANDLE = IO::Handle->new();  #Create anonymous filehandle
    my $nrCores = ${$scriptParameterHashRef}{'maximumCores'};

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					   {'directoryID' => $sampleID,
					    'programName' => $programName,
					    'programDirectory' => lc($aligner."/gatk"),
					    'nrofCores' => $nrCores,
					    'processTime' => 50,
					     'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}."/gatk/intermediary"
					   });

    ## Assign directories
    my $inSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner."/gatk";
    my $intermediarySampleDirectory = ${$scriptParameterHashRef}{'tempDirectory'}."/gatk/intermediary";
    my $outSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner."/gatk";

    my $infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pGATKRealigner'}{'fileEnding'};
    my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{"p".$programName}{'fileEnding'};

    ## Check if any files for this sampleID were merged previously to set infile and PicardToolsMergeSwitch to enable correct handling of number of infiles to process
    my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles(\%{$scriptParameterHashRef}, \%sampleInfo, \%lane, \%infilesLaneNoEnding, $sampleID);
    
    if ($PicardToolsMergeSwitch == 1) {  #Files were merged previously

	## Copy file(s) to temporary directory
	print $FILEHANDLE "## Copy file(s) to temporary directory\n"; 
	&MigrateFileToTemp($FILEHANDLE, 
			    {'path' => $inSampleDirectory."/".$infile.$infileEnding.".b*",
			     'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
			    });
	print $FILEHANDLE "wait", "\n\n";

	## GATK BaseRecalibrator
	print $FILEHANDLE "## GATK BaseRecalibrator","\n";

	## Writes java core commands to filehandle.
	&JavaCore($FILEHANDLE,
		  {'memoryAllocation' => "Xmx24g",
		   'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
		   'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		  });

	print $FILEHANDLE "-Djava.io.tmpdir=".${$scriptParameterHashRef}{'tempDirectory'}." ";  #Temporary Directory per chr
	print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	print $FILEHANDLE "-l INFO ";  #Set the minimum level of logging
	print $FILEHANDLE "-T BaseRecalibrator ";  #Type of analysis to run

	## Covariates to be used in the recalibration
	print $FILEHANDLE "-cov ReadGroupCovariate ";
	print $FILEHANDLE "-cov ContextCovariate ";
	print $FILEHANDLE "-cov CycleCovariate ";
	print $FILEHANDLE "-cov QualityScoreCovariate ";
	print $FILEHANDLE "-cov ReadGroupCovariate ";

	print $FILEHANDLE "-R ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference file
	print $FILEHANDLE "-knownSites ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'GATKBaseReCalibrationSNPKnownSet'}." ";
	print $FILEHANDLE "-nct ".${$scriptParameterHashRef}{'maximumCores'}." ";  #How many CPU threads should be allocated per data thread to running this analysis
	print $FILEHANDLE "-dcov ".${$scriptParameterHashRef}{'GATKDownSampleToCoverage'}." ";  #Coverage to downsample to at any given locus
	print $FILEHANDLE "-I ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$infileEnding.".bam ";  #InFile
	print $FILEHANDLE "-o ".$intermediarySampleDirectory."/".$infile.$infileEnding.".grp ", "\n\n";  #Recalibration table file
	
	## GATK PrintReads
	print $FILEHANDLE "## GATK PrintReads","\n";

	## Writes java core commands to filehandle.
	&JavaCore($FILEHANDLE,
		  {'memoryAllocation' => "Xmx24g",
		   'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
		   'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		  });
	
	print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	print $FILEHANDLE "-l INFO ";  #Set the minimum level of logging"-jar $gatk_path/GenomeAnalysisTK.
	print $FILEHANDLE "-T PrintReads ";  #Type of analysis to run
	print $FILEHANDLE "-R ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference file
	print $FILEHANDLE "-nct ".${$scriptParameterHashRef}{'maximumCores'}." ";  #How many CPU threads should be allocated per data thread to running this analysis	  
	print $FILEHANDLE "-dcov ".${$scriptParameterHashRef}{'GATKDownSampleToCoverage'}." ";  #Coverage to downsample to at any given locus  
	print $FILEHANDLE "-I ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$infileEnding.".bam ";  #InFile
	print $FILEHANDLE "-o ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.".bam ";  #OutFile
	print $FILEHANDLE "-BQSR ".$intermediarySampleDirectory."/".$infile.$infileEnding.".grp ", "\n\n";  #Recalibration table file

	## Copies file from temporary directory.
	print $FILEHANDLE "## Copy file from temporary directory\n"; 
	&MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.".b*", $outSampleDirectory."/", $FILEHANDLE);
	print $FILEHANDLE "wait", "\n\n";

	if ( (${$scriptParameterHashRef}{'pGATKBaseRecalibration'} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) { 
	    
	    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'} = $outSampleDirectory."/".$infile.$outfileEnding.".bam";
	}
    }
    else {  #no previous merge

	## Copies files from source to temporary folder. Loop over files specified by $arrayRef and collects files from $extractArrayRef.
	&MigrateFilesToTemp(\%{$scriptParameterHashRef}, \@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, \@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, $FILEHANDLE,
			    {'inSampleDirectory' => $inSampleDirectory,
			     'nrCores' => $nrCores,
			     'fileEnding' => $infileEnding.".b*"
			    });

	for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} });$infileCounter++) {  #For all infiles per lane
	    
	    my $infile = ${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter];

	    ## GATK BaseRecalibrator
	    print $FILEHANDLE "## GATK BaseRecalibrator","\n";

	    ## Writes java core commands to filehandle.
	    &JavaCore($FILEHANDLE,
		      {'memoryAllocation' => "Xmx24g",
		       'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
		       'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		      });
	    
	    print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	    print $FILEHANDLE "-l INFO ";  #Set the minimum level of logging
	    print $FILEHANDLE "-T BaseRecalibrator ";  #Type of analysis to run

	    ## Covariates to be used in the recalibration
	    print $FILEHANDLE "-cov ReadGroupCovariate "; 
	    print $FILEHANDLE "-cov ContextCovariate "; 
	    print $FILEHANDLE "-cov CycleCovariate ";
	    print $FILEHANDLE "-cov QualityScoreCovariate ";
	    print $FILEHANDLE "-cov ReadGroupCovariate ";

	    print $FILEHANDLE "-R ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference file
	    print $FILEHANDLE "-knownSites ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'GATKBaseReCalibrationSNPKnownSet'}." ";
	    print $FILEHANDLE "-nct ".${$scriptParameterHashRef}{'maximumCores'}." ";  #How many CPU threads should be allocated per data thread to running this analysis
	    print $FILEHANDLE "-dcov ".${$scriptParameterHashRef}{'GATKDownSampleToCoverage'}." ";  #Coverage to downsample to at any given locus	
	    print $FILEHANDLE "-I ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$infileEnding.".bam ";  #InFile
	    print $FILEHANDLE "-o ".$intermediarySampleDirectory."/".$infile.$infileEnding.".grp ", "\n\n";  #Recalibration table file

	    ## GATK PrintReads
	    print $FILEHANDLE "## GATK PrintReads","\n";
	    
	    ## Writes java core commands to filehandle.
	    &JavaCore($FILEHANDLE,
		      {'memoryAllocation' => "Xmx24g",
		       'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
		       'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		      });
	    
	    print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	    print $FILEHANDLE "-l INFO ";  #Set the minimum level of logging"-jar $gatk_path/GenomeAnalysisTK.
	    print $FILEHANDLE "-T PrintReads ";  #Type of analysis to run
	    print $FILEHANDLE "-R ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference file
	    print $FILEHANDLE "-nct ".${$scriptParameterHashRef}{'maximumCores'}." ";  #How many CPU threads should be allocated per data thread to running this analysis
	    print $FILEHANDLE "-dcov ".${$scriptParameterHashRef}{'GATKDownSampleToCoverage'}." ";  #Coverage to downsample to at any given locus
	    print $FILEHANDLE "-I ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$infileEnding.".bam ";  #InFile
	    print $FILEHANDLE "-o ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.".bam ";  #OutFile
	    print $FILEHANDLE "-BQSR ".$intermediarySampleDirectory."/".$infile.$infileEnding.".grp ", "\n\n";  #Recalibration table file

	    if ( (${$scriptParameterHashRef}{'pGATKBaseRecalibration'} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) { 

		${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'} = $outSampleDirectory."/".$infile.$outfileEnding.".bam";
	    }
	}

	## Copies files from temporary folder to source. Loop over files specified by $arrayRef and collects files from $extractArrayRef.
	&MigrateFilesFromTemp(\@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, \@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, $outSampleDirectory, ${$scriptParameterHashRef}{'tempDirectory'}, $nrCores, $outfileEnding.".b*", $FILEHANDLE);
    }

    &RemoveDirectory(\${$scriptParameterHashRef}{'tempDirectory'}, $FILEHANDLE);
    
    close($FILEHANDLE);

    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) { 

	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'sampleID' => $sampleID,
			'dependencies' => 1, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}


sub GATKReAligner { 

##GATKReAligner
    
##Function : GATK ReAlignerTargetCreator/IndelRealigner to rearrange reads around INDELs. Both ReAlignerTargetCreator and IndelRealigner will be executed within the same sbatch script.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $infilesLaneNoEndingHashRef, $sampleID, $aligner, $programName 
##         : $parameterHashRef           => The parameter hash {REF}
##         : $scriptParameterHashRef     => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef          => Info on samples and family hash {REF}
##         : $infilesLaneNoEndingHashRef => The infile(s) without the ".ending" {REF}
##         : $sampleID                   => The sampleID
##         : $aligner                    => The aligner used in the analysis
##         : $programName                => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $infilesLaneNoEndingHashRef = $_[3];
    my $sampleID = $_[4];
    my $aligner = $_[5];
    my $programName = $_[6];

    my $FILEHANDLE = IO::Handle->new(); #Create anonymous filehandle
    my $nrCores = ${$scriptParameterHashRef}{'maximumCores'};

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					   {'directoryID' => $sampleID,
					    'programName' => $programName,
					    'programDirectory' => lc($aligner."/gatk"),
					    'nrofCores' => $nrCores,
					    'processTime' => 40,
					    'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}."/gatk/intermediary"
					   });

    ## Assign directories
    my $inSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner;
    my $intermediarySampleDirectory = ${$scriptParameterHashRef}{'tempDirectory'}."/gatk/intermediary";
    my $outSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner."/gatk";

    my $infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pPicardToolsMarkduplicates'}{'fileEnding'};
    my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{"p".$programName}{'fileEnding'};

    ## Check if any files for this sampleID were merged previously to set infile and PicardToolsMergeSwitch to enable correct handling of number of infiles to process
    my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles(\%{$scriptParameterHashRef}, \%sampleInfo, \%lane, \%infilesLaneNoEnding, $sampleID);
    
    if ($PicardToolsMergeSwitch == 1) {  #Files were merged previously

	## Copy file(s) to temporary directory
	print $FILEHANDLE "## Copy file(s) to temporary directory\n"; 
	&MigrateFileToTemp($FILEHANDLE, 
			    {'path' => $inSampleDirectory."/".$infile.$infileEnding.".b*",
			     'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
			    });
	print $FILEHANDLE "wait", "\n\n";

	## GATK ReAlignerTargetCreator
	print $FILEHANDLE "## GATK ReAlignerTargetCreator","\n";

	## Writes java core commands to filehandle.
	&JavaCore($FILEHANDLE,
		  {'memoryAllocation' => "Xmx24g",
		   'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
		   'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		  });

	print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	print $FILEHANDLE "-l INFO ";  #Set the minimum level of logging
	print $FILEHANDLE "-T RealignerTargetCreator ";  #Type of analysis to run
	print $FILEHANDLE "-R ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference file 
	print $FILEHANDLE "-known ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'GATKReAlignerINDELKnownSet1'}." ";  #Input VCF file with known indels
	print $FILEHANDLE "-known ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'GATKReAlignerINDELKnownSet2'}." ";  #Input VCF file with known indels
	print $FILEHANDLE "-nt ".$nrCores." ";  #How many data threads should be allocated to running this analysis.
	print $FILEHANDLE "-dcov ".${$scriptParameterHashRef}{'GATKDownSampleToCoverage'}." ";  #Coverage to downsample to at any given locus
	print $FILEHANDLE "-I ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$infileEnding.".bam ";  #InFile	    
	print $FILEHANDLE "-o ".$intermediarySampleDirectory."/".$infile.$outfileEnding.".intervals ", "\n\n";  #Interval outFile
	
	## GATK IndelRealigner
	print $FILEHANDLE "## GATK IndelRealigner","\n";
	
	## Writes java core commands to filehandle.
	&JavaCore($FILEHANDLE,
		  {'memoryAllocation' => "Xmx24g",
		   'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
		   'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		  });
	
	print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	print $FILEHANDLE "-l INFO ";
	print $FILEHANDLE "-T IndelRealigner ";
	print $FILEHANDLE "-R ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference file
	print $FILEHANDLE "-known ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'GATKReAlignerINDELKnownSet1'}." ";  #Input VCF file with known indels
	print $FILEHANDLE "-known ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'GATKReAlignerINDELKnownSet2'}." ";  #Input VCF file with known indels	 
	print $FILEHANDLE "-dcov ".${$scriptParameterHashRef}{'GATKDownSampleToCoverage'}." ";  #Coverage to downsample to at any given locus
	print $FILEHANDLE "-I ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$infileEnding.".bam ";  #InFile	
	print $FILEHANDLE "-o ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.".bam ";  #OutFile
	print $FILEHANDLE "-targetIntervals ".$intermediarySampleDirectory."/".$infile.$outfileEnding.".intervals ", "\n\n";

	## Copies file from temporary directory.
	print $FILEHANDLE "## Copy file from temporary directory\n"; 
	&MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.".b*", $outSampleDirectory."/", $FILEHANDLE);
	print $FILEHANDLE "wait", "\n\n";

	if ( (${$scriptParameterHashRef}{'pGATKRealigner'} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {
	    
	    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'} = $outSampleDirectory."/".$infile.$outfileEnding.".bam";	    
	}	
    }
    else  {  #No previous merge

	## Copies files from source to temporary folder. Loop over files specified by $arrayRef and collects files from $extractArrayRef.
	&MigrateFilesToTemp(\%{$scriptParameterHashRef}, \@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, \@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, $FILEHANDLE,
			    {'inSampleDirectory' => $inSampleDirectory,
			     'nrCores' => $nrCores,
			     'fileEnding' => $infileEnding.".b*"
			    });
	
	for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} });$infileCounter++) {  #For all infiles per lane
	    
	    my $infile = ${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter];

	    ## GATK ReAlignerTargetCreator
	    print $FILEHANDLE "## GATK ReAlignerTargetCreator","\n";

	    ## Writes java core commands to filehandle.
	    &JavaCore($FILEHANDLE,
		      {'memoryAllocation' => "Xmx24g",
		       'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
		       'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		      });
	    
	    print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	    print $FILEHANDLE "-l INFO ";  #Set the minimum level of logging
	    print $FILEHANDLE "-T RealignerTargetCreator ";  #Type of analysis to run
	    print $FILEHANDLE "-R ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference file 
	    print $FILEHANDLE "-known ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'GATKReAlignerINDELKnownSet1'}." ";  #Input VCF file with known indels
	    print $FILEHANDLE "-known ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'GATKReAlignerINDELKnownSet2'}." ";  #Input VCF file with known indels
	    print $FILEHANDLE "-nt ".$nrCores." ";  #How many data threads should be allocated to running this analysis.
	    print $FILEHANDLE "-dcov ".${$scriptParameterHashRef}{'GATKDownSampleToCoverage'}." ";  #Coverage to downsample to at any given locus	 
	    print $FILEHANDLE "-I ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$infileEnding.".bam ";  #InFile
	    print $FILEHANDLE "-o ".$intermediarySampleDirectory."/".$infile.$outfileEnding.".intervals ", "\n\n";  #Interval outFile
	    
	    ## GATK IndelRealigner
	    print $FILEHANDLE "## GATK IndelRealigner","\n";

	    ## Writes java core commands to filehandle.
	    &JavaCore($FILEHANDLE,
		      {'memoryAllocation' => "Xmx24g",
		       'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
		       'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		      });

	    print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
	    print $FILEHANDLE "-l INFO ";
	    print $FILEHANDLE "-T IndelRealigner ";
	    print $FILEHANDLE "-R ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference file
	    print $FILEHANDLE "-known ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'GATKReAlignerINDELKnownSet1'}." ";  #Input VCF file with known indels
	    print $FILEHANDLE "-known ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'GATKReAlignerINDELKnownSet2'}." ";  #Input VCF file with known indels
	    print $FILEHANDLE "-dcov ".${$scriptParameterHashRef}{'GATKDownSampleToCoverage'}." ";  #Coverage to downsample to at any given locus
	    print $FILEHANDLE "-I ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$infileEnding.".bam ";  #InFile		
	    print $FILEHANDLE "-o ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.".bam ";  #OutFile
	    print $FILEHANDLE "-targetIntervals ".$intermediarySampleDirectory."/".$infile.$outfileEnding.".intervals ", "\n\n";

	    if ( (${$scriptParameterHashRef}{'pGATKRealigner'} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {
		
		${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'} = $outSampleDirectory."/".$infile.$outfileEnding.".bam";	    
	    }
	}

	## Copies files from temporary folder to source. Loop over files specified by $arrayRef and collects files from $extractArrayRef.
	&MigrateFilesFromTemp(\@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, \@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, $outSampleDirectory, ${$scriptParameterHashRef}{'tempDirectory'}, $nrCores, $outfileEnding.".b*", $FILEHANDLE);
    }
    
    &RemoveDirectory(\${$scriptParameterHashRef}{'tempDirectory'}, $FILEHANDLE);

    close($FILEHANDLE);

    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'sampleID' => $sampleID,
			'dependencies' => 1, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}


sub RCoveragePlots { 

##RCoveragePlots
    
##Function : Generates sbatch scripts for R scripts: 1. covplots_genome.R 2. covplots_exome.R; on files generated from calculateCoverage genomeCoverageBED.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $infilesLaneNoEndingHashRef, $sampleID, $aligner, $programName
##         : $parameterHashRef           => The parameter hash {REF}
##         : $scriptParameterHashRef     => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef          => Info on samples and family hash {REF}
##         : $infilesLaneNoEndingHashRef => The infile(s) without the ".ending" {REF}
##         : $sampleID                   => The sampleID
##         : $aligner                    => The aligner used in the analysis
##         : $programName                => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $infilesLaneNoEndingHashRef = $_[3];
    my $sampleID = $_[4];
    my $aligner = $_[5];
    my $programName = $_[6];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					   {'directoryID' => $sampleID,
					    'programName' => $programName,
					    'programDirectory' => lc($aligner."/coveragereport"),
					    'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
					   });
    
    ## Assign directories
    my $inSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner."/coveragereport";
    my $outSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner."/coveragereport";

    my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pPicardToolsMarkduplicates'}{'fileEnding'};

    ## Check if any files for this sampleID were merged previously to set infile and PicardToolsMergeSwitch to enable correct handling of number of infiles to process
    my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles(\%{$scriptParameterHashRef}, \%sampleInfo, \%lane, \%infilesLaneNoEnding, $sampleID);    
    
    if ($PicardToolsMergeSwitch == 1) {  #Files were merged previously
	
	if ( defined(${$scriptParameterHashRef}{'pGenomeCoverageBED'}) && (${$scriptParameterHashRef}{'pGenomeCoverageBED'} > 0) ) {

	    my $infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pGenomeCoverageBED'}{'fileEnding'};

	    print $FILEHANDLE "Rscript ";
	    print $FILEHANDLE ${$scriptParameterHashRef}{'inScriptDir'}."/covplots_genome.R ";
	    print $FILEHANDLE $inSampleDirectory."/".$infile.$infileEnding." ";  #InFile
	    print $FILEHANDLE $infile." ";  #Sample name
	    print $FILEHANDLE ${$scriptParameterHashRef}{'GenomeCoverageBEDMaxCoverage'}." ";  #X-axis max scale
	    print $FILEHANDLE $outSampleDirectory, " &","\n\n";  #OutFile
	}
    }
    else {  #No previous merge

	for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} });$infileCounter++) {  #For all infiles per lane
	    
	    my $infile = ${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter];
	    
	    if ( defined(${$scriptParameterHashRef}{'pGenomeCoverageBED'}) && (${$scriptParameterHashRef}{'pGenomeCoverageBED'} > 0) ) {

		my $infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pGenomeCoverageBED'}{'fileEnding'};
		
		print $FILEHANDLE "Rscript ";
		print $FILEHANDLE ${$scriptParameterHashRef}{'inScriptDir'}."/covplots_genome.R ";
		print $FILEHANDLE $inSampleDirectory."/".$infile.$infileEnding." ";  #InFile
		print $FILEHANDLE $infile." ";  #Sample name
		print $FILEHANDLE ${$scriptParameterHashRef}{'GenomeCoverageBEDMaxCoverage'}." ";  #X-axis max scale
		print $FILEHANDLE $outSampleDirectory, " &", "\n\n";  #OutFile
	    }
	}
    }
    print $FILEHANDLE "wait", "\n\n";
    close($FILEHANDLE);

    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {
	
	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'sampleID' => $sampleID,
			'dependencies' => 2, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
    return;
}


sub GenomeCoverageBED { 

##GenomeCoverageBED
    
##Function : Calculates coverage on BAM files.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $infilesLaneNoEndingHashRef, $laneHashRef, $sampleID, $aligner, $programName
##         : $parameterHashRef           => The parameter hash {REF}
##         : $scriptParameterHashRef     => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef          => Info on samples and family hash {REF}
##         : $infilesLaneNoEndingHashRef => The infile(s) without the ".ending" {REF}
##         : $laneHashRef                => The lane info hash {REF}
##         : $sampleID                   => The sampleID
##         : $aligner                    => The aligner used in the analysis
##         : $programName                => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $infilesLaneNoEndingHashRef = $_[3];
    my $laneHashRef = $_[4];
    my $sampleID = $_[5];
    my $aligner = $_[6];
    my $programName = $_[7];
    
    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $fileName;
    
    ## Assign directories
    my $inSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner;
    my $outSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner."/coveragereport";

    my $infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pPicardToolsMarkduplicates'}{'fileEnding'};
    my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{"p".$programName}{'fileEnding'};

    ## Check if any files for this sampleID were merged previously to set infile and PicardToolsMergeSwitch to enable correct handling of number of infiles to process
    my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles(\%{$scriptParameterHashRef}, \%sampleInfo, \%lane, \%infilesLaneNoEnding, $sampleID);
    my $coreCounter=1;
    my $nrCores=1;

    if ($PicardToolsMergeSwitch == 1) {  #Files were merged previously
	
	## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
	($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					    {'directoryID' => $sampleID,
					     'programName' => $programName,
					     'programDirectory' => lc($aligner."/coveragereport"),
					     'nrofCores' => $nrCores,
					     'processTime' => 4,
					     'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
					    });
	
	## Copy file(s) to temporary directory
	print $FILEHANDLE "## Copy file(s) to temporary directory\n";
	&MigrateFileToTemp($FILEHANDLE, 
			    {'path' => $inSampleDirectory."/".$infile.$infileEnding.".b*",
			     'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
			    });
	print $FILEHANDLE "wait", "\n\n";
	
	## GenomeCoverageBed
	print $FILEHANDLE "## Calculate coverage metrics on alignment\n";
	print $FILEHANDLE "genomeCoverageBed ";
	print $FILEHANDLE "-max ".${$scriptParameterHashRef}{'GenomeCoverageBEDMaxCoverage'}." ";  #Combine all positions with a depth >= max into a single bin in the histogram.
	print $FILEHANDLE "-ibam ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$infileEnding.".bam ";  #InFile
	print $FILEHANDLE "> ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding." ", "\n\n";  #OutFile
	
	## Copies file from temporary directory.
	print $FILEHANDLE "## Copy file from temporary directory\n"; 
	&MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding, $outSampleDirectory."/", $FILEHANDLE);
	print $FILEHANDLE "wait", "\n\n";
    }
    else {  #No merged files
	
	## Set the number of cores to allocate per sbatch job.
	my $nrCores = &NrofCoresPerSbatch(\%{$scriptParameterHashRef}, scalar( @{${$laneHashRef}{$sampleID}} ) );  #Detect the number of cores to from lanes	
	
	## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
	($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					    {'directoryID' => $sampleID,
					     'programName' => $programName,
					     'programDirectory' => lc($aligner."/coveragereport"),
					     'nrofCores' => $nrCores,
					     'processTime' => 4,
					     'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
					    });
	
	## Copies files from source to temporary folder. Loop over files specified by $arrayRef and collects files from $extractArrayRef.
	&MigrateFilesToTemp(\%{$scriptParameterHashRef}, \@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, \@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, $FILEHANDLE,
			     {'inSampleDirectory' => $inSampleDirectory,
			      'nrCores' => $nrCores,
			      'fileEnding' => $infileEnding.".b*"});
	
	## GenomeCoverageBed
	print $FILEHANDLE "## Calculate coverage metrics on alignment\n";
	for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} });$infileCounter++) {  #For all files from MosaikAlign or BWA_Sampe

	    &PrintWait(\$infileCounter, \$nrCores, \$coreCounter, $FILEHANDLE);	    

	    my $infile = ${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter];	    
	    
	    print $FILEHANDLE "genomeCoverageBed ";
	    print $FILEHANDLE "-max ".${$scriptParameterHashRef}{'GenomeCoverageBEDMaxCoverage'}." ";  #Combine all positions with a depth >= max into a single bin in the histogram.
	    print $FILEHANDLE "-ibam ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$infileEnding.".bam ";  #InFile
	    print $FILEHANDLE "> ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding." &", "\n\n";  #outFile
	}
	print $FILEHANDLE "wait", "\n\n";

	## Copies files from temporary folder to source. Loop over files specified by $arrayRef and collects files from $extractArrayRef.
	&MigrateFilesFromTemp(\@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, \@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, $outSampleDirectory, ${$scriptParameterHashRef}{'tempDirectory'}, $nrCores, $outfileEnding, $FILEHANDLE);
    }
    
    &RemoveDirectory(\${$scriptParameterHashRef}{'tempDirectory'}, $FILEHANDLE);

    close($FILEHANDLE);

    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'sampleID' => $sampleID,
			'dependencies' => 1, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
    return;
}


sub PicardToolsCalculateHSMetrics { 
 
##PicardToolsCalculateHSMetrics
    
##Function : Calculates coverage on exonic part of BAM files.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $infilesLaneNoEndingHashRef, $laneHashRef, $sampleID, $aligner, $programName
##         : $parameterHashRef           => The parameter hash {REF}
##         : $scriptParameterHashRef     => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef          => Info on samples and family hash {REF}
##         : $infilesLaneNoEndingHashRef => The infile(s) without the ".ending" {REF}
##         : $laneHashRef                => The lane info hash {REF}
##         : $sampleID                   => The sampleID
##         : $aligner                    => The aligner used in the analysis
##         : $programName                => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $infilesLaneNoEndingHashRef = $_[3];
    my $laneHashRef = $_[4];
    my $sampleID = $_[5];
    my $aligner = $_[6];
    my $programName = $_[7];
    
    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $fileName;

    ## Assign directories
    my $inSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner;
    my $outSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner."/coveragereport";

    my $infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pPicardToolsMarkduplicates'}{'fileEnding'};
    my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pPicardToolsMarkduplicates'}{'fileEnding'};

    ## Check if any files for this sampleID were merged previously to set infile and PicardToolsMergeSwitch to enable correct handling of number of infiles to process
    my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles(\%{$scriptParameterHashRef}, \%sampleInfo, \%lane, \%infilesLaneNoEnding, $sampleID);
    my $coreCounter=1;
    my $nrCores=1;

    if ($PicardToolsMergeSwitch == 1) {  #Files were merged previously
	
	## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
	($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					    {'directoryID' => $sampleID,
					     'programName' => $programName,
					     'programDirectory' => lc($aligner."/coveragereport"),
					     'nrofCores' => $nrCores,
					     'processTime' => 4,
					     'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
					    });
	
	## Copy file(s) to temporary directory
	print $FILEHANDLE "## Copy file(s) to temporary directory\n"; 
	&MigrateFileToTemp($FILEHANDLE, 
			    {'path' => $inSampleDirectory."/".$infile.$infileEnding.".b*",
			     'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
			    });
	print $FILEHANDLE "wait", "\n\n";

	## CalculateHsMetrics
	print $FILEHANDLE "## Calculate capture metrics on alignment\n";

	## Writes java core commands to filehandle.
	&JavaCore($FILEHANDLE,
		  {'memoryAllocation' => "Xmx4g",
		   'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
		   'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		  });

	print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'picardToolsPath'}."/CalculateHsMetrics.jar ";
	print $FILEHANDLE "INPUT=".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$infileEnding.".bam ";  #InFile
	print $FILEHANDLE "OUTPUT=".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding."_CalculateHsMetrics ";  #OutFile
	print $FILEHANDLE "REFERENCE_SEQUENCE=".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference file
	print $FILEHANDLE "BAIT_INTERVALS=".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'exomeTargetPaddedBedInfileLists'}." ";  #Capture kit padded target infile_list file
	print $FILEHANDLE "TARGET_INTERVALS=".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'exomeTargetBedInfileLists'}, "\n\n";  #Capture kit target infile_list file
	
	## Copies file from temporary directory.
	print $FILEHANDLE "## Copy file from temporary directory\n"; 
	&MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding."_CalculateHsMetrics", $outSampleDirectory, $FILEHANDLE);
	print $FILEHANDLE "wait", "\n\n";

	if ( (${$scriptParameterHashRef}{'pPicardToolsCalculateHSMetrics'} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	    ## Collect QC metadata info for later use
	    &SampleInfoQC(\%sampleInfo,
			  {'familyID' => ${$scriptParameterHashRef}{'familyID'},
			   'sampleID' => $sampleID,
			   'programName' => "CalculateHsMetrics",
			   'infile' => $infile,
			   'outDirectory' => $outSampleDirectory,
			   'outFileEnding' => $outfileEnding."_CalculateHsMetrics",
			   'outDataType' => "infileDependent"
			  });
	}
    }
    else {  #No merged files
	
	## Set the number of cores to allocate per sbatch job.
	$nrCores = &NrofCoresPerSbatch(\%{$scriptParameterHashRef}, scalar( @{${$laneHashRef}{$sampleID}} ) );  #Detect the number of cores to from lanes	
	
	## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
	($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					    {'directoryID' => $sampleID,
					     'programName' => $programName,
					     'programDirectory' => lc($aligner."/coveragereport"),
					     'nrofCores' => $nrCores,
					     'processTime' => 4,
					     'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
					    });

	## Copies files from source to temporary folder. Loop over files specified by $arrayRef and collects files from $extractArrayRef.
	&MigrateFilesToTemp(\%{$scriptParameterHashRef}, \@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, \@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, $FILEHANDLE,
			     {'inSampleDirectory' => $inSampleDirectory,
			      'nrCores' => $nrCores,
			      'fileEnding' => $infileEnding.".b*"});

	## CalculateHsMetrics
	print $FILEHANDLE "## Calculate capture metrics on alignment\n";
	for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} });$infileCounter++) {  #For all files from MosaikAlign or BWA_Sampe

	    &PrintWait(\$infileCounter, \$nrCores, \$coreCounter, $FILEHANDLE);	    
	    
	    my $infile = ${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter];	    
	    
	    ## Writes java core commands to filehandle.
	    &JavaCore($FILEHANDLE,
		      {'memoryAllocation' => "Xmx4g",
		       'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
		       'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		      });

	    print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'picardToolsPath'}."/CalculateHsMetrics.jar ";
	    print $FILEHANDLE "INPUT=".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$infileEnding.".bam ";  #InFile
	    print $FILEHANDLE "OUTPUT=".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding."_CalculateHsMetrics ";  #OutFile
	    print $FILEHANDLE "REFERENCE_SEQUENCE=".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference file
	    print $FILEHANDLE "BAIT_INTERVALS=".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'exomeTargetPaddedBedInfileLists'}." ";  #Capture kit padded target infile_list file
	    print $FILEHANDLE "TARGET_INTERVALS=".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'exomeTargetBedInfileLists'}." &", "\n\n";  #Capture kit target infile_list file 

	    if ( (${$scriptParameterHashRef}{'pPicardToolsCalculateHSMetrics'} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

		## Collect QC metadata info for later use                                                                                                 
		&SampleInfoQC(\%sampleInfo,
			       {'familyID' => ${$scriptParameterHashRef}{'familyID'},
				'sampleID' => $sampleID,
				'programName' => "CalculateHsMetrics",
				'infile' => $infile,
				'outDirectory' => $outSampleDirectory,
				'outFileEnding' => $outfileEnding."_CalculateHsMetrics",
				'outDataType' => "infileDependent"
			       });
	    }
	}
	print $FILEHANDLE "wait", "\n\n";

	## Copies files from temporary folder to source. Loop over files specified by $arrayRef and collects files from $extractArrayRef.
	&MigrateFilesFromTemp(\@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, \@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, $outSampleDirectory, ${$scriptParameterHashRef}{'tempDirectory'}, $nrCores, $outfileEnding."_CalculateHsMetrics", $FILEHANDLE);
    }
    
    &RemoveDirectory(\${$scriptParameterHashRef}{'tempDirectory'}, $FILEHANDLE);

    close($FILEHANDLE);

    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'sampleID' => $sampleID,
			'dependencies' => 1, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}


sub PicardToolsCollectMultipleMetrics { 
 
##PicardToolsCollectMultipleMetrics
    
##Function : Calculates coverage and alignment metrics on BAM files.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $infilesLaneNoEndingHashRef, $laneHashRef, $sampleID, $aligner, $programName
##         : $parameterHashRef           => The parameter hash {REF}
##         : $scriptParameterHashRef     => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef          => Info on samples and family hash {REF}
##         : $infilesLaneNoEndingHashRef => The infile(s) without the ".ending" {REF}
##         : $laneHashRef                => The lane info hash {REF}
##         : $sampleID                   => The sampleID
##         : $aligner                    => The aligner used in the analysis
##         : $programName                => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $infilesLaneNoEndingHashRef = $_[3]; 
    my $laneHashRef = $_[4];
    my $sampleID = $_[5];
    my $aligner = $_[6];
    my $programName = $_[7];
    
    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $fileName;

    ## Assign directories
    my $inSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner;
    my $outSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner."/coveragereport";

    my $infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pPicardToolsMarkduplicates'}{'fileEnding'};
    my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pPicardToolsMarkduplicates'}{'fileEnding'};

    ## Check if any files for this sampleID were merged previously to set infile and PicardToolsMergeSwitch to enable correct handling of number of infiles to process
    my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles(\%{$scriptParameterHashRef}, \%sampleInfo, \%lane, \%infilesLaneNoEnding, $sampleID);
    my $nrCores = 1;
    my $coreCounter=1;

    if ($PicardToolsMergeSwitch == 1) {  #Files were merged previously
	
	## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
	($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					    {'directoryID' => $sampleID,
					     'programName' => $programName,
					     'programDirectory' => lc($aligner."/coveragereport"),
					     'nrofCores' => $nrCores,
					     'processTime' => 4,
					     'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
					    });

	## Copy file(s) to temporary directory
	print $FILEHANDLE "## Copy file(s) to temporary directory\n"; 
	&MigrateFileToTemp($FILEHANDLE, 
			    {'path' => $inSampleDirectory."/".$infile.$infileEnding.".b*",
			     'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
			    });
	print $FILEHANDLE "wait", "\n\n";

	## CollectMultipleMetrics
	print $FILEHANDLE "## Collecting multiple metrics on alignment\n";

	## Writes java core commands to filehandle.
	&JavaCore($FILEHANDLE,
		  {'memoryAllocation' => "Xmx4g",
		   'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
		   'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		  });

	print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'picardToolsPath'}."/CollectMultipleMetrics.jar ";
	print $FILEHANDLE "INPUT=".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$infileEnding.".bam ";  #InFile
	print $FILEHANDLE "OUTPUT=".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding." ";  #OutFile
	print $FILEHANDLE "R=".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}, "\n\n";  #Reference file

	## Copies file from temporary directory.
	print $FILEHANDLE "## Copy file from temporary directory\n";
	&MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.".alignment_summary_metrics", $outSampleDirectory."/".$infile.$outfileEnding.".alignment_summary_metrics", $FILEHANDLE);
	&MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.".quality*", $outSampleDirectory, $FILEHANDLE);
	&MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.".insert*", $outSampleDirectory, $FILEHANDLE);
	print $FILEHANDLE "wait", "\n\n";

	if ( (${$scriptParameterHashRef}{'pPicardToolsCollectMultipleMetrics'} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	    ## Collect QC metadata info for later use                                                                                             
	    &SampleInfoQC(\%sampleInfo,
			  {'familyID' => ${$scriptParameterHashRef}{'familyID'},
			   'sampleID' => $sampleID,
			   'programName' => "CollectMultipleMetrics",
			   'infile' => $infile,
			   'outDirectory' => $outSampleDirectory,
			   'outFileEnding' => $outfileEnding.".alignment_summary_metrics",
			   'outDataType' => "infileDependent"
			  });
	}
    }
    else {  #No merged files
	
	## Set the number of cores to allocate per sbatch job.
	my $nrCores = &NrofCoresPerSbatch(\%{$scriptParameterHashRef}, scalar( @{${$laneHashRef}{$sampleID}} ) );  #Detect the number of cores to from lanes	
	
	## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
	($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					    {'directoryID' => $sampleID,
					     'programName' => $programName,
					     'programDirectory' => lc($aligner."/coveragereport"),
					     'nrofCores' => $nrCores,
					     'processTime' => 4,
					     'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
					    });

	## Copies files from source to temporary folder. Loop over files specified by $arrayRef and collects files from $extractArrayRef.
	&MigrateFilesToTemp(\%{$scriptParameterHashRef}, \@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, \@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, $FILEHANDLE,
			    {'inSampleDirectory' => $inSampleDirectory,
			     'nrCores' => $nrCores,
			     'fileEnding' => $infileEnding.".b*"
			    });

	## CollectMultipleMetrics
	print $FILEHANDLE "## Collecting multiple metrics on alignment\n";
	for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} });$infileCounter++) {  #For all files from MosaikAlign or BWA_Sampe

	    &PrintWait(\$infileCounter, \$nrCores, \$coreCounter, $FILEHANDLE);	

	    my $infile = ${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter];	    
	   
	    ## Writes java core commands to filehandle.
	    &JavaCore($FILEHANDLE,
		      {'memoryAllocation' => "Xmx4g",
		       'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
		       'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		      }); 

	    print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'picardToolsPath'}."/CollectMultipleMetrics.jar ";
	    print $FILEHANDLE "INPUT=".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$infileEnding.".bam ";  #InFile
	    print $FILEHANDLE "OUTPUT=".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding." ";  #outFile
	    print $FILEHANDLE "R=".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." &", "\n\n";  #Reference file

	    if ( (${$scriptParameterHashRef}{'pPicardToolsCollectMultipleMetrics'} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

		## Collect QC metadata info for later use
		&SampleInfoQC(\%sampleInfo,
			      {'familyID' => ${$scriptParameterHashRef}{'familyID'},
			       'sampleID' => $sampleID,
			       'programName' => "CollectMultipleMetrics",
			       'infile' => $infile,
			       'outDirectory' => $outSampleDirectory,
			       'outFileEnding' => $outfileEnding.".alignment_summary_metrics",
			       'outDataType' => "infileDependent"
			      });
	    }
	}
	print $FILEHANDLE "wait", "\n\n";

	## Copies files from temporary folder to source. Loop over files specified by $arrayRef and collects files from $extractArrayRef.
	&MigrateFilesFromTemp(\@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, \@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, $outSampleDirectory, ${$scriptParameterHashRef}{'tempDirectory'}, $nrCores, $outfileEnding.".alignment_summary_metrics", $FILEHANDLE);
	&MigrateFilesFromTemp(\@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, \@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, $outSampleDirectory, ${$scriptParameterHashRef}{'tempDirectory'}, $nrCores, $outfileEnding.".quality*", $FILEHANDLE);
	&MigrateFilesFromTemp(\@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, \@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, $outSampleDirectory, ${$scriptParameterHashRef}{'tempDirectory'}, $nrCores, $outfileEnding.".insert*", $FILEHANDLE);
    }
    
    &RemoveDirectory(\${$scriptParameterHashRef}{'tempDirectory'}, $FILEHANDLE);

    close($FILEHANDLE);

    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'sampleID' => $sampleID,
			'dependencies' => 1, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}


sub ChanjoImport { 
 
##ChanjoImport
    
##Function : Loads the calculated coverage to family database
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $infilesLaneNoEndingHashRef, $familyID, $aligner, $programName
##         : $parameterHashRef           => The parameter hash {REF}
##         : $scriptParameterHashRef     => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef          => Info on samples and family hash {REF}
##         : $infilesLaneNoEndingHashRef => The infile(s) without the ".ending" {REF}
##         : $familyID                   => The familyID
##         : $aligner                    => The aligner used in the analysis
##         : $programName                => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $infilesLaneNoEndingHashRef = $_[3];
    my $familyID = $_[4];
    my $aligner = $_[5];
    my $programName = $_[6];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    
    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					   {'directoryID' => $familyID,
					    'programName' => $programName,
					    'programDirectory' => lc($programName),
					    'processTime' => 3,
					    'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
					   });

    my $outFamilyDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID;

    my $coreCounter=1;

    print $FILEHANDLE "workon ".${$scriptParameterHashRef}{'pythonVirtualEnvironment'}, "\n\n";  #Activate python environment
    
    ##Build family database for coverage report

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{${$scriptParameterHashRef}{'sampleIDs'}});$sampleIDCounter++) {   
	
	my $sampleID = ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter];
	my $inSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner."/coveragereport";
	my $infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pChanjoAnnotate'}{'fileEnding'};

	## Check if any files for this sampleID were merged previously to set infile and PicardToolsMergeSwitch to enable correct handling of number of infiles to process
	my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles(\%{$scriptParameterHashRef}, \%sampleInfo, \%lane, \%infilesLaneNoEnding, $sampleID);	

	if ($PicardToolsMergeSwitch == 1) {  #Files were merged previously
	    print $FILEHANDLE "chanjo ";
	    print $FILEHANDLE "import ";
	    print $FILEHANDLE $inSampleDirectory."/".$infile.$infileEnding.".bed ";
	    print $FILEHANDLE "--db=".$outFamilyDirectory."/".$familyID.".sqlite ","\n\n";  #Central Db for family
      	
	}
	else {
	    
	    for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} });$infileCounter++) {  #For all files from independent of merged or not

		&PrintWait(\$infileCounter, \${$scriptParameterHashRef}{'maximumCores'}, \$coreCounter, $FILEHANDLE);		
		
		my $infile = ${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter];
		print $FILEHANDLE "chanjo ";
		print $FILEHANDLE "import ";
		print $FILEHANDLE $inSampleDirectory."/".$infile.$infileEnding.".bed ";
		print $FILEHANDLE "--db=".$outFamilyDirectory."/".$familyID.".sqlite ","\n\n";  #Central Db for family
		
	    }
	}
    }
    print $FILEHANDLE "\n\ndeactivate ", "\n\n";  #Deactivate python environment

    close($FILEHANDLE); 

    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'dependencies' => 5, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}


sub ChanjoSexCheck {

##ChanjoSexCheck
    
##Function : Predict gender from BAM files.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $infilesLaneNoEndingHashRef, $laneHashRef, $sampleID, $aligner, $programName
##         : $parameterHashRef           => The parameter hash {REF}
##         : $scriptParameterHashRef     => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef          => Info on samples and family hash {REF}
##         : $infilesLaneNoEndingHashRef => The infile(s) without the ".ending" {REF}
##         : $laneHashRef                => The lane info hash {REF}
##         : $sampleID                   => The sampleID
##         : $aligner                    => The aligner used in the analysis
##         : $programName                => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $infilesLaneNoEndingHashRef = $_[3];
    my $laneHashRef = $_[4];
    my $sampleID = $_[5];
    my $aligner = $_[6];
    my $programName = $_[7];

    my $FILEHANDLE = IO::Handle->new();  #Create anonymous filehandle
    my $fileName;      

    ## Assign directories               
    my $outFamilyDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".${$scriptParameterHashRef}{'familyID'};
    my $inSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner;
    my $outSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner."/coveragereport";

    my $infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pPicardToolsMarkduplicates'}{'fileEnding'};
    my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{"p".$programName}{'fileEnding'};

    ## Check if any files for this sampleID were merged previously to set infile and PicardToolsMergeSwitch to enable correct handling of number of infiles to process
    my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles(\%{$scriptParameterHashRef}, \%sampleInfo, \%lane, \%infilesLaneNoEnding, $sampleID);
    my $coreCounter=1;	

    if ($PicardToolsMergeSwitch == 1) {  #Files were merged previously
	
	## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
	($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					    {'directoryID' => $sampleID,
					     'programName' => $programName,
					     'programDirectory' => lc($aligner."/coveragereport"),
					     'processTime' => 2,
					    });

	print $FILEHANDLE "workon ".${$scriptParameterHashRef}{'pythonVirtualEnvironment'}, "\n\n";  #Activate python environment

	## ChanjoSexCheck
	print $FILEHANDLE "## Predicting sex from alignment\n";
	print $FILEHANDLE "sex-check ";
	print $FILEHANDLE $inSampleDirectory."/".$infile.$infileEnding.".bam ";  #InFile
	print $FILEHANDLE "> ".$outSampleDirectory."/".$infile.$outfileEnding, "\n\n";  #OutFile
	
	if ( (${$scriptParameterHashRef}{'pChanjoSexCheck'} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	    ## Collect QC metadata info for later use
	    &SampleInfoQC(\%sampleInfo,
			  {'familyID' => ${$scriptParameterHashRef}{'familyID'},
			   'sampleID' => $sampleID,
			   'programName' => "ChanjoSexCheck",
			   'infile' => $infile,
			   'outDirectory' => $outSampleDirectory,
			   'outFileEnding' => $outfileEnding,
			   'outDataType' => "infileDependent"
			  });
	}
    }
    else {  #No merged files
	
	## Set the number of cores to allocate per sbatch job.
	my $nrCores = &NrofCoresPerSbatch(\%{$scriptParameterHashRef}, scalar( @{${$laneHashRef}{$sampleID}} ));  #Detect the number of cores to use from lanes
	
	## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
	($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					    {'directoryID' => $sampleID,
					     'programName' => $programName,
					     'programDirectory' => lc($aligner."/coveragereport"),
					     'nrofCores' => $nrCores,
					     'processTime' => 2,
					    });

	print $FILEHANDLE "workon ".${$scriptParameterHashRef}{'pythonVirtualEnvironment'}, "\n\n";  #Activate python environment

	## ChanjoSexCheck
	print $FILEHANDLE "## Predicting sex from alignment\n";
	for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} });$infileCounter++) {  #For all files from independent of merged or not
	    
	    &PrintWait(\$infileCounter, \${$scriptParameterHashRef}{'maximumCores'}, \$coreCounter, $FILEHANDLE);
	    
	    my $infile = ${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter];
	    
	    print $FILEHANDLE "sex-check ";
	    print $FILEHANDLE $inSampleDirectory."/".$infile.$infileEnding.".bam ";  #InFile
	    print $FILEHANDLE "> ".$outSampleDirectory."/".$infile.$outfileEnding." &", "\n\n";  #OutFile

	    if ( (${$scriptParameterHashRef}{'pChanjoSexCheck'} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {
		
		## Collect QC metadata info for later use
		&SampleInfoQC(\%sampleInfo,
			       {'familyID' => ${$scriptParameterHashRef}{'familyID'},
				'sampleID' => $sampleID,
				'programName' => "ChanjoSexCheck",
				'infile' => $infile,
				'outDirectory' => $outSampleDirectory,
				'outFileEnding' => $outfileEnding,
				'outDataType' => "infileDependent"
			       });
	    }
	}
	print $FILEHANDLE "wait", "\n\n";
    }
    print $FILEHANDLE "deactivate ", "\n\n";  #Deactivate python environment

    close($FILEHANDLE);

    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {
	
	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'sampleID' => $sampleID,
			'dependencies' => 1, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}


sub ChanjoAnnotate { 

##ChanjoAnnotate
    
##Function : Generate coverage bed outfile for each individual.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $infilesLaneNoEndingHashRef, $laneHashRef, $sampleID, $aligner, $programName
##         : $parameterHashRef           => The parameter hash {REF}
##         : $scriptParameterHashRef     => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef          => Info on samples and family hash {REF}
##         : $infilesLaneNoEndingHashRef => The infile(s) without the ".ending" {REF}
##         : $laneHashRef                => The lane info hash {REF}
##         : $sampleID                   => The sampleID
##         : $aligner                    => The aligner used in the analysis
##         : $programName                => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $infilesLaneNoEndingHashRef = $_[3];
    my $laneHashRef = $_[4];
    my $sampleID = $_[5];
    my $aligner = $_[6];
    my $programName = $_[7];

    my $FILEHANDLE = IO::Handle->new();  #Create anonymous filehandle
    my $fileName;
    my $nrCores = 1;      

    ## Assign directories
    my $outFamilyDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".${$scriptParameterHashRef}{'familyID'};
    my $inSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner;
    my $outSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner."/coveragereport";

    my $infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pPicardToolsMarkduplicates'}{'fileEnding'};
    my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{"p".$programName}{'fileEnding'};

    ## Check if any files for this sampleID were merged previously to set infile and PicardToolsMergeSwitch to enable correct handling of number of infiles to process
    my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles(\%{$scriptParameterHashRef}, \%sampleInfo, \%lane, \%infilesLaneNoEnding, $sampleID);
    my $coreCounter=1;	

    if ($PicardToolsMergeSwitch == 1) {  #Files were merged previously

	## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
	($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					    {'directoryID' => $sampleID,
					     'programName' => $programName,
					     'programDirectory' => lc($aligner."/coveragereport"),
					     'nrofCores' => $nrCores,
					     'processTime' => 2,
					     'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
					    });

	## Copy file(s) to temporary directory
	print $FILEHANDLE "## Copy file(s) to temporary directory\n"; 
	&MigrateFileToTemp($FILEHANDLE, 
			    {'path' => $inSampleDirectory."/".$infile.$infileEnding.".b*",
			     'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
			    });
	print $FILEHANDLE "wait", "\n\n";

	print $FILEHANDLE "workon ".${$scriptParameterHashRef}{'pythonVirtualEnvironment'}, "\n\n";  #Activate python environment

	## ChanjoAnnotate
	print $FILEHANDLE "## Annotating bed from alignment\n";
	&ChanjoConvert($FILEHANDLE, ${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'chanjoBuildDb'});
	print $FILEHANDLE "chanjo ";
	print $FILEHANDLE "annotate ";
	print $FILEHANDLE ${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$infileEnding.".bam ";  #InFile
	print $FILEHANDLE "--cutoff ".${$scriptParameterHashRef}{'chanjoAnnotateCutoff'}." ";  #The “cutoff” is used for the completeness calculation
	print $FILEHANDLE "--sample ".$sampleID." ";  #A unique sample Id
	print $FILEHANDLE "--extendby 2 ";  #Dynamically extend intervals symetrically
	print $FILEHANDLE "--group ".${$scriptParameterHashRef}{'familyID'}." ";  #Grouping option for samples
	print $FILEHANDLE "> ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.".bed". "\n\n";  #OutFile

	## Copies file from temporary directory.
	print $FILEHANDLE "## Copy file from temporary directory\n"; 
	&MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.".bed", $outSampleDirectory."/".$infile.$outfileEnding.".bed", $FILEHANDLE);
	
	if ( (${$scriptParameterHashRef}{'pChanjoAnnotate'} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	    ## Collect QC metadata info for later use
	    &SampleInfoQC(\%sampleInfo,
			  {'familyID' => ${$scriptParameterHashRef}{'familyID'},
			   'sampleID' => $sampleID,
			   'programName' => "ChanjoAnnotate",
			   'infile' => $infile,
			   'outDirectory' => $outSampleDirectory,
			   'outFileEnding' => $outfileEnding.".bed",
			   'outDataType' => "infileDependent"
			  });
	}
    }
    else {  #No merged files
	
	## Set the number of cores to allocate per sbatch job.
	my $nrCores = &NrofCoresPerSbatch(\%{$scriptParameterHashRef}, scalar( @{${$laneHashRef}{$sampleID}} ));  #Detect the number of cores to use from lanes

	## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
	($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					    {'directoryID' => $sampleID,
					     'programName' => $programName,
					     'programDirectory' => lc($aligner."/coveragereport"),
					     'nrofCores' => $nrCores,
					     'processTime' => 2,
					     'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
					    });	

	## Copies files from source to temporary folder. Loop over files specified by $arrayRef and collects files from $extractArrayRef.
	&MigrateFilesToTemp(\%{$scriptParameterHashRef}, \@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, \@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, $FILEHANDLE,
			     {'inSampleDirectory' => $inSampleDirectory,
			      'nrCores' => $nrCores,
			      'fileEnding' => $infileEnding.".b*"});

	print $FILEHANDLE "workon ".${$scriptParameterHashRef}{'pythonVirtualEnvironment'}, "\n\n";  #Activate python environment

	## ChanjoAnnotate
	print $FILEHANDLE "## Annotating bed from alignment\n";
	for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} });$infileCounter++) {  #For all files from independent of merged or not
	    
	    &PrintWait(\$infileCounter, \$nrCores, \$coreCounter, $FILEHANDLE);
	    
	    my $infile = ${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter];
	    
	    &ChanjoConvert($FILEHANDLE, ${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'chanjoBuildDb'});
	    print $FILEHANDLE "chanjo ";
	    print $FILEHANDLE "annotate ";
	    print $FILEHANDLE $inSampleDirectory."/".$infile.$infileEnding.".bam ";  #InFile
	    print $FILEHANDLE "--cutoff ".${$scriptParameterHashRef}{'chanjoAnnotateCutoff'}." ";  #The “cutoff” is used for the completeness calculation
	    print $FILEHANDLE "--sample ".$sampleID." ";  #A unique sample Id
	    print $FILEHANDLE "--extendby 2 ";  #Dynamically extend intervals symetrically
	    print $FILEHANDLE "--group ".${$scriptParameterHashRef}{'familyID'}." ";  #Grouping option for samples
	    print $FILEHANDLE "> ".$outSampleDirectory."/".$infile.$outfileEnding.".bed &". "\n\n";  #OutFile

	    if ( (${$scriptParameterHashRef}{'pChanjoAnnotate'} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

		## Collect QC metadata info for later use
		&SampleInfoQC(\%sampleInfo,
			      {'familyID' => ${$scriptParameterHashRef}{'familyID'},
			       'sampleID' => $sampleID,
			       'programName' => "ChanjoAnnotate",
			       'infile' => $infile,
			       'outDirectory' => $outSampleDirectory,
			       'outFileEnding' => $outfileEnding.".bed",
			       'outDataType' => "infileDependent"
			      });
	    }
	}
	print $FILEHANDLE "wait", "\n\n";

	## Copies files from temporary folder to source. Loop over files specified by $arrayRef and collects files from $extractArrayRef.
	&MigrateFilesFromTemp(\@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, \@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, $outSampleDirectory, ${$scriptParameterHashRef}{'tempDirectory'}, $nrCores, $outfileEnding.".bed", $FILEHANDLE);
    }
    print $FILEHANDLE "deactivate ", "\n\n";  #Deactivate python environment

    &RemoveDirectory(\${$scriptParameterHashRef}{'tempDirectory'}, $FILEHANDLE);

    close($FILEHANDLE);

    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {
	
	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'sampleID' => $sampleID,
			'dependencies' => 5, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}


sub ChanjoBuild { 

##ChanjoBuild
    
##Function : Build database for downstream coverage analsysis.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $familyID, $programName
##         : $parameterHashRef       => The parameter hash {REF}
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef      => Info on samples and family hash {REF}
##         : $familyID               => The familyID
##         : $programName            => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $familyID = $_[3];
    my $programName = $_[4];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    
    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					   {'directoryID' => $familyID,
					    'programName' => $programName,
					    'programDirectory' => lc($programName),
					   });

    ## Assign directories
    my $outFamilyDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$familyID;

    print $FILEHANDLE "workon ".${$scriptParameterHashRef}{'pythonVirtualEnvironment'}, "\n\n";  #Activate python environment
    
    ## ChanjoBuild
    print $FILEHANDLE "## Build new coverage database\n";
    &ChanjoConvert($FILEHANDLE, ${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'chanjoBuildDb'});
    print $FILEHANDLE "chanjo ";
    print $FILEHANDLE "--db ".$outFamilyDirectory."/".$familyID.".sqlite ";  #Path/URI of the SQL database
    print $FILEHANDLE "--dialect sqlite ";  #Type of SQL database
    print $FILEHANDLE "build ";  #Chanjo sub program argument
    print $FILEHANDLE "--force", "\n\n";  #Overwrite existing assets without warning

    print $FILEHANDLE "deactivate ", "\n\n";  #Deactivate python environment

    close($FILEHANDLE); 

    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	## Collect QC metadata info for later use
	&SampleInfoQC(\%{$sampleInfoHashRef},
		      {'familyID' => ${$scriptParameterHashRef}{'familyID'},
		       'programName' => "ChanjoBuild",
		       'outDirectory' => $outFamilyDirectory,
		       'outFileEnding' => $familyID.".sqlite",
		       'outDataType' => "infileDependent"
		      });

	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'dependencies' => 5, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}


sub PicardToolsMarkDuplicates { 

##PicardToolsMarkDuplicates
    
##Function : Mark duplicated reads using PicardTools MarkDuplicates in files generated from alignment (sorted, merged).
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $infilesLaneNoEndingHashRef, $laneHashRef, $sampleID, $aligner, $programName
##         : $parameterHashRef                  => The parameter hash {REF}
##         : $scriptParameterHashRef            => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef                 => Info on samples and family hash {REF}
##         : $infilesLaneNoEndingHashRef        => The infile(s) without the ".ending" {REF}
##         : $infilesBothStrandsNoEndingHashRef => The infile(s) without the ".ending" and strand info {REF}
##         : $laneHashRef                       => The lane info hash {REF}
##         : $sampleID                          => The sampleID
##         : $aligner                           => The aligner used in the analysis
##         : $programName                       => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $infilesLaneNoEndingHashRef = $_[3];
    my $infilesBothStrandsNoEndingHashRef = $_[4];
    my $laneHashRef = $_[5];
    my $sampleID = $_[6];
    my $aligner = $_[7];
    my $programName = $_[8];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $fileName;
    my $time;
    
    if (${$scriptParameterHashRef}{'pPicardToolsMergeSamFiles'} eq 0) {  #If No merge has been performed then time requirements goes down

	$time = 3;
    }
    else {

	$time = ceil(3*scalar( @{ ${$infilesBothStrandsNoEndingHashRef}{$sampleID} }));  #One full lane on Hiseq takes approx. 3 h to process, round up to nearest full hour.	
    }
    
    ## Assign directories
    my $inSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner;
    my $outSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner;

    my $infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pPicardToolsMergeSamFiles'}{'fileEnding'};
    my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{"p".$programName}{'fileEnding'};

    ## Check if any files for this sampleID were merged previously to set infile and PicardToolsMergeSwitch to enable correct handling of number of infiles to process
    my ($infile, $PicardToolsMergeSwitch) = &CheckIfMergedFiles(\%{$scriptParameterHashRef}, \%sampleInfo, \%lane, \%infilesLaneNoEnding, $sampleID);
    my $coreCounter=1;
    
    if ($PicardToolsMergeSwitch == 1) {  #Files were merged previously

	my $nrCores = 1;
	## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
	($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					    {'directoryID' => $sampleID,
					     'programName' => $programName,
					     'programDirectory' => lc($aligner),
					     'nrofCores' => $nrCores,
					     'processTime' => $time,
					     'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
					    });
	
	## Copy file(s) to temporary directory
	print $FILEHANDLE "## Copy file(s) to temporary directory\n"; 
	&MigrateFileToTemp($FILEHANDLE, 
			    {'path' => $inSampleDirectory."/".$infile.$infileEnding.".b*",
			     'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
			    });
	print $FILEHANDLE "wait", "\n\n";

	## PicardToolsMarkDuplicates
	print $FILEHANDLE "## Marking Duplicates\n";

	## Writes java core commands to filehandle.
	&JavaCore($FILEHANDLE,
		  {'memoryAllocation' => "Xmx4g",
		   'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
		   'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		  });
	
	print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'picardToolsPath'}."/MarkDuplicates.jar ";
	print $FILEHANDLE "TMP_DIR=".${$scriptParameterHashRef}{'tempDirectory'}." ";  #Temp Directory
	print $FILEHANDLE "ASSUME_SORTED=true ";
	print $FILEHANDLE "CREATE_INDEX=TRUE ";  #Create a BAM index when writing a coordinate-sorted BAM file.
	print $FILEHANDLE "REMOVE_DUPLICATES=false ";
	print $FILEHANDLE "VALIDATION_STRINGENCY=STRICT ";
	print $FILEHANDLE "INPUT=".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$infileEnding.".bam ";  #InFile
	print $FILEHANDLE "OUTPUT=".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.".bam ";  #OutFile
	print $FILEHANDLE "METRICS_FILE=".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding."metric ", "\n\n";  #Metric file 
	
	## Copies file from temporary directory.
	print $FILEHANDLE "## Copy file from temporary directory\n"; 
	&MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.".b*", $outSampleDirectory, $FILEHANDLE);
	&MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding."metric", $outSampleDirectory."/".$infile.$outfileEnding."metric", $FILEHANDLE);
	print $FILEHANDLE "wait", "\n\n";

	if ( (${$scriptParameterHashRef}{'pPicardToolsMarkduplicates'} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	    ## Collect QC metadata info for later use
	    &SampleInfoQC(\%sampleInfo, 
			  {'familyID' => ${$scriptParameterHashRef}{'familyID'},
			   'sampleID' => $sampleID,
			   'programName' => "MarkDuplicates",
			   'infile' => $infile,
			   'outDirectory' => $outSampleDirectory,
			   'outFileEnding' => $outfileEnding."metric",
			   'outDataType' => "infileDependent"
			  });
	    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'} = $outSampleDirectory."/".$infile.$outfileEnding.".bam";
	}
    }
    else {  #No merged files

	## Set the number of cores to allocate per sbatch job.
	my $nrCores = &NrofCoresPerSbatch(\%{$scriptParameterHashRef}, scalar( @{${$laneHashRef}{$sampleID}} ));  #Detect the number of cores to use from lanes
	
	## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
	($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					    {'directoryID' => $sampleID,
					     'programName' => $programName,
					     'programDirectory' => lc($aligner),
					     'nrofCores' => $nrCores,
					     'processTime' => $time,
					     'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
					    });

	## Copies files from source to temporary folder. Loop over files specified by $arrayRef and collects files from $extractArrayRef.
	&MigrateFilesToTemp(\%{$scriptParameterHashRef}, \@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, \@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, $FILEHANDLE,
			    {'inSampleDirectory' => $inSampleDirectory,
			     'nrCores' => $nrCores,
			     'fileEnding' => $infileEnding.".b*"
			    });

	## PicardToolsMarkDuplicates
	print $FILEHANDLE "## Marking Duplicates\n";
	for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} });$infileCounter++) {  #For all files from independent of merged or not
	    
	    &PrintWait(\$infileCounter, \$nrCores, \$coreCounter, $FILEHANDLE);
	    
	    my $infile = ${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter];
	    
	    ## Writes java core commands to filehandle.
	    &JavaCore($FILEHANDLE,
		      {'memoryAllocation' => "Xmx4g",
		       'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
		       'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		      });

	    print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'picardToolsPath'}."/MarkDuplicates.jar ";
	    print $FILEHANDLE "TMP_DIR=".${$scriptParameterHashRef}{'tempDirectory'}." ";  #Temp Directory
	    print $FILEHANDLE "ASSUME_SORTED=true ";
	    print $FILEHANDLE "CREATE_INDEX=TRUE ";  #Create a BAM index when writing a coordinate-sorted BAM file.
	    print $FILEHANDLE "REMOVE_DUPLICATES=false ";
	    print $FILEHANDLE "VALIDATION_STRINGENCY=STRICT ";
	    print $FILEHANDLE "INPUT=".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$infileEnding.".bam ";  #InFile
	    print $FILEHANDLE "OUTPUT=".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.".bam ";  #OutFile
	    print $FILEHANDLE "METRICS_FILE=".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding."metric &","\n\n";  #Metric file  

	    if ( (${$scriptParameterHashRef}{'pPicardToolsMarkduplicates'} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

		## Collect QC metadata info for later use
		&SampleInfoQC(\%sampleInfo, 
			       {'familyID' => ${$scriptParameterHashRef}{'familyID'},
				'sampleID' => $sampleID,
				'programName' => "MarkDuplicates",
				'infile' => $infile,
				'outDirectory' => $outSampleDirectory,
				'outFileEnding' => $outfileEnding."metric",
				'outDataType' => "infileDependent"
			       });
		${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'} = $outSampleDirectory."/".$infile.$outfileEnding.".bam";
	    }
	}
	print $FILEHANDLE "wait", "\n\n";

	## Copies files from temporary folder to source. Loop over files specified by $arrayRef and collects files from $extractArrayRef.
	&MigrateFilesFromTemp(\@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, \@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, $outSampleDirectory, ${$scriptParameterHashRef}{'tempDirectory'}, $nrCores, $outfileEnding.".b*", $FILEHANDLE);
	&MigrateFilesFromTemp(\@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, \@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, $outSampleDirectory, ${$scriptParameterHashRef}{'tempDirectory'}, $nrCores, $outfileEnding."metric", $FILEHANDLE);
    }

    &RemoveDirectory(\${$scriptParameterHashRef}{'tempDirectory'}, $FILEHANDLE);

    close($FILEHANDLE);

    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'sampleID' => $sampleID,
			'dependencies' => 1, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}


sub PicardToolsMerge { 
#Merges all bam files using PicardTools MergeSamFiles within each sampleid and files generated previously (option if provided with '-picardToolsMergeSamFilesPrevious'). The merged files have to be sorted before attempting to merge.
##PicardToolsMarkDuplicates
    
##Function : Mark duplicated reads using PicardTools MarkDuplicates in files generated from alignment (sorted, merged).
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $infilesLaneNoEndingHashRef, $laneHashRef, $sampleID, $aligner, $fileEnding, $programName
##         : $parameterHashRef           => The parameter hash {REF}
##         : $scriptParameterHashRef     => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef          => Info on samples and family hash {REF}
##         : $infilesLaneNoEndingHashRef => The infile(s) without the ".ending" {REF}
##         : $laneHashRef                => The lane info hash {REF}
##         : $sampleID                   => The sampleID
##         : $aligner                    => The aligner used in the analysis
##         : $fileEnding                 => The sampleID file ending to use 
##         : $programName                => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $infilesLaneNoEndingHashRef = $_[3];
    my $laneHashRef = $_[4];
    my $sampleID = $_[5];
    my $aligner = $_[6];
    my $fileEnding = $_[7];
    my $programName = $_[8];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $nrCores = 1;

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					   {'directoryID' => $sampleID,
					    'programName' => $programName,
					    'programDirectory' => lc($aligner),
					    'nrofCores' => $nrCores,
					    'processTime' => 20,
					    'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
					   });

    ## Assign directories
    my $inSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner;
    my $outSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner;
   
    my $infileEnding;

    if (${$scriptParameterHashRef}{'analysisType'} ne "rapid") {

	$infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pPicardToolsSortSam'}{'fileEnding'};
    }    
    else {  #Rapid mode used

	$infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pPicardToolsMergeRapidReads'}{'fileEnding'};
    }
    my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pPicardToolsMergeSamFiles'}{'fileEnding'};
    my $lanes = join("",@{${$laneHashRef}{$sampleID}});  #Extract lanes

    ## Check that we have something to merge and then merge current files before merging with previously merged files
    if (scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} }) > 1) {

	## Copies files from source to temporary folder. Loop over files specified by $arrayRef and collects files from $extractArrayRef.
	&MigrateFilesToTemp(\%{$scriptParameterHashRef}, \@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, \@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, $FILEHANDLE,
			    {'inSampleDirectory' => $inSampleDirectory,
			     'nrCores' => $nrCores,
			     'fileEnding' => $infileEnding.".b*"
			    });

	## PicardToolsMergeSamFiles
	print $FILEHANDLE "## Merging alignment files\n";
	for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} });$infileCounter++) {  #For all files from 

	    my $infile = ${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter];

	    if ($infileCounter eq 0) {  #First round of loop
		
		## Writes java core commands to filehandle.
		&JavaCore($FILEHANDLE,
			  {'memoryAllocation' => "Xmx4g",
			   'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
			   'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
			  });

		print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'picardToolsPath'}."/MergeSamFiles.jar ";
		print $FILEHANDLE "TMP_DIR=".${$scriptParameterHashRef}{'tempDirectory'}." ";  #Temp Directory
		print $FILEHANDLE "CREATE_INDEX=TRUE ";  #create a BAM index when writing a coordinate-sorted BAM file.
		print $FILEHANDLE "OUTPUT=".${$scriptParameterHashRef}{'tempDirectory'}."/".$sampleID."_lanes_".$lanes.$outfileEnding.".bam ";  #OutFile
	    }
	    print $FILEHANDLE "INPUT=".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$infileEnding.".bam ";  #InFile
	}

	if ( (${$scriptParameterHashRef}{'pPicardToolsMergeSamFiles'} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {
	    
	    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'} = $outSampleDirectory."/".$sampleID."_lanes_".$lanes.$outfileEnding.".bam";
	}
    }
    ## Merge previously merged files with merged files generated this run
    if ( (${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'picardToolsMergeSamFilesPrevious'} == 1) && (scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} }) > 1) ) {

	for (my $mergeFileCounter=0;$mergeFileCounter<scalar(@{${$scriptParameterHashRef}{'picardToolsMergeSamFilesPrevious'}});$mergeFileCounter++) {
	    
	    if (${$scriptParameterHashRef}{'picardToolsMergeSamFilesPrevious'}[$mergeFileCounter] =~ /$sampleID/) {  #Look for sampleID in previously generated file to be merged with current run to be able to merge correct files within sampleID

		## PicardToolsMergeSamFiles
		print $FILEHANDLE "## Merging alignment files\n";

		## Copy file(s) to temporary directory
		print $FILEHANDLE "## Copy file(s) to temporary directory\n";

		my $picardToolsMergeSamFilesPreviousFile = &MigrateFileToTemp($FILEHANDLE, 
									       {'path' => ${$scriptParameterHashRef}{'picardToolsMergeSamFilesPrevious'}[$mergeFileCounter],
										'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
									       });
		print $FILEHANDLE "wait", "\n\n";
		
		if (${$scriptParameterHashRef}{'picardToolsMergeSamFilesPrevious'}[$mergeFileCounter] =~ /lane(\d+)|s_(\d+)/) {  #Look for lanes_ or lane\d in previously generated file to be merged with current run to be able to extract previous lanes
		    
		    my $mergeLanes; if($1) {$mergeLanes = $1;} else {$mergeLanes = $2;}  #Make sure to always supply lanes from previous regexp		    

		    ## Writes java core commands to filehandle.
		    &JavaCore($FILEHANDLE,
			      {'memoryAllocation' => "Xmx4g",
			       'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
			       'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
			      });

		    print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'picardToolsPath'}."/MergeSamFiles.jar ";
		    print $FILEHANDLE "TMP_DIR=".${$scriptParameterHashRef}{'tempDirectory'}." ";  #Temp directory
		    print $FILEHANDLE "CREATE_INDEX=TRUE ";  #Create a BAM index when writing a coordinate-sorted BAM file.
		    print $FILEHANDLE "OUTPUT=".${$scriptParameterHashRef}{'tempDirectory'}."/".$sampleID."_lanes_".$mergeLanes.$lanes.$outfileEnding.".bam ";  #OutFile
		    print $FILEHANDLE "INPUT=".${$scriptParameterHashRef}{'tempDirectory'}."/".$sampleID."_lanes_".$lanes.$outfileEnding.".bam ";  #InFile from previous merge
		    print $FILEHANDLE "INPUT=".${$scriptParameterHashRef}{'tempDirectory'}."/".$picardToolsMergeSamFilesPreviousFile, "\n\n";  #$mergeLanes contains lane info on previous merge, $infilesLaneNoEnding{$sampleID}[0] uses @RG for very first .bam file to include read group for subsequent merges. Complete path. 
		
		    ## Copies file from temporary directory.
		    print $FILEHANDLE "## Copy file from temporary directory\n"; 
		    &MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$sampleID."_lanes_".$mergeLanes.$lanes.$outfileEnding.".b*", $outSampleDirectory, $FILEHANDLE);
		    print $FILEHANDLE "wait", "\n\n";

		    if ( (${$scriptParameterHashRef}{'pPicardToolsMergeSamFiles'} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

			${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'} = $outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes.$lanes.$outfileEnding.".bam";
		    }
		}
	    }
	}
    }
    ## Merge files previously merged to single file with single file generated this run
    elsif (${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'picardToolsMergeSamFilesPrevious'} == 1) {

	for (my $mergeFileCounter=0;$mergeFileCounter<scalar(@{${$scriptParameterHashRef}{'picardToolsMergeSamFilesPrevious'}});$mergeFileCounter++) {

	    ## Copy file(s) to temporary directory
	    print $FILEHANDLE "## Copy file(s) to temporary directory\n";
	    my $picardToolsMergeSamFilesPreviousFile = &MigrateFileToTemp($FILEHANDLE, 
									   {'path' => ${$scriptParameterHashRef}{'picardToolsMergeSamFilesPrevious'}[$mergeFileCounter],
									    'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
									   });
	    print $FILEHANDLE "wait", "\n\n";

	    if (${$scriptParameterHashRef}{'picardToolsMergeSamFilesPrevious'}[$mergeFileCounter] =~ /lane(\d+)|s_(\d+)/) {  #Look for lanes_ or lane\d in previously generated file to be merged with current run to be able to extract previous lanes
		
		my $mergeLanes; if($1) {$mergeLanes = $1;} else {$mergeLanes = $2;}  #Make sure to always supply lanes from previous regexp
		my $infile = ${$infilesLaneNoEndingHashRef}{$sampleID}[0];  #Can only be 1 element in array due to previous if statement		    
		
		## PicardToolsMergeSamFiles
		print $FILEHANDLE "## Merging alignment files\n";

		## Writes java core commands to filehandle.
		&JavaCore($FILEHANDLE,
			  {'memoryAllocation' => "Xmx4g",
			   'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
			   'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
			  });
		
		print $FILEHANDLE "jar ".${$scriptParameterHashRef}{'picardToolsPath'}."/MergeSamFiles.jar ";
		print $FILEHANDLE "TMP_DIR=".${$scriptParameterHashRef}{'tempDirectory'}." ";  #Temp Directory
		print $FILEHANDLE "CREATE_INDEX=TRUE ";  #create a BAM index when writing a coordinate-sorted BAM file.
		print $FILEHANDLE "OUTPUT=".${$scriptParameterHashRef}{'tempDirectory'}."/".$sampleID."_lanes_".$mergeLanes.$lanes.$outfileEnding.".bam ";  #OutFile
		print $FILEHANDLE "INPUT=".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$infileEnding.".bam ";  #InFile
		print $FILEHANDLE "INPUT=".${$scriptParameterHashRef}{'tempDirectory'}."/".$picardToolsMergeSamFilesPreviousFile,"\n\n";  #$mergeLanes contains lane info on previous merge, $infilesLaneNoEnding{$sampleID}[0] uses @RG for very first .bam file to include read group for subsequent merges. Complete path. 

		## Copies file from temporary directory.
		print $FILEHANDLE "## Copy file from temporary directory\n"; 
		&MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$sampleID."_lanes_".$mergeLanes.$lanes.$outfileEnding.".b*", $outSampleDirectory, $FILEHANDLE);
		print $FILEHANDLE "wait", "\n\n";

		if ( (${$scriptParameterHashRef}{'pPicardToolsMergeSamFiles'} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {
		    
		    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'} = $outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes.$lanes.$outfileEnding.".bam";
		}
	    }
	}
    }
    else {

	##Copies file from temporary directory (From first merge)
	print $FILEHANDLE "\n\n## Copy file from temporary directory\n"; 
	&MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$sampleID."_lanes_".$lanes.$outfileEnding.".b*", $outSampleDirectory, $FILEHANDLE);
	print $FILEHANDLE "wait", "\n\n";
    }

    &RemoveDirectory(\${$scriptParameterHashRef}{'tempDirectory'}, $FILEHANDLE);

    close($FILEHANDLE);

    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'sampleID' => $sampleID,
			'dependencies' => 1, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}


sub PicardToolsSortSamIndex { 

##PicardToolsSortSamIndex
    
##Function : Sort and indexes bam files using PicardTools sort and index.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $infileHashRef, $inDirPathHashRef, $infilesLaneNoEndingHashRef, $sampleID, $aligner, $programName
##         : $parameterHashRef           => The parameter hash {REF}
##         : $scriptParameterHashRef     => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef          => Info on samples and family hash {REF}
##         : $infileHashRef              => The infiles hash {REF}
##         : $inDirPathHashRef           => The indirectories path(s) hash {REF}
##         : $infilesLaneNoEndingHashRef => The infile(s) without the ".ending" {REF}
##         : $sampleID                   => The sampleID
##         : $aligner                    => The aligner used in the analysis
##         : $programName                => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $infileHashRef = $_[3];
    my $inDirPathHashRef = $_[4];
    my $infilesLaneNoEndingHashRef = $_[5];
    my $sampleID = $_[6];
    my $aligner = $_[7];
    my $programName = $_[8];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $sbatchScriptTracker=0;
    my $time=0;
    my $infileSize;

    for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} });$infileCounter++) {  #For all files

	if (${$scriptParameterHashRef}{'analysisType'} eq "genomes") {
	    
	    $time = 25;  
	}
	else {
	    
	    $time = 15;
	}

	## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
	my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					       {'directoryID' => $sampleID,
						'programName' => $programName,
						'programDirectory' => lc($aligner),
						'processTime' => $time,
						'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
					       });
    
	my $inSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner;
	my $outSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner;
	my $infile = ${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter];
	my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{"p".$programName}{'fileEnding'};

	## Copy file(s) to temporary directory
	print $FILEHANDLE "## Copy file(s) to temporary directory\n"; 
	&MigrateFileToTemp($FILEHANDLE, 
			    {'path' => $inSampleDirectory."/".$infile.".bam ",
			     'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
			    });
	print $FILEHANDLE "wait", "\n\n";

	## PicardTools Sort and Index
	print $FILEHANDLE "#Sorting the reads\n\n";

	## Writes java core commands to filehandle.
	&JavaCore($FILEHANDLE,
		  {'memoryAllocation' => "Xmx4g",
		   'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
		   'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
		  });
	
	print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'picardToolsPath'}."/SortSam.jar ";
	print $FILEHANDLE "TMP_DIR=".${$scriptParameterHashRef}{'tempDirectory'}." ";  #Temp Directory
	print $FILEHANDLE "SORT_ORDER=coordinate"." ";  #Sort per contig and coordinate
	print $FILEHANDLE "CREATE_INDEX=TRUE ";  #create a BAM index when writing a coordinate-sorted BAM file. 
	print $FILEHANDLE "INPUT=".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.".bam ";  #InFile
	print $FILEHANDLE "OUTPUT=".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.".bam", "\n\n";  #Outfile	

	## Copies file from temporary directory.
	print $FILEHANDLE "## Copy file from temporary directory\n"; 
	&MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.$outfileEnding.".b*", $outSampleDirectory, $FILEHANDLE);
	print $FILEHANDLE "wait", "\n\n";	

	&RemoveDirectory(\${$scriptParameterHashRef}{'tempDirectory'}, $FILEHANDLE);	

	close($FILEHANDLE);

	if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'} = $outSampleDirectory."/".$infile.$outfileEnding.".bam";

	    &FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
			   {'sampleID' => $sampleID,
			    'dependencies' => 4, 
			    'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			    'sbatchFileName' => $fileName,
			    'sbatchScriptTracker' => $sbatchScriptTracker
			   });
	}
	$sbatchScriptTracker++; 
    }
}

sub BWA_Sampe {

##BWA_Sampe
    
##Function : Perform alignment of BWA Aln index reads using BWA sampe.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $infileHashRef, $inDirPathHashRef, $infilesLaneNoEndingHashRef, $infilesBothStrandsNoEndingHashRef, $sampleID, $aligner, $programName
##         : $parameterHashRef                  => The parameter hash {REF}
##         : $scriptParameterHashRef            => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef                 => Info on samples and family hash {REF}
##         : $infileHashRef                     => The infiles hash {REF}
##         : $inDirPathHashRef                  => The indirectories path(s) hash {REF}
##         : $infilesLaneNoEndingHashRef        => The infile(s) without the ".ending" {REF}
##         : $infilesBothStrandsNoEndingHashRef => The infile(s) without the ".ending" and strand info {REF}
##         : $sampleID                          => The sampleID
##         : $aligner                           => The aligner used in the analysis
##         : $programName                       => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $infileHashRef = $_[3];
    my $inDirPathHashRef = $_[4];
    my $infilesLaneNoEndingHashRef = $_[5];
    my $infilesBothStrandsNoEndingHashRef = $_[6];
    my $sampleID = $_[7];
    my $aligner = $_[8];
    my $programName = $_[9];
    
    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $time=0;
    my $infileSize;
    my $pairedEndTracker = 0;

    for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} });$infileCounter++) {  #For all files from BWA aln but process in the same command i.e. both reads per align call
	 
	if (${$scriptParameterHashRef}{'analysisType'} eq "genomes") {
	    
	    $time = 40;  
	}
	else {
	    
	    $time = 20;
	}

	my $nrCores = 2;
	my $sequenceRunMode = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{${$infilesLaneNoEndingHashRef}{ $sampleID }[$infileCounter]}{'sequenceRunType'};  #Collect paired-end or single-end sequence run mode

	## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
	my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					       {'directoryID' => $sampleID,
						'programName' => $programName,
						'programDirectory' => lc($aligner),
						'nrofCores' => $nrCores,
						'processTime' => $time,
						'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
					       });
    
	## Assign directories
	my $FASTQinSampleDirectory = ${$inDirPathHashRef}{$sampleID};
	my $inSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/bwa";
	my $outSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/bwa";

	my $infile = $infile{$sampleID}[$pairedEndTracker]; #For required .fastq file

	## Copies files from source to temporary folder. Loop over files specified by $arrayRef and collects files from $extractArrayRef.
	&MigrateFilesToTemp(\%{$scriptParameterHashRef}, \@{ ${$infileHashRef}{$sampleID} }, \@{ ${$infileHashRef}{$sampleID} }, $FILEHANDLE,
			    {'inSampleDirectory' => $FASTQinSampleDirectory,
			     'nrCores' => $nrCores
			    });  #Fastq files
	&MigrateFilesToTemp(\%{$scriptParameterHashRef}, \@{ ${$infilesBothStrandsNoEndingHashRef}{$sampleID} }, \@{ ${$infilesBothStrandsNoEndingHashRef}{$sampleID} }, $FILEHANDLE,
			    {'inSampleDirectory' => $inSampleDirectory,
			     'nrCores' => $nrCores,
			     'fileEnding' => ".sai*"
			    });
	
	## BWA Sampe	
	print $FILEHANDLE "## Aligning reads\n";
	print $FILEHANDLE "bwa sampe ";
	print $FILEHANDLE "-r ".'"@RG\tID:'.${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter].'\tSM:'.$sampleID.'\tPL:ILLUMINA" '.${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #read group header line
	print $FILEHANDLE ${$scriptParameterHashRef}{'tempDirectory'}."/".${$infilesBothStrandsNoEndingHashRef}{$sampleID}[$pairedEndTracker].".sai ";  #Read 1

	if ( $sequenceRunMode eq "Paired-end") {

	    $pairedEndTracker = $pairedEndTracker+1;  #Increment to collect correct read 2 from %infile
	    print $FILEHANDLE ${$scriptParameterHashRef}{'tempDirectory'}."/".${$infilesBothStrandsNoEndingHashRef}{$sampleID}[$pairedEndTracker].".sai ";  #Read 2
	}

	print $FILEHANDLE ${$scriptParameterHashRef}{'tempDirectory'}."/".$infile." ";  #Fastq read 1
	
	if ( $sequenceRunMode eq "Paired-end") { 

	    print $FILEHANDLE ${$scriptParameterHashRef}{'tempDirectory'}."/".${$infileHashRef}{$sampleID}[$pairedEndTracker]." ";  #Fastq read 2
	}

	print $FILEHANDLE "> ".${$scriptParameterHashRef}{'tempDirectory'}."/".${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter].".sam", "\n\n";  #Outfile (SAM)

	## Convert SAM to BAM using samTools view
	print $FILEHANDLE "## Convert SAM to BAM\n";
	print $FILEHANDLE "samtools view -bS ".${$scriptParameterHashRef}{'tempDirectory'}."/".${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter].".sam ";  #Infile (SAM)
	print $FILEHANDLE "> ".${$scriptParameterHashRef}{'tempDirectory'}."/".${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter].".bam", "\n\n";  #Outfile (BAM)

	## Copies file from temporary directory.
	print $FILEHANDLE "## Copy file from temporary directory\n"; 
	&MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter].".bam", $outSampleDirectory, $FILEHANDLE);
	print $FILEHANDLE "wait", "\n\n";

	&RemoveDirectory(\${$scriptParameterHashRef}{'tempDirectory'}, $FILEHANDLE);

	close($FILEHANDLE);

	if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'} = $outSampleDirectory."/".${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter].".bam";

	    &FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
			   {'sampleID' => $sampleID,
			    'dependencies' => 3, 
			    'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			    'sbatchFileName' => $fileName,
			    'sbatchScriptTracker' => $infileCounter
			   });
	}
	$pairedEndTracker++;
    }
}


sub BWA_Aln {

##BWA_Aln
    
##Function : Generates BWA aln index on fastq files.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $infileHashRef, $inDirPathHashRef, $infilesLaneNoEndingHashRef, $infilesBothStrandsNoEndingHashRef, $sampleID, $aligner, $programName
##         : $parameterHashRef                  => The parameter hash {REF}
##         : $scriptParameterHashRef            => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef                 => Info on samples and family hash {REF}
##         : $infileHashRef                     => The infiles hash {REF}
##         : $inDirPathHashRef                  => The indirectories path(s) hash {REF}
##         : $infilesLaneNoEndingHashRef        => The infile(s) without the ".ending" {REF}
##         : $infilesBothStrandsNoEndingHashRef => The infile(s) without the ".ending" and strand info {REF}
##         : $sampleID                          => The sampleID
##         : $aligner                           => The aligner used in the analysis
##         : $programName                       => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $infileHashRef = $_[3];
    my $inDirPathHashRef = $_[4];
    my $infilesLaneNoEndingHashRef = $_[5];
    my $infilesBothStrandsNoEndingHashRef = $_[6];
    my $sampleID = $_[7];
    my $aligner = $_[8];
    my $programName = $_[9];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $time = ceil(2.5*scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} }));  #One full lane on Hiseq takes approx. 2,5 h for BWA_Aln to process, round up to nearest full hour.
    my $nrCores = 0;

    for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} });$infileCounter++) {  #For all files   

	## Adjust the number of cores to be used in the analysis according to sequencing mode requirements.
	&AdjustNrCoresToSeqMode(\$nrCores, \${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{${$infilesLaneNoEndingHashRef}{ $sampleID }[$infileCounter]}{'sequenceRunType'});
    }

    ## Set the number of cores to allocate per sbatch job.
    $nrCores = &NrofCoresPerSbatch(\%{$scriptParameterHashRef}, $nrCores );  #Make sure that the number of cores does not exceed maximum after incrementing above

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					   {'directoryID' => $sampleID,
					    'programName' => $programName,
					    'programDirectory' => lc($aligner),
					    'nrofCores' => $nrCores,
					    'processTime' => $time,
					    'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
					   });

    ## Assign directories
    my $inSampleDirectory =  ${$inDirPathHashRef}{$sampleID};
    my $outSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/bwa";

    my $coreCounter=1;    

    ## Copies files from source to temporary folder. Loop over files specified by $arrayRef and collects files from $extractArrayRef
    &MigrateFilesToTemp(\%{$scriptParameterHashRef}, \@{ ${$infileHashRef}{$sampleID} }, \@{ ${$infileHashRef}{$sampleID} }, $FILEHANDLE,
			 {'inSampleDirectory' => $inSampleDirectory,
			  'nrCores' => $nrCores,});

    ## BWA Aln
    print $FILEHANDLE "## Creating .sai index\n";
    for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infileHashRef}{$sampleID} });$infileCounter++) {

	&PrintWait(\$infileCounter, \$nrCores, \$coreCounter, $FILEHANDLE);

	my $infile = ${$infileHashRef}{$sampleID}[$infileCounter];

	print $FILEHANDLE "bwa aln ";
	print $FILEHANDLE "-k 1 ";  #maximum differences in the seed
	print $FILEHANDLE "-t 4 ";  #number of threads
	print $FILEHANDLE "-n 3 ";  #max diff (int) or missing prob under 0.02 err rate (float)
	print $FILEHANDLE "-q ".${$scriptParameterHashRef}{'bwaAlnQualityTrimming'}." ";  #Quality trimming
	print $FILEHANDLE ${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference
	print $FILEHANDLE ${$scriptParameterHashRef}{'tempDirectory'}."/".$infile." ";  #InFile
	print $FILEHANDLE "> ".${$scriptParameterHashRef}{'tempDirectory'}."/".${$infilesBothStrandsNoEndingHashRef}{$sampleID}[$infileCounter].".sai &", "\n\n";  #OutFile 
    }
    print $FILEHANDLE "wait", "\n\n";

    ## Copies files from temporary folder to source. Loop over files specified by $arrayRef and collects files from $extractArrayRef.
    &MigrateFilesFromTemp(\@{ ${$infilesBothStrandsNoEndingHashRef}{$sampleID} }, \@{ ${$infilesBothStrandsNoEndingHashRef}{$sampleID} }, $outSampleDirectory, ${$scriptParameterHashRef}{'tempDirectory'}, $nrCores, ".sai", $FILEHANDLE);
    
    &RemoveDirectory(\${$scriptParameterHashRef}{'tempDirectory'}, $FILEHANDLE);

    close($FILEHANDLE);

    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {   

	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
			   {'sampleID' => $sampleID,
			    'dependencies' => 1, 
			    'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			    'sbatchFileName' => $fileName
			   });
    }
}

sub PicardToolsMergeRapidReads { 

##PicardToolsMergeRapidReads
    
##Function : Merges all batch read processes to one file using PicardTools MergeSamFiles within each sampleid. The read batch proccessed files have to be sorted before attempting to merge.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $infilesLaneNoEndingHashRef, $sampleID, $aligner, $programName
##         : $parameterHashRef           => The parameter hash {REF}
##         : $scriptParameterHashRef     => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef          => Info on samples and family hash {REF}
##         : $infilesLaneNoEndingHashRef => The infile(s) without the ".ending" {REF}
##         : $sampleID                   => The sampleID
##         : $aligner                    => The aligner used in the analysis
##         : $programName                => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $infilesLaneNoEndingHashRef = $_[3];
    my $sampleID = $_[4];
    my $aligner = $_[5];
    my $programName = $_[6];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					   {'directoryID' => $sampleID,
					    'programName' => $programName,
					    'programDirectory' => lc($aligner),
					    'nrofCores' => ${$scriptParameterHashRef}{'maximumCores'},
					    'processTime' => 20,
					   });

    ## Assign directories
    my $inSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner;
    my $outSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner;

    my $infileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pBwaMem'}{'fileEnding'};
    my $outfileEnding = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{"p".$programName}{'fileEnding'};
    my $coreCounter=1;
    my $coreTracker=0;  #Required to portion out cores and files before wait and to track the MOS_BU outfiles to correct lane
    
    for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} });$infileCounter++) {  #For all files from 
	
	my $infile = ${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter];
	my $nrReadBatchProcesses = ${$sampleInfoHashRef}{${$scriptParameterHashRef}{'familyID'}}{$sampleID}{${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter]}{'pBwaMem'}{'ReadBatchProcesses'}; 

	if ($nrReadBatchProcesses > 0) {  #Check that we have read batch processes to merge

	    &PrintWait(\$coreTracker, \${$scriptParameterHashRef}{'maximumCores'}, \$coreCounter, $FILEHANDLE);

	    for (my $readBatchProcessesCount=0;$readBatchProcessesCount<$nrReadBatchProcesses;$readBatchProcessesCount++) {
		
		if ($readBatchProcessesCount eq 0) {
		    
		    print $FILEHANDLE "java -Xmx4g ";

		    &WriteUseLargePages($FILEHANDLE, \${$scriptParameterHashRef}{'javaUseLargePages'});
		    
		    print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'picardToolsPath'}."/MergeSamFiles.jar ";
		    print $FILEHANDLE "TMP_DIR=".${$scriptParameterHashRef}{'tempDirectory'}." ";  #Temp Directory
		    print $FILEHANDLE "OUTPUT=".$outSampleDirectory."/".${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter].$outfileEnding.".bam ";  #OutFile
		}
		print $FILEHANDLE "INPUT=".$inSampleDirectory."/".${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter]."_".$readBatchProcessesCount."_sorted.bam ";  #InFile(s)
	    }
	    print $FILEHANDLE "CREATE_INDEX=TRUE &";  #Create a BAM index when writing a coordinate-sorted BAM file.
	    print $FILEHANDLE "\n\n";
	    $coreTracker++;  #Track nr of merge calls for infiles so that wait can be printed at the correct intervals (dependent on ${$scriptParameterHashRef}{'maximumCores'})
	}
	else {  #Still needs to rename file to be included in potential merge of BAM files in next step
	    
	    print $FILEHANDLE "java -Xmx4g ";

	    &WriteUseLargePages($FILEHANDLE, \${$scriptParameterHashRef}{'javaUseLargePages'});
	    
	    print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'picardToolsPath'}."/MergeSamFiles.jar ";
	    print $FILEHANDLE "TMP_DIR=".${$scriptParameterHashRef}{'tempDirectory'}." ";  #Temp Directory
	    print $FILEHANDLE "OUTPUT=".$outSampleDirectory."/".${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter].$outfileEnding.".bam ";  #OutFile
	    
	    print $FILEHANDLE "INPUT=".$inSampleDirectory."/".${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter]."_0_sorted_rg.bam ";  #InFile
	    print $FILEHANDLE "CREATE_INDEX=TRUE &";  #Create a BAM index when writing a coordinate-sorted BAM file.
	    print $FILEHANDLE "\n\n";
	}
    }
    print $FILEHANDLE "wait", "\n\n";
    
    print $FILEHANDLE "#Remove Temp Directory\n\n";
    print $FILEHANDLE "rm ";
    print $FILEHANDLE "-rf ".${$scriptParameterHashRef}{'tempDirectory'}, "\n\n";  #Remove Temp Directory
    
    close($FILEHANDLE);

    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'sampleID' => $sampleID,
			'dependencies' => 1, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}


sub BWA_Mem {

##BWA_Mem
    
##Function : Performs alignment.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $infileHashRef, $inDirPathHashRef, $infilesLaneNoEndingHashRef, $infilesBothStrandsNoEndingHashRef, $sampleID, $aligner, $programName
##         : $parameterHashRef                  => The parameter hash {REF}
##         : $scriptParameterHashRef            => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef                 => Info on samples and family hash {REF}
##         : $infileHashRef                     => The infiles hash {REF}
##         : $inDirPathHashRef                  => The indirectories path(s) hash {REF}
##         : $infilesLaneNoEndingHashRef        => The infile(s) without the ".ending" {REF}
##         : $infilesBothStrandsNoEndingHashRef => The infile(s) without the ".ending" and strand info {REF}
##         : $sampleID                          => The sampleID
##         : $aligner                           => The aligner used in the analysis
##         : $programName                       => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $infileHashRef = $_[3];
    my $inDirPathHashRef = $_[4];
    my $infilesLaneNoEndingHashRef = $_[5];
    my $infilesBothStrandsNoEndingHashRef = $_[6];
    my $sampleID = $_[7];
    my $aligner = $_[8];
    my $programName = $_[9];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $fileName;
    my $infileSize;
    my $totalSbatchCounter = 0;
    my $pairedEndTracker = 0;

    ## Collect fastq file(s) size
    for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} });$infileCounter++) {  #For all infiles but process in the same command i.e. both reads per align call
	
	my $sequenceRunMode = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{${$infilesLaneNoEndingHashRef}{ $sampleID }[$infileCounter]}{'sequenceRunType'};  #Collect paired-end or single-end sequence run mode
	
	## Fastq.gz
	if (${$infileHashRef}{$sampleID}[$infileCounter] =~/.fastq.gz$/) {  #Files are already gz and presently the scalar for compression has not been investigated. Therefore no automatic time allocation can be performed.
	
	    if (${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{${$infilesLaneNoEndingHashRef}{ $sampleID }[$infileCounter]}{'sequenceRunType'} eq "Paired-end") {  #Second read direction if present
                $infileSize = -s ${$inDirPathHashRef}{$sampleID}."/".${$infileHashRef}{$sampleID}[$infileCounter+$infileCounter];
	    }
	    else {  #Single-end
                $infileSize = -s ${$inDirPathHashRef}{$sampleID}."/".${$infileHashRef}{$sampleID}[$infileCounter];
	    }
        }
        else {  #Files are in fastq format
	    
	    if (${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{${$infilesLaneNoEndingHashRef}{ $sampleID }[$infileCounter]}{'sequenceRunType'} eq "Paired-end") {  #Second read direction if present        
		$infileSize = -s ${$inDirPathHashRef}{$sampleID}."/".${$infileHashRef}{$sampleID}[$infileCounter+$infileCounter];  # collect .fastq file size to enable estimation of time required for aligning, +1 for syncing multiple infiles per sampleID. Hence, filesize will be calculated on read2 (should not matter).
	    }
	    else {  #Single-end
                $infileSize = -s ${$inDirPathHashRef}{$sampleID}."/".${$infileHashRef}{$sampleID}[$infileCounter];
	    }
        }
	## Parallelize alignment by spliting of alignmnet processes as the files are read
	if (${$scriptParameterHashRef}{'analysisType'} eq "rapid") {
	    
	    my $seqLength = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{${$infilesBothStrandsNoEndingHashRef}{ $sampleID }[$infileCounter]}{'sequenceLength'};
	    my ($numberNodes, $ReadNrofLines) = &DetermineNrofRapidNodes($seqLength, $infileSize);
	    
	    for (my $sbatchCounter=0;$sbatchCounter<$numberNodes-1;$sbatchCounter++) {  #Parallization for each file handled
		
		## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
		($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
						    {'directoryID' => $sampleID,
						     'programName' => $programName,
						     'programDirectory' => lc($aligner),
						     'nrofCores' => ${$scriptParameterHashRef}{'maximumCores'},
						     'processTime' => 5,
						    });
		
		my $readStart = $sbatchCounter *  $ReadNrofLines;  #Constant for gz files
		my $readStop = $readStart + ceil( $ReadNrofLines + 1);  #Constant for gz files	

		## Assign directories
		my $BWAinSampleDirectory = ${$inDirPathHashRef}{$sampleID};
		my $BWAoutSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner; 

		my $infile;

		if ($sequenceRunMode eq "Paired-end") {  #Second read direction if present
	
		    $infile = ${$infileHashRef}{$sampleID}[$infileCounter+$infileCounter];  #For required .fastq file
                }
                else {  #Single-end
		    
		    $infile = ${$infileHashRef}{$sampleID}[$infileCounter];  #For required .fastq file
                }
		
		## BWA Mem for each read batch	
		print $FILEHANDLE "bwa mem ";
		print $FILEHANDLE "-M ";  #Mark shorter split hits as secondary (for Picard compatibility). 
		print $FILEHANDLE "-t ".${$scriptParameterHashRef}{'maximumCores'}." ";  #Number of threads 
		print $FILEHANDLE "-R ".'"@RG\tID:'.${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter].'\tSM:'.$sampleID.'\tPL:ILLUMINA" '; #Read group header line
		print $FILEHANDLE ${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference

		print $FILEHANDLE "<( ";  #Pipe to BWA Mem (Read 1)
		print $FILEHANDLE "zcat ";  #Decompress Read 1
		print $FILEHANDLE $BWAinSampleDirectory."/".$infile." ";  #Read 1
		print $FILEHANDLE "| ";  #Pipe
		print $FILEHANDLE q?perl -ne 'if ( ($.>?.$readStart.q?) && ($.<?.$readStop.q?) ) {print $_;}' ?;  #Limit to sbatch script interval
		print $FILEHANDLE ") ";  #End Read 1

		if ($sequenceRunMode eq "Paired-end") {  #Second read direction if present
		      
		    print $FILEHANDLE "<( ";  #Pipe to BWA Mem (Read 2)
		    print $FILEHANDLE "zcat ";  #Decompress Read 2
		    print $FILEHANDLE $BWAinSampleDirectory."/".${$infileHashRef}{$sampleID}[$infileCounter+$infileCounter+1]." ";  #Read 2
		    print $FILEHANDLE "| ";  #Pipe
		    print $FILEHANDLE q?perl -ne 'if ( ($.>?.$readStart.q?) && ($.<?.$readStop.q?) ) {print $_;}' ?;  #Limit to sbatch script interval
		    print $FILEHANDLE ") ";  #End Read 2
		}

		print $FILEHANDLE "| ";  #Pipe SAM to BAM conversion of aligned reads
		print $FILEHANDLE "samtools view "; 
		print $FILEHANDLE "-S ";  #Input is SAM
		print $FILEHANDLE "-h ";  #Print header for the SAM output
		print $FILEHANDLE "-u ";  #Uncompressed BAM output
		print $FILEHANDLE "- ";  #/dev/stdin
		print $FILEHANDLE "| ";  #Pipe
		print $FILEHANDLE "intersectBed ";  #Limit output to only clinically interesting genes
		print $FILEHANDLE "-abam stdin ";  #The A input file is in BAM format.  Output will be BAM as well.
		print $FILEHANDLE "-b ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'bwaMemRapidDb'}." ";  #Db file of clinically relevant variants
		print $FILEHANDLE "> ".$BWAoutSampleDirectory."/".${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter]."_".$sbatchCounter.".bam", "\n\n";  #Outfile (BAM)
		
		print $FILEHANDLE "samtools sort ";
		print $FILEHANDLE $BWAoutSampleDirectory."/".${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter]."_".$sbatchCounter.".bam ";  #Infile
		print $FILEHANDLE $BWAoutSampleDirectory."/".${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter]."_".$sbatchCounter."_sorted", "\n\n";  #OutFile

		print $FILEHANDLE "samtools index ";
		print $FILEHANDLE $BWAoutSampleDirectory."/".${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter]."_".$sbatchCounter."_sorted.bam", "\n\n";  #OutFile

		close($FILEHANDLE);

		if ( (${$scriptParameterHashRef}{'pBwaMem'} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

		    &FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
			       {'sampleID' => $sampleID,
				'dependencies' => 3, 
				'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
				'sbatchFileName' => $fileName,
				'sbatchScriptTracker' => $totalSbatchCounter
			       });
		}
		$totalSbatchCounter++;

                ## Save sbatch Counter to track how many read batch processes we have engaged
		${$sampleInfoHashRef}{${$scriptParameterHashRef}{'familyID'}}{$sampleID}{${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter]}{'pBwaMem'}{'ReadBatchProcesses'} = $sbatchCounter+1;#Used to be  $sbatchCounter
		${$sampleInfoHashRef}{${$scriptParameterHashRef}{'familyID'}}{$sampleID}{'pBwaMem'}{'sbatchBatchProcesses'} = $totalSbatchCounter;
	    }
	}
	else {  #Not rapid mode align whole file

	    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
	    ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
						{'directoryID' => $sampleID,
						 'programName' => $programName,
						 'programDirectory' => lc($aligner),
						 'nrofCores' => ${$scriptParameterHashRef}{'maximumCores'},
						 'processTime' => 5,
						 'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
						 });
	
	    ## Assign directories
	    my $inSampleDirectory = ${$inDirPathHashRef}{$sampleID};
	    my $outSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner; 
	    
	    ## Copies file to temporary directory.
	    print $FILEHANDLE "## Copy file(s) to temporary directory\n"; 
	    &MigrateFileToTemp($FILEHANDLE, 
				{'path' => $inSampleDirectory."/".${$infileHashRef}{$sampleID}[$pairedEndTracker],
				 'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
				});  #Read 1
	    if ($sequenceRunMode eq "Paired-end") {  #Second read direction if present
		
		&MigrateFileToTemp($FILEHANDLE, 
				    {'path' => $inSampleDirectory."/".${$infileHashRef}{$sampleID}[$pairedEndTracker+1],
				     'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
				    });  #Read 2
	    }
	    print $FILEHANDLE "wait", "\n\n";

	    ## BWA MEM
	    print $FILEHANDLE "## Aligning reads and converting to BAM via samtools\n";
	    print $FILEHANDLE "bwa mem ";
	    print $FILEHANDLE "-M ";  #Mark shorter split hits as secondary (for Picard compatibility). 
	    print $FILEHANDLE "-t ".${$scriptParameterHashRef}{'maximumCores'}." ";  #Number of threads 
	    print $FILEHANDLE "-R ".'"@RG\tID:'.${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter].'\tSM:'.$sampleID.'\tPL:ILLUMINA" '; #Read group header line
	    print $FILEHANDLE ${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference
	    print $FILEHANDLE ${$scriptParameterHashRef}{'tempDirectory'}."/".${$infileHashRef}{$sampleID}[$pairedEndTracker]." ";  #Read 1

	    if ($sequenceRunMode eq "Paired-end") {  #Second read direction if present

		$pairedEndTracker = $pairedEndTracker+1;  #Increment to collect correct read 2 from %infile 
		print $FILEHANDLE ${$scriptParameterHashRef}{'tempDirectory'}."/".${$infileHashRef}{$sampleID}[$pairedEndTracker]." ";  #Read 2
	    }
	    $pairedEndTracker++;
	    print $FILEHANDLE "| ";  #Pipe SAM to BAM conversion of aligned reads
	    print $FILEHANDLE "samtools view "; 
	    print $FILEHANDLE "-S ";  #Input is SAM
	    print $FILEHANDLE "-h ";  #Print header for the SAM output
	    print $FILEHANDLE "-u ";  #Uncompressed BAM output
	    print $FILEHANDLE "- ";  #/dev/stdin
	    print $FILEHANDLE "> ".${$scriptParameterHashRef}{'tempDirectory'}."/".${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter].".bam", "\n\n";  #Outfile (BAM)

	    ## Copies file from temporary directory.
	    print $FILEHANDLE "## Copy file from temporary directory\n"; 
	    &MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter].".bam", $outSampleDirectory, $FILEHANDLE);
	    print $FILEHANDLE "wait", "\n\n";
	    
	    &RemoveDirectory(\${$scriptParameterHashRef}{'tempDirectory'}, $FILEHANDLE);

	    close($FILEHANDLE);

	    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

		${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'} = $outSampleDirectory."/".${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter].".bam";

		&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
			       {'sampleID' => $sampleID,
				'dependencies' => 3, 
				'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
				'sbatchFileName' => $fileName,
				'sbatchScriptTracker' => $infileCounter
			       });
	    }
	}
    }
}


sub MosaikAlign {

##MosaikAlign
    
##Function : Performs alignment.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $infileHashRef, $inDirPathHashRef, $infilesLaneNoEndingHashRef, $infilesBothStrandsNoEndingHashRef, $sampleID, $aligner, $programName
##         : $parameterHashRef                  => The parameter hash {REF}
##         : $scriptParameterHashRef            => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef                 => Info on samples and family hash {REF}
##         : $infileHashRef                     => The infiles hash {REF}
##         : $inDirPathHashRef                  => The indirectories path(s) hash {REF}
##         : $infilesLaneNoEndingHashRef        => The infile(s) without the ".ending" {REF}
##         : $infilesBothStrandsNoEndingHashRef => The infile(s) without the ".ending" and strand info {REF}
##         : $sampleID                          => The sampleID
##         : $aligner                           => The aligner used in the analysis
##         : $programName                       => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $infileHashRef = $_[3];
    my $inDirPathHashRef = $_[4];
    my $infilesLaneNoEndingHashRef = $_[5];
    my $infilesBothStrandsNoEndingHashRef = $_[6];
    my $sampleID = $_[7];
    my $aligner = $_[8];
    my $programName = $_[9];
    
    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $sbatchScriptTracker=0;
    my $time=0;
    my $infileSize;

    for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} });$infileCounter++) {  #For all infiles per lane
	   
	if (${$scriptParameterHashRef}{'analysisType'} eq "genomes") {
	    
	    $time = 80;  
	}
	else {
	    
	    $time = 40;
	}

	## Set parameters depending on sequence length
	my $seqLength = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{${$infilesBothStrandsNoEndingHashRef}{ $sampleID }[$infileCounter]}{'sequenceLength'};
	my $actParameter = 35;  #The alignment candidate threshold (length)
	my $bwParameter = 35;  #Specifies the Smith-Waterman bandwidth.

	if ($seqLength <= 36) {
	    
	    $actParameter = 20;
	    $bwParameter = 13;   
	}
	if ($seqLength > 36 && $seqLength <= 51) {
	    
	    $actParameter = 25;
	    $bwParameter = 21;   
	}
	if ($seqLength > 51 && $seqLength <= 76) {
	    
	    $bwParameter = 29;   
	}

	## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
	my ($fileName, $stdoutPath) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
							    {'directoryID' => $sampleID,
							     'programName' => $programName,
							     'programDirectory' => lc($aligner),
							     'nrofCores' => ${$scriptParameterHashRef}{'maximumCores'},
							     'processTime' => $time,
							     'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
							    });
	my ($volume, $directories, $stdoutFile) = File::Spec->splitpath($stdoutPath);  #Split to enable submission to &SampleInfoQC later

	## Assign directories
	my $inSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner;
	my $outSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner;
	my $infile = ${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter];
	
	## Copies file to temporary directory.
	print $FILEHANDLE "## Copy file to node\n"; 
	&MigrateFileToTemp($FILEHANDLE, 
			    {'path' => $inSampleDirectory."/".$infile.".dat",
			     'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
			    });
	print $FILEHANDLE "wait", "\n\n";
	
	## MosaikAlign
	print $FILEHANDLE "## Create node temporary MOSAIK directory\n";
	print $FILEHANDLE "mkdir -p ".${$scriptParameterHashRef}{'tempDirectory'}."/"."mosaik_tmp", "\n";
	print $FILEHANDLE "export MOSAIK_TMP=".${$scriptParameterHashRef}{'tempDirectory'}."/"."mosaik_tmp", "\n\n";

	print $FILEHANDLE "## Generating .bam file from .dat files\n";
	print $FILEHANDLE "MosaikAligner ";
	print $FILEHANDLE "-in ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.".dat ";  #Infile
	print $FILEHANDLE "-out ".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile." ";  #OutFile (MosaikAligner appends .bam to infile name)
	print $FILEHANDLE "-ia ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'mosaikAlignReference'}." ";  #Mosaik Reference
	print $FILEHANDLE "-annse ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'mosaikAlignNeuralNetworkSeFile'}." ";  #NerualNetworkSE
	print $FILEHANDLE "-hs 15 ";  #Hash size
	print $FILEHANDLE "-mm 4 ";  #The # of mismatches allowed
	print $FILEHANDLE "-mhp 100 "; #The maximum of positions stored per seed

	if (${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{${$infilesLaneNoEndingHashRef}{ $sampleID }[$infileCounter]}{'sequenceRunType'} eq "Paired-end") {  #Second read direction if present

	    print $FILEHANDLE "-annpe ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'mosaikAlignNeuralNetworkPeFile'}." ";  #NerualNetworkPE
	    print $FILEHANDLE "-ls 100 "; #Enable local alignment search for PE reads
	}
	print $FILEHANDLE "-act ".$actParameter." ";  #The alignment candidate threshold (length)
	print $FILEHANDLE "-bw ".$bwParameter." ";  #Specifies the Smith-Waterman bandwidth.
	print $FILEHANDLE "-j ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'mosaikJumpDbStub'}." ";  #JumpDatabase
	print $FILEHANDLE "-p ".${$scriptParameterHashRef}{'maximumCores'}, "\n\n";  #Nr of cores

	## BAM to SAM conversion
	print $FILEHANDLE "## BAM to SAM conversion\n";
	print $FILEHANDLE "java -Xmx4g ";

	&WriteUseLargePages($FILEHANDLE, \${$scriptParameterHashRef}{'javaUseLargePages'});

	print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'picardToolsPath'}."/SamFormatConverter.jar ";  #Make sure that the BAM file BIN field is correct (Mosaik v.2.2.3 does according to Picard not set the bin field correctly)
	print $FILEHANDLE "VALIDATION_STRINGENCY=SILENT ";  #Disable errors print 
	print $FILEHANDLE "INPUT=".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.".bam ";  #InFile
	print $FILEHANDLE "OUTPUT=".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.".sam ", "\n\n";  #OutFile
	
	## SAM to BAM conversion
	print $FILEHANDLE "## SAM to BAM conversion\n";
	print $FILEHANDLE "java -Xmx4g ";

	&WriteUseLargePages($FILEHANDLE, \${$scriptParameterHashRef}{'javaUseLargePages'});

	print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'picardToolsPath'}."/SamFormatConverter.jar ";  #Make sure that the BAM file BIN field is correct (Mosaik v.2.2.3 does according to Picard not set the bin field correctly)
	print $FILEHANDLE "INPUT=".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.".sam ";  #InFile
	print $FILEHANDLE "OUTPUT=".${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.".bam ", "\n\n";  #OutFile
	
	## Copies file from temporary directory.
	print $FILEHANDLE "## Copy file from temporary directory\n"; 
	&MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.".bam", $outSampleDirectory, $FILEHANDLE);
	&MigrateFileFromTemp(${$scriptParameterHashRef}{'tempDirectory'}."/".$infile.".stat", $outSampleDirectory, $FILEHANDLE);
	print $FILEHANDLE "wait", "\n\n";
	
	&RemoveDirectory(\${$scriptParameterHashRef}{'tempDirectory'}, $FILEHANDLE);

	close($FILEHANDLE);
	
	if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {
	    
	    ## Collect QC metadata info for later use                     	
	    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'MostCompleteBAM'}{'Path'} = $outSampleDirectory."/".$infile.".bam";	
	    &SampleInfoQC(\%sampleInfo, 
			   {'familyID' => ${$scriptParameterHashRef}{'familyID'},
			    'sampleID' => $sampleID,
			    'programName' => "MosaikAligner",
			    'infile' => $infile,
			    'outDirectory' => $directories,
			    'outFileEnding' => $stdoutFile,
			    'outDataType' => "infoDirectory"
			   });
	    
	    &FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
			       {'sampleID' => $sampleID,
				'dependencies' => 3, 
				'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
				'sbatchFileName' => $fileName,
				'sbatchScriptTracker' => $sbatchScriptTracker
			       });
	}
	$sbatchScriptTracker++;  #Tracks nr of sbatch scripts
    }
}


sub MosaikBuild {

##MosaikBuild
    
##Function : Generates Mosaik hash format on reads using MosaikBuild
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $infileHashRef, $inDirPathHashRef, $infilesLaneNoEndingHashRef, $laneHashRef, $sampleID, $aligner, $programName
##         : $parameterHashRef           => The parameter hash {REF}
##         : $scriptParameterHashRef     => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef          => Info on samples and family hash {REF}
##         : $infileHashRef              => The infiles hash {REF}
##         : $inDirPathHashRef           => The indirectories path(s) hash {REF}
##         : $infilesLaneNoEndingHashRef => The infile(s) without the ".ending" {REF}
##         : $laneHashRef                => The lane info hash {REF}
##         : $sampleID                   => The sampleID
##         : $aligner                    => The sampleID
##         : $programName                => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $infileHashRef = $_[3];
    my $inDirPathHashRef = $_[4];
    my $infilesLaneNoEndingHashRef = $_[5];
    my $laneHashRef = $_[6];
    my $sampleID = $_[7];
    my $aligner = $_[8];
    my $programName = $_[9];
    
    my $FILEHANDLE = IO::Handle->new();  #Create anonymous filehandle
    my $time = ceil(2.5*scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} }));  #Scale with the number of files to process
    
    ## Set the number of cores to allocate per sbatch job.
    my $nrCores = &NrofCoresPerSbatch(\%{$scriptParameterHashRef}, scalar( @{${$laneHashRef}{$sampleID}} ));  #Detect the number of cores to use from lanes
    
    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					   {'directoryID' => $sampleID,
					    'programName' => $programName,
					    'programDirectory' => lc($aligner),
					    'nrofCores' => $nrCores,
					    'processTime' => $time,
					    'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
					   });
    
    ## Assign directories
    my $inSampleDirectory = ${$inDirPathHashRef}{$sampleID};
    my $outSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".$aligner;

    my $coreCounter=1;
    my  $pairedEndTracker = 0;
    my $stParameter = "ILLUMINA";  #Default

    ## Copies files from source to temporary folder. Loop over files specified by $arrayRef and collects files from $extractArrayRef.
    &MigrateFilesToTemp(\%{$scriptParameterHashRef}, \@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, \@{ ${$infileHashRef}{$sampleID} }, $FILEHANDLE,
			{'sampleInfoHashRef' => \%{$sampleInfoHashRef},
			 'inSampleDirectory' => $inSampleDirectory,
			 'nrCores' => $nrCores,
			 'sampleID' => $sampleID
			});

    ## MosaikBuild
    print $FILEHANDLE "## Generating .dat file from fastq files\n";    
    for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} });$infileCounter++) {  #For all files
	
	my $sequenceRunMode = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{ ${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter] }{'sequenceRunType'};  #Collect paired-end or single-end sequence run mode
	my $coreTracker=0;  #Required to portion out cores and files before wait and to track the outfiles to correct lane

	&PrintWait(\$infileCounter, \$nrCores, \$coreCounter, $FILEHANDLE);

	print $FILEHANDLE "MosaikBuild ";
	print $FILEHANDLE "-id ".${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter]." ";  #Read group ID for BAM Header
	print $FILEHANDLE "-sam ".$sampleID." ";  #Sample name for BAM Header
	print $FILEHANDLE "-st ".$stParameter." ";  #Sequencing technology for BAM Header
	print $FILEHANDLE "-mfl ".${$scriptParameterHashRef}{'mosaikBuildMedianFragLength'}." ";  #Median Fragment Length
	print $FILEHANDLE "-q ".${$scriptParameterHashRef}{'tempDirectory'}."/".${$infileHashRef}{$sampleID}[$pairedEndTracker]." ";  #Read 1
	
	if ( $sequenceRunMode eq "Paired-end") {
    
	    $pairedEndTracker = $pairedEndTracker+1;  #Increment to collect correct read 2 from %infile
	    print $FILEHANDLE "-q2 ".${$scriptParameterHashRef}{'tempDirectory'}."/".${$infileHashRef}{$sampleID}[$pairedEndTracker]." ";  #Read 2
	} 

	$pairedEndTracker++;  #Increment to correctly track both single-end runs and paired-end runs
	print $FILEHANDLE "-out ".${$scriptParameterHashRef}{'tempDirectory'}."/".${$infilesLaneNoEndingHashRef}{$sampleID}[$infileCounter].".dat &", "\n\n";  #OutFile
    }
    print $FILEHANDLE "wait", "\n\n";

    ## Copies files from temporary folder to source. Loop over files specified by $arrayRef and collects files from $extractArrayRef.
    &MigrateFilesFromTemp(\@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, \@{ ${$infilesLaneNoEndingHashRef}{$sampleID} }, $outSampleDirectory, ${$scriptParameterHashRef}{'tempDirectory'}, $nrCores, ".dat", $FILEHANDLE);
    
    &RemoveDirectory(\${$scriptParameterHashRef}{'tempDirectory'}, $FILEHANDLE);

    close($FILEHANDLE);

    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) { 

	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%{$infilesLaneNoEndingHashRef},
		       {'sampleID' => $sampleID,
			'dependencies' => 1, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}   


sub FastQC {

##FastQC
    
##Function : Raw sequence quality analysis using FASTQC.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $infileHashRef, $inDirPathHashRef, $infilesLaneNoEndingHashRef, $infilesBothStrandsNoEndingHashRef, $sampleID, $programName
##         : $parameterHashRef                  => The parameter hash {REF}
##         : $scriptParameterHashRef            => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef                 => Info on samples and family hash {REF}
##         : $infileHashRef                     => The infiles hash {REF}
##         : $inDirPathHashRef                  => The indirectories path(s) hash {REF}
##         : $infilesLaneNoEndingHashRef        => The infile(s) without the ".ending" {REF}
##         : $infilesBothStrandsNoEndingHashRef => The infile(s) without the ".ending" and strand info {REF}
##         : $sampleID                          => The sampleID
##         : $programName                       => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $infileHashRef = $_[3];
    my $inDirPathHashRef = $_[4];
    my $infilesLaneNoEndingHashRef = $_[5];
    my $infilesBothStrandsNoEndingHashRef = $_[6];
    my $sampleID = $_[7];
    my $programName = $_[8];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $time = ceil(0.5*scalar( @{ ${$infileHashRef}{$sampleID} }));  #One full lane on Hiseq takes approx. 0.5 h for FASTQC to process, round up to nearest full hour.

    my $nrCores = 0;

    for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} });$infileCounter++) {  #For all files   

	## Adjust the number of cores to be used in the analysis according to sequencing mode requirements.
	&AdjustNrCoresToSeqMode(\$nrCores, \${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{ ${$infilesLaneNoEndingHashRef}{ $sampleID }[$infileCounter] }{'sequenceRunType'});
    }

    ## Set the number of cores to allocate per sbatch job.
    $nrCores = &NrofCoresPerSbatch(\%{$scriptParameterHashRef}, $nrCores );  #Make sure that the number of cores does not exceed maximum after incrementing above
    
    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					   {'directoryID' => $sampleID,
					    'programName' => $programName,
					    'programDirectory' => lc($programName),
					    'nrofCores' => $nrCores,
					    'processTime' => $time,
					   });

    ## Assign directories
    my $inSampleDirectory = ${$inDirPathHashRef}{$sampleID};
    my $outSampleDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".$sampleID."/".lc($programName);

    my $coreCounter=1;
    
    for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infileHashRef}{$sampleID} });$infileCounter++) {

	&PrintWait(\$infileCounter, \$nrCores, \$coreCounter, $FILEHANDLE);

	my $infile = ${$infileHashRef}{$sampleID}[$infileCounter];

	print $FILEHANDLE "fastqc ";
	print $FILEHANDLE $inSampleDirectory."/".$infile." ";  #InFile
	print $FILEHANDLE "--extract ";  #the zipped output file will be uncompressed in the same directory after it has been created.
	print $FILEHANDLE "-o ".$outSampleDirectory. " &", "\n\n";  #OutFile

##Collect QC metadata info for active program for later use
	if ( (${$scriptParameterHashRef}{'pFastQC'} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	    &SampleInfoQC(\%sampleInfo,
			   {'familyID' => ${$scriptParameterHashRef}{'familyID'},
			    'sampleID' => $sampleID,
			    'programName' => "FastQC",
			    'infile' => $infile,
			    'outDirectory' => $outSampleDirectory."/".${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{${$infilesBothStrandsNoEndingHashRef}{ $sampleID }[$infileCounter]}{'originalFileNameNoEnding'}."_fastqc",
			    'outFileEnding' => "fastqc_data.txt",
			    'outDataType' => "static"
			   });
	}
    }
    print $FILEHANDLE "wait", "\n";    
    close($FILEHANDLE);
    
    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%{$infilesLaneNoEndingHashRef},
		       {'sampleID' => $sampleID,
			'dependencies' => 2, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}


sub GZipFastq { 
 
##GZipFastq
    
##Function : Automatically gzips fastq files.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $infileHashRef, $inDirPathHashRef, $infilesLaneNoEndingHashRef, $sampleID, $programName
##         : $parameterHashRef           => The parameter hash {REF}
##         : $scriptParameterHashRef     => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef          => Info on samples and family hash {REF}
##         : $infileHashRef              => The infiles hash {REF}
##         : $inDirPathHashRef           => The indirectories path(s) hash {REF}
##         : $infilesLaneNoEndingHashRef => The infile(s) without the ".ending" {REF}
##         : $sampleID                   => The sampleID
##         : $programName                => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $infileHashRef = $_[3];
    my $inDirPathHashRef = $_[4];
    my $infilesLaneNoEndingHashRef = $_[5];
    my $sampleID = $_[6];
    my $programName = $_[7];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $time = ceil(1.5*scalar( @{ ${$infileHashRef}{$sampleID} }));  #One full lane on Hiseq takes approx. 1.5 h for gzip to process, round up to nearest full hour.

    my $nrCores = 0;

    for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} });$infileCounter++) {  #For all files   
	
	## Adjust the number of cores to be used in the analysis according to sequencing mode requirements.
	&AdjustNrCoresToSeqMode(\$nrCores, \${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{ ${$infilesLaneNoEndingHashRef}{ $sampleID }[$infileCounter] }{'sequenceRunType'});
    }

    ## Set the number of cores to allocate per sbatch job.
    $nrCores = &NrofCoresPerSbatch(\%{$scriptParameterHashRef}, $nrCores );  #Make sure that the number of cores does not exceed maximum after incrementing above
    
    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					   {'directoryID' => $sampleID,
					    'programName' => $programName,
					    'programDirectory' => lc($programName),
					    'nrofCores' => $nrCores,
					    'processTime' => $time,
					   });

    ## Assign directories
    my $inSampleDirectory = ${$inDirPathHashRef}{$sampleID};   

    my $coreCounter=1;
    my $uncompressedFileCounter = 0;  #Used to print wait at the right times since infiles cannot be used (can be a mixture of .gz and .fast files)

    print $FILEHANDLE "cd ".${$inDirPathHashRef}{$sampleID}, "\n\n";
    
    for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infileHashRef}{$sampleID} });$infileCounter++) {

	if (${$infileHashRef}{$sampleID}[$infileCounter] =~/.fastq$/) {  #For files ending with .fastq required since there can be a mixture (also .fastq.gz) within the sample dir

	    if ($uncompressedFileCounter == $coreCounter*${$scriptParameterHashRef}{'maximumCores'}) {  #Using only $scriptParameter{'maximumCores'} cores
		
		print $FILEHANDLE "wait", "\n\n";
		$coreCounter=$coreCounter+1;
	    }

	    my $infile = ${$infileHashRef}{$sampleID}[$infileCounter];

	    print $FILEHANDLE "gzip ";
	    print $FILEHANDLE $inSampleDirectory."/".$infile," &", "\n\n";  #InFile
	    $uncompressedFileCounter++;
	    ${$infileHashRef}{$sampleID}[$infileCounter] =~ s/.fastq/.fastq.gz/g;  #Replace the .fastq ending with .fastq.gz since this will execute before fastQC screen and mosaikBuild, hence changing the original file name ending from ".fastq" to ".fastq.gz". 
	}
    }
    print $FILEHANDLE "wait", "\n\n";

    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'sampleID' => $sampleID,
			'dependencies' => 0, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}


sub BuildAnnovarPreRequisites {

##BuildAnnovarPreRequisites
    
##Function : Creates the AnnovarPreRequisites.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $annovarTableHashRef, $familyID, $aligner, $programName
##         : $parameterHashRef       => The parameter hash {REF}
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $annovarTableHashRef    => annovarTableHashRef {REF}
##         : $familyID               => Family ID
##         : $aligner                => The aligner used in the analysis
##         : $programName            => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $annovarTableHashRef = $_[2];
    my $familyID = $_[3];
    my $aligner = $_[4];
    my $programName = $_[5];
    
    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    ${$parameterHashRef}{'annovarBuildReference'}{'buildFile'} = 0;  #Ensure that this subrutine is only executed once
    my $annovarTemporaryDirectory = ${$scriptParameterHashRef}{'annovarPath'}."/humandb/Db_temporary";  #Temporary download directory
    
    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					   {'directoryID' => $familyID,
					    'programName' => $programName,
					    'programDirectory' => lc($aligner),
					    'processTime' => 3,
					    'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
					   });

    $logger->warn("Will try to create required Annovar database files before executing ".$programName."\n");

    print $FILEHANDLE "#Make temporary download directory\n\n"; 
    print $FILEHANDLE "mkdir -p ".$annovarTemporaryDirectory."; ", "\n\n"; 

    print $FILEHANDLE "#Downloading Annovar Db files", "\n\n";

    for (my $tableNamesCounter=0;$tableNamesCounter<scalar(@{${$scriptParameterHashRef}{'annovarTableNames'}});$tableNamesCounter++) {  #For all specified table names
	
	if (${$parameterHashRef}{${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter]}{'buildFile'} eq 1) {
	    
	    print $FILEHANDLE "perl ".${$scriptParameterHashRef}{'annovarPath'}."/annotate_variation.pl ";  #Annovar script 
	    print $FILEHANDLE "-buildver ".${$scriptParameterHashRef}{'annovarGenomeBuildVersion'}." ";  #GenomeBuild version
	    print $FILEHANDLE "-downdb ".${$annovarTableHashRef}{${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter]}{'download'}." ";  #Db to download
	    
	    if (defined(${$annovarTableHashRef}{${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter]}{'ucscAlias'})) {
		
		print $FILEHANDLE "-webfrom ucsc ";  #Download from ucsc
	    }
	    else {
		
		print $FILEHANDLE "-webfrom annovar ";  #Download from annovar
	    }
	    print $FILEHANDLE $annovarTemporaryDirectory."/ ", "\n\n";  #Annovar/humandb directory is assumed

	    if (${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter] =~/ensGene|refGene/) {  #Special case for MT download
		
		print $FILEHANDLE "perl ".${$scriptParameterHashRef}{'annovarPath'}."/annotate_variation.pl ";  #Annovar script 
		print $FILEHANDLE "-buildver GRCh37_MT ";  #GenomeBuild version
		print $FILEHANDLE "-downdb ensGene ";  #Db to download
		print $FILEHANDLE "-webfrom annovar ";  #Download from annovar
		print $FILEHANDLE $annovarTemporaryDirectory."/ ", "\n\n";  #annovar/humandb directory is assumed
	    }
	    
	    ### Check file existance and move created file if lacking 
	    my $intendedFilePathRef;
	    my $temporaryFilePathRef;
	    
	    if (defined(${$annovarTableHashRef}{ ${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter] }{'file'})) {
		
		for (my $filesCounter=0;$filesCounter<scalar(@{${$annovarTableHashRef}{ ${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter] }{'file'}});$filesCounter++) {  #All annovarTable file(s), some tables have multiple files downloaded from the same call
		    
		    $intendedFilePathRef = \(${$scriptParameterHashRef}{'annovarPath'}."/humandb/".${$annovarTableHashRef}{ ${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter] }{'file'}[$filesCounter]);  
		    $temporaryFilePathRef = \($annovarTemporaryDirectory."/".${$annovarTableHashRef}{ ${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter] }{'file'}[$filesCounter]);
	
		    ## Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
		    &PrintCheckExistandMoveFile($FILEHANDLE, $intendedFilePathRef, $temporaryFilePathRef);
		
		    if (defined(${$annovarTableHashRef}{ ${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter] }{'indexFile'})) {
		    
			$intendedFilePathRef = \(${$scriptParameterHashRef}{'annovarPath'}."/humandb/".${$annovarTableHashRef}{ ${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter] }{'file'}[$filesCounter].".idx");  
			$temporaryFilePathRef = \($annovarTemporaryDirectory."/".${$annovarTableHashRef}{ ${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter] }{'file'}[$filesCounter].".idx");

			## Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
			&PrintCheckExistandMoveFile($FILEHANDLE, $intendedFilePathRef, $temporaryFilePathRef);	
		    }
		}		
	    }
	    elsif ((defined(${$annovarTableHashRef}{ ${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter] }{'ucscAlias'}))){
	    
		$intendedFilePathRef = \(${$scriptParameterHashRef}{'annovarPath'}."/humandb/".${$scriptParameterHashRef}{'annovarGenomeBuildVersion'}."_".${$annovarTableHashRef}{ ${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter] }{'ucscAlias'}.".txt");
		$temporaryFilePathRef = \($annovarTemporaryDirectory."/".${$scriptParameterHashRef}{'annovarGenomeBuildVersion'}."_".${$annovarTableHashRef}{ ${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter] }{'ucscAlias'}.".txt");

		## Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
		&PrintCheckExistandMoveFile($FILEHANDLE, $intendedFilePathRef, $temporaryFilePathRef);
	    
		if (defined(${$annovarTableHashRef}{ ${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter] }{'indexFile'})) {
		    
		    $intendedFilePathRef = \(${$scriptParameterHashRef}{'annovarPath'}."/humandb/".${$scriptParameterHashRef}{'annovarGenomeBuildVersion'}."_".${$annovarTableHashRef}{ ${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter] }{'ucscAlias'}.".txt.idx");
		    $temporaryFilePathRef = \($annovarTemporaryDirectory."/".${$scriptParameterHashRef}{'annovarGenomeBuildVersion'}."_".${$annovarTableHashRef}{ ${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter] }{'ucscAlias'}.".txt.idx");

		    ## Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
		    &PrintCheckExistandMoveFile($FILEHANDLE, $intendedFilePathRef, $temporaryFilePathRef);
		}
	    }
	    else {
	    
		$intendedFilePathRef = \(${$scriptParameterHashRef}{'annovarPath'}."/humandb/".${$scriptParameterHashRef}{'annovarGenomeBuildVersion'}."_".${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter].".txt");
		$temporaryFilePathRef = \($annovarTemporaryDirectory."/".${$scriptParameterHashRef}{'annovarGenomeBuildVersion'}."_".${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter].".txt");    

		## Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
		&PrintCheckExistandMoveFile($FILEHANDLE, $intendedFilePathRef, $temporaryFilePathRef);
	    
		if (defined(${$annovarTableHashRef}{ ${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter] }{'indexFile'})) {
	
		    $intendedFilePathRef = \(${$scriptParameterHashRef}{'annovarPath'}."/humandb/".${$scriptParameterHashRef}{'annovarGenomeBuildVersion'}."_".${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter].".txt.idx");
		    $temporaryFilePathRef = \($annovarTemporaryDirectory."/".${$scriptParameterHashRef}{'annovarGenomeBuildVersion'}."_".${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter].".txt.idx");    

		    ## Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
		    &PrintCheckExistandMoveFile($FILEHANDLE, $intendedFilePathRef, $temporaryFilePathRef);	
		}				
	    }
	}
        ${$parameterHashRef}{${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter]}{'buildFile'} = 0;
    }
    
    print $FILEHANDLE "rm -rf $annovarTemporaryDirectory;", "\n\n";  #Cleaning up temp directory
    close($FILEHANDLE);
    
    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'dependencies' => 6, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}


sub BuildDownLoadablePreRequisites {

##BuildDownLoadablePreRequisites

##Function : Creates the downloadable resources.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, supportedCosmidReferenceHashRef, $familyID, $aligner, $programName, $FILEHANDLE, $randomInteger
##         : $parameterHashRef                 => The parameter hash {REF}
##         : $scriptParameterHashRef           => The active parameters for this analysis hash {REF}
##         : $supportedCosmidReferenceHashRef  => The supported cosmid references hash {REF}
##         : $familyID                         => Family ID
##         : $aligner                          => The aligner used in the analysis
##         : $programName                      => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $supportedCosmidReferenceHashRef = $_[2];
    my $familyID = $_[3];
    my $aligner = $_[4];
    my $programName = $_[5];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    
    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					   {'directoryID' => $familyID,
					    'programName' => $programName,
					    'programDirectory' => lc($aligner),
					    'processTime' => 4,
					   });

    print $FILEHANDLE "cd ${$scriptParameterHashRef}{'referencesDir'}", "\n\n";  #Move to reference directory

    ## Locates and sets the cosmid directory to download to
    my $cosmidResourceDirectory = &CheckCosmidYAML(\%{$scriptParameterHashRef});

    for my $parameterName (keys %{$supportedCosmidReferenceHashRef}) {

	if (${$parameterHashRef}{$parameterName}{'associatedProgram'} =~/$programName/) {

	    if (${$parameterHashRef}{$parameterName}{'buildFile'} eq 1) {
	    
		&DownloadReference(\%parameter, \%scriptParameter, \%{$supportedCosmidReferenceHashRef}, \$cosmidResourceDirectory, \$programName, $FILEHANDLE, $parameterName);
	    }
	}
    }
    
    close($FILEHANDLE);
    
    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {
	
	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'dependencies' => 6, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}


sub BuildPTCHSMetricPreRequisites {

##BuildPTCHSMetricPreRequisites

##Function : Creates the target "infiles_list" "padded.infile_list" and interval_list files.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $referenceFileEndingsHashRef, $familyID, $aligner, $programName, $FILEHANDLE
##         : $parameterHashRef            => The parameter hash {REF}
##         : $scriptParameterHashRef      => The active parameters for this analysis hash {REF}
##         : $referenceFileEndingsHashRef => The associated reference file endings {REF}
##         : $familyID                    => Family ID
##         : $aligner                     => The aligner used in the analysis
##         : $programName                 => The program name
##         : $FILEHANDLE                  => Filehandle to write to

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $referenceFileEndingsHashRef = $_[2];
    my $familyID = $_[3];
    my $aligner = $_[4];
    my $programName = $_[5];   
    my $FILEHANDLE = $_[6];  #Decides if a new sbatch script will be generated or handled by supplied FILEHANDLE

    my $parametersToEvaluate = 0;  #The number of parameters to evaluate
    my $fileName;

    unless(defined($FILEHANDLE)) {  #No supplied FILEHANDLE i.e. create new sbatch script    
	
	$FILEHANDLE = IO::Handle->new();#Create anonymous filehandle

	## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
	($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					    {'directoryID' => $familyID,
					     'programName' => $programName,
					     'programDirectory' => lc($aligner),
					    });
    }

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{${$scriptParameterHashRef}{'sampleIDs'}});$sampleIDCounter++) {  #All sampleIDs

	my $sampleIDBuildSwitchInfile = ${$parameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }{'exomeTargetBedInfileLists'}{'buildFile'};
	my $sampleIDBuildFileInfile = ${$scriptParameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }{'exomeTargetBedInfileLists'};
	my $infileListNoEnding = &RemoveFileEnding(\$sampleIDBuildFileInfile, ${$referenceFileEndingsHashRef}{'exomeTargetBedInfileLists'});  #For comparison of identical filename.bed files, to avoid creating files twice

	my $sampleIDBuildSwitchPadded = ${$parameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }{'exomeTargetPaddedBedInfileLists'}{'buildFile'};
	my $sampleIDBuildFilePadded = ${$scriptParameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }{'exomeTargetPaddedBedInfileLists'};
	my $paddedInfileListNoEnding = &RemoveFileEnding(\$sampleIDBuildFilePadded , ${$referenceFileEndingsHashRef}{'exomeTargetPaddedBedInfileLists'});  #For comparison of identical filename.bed files, to avoid creating files twice

	my $sampleIDBuildSwitchPaddedInterval = ${$parameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }{'GATKTargetPaddedBedIntervalLists'}{'buildFile'};
	my $sampleIDBuildFilePaddedInterval = ${$scriptParameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }{'GATKTargetPaddedBedIntervalLists'};
	my $paddedIntervalListNoEnding = &RemoveFileEnding(\$sampleIDBuildFilePaddedInterval , ${$referenceFileEndingsHashRef}{'GATKTargetPaddedBedIntervalLists'});  #For comparison of identical filename.bed files, to avoid creating files twice

	if ( (defined($sampleIDBuildSwitchPadded)) && ($sampleIDBuildSwitchPadded == 1) ) {  #If identical filename.bed and padded files need creation do not build infile_list in separate part of sbatch

	    if ($infileListNoEnding eq $paddedInfileListNoEnding) {		
	
		$sampleIDBuildSwitchInfile = 0;  #Turn of separate infile_list creation
	    }	
	    $parametersToEvaluate = $sampleIDBuildSwitchPadded;  #Add to parameters to evaluate (1 or 0)
	}
	if ( (defined($sampleIDBuildSwitchPaddedInterval)) && ($sampleIDBuildSwitchPaddedInterval == 1) ) {  #If identical filename.bed and padded files need creation do not build .pad100.infile_list in separate part of sbatch

	    if ($paddedInfileListNoEnding eq $paddedIntervalListNoEnding) {
		
		$sampleIDBuildSwitchPadded = 0;  #Turn of seperate paddded infile_list creation
	    }	
	    $parametersToEvaluate = $parametersToEvaluate + $sampleIDBuildSwitchPaddedInterval;   #Add to parameters to evaluate (1 or 0)
	} 
	if (defined($sampleIDBuildSwitchInfile)) {

	    $parametersToEvaluate = $parametersToEvaluate + $sampleIDBuildSwitchInfile;   #Add to parameters to evaluate (1 or 0)
	}
	
        ##Turn of build of identical filename.bed files
	&CheckUniqueTargetFiles(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \@{${$scriptParameterHashRef}{'sampleIDs'}}, \$sampleIDCounter, \$sampleIDBuildFileInfile, "exomeTargetBedInfileLists");
	&CheckUniqueTargetFiles(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \@{${$scriptParameterHashRef}{'sampleIDs'}}, \$sampleIDCounter, \$sampleIDBuildFilePadded, "exomeTargetPaddedBedInfileLists");
	&CheckUniqueTargetFiles(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \@{${$scriptParameterHashRef}{'sampleIDs'}}, \$sampleIDCounter, \$sampleIDBuildFilePaddedInterval, "GATKTargetPaddedBedIntervalLists");
	
	for (my $parameterCounter=0;$parameterCounter<$parametersToEvaluate;$parameterCounter++) {

	    my $randomInteger = int(rand(10000));  #Generate a random integer between 0-10,000.		
	    
	    ##Initiate general build variables used for all parameters
	    my $sampleIDBuildFile;
	    my $sampleIDBuildFileNoEnding;
	    my $sampleIDBuildFileNoEndingTemp;
	    
	    if ( (defined($sampleIDBuildSwitchInfile)) && ($sampleIDBuildSwitchInfile eq 1) ) {
		
		&SetTargetFileGeneralBuildParameter(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \%referenceFileEndings, \$sampleIDBuildFileInfile, "exomeTargetBedInfileLists", \$sampleIDBuildFile, \$sampleIDBuildFileNoEnding, \${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter]);
	    }
	    elsif ( (defined($sampleIDBuildSwitchPadded)) && ($sampleIDBuildSwitchPadded eq 1) ) {

		&SetTargetFileGeneralBuildParameter(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \%referenceFileEndings, \$sampleIDBuildFilePadded, "exomeTargetPaddedBedInfileLists", \$sampleIDBuildFile, \$sampleIDBuildFileNoEnding, \${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter]);
	    }
	    elsif ( (defined($sampleIDBuildSwitchPaddedInterval)) && ($sampleIDBuildSwitchPaddedInterval == 1) ) {
		
		&SetTargetFileGeneralBuildParameter(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \%referenceFileEndings, \$sampleIDBuildFilePaddedInterval, "GATKTargetPaddedBedIntervalLists", \$sampleIDBuildFile, \$sampleIDBuildFileNoEnding, \${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter]);
	    }
	    
	    if (defined($sampleIDBuildFile)) {
		
		$sampleIDBuildFileNoEndingTemp = $sampleIDBuildFileNoEnding."_".$randomInteger;  #Add random integer	
		
		$logger->warn("Will try to create required ".$sampleIDBuildFile." file before executing ".$programName."\n");
		
		print $FILEHANDLE "#SampleID:".${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter], "\n\n";
		print $FILEHANDLE "#CreateSequenceDictionary from reference", "\n";
		print $FILEHANDLE "java -Xmx2g ";

		&WriteUseLargePages($FILEHANDLE, \${$scriptParameterHashRef}{'javaUseLargePages'});
		
		print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'picardToolsPath'}."/CreateSequenceDictionary.jar ";
		print $FILEHANDLE "R=".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference genome
		print $FILEHANDLE "OUTPUT=".${$scriptParameterHashRef}{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict ", "\n\n";  #Output sequence dictionnary
		
		print $FILEHANDLE "#Add target file to headers from sequenceDictionary", "\n";
		print $FILEHANDLE "cat ";  #Concatenate
		print $FILEHANDLE ${$scriptParameterHashRef}{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict ";  #Sequence dictionnary
		print $FILEHANDLE ${$scriptParameterHashRef}{'referencesDir'}."/".$sampleIDBuildFileNoEnding." ";  #Bed file
		print $FILEHANDLE "> ";  #Write to
		print $FILEHANDLE ${$scriptParameterHashRef}{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict_body", "\n\n";  #Add bed body to dictionnary
		
		print $FILEHANDLE "#Remove target annotations, 'track', 'browse' and keep only 5 columns", "\n";
		print $FILEHANDLE q?perl  -nae 'if ($_=~/@/) {print $_;} elsif ($_=~/^track/) {} elsif ($_=~/^browser/) {} else {print @F[0], "\t", (@F[1] + 1), "\t", @F[2], "\t", "+", "\t", "-", "\n";}' ?;
		print $FILEHANDLE ${$scriptParameterHashRef}{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict_body ";  #Infile
		print $FILEHANDLE "> ";  #Write to
		print $FILEHANDLE ${$scriptParameterHashRef}{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict_body_col_5", "\n\n";  #Remove unnecessary info and reformat 
		
		print $FILEHANDLE "#Create".${$referenceFileEndingsHashRef}{'exomeTargetBedInfileLists'}, "\n";
		print $FILEHANDLE "java -Xmx2g ";

		&WriteUseLargePages($FILEHANDLE, \${$scriptParameterHashRef}{'javaUseLargePages'});

		print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'picardToolsPath'}."/IntervalListTools.jar ";
		print $FILEHANDLE "INPUT=".${$scriptParameterHashRef}{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict_body_col_5 ";
		print $FILEHANDLE "OUTPUT=".${$scriptParameterHashRef}{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict_body_col_5_".${$referenceFileEndingsHashRef}{'exomeTargetBedInfileLists'}." ", "\n\n";
		    
		my $intendedFilePathRef = \(${$scriptParameterHashRef}{'referencesDir'}."/".$sampleIDBuildFileNoEnding.${$referenceFileEndingsHashRef}{'exomeTargetBedInfileLists'});
		my $temporaryFilePathRef = \(${$scriptParameterHashRef}{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict_body_col_5_".${$referenceFileEndingsHashRef}{'exomeTargetBedInfileLists'});

		## Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
		&PrintCheckExistandMoveFile($FILEHANDLE, $intendedFilePathRef, $temporaryFilePathRef);
	    }
	    if ( (defined($sampleIDBuildSwitchPadded) && ($sampleIDBuildSwitchPadded eq 1)) || (defined($sampleIDBuildSwitchPaddedInterval) && ($sampleIDBuildSwitchPaddedInterval eq 1)) ) {
		
		print $FILEHANDLE "#Create padded interval list", "\n";
		print $FILEHANDLE "java -Xmx2g ";

		&WriteUseLargePages($FILEHANDLE, \${$scriptParameterHashRef}{'javaUseLargePages'});
		
		print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'picardToolsPath'}."/IntervalListTools.jar ";
		print $FILEHANDLE "PADDING=100 ";  #Add 100 nt on both sides of bed entry
		print $FILEHANDLE "INPUT=".${$scriptParameterHashRef}{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict_body_col_5 ";
		print $FILEHANDLE "OUTPUT=".${$scriptParameterHashRef}{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict_body_col_5".${$referenceFileEndingsHashRef}{'exomeTargetPaddedBedInfileLists'}." ", "\n\n";
		
		my $intendedFilePathRef = \(${$scriptParameterHashRef}{'referencesDir'}."/".$sampleIDBuildFileNoEnding.${$referenceFileEndingsHashRef}{'exomeTargetPaddedBedInfileLists'});
		my $temporaryFilePathRef = \(${$scriptParameterHashRef}{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict_body_col_5".${$referenceFileEndingsHashRef}{'exomeTargetPaddedBedInfileLists'});    

		## Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
		&PrintCheckExistandMoveFile($FILEHANDLE, $intendedFilePathRef, $temporaryFilePathRef);
		
		if (defined($sampleIDBuildSwitchPaddedInterval) && ($sampleIDBuildSwitchPaddedInterval eq 1)) {
		    
		    ##Softlink '.interval_list' to padded .infile_list", "\n";
		    print $FILEHANDLE "ln -s ";  #Softlink
		    print $FILEHANDLE ${$scriptParameterHashRef}{'referencesDir'}."/".$sampleIDBuildFileNoEnding.${$referenceFileEndingsHashRef}{'exomeTargetPaddedBedInfileLists'}." ";  #Origin file
		    print $FILEHANDLE ${$scriptParameterHashRef}{'referencesDir'}."/".$sampleIDBuildFileNoEnding.${$referenceFileEndingsHashRef}{'GATKTargetPaddedBedIntervalLists'};  #interval_list file
		}
		
		print $FILEHANDLE "\n\n";
	    }
	    if (defined($sampleIDBuildFile)) {
		
		print $FILEHANDLE "#Remove temporary files", "\n";
		
		print $FILEHANDLE "rm ";
		print $FILEHANDLE ${$scriptParameterHashRef}{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict_body_col_5 ", "\n\n";
		
		print $FILEHANDLE "rm ";
		print $FILEHANDLE ${$scriptParameterHashRef}{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict_body ", "\n\n";
		
		print $FILEHANDLE "rm ";
		print $FILEHANDLE ${$scriptParameterHashRef}{'referencesDir'}."/".$sampleIDBuildFileNoEndingTemp.".dict ", "\n\n";
		
		if ( (defined($sampleIDBuildSwitchPadded)) && ($sampleIDBuildSwitchPadded == 0) ) {
		    
		    $sampleIDBuildSwitchPaddedInterval = 0;
		    &SetTargetFileGeneralBuildParameter(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \%referenceFileEndings, \$sampleIDBuildFilePaddedInterval, "GATKTargetPaddedBedIntervalLists", \$sampleIDBuildFile, \$sampleIDBuildFileNoEnding, \${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter]);
		}
		if ( (defined($sampleIDBuildSwitchInfile)) && ($sampleIDBuildSwitchInfile == 0) ) {
		    
		    $sampleIDBuildSwitchPadded = 0;
		    &SetTargetFileGeneralBuildParameter(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \%referenceFileEndings, \$sampleIDBuildFilePadded, "exomeTargetPaddedBedInfileLists", \$sampleIDBuildFile, \$sampleIDBuildFileNoEnding, \${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter]);
		}
		$sampleIDBuildSwitchInfile = 0;
		&SetTargetFileGeneralBuildParameter(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \%referenceFileEndings, \$sampleIDBuildFileInfile, "exomeTargetBedInfileLists", \$sampleIDBuildFile, \$sampleIDBuildFileNoEnding, \${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter]);
	    }
	}
    }
    unless($_[6]) {  #Unless FILEHANDLE was supplied close it and submit 
    
	close($FILEHANDLE);
    
	if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {
	    
	    &FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
			   {'dependencies' => 6, 
			    'path' => "MIP",
			    'sbatchFileName' => $fileName
			   });
	}
    }
}


sub BuildBwaPreRequisites {

##BuildBwaPreRequisites

##Function : Creates the BwaPreRequisites using scriptParameters{'humanGenomeReference'} as reference.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $bwaBuildReferenceFileEndingsArrayRef, $familyID, $aligner, $programName, $FILEHANDLE
##         : $parameterHashRef                     => The parameter hash {REF}
##         : $scriptParameterHashRef               => The active parameters for this analysis hash {REF}
##         : $bwaBuildReferenceFileEndingsArrayRef => The bwa reference associated file endings {REF}
##         : $familyID                             => Family ID
##         : $aligner                              => The aligner used in the analysis
##         : $programName                          => The program name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $bwaBuildReferenceFileEndingsArrayRef = $_[2];
    my $familyID = $_[3];
    my $aligner = $_[4];
    my $programName = $_[5];
    
    my $FILEHANDLE = IO::Handle->new();  #Create anonymous filehandle
    my $randomInteger = int(rand(10000));  #Generate a random integer between 0-10,000.

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					   {'directoryID' => $familyID,
					    'programName' => $programName,
					    'programDirectory' => lc($aligner),
					    'processTime' => 3,
					   });
 
    &BuildHumanGenomePreRequisites(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \%fileInfo, \@humanGenomeReferenceFileEndings, $familyID, $aligner, $programName, $FILEHANDLE, $randomInteger);

    if (${$parameterHashRef}{'bwaBuildReference'}{'buildFile'} eq 1) {

	$logger->warn("Will try to create required ".${$scriptParameterHashRef}{'humanGenomeReference'}." index files before executing ".$programName."\n");
	
	print $FILEHANDLE "#Building BWA index", "\n\n";
	print $FILEHANDLE "bwa index ";  #Index sequences in the FASTA format
	print $FILEHANDLE "-p ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}."_".$randomInteger." "; #Prefix of the index
	print $FILEHANDLE "-a bwtsw ";  #BWT construction algorithm
	print $FILEHANDLE ${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'},"\n\n";  #The FASTA reference sequences file
	
	for (my $fileEndingsCounter=0;$fileEndingsCounter<scalar(@{$bwaBuildReferenceFileEndingsArrayRef});$fileEndingsCounter++) {  #All fileEndings
	    
	    my $intendedFilePathRef = \(${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}.${$bwaBuildReferenceFileEndingsArrayRef}[$fileEndingsCounter]);
	    my $temporaryFilePathRef = \(${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}."_".$randomInteger.${$bwaBuildReferenceFileEndingsArrayRef}[$fileEndingsCounter]);    

	    ## Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
	    &PrintCheckExistandMoveFile($FILEHANDLE, $intendedFilePathRef, $temporaryFilePathRef);
	}
	${$parameterHashRef}{'bwaBuildReference'}{'buildFile'} = 0;  #Ensure that this subrutine is only executed once
    }
    close($FILEHANDLE);
    
    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'dependencies' => 6, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}


sub BuildMosaikAlignPreRequisites {

##BuildMosaikAlignPreRequisites
    
##Function : Creates the mosaikAlignPreRequisites using scriptParameters{'humanGenomeReference'} as reference.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $mosaikJumpDbStubFileEndingsArrayRef, $humanGenomeReferenceSourceRef, $humanGenomeReferenceVersionRef, $familyID, $aligner, $programName
##         : $parameterHashRef                    => The parameter hash {REF}
##         : $scriptParameterHashRef              => The active parameters for this analysis hash {REF}
##         : $mosaikJumpDbStubFileEndingsArrayRef => The mosaikJump database file endings
##         : $humanGenomeReferenceSourceRef       => The human genome source {REF}
##         : $humanGenomeReferenceVersionRef      => The human genome build version {REF}
##         : $familyID                            => Family ID
##         : $aligner                             => Aligner used in the analysis
##         : $programName                         => Program name 

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $mosaikJumpDbStubFileEndingsArrayRef = $_[2];
    my $humanGenomeReferenceSourceRef = $_[3];
    my $humanGenomeReferenceVersionRef = $_[4];
    my $familyID = $_[5];
    my $aligner = $_[6];
    my $programName = $_[7];

    my $FILEHANDLE = IO::Handle->new();#Create anonymous filehandle
    my $randomInteger = int(rand(10000));  #Generate a random integer between 0-10,000.

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header.
    my ($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					   {'directoryID' => $familyID,
					     'programName' => $programName,
					     'programDirectory' => lc($aligner),
					     'nrofCores' => 4,
					     'processTime' => 2,
					    });
    
    ## Creates the humanGenomePreRequisites using scriptParameters{'humanGenomeReference'} as reference.
    &BuildHumanGenomePreRequisites(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \%fileInfo, \@humanGenomeReferenceFileEndings, $familyID, $aligner, $programName, $FILEHANDLE, $randomInteger);

    if (${$parameterHashRef}{'mosaikAlignReference'}{'buildFile'} eq 1) {  ##Begin autoBuild of MosaikAlignReference
	
	$logger->warn("Will try to create required ".${$scriptParameterHashRef}{'mosaikAlignReference'}." before executing ".$programName."\n");
	
	print $FILEHANDLE "#Building MosaikAligner Reference", "\n\n";
	print $FILEHANDLE "MosaikBuild ";
	print $FILEHANDLE "-fr ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #The FASTA reference sequences file
	print $FILEHANDLE "-sn Homo_sapiens ";  #Species name
	print $FILEHANDLE "-ga ".$$humanGenomeReferenceSourceRef.$$humanGenomeReferenceVersionRef." ";  #The genome assembly ID
	print $FILEHANDLE "-oa ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'mosaikAlignReference'}."_".$randomInteger, "\n\n";  #Temporary outfile

	my $intendedFilePathRef = \(${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'mosaikAlignReference'});
	my $temporaryFilePathRef = \(${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'mosaikAlignReference'}."_".$randomInteger);    

	## Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
	&PrintCheckExistandMoveFile($FILEHANDLE, $intendedFilePathRef, $temporaryFilePathRef);
    }
    if (${$parameterHashRef}{'mosaikJumpDbStub'}{'buildFile'} eq 1) {  ##Begin autoBuild of MosaikJump Database

	$logger->warn("Will try to create required ".${$scriptParameterHashRef}{'mosaikJumpDbStub'}." before executing ".$programName."\n");

	print $FILEHANDLE "#Building MosaikAligner JumpDatabase", "\n\n";
	print $FILEHANDLE "mkdir -p /scratch/mosaik_tmp", "\n";
	print $FILEHANDLE "export MOSAIK_TMP=/scratch/mosaik_tmp", "\n\n";
	
	print $FILEHANDLE "MosaikJump ";
	print $FILEHANDLE "-ia ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'mosaikAlignReference'}." ";  #The input reference file  
	print $FILEHANDLE "-hs 15 ";  #The hash size
	print $FILEHANDLE "-mem 24 ";  #The amount memory used when sorting hashes
	print $FILEHANDLE "-out ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'mosaikJumpDbStub'}."_".$randomInteger, "\n\n";  #Mosaik JumpDbStub for the output filenames

	for (my $fileEndingsCounter=0;$fileEndingsCounter<scalar(@{$mosaikJumpDbStubFileEndingsArrayRef});$fileEndingsCounter++) {  #All MosaikJumpDb assocaiated files

	    my $intendedFilePathRef = \(${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'mosaikJumpDbStub'}.${$mosaikJumpDbStubFileEndingsArrayRef}[$fileEndingsCounter]);
	    my $temporaryFilePathRef = \(${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'mosaikJumpDbStub'}."_".$randomInteger.${$mosaikJumpDbStubFileEndingsArrayRef}[$fileEndingsCounter]);

	    ## Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
	    &PrintCheckExistandMoveFile($FILEHANDLE, $intendedFilePathRef, $temporaryFilePathRef);
	}	
	
	print $FILEHANDLE "rm -rf /scratch/mosaik_tmp", "\n\n";  #Cleaning up temp directory
    }
    close($FILEHANDLE);
    
    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	&FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
		       {'dependencies' => 6, 
			'path' => ${$parameterHashRef}{"p".$programName}{'chain'},
			'sbatchFileName' => $fileName
		       });
    }
}


sub CheckBuildHumanGenomePreRequisites {

##CheckBuildHumanGenomePreRequisites
    
##Function : Checks if the HumanGenomePreRequisites needs to be built
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $fileInfoHashRef, $humanGenomeReferenceFileEndingsArrayRef, $programName
##         : $parameterHashRef                        => The parameter hash {REF}
##         : $scriptParameterHashRef                  => The active parameters for this analysis hash {REF}
##         : $fileInfoHashRef                         => The fileInfo hash {REF}
##         : $humanGenomeReferenceFileEndingsArrayRef => The human genome associated file endings array {REF}
##         : $programName                             => Program name  

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $fileInfoHashRef = $_[2];
    my $humanGenomeReferenceFileEndingsArrayRef = $_[3];
    my $programName = $_[4];

    for (my $fileEndingsCounter=0;$fileEndingsCounter<scalar(@{$humanGenomeReferenceFileEndingsArrayRef});$fileEndingsCounter++) {  #Files assocaiated with human genome reference
	
	if ( (${$parameterHashRef}{"humanGenomeReference".${$humanGenomeReferenceFileEndingsArrayRef}[$fileEndingsCounter]}{'buildFile'} eq 1) || (${$fileInfoHashRef}{'humanGenomeCompressed'} eq "compressed") ) {
	   
	    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} != 1)) {
	
		## Creates the humanGenomePreRequisites using scriptParameters{'humanGenomeReference'} as reference.
		&BuildHumanGenomePreRequisites(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \%{$fileInfoHashRef}, \@humanGenomeReferenceFileEndings, ${$scriptParameterHashRef}{'familyID'}, ${$scriptParameterHashRef}{'aligner'}, $programName);
		last;#Will handle all metafiles build within sbatch script
	    }
	}
    }
    ##Collect sequence contigs from human reference
    #&CollectSeqContigs();  #Reloads if required NOTE:Preparation for future changes but not activated yet
}


sub CheckBuildPTCHSMetricPreRequisites {

##CheckBuildPTCHSMetricPreRequisites
    
##Function : Check if PicardToolsHSMetricsPrequisites needs to be built
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $programName, $FILEHANDLE
##         : $parameterHashRef       => The parameter hash {REF}
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $programName            => Program name
##         : $FILEHANDLE             => Filehandle to write to

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $programName = $_[2];
    my $FILEHANDLE = $_[3];

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{${$scriptParameterHashRef}{'sampleIDs'}});$sampleIDCounter++) {
	
	if (${$parameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }{'exomeTargetBedInfileLists'}{'buildFile'} eq 1) {
	    
	    &BuildPTCHSMetricPreRequisites(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \%referenceFileEndings, ${$scriptParameterHashRef}{'familyID'}, ${$scriptParameterHashRef}{'aligner'}, $programName, $FILEHANDLE);
	    last;  #Will handle all build per sampleID within sbatch script
	}
	if (${$parameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }{'exomeTargetPaddedBedInfileLists'}{'buildFile'} eq 1) {
	    
	    &BuildPTCHSMetricPreRequisites(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \%referenceFileEndings, ${$scriptParameterHashRef}{'familyID'}, ${$scriptParameterHashRef}{'aligner'}, $programName, $FILEHANDLE);
	    last;  #Will handle all build per sampleID within sbatch script
	}
	if ( (defined(${$parameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }{'GATKTargetPaddedBedIntervalLists'}{'buildFile'})) && (${$parameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }{'GATKTargetPaddedBedIntervalLists'}{'buildFile'} eq 1) ){
	    
	    &BuildPTCHSMetricPreRequisites(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \%referenceFileEndings, ${$scriptParameterHashRef}{'familyID'}, ${$scriptParameterHashRef}{'aligner'}, $programName, $FILEHANDLE);
	    last;  #Will handle all build per sampleID within sbatch script
	}
    }
}


sub DownloadReference {
     
##DownloadReference
    
##Function : Downloads reference(s) using the database download manager Cosmid.
##Returns  : ""
##Arguments: $scriptParameterHashRef, $parameterHashRef, $supportedCosmidReferenceHashRef, $cosmidResourceDirectoryRef, $programRef, $FILEHANDLE, $parameterName, $cosmidResourceDirectoryRef
##         : $parameterHashRef                 => The parameter hash {REF}
##         : $scriptParameterHashRef           => The active parameters for this analysis hash {REF}
##         : $supportedCosmidReferenceHashRef  => The supported cosmid references hash {REF}
##         : $cosmidResourceDirectoryRef       => Cosmid directory {REF}
##         : $programRef                       => Program under evaluation {REF}
##         : $FILEHANDLE                       => Filehandle to write to
##         : $parameterName                    => Parameter to use for download

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $supportedCosmidReferenceHashRef = $_[2];
    my $cosmidResourceDirectoryRef = $_[3];
    my $programRef = $_[4];
    my $FILEHANDLE = $_[5];
    my $parameterName = $_[6];

    if (${$parameterHashRef}{$parameterName}{'buildFile'} eq 1) {  #Reference need to be built a.k.a downloaded
	
	## Use $parameter instead of $scriptParameter to cater for annotation files that are arrays and not supplied as flag => value
	if (defined(${$scriptParameterHashRef}{$parameterName})) {

	    $logger->warn("Will try to download ".${$scriptParameterHashRef}{$parameterName}." before executing ".$$programRef."\n");
	}
	else {

	    $logger->warn("Will try to download ".$parameterName." before executing ".$$programRef."\n");
	}
	print $FILEHANDLE "workon ".${$scriptParameterHashRef}{'pythonVirtualEnvironment'}, "\n\n";  #Activate python environment

	print $FILEHANDLE "cosmid ";  #Database download manager
	print $FILEHANDLE "clone ";  #Clone resource
	print $FILEHANDLE ${$supportedCosmidReferenceHashRef}{$parameterName}{'cosmidName'};  #The actual reference

	unless (${$supportedCosmidReferenceHashRef}{$parameterName}{'version'} eq "latest") {  #Version to download

	    print $FILEHANDLE "#".${$supportedCosmidReferenceHashRef}{$parameterName}{'version'},
	}
	print $FILEHANDLE "\n\n"; 

	print $FILEHANDLE "deactivate ", "\n\n";  #Deactivate python environment
	
	## Check if reference comes decompressed or not
	if (${$supportedCosmidReferenceHashRef}{$parameterName}{'compressedSwitch'} eq "compressed") {

	    print $FILEHANDLE "gzip ";
	    print $FILEHANDLE "-d ";  #Decompress
	    print $FILEHANDLE $$cosmidResourceDirectoryRef."/".${$supportedCosmidReferenceHashRef}{$parameterName}{'cosmidName'}."/*.gz", "\n\n";
	}

	my $intendedFilePathRef;
	## Use $parameter instead of $scriptParameter to cater for annotation files that are arrays and not supplied as flag => value
	if (defined(${$scriptParameterHashRef}{$parameterName})) {
	  
	    $intendedFilePathRef = \(${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{$parameterName});
	}
	else {
	    
	    $intendedFilePathRef = \(${$scriptParameterHashRef}{'referencesDir'}."/".$parameterName);
	}
	my $temporaryFilePathRef = \($$cosmidResourceDirectoryRef."/".${$supportedCosmidReferenceHashRef}{$parameterName}{'cosmidName'}."/*");

	## Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
	&PrintCheckExistandMoveFile($FILEHANDLE, $intendedFilePathRef, $temporaryFilePathRef);
	
	## Remove temporary Cosmid resources directory
	print $FILEHANDLE "rm -rf ";
	print $FILEHANDLE $$cosmidResourceDirectoryRef."/".${$supportedCosmidReferenceHashRef}{$parameterName}{'cosmidName'}."/;", "\n\n";

	## Remove temporary Cosmid ".cosmid.yaml" file
	print $FILEHANDLE "rm ";
	print $FILEHANDLE $$cosmidResourceDirectoryRef."/.cosmid.yaml", "\n\n";
	
	for my $supportedParameterName (keys %supportedCosmidReference) {

	    if (${$supportedCosmidReferenceHashRef}{$supportedParameterName}{'cosmidName'} eq ${$supportedCosmidReferenceHashRef}{$parameterName}{'cosmidName'}) {  #Reset to 0 for all supportedCosmidReference that are shared between modules 
		
		${$parameterHashRef}{$supportedParameterName}{'buildFile'} = 0;  #Only need to download once per analysis call
	    }
	}
    }
}


sub BuildHumanGenomePreRequisites {

##BuildHumanGenomePreRequisites
    
##Function : Creates the humanGenomePreRequisites using scriptParameters{'humanGenomeReference'} as reference.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $fileInfoHashRef, $humanGenomeReferenceFileEndingsArrayRef, $familyID, $aligner, $program, $FILEHANDLE, $randomInteger
##         : $parameterHashRef                        => The parameter hash {REF}
##         : $scriptParameterHashRef                  => The active parameters for this analysis hash {REF}
##         : $fileInfoHashRef                         => The fileInfo hash {REF}
##         : $humanGenomeReferenceFileEndingsArrayRef => The human genome associated file endings array {REF}
##         : $familyID                                => Family ID
##         : $aligner                                 => The aligner used in the analysis
##         : $program                                 => The program under evaluation
##         : $FILEHANDLE                              => Filehandle to write to. A new sbatch script will be generated if $FILEHANDLE is lacking, else write to exising $FILEHANDLE
##         : $randomInteger                           => The random integer to create temporary file name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $fileInfoHashRef = $_[2];
    my $humanGenomeReferenceFileEndingsArrayRef = $_[3];
    my $familyID = $_[4];
    my $aligner = $_[5];
    my $program = $_[6];
    my $FILEHANDLE = $_[7];
    my $randomInteger = $_[8];

    my $fileName;

    unless(defined($FILEHANDLE)) {  #No supplied FILEHANDLE i.e. create new sbatch script

	$FILEHANDLE = IO::Handle->new();  #Create anonymous filehandle
	$randomInteger = int(rand(10000));  #Generate a random integer between 0-10,000.

	## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
	($fileName) = &ProgramPreRequisites(\%{$scriptParameterHashRef}, $FILEHANDLE,
					    {'directoryID' => $familyID,
					     'programName' => $program,
					     'programDirectory' => lc($aligner),
					    });
    }

    print $FILEHANDLE "cd ${$scriptParameterHashRef}{'referencesDir'}", "\n\n";  #Move to reference directory

    ## Locates and sets the cosmid directory to download to
    my $cosmidResourceDirectory = &CheckCosmidYAML(\%{$scriptParameterHashRef});

    &DownloadReference(\%parameter, \%{$scriptParameterHashRef}, \%supportedCosmidReference, \$cosmidResourceDirectory, \$program, $FILEHANDLE, "humanGenomeReference");

    ## Check for compressed files
    if (${$fileInfoHashRef}{'humanGenomeCompressed'} eq "compressed") {

	$logger->warn("Will try to decompres ".${$scriptParameterHashRef}{'humanGenomeReference'}." before executing ".$program."\n");

	print $FILEHANDLE "gzip ";
	print $FILEHANDLE "-d ";  #Decompress
	print $FILEHANDLE ${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}, "\n\n";
	${$scriptParameterHashRef}{'humanGenomeReference'} =~ s/.fasta.gz/.fasta/g;  #Replace the .fasta.gz ending with .fasta since this will execute before the analysis, hence changing the original file name ending from ".fastq" to ".fastq.gz".
	$logger->info("Set humanGenomeReference to: ".${$scriptParameterHashRef}{'humanGenomeReference'}, "\n");
	${$fileInfoHashRef}{'humanGenomeCompressedRef'} = "unCompressed";
    }

    &CheckBuildPTCHSMetricPreRequisites(\%parameter, \%{$scriptParameterHashRef}, $program, $FILEHANDLE);

    for (my $fileEndingsCounter=0;$fileEndingsCounter<scalar(@{$humanGenomeReferenceFileEndingsArrayRef});$fileEndingsCounter++) {  #All meta files    
	
	if (${$parameterHashRef}{"humanGenomeReference.dict"}{'buildFile'} eq 1) {  #.dict file

	   $logger->warn("Will try to create dict file for ".${$scriptParameterHashRef}{'humanGenomeReference'}." before executing ".$program."\n");
	    
	    print $FILEHANDLE "#CreateSequenceDictionary from reference", "\n";
	    print $FILEHANDLE "java -Xmx2g ";

	    &WriteUseLargePages($FILEHANDLE, \${$scriptParameterHashRef}{'javaUseLargePages'});

	    print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'picardToolsPath'}."/CreateSequenceDictionary.jar ";
	    print $FILEHANDLE "R=".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference genome
	    print $FILEHANDLE "OUTPUT=".${$scriptParameterHashRef}{'referencesDir'}."/".${$fileInfoHashRef}{'humanGenomeReferenceNameNoEnding'}."_".$randomInteger.".dict ", "\n\n";  #Output sequence dictionnary
	    
	    ## Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
	    &PrintCheckExistandMoveFile($FILEHANDLE, \(${$scriptParameterHashRef}{'referencesDir'}."/".${$fileInfoHashRef}{'humanGenomeReferenceNameNoEnding'}.".dict"), \(${$scriptParameterHashRef}{'referencesDir'}."/".${$fileInfoHashRef}{'humanGenomeReferenceNameNoEnding'}."_".$randomInteger.".dict"));
	    
	    ${$parameterHashRef}{"humanGenomeReference.dict"}{'buildFile'} = 0;  #Only create once

	}
	if (${$parameterHashRef}{"humanGenomeReference.fasta.fai"}{'buildFile'} eq 1) {

	    $logger->warn("Will try to create .fai file for ".${$scriptParameterHashRef}{'humanGenomeReference'}." before executing ".$program."\n");

	    print $FILEHANDLE "#Fai file from reference", "\n";
	    print $FILEHANDLE "ln -s ";  #Softlink
	    print $FILEHANDLE ${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference genome
	    print $FILEHANDLE ${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}."_".$randomInteger, "\n\n";  #Softlink to Reference genome
	    
	    print $FILEHANDLE "samtools faidx ";#index/extract FASTA
	    print $FILEHANDLE ${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}."_".$randomInteger, "\n\n";  #Softlink to Reference genome
	    
	    ## Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
	    &PrintCheckExistandMoveFile($FILEHANDLE, \(${$scriptParameterHashRef}{'referencesDir'}."/".${$fileInfoHashRef}{'humanGenomeReferenceNameNoEnding'}.".fasta.fai"), \(${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}."_".$randomInteger.".fai"));
	
	    print $FILEHANDLE "rm ";  #Remove softLink
	    print $FILEHANDLE ${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}."_".$randomInteger, "\n\n";  #Softlink to Reference genome
	    
	    ${$parameterHashRef}{"humanGenomeReference.fasta.fai"}{'buildFile'} = 0;  #Only create once	
	}
    }
    unless($_[7]) {  #Unless FILEHANDLE was supplied close it and submit 
	
	close($FILEHANDLE);
    
	if ( (${$scriptParameterHashRef}{"p".$program} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	    &FIDSubmitJob(\%{$scriptParameterHashRef}, \%jobID, \%infilesLaneNoEnding,
			   {'dependencies' => 6, 
			    'path' => "MIP",
			    'sbatchFileName' => $fileName
			   });
	}
    }
}


sub CheckCosmidInstallation {

##CheckCosmidInstallation
    
##Function : Check that a Cosmid installation exists
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $parameterNameRef
##         : $parameterHashRef       => The parameter hash {REF}
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $parameterNameRef       => Parameter that uses Cosmid

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $parameterNameRef = $_[2];
    
    if (${$parameterHashRef}{$$parameterNameRef}{'buildFile'} eq 1) {
	
	if (defined(${$scriptParameterHashRef}{'pythonVirtualEnvironment'})) {  #Use python virtualenv
	
	    $logger->info("Checking your Cosmid installation in preparation for download of ".${$scriptParameterHashRef}{$$parameterNameRef}."\n");
 
	    my $whichReturn = `source ~/.bash_profile; workon ${$scriptParameterHashRef}{'pythonVirtualEnvironment'};which cosmid;deactivate;`;
	    
	    if ($whichReturn eq "") {

		$logger->fatal("MIP uses cosmid to download ".${$scriptParameterHashRef}{$$parameterNameRef}." and MIP could not find a cosmid installation in your python virtualenvironment".${$scriptParameterHashRef}{'pythonVirtualEnvironment'}." ","\n"); 
		exit;
	    }
	    else {  #Test ok

		$logger->info("Found installation in ".$whichReturn);
	    }
	}
	else  {  #No python virtualenv
	
	    $logger->fatal("Cannot download".${$scriptParameterHashRef}{$$parameterNameRef}." without a '-pythonVirtualEnvironment'");
	    exit;
	}
    }
}


sub ReadPlinkPedigreeFile {

##ReadPlinkPedigreeFile
    
##Function : Reads familyID_pedigree file in PLINK format. Checks for pedigree data for allowed entries and correct format. Add data to sampleInfo depending on user info. 
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $referenceFileEndingsHashRef, $supportedCaptureKitHashRef, $filePath
##         : $parameterHashRef            => The parameter hash {REF}
##         : $scriptParameterHashRef      => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef           => Info on samples and family hash {REF}
##         : $referenceFileEndingsHashRef => The associated reference file endings {REF}
##         : $supportedCaptureKitHashRef  => The supported capture kits hash {REF}
##         : $filePath                    => The pedigree file 
###FORMAT: FamliyID\tSampleID\tFather\tMother\tSex(1=male; 2=female; other=unknown)\tPhenotype(-9 missing, 0 missing, 1 unaffected, 2 affected)..n

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $referenceFileEndingsHashRef = $_[3];
    my $supportedCaptureKitHashRef = $_[4];
    my $filePath = $_[5];
    
    my @pedigreeFileElements = ("FamilyID", "SampleID", "Father", "Mother", "Sex", "Phenotype", );
    my $familyID;
    my $sampleID;
    
    ## Determine if the user supplied info on array parameter
    my $userSampleIDsSwitch = &CheckUserInfoArrays(\%{$scriptParameterHashRef}, \@{${$parameterHashRef}{'sampleIDs'}{'value'}}, "sampleIDs");
    my $userExomeTargetBedInfileListsSwitch = &CheckUserInfoArrays(\%{$scriptParameterHashRef}, \@exomeTargetBedInfileLists, "exomeTargetBedInfileLists"); 
    my $userExomeTargetPaddedBedInfileListSwitch = &CheckUserInfoArrays(\%{$scriptParameterHashRef}, \@exomeTargetPaddedBedInfileLists, "exomeTargetPaddedBedInfileLists");
    my $userExomeTargetPaddedBedIntervalListSwitch = &CheckUserInfoArrays(\%{$scriptParameterHashRef}, \@GATKTargetPaddedBedIntervalLists, "GATKTargetPaddedBedIntervalLists");

    ## Defines which entries are allowed and links them to position.
    my %plinkPedigree = &DefinePlinkPedigree();  #Holds allowed entries and positions to be checked for Plink pedigree files

    open(my $PEDF, "<", $filePath) or $logger->logdie("Can't open '".$filePath."': ".$!."\n");    
     
    while (<$PEDF>) {
	
	chomp $_;  #Remove newline
	
	if ( ($. == 1) && ($_ =~/^\#/) ) {  #Header present overwrite @pedigreeFileElements with header info
	
	    @pedigreeFileElements = split("\t", $'); #'
	    next;
	}
	if (m/^\s+$/) {  # Avoid blank lines
            next;
        }
	if (m/^\#/) {  # Avoid "#"
            next;
        }		
	if ($_ =~/(\S+)/) {	
	    
	    my @lineInfo = split("\t",$_);  #Loads pedigree file info
	    
	    ##Need to parse familyID and sampleID separately since these have not been set yet
	    if ($lineInfo[0] =~/\S+/) {  #FamilyID

		$familyID = $lineInfo[0];
	    }
	    else {

		$logger->("File: ".$filePath." at line ".$.." cannot find FamilyID in column 1\n");
		exit;
	    }
	    if ($lineInfo[1] =~/\S+/) { #SampleID

		$sampleID = $lineInfo[1];		

		if ($userSampleIDsSwitch == 0) {
		    
		    push(@{${$scriptParameterHashRef}{'sampleIDs'}}, $lineInfo[1]);  #Save sampleid info
		}
	    }
	    else {

		$logger->fatal("File: ".$filePath." at line ".$.." cannot find SampleID in column 2\n");
		exit;
	    }
	    for (my $sampleElementsCounter=0;$sampleElementsCounter<scalar(@pedigreeFileElements);$sampleElementsCounter++) {  #All pedigreeFileElements
		
		if ( defined($lineInfo[$sampleElementsCounter]) && ($lineInfo[$sampleElementsCounter] =~/\S+/) ) {  #Check that we have an non blank entry
		    
		    ## Test element for being part of hash of array at supplied key.
		    my $foundElement =  &CheckEntryHashofArray(\%plinkPedigree, \$sampleElementsCounter, \$lineInfo[$sampleElementsCounter]);

		    if ($foundElement == 1) {  #Invalid element found in file

			$logger->fatal("Found illegal element: '".$lineInfo[$sampleElementsCounter]."' in column '".$sampleElementsCounter."' in pedigree file: '".$filePath."' at line '".$.."'\n");
			$logger->fatal("Please correct the entry before analysis.\n");
			$logger->fatal("\nMIP: Aborting run.\n\n");
			exit;
		    }
		    
		    my @elementInfo = split(";", $lineInfo[$sampleElementsCounter]);  #Split element (if required)

		    if ($sampleElementsCounter < 6) {  #Mandatory elements known to be key->value
			
			${$sampleInfoHashRef}{$familyID}{$sampleID}{$pedigreeFileElements[$sampleElementsCounter]} = $lineInfo[$sampleElementsCounter];
		    }	
		    else {  #Other elements treat as lists

			## Detects if there are elements in arrayQueryRef that are not present in scalarQueryRef or arrayToCheckRef. If unique adds the unique element to arrayToCheckRef.
			&CheckUniqueArrayElement(\@{ ${$sampleInfoHashRef}{$familyID}{$sampleID}{$pedigreeFileElements[$sampleElementsCounter]} }, \@elementInfo);  #Check if there are any new info and add it to sampleInfo if so. 
		    }
		    if (${$sampleInfoHashRef}{$familyID}{$sampleID}{'Capture_kit'} && $pedigreeFileElements[$sampleElementsCounter] eq "Capture_kit") {  #Add latest capture kit for each individual
			
			my $captureKit = ${$sampleInfoHashRef}{$familyID}{$sampleID}{$pedigreeFileElements[$sampleElementsCounter]}[-1];  #Use only the last capture kit since it should be the most interesting

			${$scriptParameterHashRef}{$familyID}{$sampleID}{'exomeTargetBedInfileLists'} = &AddCaptureKit(\%{$referenceFileEndingsHashRef}, \%{$supportedCaptureKitHashRef}, 
															{'captureKit' => $captureKit, 
															 'parameterName' => "exomeTargetBedInfileLists", 
															 'userSuppliedParameterswitch' => $userExomeTargetBedInfileListsSwitch});  #Capture kit target infile_list 
			${$scriptParameterHashRef}{$familyID}{$sampleID}{'exomeTargetPaddedBedInfileLists'} = &AddCaptureKit(\%{$referenceFileEndingsHashRef}, \%{$supportedCaptureKitHashRef}, 
															      {'captureKit' => $captureKit,
															       'parameterName' => "exomeTargetPaddedBedInfileLists",
															       'userSuppliedParameterswitch' => $userExomeTargetPaddedBedInfileListSwitch});  #Capture kit padded target infile_list	
			${$scriptParameterHashRef}{$familyID}{$sampleID}{'GATKTargetPaddedBedIntervalLists'} = &AddCaptureKit(\%{$referenceFileEndingsHashRef}, \%{$supportedCaptureKitHashRef}, 
															       {'captureKit' => $captureKit,
																'parameterName' => "GATKTargetPaddedBedIntervalLists",
																'userSuppliedParameterswitch' => $userExomeTargetPaddedBedIntervalListSwitch}); #Capture kit padded target interval_list
		    }
		}
		else {  #No entry in pedigre file element
		    
		    if ($sampleElementsCounter < 6) {  #Only check mandatory elements 

			$logger->fatal($pedigreeFileElements[$sampleElementsCounter], "\t File: ".$filePath." at line ".$.."\tcannot find '".$pedigreeFileElements[$sampleElementsCounter]."' entry in column ".$sampleElementsCounter, "\n");
			exit;
		    }  
		}
	    }
	}	
    }
    if ($userSampleIDsSwitch == 0) {

	@{${$scriptParameterHashRef}{'sampleIDs'}} = sort(@{${$scriptParameterHashRef}{'sampleIDs'}});  #Lexiographical sort to determine the correct order of ids indata
    }
    $logger->info("Read pedigree file: ".$filePath, "\n");
    close($PEDF);
}


sub DefinePlinkPedigree {

##DefinePlinkPedigree
    
##Function : Defines which entries are allowed and links them to position.
##Returns  : "%plinkPedigree"
##Arguments: 
##         : 

    my %plinkPedigree;

    $plinkPedigree{4} = [1, 2, "other"];  #Sex allowed entries
    $plinkPedigree{5} = [-9, 0, 1, 2];  #Phenotype allowed entries

    return %plinkPedigree
}


sub AddToJobID {

##AddToJobID
    
##Function : Adds all previous jobIds per familyChainKey and chainKey to jobIDs string used to set the dependency in SLURM.
##Returns  : "$jobIDs"
##Arguments: $jobIDHashRef, $familyIDChainKey, $chainKey
##         : $jobIDHashRef     => The info on jobIds hash {REF}
##         : $familyIDChainKey => Family ID chain hash key
##         : $chainKey         => The current chain hash key

    my $jobIDHashRef = $_[0];
    my $familyIDChainKey = $_[1];
    my $chainKey = $_[2];
    
    my $jobIDs = "";  #JobID string
    
    if (${$jobIDHashRef}{$familyIDChainKey}{$chainKey}) {
	
	for (my $jobCounter=0;$jobCounter<scalar( @{ ${$jobIDHashRef}{$familyIDChainKey}{$chainKey} });$jobCounter++) {   #All previous jobIDs
	    
	    if ( ($jobCounter == 0) && (scalar( @{ ${$jobIDHashRef}{$familyIDChainKey}{$chainKey} }) == 1) ) {  #Only 1 previous jobID 
		
		$jobIDs .= ":${$jobIDHashRef}{$familyIDChainKey}{$chainKey}[$jobCounter]";  #First and last jobID start with ":" and end without ":"
	    }
	    elsif ($jobCounter == 0) {  #First jobID
		
		$jobIDs .= ":${$jobIDHashRef}{$familyIDChainKey}{$chainKey}[$jobCounter]:";  #First jobID start with :
	    }
	    elsif ($jobCounter eq (scalar( @{ ${$jobIDHashRef}{$familyIDChainKey}{$chainKey} }) -1) ) {  #Last jobID
		
		$jobIDs .= "${$jobIDHashRef}{$familyIDChainKey}{$chainKey}[$jobCounter]";  #Last jobID finish without :
	    }
	    else {  #JobIDs in the middle
		
		$jobIDs .= "${$jobIDHashRef}{$familyIDChainKey}{$chainKey}[$jobCounter]:";
	    }
	}
    }
    return $jobIDs;
}


sub PushToJobID {

##PushToJobID
    
##Function : Saves JobId to the correct hash array depending on chaintype.
##Returns  : ""
##Arguments: $scriptParameterHashRef, $sampleInfoHashRef, $jobIDHashRef, $infilesLaneNoEndingHashRef, $familyIDChainKey, $sampleIDChainKey, $sampleID, $path, $chainKeyType
##         : $scriptParameterHashRef     => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef          => Info on samples and family hash {REF}
##         : $jobIDHashRef               => The info on jobIds hash {REF}
##         : $infilesLaneNoEndingHashRef => The infile(s) without the ".ending" {REF}
##         : $familyIDChainKey           => Family ID chain hash key
##         : $sampleIDChainKey           => Sample ID chain hash key
##         : $sampleID                   => Sample ID
##         : $path                       => Trunk or branch
##         : $chainKeyType               => "parallel", "merged" or "family_merged" (familyID_sampleID)

    my $scriptParameterHashRef = $_[0];
    my $sampleInfoHashRef = $_[1];
    my $jobIDHashRef = $_[2];
    my $infilesLaneNoEndingHashRef = $_[3];
    my $familyIDChainKey = $_[4]; 
    my $sampleIDChainKey = $_[5];
    my $sampleID = $_[6];
    my $path = $_[7];
    my $chainKeyType = $_[8];
    
    my $chainKey;
    
    if ($chainKeyType eq "parallel") {  #Push parallel jobs

	if (${$scriptParameterHashRef}{'analysisType'} eq "rapid" && ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pBwaMem'}{'sbatchBatchProcesses'}) {  #Rapid run

	    for (my $sbatchCounter=0;$sbatchCounter<${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pBwaMem'}{'sbatchBatchProcesses'};$sbatchCounter++) {  #Iterate over sbatch processes instead of infile(s)

		$chainKey = $sampleID."_".$chainKeyType."_".$path.$sbatchCounter;  #Set key

		if (${$jobIDHashRef}{$familyIDChainKey}{$chainKey}) {  #Job exists
		    
		    for (my $jobCounter=0;$jobCounter<scalar( @{ ${$jobIDHashRef}{$familyIDChainKey}{$chainKey} });$jobCounter++) {  #All previous jobs i.e. jobs in this case equals to infiles in number
			
			push ( @{ ${$jobIDHashRef}{$familyIDChainKey}{$sampleIDChainKey} }, ${$jobIDHashRef}{$familyIDChainKey}{$chainKey}[$jobCounter]);  #Add jobID to hash
		    }    
		}
	    }	
	    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'pBwaMem'}{'sbatchBatchProcesses'} = ();
	}
	else {

	    for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} });$infileCounter++) {  #All infiles
		
		$chainKey = $sampleID."_".$chainKeyType."_".$path.$infileCounter;  #Set key
		
		if (${$jobIDHashRef}{$familyIDChainKey}{$chainKey}) {  #Job exists
		    
		    for (my $jobCounter=0;$jobCounter<scalar( @{ ${$jobIDHashRef}{$familyIDChainKey}{$chainKey} });$jobCounter++) {  #All previous jobs i.e. jobs in this case equals to infiles in number
			
			push ( @{ ${$jobIDHashRef}{$familyIDChainKey}{$sampleIDChainKey} }, ${$jobIDHashRef}{$familyIDChainKey}{$chainKey}[$jobCounter]);  #Add jobID to hash
		    }    
		}
	    }
	}
    }
    elsif ( ($chainKeyType eq "merged") || ($chainKeyType eq "family_merged")  ) {  #Push merged jobs
	
	$chainKey = $familyIDChainKey."_".$sampleIDChainKey;  #Set key
	
	if (${$jobIDHashRef}{$familyIDChainKey}{$chainKey}) {  #Job exists
	    
	    for (my $jobCounter=0;$jobCounter<scalar( @{ ${$jobIDHashRef}{$familyIDChainKey}{$chainKey} });$jobCounter++) {  #All previous jobs i.e. jobs in this case equals to infiles in number
		
		if ($chainKeyType eq "family_merged") {  #Use $familyIDChainKey instead of $sampleIDChainKey
		    
		    push ( @{ ${$jobIDHashRef}{$familyIDChainKey}{$familyIDChainKey} }, ${$jobIDHashRef}{$familyIDChainKey}{$chainKey}[$jobCounter]);  #Add jobID hash
		}
		else {
		    push ( @{ ${$jobIDHashRef}{$familyIDChainKey}{$sampleIDChainKey} }, ${$jobIDHashRef}{$familyIDChainKey}{$chainKey}[$jobCounter]);  #Add jobID to hash
		}
	    }    
	}
    }
}


sub FIDSubmitJob {

##FIDSubmitJob
    
##Function : Submits all jobIDs to SLURM using SLURM dependencies. The trunk is the "MAIN path" and any subsequent splits into  branches "other paths" later is handled by adding relevant previous jobIDs to the new paths key in jobID{family_path_key} hash. The subroutine supports parallel job within each step and submission which do not leave any dependencies. Currently any path downstream of MAIN inherits the relevant previous jobIds, but it is not possible to merge to splited paths downstream of main to each other.
##Returns  : ""
##Arguments: $scriptParameterHashRef, $jobIDHashRef, $infilesLaneNoEndingHashRef
##         : $scriptParameterHashRef            => The active parameters for this analysis hash {REF}
##         : $jobIDHashRef                      => The info on jobIds hash {REF}
##         : $infilesLaneNoEndingHashRef        => The infile(s) without the ".ending" {REF}
##         : $argHashRef{'sampleID'}            => Sample id
##         : $argHashRef{'familyID'}            => Family id
##         : $argHashRef{'dependencies'}        => Job dependencies
##         : $argHashRef{'path'}                => Trunk or Branch part of chainkey
##         : $argHashRef{'sbatchFileName'}      => Sbatch filename to submit
##         : $argHashRef{'sbatchScriptTracker'} => Track the number of parallel processes (e.g. sbatch scripts for a module)

###dependencies
    
##-1 = Not dependent on earlier scripts, and are self cul-de-sâcs
##0 = Not dependent on earlier scripts
##1 = Dependent on earlier scripts (within sampleID_path or familyID_path)
##2 = Dependent on earlier scripts (within sampleID_path or familyID_path), but are self cul-de-sâcs. 
##3 = Dependent on earlier scripts and executed in parallel within step
##4 = Dependent on earlier scripts and parallel scripts and executed in parallel within step 
##5 = Dependent on earlier scripts both family and sample and adds to both familyID and sampleId jobs
##6 = Not dependent on earlier scripts and adds to sampleId jobs, but sbatch is processed at family level i.e. affects all sampleID jobs e.g. building a reference

    my $scriptParameterHashRef = shift;
    my $jobIDHashRef = shift;
    my $infilesLaneNoEndingHashRef = shift;
    my ($argHashRef) = @_;

    my %default = ('familyID' => ${$scriptParameterHashRef}{'familyID'},
	);
	#'sampleID' => 0,
    #'sbatchScriptTracker' => 0
    my $sampleIDChainKey;
    my $sampleIDParallelChainKey;
    my $familyIDParallelChainKey;

    &SetDefaultArg(\%{$argHashRef}, \%default);

    if (defined(${$argHashRef}{'sampleID'})) {

	$sampleIDChainKey = ${$argHashRef}{'sampleID'}."_".${$argHashRef}{'path'};  #Sample chainkey
    }
    if ( (defined(${$argHashRef}{'sbatchScriptTracker'})) ) {

	$familyIDParallelChainKey = ${$argHashRef}{'familyID'}."_parallel_".${$argHashRef}{'path'}.${$argHashRef}{'sbatchScriptTracker'};  #Family parallel chainkey

	if (defined(${$argHashRef}{'sampleID'})) {

	    $sampleIDParallelChainKey = ${$argHashRef}{'sampleID'}."_parallel_".${$argHashRef}{'path'}.${$argHashRef}{'sbatchScriptTracker'};  #Sample parallel chainkey
	}
    }
    my $jobIDs = "";  #Create string with all previous jobIDs
    my $jobIDsReturn;  #Return jobID
    #my $sampleIDChainKey = ${$argHashRef}{'sampleID'}."_".${$argHashRef}{'path'};  #Sample chainkey
    my $familyIDChainKey = ${$argHashRef}{'familyID'}."_".${$argHashRef}{'path'};  #Family chainkey
    #my $sampleIDParallelChainKey = ${$argHashRef}{'sampleID'}."_parallel_".${$argHashRef}{'path'}.${$argHashRef}{'sbatchScriptTracker'};  #Sample parallel chainkey
    #my $familyIDParallelChainKey = ${$argHashRef}{'familyID'}."_parallel_".${$argHashRef}{'path'}.${$argHashRef}{'sbatchScriptTracker'};  #Family parallel chainkey
    my $jobID;  #The jobID that is returned from submission
    
    if (${$argHashRef}{'dependencies'} == -1) {  #Initiate chain - No dependencies, lonely program "sapling"
	
	$jobIDsReturn = `sbatch ${$argHashRef}{'sbatchFileName'}`;  #No jobs have been run: submit
	($jobID) = ($jobIDsReturn =~ /Submitted batch job (\d+)/);  #Just submitted jobID
    }
    if (${$argHashRef}{'dependencies'} == 6) {  #Initiate chain - No dependencies, adds to all sampleID(s)
	
	$jobIDsReturn = `sbatch ${$argHashRef}{'sbatchFileName'}`;  #No jobs have been run: submit
	($jobID) = ($jobIDsReturn =~ /Submitted batch job (\d+)/);  #Just submitted jobID
	
	for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{${$scriptParameterHashRef}{'sampleIDs'}});$sampleIDCounter++) {

	    my $sampleIDChainKey =  ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter]."_".${$argHashRef}{'path'};
	    push ( @{ ${$jobIDHashRef}{$familyIDChainKey}{$sampleIDChainKey} }, $jobID);  #Add jobID to hash
	}
    }
    elsif (${$argHashRef}{'dependencies'} == 0) {  #Initiate chain - No dependencies, initiate Trunk (Main or other)
	
	$jobIDsReturn = `sbatch ${$argHashRef}{'sbatchFileName'}`;  #No jobs have been run: submit
	($jobID) = ($jobIDsReturn =~ /Submitted batch job (\d+)/);  #Just submitted jobID
	push ( @{ ${$jobIDHashRef}{$familyIDChainKey}{$sampleIDChainKey} }, $jobID);  #Add jobID to hash
    }
    else {  #Dependent on earlier scripts and/or parallel. JobIDs that do not leave dependencies do not get pushed to jobID hash
	
	if (defined(${$argHashRef}{'sampleID'})) {  #Check jobs within sampleID (exception if dependencies = 5) 
	    
	    if (${$argHashRef}{'dependencies'} == 5) {  #Add familyID_sampleID jobs to current sampleID chain
		
		&PushToJobID(\%{$scriptParameterHashRef}, \%sampleInfo, \%jobID, \%infilesLaneNoEnding, $familyIDChainKey, $sampleIDChainKey, ${$argHashRef}{'sampleID'}, ${$argHashRef}{'path'}, "merged");
	    }
	    if ( (${$argHashRef}{'dependencies'} == 1) || (${$argHashRef}{'dependencies'} == 2) ) {  #Not parallel jobs, but check if last job submission was parallel
		
		&PushToJobID(\%{$scriptParameterHashRef}, \%sampleInfo, \%jobID, \%infilesLaneNoEnding, $familyIDChainKey, $sampleIDChainKey, ${$argHashRef}{'sampleID'}, ${$argHashRef}{'path'}, "parallel");
	    }
	    if ( (defined(${$argHashRef}{'path'})) && (${$argHashRef}{'path'} eq "MAIN") ) {
		
		if ( (${$argHashRef}{'dependencies'} == 4) || (${$argHashRef}{'dependencies'} == 3) ) {  #Parallel jobs
		    
		    $jobIDs = &AddToJobID(\%jobID, $familyIDChainKey, $sampleIDParallelChainKey);  #Add to jobID string
		    
		    if (${$jobIDHashRef}{$familyIDChainKey}{$sampleIDChainKey}) {  #Check for previous single jobs - required to initiate broken chain with correct dependencies 
               
			$jobIDs .= &AddToJobID(\%jobID, $familyIDChainKey, $sampleIDChainKey);  #Add to jobID string
		    }
		    
		}
		else {  #Previous job was a single job
		    
		    $jobIDs = &AddToJobID(\%jobID, $familyIDChainKey, $sampleIDChainKey);  #Add to jobID string
		}
	    }
	    if ( (defined(${$argHashRef}{'path'})) && (${$argHashRef}{'path'} ne "MAIN") ) {  #Check for any previous jobIDs within path current PATH. Branch.
		
		if (${$jobIDHashRef}{$familyIDChainKey}{$sampleIDChainKey}) {  #Second or later in branch chain
		    
		    $jobIDs = &AddToJobID(\%jobID, $familyIDChainKey, $sampleIDChainKey);
		}
		elsif (${$jobIDHashRef}{${$argHashRef}{'familyID'}."_MAIN"}{${$argHashRef}{'sampleID'}."_MAIN"}) {  #Inherit from potential MAIN. Trunk
		    
		    $jobIDs = &AddToJobID(\%jobID, ${$argHashRef}{'familyID'}."_MAIN", ${$argHashRef}{'sampleID'}."_MAIN");
		}
	    }     
	    if ($jobIDs) {  #Previous jobs for chainkey exists

		$jobIDsReturn = `sbatch --dependency=afterok$jobIDs ${$argHashRef}{'sbatchFileName'}`;  #Supply with dependency of previous jobs that this one is dependent on
		($jobID) = ($jobIDsReturn =~ /Submitted batch job (\d+)/);  #Just submitted jobID
	    }
	    else {  #No previous jobs

		$jobIDsReturn = `sbatch ${$argHashRef}{'sbatchFileName'}`;  #No jobs have been run: submit
		($jobID) = ($jobIDsReturn =~ /Submitted batch job (\d+)/);  #Just submitted jobID
	    }
	    if (${$argHashRef}{'dependencies'} == 1) {  #Ordinary job push to array
		
		@{ ${$jobIDHashRef}{$familyIDChainKey}{$sampleIDChainKey} } = ();  #Clear latest familyID/sampleID chain submission
		
		##Clear all latest parallel jobs within chainkey
		for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infilesLaneNoEndingHashRef}{${$argHashRef}{'sampleID'}} });$infileCounter++) {
		    
		    my $sampleIDParallelChainKey = ${$argHashRef}{'sampleID'}."_parallel_".${$argHashRef}{'path'}.$infileCounter;  #Create key
		    
		    if (${$jobIDHashRef}{$familyIDChainKey}{$sampleIDParallelChainKey}) {  #Parallel job exists
			
			@{ ${$jobIDHashRef}{$familyIDChainKey}{$sampleIDParallelChainKey} } = ();  #Clear latest familyID/sampleID chain submission
                    }
		}
		
		push ( @{ ${$jobIDHashRef}{$familyIDChainKey}{$sampleIDChainKey} }, $jobID);  #Add jobID to hash{$sampleID}[]
	    }
	    if ( (${$argHashRef}{'dependencies'} == 3) || (${$argHashRef}{'dependencies'} == 4) ) {  #Parallel job wait to push to array until all parallel jobs are finished within step
		
		push ( @{ ${$jobIDHashRef}{$familyIDChainKey}{$sampleIDParallelChainKey} }, $jobID);  #Add jobID to hash
	    }
	    if (${$argHashRef}{'dependencies'} == 5) {  #Job dependent on both familyID and sampleID push to array
		
		@{ ${$jobIDHashRef}{$familyIDChainKey}{$familyIDChainKey."_".$sampleIDChainKey} } = ();  #Clear latest familyID_sampleID chainkey
		@{ ${$jobIDHashRef}{$familyIDChainKey}{$sampleIDChainKey} } = ();  #Clear latest sampleID chainkey
		push ( @{ ${$jobIDHashRef}{$familyIDChainKey}{$familyIDChainKey."_".$sampleIDChainKey} }, $jobID);  #Add jobID to hash
	    }
	}
	else {  #AFTER merging to familyID
	    
	    if (${$argHashRef}{'dependencies'} == 5) {  #Add familyID_sampleID jobs to current familyID chain
		
		for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{${$scriptParameterHashRef}{'sampleIDs'}});$sampleIDCounter++) {  #Check jobs for sampleID          
		    
		    my $sampleIDChainKey = ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter]."_".${$argHashRef}{'path'};  #Current chain
		    &PushToJobID(\%{$scriptParameterHashRef}, \%sampleInfo, \%jobID, \%infilesLaneNoEnding, $familyIDChainKey, $sampleIDChainKey, ${$argHashRef}{'sampleID'}, ${$argHashRef}{'path'}, "family_merged");
		}
	    }
	    if ( (${$argHashRef}{'dependencies'} == 1) || (${$argHashRef}{'dependencies'} == 2) ) {  #Not parallel jobs, but check if last job submission was parallel
		
		if ( (defined($familyIDParallelChainKey)) && (${$jobIDHashRef}{$familyIDChainKey}{$familyIDParallelChainKey})) {  #Parallel job exists
		    
		    for (my $jobCounter=0;$jobCounter<scalar( @{ ${$jobIDHashRef}{$familyIDChainKey}{$familyIDParallelChainKey} });$jobCounter++) {
			
			push ( @{ ${$jobIDHashRef}{$familyIDChainKey}{$familyIDChainKey} }, ${$jobIDHashRef}{$familyIDChainKey}{$familyIDParallelChainKey}[$jobCounter]);  #Add jobID to hash{$} 
		    }
		}
	    }
	    if ( (defined(${$argHashRef}{'path'})) && (${$argHashRef}{'path'} eq "MAIN") && (${$jobIDHashRef}{$familyIDChainKey}{$familyIDChainKey}) ) {  #Check for any previous jobIDs within path MAIN. Test for previous must be done to allow initiating from broken chain. Trunk and not first in chain
		if ( (${$argHashRef}{'dependencies'} == 4) || (${$argHashRef}{'dependencies'} == 3) ) {  #Parallel jobs
		    
		    $jobIDs = &AddToJobID(\%jobID, $familyIDChainKey, $familyIDParallelChainKey);  #Add to jobID string
		}
		else {  #Previous job was a single job 
		    
		    $jobIDs = &AddToJobID(\%jobID, $familyIDChainKey, $familyIDChainKey);  #Add to jobID string
		}
	    }
	    elsif ((defined(${$argHashRef}{'path'})) && ${$argHashRef}{'path'} eq "MAIN") {  #First familyID MAIN chain 
		
		##Add all previous jobId(s) from sampleID chainkey(s)
		for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{${$scriptParameterHashRef}{'sampleIDs'}});$sampleIDCounter++) {           
		    
		    my $sampleIDChainKey = ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter]."_".${$argHashRef}{'path'};
		    
		    if (${$jobIDHashRef}{$familyIDChainKey}{$sampleIDChainKey}) {
			
			$jobIDs .= &AddToJobID(\%jobID, $familyIDChainKey, $sampleIDChainKey);  #Add to jobID string, while keeping previous additions
			
		    }
		    for (my $infileCounter=0;$infileCounter<scalar( @{ ${$infilesLaneNoEndingHashRef}{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] } });$infileCounter++) {
			
			my $sampleIDParallelChainKey = ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter]."_parallel_".${$argHashRef}{'path'}.$infileCounter;  #Create key
			
			if (${$jobIDHashRef}{$familyIDChainKey}{$sampleIDParallelChainKey}) {  #Parallel job exists
			    
			    $jobIDs .= &AddToJobID(\%jobID, $familyIDChainKey, $sampleIDParallelChainKey);  #Add to jobID string, while keeping previous additions
			    
			}
		    }
		}
	    }
	    if ((defined(${$argHashRef}{'path'})) && ${$argHashRef}{'path'} ne "MAIN" ) {  #Check for any previous jobIDs within path current PATH. Branch
		
		if (${$jobIDHashRef}{$familyIDChainKey}{$familyIDChainKey}) {  #Second or later in branch chain
		    
		    $jobIDs = &AddToJobID(\%jobID, $familyIDChainKey, $familyIDChainKey);  #Family chain
		}
		elsif (${$jobIDHashRef}{${$argHashRef}{'familyID'}."_MAIN"}{${$argHashRef}{'familyID'}."_MAIN"}) {  #Inherit from potential MAIN. Trunk
		    
		    $jobIDs = &AddToJobID(\%jobID, ${$argHashRef}{'familyID'}."_MAIN", ${$argHashRef}{'familyID'}."_MAIN");
		}
		else {  #First job in new path and first familyID MAIN chain 
		    
		    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{${$scriptParameterHashRef}{'sampleIDs'}});$sampleIDCounter++) {           
			
			my $familyIDChainKey = ${$argHashRef}{'familyID'}."_MAIN";
			my $sampleIDChainKey = ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter]."_MAIN";
			
			if (${$jobIDHashRef}{$familyIDChainKey}{$sampleIDChainKey}) {
			    
			    $jobIDs .= &AddToJobID(\%jobID, $familyIDChainKey, $sampleIDChainKey); 
			}
		    }
		}
	    }
	    if ($jobIDs) {

		$jobIDsReturn = `sbatch --dependency=afterok$jobIDs ${$argHashRef}{'sbatchFileName'}`;  #Supply with dependency of previous jobs that this one is dependent on
		($jobID) = ($jobIDsReturn =~ /Submitted batch job (\d+)/);
	    }
	    else {

		$jobIDsReturn = `sbatch ${$argHashRef}{'sbatchFileName'}`;  #No jobs have been run: submit
		($jobID) = ($jobIDsReturn =~ /Submitted batch job (\d+)/);
	    }
	    if (${$argHashRef}{'dependencies'} == 1) {  #Ordinary job push to array
		
		@{ ${$jobIDHashRef}{$familyIDChainKey}{$familyIDChainKey} } = ();  #Clear latest familyID/sampleID chain submission

		if (defined($familyIDParallelChainKey)) {
		    
		    @{ ${$jobIDHashRef}{$familyIDChainKey}{$familyIDParallelChainKey} } = ();  #Clear latest familyID/sampleID chain submission
		}
		push ( @{ ${$jobIDHashRef}{$familyIDChainKey}{$familyIDChainKey} }, $jobID);  #Add jobID to hash{$sampleID}[]
	    }
	    if ( (${$argHashRef}{'dependencies'} == 3) || (${$argHashRef}{'dependencies'} == 4) ) {  #Parallel job wait to push to array until all parallel jobs are finished within step
		
		push ( @{ ${$jobIDHashRef}{$familyIDChainKey}{$familyIDParallelChainKey} }, $jobID);  #Add jobID to hash{$sampleID_parallel}[].
	    }    
	    if (${$argHashRef}{'dependencies'} == 5) {  #Job dependent on both familyID and sampleID push to array
		
		for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{${$scriptParameterHashRef}{'sampleIDs'}});$sampleIDCounter++) {  #Check jobs for sampleID          
		    
		    my $sampleIDChainKey = ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter]."_".${$argHashRef}{'path'};  #Current chain
		    @{ ${$jobIDHashRef}{$familyIDChainKey}{$familyIDChainKey."_".$sampleIDChainKey} } = ();
		    @{ ${$jobIDHashRef}{$familyIDChainKey}{$familyIDChainKey} } = ();  #Clear latest sampleID chainkey
		    push ( @{ ${$jobIDHashRef}{$familyIDChainKey}{$familyIDChainKey."_".$sampleIDChainKey} }, $jobID);   
		}
	    }
	}
    }
    if ($jobIDsReturn !~/\d+/) {  #Catch errors since, propper sbatch submission should only return numbers

	$logger->fatal($jobIDsReturn."\n");
	$logger->fatal("MIP: Aborting run.\n");
	exit;
    }
    $logger->info("Sbatch script submitted, job id: $jobID\n");
    $logger->info("To check status of job, please run \'squeue -j $jobID\'\n");
    $logger->info("To cancel job, please run \'scancel $jobID\'\n");
}


sub NrofCoresPerSbatch {

##NrofCoresPerSbatch
    
##Function : Set the number of cores to allocate per sbatch job.
##Returns  : "$nrCores"
##Arguments: $scriptParameterHashRef, $nrCores
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $nrCores                => The number of cores to allocate

    my $scriptParameterHashRef = $_[0];
    my $nrCores = $_[1];
    
    if ($nrCores > ${$scriptParameterHashRef}{'maximumCores'}) {  #Set number of cores depending on how many lanes to process
	
	$nrCores = ${$scriptParameterHashRef}{'maximumCores'};  #Set to max on cluster
    }
    return $nrCores;
}


sub CollectInfiles {

##CollectInfiles
    
##Function : Collects the ".fastq(.gz)" files from the supplied infiles directory. Checks if any files exist.
##Returns  : ""
##Arguments: $scriptParameterHashRef, $inDirPathHashRef, $infileHashRef
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $inDirPathHashRef       => The indirectories path(s) hash {REF}
##         : $infileHashRef          => The infiles hash {REF}

    my $scriptParameterHashRef = $_[0];
    my $inDirPathHashRef = $_[1];
    my $infileHashRef = $_[2];

    $logger->info("Reads from Platform\n");

    for (my $inputDirectoryCounter=0;$inputDirectoryCounter<scalar(@{${$scriptParameterHashRef}{'sampleIDs'}});$inputDirectoryCounter++) {  #Collects inputfiles govern by sampleIDs
	
	my @infiles = `cd ${$scriptParameterHashRef}{'inFilesDirs'}[ $inputDirectoryCounter ];ls *.fastq*;`;  #'cd' to input dir and collect fastq files and fastq.gz files
	chomp(@infiles);   #Remove newline from every entry in array

	if (scalar(@infiles) == 0) {  #No "*.fastq*" infiles
	    
	    $logger->fatal("Could not find any '.fastq' files in supplied infiles directory ".${$scriptParameterHashRef}{'inFilesDirs'}[$inputDirectoryCounter], "\n");
	    exit;
	}
	foreach my $infile (@infiles) {  #Check that inFileDirs/infile contains sampleID in filename

	    unless ( $infile =~/${$scriptParameterHashRef}{'sampleIDs'}[$inputDirectoryCounter]/) {

		$logger->fatal("Could not detect sampleID: ".${$scriptParameterHashRef}{'sampleIDs'}[$inputDirectoryCounter]." in supplied infile: ".${$scriptParameterHashRef}{'inFilesDirs'}[$inputDirectoryCounter]."/".$infile, "\n");
		$logger->fatal("Check that the order of supplied: '--sampleIDs' and '--inFilesDirs' correlate.", "\n");
		$logger->fatal("NOTE: SampleIDs read from pedigree are lexiographically sorted and for instance '--inFileDirs' supplied need to be supplied in the same order to correlate", "\n");
		exit;
	    }
	}
	$logger->info("Sample ID: ".${$scriptParameterHashRef}{'sampleIDs'}[$inputDirectoryCounter]."\n");
	$logger->info("\tInputfiles:\n");
	## Log each file from platform
	foreach my $file (@infiles) {

	    $logger->info("\t\t", $file, "\n");  #Indent for visability
	}
	${$inDirPathHashRef}{ ${$scriptParameterHashRef}{'sampleIDs'}[$inputDirectoryCounter] } = ${$scriptParameterHashRef}{'inFilesDirs'}[ $inputDirectoryCounter ];   #Catch inputdir path
	${$infileHashRef}{ ${$scriptParameterHashRef}{'sampleIDs'}[$inputDirectoryCounter] }  = [@infiles];  #Reload files into hash
    }
}


sub InfilesReFormat {

##InfilesReFormat
    
##Function : Reformat files for MIP output, which have not yet been created into, correct format so that a sbatch script can be generated with the correct filenames.
##Returns  : "$uncompressedFileCounter"
##Arguments: $infileHashRef
##         : $infileHashRef => The infiles hash {REF}

    my $infileHashRef = $_[0];

    my $uncompressedFileCounter = 0;  #Used to decide later if any inputfiles needs to be compressed before starting analysis 

    for my $sampleID (keys %{$infileHashRef}) {  #For every sampleID                                                                                       
	
        my $laneTracker=0;  #Needed to be able to track when lanes are finished                                                                  
	
        for (my $infileCounter=0;$infileCounter<scalar( @ { ${$infileHashRef}{$sampleID} });$infileCounter++) {  #All inputfiles for all fastq dir and remakes format
	    
            if (${$infileHashRef}{$sampleID}[$infileCounter] =~ /\/?([^\.\/]+\.[^\.]+)\.lane(\d+)_([12FfRr])\.fastq.gz/) {  #Parse fastq.gz 'old' format, $2="lane", $3="Read direction"

		&AddInfileInfoOld(\%scriptParameter, \%lane, \%infilesLaneNoEnding, \%sampleInfo, \%infilesBothStrandsNoEnding, \%{$infileHashRef}, \%inDirPath, $1, $2, $3, $sampleID, \$laneTracker, $infileCounter, "compressed"); 
            }
            elsif (${$infileHashRef}{$sampleID}[$infileCounter] =~ /\/?([^\.\/]+\.[^\.]+)\.lane(\d+)_([12FfRr])\.fastq/) {  #Parse 'old' format

		&AddInfileInfoOld(\%scriptParameter, \%lane, \%infilesLaneNoEnding, \%sampleInfo, \%infilesBothStrandsNoEnding, \%{$infileHashRef}, \%inDirPath, $1, $2, $3, $sampleID, \$laneTracker, $infileCounter, "uncompressed");
                $uncompressedFileCounter = 1;  #File needs compression before starting analysis                               
            }
	    elsif (${$infileHashRef}{$sampleID}[$infileCounter] =~/(\d+)_(\d+)_([^_]+)_([^_]+)_index([^_]+)_(\d).fastq/) {  #Parse 'new' no "index" format $1=lane, $2=date, $3=Flow-cell, $4=SampleID, $5=index,$6=direction                             
		
		## Check if a file is gzipped.
		my $compressedSwitch = &CheckGzipped(\${$infileHashRef}{$sampleID}[$infileCounter]);#Check gzipped or not
		
		if ($compressedSwitch eq "unCompressed") {
		    
		    $uncompressedFileCounter = "unCompressed";  #File needs compression before starting analysis. Note: All files are rechecked downstream and uncompressed ones are gzipped automatically           
		}
		## Check that the sampleID provided and sampleID in infile name match.
		&CheckSampleIDMatch(\%scriptParameter, \%{$infileHashRef}, $sampleID, $4, $infileCounter);   #$4 = SampleID from filename
		&AddInfileInfo(\%scriptParameter, \%lane, \%infilesLaneNoEnding, \%sampleInfo, \%infilesBothStrandsNoEnding, \%{$infileHashRef}, \%inDirPath, $1, $2, $3, $4, $5, $6, \$laneTracker, $infileCounter, $compressedSwitch);		
	    }
            elsif (${$infileHashRef}{$sampleID}[$infileCounter] =~/(\d+)_(\d+)_([^_]+)_([^_]+)_([^_]+)_(\d).fastq/) {  #Parse 'new' no "index" format $1=lane, $2=date, $3=Flow-cell, $4=SampleID, $5=index,$6=direction                             
		
		## Check if a file is gzipped.
		my $compressedSwitch = &CheckGzipped(\${$infileHashRef}{$sampleID}[$infileCounter]);#Check gzipped or not

		if ($compressedSwitch eq "unCompressed") {
		   
		    $uncompressedFileCounter = "unCompressed";  #File needs compression before starting analysis. Note: All files are rechecked downstream and uncompressed ones are gzipped automatically           
		}
		## Check that the sampleID provided and sampleID in infile name match.
		&CheckSampleIDMatch(\%scriptParameter, \%{$infileHashRef}, $sampleID, $4, $infileCounter);  #$4 = SampleID from filename
		&AddInfileInfo(\%scriptParameter, \%lane, \%infilesLaneNoEnding, \%sampleInfo, \%infilesBothStrandsNoEnding, \%{$infileHashRef}, \%inDirPath, $1, $2, $3, $4, $5, $6, \$laneTracker, $infileCounter, $compressedSwitch);		
	    }
	    else {  #No regexp match i.e. file does not follow filename convention 

		$logger->fatal("Could not detect MIP file name convention for file: ".${$infileHashRef}{$sampleID}[$infileCounter].". \n");
		$logger->fatal("Please check that the file name follows the specified convention.", "\n");
		exit;
	    }
        }
    }
    return $uncompressedFileCounter;
}


sub CheckSampleIDMatch {

##CheckSampleIDMatch
    
##Function : Check that the sampleID provided and sampleID in infile name match.
##Returns  : ""
##Arguments: $scriptParameterHashRef, $infileHashRef, $sampleID, $infileSampleID, $infileCounter
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $infileHashRef          => The infiles hash {REF}
##         : $sampleID               => Sample id from user
##         : $infileSampleID         => SampleID collect with regexp from infile
##         : $infileCounter          => Counts the number of infiles

    my $scriptParameterHashRef = $_[0];
    my $infileHashRef = $_[1];
    my $sampleID = $_[2];
    my $infileSampleID = $_[3];
    my $infileCounter = $_[4];
    
    my %seen;
    $seen{$infileSampleID} = 1;  #Add input as first increment
    
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{${$scriptParameterHashRef}{'sampleIDs'}});$sampleIDCounter++) {
	
	$seen{${$scriptParameterHashRef}{'sampleIDs'}[ $sampleIDCounter]}++;
    }
    unless ($seen{$infileSampleID} > 1) {

	$logger->fatal($sampleID." supplied and sampleID ".$infileSampleID." found in file : ".${$infileHashRef}{$sampleID}[$infileCounter]." does not match. Please rename file to match sampleID: ".$sampleID."\n");
	exit;
    }
}


sub AddInfileInfoOld {
  
##AddInfileInfoOld
    
##Function : Adds information derived from infile name to sampleInfo hash. Tracks the number of lanes sequenced and checks unique array elementents.
##Returns  : ""
##Arguments: $scriptParameterHashRef, $laneHashRef, $infilesLaneNoEndingHashRef, $sampleInfoHashRef, $infilesBothStrandsNoEndingHashRef, $infileHashRef, $inDirPathHashRef, $dateFlowCell, $lane, $direction, $sampleID, $laneTrackerRef, $infileCounter, $compressedInfo
##         : $scriptParameterHashRef            => The active parameters for this analysis hash {REF}
##         : $laneHashRef                       => The lane info hash {REF}
##         : $infilesLaneNoEndingHashRef        => The infile(s) without the ".ending" {REF}
##         : $sampleInfoHashRef                 => Info on samples and family hash {REF}
##         : $infilesBothStrandsNoEndingHashRef => The infile(s) without the ".ending" and strand info {REF}
##         : $infileHashRef                     => The infiles hash {REF}
##         : $inDirPathHashRef                  => The indirectories path(s) hash {REF}
##         : $dateFlowCell                      => Date and flow-cell
##         : $lane                              => Flow-cell lane
##         : $direction                         => Sequencing read direction
##         : $sampleID                          => Sample id
##         : $laneTrackerRef                    => Counts the number of lanes sequenced {REF}
##         : $infileCounter                     => Counts the number of infiles
##         : $compressedInfo                    => ".fastq.gz" or ".fastq" info governs zcat or cat downstream

    my $scriptParameterHashRef = $_[0];
    my $laneHashRef = $_[1];
    my $infilesLaneNoEndingHashRef = $_[2];
    my $sampleInfoHashRef = $_[3];
    my $infilesBothStrandsNoEndingHashRef = $_[4];
    my $infileHashRef = $_[5];
    my $inDirPathHashRef = $_[6];
    my $dateFlowCell = $_[7];
    my $lane = $_[8];
    my $direction = $_[9];
    my $sampleID = $_[10];
    my $laneTrackerRef = $_[11];
    my $infileCounter = $_[12];
    my $compressedInfo = $_[13]; #".fastq.gz" or ".fastq" info governs zcat or cat downstream

    my $readFile;

    if ($compressedInfo eq "compressed") {

	$readFile = "zcat";  #Read file in compressed format
    }
    else {

	$readFile = "cat";  #Read file in uncompressed format
    }
    
    my $seqLengthRegExp = q?perl -ne 'if ($_!~/@/) {chomp($_);my $seqLength = length($_);print $seqLength;last;}' ?;  #Prints sequence length and exits

    if ($direction == 1) {  #Read 1
	
	push( @{${$laneHashRef}{$sampleID}}, $lane);  #Lane
	${$infilesLaneNoEndingHashRef}{$sampleID}[$$laneTrackerRef]= $dateFlowCell.".".$lane;  #Save old format in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry into hash for both strands and .ending is removed (.fastq).

	${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{ ${$infilesLaneNoEndingHashRef}{ $sampleID }[$$laneTrackerRef] }{'sequenceRunType'} = "Single-end";  #Single-end until proven otherwise
	$$laneTrackerRef++;
    }
    if ($direction == 2) {  #2nd read direction

	${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{ ${$infilesLaneNoEndingHashRef}{ $sampleID }[$$laneTrackerRef-1] }{'sequenceRunType'} = "Paired-end";  #$laneTracker -1 since it gets incremented after direction eq 1. 
    }
    
    ${$infilesBothStrandsNoEndingHashRef}{ $sampleID }[$infileCounter]= $dateFlowCell.".".$lane."_".$direction;  #Save new format in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry per strand and .ending is removed (.fastq).

    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{ ${$infilesBothStrandsNoEndingHashRef}{ $sampleID }[$infileCounter] }{'originalFileName'} = ${$infileHashRef}{$sampleID}[$infileCounter];  #Original fileName

    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{ ${$infilesBothStrandsNoEndingHashRef}{ $sampleID }[$infileCounter] }{'originalFileNameNoEnding'} = $dateFlowCell."lane".$lane."_".$direction;  #Original fileName, but no ending

    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{ ${$infilesBothStrandsNoEndingHashRef}{ $sampleID }[$infileCounter] }{'sampleBarcode'} = "X";  #Save barcode, but not defined

    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{ ${$infilesBothStrandsNoEndingHashRef}{ $sampleID }[$infileCounter] }{'runBarcode'} = $lane."_".$dateFlowCell;  #Save run barcode

    &CheckUniqueArrayElement(\@{ ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{ ${$infilesBothStrandsNoEndingHashRef}{ $sampleID }[$infileCounter] }{'readDirection'} }, \$direction);  #Check if there are any new info and add it to sampleInfo if so. 

    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{ ${$infilesBothStrandsNoEndingHashRef}{ $sampleID }[$infileCounter] }{'sequenceLength'} = `cd ${$inDirPathHashRef}{$sampleID};$readFile ${$infileHashRef}{$sampleID}[$infileCounter] | $seqLengthRegExp;`;  #Collect sequence length
}


sub AddInfileInfo {
  
##AddInfileInfo
    
##Function : Adds information derived from infile name to sampleInfo hash. Tracks the number of lanes sequenced and checks unique array elementents.
##Returns  : ""
##Arguments: $scriptParameterHashRef, $laneHashRef, $infilesLaneNoEndingHashRef, $sampleInfoHashRef, $infilesBothStrandsNoEndingHashRef, $infileHashRef, $inDirPathHashRef, $lane, $date, $flowCell, $sampleID, $index, $direction, $laneTrackerRef, $infileCounter, $compressedInfo
##         : $scriptParameterHashRef            => The active parameters for this analysis hash {REF}
##         : $laneHashRef                       => The lane info hash {REF}
##         : $infilesLaneNoEndingHashRef        => The infile(s) without the ".ending" {REF}
##         : $sampleInfoHashRef                 => Info on samples and family hash {REF}
##         : $infilesBothStrandsNoEndingHashRef => The infile(s) without the ".ending" and strand info {REF}
##         : $infileHashRef                     => The infiles hash {REF}
##         : $inDirPathHashRef                  => The indirectories path(s) hash {REF}
##         : $lane                              => Flow-cell lane
##         : $date                              => Flow-cell sequencing date
##         : $flowCell                          => Flow-cell id
##         : $sampleID                          => Sample id
##         : $index                             => The DNA library preparation molecular barcode
##         : $direction                         => Sequencing read direction
##         : $laneTrackerRef                    => Counts the number of lanes sequenced {REF}
##         : $infileCounter                     => Counts the number of infiles
##         : $compressedInfo                    => ".fastq.gz" or ".fastq" info governs zcat or cat downstream

    my $scriptParameterHashRef = $_[0];
    my $laneHashRef = $_[1];
    my $infilesLaneNoEndingHashRef = $_[2];
    my $sampleInfoHashRef = $_[3];
    my $infilesBothStrandsNoEndingHashRef = $_[4];
    my $infileHashRef = $_[5];
    my $inDirPathHashRef = $_[6];
    my $lane = $_[7];
    my $date = $_[8];
    my $flowCell = $_[9];
    my $sampleID = $_[10];
    my $index = $_[11];
    my $direction = $_[12];
    my $laneTrackerRef = $_[13];
    my $infileCounter = $_[14];
    my $compressedInfo = $_[15];

    my $readFile;

    if ($compressedInfo eq "compressed") {

	$readFile = "zcat";  #Read file in compressed format
    }
    else {

	$readFile = "cat";  #Read file in uncompressed format
    }
    
    my $seqLengthRegExp = q?perl -ne 'if ($_!~/@/) {chomp($_);my $seqLength = length($_);print $seqLength;last;}' ?;  #Prints sequence length and exits

    if ($direction == 1) {  #Read 1

	push( @{${$laneHashRef}{$sampleID}}, $lane);  #Lane
	${$infilesLaneNoEndingHashRef}{$sampleID}[$$laneTrackerRef]= $sampleID.".".$date."_".$flowCell."_".$index.".lane".$1;  #Save new format (sampleID_date_flow-cell_index_lane) in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry into hash for both strands and .ending is removed (.fastq).

	${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{ ${$infilesLaneNoEndingHashRef}{ $sampleID }[$$laneTrackerRef] }{'sequenceRunType'} = "Single-end";  #Single-end until proven otherwise
	$$laneTrackerRef++;
    }
    if ($direction == 2) {  #2nd read direction

	${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{ ${$infilesLaneNoEndingHashRef}{ $sampleID }[$$laneTrackerRef-1] }{'sequenceRunType'} = "Paired-end";  #$laneTracker -1 since it gets incremented after direction eq 1. 
    }
    
    ${$infilesBothStrandsNoEndingHashRef}{ $sampleID }[$infileCounter]= $sampleID.".".$date."_".$flowCell."_".$index.".lane".$1."_".$direction;  #Save new format in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry per strand and .ending is removed (.fastq).

    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{ ${$infilesBothStrandsNoEndingHashRef}{ $sampleID }[$infileCounter] }{'originalFileName'} = ${$infileHashRef}{$sampleID}[$infileCounter];  #Original fileName

    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{ ${$infilesBothStrandsNoEndingHashRef}{ $sampleID }[$infileCounter] }{'originalFileNameNoEnding'} = $1."_".$date."_".$flowCell."_".$sampleID."_".$index."_".$direction;  #Original fileName, but no ending

    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{ ${$infilesBothStrandsNoEndingHashRef}{ $sampleID }[$infileCounter] }{'lane'} = $1;  #Save sample lane                  

    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{ ${$infilesBothStrandsNoEndingHashRef}{ $sampleID }[$infileCounter] }{'date'} = $date;  #Save Sequence run date          

    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{ ${$infilesBothStrandsNoEndingHashRef}{ $sampleID }[$infileCounter] }{'flow-cell'} = $flowCell;  #Save Sequence flow-cell        

    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{ ${$infilesBothStrandsNoEndingHashRef}{ $sampleID }[$infileCounter] }{'sampleBarcode'} = $index;  #Save sample barcode

    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{ ${$infilesBothStrandsNoEndingHashRef}{ $sampleID }[$infileCounter] }{'runBarcode'} = $date."_".$flowCell."_".$1."_".$index;  #Save run barcode
    
    &CheckUniqueArrayElement(\@{  ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{ ${$infilesBothStrandsNoEndingHashRef}{ $sampleID }[$infileCounter] }{'ReadDirection'} }, \$direction);  #Check if there are any new info and add it to sampleInfo if so.   

    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'file'}{ ${$infilesBothStrandsNoEndingHashRef}{ $sampleID }[$infileCounter] }{'sequenceLength'} = `cd ${$inDirPathHashRef}{$sampleID};$readFile ${$infileHashRef}{$sampleID}[$infileCounter] | $seqLengthRegExp;`;  #Collect sequence length
}


sub CheckFileNameExists {

##CheckFileNameExists
    
##Function : Check if a file with with a filename consisting of $filePathRef.$fileCounter.$fileEndingRef exist. If so bumps the version number and return new filename and sbatch version number.
##Returns  : "$fileName, $fileNameTracker"
##Arguments: $filePathRef, $fileEndingRef
##         : $filePathRef   => The file path {REF}
##         : $fileEndingRef => The file ending {REF}

    my $filePathRef = $_[0];
    my $fileEndingRef = $_[1];

    my $fileName;  #Temp filename
    my $fileNameTracker = 0;  #Nr of sbatch scripts with identical filenames i.e. version number
  
    for (my $fileCounter=0;$fileCounter<9999;$fileCounter++) {  #Number of possible files with the same name
	
	$fileName = $$filePathRef.$fileCounter.$$fileEndingRef;  #Filename, filenr and fileending
	$fileNameTracker = $fileCounter;  #Nr of sbatch scripts with identical filenames

	unless (-f $fileName) {  #File exists

	    last;  #No file exists
	}	
    }
    return ($fileName, $fileNameTracker);
}


sub DefineHashParameters {

##DefineHashParameters
    
##Function : Check if user supplied cmd info and supplies HashParameters to scriptParameters
##Returns : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $orderParametersArrayRef, $broadcastsArrayRef, $parameterName, $parameterType, $parameterDefault, $associatedPrograms, $parameterExistsCheck
##         : $parameterHashRef        => The parameters hash {REF}
##         : $scriptParameterHashRef  => The active parameters for this analysis hash
##         : $orderParametersArrayRef => Order of addition to parameter array {REF}
##         : $broadcastsArrayRef      => Holds the parameters info for broadcasting later
##         : $parameterName           => MIP parameter to evaluate
##         : $parameterType           => Type of MIP parameter
##         : $parameterDefault        => The parameter default value
##         : $associatedPrograms      => Programs that use the parameter. Comma separated string
##         : $parameterExistsCheck    => Check if intendent file exists in reference directory

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $orderParametersArrayRef = $_[2];
    my $broadcastsArrayRef = $_[3];
    my $parameterName = $_[4];
    my $parameterType = $_[5];
    my $parameterDefault = $_[6];
    my $associatedPrograms = $_[7];
    my $parameterExistsCheck = $_[8];

    ${$parameterHashRef}{$parameterName}{'hash'} = "yes";  #To separate scalars and arrays from hashes

    if (defined(${$parameterHashRef}{$parameterName}{'value'})) { #Input from cmd 

	${$scriptParameterHashRef}{$parameterName} = ${$parameterHashRef}{$parameterName}{'value'};  #Save cmd supplied info	
    }
    else {

	${$parameterHashRef}{$parameterName}{'value'} = "nocmdinput"; #To enable use of subroutine &AddToScriptParameter
    }
    push(@{$orderParametersArrayRef}, $parameterName); #Add to enable later evaluation of parameters in proper order & write to MIP log file

    &AddToScriptParameter(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \%sampleInfo, \%fileInfo, \%referenceFileEndings, \@{$broadcastsArrayRef},
			  {'parameterName' => $parameterName,
			   'parameterValue' => ${$parameterHashRef}{$parameterName}{'value'},
			   'parameterType' => $parameterType,
			   'parameterDefault' => $parameterDefault,
			   'associatedPrograms' => $associatedPrograms,
			   'parameterExistsCheck' => $parameterExistsCheck
			  });
}


sub DefineArrayParameters {

##DefineArrayParameters
    
##Function : Check if user supplied cmd info and supplies arrayParameters to scriptParameters
##Returns : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $arrayRef, $orderParametersArrayRef, $broadcastsArrayRef, $parameterName, $parameterType, $parameterDefault, $associatedPrograms, $parameterExistsCheck
##         : $parameterHashRef        => The parameters hash {REF}
##         : $scriptParameterHashRef  => The active parameters for this analysis hash
##         : $arrayRef                => Array to loop in for parameter {REF}
##         : $orderParametersArrayRef => Order of addition to parameter array {REF}
##         : $broadcastsArrayRef      => Holds the parameters info for broadcasting later
##         : $parameterName           => MIP parameter to evaluate
##         : $parameterType           => Type of MIP parameter
##         : $parameterDefault        => The parameter default value
##         : $associatedPrograms      => Programs that use the parameter. Comma separated string
##         : $parameterExistsCheck    => Check if intendent file exists in reference directory

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $arrayRef = $_[2];
    my $orderParametersArrayRef = $_[3];
    my $broadcastsArrayRef = $_[4];
    my $parameterName = $_[5];
    my $parameterType = $_[6];
    my $parameterDefault = $_[7];
    my $associatedPrograms = $_[8];
    my $parameterExistsCheck = $_[9];

    ${$parameterHashRef}{$parameterName}{'array'} = "yes";  #To separate scalars from arrays

    if (scalar(@{$arrayRef}) == 0) { #No input from cmd 

	${$parameterHashRef}{$parameterName}{'value'} = ["nocmdinput"]; #To enable use of subroutine &AddToScriptParameter
    }
    else {

	@{${$scriptParameterHashRef}{$parameterName}} = split(',', join(',', @{$arrayRef})); #If user supplied parameter a comma separated list
    }
    push(@{$orderParametersArrayRef}, $parameterName); #Add to enable later evaluation of parameters in proper order & write to MIP log file
    
    &AddToScriptParameter(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \%sampleInfo, \%fileInfo, \%referenceFileEndings, \@{$broadcastsArrayRef},
			  {'parameterName' => $parameterName,
			   'parameterValue' => ${$parameterHashRef}{$parameterName}{'value'}[0],
			   'parameterType' => $parameterType,
			   'parameterDefault' => $parameterDefault,
			   'associatedPrograms' => $associatedPrograms,
			   'parameterExistsCheck' => $parameterExistsCheck
			  });
}


sub DefineParametersPath {

##DefineParametersPath
    
##Function : Defines all attributes of a parameter, so that the correct value can be set and added to %scriptparameter later.
##Arguments: $parameterHashRef, $orderParametersArrayRef, $parameterName, $parameterDefault, $associatedProgram, $existsCheck
##         : $parameterHashRef        => Parameter Hash {REF}
##         : $orderParametersArrayRef => Order of addition to parameter array {REF}
##         : $parameterName     => Parameter name
##         : $parameterDefault  => Default setting
##         : $associatedProgram => The parameters program
##         : $existsCheck       => Check if intendent file exists in reference directory
##         : $buildFile         => Autovivication of file if it does not exists (yes or no)

    my $parameterHashRef = $_[0];
    my $orderParametersArrayRef = $_[1];
    my $parameterName = $_[2];
    my $parameterDefault = $_[3];
    my $associatedProgram = $_[4];
    my $existsCheck = $_[5];
    my $buildFile = $_[6];
    
    ${$parameterHashRef}{$parameterName} = {
	'type' => "path",
	'value' => "nocmdinput",
	'default' => $parameterDefault,
	'associatedProgram' => $associatedProgram,
	'existsCheck' => $existsCheck,
	'buildFile' => $buildFile,
    };
    
    push(@{$orderParametersArrayRef}, $parameterName);  #Add to enable later evaluation of parameters in proper order & write to master file
}


sub DefineParameters {

##DefineParameters
    
##Function : Defines all attributes of a parameter, so that the correct value can be set and added to %scriptparameter later.
##Arguments: $parameterHashRef, $orderParametersArrayRef, $parameterName, $parameterType, $parameterDefault, $associatedProgram, $fileEnding, @programNamePath
##         : $parameterHashRef        => Parameter Hash {REF}
##         : $orderParametersArrayRef => Order of addition to parameter array {REF}
##         : $parameterName           => Parameter name
##         : $parameterType           => MIP or program
##         : $parameterDefault        => Default setting
##         : $associatedProgram       => The parameters program
##         : $fileEnding              => The filending after the module has been run
##         : $parameterChain          => The chain to which the program belongs to
##         : @programNamePath         => The path name of the program(s) for each sbatch script

    my $parameterHashRef = $_[0];
    my $orderParametersArrayRef = $_[1];
    my $parameterName = $_[2];
    my $parameterType = $_[3];
    my $parameterDefault = $_[4];
    my $associatedProgram = $_[5];
    my $fileEnding = $_[6];
    my $parameterChain = $_[7];
    my @programNamePath;

    if (defined($_[8])) {

	@programNamePath = split(":", $_[8]);
    }
    if (defined($programNamePath[0])) {
	
	${$parameterHashRef}{$parameterName} = {
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
	
	${$parameterHashRef}{$parameterName} = {
	    'type' => $parameterType,
	    'value' => "nocmdinput",
	    'default' => $parameterDefault,
	    'associatedProgram' => $associatedProgram,
	    'fileEnding' => $fileEnding,
	    'chain' => $parameterChain,
	};
    }

    push(@{$orderParametersArrayRef}, $parameterName);  #Add to enable later evaluation of parameters in proper order & write to MIP log file
}


sub AddToScriptParameter {

##AddToScriptParameter
    
##Function : Checks and sets user input or default values to scriptParameters
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleInfoHashRef, $fileInfoRef, $referenceFileEndingsHashRef, $broadcastsArrayRef, 
##         : $parameterHashRef                   => Holds all parameters
##         : $scriptParameterHashRef             => Holds all set parameter for analysis
##         : $sampleInfoHashRef                  => Info on samples and family hash {REF}
##         : $fileInfoHashRef                    => The fileInfo hash {REF}
##         : $referenceFileEndingsHashRef        => The associated reference file endings {REF}
##         : $broadcastsArrayRef                 => Holds the parameters info for broadcasting later
##         : $argHashRef{'parameterName'}        => Parameter name
##         : $argHashRef{'parameterValue'}       => Parameter value to evaluate
##         : $argHashRef{'parameterType'}        => Path, MIP or program
##         : $argHashRef{'parameterDefault'}     => Default setting
##         : $argHashRef{'associatedPrograms'}   => The parameters program(s) {array, REF}
##         : $argHashRef{'parameterExistsCheck'} => Check if intendent file exists in reference directory
##         : $argHashRef{'programNamePath'}      => Program name in system path
    
    my $parameterHashRef = shift;
    my $scriptParameterHashRef = shift;
    my $sampleInfoHashRef = shift;
    my $fileInfoHashRef = shift;
    my $referenceFileEndingsHashRef = shift;
    my $broadcastsArrayRef = shift;
    my ($argHashRef) = @_;

    my @associatedPrograms;

    if (defined(${$argHashRef}{'associatedPrograms'})) {  #Restore as array
	
	@associatedPrograms = split(/,/, ${$argHashRef}{'associatedPrograms'});
    }

    foreach my $associatedProgram (@associatedPrograms) {  #Check all programs that use parameter

	my $parameterSetSwitch = 0;
	
	if (defined(${$scriptParameterHashRef}{$associatedProgram}) && (${$scriptParameterHashRef}{$associatedProgram} > 0) ) {  #Only add active programs parameters	    
	    
	    $parameterSetSwitch = 1;

	    if (${$argHashRef}{'parameterType'} eq "path") {  #Evaluate "Path" parameters
		
		if (${$argHashRef}{'parameterValue'} eq "nocmdinput") {  #No input from cmd
		    
		    if (defined(${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}})) {  #Input from config file
			
			if (${$argHashRef}{'parameterName'} eq "exomeTargetBedInfileLists") {  #ExomeTargetBedInfileListss is a comma separated list 
			    
			    &SetTargetandAutoBuild(\%parameter, \%{$scriptParameterHashRef}, \%{$fileInfoHashRef}, \@{${$scriptParameterHashRef}{'sampleIDs'}}, \${$argHashRef}{'parameterName'}, \${$referenceFileEndingsHashRef}{'exomeTargetBedInfileLists'});
			}
			if (${$argHashRef}{'parameterName'} eq "exomeTargetPaddedBedInfileLists") {  #ExomeTargetPaddedBedInfileLists is a comma separated list 
			    
			    &SetTargetandAutoBuild(\%parameter, \%{$scriptParameterHashRef}, \%{$fileInfoHashRef}, \@{${$scriptParameterHashRef}{'sampleIDs'}}, \${$argHashRef}{'parameterName'}, \${$referenceFileEndingsHashRef}{'exomeTargetPaddedBedInfileLists'});
			}
			if ( (${$argHashRef}{'parameterName'} eq "GATKTargetPaddedBedIntervalLists") && (${$scriptParameterHashRef}{'analysisType'} ne "genomes") ) {  #GATKTargetPaddedBedIntervalLists is a comma separated list 
			    
			    &SetTargetandAutoBuild(\%parameter, \%{$scriptParameterHashRef}, \%{$fileInfoHashRef}, \@{${$scriptParameterHashRef}{'sampleIDs'}}, \${$argHashRef}{'parameterName'}, \${$referenceFileEndingsHashRef}{'GATKTargetPaddedBedIntervalLists'});
			}
			if (${$argHashRef}{'parameterName'} eq "snpSiftAnnotationFiles") {
			    
			    &BreakString(\%fileInfo, \${$scriptParameterHashRef}{'snpSiftAnnotationFiles'}, \${$argHashRef}{'parameterName'}, \$associatedProgram);
			}
			if (${$argHashRef}{'parameterName'} eq "humanGenomeReference") {

			    &ParseHumanGenomeReference(\%{$fileInfoHashRef}, \${$scriptParameterHashRef}{'humanGenomeReference'});
			}
			if (${$argHashRef}{'parameterName'} eq "pedigreeFile") {
			    
			    &ReadPlinkPedigreeFile(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \%sampleInfo, \%referenceFileEndings, \%supportedCaptureKit, ${$scriptParameterHashRef}{'pedigreeFile'});
			}
		    }
		    elsif (${$argHashRef}{'parameterDefault'} ne "nodefault") {  #Add default value
			
			if (${$argHashRef}{'parameterName'} eq "inFilesDirs") {
			    
			    for (my $indirectoryCount=0;$indirectoryCount<scalar(@{${$scriptParameterHashRef}{'sampleIDs'}});$indirectoryCount++) {
				
				push(@{${$scriptParameterHashRef}{'inFilesDirs'}}, ${$scriptParameterHashRef}{'clusterConstantPath'}."/".${$scriptParameterHashRef}{'analysisType'}."/".${$scriptParameterHashRef}{'sampleIDs'}[$indirectoryCount]."/fastq");					
			    }
			}
			elsif (${$argHashRef}{'parameterName'} eq "exomeTargetBedInfileLists") {
			    
			    &SetTargetandAutoBuild(\%parameter, \%{$scriptParameterHashRef}, \%{$fileInfoHashRef}, \@{${$scriptParameterHashRef}{'sampleIDs'}}, \${$argHashRef}{'parameterName'}, \${$referenceFileEndingsHashRef}{'exomeTargetBedInfileLists'});
			}
			elsif (${$argHashRef}{'parameterName'} eq "exomeTargetPaddedBedInfileLists") {  #Note that potential pedigree files entries will be updated with GenomeReferenceSource and version here
			    
			    &SetTargetandAutoBuild(\%parameter, \%{$scriptParameterHashRef}, \%{$fileInfoHashRef}, \@{${$scriptParameterHashRef}{'sampleIDs'}}, \${$argHashRef}{'parameterName'}, \${$referenceFileEndingsHashRef}{'exomeTargetPaddedBedInfileLists'});
			}
			elsif ( (${$argHashRef}{'parameterName'} eq "GATKTargetPaddedBedIntervalLists") && (${$scriptParameterHashRef}{'analysisType'} ne "genomes") ) {  #Note that potential pedigree files entries will be updated with GenomeReferenceSource and version here 
			    
			    &SetTargetandAutoBuild(\%parameter, \%{$scriptParameterHashRef}, \%{$fileInfoHashRef}, \@{${$scriptParameterHashRef}{'sampleIDs'}}, \${$argHashRef}{'parameterName'}, \${$referenceFileEndingsHashRef}{'GATKTargetPaddedBedIntervalLists'});
			}
			elsif (${$argHashRef}{'parameterName'} eq "vepFeatures") {
			    
			    @{${$scriptParameterHashRef}{'vepFeatures'}} = ("refseq", "hgvs", "symbol", "numbers", "sift", "polyphen", "humdiv");  #Set default vep features
			}
			elsif (${$argHashRef}{'parameterName'} eq "snpSiftAnnotationFiles") {
			    
			    ${$scriptParameterHashRef}{'snpSiftAnnotationFiles'} = "dbsnp_138.b37.excluding_sites_after_129.vcf.gz:CAF,ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz:AF,ExAC.r0.1.sites.vep.vcf:AF";  #Set default snpSiftAnnotationFiles
			    &BreakString(\%fileInfo, \${$scriptParameterHashRef}{'snpSiftAnnotationFiles'}, \${$argHashRef}{'parameterName'}, \$associatedProgram);
			}
			elsif (${$argHashRef}{'parameterName'} eq "snpSiftDbNSFPAnnotations") {
    
			    @{${$scriptParameterHashRef}{'snpSiftDbNSFPAnnotations'}} = ("SIFT_pred", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_pred", "LRT_pred", "MutationTaster_pred", "GERP++_NR", "GERP++_RS", "phastCons100way_vertebrate", "1000Gp1_AF", "ESP6500_AA_AF");  #Set default snpSiftDbNSFPAnnotations
			}
			elsif (${$argHashRef}{'parameterName'} eq "annovarTableNames") {
			    
			    @{${$scriptParameterHashRef}{'annovarTableNames'}} = ("refGene", "mce46way", "gerp++elem", "segdup", "tfbs", "mirna", "snp137NonFlagged", "1000g2012apr_all", "esp6500si_all", "ljb2_sift", "ljb2_pp2hdiv", "ljb2_pp2hvar", "ljb2_mt", "ljb2_lrt", "ljb2_gerp++", "ljb2_phylop");  #Set default annovar table names
			}
			else {
			    
			    if (${$argHashRef}{'parameterName'} eq "humanGenomeReference") {

				&ParseHumanGenomeReference(\%{$fileInfoHashRef}, \${$scriptParameterHashRef}{'humanGenomeReference'});
			    }			    
			    ${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}} = ${$argHashRef}{'parameterDefault'};  #Set default value
			}
		    }
		    else {  #No default

			if (${$argHashRef}{'parameterName'} eq "mosaikAlignReference") {  #Special case - do nothing, since file can be created by MIP from the humanGenomeReference if required
			}
			elsif ( (${$argHashRef}{'parameterName'} eq "bwaMemRapidDb") && (${$scriptParameterHashRef}{'analysisType'} ne "rapid")) {  #Do nothing since file is not required unless rapid mode is enabled
			}
			elsif ( (${$argHashRef}{'parameterName'} eq "GATKGenoTypeGVCFsRefGVCF") && (${$scriptParameterHashRef}{'analysisType'} =~/genomes/) ) {  #Do nothing since file is not required unless exome or rapid mode is enabled
			}
			elsif ( (${$argHashRef}{'parameterName'} eq "vcfParserRangeFeatureAnnotationColumns") && ( ${$scriptParameterHashRef}{'vcfParserRangeFeatureFile'} eq "noUserInfo") ) {  #Do nothing since no SelectFile was given
			} 
			elsif ( (${$argHashRef}{'parameterName'} eq "vcfParserSelectFeatureAnnotationColumns") && ( ${$scriptParameterHashRef}{'vcfParserSelectFile'} eq "noUserInfo") ) {  #Do nothing since no SelectFile was given
			}
			elsif ( (${$argHashRef}{'parameterName'} eq "vcfParserSelectFileMatchingColumn") && ( ${$scriptParameterHashRef}{'vcfParserSelectFile'} eq "noUserInfo") ) {  #Do nothing since no SelectFile was given
			}
			elsif ( (${$argHashRef}{'parameterName'} eq "geneFile") && (${$scriptParameterHashRef}{'pVariantEffectPredictor'} > 0) ) {  #Do nothing since VEP annotations can be used
			}
			elsif ( (${$argHashRef}{'parameterName'} eq "caddWGSSNVsFile") && ( ${$scriptParameterHashRef}{'caddWGSSNVs'} == 0) ) {  #Do nothing since no CADD annotation should be performed
			}
			elsif ( (${$argHashRef}{'parameterName'} eq "cadd1000GenomesFile") && ( ${$scriptParameterHashRef}{'cadd1000Genomes'} == 0) ) {  #Do nothing since no CADD annotation should be performed
			}
			elsif ( (${$argHashRef}{'parameterName'} eq "rankModelFile") && ( ${$scriptParameterHashRef}{'rankModelFile'} eq "noUserInfo") ) {  #Do nothing since no rank model was given i.e. use rank scripts deafult supplied with distribution
			}
			else {
			    
			    $logger->fatal($USAGE, "\n");
			    $logger->fatal("Supply '-".${$argHashRef}{'parameterName'}."' if you want to run ".$associatedProgram, "\n");
			    exit;
			}
		    }
		}
		else {  #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
		    
		    if (${$argHashRef}{'parameterName'} eq "exomeTargetBedInfileLists") {	    
			
			&EnableArrayParameter(\%{$scriptParameterHashRef}, \@exomeTargetBedInfileLists, \${$argHashRef}{'parameterName'});
			&CompareArrayElements(\@{${$scriptParameterHashRef}{'sampleIDs'}}, \@exomeTargetBedInfileLists, "sampleIDs", ${$argHashRef}{'parameterName'});
			&SetAutoBuildAndScriptParameterPerSample(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \@{${$scriptParameterHashRef}{'sampleIDs'}}, \@exomeTargetBedInfileLists, \${$argHashRef}{'parameterName'});
		    }
		    elsif (${$argHashRef}{'parameterName'} eq "exomeTargetPaddedBedInfileLists") {	    
			
			&EnableArrayParameter(\%{$scriptParameterHashRef}, \@exomeTargetPaddedBedInfileLists, \${$argHashRef}{'parameterName'});
			&CompareArrayElements(\@{${$scriptParameterHashRef}{'sampleIDs'}}, \@exomeTargetPaddedBedInfileLists, "sampleIDs", ${$argHashRef}{'parameterName'});
			&SetAutoBuildAndScriptParameterPerSample(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \@{${$scriptParameterHashRef}{'sampleIDs'}}, \@exomeTargetPaddedBedInfileLists, \${$argHashRef}{'parameterName'});
		    }
		    elsif ( (${$argHashRef}{'parameterName'} eq "GATKTargetPaddedBedIntervalLists") && (${$scriptParameterHashRef}{'analysisType'} ne "genomes") ) {	    
			
			&EnableArrayParameter(\%{$scriptParameterHashRef}, \@GATKTargetPaddedBedIntervalLists, \${$argHashRef}{'parameterName'});
			&CompareArrayElements(\@{${$scriptParameterHashRef}{'sampleIDs'}}, \@GATKTargetPaddedBedIntervalLists, "sampleIDs", ${$argHashRef}{'parameterName'});
			&SetAutoBuildAndScriptParameterPerSample(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \@{${$scriptParameterHashRef}{'sampleIDs'}}, \@GATKTargetPaddedBedIntervalLists, \${$argHashRef}{'parameterName'});
		    }
		    elsif (${$argHashRef}{'parameterName'} eq "pedigreeFile") {  #Must come after arrays that can be populated from pedigree file to not overwrite user cmd input 

			${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}} = ${$argHashRef}{'parameterValue'};
			&ReadPlinkPedigreeFile(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \%sampleInfo, \%referenceFileEndings, \%supportedCaptureKit, ${$scriptParameterHashRef}{'pedigreeFile'});
		    }
		    elsif (${$argHashRef}{'parameterName'} eq "snpSiftAnnotationFiles") {

			&BreakString(\%fileInfo, \${$parameterHashRef}{${$argHashRef}{'parameterName'}}{'value'}, \${$argHashRef}{'parameterName'}, \$associatedProgram);
		    }
		    else {
			
			if (${$argHashRef}{'parameterName'} eq "humanGenomeReference") {
			    
			    &ParseHumanGenomeReference(\%{$fileInfoHashRef}, \${$scriptParameterHashRef}{'humanGenomeReference'});
			}
			if (defined(${$parameterHashRef}{${$argHashRef}{'parameterName'}}{'array'})) {  #Do nothing
			}
			else {
			    
			    ${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}} = ${$argHashRef}{'parameterValue'};
			}
		    }
		}
		if ( ${$argHashRef}{'parameterExistsCheck'} && (${$argHashRef}{'parameterExistsCheck'} eq "directory") ) {  #Check dir existence
		    
		    if (${$argHashRef}{'parameterName'} eq "inFilesDirs") {
			
			for (my $indirectoryCount=0;$indirectoryCount<scalar(@{${$scriptParameterHashRef}{'inFilesDirs'}});$indirectoryCount++) {
			    
			    &CheckExistance(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \${$scriptParameterHashRef}{'inFilesDirs'}[$indirectoryCount], \${$argHashRef}{'parameterName'}, "d");
			}
		    }
		    else {
			
			&CheckExistance(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}}, \${$argHashRef}{'parameterName'}, "d");

			if (${$argHashRef}{'parameterName'} eq "genomeAnalysisToolKitPath") {  #To enable addition of version to sampleInfo
			    
			    if (${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}}=~/GenomeAnalysisTK-([^,]+)/) {
				
				${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'program'}{"GATK"}{'Version'} = $1;
			    }
			}
			if (${$argHashRef}{'parameterName'} eq "picardToolsPath") {  #To enable addition of version to sampleInfo

                            if (${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}}=~/picard-tools-([^,]+)/) {
                                
				${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'program'}{"PicardTools"}{'Version'} = $1;
                            }
                        }
		    }
		}
		elsif ( (${$argHashRef}{'parameterExistsCheck'}) && (${$argHashRef}{'parameterExistsCheck'} eq "file") && (defined(${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}})) ) {  #Check file existence in reference directory
		    
		    if (${$argHashRef}{'parameterName'} eq "mosaikJumpDbStub") {
			
			&CheckFileEndingsToBeBuilt(\%{$scriptParameterHashRef}, \@mosaikJumpDbStubFileEndings, ${$argHashRef}{'parameterName'}); 
		    }
		    elsif (${$argHashRef}{'parameterName'} eq "bwaBuildReference") {
			
			&CheckFileEndingsToBeBuilt(\%{$scriptParameterHashRef}, \@bwaBuildReferenceFileEndings, ${$argHashRef}{'parameterName'});
		    }
		    elsif (${$argHashRef}{'parameterName'} eq "picardToolsMergeSamFilesPrevious") {

			for (my $fileCounter=0;$fileCounter<scalar(@{${$scriptParameterHashRef}{'picardToolsMergeSamFilesPrevious'}});$fileCounter++) {  #All AnnovarTables

			    &CheckExistance(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \${$scriptParameterHashRef}{'picardToolsMergeSamFilesPrevious'}[$fileCounter], \${$argHashRef}{'parameterName'}, "f");
			}
		    }
		    elsif (${$argHashRef}{'parameterName'} eq "humanGenomeReference") {
			
			&CheckExistance(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \(${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}}), \${$argHashRef}{'parameterName'}, "f");#Check reference genome
			${$sampleInfoHashRef}{${$scriptParameterHashRef}{'familyID'}}{${$scriptParameterHashRef}{'familyID'}}{"HumanGenomeBuild"}{'Path'} = ${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}};
			${$sampleInfoHashRef}{${$scriptParameterHashRef}{'familyID'}}{${$scriptParameterHashRef}{'familyID'}}{"HumanGenomeBuild"}{'Source'} = ${$fileInfoHashRef}{'humanGenomeReferenceSource'};
			${$sampleInfoHashRef}{${$scriptParameterHashRef}{'familyID'}}{${$scriptParameterHashRef}{'familyID'}}{"HumanGenomeBuild"}{'Version'} = ${$fileInfoHashRef}{'humanGenomeReferenceVersion'};

			&CheckHumanGenomeFileEndings(\%parameter, \${$scriptParameterHashRef}{'referencesDir'}, \@humanGenomeReferenceFileEndings, \${$fileInfoHashRef}{'humanGenomeReferenceNameNoEnding'}, \${$argHashRef}{'parameterName'});
		    }
		    elsif ( (${$argHashRef}{'parameterName'} eq "exomeTargetBedInfileLists") || (${$argHashRef}{'parameterName'} eq "exomeTargetPaddedBedInfileLists") || (${$argHashRef}{'parameterName'} eq "GATKTargetPaddedBedIntervalLists") ) {
			
			if ( (${$argHashRef}{'parameterName'} eq "GATKTargetPaddedBedIntervalLists") && (${$scriptParameterHashRef}{'analysisType'} eq "genomes") ) {  #No need to check since genomes does not use GATKTargetPaddedBedIntervalLists
			}
			else {
			    
			    for (my $sampleIDsCounter=0;$sampleIDsCounter<scalar(@{${$scriptParameterHashRef}{'sampleIDs'}});$sampleIDsCounter++) {  #All sampleIDs
				
				unless (defined(${$scriptParameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDsCounter] }{${$argHashRef}{'parameterName'}})) {  #No capture kit supplied

				    my $captureKit = &AddCaptureKit(\%{$referenceFileEndingsHashRef}, \%supportedCaptureKit, 
								     {'captureKit' => "Latest", 
								      'parameterName' => ${$argHashRef}{'parameterName'},
								     });
				    $logger->warn("Could not detect a supplied capture kit. Will Try to use 'Latest' capture kit: ".$captureKit, "\n");
				    ${$scriptParameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDsCounter] }{${$argHashRef}{'parameterName'}} = $captureKit;
				}
				 &CheckSupportedFileEnding(\(${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDsCounter] }{ ${$argHashRef}{'parameterName'} }), \${$referenceFileEndingsHashRef}{${$argHashRef}{'parameterName'}}, \${$argHashRef}{'parameterName'});
				 &CheckExistance(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \(${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDsCounter]}{${$argHashRef}{'parameterName'}}), \${$argHashRef}{'parameterName'}, "f", \${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDsCounter]);
			    
				my $exomeTargetBedFileNoEnding = &RemoveFileEnding(\${$scriptParameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDsCounter]}{${$argHashRef}{'parameterName'}} , ${$referenceFileEndingsHashRef}{${$argHashRef}{'parameterName'}});  #Remove ".fileending" from reference filename
				 &CheckTargetExistFileBed(\%{$scriptParameterHashRef}, \$exomeTargetBedFileNoEnding, ${$argHashRef}{'parameterName'});
			    }
			    undef(${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}});  #Remove parameter to avoid unnecessary print to STDOUT and config
			}
		    }
                    elsif (${$argHashRef}{'parameterName'} eq "snpSiftAnnotationFiles"){
			
			my %snpEffFile = &DefineSnpEffFiles(\%{$parameterHashRef});
			my $intendedFilePathRef;
			
			for my $file (keys %{${$fileInfoHashRef}{$associatedProgram}{'snpSiftAnnotationFiles'}}) {
			    
			    my $intendedFilePathRef = \(${$scriptParameterHashRef}{'referencesDir'}."/".$file);
			    &CheckExistance(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, $intendedFilePathRef, \$file, "f");

			    if ($file =~/\.gz$/) {  #Check for tabix index as well

				my $fileIndex = $file.".tbi";
				my $intendedFilePathRef = \(${$scriptParameterHashRef}{'referencesDir'}."/".$fileIndex);
				&CheckExistance(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, $intendedFilePathRef, \$fileIndex, "f");
			    }
			}
		    }
		    elsif (${$argHashRef}{'parameterName'} eq "annovarTableNames") {
			
			&CheckAnnovarTables(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \%annovarTable);
		    }
		    elsif (${$argHashRef}{'parameterName'} eq "configFile") {  #Do nothing since file existence is checked by &LoadYAML
		    }
		    elsif (${$argHashRef}{'parameterName'} eq "pedigreeFile") {  #Do nothing since file existence is checked by ReadPlinkPedigreeFile
		    }
		    elsif (${$argHashRef}{'parameterName'} eq "sampleInfoFile") {

			if (defined(${$scriptParameterHashRef}{'sampleInfoFile'})) {

			    if (-f ${$scriptParameterHashRef}{'sampleInfoFile'}) {

				%sampleInfo = &LoadYAML(\%scriptParameter, ${$scriptParameterHashRef}{'sampleInfoFile'});  #Load parameters from previous run from sampleInfoFile	  
				
			    }
			    if (defined(${$scriptParameterHashRef}{'pedigreeFile'}) ) {
				
				${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'pedigreeFile'}{'Path'} = ${$scriptParameterHashRef}{'pedigreeFile'};  #Add pedigreeFile to sampleInfo
				${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'pedigreeFileAnalysis'}{'Path'} = ${$scriptParameterHashRef}{'outDataDir'}."/".${$scriptParameterHashRef}{'familyID'}."/qc_pedigree.yaml";  #Add pedigreeFile info used in this analysis to SampleInfoFile
			    }
			} 
		    }
		    elsif (${$argHashRef}{'parameterName'} eq "logFile") {  #Do nothing since file is to be created
		    }
		    elsif ( (${$argHashRef}{'parameterName'} eq "bwaMemRapidDb") && (${$scriptParameterHashRef}{'analysisType'} ne "rapid")) {  #Do nothing since file is not required unless rapid mode is enabled
		    }
		    elsif ( (${$argHashRef}{'parameterName'} eq "GATKGenoTypeGVCFsRefGVCF") && (${$scriptParameterHashRef}{'analysisType'} =~/genomes/) ) {  #Do nothing since file is not required unless exome mode is enabled
		    }
                    elsif ( (${$argHashRef}{'parameterName'} eq "vcfParserRangeFeatureFile") && ( ${$scriptParameterHashRef}{'vcfParserRangeFeatureFile'} eq "noUserInfo") ) {  #Do nothing since no RangeFile was given
		    }
		    elsif (${$argHashRef}{'parameterName'} eq "vcfParserSelectFile") {
			
			if (${$scriptParameterHashRef}{'vcfParserSelectFile'} eq "noUserInfo") {  #Do nothing since no SelectFile was given
			    
			}
			else {  #To enable addition of selectFile to sampleInfo                                                                       
			    
			    &CheckExistance(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \(${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}}), \${$argHashRef}{'parameterName'}, "f");

			    if (${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}}=~/v(\d+\.\d+.\d+|\d+\.\d+)/) {
				
				${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'database'}{"SelectFile"}{'Version'} = $1;
			    }
			    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'database'}{"SelectFile"}{'File'} = ${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}};
			}
		    }
                    elsif ( (${$argHashRef}{'parameterName'} eq "geneFile") && (${$scriptParameterHashRef}{'pVariantEffectPredictor'} > 0) ) {  #Do nothing since VEP annotations can be used			    
		    }
		    elsif ( (${$argHashRef}{'parameterName'} eq "caddWGSSNVsFile") && ( ${$scriptParameterHashRef}{'caddWGSSNVs'} == 0) ) {  #Do nothing since no CADD annotation should be performed
		    }
		    elsif ( (${$argHashRef}{'parameterName'} eq "cadd1000GenomesFile") && ( ${$scriptParameterHashRef}{'cadd1000Genomes'} == 0) ) {  #Do nothing since no CADD annotation should be performed
		    }
		    elsif (${$argHashRef}{'parameterName'} eq "rankModelFile") {  
			
			if (${$scriptParameterHashRef}{'rankModelFile'} eq "noUserInfo") {  #Do nothing since no rank model config file was given. Usse default supplied by ranking script
			}
			else {  #To enable addition of rankModel file and version to sampleInfo                                                                       
			    
			    &CheckExistance(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \(${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}}), \${$argHashRef}{'parameterName'}, "f");
			    
			    if (${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}}=~/v(\d+\.\d+.\d+|\d+\.\d+)/) {
				
				${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'program'}{"RankVariants"}{'RankModel'}{'Version'} = $1;
			    }
			    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{'program'}{"RankVariants"}{'RankModel'}{'File'} = ${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}};
			}
		    }
		    else {
			
			 &CheckExistance(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \(${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}}), \${$argHashRef}{'parameterName'}, "f");
                    }		    
		}
	    }
	    
	    if (${$argHashRef}{'parameterType'} eq "MIP") {  #Evaluate "MIP" parameters
		
		if (${$argHashRef}{'parameterValue'} eq "nocmdinput") {  #No input from cmd
		    
		    if (defined(${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}})) {  #Input from config file - do nothing
		    }
		    elsif (${$argHashRef}{'parameterDefault'} ne "nodefault") {
		
			${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}} = ${$argHashRef}{'parameterDefault'};  #Set default value
		    }
		    else {
			
			if (${$argHashRef}{'parameterName'} eq "aligner") {  #Set to "nocmdinput"
			
			    ${$scriptParameterHashRef}{'aligner'} = "nocmdinput";
			}
			else {
			    
			    $logger->fatal($USAGE, "\n");
			    $logger->fatal("Supply '-".${$argHashRef}{'parameterName'}."' if you want to run ".$associatedProgram, "\n");
			    exit;
			}
		    }
		}
		else {  #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
		    
		    ${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}} = ${$argHashRef}{'parameterValue'}; 
		}
                if (${$argHashRef}{'parameterName'} eq "email") {

		    if (${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}} ne 0) {  #Allow no supplied email info
			
			&CheckEmailAddress(\${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}});
                    }
                }
		if (${$argHashRef}{'parameterName'} eq "instanceTag") {

		    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{${$argHashRef}{'parameterName'}} = ${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}};
		}
		if (${$argHashRef}{'parameterName'} eq "researchEthicalApproval") {

		    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{${$argHashRef}{'parameterName'}} = ${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}};
		}
	    }
	    
	    if ( ${$argHashRef}{'parameterType'} eq "program") {  #Evaluate "program" parameters

		if(${$argHashRef}{'parameterValue'} eq "nocmdinput") {  #No input from cmd
		    
		    if (defined(${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}})) {  #Input from config file - do nothing
			
		    }
		    elsif (${$argHashRef}{'parameterDefault'} ne "nodefault") {
			    
			${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}} = ${$argHashRef}{'parameterDefault'};  #Set default value
		    }
		}
		else {

		    ${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}} = ${$argHashRef}{'parameterValue'};
		}
		###Code for checking commands in your path and executable
		if ( (defined(${$argHashRef}{'programNamePath'})) && (defined(${$parameterHashRef}{${$argHashRef}{'parameterName'}}{ ${$argHashRef}{'programNamePath'} }[0])) ) {
		
		    if (${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}} > 0) {  #Only check path(s) for active programs
		
			for (my $programNamePathCounter=0;$programNamePathCounter<scalar(@{${$parameterHashRef}{${$argHashRef}{'parameterName'}}{ ${$argHashRef}{ ${$argHashRef}{'programNamePath'} } }});$programNamePathCounter++) {  #Check all binaries for sbatch program
			    
			    if ( grep { -x "$_/".${$parameterHashRef}{${$argHashRef}{'parameterName'}}{ ${$argHashRef}{'programNamePath'} }[$programNamePathCounter]} split(/:/,$ENV{PATH}) ) {
				
				$logger->info("ProgramCheck: ".${$parameterHashRef}{${$argHashRef}{'parameterName'}}{ ${$argHashRef}{'programNamePath'} }[$programNamePathCounter]." installed\n"); 
			    }
			    else {
				$logger->fatal("Could not detect ".${$parameterHashRef}{${$argHashRef}{'parameterName'}}{ ${$argHashRef}{'programNamePath'} }[$programNamePathCounter]." in your Path\n");
				exit;
			    }
			}
		    }
		}
	    }
	    
	    if (${$argHashRef}{'parameterName'} eq "aligner") {
		
		if ( (${$scriptParameterHashRef}{'pMosaikBuild'} > 0) || (${$scriptParameterHashRef}{'pMosaikAlign'} > 0)) {  #Mosaik track
		    
		    if ( (${$scriptParameterHashRef}{'pBwaAln'} == 0) && (${$scriptParameterHashRef}{'pBwaSampe'} == 0) && (${$scriptParameterHashRef}{'pBwaMem'} == 0) ) {
			
			if (${$scriptParameterHashRef}{'aligner'} eq "bwa") {
			    
			    ${$scriptParameterHashRef}{'aligner'} = "mosaik";
			}
		    }
		    else {
		
			$logger->fatal($USAGE, "\n");
			$logger->fatal("You have to choose either mosaik or bwa to perform alignments or specify which aligner (-aligner 'mosaik' or 'bwa') was used if you want to only run programs after alignment.", "\n");
			exit;
		    }
		}
		elsif ( (${$scriptParameterHashRef}{'pBwaAln'} > 0) || (${$scriptParameterHashRef}{'pBwaSampe'} > 0) || (${$scriptParameterHashRef}{'pBwaMem'} > 0)) {  #BWA track
		    
		    if ( (${$scriptParameterHashRef}{'aligner'} eq "mosaik") || (${$scriptParameterHashRef}{'aligner'} =~ /bwa/i) ) {

			${$scriptParameterHashRef}{'aligner'} = "bwa";
		    }
		    else {

			$logger->fatal($USAGE, "\n");
			$logger->fatal("You have to choose either mosaik or bwa to perform alignments or specify which aligner (-aligner 'mosaik' or 'bwa') was used if you want to only run programs after alignment.", "\n");
			exit;
		    }
		}
		elsif (${$scriptParameterHashRef}{'aligner'} eq "nocmdinput") {

		    $logger->fatal($USAGE, "\n");
		    $logger->fatal("You have to choose either mosaik or bwa to perform alignments or specify which aligner (-aligner 'mosaik' or 'bwa') was used if you want to only run programs after alignment.", "\n");
		    exit;
		}
	    }
	}
	if ($parameterSetSwitch eq 1) {  #No need to set parameter more than once
	    last;
	}
    }	
##Parameter set
    if (defined(${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}})) {
	
	my $info = "";  #Hold parameters info

	if (defined(${$parameterHashRef}{${$argHashRef}{'parameterName'}}{'array'})) {
	 
	    $info = "Set ".${$argHashRef}{'parameterName'}." to: ".join(',',@{${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}}});
	    push(@{$broadcastsArrayRef}, $info);  #Add info to broadcasts
	}
	else {

	    $info = "Set ".${$argHashRef}{'parameterName'}." to: ".${$scriptParameterHashRef}{${$argHashRef}{'parameterName'}};
	    push(@{$broadcastsArrayRef}, $info);  #Add info to broadcasts
	}
    }
}


sub CreateFileEndings {

##CreateFileEndings
    
##Function : Creates the fileEndings depending on which modules are used by the user to relevant chain.
##Returns  : ""
##Arguments: $scriptParameterHashRef, $orderParametersArrayRef
##         : $parameterHashRef           => The parameter hash {REF}
##         : $scriptParameterHashRef     => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef          => Info on samples and family hash {REF}
##         : $infilesLaneNoEndingHashRef => The infile(s) without the ".ending" {REF}
##         : $orderParametersArrayRef    => Order of addition to parameter array {REF}

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleInfoHashRef = $_[2];
    my $infilesLaneNoEndingHashRef = $_[3];
    my $orderParametersArrayRef = $_[4];

    my %tempFileEnding;  #Used to enable seqential build-up of fileEndings between modules
    
    foreach my $orderParameterElement (@{$orderParametersArrayRef}) {
	
	if (defined(${$scriptParameterHashRef}{$orderParameterElement})) {  #Only active parameters

	    if ( ($orderParameterElement =~ /^p[A-Z]/) && (${$parameterHashRef}{$orderParameterElement}{'associatedProgram'}) ) {  #Only process programs

		if (${$parameterHashRef}{$orderParameterElement}{'chain'} eq "MAIN") {  #MAIN chain
		    
		    if (${$parameterHashRef}{$orderParameterElement}{'fileEnding'} ne "nofileEnding") {  #FileEnding exist
			
###MAIN/Per sampleID
			for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{${$scriptParameterHashRef}{'sampleIDs'}});$sampleIDCounter++) {
			    
			    if (${$scriptParameterHashRef}{$orderParameterElement} > 0) {  #Fileending should be added    
				
				if ($orderParameterElement eq "pPicardToolsMergeSamFiles") {  #Special case
				    
				    if ( (defined(${$scriptParameterHashRef}{'picardToolsMergeSamFilesPrevious'})) || (scalar( @{ ${$infilesLaneNoEndingHashRef}{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] } }) > 1) ) {  #Sanity check that we have something to merge and hence to fileEnding should be added
					
					${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }{'pPicardToolsMergeSamFiles'}{'fileEnding'} = $tempFileEnding{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }.${$parameterHashRef}{$orderParameterElement}{'fileEnding'};  #Adds from previous entry 
				    }
				    else {
					
					${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }{'pPicardToolsMergeSamFiles'}{'fileEnding'} = $tempFileEnding{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }."";
				    }
				}
				else {

				    if (defined($tempFileEnding{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] })) {
					
					${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }{$orderParameterElement}{'fileEnding'} = $tempFileEnding{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }.${$parameterHashRef}{$orderParameterElement}{'fileEnding'};
				    }
				    else  {  #First module that should add filending

					${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }{$orderParameterElement}{'fileEnding'} = ${$parameterHashRef}{$orderParameterElement}{'fileEnding'};
				    } 
				}
			    }
			    else {  #Do not add new module fileEnding

				${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }{$orderParameterElement}{'fileEnding'} = $tempFileEnding{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] };
			    }
			    $tempFileEnding{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] } = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }{$orderParameterElement}{'fileEnding'};  #To enable sequential build-up of fileending
			}
			
###MAIN/Per familyID
			if (${$scriptParameterHashRef}{$orderParameterElement} > 0) {  #Fileending should be added 

			    if ($orderParameterElement eq "pPicardToolsMergeSamFiles") {  #Special case - do nothing
			    }
			    elsif ( ($orderParameterElement eq "pPicardToolsSortSam") && (${$scriptParameterHashRef}{'analysisType'} eq "rapid") ) {  #Special case - do nothing
			    }
			    elsif ( ($orderParameterElement eq "pPicardToolsMergeRapidReads") && (${$scriptParameterHashRef}{'analysisType'} ne "rapid") ) {  #Special case - do nothing
			    }
			    else {
				
				if (defined($tempFileEnding{${$scriptParameterHashRef}{'familyID'}})) {
			
				    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{$orderParameterElement}{'fileEnding'} = $tempFileEnding{${$scriptParameterHashRef}{'familyID'}}.${$parameterHashRef}{$orderParameterElement}{'fileEnding'};
				}
				else  {  #First module that should add filending
				    
				    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{$orderParameterElement}{'fileEnding'} = ${$parameterHashRef}{$orderParameterElement}{'fileEnding'};
				}
				$tempFileEnding{ ${$scriptParameterHashRef}{'familyID'} } = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{$orderParameterElement}{'fileEnding'};  #To enable sequential build-up of fileending 
			    }		
			}
			else {  #Do not add new module fileEnding
			 
			    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{  ${$scriptParameterHashRef}{'familyID'} }{$orderParameterElement}{'fileEnding'} = $tempFileEnding{ ${$scriptParameterHashRef}{'familyID'} };
			}
		    }
		}
		if (${$parameterHashRef}{$orderParameterElement}{'chain'} ne "MAIN") {  #Other chain(s)
		    
		    my $chainfork = ${$parameterHashRef}{$orderParameterElement}{'chain'}; 

		    if (${$parameterHashRef}{$orderParameterElement}{'fileEnding'} ne "nofileEnding") {  #FileEnding exist
			
###OTHER/Per sampleID
			for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{${$scriptParameterHashRef}{'sampleIDs'}});$sampleIDCounter++) {
			    
			    if (${$scriptParameterHashRef}{$orderParameterElement} > 0) {  #Fileending should be added    
			
				unless (defined($tempFileEnding{$chainfork}{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] })) {	

				    $tempFileEnding{$chainfork}{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] } = $tempFileEnding{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] };  #Inherit current MAIN chain. 
				}
				if (defined($tempFileEnding{$chainfork}{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] })) {

				    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }{$orderParameterElement}{'fileEnding'} = $tempFileEnding{$chainfork}{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }.${$parameterHashRef}{$orderParameterElement}{'fileEnding'};
				}
				else  {  #First module that should add filending

				    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }{$orderParameterElement}{'fileEnding'} = ${$parameterHashRef}{$orderParameterElement}{'fileEnding'};
				} 
			    }
			    else {  #Do not add new module fileEnding

				${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }{$orderParameterElement}{'fileEnding'} = $tempFileEnding{$chainfork}{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] };
			    }
			    ##NOTE: No sequential build-up of fileending
			}
###Other/Per familyID

			if (${$scriptParameterHashRef}{$orderParameterElement} > 0) {  #File ending should be added
			    
			    unless (defined($tempFileEnding{$chainfork}{${$scriptParameterHashRef}{'familyID'}})) {	

				$tempFileEnding{$chainfork}{${$scriptParameterHashRef}{'familyID'}} =  $tempFileEnding{${$scriptParameterHashRef}{'familyID'}};  #Inherit current MAIN chain. 
			    }
			    if (defined($tempFileEnding{$chainfork}{${$scriptParameterHashRef}{'familyID'}})) {

				${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{$orderParameterElement}{'fileEnding'} = $tempFileEnding{$chainfork}{${$scriptParameterHashRef}{'familyID'}}.${$parameterHashRef}{$orderParameterElement}{'fileEnding'};
			    }
			    else  {  #First module that should add filending

				${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{$orderParameterElement}{'fileEnding'} = ${$parameterHashRef}{$orderParameterElement}{'fileEnding'};
			    }
			    $tempFileEnding{$chainfork}{${$scriptParameterHashRef}{'familyID'}} = ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'familyID'} }{$orderParameterElement}{'fileEnding'};  #To enable sequential build-up of fileending 
			}
		    }
		}
	    }
	}
    }
}


sub ProgramPreRequisites {

##ProgramPreRequisites
    
##Function : Creates program directories (info & programData & programScript), program script filenames and writes sbatch header.
##Returns  : Path to stdout
##Arguments: $scriptParameterHashRef, $FILEHANDLE, $argHashRef
##         : $scriptParameterHashRef         => The active parameters for this analysis hash {REF}
##         : $FILEHANDLE                     => FILEHANDLE to write to
##         : $argHashRef{'directoryID'}      => $samplID|$familyID
##         : $argHashRef{'programName'}      => Assigns filename to sbatch script
##         : $argHashRef{'programDirectory'} => Builds from $directoryID/$aligner
##         : $argHashRef{'callType'}         => SNV,INDEL or BOTH
##         : $argHashRef{'nrofCores'}        => The number of cores to allocate
##         : $argHashRef{'processTime'}      => Hours
##         : $argHashRef{'tempDirectory'}    => Temporary directory for program {Optional}

    my $scriptParameterHashRef = shift;
    my $FILEHANDLE = shift;
    my ($argHashRef) = @_;


    my %default = ('callType' => "",
		   'nrofCores' => 1,
		   'processTime' => 1
	);
    
    if (defined(${$argHashRef}{'callType'})) {
	
	${$argHashRef}{'callType'} = "_".${$argHashRef}{'callType'};
    }
	
    &SetDefaultArg(\%{$argHashRef}, \%default);

    my $fileNameEnd = ".sh";
    my $fileName;  #The sbatch script - to be created filename
    my $fileNamePath;
    my $dryRunFilenamePath;
    my $programDataDirectory;
    my $fileInfoPath;
    my $dryRunFileInfoPath;
    my $fileNameTracker;

###Sbatch script names and directory creation

    $programDataDirectory = ${$scriptParameterHashRef}{'outDataDir'}."/".${$argHashRef}{'directoryID'}."/".${$argHashRef}{'programDirectory'};
    $fileNamePath = ${$scriptParameterHashRef}{'outScriptDir'}."/".${$argHashRef}{'directoryID'}."/".${$argHashRef}{'programDirectory'}."/".${$argHashRef}{'programName'}."_".${$argHashRef}{'directoryID'};
    $dryRunFilenamePath = ${$scriptParameterHashRef}{'outScriptDir'}."/".${$argHashRef}{'directoryID'}."/".${$argHashRef}{'programDirectory'}."/dry_run_".${$argHashRef}{'programName'}."_".${$argHashRef}{'directoryID'};
    $fileInfoPath = ${$scriptParameterHashRef}{'outDataDir'}."/".${$argHashRef}{'directoryID'}."/".${$argHashRef}{'programDirectory'}."/info/".${$argHashRef}{'programName'}."_".${$argHashRef}{'directoryID'};
    $dryRunFileInfoPath = ${$scriptParameterHashRef}{'outDataDir'}."/".${$argHashRef}{'directoryID'}."/".${$argHashRef}{'programDirectory'}."/info/dry_run_".${$argHashRef}{'programName'}."_".${$argHashRef}{'directoryID'};

    $fileNamePath .= ${$argHashRef}{'callType'}.".";
    $dryRunFilenamePath .= ${$argHashRef}{'callType'}.".";
    $fileInfoPath .= ${$argHashRef}{'callType'}.".";
    $dryRunFileInfoPath .= ${$argHashRef}{'callType'}.".";

    ## Create directories
    `mkdir -p ${$scriptParameterHashRef}{'outDataDir'}/${$argHashRef}{'directoryID'}/${$argHashRef}{'programDirectory'}/info;`;  #Creates the aligner folder and info data file directory
    `mkdir -p $programDataDirectory`;  #Creates the aligner folder and if supplied the program data file directory
    `mkdir -p ${$scriptParameterHashRef}{'outScriptDir'}/${$argHashRef}{'directoryID'}/${$argHashRef}{'programDirectory'}`;  #Creates the aligner folder script file directory

    if ( (${$scriptParameterHashRef}{"p".${$argHashRef}{'programName'}} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	$fileName = $fileNamePath; 
    }
    elsif (${$scriptParameterHashRef}{"p".${$argHashRef}{'programName'}} == 2) {  #Dry run single program

	$fileName = $dryRunFilenamePath; 
	$logger->info("Dry Run:\n");
    }
    else {  #Dry run

	$fileName = $dryRunFilenamePath;
	$logger->info("Dry Run:\n");
    }

    ($fileName, $fileNameTracker) = &CheckFileNameExists(\$fileName, \$fileNameEnd);

###Info and Log
    $logger->info("Creating sbatch script for ".${$argHashRef}{'programName'}." and writing script file(s) to: ".$fileName."\n");
    $logger->info("Sbatch script ".${$argHashRef}{'programName'}." data files will be written to: ".$programDataDirectory."\n");

###Sbatch header
    open ($FILEHANDLE, ">",$fileName) or $logger->logdie("Can't write to '".$fileName."' :".$!."\n");
    
    print $FILEHANDLE "#! /bin/bash -l", "\n";
    print $FILEHANDLE "#SBATCH -A ".${$scriptParameterHashRef}{'projectID'}, "\n";
    print $FILEHANDLE "#SBATCH -n ".${$argHashRef}{'nrofCores'}, "\n";
    print $FILEHANDLE "#SBATCH -t ".${$argHashRef}{'processTime'}.":00:00", "\n";	
    print $FILEHANDLE "#SBATCH -J ".${$argHashRef}{'programName'}."_".${$argHashRef}{'directoryID'}.${$argHashRef}{'callType'}, "\n";

    if ( (${$scriptParameterHashRef}{"p".${$argHashRef}{'programName'}} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} == 0) ) {

	print $FILEHANDLE "#SBATCH -e ".$fileInfoPath.$fileNameTracker.".stderr.txt", "\n";
	print $FILEHANDLE "#SBATCH -o ".$fileInfoPath.$fileNameTracker.".stdout.txt", "\n";
    }
    elsif (${$scriptParameterHashRef}{'pSampleCheck'} == 2) {  #Single program dry run

	print $FILEHANDLE "#SBATCH -e ".$dryRunFileInfoPath.$fileNameTracker.".stderr.txt", "\n";
	print $FILEHANDLE "#SBATCH -o ".$dryRunFileInfoPath.$fileNameTracker.".stdout.txt", "\n";
    }
    else {  #Dry run

	print $FILEHANDLE "#SBATCH -e ".$dryRunFileInfoPath.$fileNameTracker.".stderr.txt", "\n";
	print $FILEHANDLE "#SBATCH -o ".$dryRunFileInfoPath.$fileNameTracker.".stdout.txt", "\n";
    }
    
    unless (${$scriptParameterHashRef}{'email'} eq 0) {
	
	if (${$scriptParameterHashRef}{'emailType'} =~/B/i) {

	    print $FILEHANDLE "#SBATCH --mail-type=BEGIN", "\n";
	}
	if (${$scriptParameterHashRef}{'emailType'} =~/E/i) {
	 
	    print $FILEHANDLE "#SBATCH --mail-type=END", "\n";
	}
	if (${$scriptParameterHashRef}{'emailType'} =~/F/i) {
	    
	    print $FILEHANDLE "#SBATCH --mail-type=FAIL", "\n";
	}
	print $FILEHANDLE "#SBATCH --mail-user=".${$scriptParameterHashRef}{'email'}, "\n\n";	
    }
    
    print $FILEHANDLE q?echo "Running on: $(hostname)"?,"\n\n";

    if (defined(${$argHashRef}{'tempDirectory'})) {  #Not all programs need a temporary directory

	print $FILEHANDLE "## Create temporary directory\n";
	print $FILEHANDLE "mkdir -p ".${$argHashRef}{'tempDirectory'}, "\n\n";
    }
    return ($fileName, $fileInfoPath.$fileNameTracker.".stdout.txt", $fileInfoPath.$fileNameTracker.".stderr.txt");  #Return filen ame, stdout, stderr path for QC check later
}



sub CheckIfMergedFiles {

##CheckIfMergedFiles
    
##Function : Check if any files for this sampleID were merged previously to set infile and PicardToolsMergeSwitch to enable correct handling of number of infiles to process
##Returns  : "$infile, $PicardToolsMergeSwitch"
##Arguments: $scriptParameterHashRef, $sampleInfoHashRef, $laneHashRef, $infilesLaneNoEndingHashRef, $sampleID
##         : $scriptParameterHashRef     => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef          => Info on samples and family hash {REF}
##         : $laneHashRef                => The lane info hash {REF}
##         : $infilesLaneNoEndingHashRef => The infile(s) without the ".ending" {REF}
##         : $sampleID                   => The sampleID

    my $scriptParameterHashRef = $_[0];
    my $sampleInfoHashRef = $_[1];
    my $laneHashRef = $_[2];
    my $infilesLaneNoEndingHashRef = $_[3];
    my $sampleID = $_[4];

    my $infile;
    my $mergeLanes;  #To pick up merged lanes later 
    my $PicardToolsMergeSwitch = 0;  #0=no merge was previously performed

    if (${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$sampleID}{'picardToolsMergeSamFilesPrevious'} == 1) {  # Files merged this round with merged file from previous round
	
	for (my $mergeFileCounter=0;$mergeFileCounter<scalar(@{${$scriptParameterHashRef}{'picardToolsMergeSamFilesPrevious'}});$mergeFileCounter++) {
	    
	    if (${$scriptParameterHashRef}{'picardToolsMergeSamFilesPrevious'}[$mergeFileCounter] =~ /lane(\d+)|s_(\d+)/) {  #Look for lanes_ or lane\d in previously generated file to be merged with current run to be able to extract previous lanes
		
		##Make sure to always supply lanes from previous regexp 
		if($1) {
		
		    $mergeLanes = $1;
		} 
		else {

		    $mergeLanes = $2;
		}  
		$infile = $sampleID."_lanes_".$mergeLanes;

		for (my $laneCounter=0;$laneCounter<scalar(@ { ${$laneHashRef}{$sampleID} });$laneCounter++) {
		
		    $infile .= ${$laneHashRef}{$sampleID}[$laneCounter];  #Extract lanes per sampleID
		}
		$PicardToolsMergeSwitch = 1;
	    }
	}
    }
    elsif ( (${$scriptParameterHashRef}{'pPicardToolsMergeSamFiles'} > 0) && (scalar( @{ ${$infilesLaneNoEndingHashRef}{$sampleID} }) > 1) ) {  #But only if there is more than one mosaikBuild/BWA_Aln file per sample ID (Sanity check)

	$infile = $sampleID."_lanes_";

	for (my $laneCounter=0;$laneCounter<scalar(@ { ${$laneHashRef}{$sampleID} });$laneCounter++) {
	   
	    $infile .= ${$laneHashRef}{$sampleID}[$laneCounter];  #Extract lanes per sampleID
	}
	$PicardToolsMergeSwitch = 1;
    }
    else {

	$PicardToolsMergeSwitch = 0;
    }
    return ($infile, $PicardToolsMergeSwitch);
}


sub SampleInfoQC {

##SampleInfoQC
    
##Function : Adds outDirectory and outFile to sampleInfo to track all files that QC metrics are to be extracted from later
##Returns  : ""
##Arguments: $sampleInfoHashRef, $argHashRef
##         : $sampleInfoHashRef           => Info on samples and family hash {REF}
##         : $argHashRef{'familyID'}      => The familyID
##         : $argHashRef{'sampleID'}      => SampleID or "noSampleID" for family level data
##         : $argHashRef{'programName'}   => The program
##         : $argHashRef{'infile'}        => Infile or "noInFile" for family level data
##         : $argHashRef{'outDirectory'}  => The outdirectory of the QC file
##         : $argHashRef{'outFileEnding'} => The outfile ending. Actually complete outfile for "static" & "infoDirectory"
##         : $argHashRef{'outDataType'}   => Type of data produced by program (infoDirectory|infileDependent|static) 

    my $sampleInfoHashRef = shift;
    my ($argHashRef) = @_;

    unless (defined(${$argHashRef}{'sampleID'})) {

	${$sampleInfoHashRef}{ ${$argHashRef}{'familyID'} }{ ${$argHashRef}{'familyID'} }{'program'}{ ${$argHashRef}{'programName'} }{'OutDirectory'} = ${$argHashRef}{'outDirectory'};  #OutDirectory of QC file
                                                            
	if (${$argHashRef}{'outDataType'} eq "static") {  #Programs which add a static file in its own directory                                                                                                 

	    ${$sampleInfoHashRef}{ ${$argHashRef}{'familyID'} }{ ${$argHashRef}{'familyID'} }{'program'}{ ${$argHashRef}{'programName'} }{'OutFile'} = ${$argHashRef}{'outFileEnding'};  #Static QC outFile                                                                     
	}
	if (${$argHashRef}{'outDataType'} eq "infoDirectory") {  #QC metrics are sent to info files                                                                                                                   

	    ${$sampleInfoHashRef}{ ${$argHashRef}{'familyID'} }{ ${$argHashRef}{'familyID'} }{'program'}{ ${$argHashRef}{'programName'} }{'OutFile'} = ${$argHashRef}{'outFileEnding'};  #Info stdout file                                                                      
	}
	if (${$argHashRef}{'outDataType'} eq "infileDependent") {  #Programs which add a filending to infile                                                                                                          

	    ${$sampleInfoHashRef}{ ${$argHashRef}{'familyID'} }{ ${$argHashRef}{'familyID'} }{'program'}{ ${$argHashRef}{'programName'} }{'OutFile'} = ${$argHashRef}{'outFileEnding'};  #Infile dependent QC outFile                                                                                                                                                                                       
	}

    }
    else {
	
	${$sampleInfoHashRef}{ ${$argHashRef}{'familyID'} }{ ${$argHashRef}{'sampleID'} }{'program'}{ ${$argHashRef}{'programName'} }{ ${$argHashRef}{'infile'} }{'OutDirectory'} = ${$argHashRef}{'outDirectory'};  #OutDirectory of QC file                                                              

	if (${$argHashRef}{'outDataType'} eq "static") {  #Programs which add a static file in its own directory 
	    
	    ${$sampleInfoHashRef}{ ${$argHashRef}{'familyID'} }{ ${$argHashRef}{'sampleID'} }{'program'}{ ${$argHashRef}{'programName'} }{ ${$argHashRef}{'infile'} }{'OutFile'} = ${$argHashRef}{'outFileEnding'};  #Static QC outFile
	}
	if (${$argHashRef}{'outDataType'} eq "infoDirectory") {  #QC metrics are sent to info files
	    
	    ${$sampleInfoHashRef}{ ${$argHashRef}{'familyID'} }{ ${$argHashRef}{'sampleID'} }{'program'}{ ${$argHashRef}{'programName'} }{ ${$argHashRef}{'infile'} }{'OutFile'} = ${$argHashRef}{'outFileEnding'};  #Info stdout file
	}
	if (${$argHashRef}{'outDataType'} eq "infileDependent") {  #Programs which add a filending to infile
	    
	    ${$sampleInfoHashRef}{ ${$argHashRef}{'familyID'} }{ ${$argHashRef}{'sampleID'} }{'program'}{ ${$argHashRef}{'programName'} }{ ${$argHashRef}{'infile'} }{'OutFile'} = ${$argHashRef}{'infile'}.${$argHashRef}{'outFileEnding'};  #Infile dependent QC outFile                                                                      
	}
    }
}


sub GATKTargetListFlag {

##GATKTargetListFlag
    
##Function : Detects if there are different capture kits across sampleIDs. Creates a temporary merged interval_list for all interval_list that have been supplied and returns temporary list. Will also extract specific contigs if requested and return that list if enabled.
##Returns  : "Filepath"
##Arguments: $scriptParameterHashRef, $FILEHANDLE, $contigRef
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $FILEHANDLE             => FILEHANDLE to write to
##         : $contigRef              => The contig to extract {REF}

    my $scriptParameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];
    my $contigRef = $_[2];

    my %GATKTargetPaddedBedIntervalListTracker;
    my @GATKTargetPaddedBedIntervalListFiles;

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{${$scriptParameterHashRef}{'sampleIDs'}});$sampleIDCounter++) {  #Collect infiles for all sampleIDs
	
	if (defined(${$scriptParameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }{'GATKTargetPaddedBedIntervalLists'})) {

	    ${$scriptParameterHashRef}{'GATKTargetPaddedBedIntervalLists'} = ${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }{'GATKTargetPaddedBedIntervalLists'};  #Transfer to scriptParameter top level
       
	    $GATKTargetPaddedBedIntervalListTracker{ ${$scriptParameterHashRef}{'GATKTargetPaddedBedIntervalLists'} }++;  #Increment to track file record
	    
	    if ($GATKTargetPaddedBedIntervalListTracker{ ${$scriptParameterHashRef}{'GATKTargetPaddedBedIntervalLists'} } == 1) {  #Not detected previously
		
		push(@GATKTargetPaddedBedIntervalListFiles, ${$scriptParameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }{'GATKTargetPaddedBedIntervalLists'});
	    }
	}
    }
    
    ##Determine file to print to module (untouched/merged and/or splited)
    my $outDirectory = ${$scriptParameterHashRef}{'tempDirectory'};  #For merged and/or splitet

    if (scalar(@GATKTargetPaddedBedIntervalListFiles) > 1) {  #Merge files
      
	print $FILEHANDLE "\n#Generate merged interval_list\n\n"; 
	print $FILEHANDLE "java -Xmx2g ";

	&WriteUseLargePages($FILEHANDLE, \${$scriptParameterHashRef}{'javaUseLargePages'});
	
	print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'picardToolsPath'}."/IntervalListTools.jar ";
	print $FILEHANDLE "UNIQUE=TRUE ";  #Merge overlapping and adjacent intervals to create a list of unique intervals
    
	for (my $fileCounter=0;$fileCounter<scalar(@GATKTargetPaddedBedIntervalListFiles);$fileCounter++) {
	
	    print $FILEHANDLE "INPUT=".${$scriptParameterHashRef}{'referencesDir'}."/".$GATKTargetPaddedBedIntervalListFiles[$fileCounter]." ";
	}
	print $FILEHANDLE "OUTPUT=".$outDirectory."/merged.interval_list", "\n\n";  #Merged outfile

	if (defined($$contigRef)) {
	    
	    my $inDirectory = ${$scriptParameterHashRef}{'tempDirectory'};
	    return &SplitTargetFile(*$FILEHANDLE, \$inDirectory, \$outDirectory, \("merged.interval_list"), \$$contigRef);  #Split
	}
        return $outDirectory."/merged.interval_list";  #No split
    }
    elsif (defined($$contigRef)) {  #Supply original file but create splitted temp file

	return &SplitTargetFile(*$FILEHANDLE, \${$scriptParameterHashRef}{'referencesDir'}, \$outDirectory, \$GATKTargetPaddedBedIntervalListFiles[0], \$$contigRef);  #Only 1 file for all samples
    }
    else {#No merge and no split. return original and only file

	return  ${$scriptParameterHashRef}{'referencesDir'}."/".$GATKTargetPaddedBedIntervalListFiles[0];
    }
}


sub SplitTargetFile {

##SplitTargetFile
    
##Function : Splits a target file into new contig specific target file
##Returns  : "Filepath to splitted file"
##Arguments: $FILEHANDLE, $inDirectoryRef, $outDirectoryRef, $infileRef, $contigRef
##         : $FILEHANDLE      => FILEHANDLE to write to
##         : $inDirectoryRef  => Indirectory {REF}
##         : $outDirectoryRef => Analysis outdirectory {REF}
##         : $infileRef       => Target file {REF}
##         : $contigRef       => The contig to extract {REF}

    my $FILEHANDLE = $_[0];
    my $inDirectoryRef = $_[1];
    my $outDirectoryRef = $_[2];
    my $infileRef = $_[3];
    my $contigRef = $_[4];
    
    if (defined($$contigRef)) {  #The contig to split
	
	print $FILEHANDLE "\n#Generate contig specific interval_list\n\n"; 
	print $FILEHANDLE q?perl -nae 'if($_=~/^\@/) {print $_;} elsif($_=~/^?.$$contigRef.q?\s+/) {print $_;}' ?;  #Select header and contig
	print $FILEHANDLE $$inDirectoryRef."/".$$infileRef." ";  #Infile
	print $FILEHANDLE "> ".$$outDirectoryRef."/".$$contigRef."_".$$infileRef, "\n\n";#Extract genomic interval info
	
	return $$outDirectoryRef."/".$$contigRef."_".$$infileRef;
    }
}


sub GATKPedigreeFlag {

##GATKPedigreeFlag
    
##Function : Check if "--pedigree" and "--pedigreeValidationType" should be included in analysis
##Returns  : ""
##Arguments: $scriptParameterHashRef, $FILEHANDLE, $outFamilyFileDirectory, $pedigreeValidationType, $program
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $FILEHANDLE             => FILEHANDLE to write to
##         : $outFamilyFileDirectory => The family data analysis directory 
##         : $pedigreeValidationType => The pedigree validation strictness level
##         : $program                => The program to use the pedigree file

    my $scriptParameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];
    my $outFamilyFileDirectory = $_[2];
    my $pedigreeValidationType = $_[3];
    my $program = $_[4];
    
    my $famFile = $outFamilyFileDirectory."/".${$scriptParameterHashRef}{'familyID'}.".fam";
    my $parentCounter;
    my $pqParentCounter = q?perl -ne 'my $parentCounter=0; while (<>) { my @line = split(/\t/, $_); unless ($_=~/^#/) { if ( ($line[2] eq 0) || ($line[3] eq 0) ) { $parentCounter++} } } print $parentCounter; last;'?;
    my $childCounter;
    my $pqChildCounter = q?perl -ne 'my $childCounter=0; while (<>) { my @line = split(/\t/, $_); unless ($_=~/^#/) { if ( ($line[2] ne 0) || ($line[3] ne 0) ) { $childCounter++} } } print $childCounter; last;'?;
    
    $parentCounter = `$pqParentCounter $famFile`;  #Count the number of parents
    $childCounter = `$pqChildCounter $famFile`;  #Count the number of children
    
    if ($program ne "GATKPhaseByTransmission") {
	
	if ($parentCounter > 0) {  #Parents present
	    
	    print $FILEHANDLE "--pedigreeValidationType ".$pedigreeValidationType." --pedigree ".$outFamilyFileDirectory."/".${$scriptParameterHashRef}{'familyID'}.".fam ";  #Pedigree files for samples		
	}
    }
    else {
	
	&CheckPedigreeMembers(\%{$scriptParameterHashRef}, $FILEHANDLE, \$outFamilyFileDirectory, \$pedigreeValidationType, \$parentCounter, \$childCounter);  #Special case - GATK PhaseByTransmission needs parent/child or trio 
    }
}


sub CheckPedigreeMembers {

##CheckPedigreeMembers
    
##Function : Detect if the pedigree file contains a valid parent/child or trio
##Returns  : ""
##Arguments: $scriptParameterHashRef, $FILEHANDLE, $outFamilyFileDirectoryRef, $pedigreeValidationTypeRef, $parentCounterRef, $childCounterRef
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $FILEHANDLE             => FILEHANDLE to write to
##         : $outFamilyFileDirectoryRef => The family data analysis directory {REF}
##         : $pedigreeValidationTypeRef => The pedigree validation strictness level {REF}
##         : $parentCounterRef       => The number of parent(s) {REF}
##         : $childCounterRef        => The number of children(s) {REF}

    my $scriptParameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];
    my $outFamilyFileDirectoryRef = $_[2];
    my $pedigreeValidationTypeRef = $_[3];
    my $parentCounterRef = $_[4];
    my $childCounterRef = $_[5];
	    
    if (scalar(@{${$scriptParameterHashRef}{'sampleIDs'}}) < 4) {  #i.e.1-3 individuals in pedigree		    
		
	if ( ($$childCounterRef == 1) && ($$parentCounterRef > 0) ) {  #Parent/child or trio

	    print $FILEHANDLE "--pedigreeValidationType ".$$pedigreeValidationTypeRef." --pedigree ".$$outFamilyFileDirectoryRef."/".${$scriptParameterHashRef}{'familyID'}.".fam ";  #Pedigree files for samples
	}
	else {

	    ${$scriptParameterHashRef}{'pGATKPhaseByTransmission'} = 0;  #Override input since pedigree is not valid for analysis
	    $logger->info("Switched GATK PhaseByTransmission to 'no run' mode since MIP did not detect a valid pedigree for this type of analysis.");
	    
	    if (${$scriptParameterHashRef}{'pGATKReadBackedPhasing'} > 0) {  #Broadcast
		
		$logger->info("MIP will still try to run GATK ReadBackedPhasing, but with the '-respectPhaseInInput' flag set to false\n");
	    }
	}
    }
    else {
	
	${$scriptParameterHashRef}{'pGATKPhaseByTransmission'} = 0;  #Override input since pedigree is not valid for analysis
	$logger->info("Switched GATK PhaseByTransmission to 'no run' mode since MIP did not detect a valid pedigree for this type of analysis.");
	
	if (${$scriptParameterHashRef}{'pGATKReadBackedPhasing'} > 0) {  #Broadcast
	    
	    $logger->info("MIP will still try to run GATK ReadBackedPhasing, but with the '-respectPhaseInInput' flag set to false\n");
	}
    }
}


sub WriteCMDMipLog {

##WriteCMDMipLog
    
##Function : Write CMD to MIP log file
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $orderParametersArrayRef, $scriptRef, $logFileRef
##         : $parameterHashRef        => The parameter hash {REF}
##         : $scriptParameterHashRef  => The active parameters for this analysis hash {REF}
##         : $orderParametersArrayRef => Order of addition to parameter array {REF}
##         : $scriptRef               => The script that is being executed {REF}
##         : $logFileRef              => The log file
##         : $mipVersionRef           => The MIP version
    
    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $orderParametersArrayRef = $_[2];
    my $scriptRef = $_[3];
    my $logFileRef = $_[4];
    my $mipVersionRef = $_[5];

    my $cmdLine = $$scriptRef." ";

    foreach my $orderParameterElement (@{$orderParametersArrayRef}) {
	
	if (defined(${$scriptParameterHashRef}{$orderParameterElement}) ) {

	    if ( ($orderParameterElement eq "configFile") && (${$scriptParameterHashRef}{'configFile'} eq 0) ) {  #Do not print
	    }
	    else {

		if (defined(${$parameterHashRef}{$orderParameterElement}{'array'})) {  #Array parameters need to be comma sep 

		    $cmdLine .= "-".$orderParameterElement." ".join(',', @{${$scriptParameterHashRef}{$orderParameterElement}})." ";
		}
		else {
		    $cmdLine .="-".$orderParameterElement." ".${$scriptParameterHashRef}{$orderParameterElement}." ";
		}
	    }
	}
    }
    $logger->info($cmdLine,"\n");
    $logger->info("MIP Version: ".$$mipVersionRef, "\n");
    $logger->info("Script parameters and info from ".$$scriptRef." are saved in file: ".$$logFileRef, "\n");
}


sub WriteYAML {
 
##WriteYAML
    
##Function : Writes a YAML hash to file
##Returns  : ""
##Arguments: $yamlHashRef, $yamlFileRef
##         : $yamlHashRef => The hash to dump {REF}
##         : $yamlFileRef    => The yaml file to write to {REF}

    my $yamlHashRef = $_[0];  #Hash reference to write to file
    my $yamlFileRef = $_[1];  #Filename
    
    open (my $YAML, ">", $$yamlFileRef) or $logger->logdie("Can't open '".$$yamlFileRef."':".$!."\n");
    print $YAML Dump( $yamlHashRef ), "\n";
    close($YAML);

    $logger->info("Wrote: ".$$yamlFileRef, "\n");
}


sub LoadYAML {
 
##LoadYAML
    
##Function : Loads a YAML file into an arbitrary hash and returns it. Note: Currently only supports hashreferences and hashes and no mixed entries.
##Returns  : %yamlHash
##Arguments: $scriptParameterHashRef, $yamlFile
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $yamlFile => The yaml file to load

    my $scriptParameterHashRef = $_[0];
    my $yamlFile = $_[1];

    my %yamlHash;

    open (my $YAML, "<", $yamlFile) or die "can't open ".$yamlFile.":".$!, "\n";  #Log4perl not initialised yet, hence no logdie
    %yamlHash = %{ YAML::LoadFile($yamlFile) };  #Load hashreference as hash
        
    close($YAML);
    
    if (defined(${$scriptParameterHashRef}{'logFile'})) {

	$logger->info("Read Yaml file: ". $yamlFile, "\n");
    }
    return %yamlHash;
}


sub CheckUniqueArrayElement {

##CheckUniqueArrayElement
    
##Function : Detects if there are elements in arrayQueryRef that are not present in scalarQueryRef or arrayToCheckRef. If unique adds the unique element to arrayToCheckRef.
##Returns  : ""
##Arguments: $arrayToCheckRef, $arrayQueryRef, scalarQueryRef
##         : $arrayToCheckRef => The arrayref to be queried {REF}
##         : $arrayQueryRef   => The query array {REF}
##         : $scalarQueryRef  => The query scalar {REF}

    my $arrayToCheckRef = $_[0];
    my $arrayQueryRef;
    my $scalarQueryRef;

    if (ref($_[1]) eq "ARRAY") {  #Array reference

	$arrayQueryRef = $_[1];
	
	##For each arrayQueryRef element, loop through corresponding arrayToCheckRef element(s), add if there are none or an updated/unique entry.
	for (my $elementsInfoCounter=0;$elementsInfoCounter<scalar(@{$arrayQueryRef});$elementsInfoCounter++) {  #All element(s)
    
	    my $elementFound = 0; #Track if there element is present in arrayToCheckRef
	    
	    for (my $elementsCounter=0;$elementsCounter<scalar( @{$arrayToCheckRef});$elementsCounter++) {  #All arrayToCheckRef elements
		
		if (${$arrayToCheckRef}[$elementsCounter] eq ${$arrayQueryRef}[$elementsInfoCounter]) {  #Check presence

		    $elementFound = 1;   #Entry is present in both arrays
		}
	    }
	    if ($elementFound == 0) {  #Not seen in arrayToCheckRef

		push( @{$arrayToCheckRef}, ${$arrayQueryRef}[$elementsInfoCounter]);  #Go ahead and add	
	    }
	}
    }
    if (ref($_[1]) eq "SCALAR") {  #Scalar reference

	$scalarQueryRef = $_[1]; 

	my $elementFound = 0;  #Track if there element is present in arrayToCheckRef
	
	for (my $elementsCounter=0;$elementsCounter<scalar( @{$arrayToCheckRef});$elementsCounter++) {  #All arrayToCheckRef elements
	    
	    if (${$arrayToCheckRef}[$elementsCounter] eq $$scalarQueryRef) {  #Check presence
		
		$elementFound = 1;  #Entry is present in both arrays
	    }
	}
	if ($elementFound == 0) {  #Not seen in arrayToCheckRef

	    push( @{$arrayToCheckRef}, $$scalarQueryRef);  #Go ahead and add
	}
    }
}

sub DetermineNrofRapidNodes {

##DetermineNrofRapidNodes
    
##Function : Determines the number of nodes to allocate depending on the sequence read length, which affects the infile size.
##Returns  : $numberNodes, $ReadNrofLines
##Arguments: $seqLength, $infileSize
##         : $seqLength     => Length of sequence reads
##         : $infileSize    => Size of the infile

    my $seqLength = $_[0]; 
    my $infileSize = $_[1];
    
    my $numberNodes = 0;  #Nodes to allocate
    my $readPositionWeight = 1;  #Scales the readStart and readStop position
    my $ReadNrofLines;    

    if ($seqLength > 75 && $seqLength <= 101) {

	$ReadNrofLines = 190000000;  #Read batch size
	$numberNodes = floor($infileSize / (12 * $ReadNrofLines) );  #Determines the number of nodes to use, 150000000 ~ 37,5 million reads, 13 = 2 sdtdev from sample population - currently poor estimate with compression confunding calculation.	
	$logger->info("Number of Nodes: ".$numberNodes, "\n");
    }
    if ($seqLength > 50 && $seqLength <= 75) {

	$ReadNrofLines = 190000000;  #Read batch size
	$numberNodes = floor($infileSize / (9.75 * $ReadNrofLines) );  #Determines the number of nodes to use, 150000000 ~ 37,5 million reads, 13 = 2 sdtdev from sample population - currently poor estimate with compression confunding calculation.	
	$logger->info("Number of Nodes: ".$numberNodes, "\n");
    }
    if ($seqLength >= 50 && $seqLength < 75) {

	$ReadNrofLines = 130000000;  #Read batch size
	$numberNodes = floor($infileSize / (7 * $ReadNrofLines) );  #Determines the number of nodes to use, 150000000 ~ 37,5 million reads, 13 = 2 sdtdev from sample population - currently poor estimate with compression confunding calculation.
	$logger->info("Number of Nodes: ".$numberNodes, "\n");
    }
    if ($seqLength >= 35 && $seqLength < 50) {

	$ReadNrofLines = 95000000;  #Read batch size
	$numberNodes = floor($infileSize / (6 * $ReadNrofLines) );  #Determines the number of nodes to use, 150000000 ~ 37,5 million reads, 13 = 2 sdtdev from sample population - currently poor estimate with compression confunding calculation.
	$logger->info("Number of Nodes: ".$numberNodes, "\n");
    }
    if ($numberNodes <= 1) {
	
	$numberNodes = 2;  #Ensure that at least 1 readbatch is processed
    }
    return $numberNodes, $ReadNrofLines;
}


sub CheckUniqueIDNs {

##CheckUniqueIDNs
    
##Function : Test that the familyID and the sampleID(s) exists and are unique. Check if id sampleID contains "_".
##Returns  : "" 
##Arguments: $scriptParameterHashRef, $sampleIdArrayRef
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $sampleIDArrayRef => Array to loop in for parameter {REF}

    my $scriptParameterHashRef = $_[0];
    my $sampleIdArrayRef = $_[1];

    my %seen;  #Hash to test duplicate sampleIDs later

    if (scalar(@{$sampleIdArrayRef}) == 0) {

	$logger->fatal("Please provide sampleID(s)\n");
	exit;
    }

    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{$sampleIdArrayRef});$sampleIDCounter++) {

	$seen{ ${$sampleIdArrayRef}[$sampleIDCounter] }++;  #Increment instance to check duplicates later
	
	if (${$scriptParameterHashRef}{'familyID'} eq ${$sampleIdArrayRef}[$sampleIDCounter]) {  #FamilyID cannot be the same as sampleID
	    
	    $logger->fatal("FamilyID: ".${$scriptParameterHashRef}{'familyID'}." equals sampleID: ".${$sampleIdArrayRef}[$sampleIDCounter].". Please make sure that the familyID and sampleID(s) are unique.\n");
	    exit;
	}
	if ($seen{ ${$sampleIdArrayRef}[$sampleIDCounter] } > 1) {  #Check sampleID are unique
	
	    $logger->fatal("SampleID: ".${$sampleIdArrayRef}[$sampleIDCounter]." is not uniqe.\n");
	    exit;
	}
	if (${$sampleIdArrayRef}[$sampleIDCounter] =~/_/) {  #SampleID contains "_", which is not allowed accrding to filename conventions

	    $logger->fatal("SampleID: ".${$sampleIdArrayRef}[$sampleIDCounter]." contains '_'. Please rename sampleID according to MIP's filename convention, removing the '_'.\n");
	    exit;
	}
    }
}


sub UpdateYAML {

##UpdateYAML
    
##Function : Updates the config file to particular user/cluster for entries following specifications. Leaves other entries untouched.
##Returns  : "" 
##Arguments: $scriptParameterHashRef, $parameterNameRef, $familyID
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $parameterNameRef       => MIP Parameter to update {REF}
##         : $familyID               => Sets the familyID

    my $scriptParameterHashRef = $_[0];
    my $parameterNameRef = $_[1]; 
    my $familyIDRef = $_[2]; 
    
    if (${$scriptParameterHashRef}{$$parameterNameRef}) {  #Active parameter
	
	if (defined(${$scriptParameterHashRef}{'clusterConstantPath'})) {  #Set the project specific path for this cluster
	 
	    ${$scriptParameterHashRef}{$$parameterNameRef} =~ s/CLUSTERCONSTANTPATH!/${$scriptParameterHashRef}{'clusterConstantPath'}/gi;  #Exchange CLUSTERCONSTANTPATH! for current cluster path
	}
	if (defined(${$scriptParameterHashRef}{'analysisConstantPath'})) { #Set the project specific path for this cluster
	 
	    ${$scriptParameterHashRef}{$$parameterNameRef} =~ s/ANALYSISCONSTANTPATH!/${$scriptParameterHashRef}{'analysisConstantPath'}/gi;  #Exchange ANALYSISCONSTANTPATH! for the current analysis path
	}
	if (defined(${$scriptParameterHashRef}{'analysisType'})) {  #Set the analysis run type e.g., "exomes", "genomes", "rapid"

	    ${$scriptParameterHashRef}{$$parameterNameRef} =~ s/ANALYSISTYPE!/${$scriptParameterHashRef}{'analysisType'}/gi;  #Exchange ANALYSISTYPE! for the current analysis type
	}
	if (defined($$familyIDRef)) {  #Set the familyID

	    ${$scriptParameterHashRef}{$$parameterNameRef} =~ s/FDN!/$$familyIDRef/gi;  #Exchange FND! for the current familyID
	}
	if (defined(${$scriptParameterHashRef}{'aligner'})) {  #Set the aligner used

	    ${$scriptParameterHashRef}{$$parameterNameRef} =~ s/ALIGNER!/${$scriptParameterHashRef}{'aligner'}/gi;  #Exchange ALIGNER! for the current aligner
	}
    }
}


sub CheckAutoBuild {

##CheckAutoBuild
    
##Function : Checks if autobuild is on and returns "1" if enabled or "0" if not
##Returns  : "0|1" 
##Arguments: $parameterHashRef, $scriptParameterHashRef, $parameterNameRef, $sampleIDRef
##         : $parameterHashRef       => The parameter hash {REF}
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $parameterNameRef       => MIP parameter name {REF}
##         : $sampleIDRef            => SampleId {REF}

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $parameterNameRef = $_[2];
    my $sampleIDRef = $_[3];

    if (defined($sampleIDRef)) {
	
	if ( (${$parameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$$sampleIDRef}{$$parameterNameRef}{'buildFile'} eq "yesAutoBuild") || (${$parameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$$sampleIDRef}{$$parameterNameRef}{'buildFile'} eq "yesAutoDownLoad") || (${$parameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$$sampleIDRef}{$$parameterNameRef}{'buildFile'} eq 1) ) {
	    
	    return "1";  #Flag that autobuild is needed
	}
	else {
	    return "0";  #No autobuild is needed   
	}
    }
    else {
	
	if ( (${$parameterHashRef}{$$parameterNameRef}{'buildFile'} eq "yesAutoBuild") || (${$parameterHashRef}{$$parameterNameRef}{'buildFile'} eq 1) ) {  #1 for arrays

	    return "1";  #Flag that autobuild is needed
	}
	elsif ( (${$parameterHashRef}{$$parameterNameRef}{'buildFile'} eq "yesAutoDownLoad") || (${$parameterHashRef}{$$parameterNameRef}{'buildFile'} eq 1) ) {
	    
	    if (${$parameterHashRef}{$$parameterNameRef}{'default'} eq ${$scriptParameterHashRef}{$$parameterNameRef}) {
		
		return "1";  #Flag that autobuild is needed
	    }
	    else {
		
		$logger->fatal("Could not find file ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{$$parameterNameRef}, "\n");
		$logger->fatal("Make sure that file exists or use the default for this parameter to enable automatic download via Cosmid", "\n");
		exit;
	    }
	}
	else {

	    return "0";  #No autobuild is needed   
	}
    }
}


sub ParseHumanGenomeReference {

##ParseHumanGenomeReference
    
##Function : Detect version and source of the humanGenomeReference: Source (hg19 or GRCh).
##Returns  : ""
##Arguments: $fileInfoHashRef, $humanGenomeReferenceRef
##         : $fileInfoHashRef         => The fileInfo hash {REF}
##         : $humanGenomeReferenceRef => The human genome {REF}
    
    my $fileInfoHashRef = $_[0];
    my $humanGenomeReferenceRef = $_[1];
    
    if ($$humanGenomeReferenceRef =~/^Homo_sapiens.GRCh(\d+\.\d+|\d+)/) {  #Used to change capture kit genome reference version later

	${$fileInfoHashRef}{'humanGenomeReferenceVersion'} = $1;
	${$fileInfoHashRef}{'humanGenomeReferenceSource'} = "GRCh";  #Ensembl
    }
    elsif ($$humanGenomeReferenceRef =~/^Homo_sapiens.hg(\d+)/) {  #Used to change capture kit genome reference version later

	${$fileInfoHashRef}{'humanGenomeReferenceVersion'} = $1;
	${$fileInfoHashRef}{'humanGenomeReferenceSource'} = "hg";  #Refseq
    }
    else {

	$logger->warn("MIP cannot detect what kind of humanGenomeReference you have supplied. If you want to automatically set the capture kits used please supply the refrence on this format: [Species].[Source][Version].", "\n");
    }
    ${$fileInfoHashRef}{'humanGenomeReferenceNameNoEnding'} = &RemoveFileEnding(\$$humanGenomeReferenceRef, ".fasta");
    ${$fileInfoHashRef}{'humanGenomeCompressed'} = &CheckGzipped(\$$humanGenomeReferenceRef);
}


sub CheckFileEndingsToBeBuilt {

##CheckFileEndingsToBeBuilt
    
##Function : Checks files to be built by combining filename stub with fileendings.
##Returns  : "" 
##Arguments: $scriptParameterHashRef, $fileEndingsRef
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $fileEndingsRef         => Reference to the fileEndings to be added to the filename stub {REF}
##         : $parameterName          => MIP parameter name
    
    my $scriptParameterHashRef = $_[0];
    my $fileEndingsRef = $_[1];
    my $parameterName = $_[2]; 
    
    for (my $fileEndingsRefCounter=0;$fileEndingsRefCounter<scalar(@{$fileEndingsRef});$fileEndingsRefCounter++) {  #All fileEndings
	
	&CheckExistance(\%parameter, \%{$scriptParameterHashRef}, \(${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{$parameterName}.${$fileEndingsRef}[$fileEndingsRefCounter]), \$parameterName, "f");
    }
}


sub CheckExistance {

##CheckExistance
    
##Function : Checks if a file/directory exists and if autoBuild is on or not. If file/directory does not extis and there is no autobuild, croaks and exists.
##Returns  : "" 
##Arguments: $parameterHashRef, $scriptParameterHashRef, $itemNameRef, $parameterNameRef, $itemTypeToCheck, $sampleIDRef
##         : $parameterHashRef       => The parameters hash
##         : $scriptParameterHashRef => The active parameter for this analysis hash
##         : $itemNameRef            => Item to check for existance {REF}
##         : $parameterNameRef       => MIP parameter name {REF}
##         : $itemTypeToCheck        => The type of item to check
##         : $sampleIDRef            => SampleId {REF}

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $itemNameRef = $_[2];
    my $parameterNameRef = $_[3];
    my $itemTypeToCheck = $_[4];    
    my $sampleIDRef = $_[5];
    
    if ($itemTypeToCheck eq "d") {

	unless (-d $$itemNameRef) {  #Check existence of supplied directory
	    
	    $logger->fatal($USAGE, "\n");
	    $logger->fatal("Could not find intended ".$$parameterNameRef." directory: ".$$itemNameRef, "\n");
	    exit;		
	}
    }
    elsif ($itemTypeToCheck eq "f") {
	
	unless (-f $$itemNameRef) {  #Check existence of supplied file in supplied reference dir
	    
	    if (defined($sampleIDRef)) {  #Individual files per sampleID
		
		${$parameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$$sampleIDRef}{$$parameterNameRef}{'buildFile'} = &CheckAutoBuild(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \$$parameterNameRef, \$$sampleIDRef);  #Check autoBuild or not and return value
		
		if (${$parameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$$sampleIDRef}{$$parameterNameRef}{'buildFile'} == 0) {  #No autobuild
		    
		    $logger->fatal($USAGE, "\n");
		    $logger->fatal("Could not find intended ".$$parameterNameRef." file: ".$$itemNameRef, "\n");
		    exit;		
		}
	    }
	    else {
		
		${$parameterHashRef}{$$parameterNameRef}{'buildFile'} = &CheckAutoBuild(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \$$parameterNameRef);  #Check autoBuild or not and return value
		 
		if (${$parameterHashRef}{$$parameterNameRef}{'buildFile'} == 0) {  #No autobuild
		    
		    $logger->fatal($USAGE, "\n");
		    $logger->fatal("Could not find intended ".$$parameterNameRef." file: ".$$itemNameRef, "\n");
		    exit;		
		}
	    }
	}
	else {
	    
	    if (defined($sampleIDRef)) {
		
		${$parameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$$sampleIDRef}{$$parameterNameRef}{'buildFile'} = 0;  #File exist in this check
	    }
	    else {

		${$parameterHashRef}{$$parameterNameRef}{'buildFile'} =  0;  #File exist in this check
	    }
	}
    }
}


sub SetAutoBuildFeature {

##SetAutoBuildFeature
    
##Function : Sets parameters with autoBuild enabled to the new value dependent on $referenceFileNameRef
##Returns  : "" 
##Arguments: $scriptParameterHashRef, $parameterName, $referenceFileEndingRef, $referenceFileNameRef, $printSwitch
##         : $scriptParameterHashRef => The activa parameters for this analysis hash {REF}
##         : $parameterName          => MIP parameter name 
##         : $referenceFileEndingRef => Reference file name ending {REF}
##         : $referenceFileNameRef   => Reference file name {REF}
##         : $printSwitch            => To print or not

    my $scriptParameterHashRef = $_[0];
    my $parameterName = $_[1];
    my $referenceFileEndingRef = $_[2];
    my $referenceFileNameRef = $_[3];
    my $printSwitch = $_[4];
     
     if( defined(${$scriptParameterHashRef}{$parameterName}) && (${$scriptParameterHashRef}{$parameterName} eq "notSetYet") ) { 

	 ${$scriptParameterHashRef}{$parameterName} = $$referenceFileNameRef.$$referenceFileEndingRef;

	 if ( (defined($printSwitch)) && ($printSwitch ne "noPrint") ) {

	     $logger->info("Set ".$parameterName." to: ".${$scriptParameterHashRef}{$parameterName}, "\n");
	 }
	 if ($parameterName eq "bwaBuildReference") {

	     &CheckFileEndingsToBeBuilt(\%{$scriptParameterHashRef}, \@bwaBuildReferenceFileEndings, "bwaBuildReference");
	 }
	 elsif ($parameterName eq "mosaikJumpDbStub") {

	     &CheckFileEndingsToBeBuilt(\%{$scriptParameterHashRef}, \@mosaikJumpDbStubFileEndings, "mosaikJumpDbStub");
	 }
	 else {  #Complete fileName - No stubs
	    
	     &CheckExistance(\%parameter, \%{$scriptParameterHashRef}, \(${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{$parameterName}), \$parameterName, "f");
         }
    }
}


sub MoveMosaikNN {

##MoveMosaikNN

##Function : Locate MOSAIK path and move neural network files in place if lacking
##Returns  : "" 
##Arguments: $scriptParameterHashRef
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}

    my $scriptParameterHashRef = $_[0];

    my @paths = split(/:/,$ENV{PATH});  #Find Mosaik installation path

    for (my $pathsCounter=0;$pathsCounter<scalar(@paths);$pathsCounter++) {

	if ($paths[$pathsCounter] =~/MOSAIK/) {  #Select MOSAIK path
	    
	   $paths[$pathsCounter] =~ s/bin\//src\/networkFile/g;  #Location of NN files

	   $logger->warn("Could not find Mosaik Network Files in ".${$scriptParameterHashRef}{'referencesDir'},"\n");
	   $logger->info("Copying Mosaik Network Files ".${$scriptParameterHashRef}{'mosaikAlignNeuralNetworkSeFile'}." and ".${$scriptParameterHashRef}{'mosaikAlignNeuralNetworkPeFile'}." to ".${$scriptParameterHashRef}{'referencesDir'}." from ".$paths[$pathsCounter], "\n");
	   `cp $paths[$pathsCounter]/${$scriptParameterHashRef}{'mosaikAlignNeuralNetworkSeFile'} ${$scriptParameterHashRef}{'referencesDir'}/`;  #Copying files in place
	   `cp $paths[$pathsCounter]/${$scriptParameterHashRef}{'mosaikAlignNeuralNetworkPeFile'} ${$scriptParameterHashRef}{'referencesDir'}/`;  #Copying files in place
	   last;
	}
    }
}


sub CheckUserInfoArrays {

##CheckUserInfoArrays
    
##Function : Determine if the user supplied info on array parameter
##Returns  : "0|1" 
##Arguments: $scriptParameterHashRef, $arrayRef, $parameterName
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $arrayRef               => Array to loop in for parameter {REF}
##         : $parameterName          => MIP parameter to evaluate
    
    
    my $scriptParameterHashRef = $_[0];
    my $arrayRef = $_[1];
    my $parameterName = $_[2];
    
    my $userSuppliedInfoSwitch;
    
    if (scalar(@{$arrayRef}) == 0) {  #No user supplied sample info
	
	if (defined(${$scriptParameterHashRef}{$parameterName})) {  #sampleIDs info in config file
	    
	    $userSuppliedInfoSwitch = 1;  #No user supplied sample info, but present in config file do NOT overwrite using info from pedigree file
	}
	else {  #No sampleIDs info in config file
	    
	    $userSuppliedInfoSwitch = 0;  #No user supplied sample info, not defined $scriptParameter{'sampleIDs'} in config file, add it from pedigree file
	}
    }
    else {
	$userSuppliedInfoSwitch = 1;  #User supplied sample info, do NOT overwrite using info from pedigree file	
    }
    return $userSuppliedInfoSwitch;
}


sub SetTargetFiles {
    
##SetTargetFiles
    
##Function : Sets the target files and replaces constant genomic information with the info used in present analysis. Adds target file to sampleInfo for qc print later.
##Returns  : "" 
##Arguments: $scriptParameterHashRef, $sampleInfoHashRef, $supportedCaptureKitHashRef, $humanGenomeReferenceSourceRef, $humanGenomeReferenceVersionRef, $familyIDRef, $sampleIDRef, $parameterNameRef, $referenceFileEndingRef
##         : $scriptParameterHashRef         => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef              => Info on samples and family hash {REF}
##         : $supportedCaptureKitHashRef     => The supported capture kits hash {REF}
##         : $humanGenomeReferenceSourceRef  => The human genome source {REF}
##         : $humanGenomeReferenceVersionRef => The human genome build version {REF}
##         : $familyIDRef                    => Family ID {REF}
##         : $sampleIDRef                    => Sample ID  {REF}
##         : $parameterName                  => MIP parameter name
##         : $referenceFileEndingRef         => File name ending {REF}
    
    my $scriptParameterHashRef = $_[0];
    my $sampleInfoHashRef = $_[1];
    my $supportedCaptureKitHashRef = $_[2];
    my $humanGenomeReferenceSourceRef = $_[3];
    my $humanGenomeReferenceVersionRef = $_[4];
    my $familyIDRef = $_[5];
    my $sampleIDRef = $_[6];
    my $parameterNameRef = $_[7];
    my $referenceFileEndingRef = $_[8];
    
    if (defined(${$scriptParameterHashRef}{$$familyIDRef}{$$sampleIDRef}{$$parameterNameRef})) {  #Capture kit check

	${$scriptParameterHashRef}{$$familyIDRef}{$$sampleIDRef}{$$parameterNameRef} =~ s/GenomeReferenceSource/$$humanGenomeReferenceSourceRef/;  #Replace with Refseq genome or Ensembl genome
	${$scriptParameterHashRef}{$$familyIDRef}{$$sampleIDRef}{$$parameterNameRef} =~ s/Version/$$humanGenomeReferenceVersionRef/;  #Replace with actual version 
	
	${$sampleInfoHashRef}{$$familyIDRef}{$$sampleIDRef}{$$parameterNameRef} = ${$scriptParameterHashRef}{$$familyIDRef}{$$sampleIDRef}{$$parameterNameRef};  #Add to sampleInfo for qc print later
	
	&CheckExistance(\%parameter, \%{$scriptParameterHashRef}, \(${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{$$familyIDRef}{$$sampleIDRef}{$$parameterNameRef}), \$$parameterNameRef, "f", \$$sampleIDRef);

	$logger->info("Set ".$$parameterNameRef." to: ".${$scriptParameterHashRef}{$$familyIDRef}{$$sampleIDRef}{$$parameterNameRef}, "\n");
    }
    else {
	
	${$supportedCaptureKitHashRef}{'Latest'} =~ s/GenomeReferenceSource/$$humanGenomeReferenceSourceRef/;  #Replace with Refseq genome or Ensembl genome
	${$supportedCaptureKitHashRef}{'Latest'} =~ s/Version/$$humanGenomeReferenceVersionRef/;  #Replace with actual version
	
	${$scriptParameterHashRef}{$$parameterNameRef} = "notSetYet";  #Required for autobuild
	
	&SetAutoBuildFeature(\%{$scriptParameterHashRef}, $$parameterNameRef, \$$referenceFileEndingRef, \${$supportedCaptureKitHashRef}{'Latest'}, "noPrint");  #Always use the most updated capture kit when building target list
    }
}


sub PrepareArrayParameters {

##PrepareArrayParameters
    
##Function : Check if user supplied cmd info and supplies arrayParameters to scriptParameters
##Returns  : "" 
##Arguments: $parameterHashRef, $arrayRef, $orderParametersArrayRef, $broadcastsArrayRef, $parameterName, $parameterType, $parameterDefault, $associatedPrograms, $parameterExistsCheck
##         : $parameterHashRef        => The parameter hash {REF}
##         : $arrayRef                => Array to loop in for parameter {REF}
##         : $orderParametersArrayRef => Order of addition to parameter array {REF}
##         : $broadcastsArrayRef      => Holds the parameters info for broadcasting later {REF}
##         : $parameterName           => MIP parameter to evaluate
##         : $parameterType           => Type of MIP parameter 
##         : $parameterDefault        => The parameter default value
##         : $associatedPrograms      => Programs that use the parameter. Comma separated string
##         : $parameterExistsCheck    => Check if intendent file exists in reference directory

    my $parameterHashRef = $_[0];
    my $arrayRef = $_[1];
    my $orderParametersArrayRef = $_[2];
    my $broadcastsArrayRef = $_[3];
    my $parameterName = $_[4];
    my $parameterType = $_[5];
    my $parameterDefault = $_[6];
    my $associatedPrograms = $_[7]; 
    my $parameterExistsCheck = $_[8];

    if (scalar(@{$arrayRef}) == 0) {  #No input from cmd
	
	${$parameterHashRef}{$parameterName}{'value'} = "nocmdinput";  #To enable use of subroutine &AddToScriptParameter
    }
    else {

	${$parameterHashRef}{$parameterName}{'value'} = "SetbyUser";
	@{$arrayRef} = join(',',@{$arrayRef});  #If user supplied parameter a comma separated list
    }
    push(@{$orderParametersArrayRef}, $parameterName);  #Add to enable later evaluation of parameters in proper order & write to master file
    &AddToScriptParameter(\%{$parameterHashRef}, \%scriptParameter, \%sampleInfo, \%fileInfo, \%referenceFileEndings, \@{$broadcastsArrayRef},
			  {'parameterName' => $parameterName,
			   'parameterValue' => ${$parameterHashRef}{$parameterName}{'value'},
			   'parameterType' => $parameterType,
			   'parameterDefault' => $parameterDefault,
			   'associatedPrograms' => $associatedPrograms,
			   'parameterExistsCheck' => $parameterExistsCheck
			  });
}


sub CheckUniqueTargetFiles {

##CheckUniqueTargetFiles
    
##Function : Checks target files within parameters for identical entries
##Returns  : "" 
##Arguments: $parameterHashRef, $scriptParameterHashRef, $arrayRef, $countRef, $fileToCompareRef, $parameterName
##         : $parameterHashRef       => The parameter hash {REF}
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $arrayRef               => Array to loop in for parameter (e.g. sampleID) {REF}
##         : $countRef               => Offset in array {REF}
##         : $fileToCompareRef       => The file to compare against rest of array {REF}
##         : $parameterName          => MIP parameter to evaluate
    
    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $arrayRef = $_[2]; 
    my $countRef = $_[3]; 
    my $fileToCompareRef = $_[4];
    my $parameterName = $_[5];
    
    for (my $compareCounter=($$countRef + 1);$compareCounter<scalar(@{$arrayRef});$compareCounter++) {  #Compare all target files to remove autoBuild if file names are identical for each flag
	
	if ( (defined($$fileToCompareRef)) && (defined(${$scriptParameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$arrayRef}[$compareCounter] }{$parameterName})) ) {

	    if ($$fileToCompareRef eq ${$scriptParameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$arrayRef}[$compareCounter] }{$parameterName}) {  #Identical target files
	    
		${$parameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$arrayRef}[$compareCounter] }{$parameterName}{'buildFile'} = 0;  #Turn off autoBuild
	    }
	}
    }
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
##Arguments: $scriptParameterHashRef, $familyIDRef, $sampleIDRef, $parameterName
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $familyIDRef            => Family ID {REF}
##         : $sampleIDRef            => Sample ID  {REF}
##         : $parameterName          => MIP parameter name

    my $scriptParameterHashRef = $_[0];
    my $familyIDRef = $_[1];
    my $sampleIDRef = $_[2];
    my $parameterName = $_[3];

    if (defined(${$scriptParameterHashRef}{$$familyIDRef}{$$sampleIDRef}{$parameterName})) {
	
	${$scriptParameterHashRef}{$parameterName} = 1;  #Define in scriptParameter so that we now that parameter is present per SampleID
    }
}


sub EnableArrayParameter {

##EnableArrayParameter
    
##Function : Adds arrayRef to scriptParameters for recreation of cmd in log and seperated input parameter string into array elements
##Returns  : ""
##Arguments: $scriptParameterHashRef, $arrayRef, $parameterNameRef
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $arrayRef               => Array to loop through {REF}
##         : $parameterNameRef       => MIP parameter name {REF}

    my $scriptParameterHashRef = $_[0];
    my $arrayRef = $_[1];
    my $parameterNameRef = $_[2];
    
    ${$scriptParameterHashRef}{$$parameterNameRef} = join(',',@{$arrayRef});  #Add to enable recreation of cmd line later
    @{$arrayRef} = split(/,/,join(',', @{$arrayRef}));  #Enables comma separated list of sample IDs from user supplied cmd info
}


sub SetTargetandAutoBuild {

##SetTargetandAutoBuild
    
##Function : Set autoBuild for target files and calls SetTargetFile subroutine
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $fileInfoHashRef, $arrayRef, $parameterNameRef, $fileEndingRef
##         : $parameterHashRef       => The parameter hash {REF}
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $fileInfoHashRef        => The fileInfo hash {REF}
##         : $arrayRef               => Array to loop through {REF}
##         : $parameterNameRef       => MIP parameter name {REF}
##         : $fileEndingRef          => File ending {REF}

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $fileInfoHashRef = $_[2];
    my $arrayRef = $_[3];
    my $parameterNameRef = $_[4];
    my $fileEndingRef = $_[5];
    
    for (my $elementsCounter=0;$elementsCounter<scalar(@{$arrayRef});$elementsCounter++) {

	${$parameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{${$arrayRef}[$elementsCounter]}{$$parameterNameRef}{'buildFile'} = "yesAutoBuild";
	&SetTargetFiles(\%{$scriptParameterHashRef}, \%sampleInfo, \%supportedCaptureKit, \${$fileInfoHashRef}{'humanGenomeReferenceSource'}, \${$fileInfoHashRef}{'humanGenomeReferenceVersion'}, \${$scriptParameterHashRef}{'familyID'}, \${$arrayRef}[$elementsCounter], \$$parameterNameRef, \$$fileEndingRef);
    }   
}


sub CheckTargetExistFileBed {

##CheckTargetExistFileBed
    
##Function : Check that supplied target file ends with ".bed" and exists.
##Returns  : ""
##Arguments: $scriptParameterHashRef, $fileRef, $parameterName
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $fileRef                => File to check for existance and file ending {REF}
##         : $parameterName          => MIP parameter name

    my $scriptParameterHashRef = $_[0];
    my $fileRef = $_[1];
    my $parameterName = $_[2];

    if (defined($$fileRef)) {

	if ($$fileRef !~/.bed$/) {
	
	$logger->fatal("Could not find intendended 'file ending with .bed' for target file: ".$$fileRef." in parameter '-".$parameterName."'", "\n");
	exit;
	}
	unless (-f ${$scriptParameterHashRef}{'referencesDir'}."/".$$fileRef) {

	    $logger->fatal("Could not find intendended '.bed' file for target file: ".${$scriptParameterHashRef}{'referencesDir'}."/".$$fileRef." in parameter '-".$parameterName."'", "\n");
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
	
	$logger->fatal("The number of supplied '-".$parameterNameQuery."' (=".scalar(@{$arrayQueryRef}).") do not equal the number of '-".$parameterName."' (=".scalar(@{$arrayRef})."). Please specify a equal number of elements for both parameters", "\n");
	exit;
    }
}

sub SetAutoBuildAndScriptParameterPerSample {

##SetAutoBuildAndScriptParameterPerSample
    
##Function : Sets autoBuild and populates scriptParameter hash with array elements per sampleID.
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $sampleIDArrayRef, $parameterArrayRef, $parameterNameRef
##         : $parameterHashRef             => The parameter hash {REF}
##         : $scriptParameterHashRef       => The active parameters for this analysis hash {REF}
##         : $sampleIDArrayRef             => SampleID array {REF}
##         : $parameterArrayRef            => Parameter array {REF}
##         : $parameterNameRef             => MIP parameter name {REF}

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $sampleIDArrayRef = $_[2];
    my $parameterArrayRef = $_[3];
    my $parameterNameRef = $_[4];

    for (my $sampleIDsCounter=0;$sampleIDsCounter<scalar(@{$sampleIDArrayRef});$sampleIDsCounter++) {  #All sampleIDs
	
	${$parameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{${$sampleIDArrayRef}[$sampleIDsCounter]}{$$parameterNameRef}{'buildFile'} = "yesAutoBuild";  #Turn on autoBuild
	${$scriptParameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$sampleIDArrayRef}[$sampleIDsCounter] }{$$parameterNameRef} = ${$parameterArrayRef}[$sampleIDsCounter];  #Populate hash that is used in modules
    }
}


sub SetTargetFileGeneralBuildParameter {

##SetTargetFileGeneralBuildParameter 
    
##Function : Sets the general build parameters $$sampleIDBuildFileRef and $$sampleIDBuildFileNoEndingRef and sets buildfile key to "0".
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $referenceFileEndingsHashRef, $targetfileRef, $parameterName, $sampleIDBuildFileRef, $sampleIDBuildFileNoEndingRef, $sampleIDRef
##         : $parameterHashRef             => The parameter hash {REF}
##         : $scriptParameterHashRef       => The active parameters for this analysis hash {REF}
##         : $referenceFileEndingsHashRef  => The associated reference file endings {REF}
##         : $targetfileRef                => Final file {REF}
##         : $parameterName                => MIP parameter
##         : $sampleIDBuildFileRef         => File that will be created {REF}
##         : $sampleIDBuildFileNoEndingRef => File that will be created with file ending removed {REF}
##         : $sampleIDRef                  => SampleID {REF}

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $referenceFileEndingsHashRef = $_[2];
    my $targetfileRef = $_[3];
    my $parameterName = $_[4];
    my $sampleIDBuildFileRef = $_[5];
    my $sampleIDBuildFileNoEndingRef = $_[6];
    my $sampleIDRef = $_[7];
    
    $$sampleIDBuildFileNoEndingRef = &RemoveFileEnding(\$$targetfileRef, ${$referenceFileEndingsHashRef}{$parameterName});  #Remove ".fileending" from reference filename
    ${$parameterHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{$$sampleIDRef}{$parameterName}{'buildFile'} = 0;  #Build once then done
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
    
    print $FILEHANDLE "[ -s ".$$intendedFilePathRef." ] ";  #Check file exists and is larger than 0
    print $FILEHANDLE "&& rm ".$$temporaryFilePathRef." ";  #If other processes already has created file, remove temp file
    print $FILEHANDLE "|| ";  #File has not been created by other processes
    print $FILEHANDLE "mv ".$$temporaryFilePathRef." ".$$intendedFilePathRef,"\n\n";  #Move file in place
    
}


sub DefineSnpEffFiles {

##DefineSnpEffFiles

##Function : Defines and adds snpEff/snpSift files and features to hash
##Returns  : ""
##Arguments: $parameterHashRef
##         : $parameterHashRef        => The parameter hash {REF}

    my $parameterHashRef = $_[0];

    my %snpEffFile;
    my @snpSiftDownloadableFiles = ("dbsnp_138.b37.excluding_sites_after_129.vcf", "dbsnp_138.b37.vcf", "1000G_phase1.indels.b37.vcf", "1000G_phase1.snps.high_confidence.b37.vcf");

    foreach my $file (@snpSiftDownloadableFiles) {

	$snpEffFile{'snpSift'}{$file} = {  #Files that are downloadable via Cosmid
	    'downloadable' => "yes",
	};
	${$parameterHashRef}{$file}{'associatedProgram'} = "pSnpEff";
	${$parameterHashRef}{$file}{'buildFile'} = "yesAutoBuild";  #Allow autoDownLoad, but yesAutoBuild is set since the file is its own default, so no extra check is required (compared with yesAutoDownLoad)
    }
    return %snpEffFile;
}

sub DefineAnnovarTables {
    
##DefineAnnovarTables
    
##Function : Defines and adds annovar tables parameters to hash
##Returns  : ""
##Arguments: $parameterHashRef, $annovarGenomeBuildVersionRef
##         : $parameterHashRef             => The parameter hash {REF}
##         : $annovarGenomeBuildVersionRef => The current annovar genome build

    my $parameterHashRef = $_[0];
    my $annovarGenomeBuildVersionRef = $_[1];

    my %annovarTable;  #Holds all annovar tables parameters and features

    ##Define annovar tables
    my @annovarTablesGeneAnno = ("refGene", "knownGene", "ensGene");  #Tables using annotation option "geneanno"
    my @annovarTablesRegionAnno = ("mce46way", "gerp++elem", "segdup", "tfbs", "mirna");  #Tables using annotation option "regionanno"
    my @annovarTablesFilter = ("snp137", "snp135", "snp132", "snp131", "snp130", "snp129", "snp137NonFlagged", "snp135NonFlagged", "snp132NonFlagged", "snp131NonFlagged", "snp130NonFlagged", "1000g2012apr_all", "1000g2012apr_amr", "1000g2012apr_eur", "1000g2012apr_asn", "1000g2012apr_afr", "1000g2012feb_all", "esp6500si_all", "esp6500_all", "esp6500_aa", "esp6500_ea", "esp5400_all", "esp5400_aa", "esp5400_ea","clinvar_20131105", "ljb2_sift", "ljb2_pp2hdiv", "ljb2_pp2hvar", "ljb2_mt", "ljb2_ma", "ljb2_fathmm", "ljb2_siphy", "ljb2_lrt", "ljb_all", "ljb2_gerp++", "ljb2_phylop", "caddgt20", "caddgt10");  #Tables using annotation option "filter"
    my @annovarTablesUrlUcsc = ("mce46way", "segdup", "tfbs", "mirna");  #Tables using urlAlias "ucsc"
    my @annovarGenericFiltering = ("esp6500si_all", "esp6500_all", "esp6500_aa", "esp6500_ea", "esp5400_all", "esp5400_aa", "esp5400_ea","clinvar_20131105");  #Tables using generic option
    my @annovarGenericFiles = ($$annovarGenomeBuildVersionRef."_esp6500si_all.txt", $$annovarGenomeBuildVersionRef."_esp6500_all.txt", $$annovarGenomeBuildVersionRef."_esp6500_aa.txt", $$annovarGenomeBuildVersionRef."_esp6500_ea.txt", $$annovarGenomeBuildVersionRef."_esp5400_all.txt", $$annovarGenomeBuildVersionRef."_esp5400_aa.txt", $$annovarGenomeBuildVersionRef."_esp5400_ea.txt", $$annovarGenomeBuildVersionRef."_clinvar_20131105.txt");  #Generic table files
    my @annovarRefgeneFiles = ($$annovarGenomeBuildVersionRef."_refGene.txt", $$annovarGenomeBuildVersionRef."_refGeneMrna.fa", $$annovarGenomeBuildVersionRef."_refLink.txt", "GRCh37_MT_ensGene.txt", "GRCh37_MT_ensGeneMrna.fa");  #Cater for multiple download
    my @annovarKnownGeneFiles = ($$annovarGenomeBuildVersionRef."_knownGene.txt", $$annovarGenomeBuildVersionRef."_kgXref.txt", $$annovarGenomeBuildVersionRef."_knownGeneMrna.fa");  #Cater for multiple download
    my @annovarEnsGeneFiles = ($$annovarGenomeBuildVersionRef."_ensGene.txt", $$annovarGenomeBuildVersionRef."_ensGeneMrna.fa", "GRCh37_MT_ensGene.txt", "GRCh37_MT_ensGeneMrna.fa");  #Cater for multiple download

    #Set UCSC alias for download from UCSC
    $annovarTable{'mce46way'}{'ucscAlias'} = "phastConsElements46way";
    $annovarTable{'segdup'}{'ucscAlias'} = "genomicSuperDups";
    $annovarTable{'tfbs'}{'ucscAlias'} = "tfbsConsSites";
    $annovarTable{'mirna'}{'ucscAlias'} = "wgRna";

    #Set GeneAnno files
    push(@{$annovarTable{'refGene'}{'file'}}, @annovarRefgeneFiles);
    push(@{$annovarTable{'knownGene'}{'file'}}, @annovarKnownGeneFiles);
    push(@{$annovarTable{'ensGene'}{'file'}}, @annovarEnsGeneFiles); 

    for (my $tablesCounter=0;$tablesCounter<scalar(@annovarSupportedTableNames);$tablesCounter++) {
	
	&AnnovarTableParameters(\%annovarTable, \$annovarSupportedTableNames[$tablesCounter], \@annovarSupportedTableNames, "dbtype", $annovarSupportedTableNames[$tablesCounter]);
	&AnnovarTableParameters(\%annovarTable, \$annovarSupportedTableNames[$tablesCounter], \@annovarSupportedTableNames, "download", $annovarSupportedTableNames[$tablesCounter]);
	${$parameterHashRef}{$annovarSupportedTableNames[$tablesCounter]}{'buildFile'} = "yesAutoBuild";  #Allow autobuild
    }


    #Tables using different download call from dbtype call
    $annovarTable{'1000g2012apr_all'}{'download'} = "ALL.sites.2012_04";
    $annovarTable{'1000g2012feb_all'}{'download'} = "ALL.sites.2012_02";
    $annovarTable{'1000g2012apr_afr'}{'download'} = "AFR.sites.2012_04";
    $annovarTable{'1000g2012apr_amr'}{'download'} = "AMR.sites.2012_04";
    $annovarTable{'1000g2012apr_eur'}{'download'} = "EUR.sites.2012_04";
    $annovarTable{'1000g2012apr_asn'}{'download'} = "ASN.sites.2012_04";

    #Set 1000G Table filename
    $annovarTable{'1000g2012apr_all'}{'file'}[0] = $$annovarGenomeBuildVersionRef."_ALL.sites.2012_04.txt";
    $annovarTable{'1000g2012feb_all'}{'file'}[0] = $$annovarGenomeBuildVersionRef."_ALL.sites.2012_02.txt";
    $annovarTable{'1000g2012apr_afr'}{'file'}[0] = $$annovarGenomeBuildVersionRef."_AFR.sites.2012_04.txt";
    $annovarTable{'1000g2012apr_amr'}{'file'}[0] = $$annovarGenomeBuildVersionRef."_AMR.sites.2012_04.txt";
    $annovarTable{'1000g2012apr_eur'}{'file'}[0] = $$annovarGenomeBuildVersionRef."_EUR.sites.2012_04.txt";
    $annovarTable{'1000g2012apr_asn'}{'file'}[0] = $$annovarGenomeBuildVersionRef."_ASN.sites.2012_04.txt";

    for (my $tablesCounter=0;$tablesCounter<scalar(@annovarTablesGeneAnno);$tablesCounter++) {

	&AnnovarTableParameters(\%annovarTable, \$annovarTablesGeneAnno[$tablesCounter], \@annovarTablesGeneAnno, "annotation", "geneanno");
	&AnnovarTableParameters(\%annovarTable, \$annovarTablesGeneAnno[$tablesCounter], \@annovarTablesGeneAnno, "urlAlias", "annovar");
    }
    for (my $tablesCounter=0;$tablesCounter<scalar(@annovarTablesRegionAnno);$tablesCounter++) {

	&AnnovarTableParameters(\%annovarTable, \$annovarTablesRegionAnno[$tablesCounter], \@annovarTablesRegionAnno, "annotation", "regionanno");
	&AnnovarTableParameters(\%annovarTable, \$annovarTablesRegionAnno[$tablesCounter], \@annovarTablesRegionAnno, "urlAlias", "annovar");
	&AnnovarTableParameters(\%annovarTable, \$annovarTablesRegionAnno[$tablesCounter], \@annovarTablesUrlUcsc, "urlAlias", "ucsc");  #Replace for ucsc tables NOTE: not all in RegionAnno
    }
    for (my $tablesCounter=0;$tablesCounter<scalar(@annovarTablesFilter);$tablesCounter++) {

	&AnnovarTableParameters(\%annovarTable, \$annovarTablesFilter[$tablesCounter], \@annovarTablesFilter, "annotation", "filter");
	&AnnovarTableParameters(\%annovarTable, \$annovarTablesFilter[$tablesCounter], \@annovarTablesFilter, "urlAlias", "annovar");
	&AnnovarTableParameters(\%annovarTable, \$annovarTablesFilter[$tablesCounter], \@annovarTablesFilter, "indexFile", ".idx"); 
    }	
    for (my $tablesCounter=0;$tablesCounter<scalar(@annovarGenericFiltering);$tablesCounter++) {
	
	&AnnovarTableParameters(\%annovarTable, \$annovarGenericFiltering[$tablesCounter], \@annovarGenericFiltering, "dbtype", "generic");
	&AnnovarTableParameters(\%annovarTable, \$annovarGenericFiltering[$tablesCounter], \@annovarGenericFiltering, "file", $annovarGenericFiles[$tablesCounter]);
    }
    return %annovarTable;
}


sub AnnovarTableParameters {
   
##AnnovarTableParameters
    
##Function : Populates annovarTable hash by looping through array and adding table with identified membership and associated parameters into hash.
##           Parameters of type "file" are added as a hash of array since multiple files can be downloaded

##Returns  : ""
##Arguments: $annovarTableHashRef, $tableNameRef, $arrayRef, $parameterType, $parameterValue
##         : $annovarTableHashRef => annovarTableHashRef {REF}
##         : $tableNameRef         => Annovar table name {REF}
##         : $arrayRef             => Array to search for membership {REF}
##         : $parameterType        => Type of table parameter
##         : $parameterValue       => Parameter value
    
    my $annovarTableHashRef = $_[0];
    my $tableNameRef = $_[1];
    my $arrayRef = $_[2];
    my $parameterType = $_[3];
    my $parameterValue = $_[4];
    
    for (my $tablesCounter=0;$tablesCounter<scalar(@{$arrayRef});$tablesCounter++) {
    
	if (${$arrayRef}[$tablesCounter] eq $$tableNameRef) {  #Membership test
	    
	    if ($parameterType eq "file") {  #Add as array instead, since some annovar tables have multiple files downloaded for the same call
		
		push(@{${$annovarTableHashRef}{$$tableNameRef}{$parameterType}}, $parameterValue);  #Add as array to hashRef
	    }
	    else {
		
		${$annovarTableHashRef}{$$tableNameRef}{$parameterType} = $parameterValue;
	    }
	    last;  #No table should be represented twice within the same array
	}
    }
}


sub CollectSeqContigs {

##CollectSeqContigs
    
##Function : Collects sequences contigs used in analysis from human genome sequence dictionnary associated with $humanGenomeReference
##Returns  : ""
##Arguments: $contigsArrayRef, $referencesDirRef, $humanGenomeReferenceNameNoEndingRef
##         : $contigsArrayRef                         => Contig array {REF}
##         : $referencesDirRef                        => The MIP reference directory
##         : $humanGenomeReferenceNameNoEndingRef     => The associated human genome file without file ending                        

    my $contigsArrayRef = $_[0];
    my $referencesDirRef = $_[1];
    my $humanGenomeReferenceNameNoEndingRef = $_[2];

    my $pqSeqDict = q?perl -nae 'if($F[0]=~/^\@SQ/) { if($F[1]=~/SN\:(\S+)/) {print $1, ",";} }' ?; 
    my $SeqDictLocation = $$referencesDirRef."/".$$humanGenomeReferenceNameNoEndingRef.".dict";
    @{$contigsArrayRef} = `$pqSeqDict $SeqDictLocation `;  #Returns a comma seperated string of sequence contigs from dict file
    @{$contigsArrayRef} = split(/,/,join(',', @{$contigsArrayRef}));
}


sub ReplaceConfigParamWithCMDInfo {

##ReplaceConfigParamWithCMDInfo
    
##Function : Replace config parameter with cmd info for active parameter
##Returns  : Path to Cosmid Resource directory for current analysis
##Arguments: $parameterHashRef, $scriptParameterHashRef, $parameter, $parameterName
##         : $parameterHashRef       => The parameter hash {REF}
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $parameter              => The parameter hash {REF}
##         : $parameterName          => MIP parameter name

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $parameterName = $_[2];

    if (${$parameterHashRef}{$parameterName}{'value'} ne "nocmdinput") {  #Replace config parameter with cmd info for parameter
	
	${$scriptParameterHashRef}{$parameterName} = ${$parameterHashRef}{$parameterName}{'value'};  #Transfer to active parameter
    }
}


sub DefineSupportedCosmidReferences {

##DefineSupportedCosmidReferences
    
##Function : Defines the Cosmid manager hash keys and populates it from arguments 
##Returns  : ""
##Arguments: $supportedCosmidReferenceHashRef, $parameterName, $cosmidResourceName, $CosmidResourceVersion, 
##         : $supportedCosmidReferenceHashRef => The supported cosmid references hash {REF} 
##         : $parameterName                    => MIP parameter name
##         : $cosmidResourceName               => Cosmid Resource name
##         : $CosmidResourceVersion            => Version of the cosmid Resource to download
##         : $humanGenomeReferenceVersionRef   => The human genome build used in the analysis
##         : $compressedSwitch                 => If files after download are compressed or not

    my $supportedCosmidReferenceHashRef = $_[0];
    my $parameterName = $_[1];
    my $cosmidResourceName = $_[2];
    my $CosmidResourceVersion = $_[3];
    my $humanGenomeReferenceVersionRef = $_[4];
    my $compressedSwitch = $_[5];

    ${$supportedCosmidReferenceHashRef}{$parameterName} = {

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
##Arguments: $scriptParameterHashRef
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}

    my $scriptParameterHashRef = $_[0];

    my %cosmidResources;  #Hash to load cosmid info to

    if (-f ${$scriptParameterHashRef}{'referencesDir'}."/cosmid.yaml") {  #Cosmid.yaml file exists in reference directory
	
	%cosmidResources = &LoadYAML(\%scriptParameter, ${$scriptParameterHashRef}{'referencesDir'}."/cosmid.yaml");  #Load yaml file
	
	unless (defined($cosmidResources{'directory'})) {  #Set Directory entry if not defined
	    
	    $cosmidResources{'directory'} = ${$scriptParameterHashRef}{'referencesDir'}."/resources";  #Set the Cosmid default directory
	}
    }
    else {  #No cosmid.yaml exist in reference directory
	
	$cosmidResources{'directory'} = ${$scriptParameterHashRef}{'referencesDir'}."/resources";  #Set the Cosmid default directory
    }
    return $cosmidResources{'directory'};
}


sub AdjustNrCoresToSeqMode {

##AdjustNrCoresToSeqMode
    
##Function : Adjust the number of cores to be used in the analysis according to sequencing mode requirements.
##Returns  : ""
##Arguments: $nrCoresRef, $sequenceRunTypeRef
##         : $nrCoresRef         => The maximum number of cores to be use before printing "wait" statement {REF}
##         : $sequenceRunTypeRef => Type of sequencing [Paired-end|Single-end] {REF}

    my $nrCoresRef = $_[0];
    my $sequenceRunTypeRef = $_[1];
    
    if ($$sequenceRunTypeRef eq "Paired-end") {  #Second read direction if present

	$$nrCoresRef =  $$nrCoresRef + 2;  #2 processes per file
    }
    else {  #Single-end

	$$nrCoresRef = $$nrCoresRef + 1;  #Only 1 file and one process
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
    
    if ($$counterRef == $$coreCounterRef * $$nrCoresRef) {  #Using only nr of cores eq to lanes or maximumCores
	
	print $FILEHANDLE "wait", "\n\n";
	$$coreCounterRef=$$coreCounterRef+1;  #Increase the maximum number of cores allowed to be used since "wait" was just printed 
    }    
}


sub CheckBuildDownLoadPreRequisites {

##CheckBuildDownLoadPreRequisites
    
##Function : Checks if some of the active prerequisites needs to be downloaded and calls subroutine to download them if required
##Returns  : ""
##Arguments: $parameterHashRef, $supportedCosmidReferenceHashRef, $programName
##         : $parameterHashRef                => The parameter hash {REF}
##         : $scriptParameterHashRef          => The active parameters for this analysis hash {REF}
##         : $supportedCosmidReferenceHashRef => The supported cosmid references hash {REF}
##         : $programName                     => Active program

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $supportedCosmidReferenceHashRef = $_[2];
    my $programName = $_[3];
    
    for my $parameterName (keys %{$supportedCosmidReferenceHashRef}) {  #Supported cosmid references for MIP parameters
	
	if ( (defined(${$parameterHashRef}{$parameterName}{'associatedProgram'})) && (${$parameterHashRef}{$parameterName}{'associatedProgram'} =~/$programName/) ) {  #If the cosmid supported parameter is associated with the MIP program
	    
	    if ( (${$scriptParameterHashRef}{"p".$programName} == 1) && (${$scriptParameterHashRef}{'dryRunAll'} != 1) ) {  #Only enable autoDownload for active programs

		if (${$parameterHashRef}{$parameterName}{'buildFile'} eq 1) {  #Enable autoBuild
		    
		    &BuildDownLoadablePreRequisites(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, \%{$supportedCosmidReferenceHashRef}, ${$scriptParameterHashRef}{'familyID'}, ${$scriptParameterHashRef}{'aligner'}, $programName);
		    last;  #Perform once
		}
	    }
	}
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
	
	$logger->fatal("The supplied file: ".$$fileNameRef." for parameter '".$$parameterNameRef."' does not have the supported file ending '".$$fileEndingRef."'.", "\n");
	exit;
    }
}


sub PrintSupportedAnnovarTableNames {

##PrintSupportedAnnovarTableNames
    
##Function : Print the by MIP supported Annovar Table names to STDOUT and exists
##Returns  : ""
##Arguments: $scriptParameterHashRef, $annovarSupportedTableNames
##         : $scriptParameterHashRef     => The active parameters for this analysis
##         : $annovarSupportedTableNames => The supported Annovar tables
    
    my $scriptParameterHashRef = $_[0];
    my $annovarSupportedTableNamesArrayRef = $_[1];

    if (${$scriptParameterHashRef}{'logFile'}) {

	$logger->info("These Annovar databases are supported by MIP:\n");

	foreach my $annovarSupportedTableName (@{$annovarSupportedTableNamesArrayRef}) {
	    
	    $logger->info($annovarSupportedTableName, "\n");
	}
    }
    else {

	print STDOUT "These Annovar databases are supported by MIP:\n";

	foreach my $annovarSupportedTableName (@{$annovarSupportedTableNamesArrayRef}) {
	    
	    print STDOUT $annovarSupportedTableName, "\n";
	}
    }
    exit;
}


sub CheckEntryHashofArray {

##CheckEntryHashofArray
    
##Function : Test element for being part of hash of array at supplied key. 
##Returns  : Return "1" if element is not part of array
##Arguments: $hashRef, $keyRef, $elementRef
##         : $hashRef    => Hash {REF} 
##         : $keyRef     => The key pointing to the array in the $hashRef {REF}
##         : $elementRef => Element to look for in hash of array {REF}

    my $hashRef = $_[0];
    my $keyRef = $_[1];
    my $elementRef = $_[2];

    if (defined($$hashRef{$$keyRef})) {  #Information on entry present

	if ( ! ( grep /$$elementRef/, @{$$hashRef{$$keyRef}} ) ) {  #If element is not part of array

	    return 1;
	}
    }
}


sub CheckMostCompleteAndRemoveFile {
    
##CheckMostCompleteAndRemoveFile
	
##Function  : Checks if the file is recorded as the "MostCompleteBAM|VCF". If false writes removal of file(s) to supplied filehandle
##Returns   : ""
##Arguments : $FILEHANDLE, $mostCompleteRef, $fileRef, $fileEnding
##          : $FILEHANDLE      => SBATCH script FILEHANDLE to print to 
##          : $mostCompleteRef => The mostComplete file (BAM|VCF) {REF}
##          : $fileRef         => Current file {REF}
##          : $fileEnding      => File ending of $fileRef
    
    my $FILEHANDLE = $_[0];
    my $mostCompleteRef = $_[1];
    my $fileRef = $_[2];
    my $fileEnding = $_[3];
    
    if ( (defined($$mostCompleteRef)) && (defined($$fileRef)) ) {  #Not to disturb first dry_run of analysis

	unless ($$mostCompleteRef eq $$fileRef) {  #Do not remove mostCompleteBAM|VCF
		
	    my $fileName = &RemoveFileEnding(\$$fileRef, $fileEnding);

	    if (defined($fileName)) {  #Successfully removed file ending using &RemoveFileEnding
	    
		my $end = ".*";  #Remove all files with ending with ".*"

		if ($fileEnding eq ".bam") {  #For BAM files 
		
		    $end = ".ba*";  #Removes both .bam and .bai
		}
		if ($fileEnding eq ".vcf") {  #For VCF files
		
		    $end = ".vcf*";  #Removes both .vcf and .vcf.idx
		}
		##Print removal of file to sbatch script 
		print $FILEHANDLE "rm ";
		print $FILEHANDLE $fileName.$end, "\n\n";  #Remove file(s)
	    }
	}
    }
}


sub ConcatenateVCFs {

##ConcatenateVCFs
    
##Function : Concatenate VCFs
##Returns  : ""
##Arguments: $scriptParameterHashRef, $FILEHANDLE, $arrayRef, $infilePrefix, $infilePostfix, $outfile, $reorderSwith
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $FILEHANDLE             => SBATCH script FILEHANDLE to print to
##         : $arrayRef               => Holding the number and part of file names to be combined
##         : $infilePrefix           => Will be combined with the each array element
##         : $infilePostfix          => Will be combined with the each array element
##         : $outfile                => The combined outfile
##         : $reorderSwith           => Reorder header
    
    my $scriptParameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];
    my $arrayRef = $_[2];
    my $infilePrefix = $_[3];
    my $infilePostfix = $_[4];
    my $outfile = $_[5];
    my $reorderSwith = $_[6];
    
    unless (defined($infilePostfix)) {
	
	$infilePostfix = "";  #No postfix
    }
    unless (defined($outfile)) {
	
	$outfile = $infilePrefix.".vcf";  
    }
    print $FILEHANDLE "## SNPSift Split -j","\n";

    ## Writes java core commands to filehandle.
    &JavaCore($FILEHANDLE,
	      {'memoryAllocation' => "Xmx4g",
	       'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
	       'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
	      });
    
    print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'snpEffPath'}."/SnpSift.jar ";
    print $FILEHANDLE "split -j ";  #Joinf VCFs together
    
    for (my $elementCounter=0;$elementCounter<scalar(@{$arrayRef});$elementCounter++) {
	
	print $FILEHANDLE $infilePrefix.${$arrayRef}[$elementCounter].$infilePostfix." ";  #files to combined
    }
    if ( (defined($_[6])) && $reorderSwith eq "reOrderHeader") {

	print $FILEHANDLE "| ";  #Pipe
	print $FILEHANDLE "perl ".${$scriptParameterHashRef}{'inScriptDir'}."/vcfParser.pl ";  #Parses the vcf output	
    }
    
    print $FILEHANDLE "> ".$outfile, "\n\n";  #OutFile
}


sub CombineVariants {

##CombineVariants
    
##Function : Writes sbatch code to supplied filehandle to combine variants in vcf format. Each array element is combined with the infilePre and Postfix.
##Returns  : ""
##Arguments: $scriptParameterHashRef, $FILEHANDLE, $arrayRef, $infilePrefix, $infilePostfix, $outfile
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $FILEHANDLE             => SBATCH script FILEHANDLE to print to
##         : $arrayRef               => Holding the number and part of file names to be combined
##         : $infilePrefix           => Will be combined with the each array element
##         : $infilePostfix          => Will be combined with the each array element
##         : $outfile                => The combined outfile

    my $scriptParameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];
    my $arrayRef = $_[2];
    my $infilePrefix = $_[3];
    my $infilePostfix = $_[4];
    my $outfile = $_[5];
    
    unless (defined($infilePostfix)) {
	
	$infilePostfix = "";  #No postfix
    }
    unless (defined($outfile)) {
	
	$outfile = $infilePrefix.".vcf";  
    }
    print $FILEHANDLE "## GATK CombineVariants","\n";

    ## Writes java core commands to filehandle.
    &JavaCore($FILEHANDLE,
	      {'memoryAllocation' => "Xmx4g",
	       'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
	       'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
	      });

    print $FILEHANDLE "-jar ".${$scriptParameterHashRef}{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar ";
    print $FILEHANDLE "-l INFO ";  #Set the minimum level of logging
    print $FILEHANDLE "-T CombineVariants ";  #Type of analysis to run
    print $FILEHANDLE "-R ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference file
    
    for (my $elementCounter=0;$elementCounter<scalar(@{$arrayRef});$elementCounter++) {
	
	print $FILEHANDLE "-V: ".$infilePrefix.${$arrayRef}[$elementCounter].$infilePostfix." ";  #files to combined
    }
    print $FILEHANDLE "-genotypeMergeOptions UNSORTED ";  #Take the genotypes in any order. Should be fine since the same order and number of samples exists in all files

    print $FILEHANDLE "-o ".$outfile, "\n\n";  #OutFile
}


sub ConcatenateVariants {

##ConcatenateVariants
    
##Function : Writes sbatch code to supplied filehandle to concatenate variants in vcf format. Each array element is combined with the infilePre and Postfix.
##Returns  : ""
##Arguments: $scriptParameterHashRef, $FILEHANDLE, $arrayRef, $infilePrefix, $infilePostfix, $outfile
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $FILEHANDLE             => SBATCH script FILEHANDLE to print to
##         : $arrayRef               => Holding the number and part of file names to be combined
##         : $infilePrefix           => Will be combined with the each array element
##         : $infilePostfix          => Will be combined with the each array element
##         : $outfile                => The combined outfile

    my $scriptParameterHashRef = $_[0];
    my $FILEHANDLE = $_[1];
    my $arrayRef = $_[2];
    my $infilePrefix = $_[3];
    my $infilePostfix = $_[4];
    my $outfile = $_[5];
    
    unless (defined($infilePostfix)) {
	
	$infilePostfix = "";  #No postfix
    }
    unless (defined($outfile)) {
	
	$outfile = $infilePrefix.".vcf";  
    }
    print $FILEHANDLE "## GATK CatVariants","\n";

    ## Writes java core commands to filehandle.
    &JavaCore($FILEHANDLE,
	      {'memoryAllocation' => "Xmx4g",
	       'javaUseLargePagesRef' => \${$scriptParameterHashRef}{'javaUseLargePages'},
	       'javaTemporaryDirectory' => ${$scriptParameterHashRef}{'tempDirectory'}
	      });

    print $FILEHANDLE "-cp ".${$scriptParameterHashRef}{'genomeAnalysisToolKitPath'}."/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants ";  #Type of analysis to run
    print $FILEHANDLE "-l INFO ";  #Set the minimum level of logging
    print $FILEHANDLE "-R ".${$scriptParameterHashRef}{'referencesDir'}."/".${$scriptParameterHashRef}{'humanGenomeReference'}." ";  #Reference file
    print $FILEHANDLE "-assumeSorted ";  #assumeSorted should be true if the input files are already sorted

    for (my $elementCounter=0;$elementCounter<scalar(@{$arrayRef});$elementCounter++) {
	
	print $FILEHANDLE "-V: ".$infilePrefix.${$arrayRef}[$elementCounter].$infilePostfix." ";  #files to combined
    }
    print $FILEHANDLE "-out ".$outfile, "\n\n";  #OutFile
}


sub ChanjoConvert {
    
##ChanjoConvert
    
##Function : Writes sbatch code to supplied filehandle to convert a database file into a sorted bed file and piping it
##Returns  : ""
##Arguments: $FILEHANDLE, $chanjoDatabase
##         : $FILEHANDLE => SBATCH script FILEHANDLE to print to
##         : $chanjoDatabase => The database to convert and pipe
    
    my $FILEHANDLE = $_[0];
    my $chanjoDatabase = $_[1];
    
    print $FILEHANDLE "sort ";  #Unix sort
    print $FILEHANDLE "-k1,1 -k2,2n ";  #Sort numerically 
    print $FILEHANDLE $chanjoDatabase." ";  #Chanjo Db file to be sorted
    print $FILEHANDLE "| ";  #Pipe
    print $FILEHANDLE "chanjo convert ";  #Convert ccds to bed
    print $FILEHANDLE "| ";  #Pipe
    
}


sub RemovePedigreeElements {

##RemovePedigreeElements
    
##Function : Removes ALL keys at third level except 'Capture_kit'. 
##Returns  : ""
##Arguments: $hashRef
##         : $hashRef => Hash {REF}
    
    my $hashRef = $_[0];
    
    for my $familyID (keys %{$hashRef}) {
	
	for my $sampleID (keys %{ ${$hashRef}{$familyID} })  {
	    
	    for my $pedigreeElements (keys %{ ${$hashRef}{$familyID}{$sampleID} })  {
		
		unless ($pedigreeElements eq 'Capture_kit') {
		    
		    delete(${$hashRef}{$familyID}{$sampleID}{$pedigreeElements})
		}
		
	    }
	}
    }
}


sub WriteUseLargePages {

##WriteUseLargePages
    
##Function : Write useLargePages java flag to sbatch script if enabled
##Returns  : ""
##Arguments: $FILEHANDLE, $arrayRef, $infilePrefix, $infilePostfix, $outfile
##         : $FILEHANDLE       => SBATCH script FILEHANDLE to print to
##         : $useLargePagesRef => UseLargePages for requiring large memory pages (cross-platform flag) {REF}
	
    my $FILEHANDLE = $_[0];
    my $useLargePagesRef = $_[1];
    
    if ($$useLargePagesRef ne "no") {
	
	print $FILEHANDLE "-XX:-UseLargePages ";  #UseLargePages for requiring large memory pages (cross-platform flag)
    }
}


sub CreateLog4perlCongfig {

##CreateLog4perlCongfig
    
##Function : Create log4perl config file. 
##Returns  : "$config"
##Arguments: $fileName
##         : $fileName => log4perl config file {REF}

    my $fileNameRef = $_[0];

    my $conf = q?
        log4perl.category.MIPLogger = TRACE, LogFile, ScreenApp
        log4perl.appender.LogFile = Log::Log4perl::Appender::File
        log4perl.appender.LogFile.filename = ?.$$fileNameRef.q?
        log4perl.appender.LogFile.layout=PatternLayout
        log4perl.appender.LogFile.layout.ConversionPattern = [%p] %d %c - %m%n

        log4perl.appender.ScreenApp = Log::Log4perl::Appender::Screen
        log4perl.appender.ScreenApp.layout = PatternLayout
        log4perl.appender.ScreenApp.layout.ConversionPattern = [%p] %d %c - %m%n
        ?;
    return $conf;
}


sub DeafultLog4perlFile {

##DeafultLog4perlFile
    
##Function : Set the default Log4perl file using supplied dynamic parameters.
##Returns  : "$LogFile"
##Arguments: $scriptParameterHashRef, $cmdInputRef, $scriptRef, $dateRef, $dateTimeStampRef
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $cmdInputRef            => User supplied info on cmd for logFile option {REF}
##         : $scriptRef              => The script that is executed {REF}
##         : $dateRef                => The date {REF}
##         : $dateTimeStampRef       => The date and time {REF}
   
    my $scriptParameterHashRef = $_[0];
    my $cmdInputRef = $_[1];
    my $scriptRef = $_[2];
    my $dateRef = $_[3];
    my $dateTimeStampRef = $_[4];
    
    if ($$cmdInputRef eq "nocmdinput") {  #No input from cmd i.e. do not create default logging directory or set default

	`mkdir -p ${$scriptParameterHashRef}{'outDataDir'}/${$scriptParameterHashRef}{'familyID'}/mip_log/$$dateRef;`;  #Creates the default log dir
	my $LogFile = ${$scriptParameterHashRef}{'outDataDir'}."/".${$scriptParameterHashRef}{'familyID'}."/mip_log/".$$dateRef."/".$$scriptRef."_".$$dateTimeStampRef.".log";  #concatenates log filename	
	return $LogFile;
    }
}


sub CheckHumanGenomeFileEndings {

##CheckHumanGenomeFileEndings
    
##Function : Check the existance of associated Human genome files.
##Returns  : ""
##Arguments: $parameterHashRef, $referencesDirRef
##         : $parameterHashRef                        => The parameter hash {REF}
##         : $referencesDirRef                        => The MIP reference directory {REF}
##         : $humanGenomeReferenceFileEndingsArrayRef => The associated human genome file endings {REF}
##         : $humanGenomeReferenceNameNoEndingRef     => The associated human genome file without file ending {REF}
##         : $parameterNameRef                        => The parameter under evaluation {REF}                   

    my $parameterHashRef = $_[0];
    my $referencesDirRef = $_[1];
    my $humanGenomeReferenceFileEndingsArrayRef = $_[2];
    my $humanGenomeReferenceNameNoEndingRef = $_[3];
    my $parameterNameRef = $_[4];
    
##Enable autoBuild of metafiles 	       
    ${$parameterHashRef}{$$parameterNameRef.".dict"}{'buildFile'} = "yesAutoBuild";
    ${$parameterHashRef}{$$parameterNameRef.".fasta.fai"}{'buildFile'} = "yesAutoBuild";

    for (my $fileEndingsCounter=0;$fileEndingsCounter<scalar(@{$humanGenomeReferenceFileEndingsArrayRef});$fileEndingsCounter++) {
	
	my $intendedFilePathRef = \($$referencesDirRef."/".$$humanGenomeReferenceNameNoEndingRef.${$humanGenomeReferenceFileEndingsArrayRef}[$fileEndingsCounter]);
    
	&CheckExistance(\%{$parameterHashRef}, \%scriptParameter, $intendedFilePathRef, \($$parameterNameRef.${$humanGenomeReferenceFileEndingsArrayRef}[$fileEndingsCounter]), "f");
    }
    if (${$parameterHashRef}{$$parameterNameRef.".dict"}{'buildFile'} eq 0) {
	
	##Collect sequence contigs from human reference ".dict" file since it exists
	&CollectSeqContigs(\@contigs, \$$referencesDirRef, \$$humanGenomeReferenceNameNoEndingRef);  #Preparation for future changes but not active yet
    }
}


sub CheckMergePicardToolsMergeSamFilesPrevious {
    
##CheckMergePicardToolsMergeSamFilesPrevious
    
##Function : Checks if previous alignments have been supplied for each sampleID. Saves merge info in sampleInfo hash.
##Returns  : "" 
##Arguments: $scriptParameterHashRef, $sampleInfoHashRef
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}
##         : $sampleInfoHashRef      => Info on samples and family hash {REF}
    
    my $scriptParameterHashRef = $_[0];
    my $sampleInfoHashRef = $_[1];
    
    for (my $sampleIDCounter=0;$sampleIDCounter<scalar(@{${$scriptParameterHashRef}{'sampleIDs'}});$sampleIDCounter++) {  #Check all samples to check, which are to be merged with previous files later
	
	if (scalar(@{${$scriptParameterHashRef}{'picardToolsMergeSamFilesPrevious'}}) > 0) {  #Supplied info - check for which sampleID(s)  	
	    
	    for (my $mergeFileCounter=0;$mergeFileCounter<scalar(@{${$scriptParameterHashRef}{'picardToolsMergeSamFilesPrevious'}});$mergeFileCounter++) {
		
		if (${$scriptParameterHashRef}{'picardToolsMergeSamFilesPrevious'}[$mergeFileCounter] =~ /${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter]/) {  #Look for sampleID in previously generated file to be merged with current run to be able to merge correct files
		    
		    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }{'picardToolsMergeSamFilesPrevious'} = 1;
		}
		else {
		    
		    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }{'picardToolsMergeSamFilesPrevious'} = 0;
		}
	    }
	}
	else {  #Not supplied - Set to 0 
	    
	    ${$sampleInfoHashRef}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$scriptParameterHashRef}{'sampleIDs'}[$sampleIDCounter] }{'picardToolsMergeSamFilesPrevious'} = 0;
	}
    }
}


sub CreateFamFile {

##CreateFamFile

##Function : Create .fam file to be used in variant calling analyses.
##Returns  : "" 
##Arguments: $scriptParameterHashRef
##         : $scriptParameterHashRef => The active parameters for this analysis hash {REF}

    my $scriptParameterHashRef = $_[0];

    my $pqFamFile = q?perl -nae 'print $F[0], "\t", $F[1], "\t", $F[2], "\t", $F[3], "\t", $F[4], "\t", $F[5], "\n";'?;
    my $famFile = ${$scriptParameterHashRef}{'outDataDir'}."/".${$scriptParameterHashRef}{'familyID'}."/".${$scriptParameterHashRef}{'familyID'}.".fam";
    `$pqFamFile ${$scriptParameterHashRef}{'pedigreeFile'} > $famFile;`;
}


sub CheckAnnovarTables {

##CheckAnnovarTables
    
##Function : Checks supported Annovar Tables and that the supplied supported ones exists 
##Returns  : ""
##Arguments: $parameterHashRef, $scriptParameterHashRef, $annovarTableHashRef
##         : $parameterHashRef       => Holds all parameters
##         : $scriptParameterHashRef => Holds all set parameter for analysis
##         : $annovarTableHashRef    => annovarTableHashRef {REF}

    my $parameterHashRef = $_[0];
    my $scriptParameterHashRef = $_[1];
    my $annovarTableHashRef = $_[2];
    
    my $intendedFilePathRef;
    
    for (my $tableNamesCounter=0;$tableNamesCounter<scalar(@{${$scriptParameterHashRef}{'annovarTableNames'}});$tableNamesCounter++) {  #All AnnovarTables    
	
	if (defined(${$annovarTableHashRef}{${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter]})) {  #Supported Annovar database
	    
	    if (defined(${$annovarTableHashRef}{ ${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter] }{'file'})) {
		
		for (my $filesCounter=0;$filesCounter<scalar(@{${$annovarTableHashRef}{ ${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter] }{'file'}});$filesCounter++) {  #All annovarTable file(s), some tables have multiple files downloaded from the same call
		    
		    $intendedFilePathRef = \(${$scriptParameterHashRef}{'annovarPath'}."/humandb/".${$annovarTableHashRef}{ ${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter] }{'file'}[$filesCounter]);
		    &CheckExistance(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, $intendedFilePathRef, \${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter], "f");
		}
	    }
	    elsif (defined(${$annovarTableHashRef}{ ${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter] }{'ucscAlias'})){
		
		$intendedFilePathRef = \(${$scriptParameterHashRef}{'annovarPath'}."/humandb/".${$scriptParameterHashRef}{'annovarGenomeBuildVersion'}."_".${$annovarTableHashRef}{ ${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter] }{'ucscAlias'}.".txt");
		&CheckExistance(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, $intendedFilePathRef, \${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter], "f");
	    }
	    else {
		
		$intendedFilePathRef = \(${$scriptParameterHashRef}{'annovarPath'}."/humandb/".${$scriptParameterHashRef}{'annovarGenomeBuildVersion'}."_".${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter].".txt");
		&CheckExistance(\%{$parameterHashRef}, \%{$scriptParameterHashRef}, $intendedFilePathRef, \${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter], "f");
	    }
	}
	else {  #Annovar Table not supported by MIP
	    
	    $logger->error("You supplied Annovar database: ".${$scriptParameterHashRef}{'annovarTableNames'}[$tableNamesCounter]." which is not supported by MIP. MIP can only process supported annovar databases\n");
	    &PrintSupportedAnnovarTableNames(\%{$scriptParameterHashRef}, \@annovarSupportedTableNames);
	}
    }
}


sub CollectOutDataPathsEntries {
    
##CollectOutDataPathsEntries
    
##Function : Collects all programs outfile path(s) created by MIP as OutDirectory->value and outFile->value located in %sampleInfo.  
##Returns  : ""
##Arguments: $sampleInfoHashRef, $pathsArrayRef
##         : $sampleInfoHashRef => Info on samples and family hash {REF}
##         : $pathsArrayRef     => Holds the collected paths {REF}
    
    my $sampleInfoHashRef = $_[0];
    my $pathsArrayRef = $_[1];
    
    for my $familyID ( keys %{$sampleInfoHashRef} ) {  #For every family id 
	
	for my $member ( keys %{ ${$sampleInfoHashRef}{$familyID} }) {  #For every familyID and sampleID
	    
	    if (${$sampleInfoHashRef}{$familyID}{$member}{'program'}) {  #Only examine programs     
		
		for my $program ( keys %{ ${$sampleInfoHashRef}{$familyID}{$member}{'program'} } ) {  #For every programs           
		    
		    my @outDirectoryArray;  #Temporary array for collecting outDirectories within the same program
		    my @outFileArray;  #Temporary array for collecting outFile within the same program
		    
		    for my $key ( keys %{ ${$sampleInfoHashRef}{$familyID}{$member}{'program'}{$program} } ) { #For every key within program
			
			## Check if KeyName is "OutDirectory" or "OutFile"  and adds to @pathsArrayRef if true.
			&CollectOutFile(\@{$pathsArrayRef}, \@outDirectoryArray, \@outFileArray, ${$sampleInfoHashRef}{$familyID}{$member}{'program'}{$program}{$key}, $key);
			
			if (ref(${$sampleInfoHashRef}{$familyID}{$member}{'program'}{$program}{$key}) eq "HASH" ) { #HASH reference indicating more levels
			    
			    for my $secondKey ( keys %{ ${$sampleInfoHashRef}{$familyID}{$member}{'program'}{$program}{$key} } ) { #For every programs
				
				## Check if KeyName is "OutDirectory" or "OutFile"  and adds to @pathsArrayRef if true.
				&CollectOutFile(\@{$pathsArrayRef}, \@outDirectoryArray, \@outFileArray, ${$sampleInfoHashRef}{$familyID}{$member}{'program'}{$program}{$key}{$secondKey}, $secondKey);
			    }
			}
		    }
		}
	    }
	}
    }
}


sub CollectOutFile {
    
##CollectOutFile
    
##Function  : Check if KeyName is "OutDirectory" or "OutFile"  and adds to @pathsArrayRef if true.
##Returns   : ""
##Arguments: $pathsArrayRef, $outDirectoryArrayRef, $outFileArrayRef
##         : $pathsArrayRef        => Holds the collected paths {REF}
##         : $outDirectoryArrayRef => Holds temporary outDirectory path(s) {REF}
##         : $outFileArrayRef      => Holds temporary outDirectory path(s) {REF}
##         : $key                  => The hash key
##         : $keyName              => The actual key  
    
    my $pathsArrayRef = $_[0];
    my $outDirectoryArrayRef = $_[1];
    my $outFileArrayRef = $_[2];
    my $key = $_[3];
    my $keyName = $_[4];	
    
    if ($keyName eq "OutDirectory") {
	
	push(@{$outDirectoryArrayRef}, $key);
    }
    if ($keyName eq "OutFile") {
	
	push(@{$outFileArrayRef}, $key);
    }
    if (scalar(@{$outDirectoryArrayRef}) == 1 && (scalar(@{$outFileArrayRef}) == 1) ) {  #Both outDirectory and outFile have been collected, time to join
	
	push(@{$pathsArrayRef}, ${$outDirectoryArrayRef}[0]."/".${$outFileArrayRef}[0]);
	@{$outDirectoryArrayRef} = ();  #Restart
	@{$outFileArrayRef} = ();  #Restart
    }
}


sub CollectPathEntries {
    
##CollectPathEntries
    
##Function  : Collects all programs outfile path(s) created by MIP as Path->value located in %sampleInfo.
##Returns   : ""
##Arguments: $sampleInfoHashRef, $pathsArrayRef
##         : $sampleInfoHashRef => Info on samples and family hash {REF}
##         : $pathsArrayRef     => Holds the collected paths {REF}
    
    my $sampleInfoHashRef = $_[0];
    my $pathsArrayRef = $_[1];
    
    for my $familyID ( keys %{$sampleInfoHashRef} ) {  #For every familyID 
	
	for my $member ( keys %{ ${$sampleInfoHashRef}{$familyID} }) {  #For every familyID and sampleID     
	    
	    for my $key ( keys %{ ${$sampleInfoHashRef}{$familyID}{$member} } ) {  #For every key within member
		
		## Check if KeyName is "PATH" and adds to @pathsArrayRef if true.
		&CheckAndAddToArray(\@{$pathsArrayRef}, ${$sampleInfoHashRef}{$familyID}{$member}{$key}, $key);
		
		if (ref(${$sampleInfoHashRef}{$familyID}{$member}{$key}) eq "HASH" ) {   #HASH reference indicating more levels
		    
		    for my $secondKey ( keys %{ ${$sampleInfoHashRef}{$familyID}{$member}{$key} } ) { #For every secondkey with program
			
			## Check if KeyName is "PATH" and adds to @pathsArrayRef if true.
			&CheckAndAddToArray(\@{$pathsArrayRef}, ${$sampleInfoHashRef}{$familyID}{$member}{$key}{$secondKey}, $secondKey);
		    }
		}
	    }
	}
    } 
}


sub CheckAndAddToArray {
    
##CheckAndAddToArray
    
##Function  : Check if KeyName is "PATH" and adds to @pathsArrayRef if true.
##Returns   : ""
##Arguments: $pathsArrayRef, $key, $keyName
##         : $pathsArrayRef => Holds the collected paths {REF}
##         : $key           => The hash key
##         : $keyName       => The actual key
    
    my $pathsArrayRef = $_[0];
    my $key = $_[1];
    my $keyName = $_[2];
    
    if ($keyName eq "Path") {
	
	push(@{$pathsArrayRef}, $key);
    }
}

sub MigrateFilesToTemp {

##MigrateFilesToTemp
    
##Function : Copies files from source to temporary folder. Loop over files specified by $arrayRef and collects files from $extractArrayRef.
##Returns  : ""
##Arguments: $scriptParameterHashRef, $sampleInfoHashRef, $arrayRef, $extractArrayRef, $FILEHANDLE, $argHashRef
##         : $scriptParameterHashRef               => The active parameters for this analysis hash {REF}
##         : $arrayRef                             => The array of files to copy
##         : $extractArrayRef                      => The array to extract files from
##         : $FILEHANDLE                           => Filehandle to write to
##         : $argHashRef{'$sampleInfoHashRef'}     => Info on samples and family hash {REF}
##         : $argHashRef{'inSampleDirectory'}      => The directory for the file to be copied
##         : $argHashRef{'nrCores'}                => The number of cores that can be used
##         : $argHashRef{'fileEnding'}             => File ending. Set to "" to not add any file ending or omit from call {Optional, unless further arguments are given.}
##         : $argHashRef{'sampleID'}               => the sampleID {Optional}

    my $scriptParameterHashRef = shift;
    my $arrayRef = shift;
    my $extractArrayRef = shift;
    my $FILEHANDLE = shift;
    my ($argHashRef) = @_;

    my $pairedEndTracker = 0;
    my $coreCounter=1;

    print $FILEHANDLE "## Copying file(s) to temporary directory\n";
    for (my $fileCounter=0;$fileCounter<scalar( @{$arrayRef});$fileCounter++) {  #For all files
	
	my $sequenceRunMode;

	if ( (defined(${$argHashRef}{'sampleInfoHashRef'})) && (defined(${$argHashRef}{'sampleID'})) ) {
	
	    $sequenceRunMode = ${$argHashRef}{'sampleInfoHashRef'}{ ${$scriptParameterHashRef}{'familyID'} }{ ${$argHashRef}{'sampleID'} }{'file'}{ ${$arrayRef}[$fileCounter] }{'sequenceRunType'};  #Collect paired-end or single-end sequence run mode
	}
	&PrintWait(\$fileCounter, \${$argHashRef}{'nrCores'}, \$coreCounter, $FILEHANDLE);
	
	## Copies file to temporary directory.
	&MigrateFileToTemp($FILEHANDLE,
			    {'path' => ${$argHashRef}{'inSampleDirectory'}."/".${$extractArrayRef}[$pairedEndTracker],
			     'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'},
			     'fileEnding' => ${$argHashRef}{'fileEnding'}
			    });

	if ( (defined($sequenceRunMode)) && ($sequenceRunMode eq "Paired-end") ) {
	    
	    $pairedEndTracker = $pairedEndTracker+1;  #Increment to collect correct read 2 from %infile
	    &MigrateFileToTemp($FILEHANDLE,
				{'path' => ${$argHashRef}{'inSampleDirectory'}."/".${$extractArrayRef}[$pairedEndTracker],
				 'tempDirectory' => ${$scriptParameterHashRef}{'tempDirectory'},
				 'fileEnding' => ${$argHashRef}{'fileEnding'}
				});
	}
	$pairedEndTracker++;  #Increment to correctly track both single-end runs and paired-end runs
    }
    print $FILEHANDLE "wait", "\n\n";
}


sub MigrateFilesFromTemp {

##MigrateFilesFromTemp
    
##Function : Copies files from temporary folder to source. Loop over files specified by $arrayRef and collects files from $extractArrayRef.
##Returns  : ""
##Arguments: $arrayRef, $extractArrayRef, $outSampleDirectory, $tempDirectory, $nrCores, $scriptParameterHashRef, $sampleInfoHashRef, $sampleID, $FILEHANDLE
##         : $arrayRef           => The array of files to copy
##         : $extractArrayRef    => The array to extract files from
##         : $outSampleDirectory => The directory for the file to be copied
##         : $tempDirectory      => The node directory to copy to
##         : $fileEnding         => The fileending to use for the outfile
##         : $nrCores            => The number of cores that can be used
##         : $FILEHANDLE         => Filehandle to write to

    my $arrayRef = $_[0];
    my $extractArrayRef = $_[1];
    my $outSampleDirectory = $_[2];
    my $tempDirectory = $_[3];
    my $nrCores = $_[4];
    my $fileEnding = $_[5];
    my $FILEHANDLE = $_[6];

    my $coreCounter=1;

    print $FILEHANDLE "## Copying file(s) from node\n";
    for (my $fileCounter=0;$fileCounter<scalar( @{$arrayRef});$fileCounter++) {  #For all files
	
	&PrintWait(\$fileCounter, \$nrCores, \$coreCounter, $FILEHANDLE);
	
	## Copies file from temporary directory.
	&MigrateFileFromTemp($tempDirectory."/".${$extractArrayRef}[$fileCounter].$fileEnding, $outSampleDirectory, $FILEHANDLE);
    }
    print $FILEHANDLE "wait", "\n\n";

}


sub MigrateFileToTemp {

##MigrateFileToTemp
    
##Function : Copies file to temporary directory.  
##Returns  : "$fileName"
##Arguments: $FILEHANDLE, $argHashRef
##         : $FILEHANDLE                  => Filehandle to write to
##         : $argHashRef{'path'}          => The infile path
##         : $argHashRef{'tempDirectory'} => The node directory to copy to
##         : $argHashRef{'fileEnding'}    => File ending {Optional}

    my $FILEHANDLE = shift;
    my ($argHashRef) = @_;

    my $tempFilePath;

    ## Split relative path to file(s)
    my ($pathVolume, $pathDirectories, $pathFileName) = File::Spec->splitpath(${$argHashRef}{'path'});

    if (defined(${$argHashRef}{'fileEnding'})) {

	${$argHashRef}{'path'} .= ${$argHashRef}{'fileEnding'};  #Add fileEnding if supplied
    }
    print $FILEHANDLE "cp ";  #copy
    print $FILEHANDLE ${$argHashRef}{'path'}." ";  #Infile
    print $FILEHANDLE ${$argHashRef}{'tempDirectory'}, " & \n\n";  #Temp file
    return $pathFileName; 
}


sub MigrateFileFromTemp {

##MigrateFileFromTemp
    
##Function : Copies file from temporary directory.  
##Returns  : ""
##Arguments: $tempPath, $filePath, $FILEHANDLE
##         : $tempPath   => The node temp file path
##         : $filePath   => The node directory to copy to
##         : $FILEHANDLE => Filehandle to write to

    my $tempPath = $_[0];
    my $filePath = $_[1];
    my $FILEHANDLE = $_[2];

    print $FILEHANDLE "cp ";
    print $FILEHANDLE $tempPath." ";  #Infile
    print $FILEHANDLE $filePath, " & \n\n";  #Local temp file
}


sub RemoveDirectory {

##RemoveDirectory
    
##Function : Writes command to removes directory to filehandle.  
##Returns  : ""
##Arguments: $directoryRef, $FILEHANDLE
##         : $directoryRef => the directory to remove
##         : $FILEHANDLE   => Filehandle to write to

    my $tempDirectoryRef = $_[0];
    my $FILEHANDLE = $_[1];

    print $FILEHANDLE "## Remove directory\n";
    print $FILEHANDLE "rm ";  #Remove
    print $FILEHANDLE "-rf ";  #Directory
    print $FILEHANDLE $$tempDirectoryRef, "\n";  #Directory to remove
}



sub JavaCore {

##JavaCore
    
##Function : Writes java core commands to filehandle.  
##Returns  : ""
##Arguments: $FILEHANDLE, $memoryAllocation, $javaUseLargePagesRef, $javaTemporaryDirectory
##         : $FILEHANDLE             => Filehandle to write to
##         : $memoryAllocation       => Memory allocation for java 
##         : $javaUseLargePagesRef   => Use java large pages {REF}
##         : $javaTemporaryDirectory => Redirect tmp files to java temp {Optional}

    my $FILEHANDLE = shift;
    my ($argHashRef) = @_;

    #my $memoryAllocation = $_[0];
    #my $javaUseLargePagesRef = $_[1];
    #my $FILEHANDLE = $_[2];
    #my $javaTemporaryDirectory = $_[3];
    
    print $FILEHANDLE "java ";
    print $FILEHANDLE "-".${$argHashRef}{'memoryAllocation'}." "; 
    
    &WriteUseLargePages($FILEHANDLE, \${$argHashRef}{'javaUseLargePagesRef'});

    if (defined(${$argHashRef}{'javaTemporaryDirectory'})) {

	print $FILEHANDLE "-Djava.io.tmpdir=".${$argHashRef}{'javaTemporaryDirectory'}." ";  #Temporary Directory
    }
}


sub CheckEmailAddress { 
    
##CheckEmailAddress
    
##Function : Check the syntax of the email adress is valid not that it is actually exists.  
##Returns  : ""
##Arguments: $emailRef
##         : $emailRef => The email adress

    my $emailRef = $_[0];

    $$emailRef =~ /[ |\t|\r|\n]*\"?([^\"]+\"?@[^ <>\t]+\.[^ <>\t][^ <>\t]+)[ |\t|\r|\n]*/;

    unless (defined($1)) {
	
	$logger->fatal("The supplied email: ".$$emailRef." seem to be malformed. ", "\n");
	exit;
    }
}


sub BreakString {

##BreakString
    
##Function : Breaks the string supplied on format key1:value1_value2,keyN:valueN_valueN,..n . Add key to %fileInfoHashRef and values as array. This enables association of values to key supplied in config or cmd.   
##Returns  : ""
##Arguments: $fileInfoHashRef, $parameterValueRef, $parameterName, $associatedProgram
##         : $fileInfoHashRef   => The fileInfo hash {REF}
##         : $parameterValueRef => MIP parameter value {REF}
##         : $parameterName     => MIP parameter name {REF}
##         : $associatedProgram => The parameters program {REF}

    my $fileInfoHashRef = $_[0];
    my $parameterValueRef = $_[1];
    my $parameterNameRef = $_[2];
    my $associatedProgramRef = $_[3];
    
    ## Break string into key value pairs
    my @fileKeys = split(',', join(',', $$parameterValueRef));

    foreach my $element (@fileKeys) {
	
	my @tempArray = split(/:/, $element);
	@{${$fileInfoHashRef}{$$associatedProgramRef}{$$parameterNameRef}{$tempArray[0]}} = split('_', $tempArray[1]);  #Save infoKey associated with fileName
    }
}

sub AddCaptureKit {

##AddCaptureKit
    
##Function : Adds a capture kit to the scriptParameterHash. If arg->{userSuppliedParameterswitchRef} is set, go a head and add capture kit no matter what the switch was.
##Returns  : "Set capture kit or ''"
##Arguments: $referenceFileEndingsHashRef, $supportedCaptureKitHashRef, $argHashRef
##         : $referenceFileEndingsHashRef               => The associated reference file endings {REF}
##         : $supportedCaptureKitHashRef                => The supported capture kits hash {REF}
##         : $argHashRef{'captureKit'}                  => The capture kit to add
##         : $argHashRef{'parameterName'}               => The parameter name
##         : $argHashRef{'userSuppliedParameterswitch'} => Has user supplied parameter {OPTIONAL}
    
    my $referenceFileEndingsHashRef = shift;
    my $supportedCaptureKitHashRef = shift;
    my ($argHashRef) = @_;
    
    unless (defined(${$argHashRef}{'userSuppliedParameterswitch'})) {  #No detect supplied capture kit
	
	if ( defined(${$supportedCaptureKitHashRef}{ ${$argHashRef}{'captureKit'} }) ) {  #Supported capture kit alias
	    
	    return  ${$supportedCaptureKitHashRef}{ ${$argHashRef}{'captureKit'} }.${$referenceFileEndingsHashRef}{ ${$argHashRef}{'parameterName'} };
	}
	else {  #Return unchanged capture_kit string
	    
	    return ${$argHashRef}{'captureKit'}.${$referenceFileEndingsHashRef}{ ${$argHashRef}{'parameterName'} };
	}
    }
    if ( (defined(${$argHashRef}{'userSuppliedParameterswitch'})) && (${$argHashRef}{'userSuppliedParameterswitch'} == 0) ) {  #Only add if user supplied no info on parameter
	
	if ( defined(${$supportedCaptureKitHashRef}{ ${$argHashRef}{'captureKit'} }) ) {  #Supported capture kit alias
	    
	    return  ${$supportedCaptureKitHashRef}{ ${$argHashRef}{'captureKit'} }.${$referenceFileEndingsHashRef}{ ${$argHashRef}{'parameterName'} };
	} 
	else {  #Return unchanged capture_kit string
	    
	    return ${$argHashRef}{'captureKit'}.${$referenceFileEndingsHashRef}{ ${$argHashRef}{'parameterName'} };
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

####
#Decommissioned
####

sub CheckTemplateFilesPaths {

##CheckTemplateFilesPaths
    
##Function : Checks that file paths in template files exist
##Returns  : "" 
##Arguments: $fileNameRef, $parameterName
##         : $fileNameRef   => File name {REF}
##         : $parameterName => MIP parameter name

    my $fileNameRef = $_[0];
    my $parameterName = $_[1];

    open(my $TF, "<", $$fileNameRef) or $logger->logdie("Can't open '".$$fileNameRef."':".$!."\n");  

    while (<TF>) {

	chomp $_;

	if (m/^\s+$/) {	 # Avoid blank lines
            next;
        }
	if (m/^\#/) {  # Avoid "#"
            next;
        }	
	if ($_ =~/(\S+)/) {	
	    
	    my $filePath = $_;

	    if ($filePath=~/^(RD!)/) {  #intersectCollect file
		
		my @filePath = split('\t', $filePath);
		$filePath[0] =~ s/^RD!/$scriptParameter{'referencesDir'}/g;

		&CheckExistance(\%parameter, \%scriptParameter, \$filePath[0], \$parameterName, "f");  #Only check paths that pointing to reference directory
	    }
	    if ($parameterName eq "GATKHaploTypeCallerRefBAMInfile") {  #Only Paths should be present i.e. check all lines
	       
		&CheckExistance(\%parameter, \%scriptParameter, \$filePath, \$parameterName, "f");
	    }
	}
    }
    close(TF);
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
