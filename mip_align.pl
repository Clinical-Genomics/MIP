#!/usr/bin/perl -w

use strict;
use warnings;

#Master script for analysing paired end reads from the Illumina plattform in fastq(.gz) format to sorted, dedupped and merged bam files. The program performs QC, aligns reads using Mosaik or BWA and generates a coverage report.
 
#Copyright 2011 Henrik Stranneheim

=head1 SYNOPSIS
    
mip_align.pl  -id [inFilesDir,.,.,.,n] -ids [inScriptDir,.,.,.,n] -rd [reference dir] -a [project ID] -s [sampleIDs...n] -em [e-mail] -ods [outScriptDir] -odf [outDataDir] -f [familyID] -p[program]
    
=head2 COMMANDS AND OPTIONS

-ifd/--inFilesDir Infile directory(s). Comma sep (Mandatory: Supply whole path)

-isd/--inScriptDir The pipeline script in directory (Mandatory: Supply whole path)

-rd/--referencesDir Reference(s) directory (Mandatory: Supply whole path)

-p/--projectID The project ID (Mandatory)

-s/--sampleIDs The sample ID(s) (Mandatory)

-em/--email

-odd/--outDataDir The data files output directory (Mandatory: Supply whole path)

-osd/--outScriptDir The script files (.sh) output directory (Mandatory: Supply whole path)

-f/--familyID Group id of samples to be compared (defaults to "0" (=no), (Ex: 1 for IDN 1-1-1A))

-pedigree/--pedigreeFile (Supply whole path, defaults to "")

-huref/--humanGenomeReference Fasta file for the human genome reference (defaults to "")

-pFQC/--pFastQC Sequence quality analysis using FastQC (defaults to "1" (=yes))

-pMoB/--pMosaikBuild Convert reads to Mosaik format using MosaikBuild (defaults to "1" (=yes))

-mobmfl/--mosaikBuildMedianFragLength Flag for setting the mean fragment length, mfl, (defaults to (=375) bp)

-pMoA/--pMosaikAlign Align reads using MosaikAlign (defaults to "1" (=yes))

-moaref/--mosaikAlignReference MosaikAlign reference (defaults to "")

-moaannpe/--mosaikAlignNeuralNetworkPeFile MosaikAlign Neural Network PE File (defaults to "")

-moaannse/--mosaikAlignNeuralNetworkSeFile MosaikAlign Neural Network SE File (defaults to "")

-mojdb/--mosaikJumpDbStub MosaikJump stub (defaults to "")

-pBWA_aln/--pBwaAln Index reads using BWA Aln (defaults to "0" (=no))

-bwaalnq/--bwaAlnQualityTrimming BWA Aln quality threshold for read trimming (defaults to "20")

-pBWA_sampe/--pBwaSampe Align reads using BWA Sampe (defaults to "0" (=no))

-pSamT_sort/--pSamToolsSort Sort & index aligned reads using SamTools sort & index (defaults to "1" (=yes))

-pPicT_merge/--pPicardToolsMerge Merge (BAM file(s)) using PicardTools MergeSamFiles (defaults to "0" (=no))

-picT_mergeprev/--picardToolsMergePreviousFiles Flag running picardTools MergeSamFiles on merged current files and previous BAM-file(s) (Supply whole path and name, name must contain sample id, and lanes_Xn info)

-pPicT_MarkDup/--pPicardToolsMarkduplicates Markduplicates using PicardTools MarkDuplicates (defaults to "1" (=yes))

-pic_path/--picardToolsPath Path to PicardTools. Mandatory for use of PicardTools (defaults to "")

-pCC/--pCalculateCoverage Use coverage calculation tools: qaCompute, genomeCoverageBED and PicardTools (MultipleMetrics & HSmetrics) (defaults to "1" (=yes))

-pCC_Bedgc/--pGenomeCoverageBED Genome coverage calculation using genomeCoverageBED under '-pCC' (defaults to "1" (=yes))

-pCC_Bedc/--pCoverageBED BED file coverage calculation using coverageBED under '-pCC' (defaults to "1" (=yes))

-extb/--exomeTargetBed Target BED file of exome capture for coverageBed '-pCC_Bedc' (defaults to "")

-pCC_Qac/--pQaCompute Genome coverage calculation using qaCompute under '-pCC' (defaults to "1" (=yes))

-xcov/--xCoverage Max coverage depth when using '-genomeCoverageBED', '-qaCompute' (defaults to "30")

-pCC_PicMM/--pPicardToolsCollectMultipleMetrics Metrics calculation using PicardTools collectMultipleMetrics under '-pCC' (defaults to "1" (=yes))

-pCCE_PicHS/--pPicardToolsCalculateHSMetrics Capture calculation using PicardTools CalculateHSmetrics under '-pCC' (defaults to "1" (=yes))

-extbl/--exomeTargetBedInfileList Prepared target BED file for PicardTools CalculateHSMetrics (defaults to "". File ending should be ".infile_list")
              
-extpbl/--exomeTargetPaddedBedInfileList Prepared padded target BED file for PicardTools CalculateHSMetrics (defaults to "". File should be ".padXXX.infile_list")

-pRCP/--pRCovPlots Plots of genome coverage using rCovPlots (defaults to "1" (=yes))

-pGZ/--pGZip GZip fastq files (defaults to "1" (=yes))

-pREM/--pRemovalRedundantFiles Generating sbatch script for deletion of redundant files (defaults to "1" (=yes);Note: Must be submitted manually)

-al/--aligner Setting which aligner was used for alignment in previous analysis (defaults to "")

-wgs/--wholeGenomeSequencing Analysis to perform are based on whole genome sequencing data or not (defaults to "0" (=no))

-mc/--maximumCores The maximum number of cores per node used in the analysis (defaults to "8")

-env_up/--environmentUppmax Sets the environment to UPPMAX (defaults to "0" (=no))

-c/--configFile YAML config file for script parameters (defaults to "")

-wc/--writeConfigFile Write YAML config file with used script parameters. (defaults to "";Supply whole path)

=head3 I/O

Input format ( dir/infile.fastq or dir/infile.fastq.gz)

Output format

1. FastQC files/Gzipped fastq files
2. Mosaik.dat 
3. Mosaik.bam
4. Mosaik_sorted.bam(.bai)
5. Mosaik_lanes_sorted_merged.bam(.bai)
6. calculateCoverage_results
7. Coverage plots
8. Bwa.sai
9. Bwa.sam
10. Bwa.bam
11. Bwa_sorted.bam(.bai)
12. Bwa_lanes_sorted_merged.bam

=head4 Dependencies

Local installation of:
Mosaik
BWA
SamTools
BedTools
PicardTools
qaCompute

Located in -rd, reference dir
Genome reference

Mosaik
1. .dat files of genome reference
2. Jump database keys, positions, meta
3. Neural network .ann (PE & SE)

BWA
1. BWA index files (amb,ann,bwt etc)

calculateCoverage
1. Target file
2. Genome reference file
3. Target.infile_list
4. Padded bed.infile_list

Located in -ids, inScriptDir
R scripts
1. covplots_genome.R
2. covplots_exome.R 

=cut
    
use Pod::Usage;
use Pod::Text;
use Getopt::Long;
use POSIX;

use vars qw($USAGE);

BEGIN {
    $USAGE =
	qq{mip_align.pl  -id [inFilesDir,.,.,.,n] -ids [inScriptDir,.,.,.,n] -rd [refdir] -p [project ID] -s [sample ID...n] -em [e-mail] -ods [outdirscripts] -odf [outDataDir] -f [familyID] -p[program]
	       -ifd/--inFilesDir Infile directory(s), comma sep (Mandatory: Supply whole path,)
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
	       -pFQC/--pFastQC Sequence quality analysis using FastQC (defaults to "1" (=1))
	       -pMoB/--pMosaikBuild  Convert reads to Mosaik format using MosaikBuild (defaults to "1" (=yes))
               -mobmfl/--mosaikBuildMedianFragLength Flag for setting the mean fragment length, mfl, (defaults to (=375) bp)
	       -pMoA/--pMosaikAlign Align reads using MosaikAlign (defaults to "1" (=yes))
               -moaref/--mosaikAlignReference MosaikAlign reference (defaults to "")
               -moaannpe/--mosaikAlignNeuralNetworkPeFile MosaikAlign Neural Network PE File (defaults to "")
               -moaannse/--mosaikAlignNeuralNetworkSeFile MosaikAlign Neural Network SE File (defaults to "")
               -mojdb/--mosaikJumpDbStub MosaikJump stub (defaults to "")
               -pBWA_aln/--pBwaAln Index reads using BWA Aln (defaults to "0" (=no))
               -bwaalnq/--bwaAlnQualityTrimming BWA Aln quality threshold for read trimming (defaults to "20")
               -pBWA_sampe/--pBwaSampe Align reads using BWA Sampe (defaults to "0" (=no))
               -pSamT_sort/--pSamToolsSort Sort & index aligned reads using SamTools sort & index (defaults to "1" (=yes))
               -pPicT_merge/--pPicardToolsMerge Merge (BAM file(s) ) using PicardTools MergeSamFiles (defaults to "0" (=no))
               -picT_mergeprev/--picardToolsMergePreviousFiles PicardTools MergeSamFiles on merged current files and previous BAM-file(s) (Supply whole path and name, name must contain sample id, and lanes_Xn info)
               -pPicT_MarkDup/--pPicardToolsMarkduplicates Markduplicates using PicardTools MarkDuplicates (defaults to "1" (=yes))
               -pic_path/--picardToolsPath Path to PicardTools. Mandatory for use of PicardTools (defaults to "")
               -pCC/--pCalculateCoverage Use coverage calculation tools: qaCompute, genomeCoverageBED and PicardTools (MultipleMetrics & HSmetrics) (defaults to "1" (=yes))
               -pCC_Bedgc/--pGenomeCoverageBED Genome coverage calculation using genomeCoverageBED under '-pCC' (defaults to "1" (=yes))
               -pCC_Bedc/--pCoverageBED BED file coverage calculation using coverageBED under '-pCC' (defaults to "1" (=yes))
               -extb/--exomeTargetBed Target BED file of exome capture for coverageBed '-pCC_Bedc' (defaults to "")
               -pCC_Qac/--pQaCompute Genome coverage calculation using qaCompute under '-pCC' (defaults to "1" (=yes))
               -xcov/--xCoverage Max coverage depth when using '-genomeCoverageBED', '-qaCompute' (defaults to "30")
               -pCC_PicMM/--pPicardToolsCollectMultipleMetrics Metrics calculation using PicardTools collectMultipleMetrics under '-pCC' (defaults to "1" (=yes))
               -pCCE_PicHS/--pPicardToolsCalculateHSMetrics Capture calculation using PicardTools CalculateHSmetrics under '-pCC' (defaults to "1" (=yes))
               -extbl/--exomeTargetBedInfileList Prepared target BED file for PicardTools CalculateHSMetrics (defaults to "". File ending should be ".infile_list") 
               -extpbl/--exomeTargetPaddedBedInfileList Prepared padded target BED file for PicardTools CalculateHSMetrics (defaults to "". File ending should be ".padXXX.infile_list")
	       -pRCP/--pRCovPlots Plots of genome coverage using rCovPlots (defaults to "1" (=yes))
	       -pGZ/--pGZip GZip fastq files (defaults to "1" (=yes))
               -pREM/--pRemovalRedundantFiles Generating sbatch script for deletion of redundant files (defaults to "1" (=yes);Note: Must be submitted manually)
               -al/--aligner Setting which aligner was used for alignment in previous analysis (defaults to "")
               -wgs/--wholeGenomeSequencing Analysis to perform are based on whole genome sequencing data or not (defaults to "0" (=no))
               -mc/--maximumCores The maximum number of cores per node used in the analysis (defaults to "8")
               -env_up/--environmentUppmax Sets the environment to UPPMAX (defaults to "0" (=no))
               -c/--configFile YAML config file for script parameters (defaults to "")
               -wc/--writeConfigFile Write YAML configuration file for script parameters (defaults to "";Supply whole path)
	   };
}

###
#Program parameters
###

#Project specific
my ($projectID,$email, $wholeGenomeSequencing, $inScriptDir, $referencesDir, $outDataDir, $outScriptDir, $familyID, $pedigreeFile, $configFile, $writeConfigFile) = (0,0,-1,0,0,0,0,0,0,0,0);
my (@inFilesDir,@sampleIDs); #Arrays for input file directorys,sampleIDs

#GZip
my ($pGZ) = (-1);
#FastQC
my ($pFQC) = (-1);

#Mosaik
my ($pMoB, $pMoA) = (-1,-1);
my ($mosaikBuildMedianFragLength, $mosaikAlignReference, $mosaikAlignNeuralNetworkPeFile, $mosaikAlignNeuralNetworkSeFile, $mosaikJumpDbStub) = (-1,0,0,0,0);

#BWA
my ($pBWA_aln, $pBWA_sampe) = (-1,-1);
my ($bwaAlnQualityTrimming) = (-1);

#SamTools
my ($pSamT_sort) = (-1);

#PicardTools
my ($pPicT_merge, $pPicT_MarkDup, $picardToolsPath) = (-1,-1,0);
my (@picardToolsMergePreviousFiles);

#Coverage
my ($pCC, $pCC_Bedgc, $pCC_Bedc, $pCC_Qac, $pCC_PicMM, $pCCE_PicHS, $pGdb, $pRCP) = (-1,-1,-1,-1,-1,-1,-1,-1);
my ($exomeTargetBed, $exomeTargetBedInfileList, $exomeTargetPaddedBedInfileList, $xCoverage, $identicalCaptureBedCounter, $identicalCaptureBedIntervalCounter, $identicalCaptureBedPaddedIntervalCounter) = (0,0,0,-1,0,0,0);

#Pipe
my ($pREM) = (-1);
my ($humanGenomeReference, $humanGenomeReferenceSource, $humanGenomeRefereceChromosomePrefix, $humanGenomeReferenceVersion, $fnend, $aligner, $maximumCores, $environmentUppmax, $filename, $fnt, $fnt2, $help) = (0, 0, 0, 0, ".sh", 0,0,-1); #Arguments for project
my (@chromosomes);
my (%infiles, %indirpath, %InfilesLaneNoEnding, %lanes, %InfilesBothStrandsNoEnding, %jobID, %paralleljobID, %allsampleIDjobID, %sampleInfo, %scriptParameters); 
#%infiles=from platform (Illumina), %indirpath for the path to infiles, %InfilesLaneNoEnding for MosaikBuild (one entry for both strands), %lanes for sample lanes, InfilesBothStrandsNoEnding for bwa_aln (one entry per strand)

###
#Staging Area
###

#Capture kits supported from pedigree file.
my %supportedCaptureKits = (
    'Nimblegen_SeqCapEZExome.V2' => "Nimblegen_SeqCapEZExome.V2.GenomeReferenceSourceVersion_targets_ChromosomePrefix.bed",
    'Agilent_SureSelect.V2' => "Agilent_SureSelect.V2.GenomeReferenceSourceVersion_targets_ChromosomePrefix.bed",
    'Agilent_SureSelect.V3' => "Agilent_SureSelect.V3.GenomeReferenceSourceVersion_targets_ChromosomePrefix.bed",
    'Agilent_SureSelect.V4' => "Agilent_SureSelect.V4.GenomeReferenceSourceVersion_targets_ChromosomePrefix.bed",
    'Agilent_SureSelect.V5' => "Agilent_SureSelect.V5.GenomeReferenceSourceVersion_targets_ChromosomePrefix.bed",
    );

#
#User Options
#

GetOptions('ifd|inFilesDir:s'  => \@inFilesDir, #Comma separated list
	   'isd|inScriptDir:s'  => \$inScriptDir, #Directory for custom scripts required by the pipeline
	   'rd|referencesDir:s'  => \$referencesDir, #directory containing references
	   'p|projectID:s'  => \$projectID,
	   's|sampleIDs:s'  => \@sampleIDs, #Comma separated list, one below outDataDir
	   'em|email:s'  => \$email,
	   'odd|outDataDir:s'  => \$outDataDir, #One dir above sample id, must supply whole path i.e. /proj/...
	   'osd|outScriptDir:s'  => \$outScriptDir,  #One dir above sample id, must supply whole path i.e. /proj/...
	   'f|familyID:s' => \$familyID, #Family group ID (Merged to same vcf file after GATK Base Recalibration)
	   'pedigree|pedigreeFile:s' => \$pedigreeFile, #Pedigree file
	   'huref|humanGenomeReference:s' => \$humanGenomeReference, #Human genome reference
	   'pFQC|pFastQC:n' => \$pFQC,
	   'pMoB|pMosaikBuild:n' => \$pMoB,
	   'mobmfl|mosaikBuildMedianFragLength:n' => \$mosaikBuildMedianFragLength, #for fragment length estimation and local search
	   'pMoA|pMosaikAlign:n' => \$pMoA,
	   'moaref|mosaikAlignReference:s' => \$mosaikAlignReference, #MosaikAlign reference file assumes existance of jump database files in same dir
	   'moaannpe|mosaikAlignNeuralNetworkPeFile:s' => \$mosaikAlignNeuralNetworkPeFile,
	   'moaannse|mosaikAlignNeuralNetworkSeFile:s' => \$mosaikAlignNeuralNetworkSeFile, 
	   'mojdb|mosaikJumpDbStub:s' => \$mosaikJumpDbStub, #Stub for MosaikJump database
	   'pBWA_aln|pBwaAln:n' => \$pBWA_aln,
	   'bwaalnq|bwaAlnQualityTrimming:n' => \$bwaAlnQualityTrimming, #BWA aln quality threshold for read trimming down to 35bp
	   'pBWA_sampe|pBwaSampe:n' => \$pBWA_sampe,
	   'pSamT_sort|pSamToolsSort:n' => \$pSamT_sort,
	   'pPicT_merge|pPicardToolsMerge:n' => \$pPicT_merge, #PicardTools MergeSamFiles
	   'pict_mergeprev|picardToolsMergePreviousFiles:s' => \@picardToolsMergePreviousFiles, #Comma separated list
	   'pPicT_MarkDup|pPicardToolsMarkduplicates:s' => \$pPicT_MarkDup, #PicardTools MarkDuplicates
	   'pic_path|picardToolsPath:s' => \$picardToolsPath, #Path to picardtools
	   'pCC|pCalculateCoverage:n' => \$pCC,
	   'pCC_Bedgc|pGenomeCoverageBED:n' => \$pCC_Bedgc,
	   'pCC_Bedc|pCoverageBED:n' => \$pCC_Bedc,
	   'extb|exomeTargetBed:s' => \$exomeTargetBed, #target file for coverageBed
	   'pCC_Qac|pQaCompute:n' => \$pCC_Qac,
	   'xcov|xCoverage:n' => \$xCoverage, #Sets max depth to calculate coverage
	   'pCC_PicMM|pPicardToolsCollectMultipleMetrics:n' => \$pCC_PicMM,
	   'pCCE_PicHS|pPicardToolsCalculateHSMetrics:n' => \$pCCE_PicHS,
	   'extbl|exomeTargetBedInfileList:s' => \$exomeTargetBedInfileList, #target file for CalculateHsMetrics
	   'extpbl|exomeTargetPaddedBedInfileList:s' => \$exomeTargetPaddedBedInfileList, #Padded target file for CalculateHsMetrics
	   'pRCP|pRCovPlots:n' => \$pRCP,
	   'pGZ|pGZip:n' => \$pGZ,
	   'pREM|pRemovalRedundantFiles:n' => \$pREM,
	   'al|aligner:s' => \$aligner, #determining which aligner was used previously (if not specified)
	   'wgs|wholeGenomeSequencing:n' => \$wholeGenomeSequencing,
	   'mc|maximumCores:n' => \$maximumCores, #Per node
	   'env_up|environmentUppmax:n' => \$environmentUppmax, #Sets several default paths, so that they do not have to be supplied
	   'c|configFile:s' => \$configFile,
	   'wc|writeConfigFile:s' => \$writeConfigFile,
	   'h|help' => \$help,
	   );

die $USAGE if( $help );

if ($configFile ne 0) {
    
#    if ($environmentUppmax == 1) {
	use lib '/bubo/home/h12/henriks/lib/'; #YAML not installed at @UPPMAX
 #   }
    #else {
#	use lib '/bubo/home/h12/henriks/lib/'; #Remove later 
 #   }
    use YAML;
    open (YAML, "<".$configFile) or die "can't open ".$configFile.": $!\n";
    %scriptParameters = YAML::LoadFile($configFile);
    close(YAML);
}

if ($projectID eq 0) { #No input from cmd line   
    unless ( defined($scriptParameters{'projectID'}) ) { #No input from config file   
	print STDERR "\n";
	print STDERR "Must supply a project ID. Use '-projectID'", "\n\n";
	die $USAGE;
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $scriptParameters{'projectID'} = $projectID;   
}

if ($email eq 0) { #No input from cmd line
    if ( defined($scriptParameters{'email'}) ) { #Input from config file - Do nothing
    }
    else {
	$scriptParameters{'email'} = 0; #Default
    }  
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $scriptParameters{'email'} = $email;   
}

if ( $familyID eq 0 ) {
     if ( defined($scriptParameters{'familyID'}) ) { #Input from config file - Do nothing
     }
     else {
	 print STDERR "\n";
	 print STDERR "Must supply a family id. Use '-familyID'. -If not applicable supply the same familyID as the sampleid ", "\n\n";
	 die $USAGE;
     }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $scriptParameters{'familyID'} = $familyID;   
}

#wholeGenomeSequencing
if( $wholeGenomeSequencing == -1) { #No input from cmd
    if ( defined($scriptParameters{'wholeGenomeSequencing'}) ) { #Input from config file - do nothing
    }
    else {
	$scriptParameters{'wholeGenomeSequencing'} = 0; #Default
    }
}
else {
    $scriptParameters{'wholeGenomeSequencing'} = $wholeGenomeSequencing;
}

#$environmentUppmax
if( $environmentUppmax == -1) { #No input from cmd
    if ( defined($scriptParameters{'environmentUppmax'}) ) { #Input from config file - do nothing
    }
    else {
	$scriptParameters{'environmentUppmax'} = 0; #Default
    }
}
else {
    $scriptParameters{'environmentUppmax'} = $environmentUppmax;
}

if ($outDataDir eq 0) {
    
    if ( defined($scriptParameters{'outDataDir'}) ) { #Input from config file - do nothing	
    }
    elsif ($scriptParameters{'environmentUppmax'} == 1) { #Use Uppmax default
	print STDOUT "\n";
	if ($scriptParameters{'wholeGenomeSequencing'} == 1) {
	    $outDataDir = "/proj/".$scriptParameters{'projectID'}."/private/nobackup/genomes";
	}
	else {
	    $outDataDir = "/proj/".$scriptParameters{'projectID'}."/private/nobackup/exomes";
	}
	print STDOUT "Setting MIP output data dir to: ".$outDataDir, "\n\n";
	$scriptParameters{'outDataDir'} = $outDataDir;
    }
    else {
	print STDERR "\n";
	print STDERR "Must supply a MIP output data dir", "\n\n";
	die $USAGE;
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $scriptParameters{'outDataDir'} = $outDataDir;   
}

if ($pedigreeFile eq 0) {
    
    if ( defined($scriptParameters{'pedigreeFile'}) ) { #Input from config file - ReadPedigreefile
	if ( scalar(@sampleIDs) == 0 ) { #No user supplied sample info
	    if ( defined($scriptParameters{'sampleIDs'}) ) { #sampleIDs info in config file
		ReadPedigreeFile($scriptParameters{'pedigreeFile'}, 1);  # scalar(@sampleIDs) = 0:No user supplied sample info, but present in config file do NOT overwrite using info from pedigree file
	    }
	    else { #No sampleIDs info in config file
		ReadPedigreeFile($scriptParameters{'pedigreeFile'}, scalar(@sampleIDs) );  # scalar(@sampleIDs) = 0:No user supplied sample info, not defined $scriptParameters{'sampleIDs'} in config file, add it from pedigree file
	    }
	}
	else {
	    ReadPedigreeFile($scriptParameters{'pedigreeFile'}, scalar(@sampleIDs) );  # User supplied sample info, do NOT overwrite using info from pedigree file
	}
    }
    elsif ($scriptParameters{'environmentUppmax'} == 1) {
	print STDOUT "\n";
	if ( $scriptParameters{'wholeGenomeSequencing'} ==0) {
	    $pedigreeFile = "/proj/".$scriptParameters{'projectID'}."/private/exomes/".$scriptParameters{'familyID'}."/".$scriptParameters{'familyID'}."_pedigree.txt";
	}
	else {
	    $pedigreeFile = "/proj/".$scriptParameters{'projectID'}."/private/genomes/".$scriptParameters{'familyID'}."/".$scriptParameters{'familyID'}."_pedigree.txt";
	}
#	$pedigreeFile = $scriptParameters{'outDataDir'}."/".$scriptParameters{'familyID'}."/".$scriptParameters{'familyID'}."_pedigree.txt";
	print STDOUT "Assuming location of pedigree file to be: ".$pedigreeFile, "\n\n";
	if (-e $pedigreeFile) { #if file exists 
	    print STDOUT "Found pedigree file at: ".$pedigreeFile, "\n\n";
	    $scriptParameters{'pedigreeFile'} = $pedigreeFile; #Add to enable recreation of cmd line later
	    ReadPedigreeFile($pedigreeFile, scalar(@sampleIDs) ); #  scalar(@sampleIDs)= 0:No user supplied sample info, add it from pedigree file
	}
	else { 
	    print STDERR "Could not find pedigree file at: ".$pedigreeFile, "\n";
	    die $USAGE;
	} 
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    ReadPedigreeFile($pedigreeFile, scalar(@sampleIDs) );  # scalar(@sampleIDs) = 0:No user supplied sample info, add it from pedigree file
    $scriptParameters{'pedigreeFile'} = $pedigreeFile;   
}

if (scalar (@inFilesDir) == 0) {
    if ( defined($scriptParameters{'inFilesDir'}) ) { #Input from config file - do nothing
    }
    elsif ($scriptParameters{'environmentUppmax'} == 1 && $scriptParameters{'pedigreeFile'} ne 0) { #Locate infile directory. Order dictated by @sampleIDs (lexiographical)
	for (my $indirectoryCount=0;$indirectoryCount<scalar(@sampleIDs);$indirectoryCount++) {
	    if ($scriptParameters{'wholeGenomeSequencing'} == 0) {
		push(@inFilesDir, "/proj/".$scriptParameters{'projectID'}."/private/exomes/".$sampleIDs[$indirectoryCount]."/fastq");
	    }
	    else {
		push(@inFilesDir, "/proj/".$scriptParameters{'projectID'}."/private/genomes/".$sampleIDs[$indirectoryCount]."/fastq");
	    }
	}
	$scriptParameters{'inFilesDir'} = join(',',@inFilesDir); #Add to enable recreation of cmd line later
    }
    else {
	my $verbosity = 2;
	print"\n";
	pod2usage({-message => "Must supply an infile directory as comma separated list.\n",
		   -verbose => $verbosity
	 });
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $scriptParameters{'inFilesDir'} = join(',',@inFilesDir);   
}

if ( scalar(@sampleIDs) == 0) {
    if ( defined($scriptParameters{'sampleIDs'}) ) { #Input from config file - do nothing
	@sampleIDs = split(/,/, $scriptParameters{'sampleIDs'}); #Transfer to array
    }
    else {
	print STDERR "\n";
	print STDERR "Must supply a sample ID as a comma separated list", "\n\n";
	die $USAGE;
    }
}
else {
    $scriptParameters{'sampleIDs'} = join(',',@sampleIDs); #Add to enable recreation of cmd line later
    @sampleIDs = split(/,/,join(',',@sampleIDs)); #Enables comma separated list of sample IDs from user supplied cmd info
}

if ( $inScriptDir eq 0) {
    if ( defined($scriptParameters{'inScriptDir'}) ) { #Input from config file - do nothing
	
    }
    elsif ($scriptParameters{'environmentUppmax'} == 1) {
	print STDOUT "\n";
	$inScriptDir = "/proj/".$scriptParameters{'projectID'}."/private/mip_scripts_master";
	print STDOUT "Setting the MIP scripts dir to: ".$inScriptDir, "\n\n";
	$scriptParameters{'inScriptDir'} = $inScriptDir; #Add to enable recreation of cmd line later
    }
    else {
	print STDERR "\n";
	print STDERR "Must supply the MIP scripts dir", "\n\n";
	die $USAGE;
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $scriptParameters{'inScriptDir'} = $inScriptDir;   
}

if ( $referencesDir eq 0) {
    if ( defined($scriptParameters{'referencesDir'}) ) { #Input from config file - do nothing
	
    }
    elsif ($scriptParameters{'environmentUppmax'} == 1) {
	print STDOUT "\n";
	$referencesDir = "/proj/".$scriptParameters{'projectID'}."/private/mip_references";
	print STDOUT "Setting MIP reference dir to: ".$referencesDir, "\n\n";
	$scriptParameters{'referencesDir'} = $referencesDir; #Add to enable recreation of cmd line later
    }
    else {
	print STDERR "\n";
	print STDERR "Must supply a MIP reference dir", "\n\n";
	die $USAGE;
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $scriptParameters{'referencesDir'} = $referencesDir; #Add to enable recreation of cmd line later
}

if ($outScriptDir eq 0) {
    if ( defined($scriptParameters{'outScriptDir'}) ) { #Input from config file - do nothing	
    }
    elsif ($scriptParameters{'environmentUppmax'} == 1) {
	print STDOUT "\n";
	if ($scriptParameters{'wholeGenomeSequencing'} == 1) {
	    $outScriptDir = "/proj/".$scriptParameters{'projectID'}."/private/genomes_scripts";
	}
	else {
	    $outScriptDir = "/proj/".$scriptParameters{'projectID'}."/private/exomes_scripts";
	}
	print STDOUT "Setting MIP output scripts dir to: ".$outScriptDir, "\n\n";
	$scriptParameters{'outScriptDir'} = $outScriptDir; #Add to enable recreation of cmd line later
    }
    else {
	print STDERR "\n";
	print STDERR "Must supply a MIP output script dir", "\n\n";
	die $USAGE;
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $scriptParameters{'outScriptDir'} = $outScriptDir; #Add to enable recreation of cmd line later
}

#
#GZip
#
if ( $pGZ == -1) { #Not set by cmd 
    if ( defined($scriptParameters{'pGZip'}) ) { #Input from config file - do nothing	
    }
    else {
	$scriptParameters{'pGZip'} = 1; #Default
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $scriptParameters{'pGZip'} = $pGZ;
}

#
#FASTQC
#
if ( $pFQC == -1) { #Not set by cmd 
    if ( defined($scriptParameters{'pFastQC'}) ) { #Input from config file - do nothing	
    }
    else {
	$scriptParameters{'pFastQC'} = 1; #Default
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $scriptParameters{'pFastQC'} = $pFQC;
}

#
#MosaikBuild
#
if ( $pMoB == -1) { #Not set by cmd 
    if ( defined($scriptParameters{'pMosaikBuild'}) ) { #Input from config file - do nothing	
    }
    else {
	$scriptParameters{'pMosaikBuild'} = 1; #Default
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $scriptParameters{'pMosaikBuild'} = $pMoB;
}
if ( $scriptParameters{'pMosaikBuild'} > 0 ) { #MosaikBuild is to be used
    if ( $mosaikBuildMedianFragLength == -1) {
	if ( defined($scriptParameters{'mosaikBuildMedianFragLength'}) ) { #Input from config file - do nothing	
	}
	else {
	    $scriptParameters{'mosaikBuildMedianFragLength'} = 375; #Default
	}
    }
    else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
	$scriptParameters{'mosaikBuildMedianFragLength'} = $mosaikBuildMedianFragLength;
    }
}

#
#MosaikAlign
#
if ( $pMoA == -1) { #Not set by cmd
    if ( defined($scriptParameters{'pMosaikAlign'}) ) { #Input from config file - do nothing	
    }
    else {
	$scriptParameters{'pMosaikAlign'} = 1; #Default
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $scriptParameters{'pMosaikAlign'} = $pMoA;
}
if ( $scriptParameters{'pMosaikAlign'} > 0 ) { #MosaikAlign is to be used - check prerequisets
    
    if ($mosaikAlignReference eq 0) { #No input from cmd
	if ( defined($scriptParameters{'mosaikAlignReference'}) ) { #Input from config file - do nothing	
	}
	elsif ($scriptParameters{'environmentUppmax'} == 1) { #Use default 
	    $mosaikAlignReference = "Homo_sapiens.GRCh37.70_nochr.dat";
	    $scriptParameters{'mosaikAlignReference'} = $mosaikAlignReference;
	}
	else {
	    print STDERR "\nSupply '-mosaikAlignReference' if you want to run MosaikAlign\n\n";
	    die $USAGE;
	}
    }
    else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
	$scriptParameters{'mosaikAlignReference'} = $mosaikAlignReference;
    }
    if ($mosaikAlignNeuralNetworkPeFile eq 0) { #No input from cmd
	if ( defined($scriptParameters{'mosaikAlignNeuralNetworkPeFile'}) ) { #Input from config file - do nothing	
	}
	elsif ($scriptParameters{'environmentUppmax'} == 1) {
	    $mosaikAlignNeuralNetworkPeFile = "2.1.78.pe.ann";
	    $scriptParameters{'mosaikAlignNeuralNetworkPeFile'} = $mosaikAlignNeuralNetworkPeFile;
	}
	else {
	    print STDERR "\nSupply '-mosaikAlignNeuralNetworkPeFile' if you want to run MosaikAlign\n\n";
	    die $USAGE;
	}
    }
    else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
	$scriptParameters{'mosaikAlignNeuralNetworkPeFile'} = $mosaikAlignNeuralNetworkPeFile;
    }
    if ($mosaikAlignNeuralNetworkSeFile eq 0) { #No input from cmd
	if ( defined($scriptParameters{'mosaikAlignNeuralNetworkSeFile'}) ) { #Input from config file - do nothing	
	}
	elsif ($scriptParameters{'environmentUppmax'} == 1) {
	    $mosaikAlignNeuralNetworkSeFile = "2.1.78.se.ann";
	    $scriptParameters{'mosaikAlignNeuralNetworkSeFile'} = $mosaikAlignNeuralNetworkSeFile;
	}
	else {
	    print STDERR "\nSupply '-mosaikAlignNeuralNetworkSeFile' if you want to run MosaikAlign\n\n";
	    die $USAGE;
	}
    }
    else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
	$scriptParameters{'mosaikAlignNeuralNetworkSeFile'} = $mosaikAlignNeuralNetworkSeFile;
    }
    if ($mosaikJumpDbStub eq 0) { #No input from cmd
	if ( defined($scriptParameters{'mosaikJumpDbStub'}) ) { #Input from config file - do nothing	
	}
	elsif ($scriptParameters{'environmentUppmax'} == 1) {
	    $mosaikJumpDbStub = "Homo_sapiens.GRCh37.70_nochr_jdb_15";
	    $scriptParameters{'mosaikJumpDbStub'} = $mosaikJumpDbStub;
	}
	else {
	    print STDERR "\nSupply '-mosaikJumpDbStub' if you want to run MosaikAlign\n\n";
	    die $USAGE;
	}
    }
    else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
	$scriptParameters{'mosaikJumpDbStub'} = $mosaikJumpDbStub;
    }

#File existence checks
    unless (-e $scriptParameters{'referencesDir'}."/".$scriptParameters{'mosaikAlignReference'} ) { #Check existence of mosaik hashed reference in supplied reference dir
	print STDERR "\nCould not find intended Mosaik Reference file: ".$scriptParameters{'referencesDir'}."/".$scriptParameters{'mosaikAlignReference'}, "\n\n";
	die $USAGE;		
    }
    unless (-e $scriptParameters{'referencesDir'}."/".$scriptParameters{'mosaikAlignNeuralNetworkPeFile'} ) { #Check existence of mosaik neural network file in supplied reference dir
	print STDERR "\nCould not find intended Mosaik Neural Network file: ".$scriptParameters{'referencesDir'}."/".$scriptParameters{'mosaikAlignNeuralNetworkPeFile'}, "\n\n";
	die $USAGE;		
    }
    unless (-e $scriptParameters{'referencesDir'}."/".$scriptParameters{'mosaikAlignNeuralNetworkSeFile'} ) { #Check existence of mosaik neural network file in supplied reference dir
	print STDERR "\nCould not find intended Mosaik Neural Network file: ".$scriptParameters{'referencesDir'}."/".$scriptParameters{'mosaikAlignNeuralNetworkSeFile'}, "\n\n";
	die $USAGE;		
    }
    unless (-e $scriptParameters{'referencesDir'}."/".$scriptParameters{'mosaikJumpDbStub'}."_keys.jmp" ) { #Check existence of Mosaik jumpDb in supplied reference dir
	print STDERR "\nCould not find intended Mosaik JumpDb: ".$scriptParameters{'referencesDir'}."/".$scriptParameters{'mosaikJumpDbStub'}."_keys.jmp\n\n";
	die $USAGE;		
    }
    unless (-e $scriptParameters{'referencesDir'}."/".$scriptParameters{'mosaikJumpDbStub'}."_meta.jmp" ) { #Check existence of Mosaik jumpDb in supplied reference dir
	print STDERR "\nCould not find intended Mosaik JumpDb: ".$scriptParameters{'referencesDir'}."/".$scriptParameters{'mosaikJumpDbStub'}."_meta.jmp\n\n";
	die $USAGE;		
    }
    unless (-e $scriptParameters{'referencesDir'}."/".$scriptParameters{'mosaikJumpDbStub'}."_positions.jmp" ) { #Check existence of Mosaik jumpDb in supplied reference dir
	print STDERR "\nCould not find intended Mosaik JumpDb: ".$scriptParameters{'referencesDir'}."/".$scriptParameters{'mosaikJumpDbStub'}."_positions.jmp\n\n";
	die $USAGE;		
    }
}


#
#BWA_aln
#
if ( $pBWA_aln == -1) { #No input from cmd
    if ( defined($scriptParameters{'pBwaAln'}) ) { #Input from config file - do nothing	
    }
    else {
	$scriptParameters{'pBwaAln'} = 0; #Default
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $scriptParameters{'pBwaAln'} = $pBWA_aln;
}
if ( $scriptParameters{'pBwaAln'} > 0 ) {
    if ( $bwaAlnQualityTrimming == -1) { #No input from cmd
	if ( defined($scriptParameters{'bwaAlnQualityTrimming'}) ) { #Input from config file - do nothing	
	}
	else {
	    $scriptParameters{'bwaAlnQualityTrimming'} = 20; #Default
	}
    }
    else {
	$scriptParameters{'bwaAlnQualityTrimming'} = $bwaAlnQualityTrimming;
    }
}

#
#BWA_sampe
#
if ( $pBWA_sampe == -1) { #No input from cmd
    if ( defined($scriptParameters{'pBwaSampe'}) ) { #Input from config file - do nothing	
    }
    else {
	$scriptParameters{'pBwaSampe'} = 0; #Default
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $scriptParameters{'pBwaSampe'} = $pBWA_sampe;
}

if ( $aligner eq 0 ) { #No input from cmd
    unless ( defined($scriptParameters{'aligner'}) ) { #No input from config file
	if ( ($scriptParameters{'pMosaikBuild'} ==1) || ($scriptParameters{'pMosaikAlign'} ==1) ) { #Mosaik track
	    if ( ($scriptParameters{'pBwaAln'} ==0) && ($scriptParameters{'pBwaSampe'} ==0) ) {
		$scriptParameters{'aligner'} = "mosaik";
	    }
	    else {
		print STDERR "\n";
		print STDERR "You have to choose either mosaik or bwa to perform alignments or specify which aligner (-aligner 'mosaik' or 'bwa') was used if you want to only run programs after alignment.", "\n\n";
		die $USAGE;
	    }
	}
	if ( ($scriptParameters{'pBwaAln'} ==1) || ($scriptParameters{'pBwaSampe'} ==1) ) { #BWA track
	    if ( ($scriptParameters{'pMosaikBuild'} ==0) && ($scriptParameters{'pMosaikAlign'} ==0) ) {
		$scriptParameters{'aligner'} = "bwa";
	    }
	    else {
		print STDERR "\n";
		print STDERR "You have to choose either mosaik or bwa to perform alignments or specify which aligner (-aligner 'mosaik' or 'bwa') was used if you want to only run programs after alignment.", "\n\n";
		die $USAGE;
	    }
	}
    }
}
else {
    $scriptParameters{'aligner'} = $aligner;
}

#
#SamTools Sort/Index
#
if ( $pSamT_sort == -1) { #No input from cmd
    if ( defined($scriptParameters{'pSamToolsSort'}) ) { #Input from config file - do nothing	
    }
    else {
	$scriptParameters{'pSamToolsSort'} = 1; #Default
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $scriptParameters{'pSamToolsSort'} = $pSamT_sort;
}

#
#PicardToolsMerge
#
if ( $pPicT_merge == -1) { #No input from cmd
    if ( defined($scriptParameters{'pPicardToolsMerge'}) ) { #Input from config file - do nothing	
    }
    else {
	$scriptParameters{'pPicardToolsMerge'} = 0; #Default
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $scriptParameters{'pPicardToolsMerge'} = $pPicT_merge;
}
if ( $scriptParameters{'pPicardToolsMerge'} > 0) {
    if (scalar(@picardToolsMergePreviousFiles) == 0 ) {
	if ( defined($scriptParameters{'picardToolsMergePreviousFiles'}) ) { #Input from config file - transfer to array
	    @picardToolsMergePreviousFiles = split(/,/, $scriptParameters{'picardToolsMergePreviousFiles'}); #Transfer to array
	}
    }
    else {
	$scriptParameters{'picardToolsMergePreviousFiles'} = join(',',@picardToolsMergePreviousFiles);
    }
}

#
#PicardToolsMarkduplicates
#
if ( $pPicT_MarkDup == -1) { #No input from cmd
    if ( defined($scriptParameters{'pPicardToolsMarkduplicates'}) ) { #Input from config file - do nothing	
    }
    else {
	$scriptParameters{'pPicardToolsMarkduplicates'} = 1; #Default
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $scriptParameters{'pPicardToolsMarkduplicates'} = $pPicT_MarkDup;
}

#
#Coverage Calculations
#
if ( $pCC == -1) { #No input from cmd
    if ( defined($scriptParameters{'pCalculateCoverage'}) ) { #Input from config file - do nothing	
    }
    else {
	$scriptParameters{'pCalculateCoverage'} = 1; #Default
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $scriptParameters{'pCalculateCoverage'} = $pCC;
}
if ( $scriptParameters{'pCalculateCoverage'} > 0) {
    if ( $pCC_Bedgc == -1) { #No input from cmd
	if ( defined($scriptParameters{'pGenomeCoverageBED'}) ) { #Input from config file - do nothing	
	}
	else {
	    $scriptParameters{'pGenomeCoverageBED'} = 1; #Default
	}
    }
    else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
	$scriptParameters{'pGenomeCoverageBED'} = $pCC_Bedgc;
    }
    if ( $pCC_Bedc == -1) { #No input from cmd
	if ( defined($scriptParameters{'pCoverageBED'}) ) { #Input from config file - do nothing	
	}
	else {
	    $scriptParameters{'pCoverageBED'} = 1; #Default
	}
    }
    else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
	$scriptParameters{'pCoverageBED'} = $pCC_Bedc;
    }
    if ( $pCC_Qac == -1) { #No input from cmd
	if ( defined($scriptParameters{'pQaCompute'}) ) { #Input from config file - do nothing	
	}
	else {
	    $scriptParameters{'pQaCompute'} = 1; #Default
	}
    }
    else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
	$scriptParameters{'pQaCompute'} = $pCC_Qac;
    }
    if ( $pCC_PicMM == -1) { #No input from cmd
	if ( defined($scriptParameters{'pPicardToolsCollectMultipleMetrics'}) ) { #Input from config file - do nothing	
	}
	else {
	    $scriptParameters{'pPicardToolsCollectMultipleMetrics'} = 1; #Default
	}
    }
    else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
	$scriptParameters{'pPicardToolsCollectMultipleMetrics'} = $pCC_PicMM;
    }
    if ( $pCCE_PicHS == -1) { #No input from cmd
	if ( defined($scriptParameters{'pPicardToolsCalculateHSMetrics'}) ) { #Input from config file - do nothing	
	}
	else {
	    $scriptParameters{'pPicardToolsCalculateHSMetrics'} = 1; #Default
	}
    }
    else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
	$scriptParameters{'pPicardToolsCalculateHSMetrics'} = $pCCE_PicHS;
    }
}
if ( $scriptParameters{'pCalculateCoverage'} > 0) {
    if ( ($scriptParameters{'pQaCompute'} > 0) || ($scriptParameters{'pGenomeCoverageBED'} > 0) ) {
#xCoverage
	if( $xCoverage == -1) { #No input from cmd
	    if ( defined($scriptParameters{'xCoverage'}) ) { #Input from config file - do nothing
	    }
	    else {
		$scriptParameters{'xCoverage'} = 30; #Default
	    }
	}
	else {
	    $scriptParameters{'xCoverage'} = $xCoverage;
	}
    }
}

#
#RCovPlots
#
if ( $pRCP == -1) { #No input from cmd
    if ( defined($scriptParameters{'pRCovPlots'}) ) { #Input from config file - do nothing	
    }
    else {
	$scriptParameters{'pRCovPlots'} = 1; #Default
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $scriptParameters{'pRCovPlots'} = $pRCP;
}

#
#RemovalRedundantFiles
#
if ( $pREM == -1) { #No input from cmd
    if ( defined($scriptParameters{'pRemovalRedundantFiles'}) ) { #Input from config file - do nothing	
    }
    else {
	$scriptParameters{'pRemovalRedundantFiles'} = 1; #Default
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $scriptParameters{'pRemovalRedundantFiles'} = $pREM;
}

if ($humanGenomeReference eq 0) { #No input from cmd

    if ( ($scriptParameters{'pPicardToolsCalculateHSMetrics'} > 0) || ($scriptParameters{'pPicardToolsCollectMultipleMetrics'} > 0) || ($scriptParameters{'pBwaAln'} > 0)|| ($scriptParameters{'pBwaSampe'} > 0) ) { #Programs that uses $humanGenomeReference
	if ( defined($scriptParameters{'humanGenomeReference'}) ) { #Input from config file - do nothing
	    if ( $scriptParameters{'humanGenomeReference'} =~/^Homo_sapiens.GRCh(\d+\.\d+)/) { #Used to change capture kit genome reference version later
		$humanGenomeReferenceVersion = $1;
		$humanGenomeReferenceSource = "GRCh"; #Ensembl
		$humanGenomeRefereceChromosomePrefix = "nochr";
	    }
	    elsif ( $scriptParameters{'humanGenomeReference'} =~/^Homo_sapiens.hg(\d+)/) { #Used to change capture kit genome reference version later
		$humanGenomeReferenceVersion = $1;
		$humanGenomeReferenceSource = "hg"; #Refseq
		$humanGenomeRefereceChromosomePrefix = "chr";
	    }	
	}
	elsif ($scriptParameters{'environmentUppmax'} == 1) {
	    $humanGenomeReference = "Homo_sapiens.GRCh37.70_nochr.fasta";
	    if ( $humanGenomeReference =~/Homo_sapiens.GRCh(\d+\.\d+)/) {
		$humanGenomeReferenceVersion = $1;
		$humanGenomeReferenceSource = "GRCh"; #Ensembl
		$humanGenomeRefereceChromosomePrefix = "nochr";
	    }   
	    $scriptParameters{'humanGenomeReference'} = $humanGenomeReference; #Add to enable recreation of cmd line later
	}
	else {
	    print STDERR "\nSupply human genome reference fasta file to run 'pCCE_PicHS' and/or 'pCC_PicMM'\n\n";
	    die $USAGE
	}
    }
    else {
	$scriptParameters{'humanGenomeReference'} = $humanGenomeReference; #Add to enable recreation of cmd line later
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    if ( $humanGenomeReference =~/^Homo_sapiens.GRCh(\d+\.\d+)/) { #Used to change capture kit genome reference version later
	$humanGenomeReferenceVersion = $1;
	$humanGenomeReferenceSource = "GRCh"; #Ensembl
	$humanGenomeRefereceChromosomePrefix = "nochr";
    }
    elsif ( $humanGenomeReference =~/^Homo_sapiens.hg(\d+)/) { #Used to change capture kit genome reference version later
	$humanGenomeReferenceVersion = $1;
	$humanGenomeReferenceSource = "hg"; #Refseq
	$humanGenomeRefereceChromosomePrefix = "chr";
    }
    $scriptParameters{'humanGenomeReference'} = $humanGenomeReference; 
}
unless (-e $scriptParameters{'referencesDir'}."/".$scriptParameters{'humanGenomeReference'} ) { #Check for human genome reference in supplied reference dir
    print STDERR "\nCould not find human genome reference fasta file: ".$scriptParameters{'referencesDir'}."/".$scriptParameters{'humanGenomeReference'}, "\n\n";
    die $USAGE;
}

if ($picardToolsPath eq 0 ) { #No input from cmd
    if ( ($scriptParameters{'pPicardToolsMerge'} > 0) || ($scriptParameters{'pPicardToolsMarkduplicates'} > 0) || ($scriptParameters{'pPicardToolsCalculateHSMetrics'} > 0) || ($scriptParameters{'pPicardToolsCollectMultipleMetrics'} > 0) ) {
	if ( defined($scriptParameters{'picardToolsPath'}) ) { #Input from config file - do nothing	
	}
	elsif ($scriptParameters{'environmentUppmax'} == 1) {
	    $scriptParameters{'picardToolsPath'} = "/bubo/home/h12/henriks/programs/picard-tools-1.59"; #Add to enable recreation of cmd line later
	    print STDOUT "\nSet picard path to: ".$scriptParameters{'picardToolsPath'}, "\n"; 
	}
	else {
	    print STDERR "\nMust supply picardTools path to use PicardTools", "\n";
	    die $USAGE;
	}
    }
}
else { #Add to enable recreation of cmd line later
    $scriptParameters{'picardToolsPath'} = $picardToolsPath;
}

#
#exomeTargetBed
#

if ( $exomeTargetBed eq 0 ) { #No input from cmd
    if ( ($scriptParameters{'pCalculateCoverage'} > 0) && ($scriptParameters{'pCoverageBED'} > 0) ) {
	my $uncorrectCaptureCounter = 0; #Track no entries or wrong format entry in pedigree file
	for (my $sampleidCounter=0;$sampleidCounter<scalar(@sampleIDs);$sampleidCounter++) { #Check all samples
	    if ( defined($scriptParameters{ $sampleIDs[$sampleidCounter] }{'exomeTargetBed'}) ) { #Input from config file - transfer to sampleInfo
		$sampleInfo{$scriptParameters{'familyID'}}{$sampleIDs[$sampleidCounter]}{'exomeTargetBed'} = $scriptParameters{ $sampleIDs[$sampleidCounter] }{'exomeTargetBed'};
	    }	
	    elsif ($scriptParameters{'environmentUppmax'} == 1) {
		
		if ( defined( $sampleInfo{ $scriptParameters{'familyID'} }{$sampleIDs[$sampleidCounter]}{'exomeTargetBed'} ) ) { #Capture kit check
		    $sampleInfo{ $scriptParameters{'familyID'} }{ $sampleIDs[$sampleidCounter] }{'exomeTargetBed'} =~ s/GenomeReferenceSource/$humanGenomeReferenceSource/; #Replace with Refseq genome or Ensembl genome
		    $sampleInfo{ $scriptParameters{'familyID'} }{ $sampleIDs[$sampleidCounter] }{'exomeTargetBed'} =~ s/Version/$humanGenomeReferenceVersion/; #Replace with actual version 
		    $sampleInfo{ $scriptParameters{'familyID'} }{ $sampleIDs[$sampleidCounter] }{'exomeTargetBed'} =~ s/ChromosomePrefix/$humanGenomeRefereceChromosomePrefix/; #Replace with chromosome prefix
		    $scriptParameters{ $sampleIDs[$sampleidCounter] }{'exomeTargetBed'} = $sampleInfo{ $scriptParameters{'familyID'} }{ $sampleIDs[$sampleidCounter] }{'exomeTargetBed'}; #Add to enable recreation of cmd line later
		}
		else {
		    print STDERR "\nCould not find a capture kit entry for sample: ".$sampleIDs[$sampleidCounter], "\n";
		    print STDERR "\nYou must supply a capture kit bed file when running 'pCC_Bedc'\n";
		    $uncorrectCaptureCounter++;
		}
	    }
	    else { #No capture kit information   
		print STDERR "\nYou must supply a capture kit bed file when running 'pCC_Bedc'\n";
		print STDERR "\n";
		die $USAGE;
	    }
	}
	if ($uncorrectCaptureCounter > 0) { #If lacking or not supported
	    print STDERR "\nChange/add capture kit record in pedigree file: ".$scriptParameters{'pedigreeFile'}, "\n";
	    print STDERR "List of pedigree supported capture kits records:\n\n";
	    print STDERR "Pedigree record", "\t", "Capture kit BED-file\n";
	    for my $supportedCaptureKit (keys %supportedCaptureKits) {
		print STDERR $supportedCaptureKit, "\t", $supportedCaptureKits{$supportedCaptureKit}, "\n";
	    }	    
	    print STDERR "\n";
	    die $USAGE;
	}
	for (my $sampleidCounter=0;$sampleidCounter<scalar(@sampleIDs);$sampleidCounter++) { #Check that all samples in pedigree have a capture kit file
	    unless (-e $scriptParameters{'referencesDir'}."/".$sampleInfo{$scriptParameters{'familyID'}}{$sampleIDs[$sampleidCounter]}{'exomeTargetBed'} ) { #Check existence of capture file in supplied reference dir
		print STDERR "\nCould not find capture kit target BED-file: ".$scriptParameters{'referencesDir'}."/".$sampleInfo{$scriptParameters{'familyID'}}{$sampleIDs[$sampleidCounter]}{'exomeTargetBed'}, "\n\n";
		die $USAGE;
	    }	
	    if ($sampleInfo{$scriptParameters{'familyID'}}{$sampleIDs[0]}{'exomeTargetBed'} eq $sampleInfo{$scriptParameters{'familyID'}}{$sampleIDs[$sampleidCounter]}{'exomeTargetBed'}) { #Enable later print of cmd to be executable if all IDN have been sequenced using the same capture kit
		$identicalCaptureBedCounter++;
	    }
	} 	
    }
}
else {
    for (my $sampleidCounter=0;$sampleidCounter<scalar(@sampleIDs);$sampleidCounter++) { #Add capture kit to all samples
	$sampleInfo{ $scriptParameters{'familyID'} }{ $sampleIDs[$sampleidCounter] }{'exomeTargetBed'} = $exomeTargetBed; #Add capture target BED-file to sampleInfo info to enable individal adjusted capture calculation for each family member
#Check for file existence
	unless (-e $scriptParameters{'referencesDir'}."/".$sampleInfo{ $scriptParameters{'familyID'} }{ $sampleIDs[$sampleidCounter] }{'exomeTargetBed'} ) { #Check for capture file in supplied reference dir
	    print STDERR "\nCould not find capture kit target BED-file: ".$scriptParameters{'referencesDir'}."/".$sampleInfo{ $scriptParameters{'familyID'} }{ $sampleIDs[$sampleidCounter] }{'exomeTargetBed'}, "\n\n";
	    die $USAGE;		
	}
	$scriptParameters{ $sampleIDs[$sampleidCounter] }{'exomeTargetBed'} = $exomeTargetBed; #Add to enable recreation of cmd line later
    }
}

#
#exomeTargetBedInfileList
#
if ( $exomeTargetBedInfileList eq 0 ) {
    
    if ( ($scriptParameters{'pCalculateCoverage'} > 0 ) && ($scriptParameters{'pPicardToolsCalculateHSMetrics'} > 0) ) {
	my $uncorrectCaptureCounter = 0; #Track no entries or wrong format entry in pedigree file
	for (my $sampleidCounter=0;$sampleidCounter<scalar(@sampleIDs);$sampleidCounter++) { #Check all samples
	    if ( defined($scriptParameters{ $sampleIDs[$sampleidCounter] }{'exomeTargetBedInfileList'}) ) { #Input from config file - transfer to sampleInfo
		$sampleInfo{$scriptParameters{'familyID'}}{$sampleIDs[$sampleidCounter]}{'exomeTargetBedInfileList'} = $scriptParameters{ $sampleIDs[$sampleidCounter] }{'exomeTargetBedInfileList'};
	    }	
	    elsif ($scriptParameters{'environmentUppmax'} == 1) {		
		
		if ( defined( $sampleInfo{$scriptParameters{'familyID'} }{$sampleIDs[$sampleidCounter]}{'exomeTargetBedInfileList'} ) ) { #Capture kit check
		    #Add actual run info to the capture kit file names
		    $sampleInfo{ $scriptParameters{'familyID'} }{ $sampleIDs[$sampleidCounter] }{'exomeTargetBedInfileList'} =~ s/GenomeReferenceSource/$humanGenomeReferenceSource/; #Replace with Refseq genome or Ensembl genome
		    $sampleInfo{ $scriptParameters{'familyID'} }{ $sampleIDs[$sampleidCounter] }{'exomeTargetBedInfileList'} =~ s/Version/$humanGenomeReferenceVersion/; #Replace with actual version
		    $sampleInfo{ $scriptParameters{'familyID'} }{ $sampleIDs[$sampleidCounter] }{'exomeTargetBedInfileList'} =~ s/ChromosomePrefix/$humanGenomeRefereceChromosomePrefix/; #Replace with chromosome prefix
		    $scriptParameters{ $sampleIDs[$sampleidCounter] }{'exomeTargetBedInfileList'} = $sampleInfo{ $scriptParameters{'familyID'} }{ $sampleIDs[$sampleidCounter] }{'exomeTargetBedInfileList'}; #Add to enable recreation of cmd line later
		}
		else {
		    print STDERR "\nCould not find a capture kit entry for sample: ".$sampleIDs[$sampleidCounter], "\n";
		    print STDERR "\nYou must supply a capture kit target bed infile_list when running 'pCCE_PicHS'\n";
		    $uncorrectCaptureCounter++;
		}
	    }
	    else { #Lacking capture kit information   
		print STDERR "\nYou must supply a capture kit target bed infile_list when running 'pCCE_PicHS'\n";
		print STDERR "\n";
		die $USAGE;
	    }
	}
	if ($uncorrectCaptureCounter > 0) { #If lacking or not supported
	    print STDERR "\nChange/add capture kit record in pedigree file: ".$scriptParameters{'pedigreeFile'}, "\n";
	    print STDERR "List of pedigree supported capture kits records:\n\n";
	    print STDERR "Pedigree record", "\t", "Capture kit BED-file\n";
	    for my $supportedCaptureKit (keys %supportedCaptureKits) {
		print STDERR $supportedCaptureKit, "\t", $supportedCaptureKits{$supportedCaptureKit}, "\n";
	    }	 
	    print "\n";
	    die $USAGE;
	}
	for (my $sampleidCounter=0;$sampleidCounter<scalar(@sampleIDs);$sampleidCounter++) { #Check all samples
	    unless (-e $scriptParameters{'referencesDir'}."/".$sampleInfo{ $scriptParameters{'familyID'} }{ $sampleIDs[$sampleidCounter] }{'exomeTargetBedInfileList'} ) { #Check existence of capture file in supplied reference dir
		print STDERR "\nCould not find capture kit target bed infile_list file: ".$scriptParameters{'referencesDir'}."/".$sampleInfo{ $scriptParameters{'familyID'} }{ $sampleIDs[$sampleidCounter] }{'exomeTargetBedInfileList'},"\n\n";
		die $USAGE;		
	    }
	    if ($sampleInfo{$scriptParameters{'familyID'}}{$sampleIDs[0]}{'exomeTargetBedInfileList'} eq $sampleInfo{ $scriptParameters{'familyID'} }{ $sampleIDs[$sampleidCounter] }{'exomeTargetBedInfileList'}) { #Enable later print of cmd to be executable if all IDN have been sequenced using the same capture kit
		$identicalCaptureBedIntervalCounter++;
	    }
	}
    }
}   
else {
    for (my $sampleidCounter=0;$sampleidCounter<scalar(@sampleIDs);$sampleidCounter++) { #Add capture kit to all samples using user supplied info
	$sampleInfo{ $scriptParameters{'familyID'} }{ $sampleIDs[$sampleidCounter] }{'exomeTargetBedInfileList'} = $exomeTargetBedInfileList; #Add capture target file to sampleInfo  to enable individal adjusted capture calculation for each family member
	
#Check for file existence 
	unless (-e $scriptParameters{'referencesDir'}."/".$sampleInfo{ $scriptParameters{'familyID'} }{ $sampleIDs[$sampleidCounter] }{'exomeTargetBedInfileList'} ) { #Check for capture file in supplied reference dir
	    print STDERR "\nCould not find capture kit target file: ".$scriptParameters{'referencesDir'}."/".$sampleInfo{ $scriptParameters{'familyID'} }{ $sampleIDs[$sampleidCounter] }{'exomeTargetBedInfileList'}, "\n\n";
	    die $USAGE;		
	}
	$scriptParameters{$sampleIDs[$sampleidCounter]}{'exomeTargetBedInfileList'} = $sampleInfo{ $scriptParameters{'familyID'} }{ $sampleIDs[$sampleidCounter] }{'exomeTargetBedInfileList'}; #Add to enable recreation of cmd line later
    }
}

#
#exomeTargetPaddedBedInfileList
#
if ( $exomeTargetPaddedBedInfileList eq 0 ) {
    if ( ($scriptParameters{'pCalculateCoverage'} > 0) && ($scriptParameters{'pPicardToolsCalculateHSMetrics'} > 0) ) {
	my $uncorrectCaptureCounter = 0; #Track no entries or wrong format entry in pedigree file
	for (my $sampleidCounter=0;$sampleidCounter<scalar(@sampleIDs);$sampleidCounter++) { #Check all samples
	    if ( defined($scriptParameters{ $sampleIDs[$sampleidCounter] }{'exomeTargetPaddedBedInfileList'}) ) { #Input from config file - transfer to sampleInfo
		$sampleInfo{$scriptParameters{'familyID'}}{$sampleIDs[$sampleidCounter]}{'exomeTargetPaddedBedInfileList'} = $scriptParameters{ $sampleIDs[$sampleidCounter] }{'exomeTargetPaddedBedInfileList'};
	    }	
	    elsif ($scriptParameters{'environmentUppmax'} == 1) {		
		if ( defined( $sampleInfo{$scriptParameters{'familyID'} }{$sampleIDs[$sampleidCounter]}{'exomeTargetPaddedBedInfileList'} ) ) { #Capture kit check
		    #Add actual run info to the capture kit file names
		    $sampleInfo{ $scriptParameters{'familyID'} }{ $sampleIDs[$sampleidCounter] }{'exomeTargetPaddedBedInfileList'} =~ s/GenomeReferenceSource/$humanGenomeReferenceSource/; #Replace with Refseq genome or Ensembl genome
		    $sampleInfo{ $scriptParameters{'familyID'} }{ $sampleIDs[$sampleidCounter] }{'exomeTargetPaddedBedInfileList'} =~ s/Version/$humanGenomeReferenceVersion/; #Replace with actual version
		    $sampleInfo{ $scriptParameters{'familyID'} }{ $sampleIDs[$sampleidCounter] }{'exomeTargetPaddedBedInfileList'} =~ s/ChromosomePrefix/$humanGenomeRefereceChromosomePrefix/; #Replace with chromosome prefix
		    $scriptParameters{ $sampleIDs[$sampleidCounter] }{'exomeTargetPaddedBedInfileList'} = $sampleInfo{ $scriptParameters{'familyID'} }{ $sampleIDs[$sampleidCounter] }{'exomeTargetPaddedBedInfileList'}; #Add to enable recreation of cmd line later
		}
		else {
		    print STDERR "\nCould not find a capture kit target bed padded infile_list entry for sample: ".$sampleIDs[$sampleidCounter], "\n";
		    print STDERR "\nYou must supply a capture kit target bed padded infile_list when running 'pCCE_PicHS'\n";
		    $uncorrectCaptureCounter++;
		}
	    }
	    else { #Lacking capture kit information   
		print STDERR "\nYou must supply a capture kit target bed padded infile_list when running 'pCCE_PicHS'\n";
		print STDERR "\n";
		die $USAGE;
	    }
	}
	if ($uncorrectCaptureCounter > 0) { #If lacking or not supported
	    print STDERR "\nChange/add capture kit record in pedigree file: ".$scriptParameters{'pedigreeFile'}, "\n";
	    print STDERR "List of pedigree supported capture kits records:\n\n";
	    print STDERR "Pedigree record", "\t", "Capture kit BED-file\n";
	    for my $supportedCaptureKit (keys %supportedCaptureKits) {
		print STDERR $supportedCaptureKit, "\t", $supportedCaptureKits{$supportedCaptureKit}, "\n";
	    }	 
	    print "\n";
	    die $USAGE;
	}
	for (my $sampleidCounter=0;$sampleidCounter<scalar(@sampleIDs);$sampleidCounter++) { #Check all samples
	    unless (-e $scriptParameters{'referencesDir'}."/".$sampleInfo{ $scriptParameters{'familyID'} }{ $sampleIDs[$sampleidCounter] }{'exomeTargetPaddedBedInfileList'} ) { #Check existence of capture file in supplied reference dir
		print STDERR "\nCould not find padded capture kit target file: ".$scriptParameters{'referencesDir'}."/".$sampleInfo{ $scriptParameters{'familyID'} }{ $sampleIDs[$sampleidCounter] }{'exomeTargetPaddedBedInfileList'},"\n\n";
		die $USAGE;		
	    }
	    if ($sampleInfo{$scriptParameters{'familyID'}}{$sampleIDs[0]}{'exomeTargetPaddedBedInfileList'} eq $sampleInfo{ $scriptParameters{'familyID'} }{ $sampleIDs[$sampleidCounter] }{'exomeTargetPaddedBedInfileList'}) { #Enable later print of cmd to be executable if all IDN have been sequenced using the same capture kit
		$identicalCaptureBedPaddedIntervalCounter++;
	    }
	}
    }
}   
else {
    
    for (my $sampleidCounter=0;$sampleidCounter<scalar(@sampleIDs);$sampleidCounter++) { #Add capture kit to all samples using user supplied info
	$sampleInfo{ $scriptParameters{'familyID'} }{ $sampleIDs[$sampleidCounter] }{'exomeTargetPaddedBedInfileList'} = $exomeTargetPaddedBedInfileList; #Add padded capture target file to sampleInfo  to enable individal adjusted capture calculation for each family member
	
#Check for file existence 
	unless (-e $scriptParameters{'referencesDir'}."/".$sampleInfo{ $scriptParameters{'familyID'} }{ $sampleIDs[$sampleidCounter] }{'exomeTargetPaddedBedInfileList'} ) { #Check for capture file in supplied reference dir
	    print STDERR "\nCould not find padded capture kit target file: ".$scriptParameters{'referencesDir'}."/".$sampleInfo{ $scriptParameters{'familyID'} }{ $sampleIDs[$sampleidCounter] }{'exomeTargetPaddedBedInfileList'}, "\n\n";
	    die $USAGE;		
	}
	$scriptParameters{$sampleIDs[$sampleidCounter]}{'exomeTargetPaddedBedInfileList'} = $sampleInfo{ $scriptParameters{'familyID'} }{ $sampleIDs[$sampleidCounter] }{'exomeTargetPaddedBedInfileList'}; #Add to enable recreation of cmd line later
    }
}

#maximumCores
if( $maximumCores == 0) { #No input from cmd
    if ( defined($scriptParameters{'maximumCores'}) ) { #Input from config file - do nothing
    }
    else {
	$scriptParameters{'maximumCores'} = 8; #Default
    }
}
else {
   $scriptParameters{'maximumCores'} = $maximumCores;
}

#configFile
if( $configFile eq 0) { #No input from cmd
    if ( defined($scriptParameters{'configFile'}) ) { #Input from config file - do nothing
    }
    else {
	$scriptParameters{'configFile'} = 0; #Default
    }
}
else {
    $scriptParameters{'configFile'} = $configFile;
}

if( $writeConfigFile eq 0) { #No input from cmd
    if ( defined($scriptParameters{'writeConfigFile'}) ) { #Input from config file - do nothing
    }
    else {
	$scriptParameters{'writeConfigFile'} = 0; #Default
    }
}
else {
    $scriptParameters{'writeConfigFile'} = $writeConfigFile;
}
if ( $scriptParameters{'writeConfigFile'} ne 0) { #Write config file for family
    open (YAML, '>', $scriptParameters{'writeConfigFile'}) or die "can't open ".$scriptParameters{'writeConfigFile'}.": $!\n";
    print YAML Dump(%scriptParameters), "\n";
    close (YAML);
}

###
#Creates master_logg for the master script 
###
`mkdir -p $scriptParameters{'outDataDir'}/$scriptParameters{'familyID'}/master_logg;`; #Creates the master_logg dir
my ($base,$script) = (`date +%Y%m%d`,`basename $0`); #Catches current date and script name
chomp($base,$script); #Remove \n;
my $masterLoggName = $scriptParameters{'outDataDir'}."/".$scriptParameters{'familyID'}."/master_logg/".$script."_".$base.".txt"; #concatenates master_logg filename
open (MASTERL, ">>".$masterLoggName) or die "Can't write to ".$masterLoggName.": $!\n"; #Open file masterLogg
#Add parameters
print MASTERL "\n".$script." "; #Adds script name to recontruct command line
WriteCMDMasterLogg();

print STDOUT "\nScript parameters and info from $script are saved in file: ".$masterLoggName, "\n";

@inFilesDir = split(/,/,$scriptParameters{'inFilesDir'}); #Enables comma separated indir(s)
@picardToolsMergePreviousFiles = split(/,/,join(',',@picardToolsMergePreviousFiles)); #Enables comma separated list of previously generated _sorted.bam files

for (my $inputDirectoryCounter=0;$inputDirectoryCounter<scalar(@inFilesDir);$inputDirectoryCounter++) { #Collects inputfiles
    
    my @infiles = `cd $inFilesDir[ $inputDirectoryCounter ];ls *.fastq*;`; #cd to input dir and collect fastq files and fastq.gz files
   
    print STDOUT "\nReads from Platform", "\n";print MASTERL "\nReads from Platform", "\n";
    print STDOUT "\nSample ID\t".$sampleIDs[$inputDirectoryCounter],"\n";print MASTERL "\nSample ID\t".$sampleIDs[$inputDirectoryCounter],"\n";
    print STDOUT "Inputfiles\n",@ { $infiles{ $sampleIDs[$inputDirectoryCounter] }  =[@infiles] }, "\n"; #hash with sample id as key and inputfiles in dir as array 
    print MASTERL "Inputfiles\n",@ { $infiles{ $sampleIDs[$inputDirectoryCounter] }  =[@infiles] }, "\n";
    $indirpath{$sampleIDs[$inputDirectoryCounter]} = $inFilesDir[ $inputDirectoryCounter ];  #Catch inputdir path
    chomp(@infiles);    #Remove newline from every entry in array
    $infiles{ $sampleIDs[$inputDirectoryCounter] }  =[@infiles]; #Reload files into hash (kept above newline just for print STDOUT)
}
close(MASTERL);

#Set chr prefix and chromosome names depending on reference used
if ($scriptParameters{'humanGenomeReference'}=~/hg\d+/) { #Refseq - prefix and M
    @chromosomes = ("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM"); #Chr for filtering of bam file
}
elsif ($scriptParameters{'humanGenomeReference'}=~/GRCh\d+/) { #Ensembl - no prefix and MT
    @chromosomes = ("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT"); #Chr for filtering of bam file
}

my $uncompressedFileSwitch = InfilesReFormat(); #Required to format infiles correctly for subsequent input into aligners

###
#MAIN
###

open (MASTERL, ">>".$masterLoggName) or die "Can't write to ".$masterLoggName.": $!\n"; #Open file run logg

if ( ($scriptParameters{'pGZip'} > 0) && ($uncompressedFileSwitch eq 1) ) { #GZip of fastq files

    print STDOUT "\nGZip for fastq files", "\n";print MASTERL "\nGZip for fastq files", "\n";

    for (my $sampleidCounter=0;$sampleidCounter<scalar(@sampleIDs);$sampleidCounter++) {  

	for (my $infileCounter=0;$infileCounter<scalar( @{ $infiles{$sampleIDs[$sampleidCounter]} });$infileCounter++) { #To determine which sampleID had the uncompressed files
	    
	    if ($infiles{$sampleIDs[$sampleidCounter]}[$infileCounter] =~/.fastq$/) {
	
		GZipfastq($sampleIDs[$sampleidCounter]);
		last; #Return to sampleID loop i.e. only call subroutine GZipfastq once per sampleID
	    }
	}
    }
}

if ($scriptParameters{'pFastQC'} > 0) { #Run FastQC
    
    print STDOUT "\nFastQC", "\n";print MASTERL "\nFastQC", "\n";
    
    for (my $sampleidCounter=0;$sampleidCounter<scalar(@sampleIDs);$sampleidCounter++) {  
	
	FastQC($sampleIDs[$sampleidCounter]);	
    }
}

if ($scriptParameters{'pMosaikBuild'} > 0) { #Run MosaikBuild
    
    print STDOUT "\nMosaikBuild", "\n";print MASTERL "\nMosaikBuild", "\n";

    for (my $sampleidCounter=0;$sampleidCounter<scalar(@sampleIDs);$sampleidCounter++) {  
	
	MosaikBuild($sampleIDs[$sampleidCounter]);	
    }
}


if ($scriptParameters{'pMosaikAlign'} > 0) { #Run MosaikAlign
    
    print STDOUT "\nMosaikAlign", "\n"; print MASTERL "\nMosaikAlign", "\n";

    for (my $sampleidCounter=0;$sampleidCounter<scalar(@sampleIDs);$sampleidCounter++) {  
	
	MosaikAlign($sampleIDs[$sampleidCounter]);	
    }
}

if ($scriptParameters{'pBwaAln'} > 0) { #Run BWA Aln
    
    print STDOUT "\nBWA Aln", "\n";print MASTERL "\nBWA Aln", "\n";
    
    for (my $sampleidCounter=0;$sampleidCounter<scalar(@sampleIDs);$sampleidCounter++) {  
	
	BWA_Aln($sampleIDs[$sampleidCounter]);	
    }    
}

if ($scriptParameters{'pBwaSampe'} > 0) { #Run BWA Sampe
    
    print STDOUT "\nBWA Sampe", "\n";print MASTERL "\nBWA Sampe", "\n";
    
    for (my $sampleidCounter=0;$sampleidCounter<scalar(@sampleIDs);$sampleidCounter++) {  
	
	BWA_Sampe($sampleIDs[$sampleidCounter]);
    }
}

if ($scriptParameters{'pSamToolsSort'} > 0) { #Run samtools Sort and Index

    print STDOUT "\nSamTools sort & index", "\n";

    for (my $sampleidCounter=0;$sampleidCounter<scalar(@sampleIDs);$sampleidCounter++) {  
    
	SamToolsSortIndex($sampleIDs[$sampleidCounter], $scriptParameters{'aligner'});
	
    }
}

if ($scriptParameters{'pPicardToolsMerge'} > 0) { #Run picardtools merge (Requires sorted bam files)

    print STDOUT "\nPicardTool MergeSamFiles", "\n";

    for (my $sampleidCounter=0;$sampleidCounter<scalar(@sampleIDs);$sampleidCounter++) {  

	PicardToolsMerge($sampleIDs[$sampleidCounter], $scriptParameters{'aligner'});	
    }
}

if ( $scriptParameters{'pPicardToolsMarkduplicates'} > 0 ) { #PicardTools MarkDuplicates

    print STDOUT "\nPicardTools MarkDuplicates", "\n";print MASTERL "\nPicardTools MarkDuplicates", "\n";

    for (my $sampleidCounter=0;$sampleidCounter<scalar(@sampleIDs);$sampleidCounter++) {  
    
	PicardToolsMarkDuplicates($sampleIDs[$sampleidCounter], $scriptParameters{'aligner'});	
    }
}

if ( $scriptParameters{'pCalculateCoverage'} > 0 ) { #Run GenomeCoverageBED, qaCompute (Paul Costea), Picard (CollectAlignmentSummaryMetrics, CalculateHsMetrics)
    
    print STDOUT "\nCalculate Coverage", "\n";print MASTERL "\nCalculate Coverage", "\n";    
    
    for (my $sampleidCounter=0;$sampleidCounter<scalar(@sampleIDs);$sampleidCounter++) {  

	CalculateCoverage($sampleIDs[$sampleidCounter], $scriptParameters{'aligner'});
    }
}

if ( $scriptParameters{'pRCovPlots'} > 1 ) { #Run Rcovplot scripts   
    print STDOUT "\nRCovPlots", "\n";print MASTERL "\nRCovPlots", "\n";	

    for (my $sampleidCounter=0;$sampleidCounter<scalar(@sampleIDs);$sampleidCounter++) {  
	
	RCoveragePlots($sampleIDs[$sampleidCounter], $scriptParameters{'aligner'});	
    }
}

if ($scriptParameters{'pRemovalRedundantFiles'} > 1) { #Sbatch generation of removal of alignment files
    
    print STDOUT "\nRemoval of alignment files", "\n"; print MASTERL "\nRemoval of alignment files", "\n";
    
    for (my $sampleidCounter=0;$sampleidCounter<scalar(@sampleIDs);$sampleidCounter++) {  
	
	RemoveRedundantFiles($sampleIDs[$sampleidCounter], $scriptParameters{'aligner'});	
    }
}

close(MASTERL); #Close Master_logg file

######################
###Sub Routines#######
######################

sub RemoveRedundantFiles {
#Generates a sbatch script, which removes some alignment files.
    
    my $sampleID = $_[0]; 
    my $aligner = $_[1];
    
    `mkdir -p $scriptParameters{'outDataDir'}/$sampleID/$aligner/info;`; #Creates the aligner and info data file directory
    `mkdir -p $scriptParameters{'outScriptDir'}/$sampleID/$aligner;`; #Creates the aligner script directory
    if ( $scriptParameters{'pRemovalRedundantFiles'} == 1 ) {
	$filename = $scriptParameters{'outScriptDir'}."/".$sampleID."/".$aligner."/removeRedundantFiles_".$sampleID.".";
    }
    elsif ( $scriptParameters{'pRemovalRedundantFiles'} == 2 ) { #Dry run
	$filename = $scriptParameters{'outScriptDir'}."/".$sampleID."/".$aligner."/dry_run_removeRedundantFiles_".$sampleID.".";
	print STDOUT "Dry Run:\n";print MASTERL  "Dry Run:\n";
    }

    Checkfnexists($filename, $fnend);

#Info and Logg
    print STDOUT "Creating sbatch script RemoveRedundantFiles and writing script file(s) to: ".$filename, "\n";print MASTERL "Creating sbatch script RemoveRedundantFiles and writing script file(s) to: ".$filename, "\n";
    print STDOUT "Sbatch script RemoveRedundantFiles data files will be removed in: ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner, "\n";print MASTERL "Sbatch script RemoveRedundantFiles data files will be removed in: ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner, "\n";

    open (REM, ">".$filename) or die "Can't write to ".$filename.": $!\n";
    
    print REM "#! /bin/bash -l", "\n";
    print REM "#SBATCH -A ".$scriptParameters{'projectID'}, "\n";
    print REM "#SBATCH -n 1", "\n";
    print REM "#SBATCH -C thin", "\n";
    print REM "#SBATCH -t 00:15:00", "\n";
    print REM "#SBATCH -J REM_".$sampleID, "\n";
    if ( $scriptParameters{'pRemovalRedundantFiles'} == 1 ) {
	print REM "#SBATCH -e ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/info/rem_".$sampleID.".".$fnt.".stderr.txt", "\n";
	print REM "#SBATCH -o ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/info/rem_".$sampleID.".".$fnt.".stdout.txt", "\n";
    }
    elsif ( $scriptParameters{'pRemovalRedundantFiles'} == 2 ) { #Dry run
	print REM "#SBATCH -e ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/info/rem_".$sampleID.".".$fnt.".stderr.txt", "\n";
	print REM "#SBATCH -o ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/info/rem_".$sampleID.".".$fnt.".stdout.txt", "\n";
    }
    unless ($scriptParameters{'email'} eq 0) {
	print REM "#SBATCH --mail-type=END", "\n";
	print REM "#SBATCH --mail-type=FAIL", "\n";
	print REM "#SBATCH --mail-user=".$scriptParameters{'email'}, "\n\n";
	
    }
    print REM 'echo "Running on: $(hostname)"',"\n\n";

    print REM "cd ";
    print REM $scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner, "\n\n";

    #print REM "#Samples", "\n\n";
    #print REM 'inSampleDir="',"$scriptParameters{'outDataDir'}/$sampleID/$aligner", '"', "\n\n";
    my $inSampleDirectory = $scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner;
###
#Remove Mosaik files
###
    if ( ($scriptParameters{'pMosaikBuild'} > 0) || ($scriptParameters{'pMosaikAlign'} > 0) || ($scriptParameters{'aligner'} eq "mosaik") ) {
	for (my $infileCounter=0;$infileCounter < scalar( @{ $InfilesLaneNoEnding{$sampleID} } );$infileCounter++) { #MosaikBuild takes both reads at once
	    
	    my $tempInfile = $InfilesLaneNoEnding{$sampleID}[$infileCounter]; 
	    print REM "rm ";
	    print REM $inSampleDirectory."/".$tempInfile.".dat", "\n\n"; #MosaikBuild
	    
	    print REM "rm ";
	    print REM $inSampleDirectory."/".$tempInfile.".stat", "\n\n"; #MosaikAlign Stats

	    print REM "rm ";
	    print REM $inSampleDirectory."/".$tempInfile.".bam", "\n\n"; #MosaikAlign

	    #print REM "rm ";
	    #print REM $inSampleDirectory."/".$tempInfile.".multiple.bam", "\n\n"; #MosaikAlign Multiple

	    #print REM "rm ";
	    #print REm $inSampleDirectory."/".$tempInfile."_sorted.bam", "\n\n"; #MosaikAlign/samtools

	    #print REM "rm ";
	    #print REM $inSampleDirectory."/".$tempInfile."_sorted.bam.bai", "\n\n"; #MosaikAlign/samtools index

	    print REM "rm ";
	    print REM $inSampleDirectory."/coverageReport/".$tempInfile."_sorted_pmd_coverageBed", "\n\n"; #Coverage of features in BED-file

	    print REM "rm ";
	    print REM $inSampleDirectory."/coverageReport/".$tempInfile."_sorted_pmd_coverageBed_hist", "\n\n"; #bedtools histogram of BED-file
	    
	    print REM "rm ";
	    print REM $inSampleDirectory."/coverageReport/".$tempInfile."_sorted_pmd_rmdup_coverageBed_hist", "\n\n"; #bedtools histogram of BED-file
	    
	    if ( ($scriptParameters{'pPicardToolsMerge'} == 1) && ( $infileCounter == scalar( @{ $InfilesLaneNoEnding{$sampleID} } )-1 ) && ($infileCounter >= 1) ) { #If merged, last file and that there exits at least 2 files
		print REM "rm ";
		print REM $inSampleDirectory."/coverageReport/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,"_sorted_merged_pmd_coverageBed", "\n\n"; #Coverage of features in BED-file

		print REM "rm ";
		print REM $inSampleDirectory."/coverageReport/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,"_sorted_merged_pmd_coverageBed_hist", "\n\n"; #bedtools histogram of BED-file

		print REM "rm ";
		print REM $inSampleDirectory."/coverageReport/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,"_sorted_merged_pmd_rmdup_coverageBed_hist", "\n\n"; #bedtools histogram of BED-file

		#print REM "rm ";
		#print REM $inSampleDirectory."/coverageReport/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,"_sorted_merged_pmd_coverageBed_depth_pos", "\n\n"; #Coverage of features BED-file per position
	    }
	}    
	print REM "rm -rf ", '${inSampleDir}', "/per_chr", "\n\n"; #samtools/GATK (real/recal)
    }
###
#Remove BWA files
###
    if ( ($scriptParameters{'pBwaAln'} > 0) || ($scriptParameters{'pBwaSampe'} >0) || ($scriptParameters{'aligner'} eq "bwa")) {

	for (my $infileCounter=0;$infileCounter < scalar( @{ $InfilesBothStrandsNoEnding{$sampleID} } );$infileCounter++) { #BWA_Aln takes 1 read at a time 
	    
	    my $tempInfile = $InfilesBothStrandsNoEnding{$sampleID}[$infileCounter]; 
	    print REM "rm ";
	    print REM $inSampleDirectory."/".$tempInfile.".sai", "\n\n"; #BWA_Aln
	}
	for (my $infileCounter=0;$infileCounter < scalar( @{ $InfilesLaneNoEnding{$sampleID} } );$infileCounter++) { #BWA_Sampe 
	    
	    my $tempInfile = $InfilesLaneNoEnding{$sampleID}[$infileCounter]; 
	    print REM "rm ";
	    print REM $inSampleDirectory."/".$tempInfile.".bam", "\n\n"; #BWA_Sampe
	}    
    }
    close(REM);
    return;
}

sub RCoveragePlots { 
#Generates sbatch scripts for R scripts:
#1. covplots_genome.R 
#2. covplots_exome.R
#on files generated from calculateCoverage genomeCoverageBED

    my $sampleID = $_[0]; 
    my $aligner = $_[1];
    
    `mkdir -p $scriptParameters{'outDataDir'}/$sampleID/$aligner/info;`; #Creates the aligner folder and info data file directory
    `mkdir -p $scriptParameters{'outDataDir'}/$sampleID/$aligner/coverageReport;`; #Creates the aligner and coverageReport folder
    `mkdir -p $scriptParameters{'outScriptDir'}/$sampleID/$aligner;`; #Creates the aligner script directory
    if ( $scriptParameters{'pRCovPlots'} == 1) {
	$filename = $scriptParameters{'outScriptDir'}."/".$sampleID."/".$aligner."/rCovPlots_".$sampleID.".";
    }
    elsif ( $scriptParameters{'pRCovPlots'} == 2) { #Dry run
	$filename = $scriptParameters{'outScriptDir'}."/".$sampleID."/".$aligner."/dry_run_rCovPlots_".$sampleID.".";
	print STDOUT "Dry Run:\n";print MASTERL  "Dry Run:\n";
    }

    Checkfnexists($filename, $fnend);

#Info and Logg
    print STDOUT "Creating sbatch script RCoveragePlots and writing script file(s) to: ".$filename, "\n";print MASTERL "Creating sbatch script RCoveragePlots and writing script file(s) to: ".$filename, "\n";
    print STDOUT "Sbatch script RCoveragePlots data files will be written to: ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/coverageReport", "\n";print MASTERL "Sbatch script RCoveragePlots data files will be written to: ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/coverageReport", "\n";
    
    open (RCOVP, ">".$filename) or die "Can't write to ".$filename.": $!\n";
    
    print RCOVP "#! /bin/bash -l", "\n";
    print RCOVP "#SBATCH -A ".$scriptParameters{'projectID'}, "\n";
    print RCOVP "#SBATCH -n 1 ", "\n";
    print RCOVP "#SBATCH -C thin", "\n";	
    print RCOVP "#SBATCH -t 01:00:00", "\n"; 
    print RCOVP "#SBATCH -J RCOVPlots_".$sampleID."_".$aligner, "\n";
    if ( $scriptParameters{'pRCovPlots'} == 1) {
	print RCOVP "#SBATCH -e ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/info/rCovPlots_".$sampleID.".".$fnt.".stderr.txt", "\n";
	print RCOVP "#SBATCH -o ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/info/rCovPlots_".$sampleID.".".$fnt.".stdout.txt", "\n";
    }
    elsif ( $scriptParameters{'pRCovPlots'} == 2) { #Dry run
	print RCOVP "#SBATCH -e ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/info/dry_run_rCovPlots_".$sampleID.".".$fnt.".stderr.txt", "\n";
	print RCOVP "#SBATCH -o ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/info/dry_run_rCovPlots_".$sampleID.".".$fnt.".stdout.txt", "\n";
    }

    unless ($scriptParameters{'email'} eq 0) {	
	print RCOVP "#SBATCH --mail-type=END", "\n";
	print RCOVP "#SBATCH --mail-type=FAIL", "\n";
	print RCOVP "#SBATCH --mail-user=".$scriptParameters{'email'}, "\n\n";
    }
    
    print RCOVP 'echo "Running on: $(hostname)"',"\n\n";
    print RCOVP "module load bioinfo-tools", "\n\n"; 
    print RCOVP "module load R/2.12.2", "\n\n";
    print RCOVP "#Samples", "\n";
    
    print RCOVP 'inSampleDir="',"$scriptParameters{'outDataDir'}/$sampleID/$aligner/coverageReport", '"', "\n";
    print RCOVP 'outSampleDir="', "$scriptParameters{'outDataDir'}/$sampleID/$aligner/coverageReport", '"', "\n\n"; 
    my $inSampleDirectory = $scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/coverageReport";
    my $outSampleDirectory = $scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/coverageReport";
    my $fileending;
    if ($scriptParameters{'pPicardToolsMarkduplicates'} == 0) { #Ensure correct file ending
	$fileending = "_sorted";
    }
    else {
	$fileending = "_sorted_pmd";
    }	
    for (my $infileCounter=0;$infileCounter<scalar( @{ $InfilesLaneNoEnding{$sampleID} } );$infileCounter++) { #For all files from MosaikBuild
	
	my $tempInfile = $InfilesLaneNoEnding{$sampleID}[$infileCounter];
	
	if ($scriptParameters{'pGenomeCoverageBED'} > 0) {
	    print RCOVP "Rscript ";
	    print RCOVP $scriptParameters{'inScriptDir'}."/covplots_genome.R ";
	    print RCOVP $inSampleDirectory."/".$tempInfile.$fileending."_genomeCoverageBed "; #InFile
	    print RCOVP $tempInfile." "; #Sample name
	    print RCOVP $scriptParameters{'xCoverage'}." "; #X-axis max scale
	    print RCOVP $outSampleDirectory, "\n\n"; #OutFile
	}
	if ($scriptParameters{'pCoverageBED'} > 0) {
	    print RCOVP "#Prepp indata file to contain only all features\n";
	    print RCOVP "grep ";
	    print RCOVP "^all "; #Prepp indata file to contain only all features
	    print RCOVP $inSampleDirectory."/".$tempInfile.$fileending."_coverageBed_hist "; #InFile
	    print RCOVP "> ".$outSampleDirectory."/".$tempInfile.$fileending."_coverageBed_all_hist", "\n\n"; #OutFile 

	    print RCOVP "Rscript ";
	    print RCOVP $scriptParameters{'inScriptDir'}."/covplots_exome_all.R ";
	    print RCOVP $inSampleDirectory."/".$tempInfile.$fileending."_coverageBed_all_hist "; #InFile
	    print RCOVP $tempInfile." "; #X-axis max scale
	    print RCOVP $outSampleDirectory, "\n\n"; #OutFile
#Duplicates removed
	    print RCOVP "#Duplicates removed\n\n";
	    print RCOVP "#Prepp indata file to contain only all features\n";
	    print RCOVP "grep ";
	    print RCOVP "^all "; #Prepp indata file to contain only all features 
	    print RCOVP $inSampleDirectory."/".$tempInfile.$fileending."_rmdup_coverageBed_hist "; #InFile
	    print RCOVP "> ".$outSampleDirectory."/".$tempInfile.$fileending."_rmdup_coverageBed_all_hist", "\n\n"; #OutFile
 
	    print RCOVP "Rscript ";
	    print RCOVP $scriptParameters{'inScriptDir'}."/covplots_exome_all.R ";
	    print RCOVP $inSampleDirectory."/".$tempInfile.$fileending."_rmdup_coverageBed_all_hist "; #InFile
	    print RCOVP $tempInfile."_rmdup "; #Sample name
	    print RCOVP $outSampleDirectory, "\n\n"; #OutFile	    
	}
    }
    
    if ( ($scriptParameters{'pPicardToolsMerge'} == 1) && (scalar( @{ $InfilesLaneNoEnding{$sampleID} }) > 1) ) { #but only if there is more than one mosaikBuild file
	if ($scriptParameters{'pPicardToolsMarkduplicates'} == 0) { #Ensure correct file ending
	    $fileending = "_sorted_merged";
	}
	else {
	    $fileending = "_sorted_merged_pmd";
	}
	if ($scriptParameters{'pGenomeCoverageBED'} > 0) {
	print RCOVP "Rscript ";
	print RCOVP $scriptParameters{'inScriptDir'}."/covplots_genome.R ";
	print RCOVP $inSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileending."_genomeCoverageBed "; #InFile
	print RCOVP $sampleID."_lanes_", @{ $lanes{$sampleID} } ," "; #Sample name
	print RCOVP $scriptParameters{'xCoverage'}." "; #X-axis max scale
	print RCOVP $outSampleDirectory, "\n\n"; #OutFile
	}
	if ($scriptParameters{'pCoverageBED'} > 0) {	    
	    print RCOVP "#Prepp indata file to contain only all features\n";
	    print RCOVP "grep ";
	    print RCOVP "^all "; #Prepp indata file to contain only all features
	    print RCOVP $inSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileending."_coverageBed_hist "; #InFile
	    print RCOVP "> ".$outSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileending."_coverageBed_all_hist", "\n\n"; #OutFile
	    
	    print RCOVP "Rscript ";
	    print RCOVP $scriptParameters{'inScriptDir'}."/covplots_exome_all.R ";
	    print RCOVP $inSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileending."_coverageBed_all_hist "; #InFile
	    print RCOVP $sampleID."_lanes_", @{ $lanes{$sampleID} } ," "; #Sample name
	    print RCOVP $outSampleDirectory, "\n\n"; #OutFile

	    #Duplicates removed
	    print RCOVP "#Duplicates removed\n\n";
	    print RCOVP "#Prepp indata file to contain only all features\n";
	    print RCOVP "grep ";
	    print RCOVP "^all "; #Prepp indata file to contain only all features 
	    print RCOVP $inSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileending."_rmdup_coverageBed_hist "; #InFile
	    print RCOVP "> ".$outSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileending."_rmdup_coverageBed_all_hist", "\n\n"; #OutFile

	    print RCOVP "Rscript ";
	    print RCOVP $scriptParameters{'inScriptDir'}."/covplots_exome_all.R ";
	    print RCOVP $inSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileending."_rmdup_coverageBed_all_hist "; #InFile
	    print RCOVP $sampleID."_lanes_", @{ $lanes{$sampleID} } ,"_rmdup "; #Sample name
	    print RCOVP $outSampleDirectory, "\n\n"; #OutFile
	}
    }
    if (@picardToolsMergePreviousFiles ) { # Coverage report R plots on files merged this round with merged file from previous round
	if ( ($scriptParameters{'pPicardToolsMerge'} == 0) && ($scriptParameters{'pPicardToolsMarkduplicates'} == 0) ) { #Ensure correct file ending
	    $fileending = "_sorted";
	}
	if ( ($scriptParameters{'pPicardToolsMerge'} > 0) && ($scriptParameters{'pPicardToolsMarkduplicates'} == 0 ) ) {
	    $fileending = "_sorted_merged";
	}
	if ( ($scriptParameters{'pPicardToolsMerge'} == 0) && ($scriptParameters{'pPicardToolsMarkduplicates'} > 0) ) {
	    $fileending = "_sorted_pmd";
	}
	if ( ($scriptParameters{'pPicardToolsMerge'} > 0) && ($scriptParameters{'pPicardToolsMarkduplicates'} > 0) ) {
	    $fileending = "_sorted_merged_pmd";
	}
	for (my $mergeFile=0;$mergeFile<scalar(@picardToolsMergePreviousFiles);$mergeFile++) {
	    
	    if ($picardToolsMergePreviousFiles[$mergeFile] =~ /$sampleID/) { #Look for sampleID in previously generated file to be merged with current run to be able to merge correct files
		if ($picardToolsMergePreviousFiles[$mergeFile] =~ /lane(\d+)|s_(\d+)/) { #Look for lanes_ or lane\d in previously generated file to be merged with current run to be able to extract previous lanes

		    my $mergeLanes; if($1) {$mergeLanes = $1;} else {$mergeLanes = $2;} #Make sure to always supply lanes from previous regexp

		    if ($scriptParameters{'pGenomeCoverageBED'} > 0) {
			print RCOVP "Rscript ";
			print RCOVP $scriptParameters{'inScriptDir'}."/covplots_genome.R ";
			print RCOVP $inSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileending."_genomeCoverageBed "; #InFile
			print RCOVP $sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ," "; #Sample name
			print RCOVP $scriptParameters{'xCoverage'}." "; #X-axis max scale
			print RCOVP $outSampleDirectory, "\n\n"; #OutFile
		    }
		    if ($scriptParameters{'pCoverageBED'} > 0) {
			print RCOVP "grep ";
			print RCOVP "^all "; #Prepp indata file to contain only all features
			print RCOVP $inSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileending."_coverageBed_hist "; #InFile
			print RCOVP "> ".$inSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileending."_coverageBed_all_hist", "\n\n"; #OutFile
			print RCOVP "Rscript ";
			print RCOVP $scriptParameters{'inScriptDir'}."/covplots_exome_all.R ";
			print RCOVP $inSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileending."_coverageBed_all_hist "; #InFile
			print RCOVP $sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ," "; #Sample name
			print RCOVP $outSampleDirectory, "\n\n"; #OutFile

			#Duplicates removed
			print RCOVP "#Duplicates removed\n\n";
			print RCOVP "#Prepp indata file to contain only all features\n";
			print RCOVP "grep ";
			print RCOVP "^all "; #Prepp indata file to contain only all features
			print RCOVP $inSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileending."_rmdup_coverageBed_hist "; #InFile
			print RCOVP "> ".$outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileending."_rmdup_coverageBed_all_hist", "\n\n"; #OutFile
			print RCOVP "Rscript ";
			print RCOVP $scriptParameters{'inScriptDir'}."/covplots_exome_all.R ";
			print RCOVP $inSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileending."rmdup_coverageBed_all_hist "; #InFile
			print RCOVP $sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,"_rmdup "; #Sample name
			print RCOVP $outSampleDirectory, "\n\n"; #OutFile
		    }
		}
	    }
	}
	print RCOVP "wait", "\n\n";
    }
    close(RCOVP);
    if ( $scriptParameters{'pRCovPlots'} == 1) {
	ParallelSampleIDSubmitJob($sampleID,$filename,"all");
    }
    return;
}

sub CalculateCoverage { 
#Generates sbatch scripts and calculates coverage on alignment files (sorted). 
#NOTE:Collect_info.pl collects key metric reference file from .alignment_summary_metrics. If not processed genome build will be missing in key metric file.

    my $sampleID = $_[0]; 
    my $aligner = $_[1]; 
    
    `mkdir -p $scriptParameters{'outDataDir'}/$sampleID/$aligner/info;`; #Creates the aligner folder and info data file directory
    `mkdir -p $scriptParameters{'outDataDir'}/$sampleID/$aligner/coverageReport;`; #Creates the aligner and coverageReport folder
    #`mkdir -p $scriptParameters{'outDataDir'}/$_[2]/$aligner/GATK/candidates/coverageReport;`; #Creates the aligner and coverageReport folder. NOTE: FamilyID
    `mkdir -p $scriptParameters{'outScriptDir'}/$sampleID/$aligner;`; #Creates the aligner script directory
    if ( $scriptParameters{'pCalculateCoverage'} == 1) {
	$filename = $scriptParameters{'outScriptDir'}."/".$sampleID."/".$aligner."/calculate_coverage_".$sampleID.".";
    }
    elsif ( $scriptParameters{'pCalculateCoverage'} == 2) { #Dry run
	$filename = $scriptParameters{'outScriptDir'}."/".$sampleID."/".$aligner."/dry_run_calculate_coverage_".$sampleID.".";
	print STDOUT "Dry Run:\n";print MASTERL  "Dry Run:\n";
    }
    Checkfnexists($filename, $fnend);

#Info and Logg
    print STDOUT "Creating sbatch script Calculate Coverage and writing script file(s) to: ".$filename, "\n";print MASTERL "Creating sbatch script Calculate Coverage and writing script file(s) to: ".$filename, "\n";
    print STDOUT "Sbatch script Calculate Coverage data files will be written to: ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/coverageReport", "\n";print MASTERL "Sbatch script Calculate Coverage data files will be written to: ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/coverageReport", "\n";

    my $time = ceil(3*scalar( @{ $InfilesBothStrandsNoEnding{$sampleID} })); #One full lane on Hiseq takes approx. 2 h to process, round up to nearest full hour.
    
    open (COCCAL, ">".$filename) or die "Can't write to ".$filename.": $!\n";
    
    print COCCAL "#! /bin/bash -l", "\n";
    print COCCAL "#SBATCH -A ".$scriptParameters{'projectID'}, "\n";
    print COCCAL "#SBATCH -p node -n ".$scriptParameters{'maximumCores'}, "\n";
    print COCCAL "#SBATCH -C thin", "\n";

    if ($scriptParameters{'pPicardToolsMerge'} eq 0) {	
	print COCCAL "#SBATCH -t 4:00:00", "\n";	
    }
    else{
	print COCCAL "#SBATCH -t ".$time.":00:00", "\n";	
    }	

    print COCCAL "#SBATCH -J CCov_".$sampleID."_".$aligner, "\n";
    if ( $scriptParameters{'pCalculateCoverage'} == 1) {
	print COCCAL "#SBATCH -e ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/info/calculate_coverage_".$sampleID.".".$fnt.".stderr.txt", "\n";
	print COCCAL "#SBATCH -o ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/info/calculate_coverage_".$sampleID.".".$fnt.".stdout.txt", "\n";
    }
    elsif ( $scriptParameters{'pCalculateCoverage'} == 2) {
	print COCCAL "#SBATCH -e ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/info/dry_run_calculate_coverage_".$sampleID.".".$fnt.".stderr.txt", "\n";
	print COCCAL "#SBATCH -o ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/info/dry_run_calculate_coverage_".$sampleID.".".$fnt.".stdout.txt", "\n";
    }

    unless ($scriptParameters{'email'} eq 0) {
	print COCCAL "#SBATCH --mail-type=END", "\n";
	print COCCAL "#SBATCH --mail-type=FAIL", "\n";
	print COCCAL "#SBATCH --mail-user=".$scriptParameters{'email'}, "\n\n";
    }
    
    print COCCAL 'echo "Running on: $(hostname)"',"\n\n";
   
    my $inSampleDirectory = $scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner;
    my $outSampleDirectory = $scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/coverageReport";
    my $coreCounter=1;
    my $fileEnding;
    if ($scriptParameters{'pPicardToolsMarkduplicates'} == 0) { #Ensure correct file ending
	$fileEnding = "_sorted";
    }
    else {
	$fileEnding = "_sorted_pmd";
    }
    for (my $infileCounter=0;$infileCounter<scalar( @{ $InfilesLaneNoEnding{$sampleID} } );$infileCounter++) { #For all files from MosaikAlign or BWA_Sampe
	
	my $tempInfile = $InfilesLaneNoEnding{$sampleID}[$infileCounter];
	if ($scriptParameters{'pGenomeCoverageBED'} > 0) {
	    print COCCAL "genomeCoverageBed ";
	    print COCCAL "-max ".$scriptParameters{'xCoverage'}." "; #Combine all positions with a depth >= max into a single bin in the histogram.
	    print COCCAL "-ibam ".$inSampleDirectory."/".$tempInfile.$fileEnding.".bam "; #InFile
	    print COCCAL "> ".$outSampleDirectory."/".$tempInfile.$fileEnding."_genomeCoverageBed &", "\n\n"; #outFile
	}
	if ($scriptParameters{'pQaCompute'} > 0) { #Genome coverage calculations
	    print COCCAL "qaCompute ";
	    print COCCAL "-m "; #Compute median coverage
	    print COCCAL "-d "; #Print per-chromosome histogram
	    print COCCAL "-i "; #Silent 
	    print COCCAL "-c ".$scriptParameters{'xCoverage'}." "; #Max depth to calculate coverage on
	    print COCCAL $inSampleDirectory."/".$tempInfile.$fileEnding.".bam "; #InFile
	    print COCCAL $outSampleDirectory."/".$tempInfile.$fileEnding."_qaCompute &", "\n\n"; #OutFile
	}
	if ($scriptParameters{'pPicardToolsCollectMultipleMetrics'} > 0) {
	    print COCCAL "java -Xmx4g -jar ".$scriptParameters{'picardToolsPath'}."/CollectMultipleMetrics.jar ";
	    print COCCAL "INPUT=".$inSampleDirectory."/".$tempInfile.$fileEnding.".bam "; #InFile
	    print COCCAL "OUTPUT=".$outSampleDirectory."/".$tempInfile.$fileEnding." "; #outFile
	    print COCCAL "R=".$scriptParameters{'referencesDir'}."/".$scriptParameters{'humanGenomeReference'}." &", "\n\n"; #Reference file
	}
	if ($scriptParameters{'pPicardToolsCalculateHSMetrics'} > 0) { #Run CalculateHsMetrics (Target BED-file)
	    print COCCAL "java -Xmx4g -jar ".$scriptParameters{'picardToolsPath'}."/CalculateHsMetrics.jar ";
	    print COCCAL "INPUT=".$inSampleDirectory."/".$tempInfile.$fileEnding.".bam "; #InFile
	    print COCCAL "OUTPUT=".$outSampleDirectory."/".$tempInfile.$fileEnding."_CalculateHsMetrics "; #OutFile
	    print COCCAL "REFERENCE_SEQUENCE=".$scriptParameters{'referencesDir'}."/".$scriptParameters{'humanGenomeReference'}." "; #Reference file
	    print COCCAL "BAIT_INTERVALS=".$scriptParameters{'referencesDir'}."/".$sampleInfo{ $scriptParameters{'familyID'} }{$sampleID}{'exomeTargetPaddedBedInfileList'}." "; #Capture kit padded target infile_list file
	    print COCCAL "TARGET_INTERVALS=".$scriptParameters{'referencesDir'}."/".$sampleInfo{ $scriptParameters{'familyID'} }{$sampleID}{'exomeTargetBedInfileList'}." &", "\n\n"; #Capture kit target infile_list file 
	}
	if ($scriptParameters{'pCoverageBED'} > 0) { #Run coverageBed (Target BED-file)
	    print COCCAL "coverageBed ";
	    print COCCAL "-hist "; #Report a histogram of coverage for each feature in B as well as a summary histogram for _all_ features in B.
	    print COCCAL "-abam ".$inSampleDirectory."/".$tempInfile.$fileEnding.".bam "; #InFile in BAM format
	    print COCCAL "-b ".$scriptParameters{'referencesDir'}."/".$sampleInfo{ $scriptParameters{'familyID'} }{$sampleID}{'exomeTargetBed'}." "; #InFile
	    print COCCAL "> ".$outSampleDirectory."/".$tempInfile.$fileEnding."_coverageBed_hist &", "\n\n"; #OutFile
	    #Remove PCR and Optical duplicates
	    print COCCAL "samtools view ";
	    print COCCAL "-F 0x400 "; #Skip alignments where read is PCR or optical duplicate
	    print COCCAL "-b ".$inSampleDirectory."/".$tempInfile.$fileEnding.".bam "; #InFile
	    print COCCAL "| coverageBed "; #Note "|"
	    print COCCAL "-hist "; #Report a histogram of coverage for each feature in B as well as a summary histogram for _all_ features in B.
	    print COCCAL "-abam stdin "; #InStream
	    print COCCAL "-b ".$scriptParameters{'referencesDir'}."/".$sampleInfo{ $scriptParameters{'familyID'} }{$sampleID}{'exomeTargetBed'}." "; #InFile
	    print COCCAL "> ".$outSampleDirectory."/".$tempInfile.$fileEnding."_rmdup_coverageBed_hist &", "\n\n"; #OutFile
	    
	}
	print COCCAL "wait", "\n\n";
    }

    if ( ($scriptParameters{'pPicardToolsMerge'} > 0) && (scalar( @{ $InfilesLaneNoEnding{$sampleID} }) > 1) ) { #but only if there is more than one mosaikBuild/BWA_Aln file per sample ID (Sanity check)
	if ($scriptParameters{'pPicardToolsMarkduplicates'} == 0) { #Ensure correct file ending
	    $fileEnding = "_sorted_merged";
	}
	else {
	    $fileEnding = "_sorted_merged_pmd";
	}	
	if ($scriptParameters{'pGenomeCoverageBED'} > 0) {
	    print COCCAL "genomeCoverageBed ";
	    print COCCAL "-max ".$scriptParameters{'xCoverage'}." ";  #Combine all positions with a depth >= max into a single bin in the histogram.
	    print COCCAL "-ibam ".$inSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding.".bam "; #InFile
	    print COCCAL "> ".$outSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding."_genomeCoverageBed &", "\n\n"; #OutFile
	}
	if ($scriptParameters{'pQaCompute'} > 0) { #Genome coverage calculations
	    print COCCAL "qaCompute ";
	    print COCCAL "-m "; #Compute median coverage
	    print COCCAL "-d "; #Print per-chromosome histogram
	    print COCCAL "-i "; #Silent 
	    print COCCAL "-c ".$scriptParameters{'xCoverage'}." "; #Max depth to calculate coverage on
	    print COCCAL $inSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding.".bam "; #InFile
	    print COCCAL $outSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding."_qaCompute &", "\n\n"; #OutFile
	}
	if ($scriptParameters{'pPicardToolsCollectMultipleMetrics'} > 0) {
	    print COCCAL "java -Xmx4g -jar ".$scriptParameters{'picardToolsPath'}."/CollectMultipleMetrics.jar ";
	    print COCCAL "INPUT=".$inSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding.".bam "; #InFile
	    print COCCAL "OUTPUT=".$outSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding." "; #OutFile
	    print COCCAL "R=".$scriptParameters{'referencesDir'}."/".$scriptParameters{'humanGenomeReference'}." &", "\n\n"; #Reference file
	}
	if ($scriptParameters{'pPicardToolsCalculateHSMetrics'} > 0) { #Run CalculateHsMetrics (Target BED-file)
	    print COCCAL "java -Xmx4g -jar ".$scriptParameters{'picardToolsPath'}."/CalculateHsMetrics.jar ";
	    print COCCAL "INPUT=".$inSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding.".bam "; #inFile
	    print COCCAL "OUTPUT=".$outSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding."_CalculateHsMetrics "; #OutFile
	    print COCCAL "REFERENCE_SEQUENCE=".$scriptParameters{'referencesDir'}."/".$scriptParameters{'humanGenomeReference'}." "; #Reference file
	    print COCCAL "BAIT_INTERVALS=".$scriptParameters{'referencesDir'}."/".$sampleInfo{ $scriptParameters{'familyID'} }{$sampleID}{'exomeTargetPaddedBedInfileList'}." "; #Capture kit padded target infile_list file
	    print COCCAL "TARGET_INTERVALS=".$scriptParameters{'referencesDir'}."/".$sampleInfo{ $scriptParameters{'familyID'} }{$sampleID}{'exomeTargetBedInfileList'}." &", "\n\n"; #Capture kit target infile_list file
	}
	if ($scriptParameters{'pCoverageBED'} > 0) { #Run coverageBed (Target BED-file)
	    
	    print COCCAL "coverageBed ";
	    print COCCAL "-hist "; #Report a histogram of coverage for each feature in B as well as a summary histogram for _all_ features in B.
	    print COCCAL "-abam ".$inSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding.".bam "; #InFile in BAM format
	    print COCCAL "-b ".$scriptParameters{'referencesDir'}."/".$sampleInfo{ $scriptParameters{'familyID'} }{$sampleID}{'exomeTargetBed'}." "; #InFile
	    print COCCAL "> ".$outSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed_hist &", "\n\n"; #OutFile
	    #Remove PCR and Optical duplicates
	    print COCCAL "samtools view ";
	    print COCCAL "-F 0x400 "; #Skip alignments where read is PCR or optical duplicate
	    print COCCAL "-b ".$inSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding.".bam "; #InFile
	    print COCCAL "| coverageBed "; #Note "|"
	    print COCCAL "-hist "; #Report a histogram of coverage for each feature in B as well as a summary histogram for _all_ features in B.
	    print COCCAL "-abam stdin "; #InStream
	    print COCCAL "-b ".$scriptParameters{'referencesDir'}."/".$sampleInfo{ $scriptParameters{'familyID'} }{$sampleID}{'exomeTargetBed'}." "; #InFile
	    print COCCAL "> ".$outSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding."_rmdup_coverageBed_hist &", "\n\n"; #OutFile
	}
	print COCCAL "wait", "\n\n"; 
    }
    
    if (@picardToolsMergePreviousFiles ) { # Coverage report on files merged this round with merged file from previous round
	if ( ($scriptParameters{'pPicardToolsMerge'} == 0) && ($scriptParameters{'pPicardToolsMarkduplicates'} == 0) ) { #Ensure correct file ending
	    $fileEnding = "_sorted";
	}
	if ( ($scriptParameters{'pPicardToolsMerge'} > 0) && ($scriptParameters{'pPicardToolsMarkduplicates'} == 0 ) ) {
	    $fileEnding = "_sorted_merged";
	}
	if ( ($scriptParameters{'pPicardToolsMerge'} == 0) && ($scriptParameters{'pPicardToolsMarkduplicates'} > 0) ) {
	    $fileEnding = "_sorted_pmd";
	}
	if ( ($scriptParameters{'pPicardToolsMerge'} > 0) && ($scriptParameters{'pPicardToolsMarkduplicates'} > 0) ) {
	    $fileEnding = "_sorted_merged_pmd";
	}
	print COCCAL "wait", "\n\n"; #For previous runs
	for (my $mergeFileCounter=0;$mergeFileCounter<scalar(@picardToolsMergePreviousFiles);$mergeFileCounter++) {
	    
	    if ($picardToolsMergePreviousFiles[$mergeFileCounter] =~ /$sampleID/) { #Look for sampleID in previously generated file to be merged with current run to be able to merge correct files
		if ($picardToolsMergePreviousFiles[$mergeFileCounter] =~ /lane(\d+)|s_(\d+)/) { #Look for lanes_ or lane\d in previously generated file to be merged with current run to be able to extract previous lanes

		    my $mergeLanes; if($1) {$mergeLanes = $1;} else {$mergeLanes = $2;} #Make sure to always supply lanes from previous regexp

		    if ($scriptParameters{'pGenomeCoverageBED'} > 0) {
			print COCCAL "genomeCoverageBed ";
			print COCCAL "-max ".$scriptParameters{'xCoverage'}." "; #Combine all positions with a depth >= max into a single bin in the histogram.
			print COCCAL "-ibam ".$inSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding.".bam "; #InFile
			print COCCAL "> ".$outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding."_genomeCoverageBed &", "\n\n"; #OutFile
		    }
		    if ($scriptParameters{'pQaCompute'} > 0) {
			print COCCAL "qaCompute ";
			print COCCAL "-m "; #Compute median coverage
			print COCCAL "-d "; #Print per-chromosome histogram
			print COCCAL "-i "; #Silent
			print COCCAL "-c ".$scriptParameters{'xCoverage'}." ";
			print COCCAL $inSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding.".bam "; #InFile
			print COCCAL $outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding."_qaCompute &", "\n\n"; #OutFile
		    }
		    if ($scriptParameters{'pPicardToolsCollectMultipleMetrics'} > 0) {
			print COCCAL "java -Xmx4g -jar ".$scriptParameters{'picardToolsPath'}."/CollectMultipleMetrics.jar ";
			print COCCAL "INPUT=".$inSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding.".bam "; #InFile
			print COCCAL "OUTPUT=".$outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding." "; #OutFile
			print COCCAL "R=".$scriptParameters{'referencesDir'}."/".$scriptParameters{'humanGenomeReference'}." &", "\n\n"; #Reference file
		    }
		    if ($scriptParameters{'pPicardToolsCalculateHSMetrics'} > 0) { #Run CalculateHsMetrics (Target BED-file)
			print COCCAL "java -Xmx4g -jar ".$scriptParameters{'picardToolsPath'}."/CalculateHsMetrics.jar ";
			print COCCAL "INPUT=".$inSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding.".bam "; #InFile
			print COCCAL "OUTPUT=".$outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding."_CalculateHsMetrics "; #OutFile
			print COCCAL "REFERENCE_SEQUENCE=".$scriptParameters{'referencesDir'}."/".$scriptParameters{'humanGenomeReference'}." "; #Reference file
			print COCCAL "BAIT_INTERVALS=".$scriptParameters{'referencesDir'}."/".$sampleInfo{ $scriptParameters{'familyID'} }{$sampleID}{'exomeTargetPaddedBedInfileList'}." "; #Capture kit padded target infile_list file
			print COCCAL "TARGET_INTERVALS=".$scriptParameters{'referencesDir'}."/".$sampleInfo{ $scriptParameters{'familyID'} }{$sampleID}{'exomeTargetBedInfileList'}." &", "\n\n"; #Capture kit target infile_list file
		    }
		    if ($scriptParameters{'pCoverageBED'} > 0) { #Run coverageBed (exome)
			
			print COCCAL "coverageBed ";
			print COCCAL "-hist "; #Report a histogram of coverage for each feature in B as well as a summary histogram for _all_ features in B.
			print COCCAL "-abam ".$inSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding.".bam "; #InFile in BAM format
			print COCCAL "-b ".$scriptParameters{'referencesDir'}."/".$sampleInfo{$scriptParameters{'familyID'}}{$sampleID}{'exomeTargetBed'}." "; #InFile
			print COCCAL "> ".$outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed_hist &", "\n\n"; #OutFile
			#Remove PCR and Optical duplicates
			print COCCAL "samtools view ";
			print COCCAL "-F 0x400 "; #Skip alignments where read is PCR or optical duplicate
			print COCCAL "-b ".$inSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding.".bam "; #InFile
			print COCCAL "| coverageBed "; #Note "|"
			print COCCAL "-hist "; #Report a histogram of coverage for each feature in B as well as a summary histogram for _all_ features in B.
			print COCCAL "-abam stdin "; #InStream
			print COCCAL "-b ".$scriptParameters{'referencesDir'}."/".$sampleInfo{ $scriptParameters{'familyID'} }{$sampleID}{'exomeTargetBed'}." "; #InFile
			print COCCAL "> ".$outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding."_rmdup_coverageBed_hist &", "\n\n"; #OutFile
		    }
		}
	    }
	}
    }
    print COCCAL "wait", "\n\n";

###
#Coverage Report
###    
#Only 1 set of files should be generated with the BAM infile that contains the most information.
    print COCCAL "#Coverage Report Generation\n\n";

    if (@picardToolsMergePreviousFiles ) { # Coverage report on files merged this round with merged file from previous round
	if ( ($scriptParameters{'pPicardToolsMerge'} == 0) && ($scriptParameters{'pPicardToolsMarkduplicates'} == 0) ) { #Ensure correct file ending
	    $fileEnding = "_sorted";
	}
	if ( ($scriptParameters{'pPicardToolsMerge'} > 0) && ($scriptParameters{'pPicardToolsMarkduplicates'} == 0) ) {
	    $fileEnding = "_sorted_merged";
	}
	if ( ($scriptParameters{'pPicardToolsMerge'} == 0) && ($scriptParameters{'pPicardToolsMarkduplicates'} > 0) ) {
	    $fileEnding = "_sorted_pmd";
	}
	if ( ($scriptParameters{'pPicardToolsMerge'} > 0) && ($scriptParameters{'pPicardToolsMarkduplicates'} > 0) ) {
	    $fileEnding = "_sorted_merged_pmd";
	}
	for (my $mergeFileCounter=0;$mergeFileCounter<scalar(@picardToolsMergePreviousFiles);$mergeFileCounter++) {
	    
	    if ($picardToolsMergePreviousFiles[$mergeFileCounter] =~ /$sampleID/) { #Look for sampleID in previously generated file to be merged with current run to be able to merge correct files
		if ($picardToolsMergePreviousFiles[$mergeFileCounter] =~ /lane(\d+)|s_(\d+)/) { #Look for lanes_ or lane\d in previously generated file to be merged with current run to be able to extract previous lanes
		    
		    my $mergeLanes; if($1) {$mergeLanes = $1;} else {$mergeLanes = $2;} #Make sure to always supply lanes from previous regexp
		    
		    print COCCAL "#Returns the depth and breadth of coverage of features from A (-abam)\n";
		    print COCCAL "coverageBed "; #Returns the depth and breadth of coverage of features from A
		    print COCCAL "-abam ".$inSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding.".bam "; #InFile in BAM format
		    print COCCAL "-b ".$scriptParameters{'referencesDir'}."/".$sampleInfo{ $scriptParameters{'familyID'} }{$sampleID}{'exomeTargetBed'}." "; #InFile
		    print COCCAL "> ".$outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed &", "\n\n"; #OutFile

		    print COCCAL "coverageBed ";
		    print COCCAL "-d "; #Report the depth at each position in each B feature.
		    print COCCAL "-abam ".$inSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding.".bam "; #InFile in BAM format
		    print COCCAL "-b ".$scriptParameters{'referencesDir'}."/".$sampleInfo{ $scriptParameters{'familyID'} }{$sampleID}{'exomeTargetBed'}." "; #InFile
		    print COCCAL "> ".$outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed_depth_pos &", "\n\n"; #OutFile
		    
		    print COCCAL "wait", "\n\n";
		    
                    #Generate the average depth per feature and the number of zero-coverage bases as well as their relative to feature position
		    print COCCAL "#Generate the average depth per feature and the number of zero-coverage bases as well as their relative to feature position", "\n";
		    print COCCAL q?perl -nae'my $prev_chr=0;my $prev_start;my $prev_end;my $average_cov=0;my $nr_positions=0;my %feature;my %zerofeature;my @temp_line; while (<>) { chomp($_);  @temp_line = split("\t",$_); if($feature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]}) { if ($temp_line[5] >= 10) { $feature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]}{'high-coverage'}++;} if ($temp_line[5] == 0) { push (@ {$zerofeature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]}{'zero-coverage'}}, $temp_line[1]+$temp_line[4]-1);} $average_cov=$average_cov+$temp_line[5];$nr_positions++;} else { if ( $prev_chr && ($prev_start != $temp_line[1]) && ($prev_end != $temp_line[2]) ) { print "$prev_chr\t","$prev_start\t","$prev_end\t",$average_cov/$nr_positions,"\t"; if ($feature{$prev_chr}{$prev_start}{$prev_end}{'high-coverage'}) { print $feature{$prev_chr}{$prev_start}{$prev_end}{'high-coverage'}/$nr_positions, "\t";} else {print "0\t";} if ($zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}) { print scalar(@ {$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'} }), "\t"; for (my $i=0;$i<scalar(@ {$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}});$i++) { print "$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}[$i];";} print "\n"; } else {print "0\t","Na\n";} } %feature = ();%zerofeature = ();$average_cov=0;$nr_positions=0;$feature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]} = $_; $prev_chr = $temp_line[0];$prev_start = $temp_line[1];$prev_end = $temp_line[2];if ($temp_line[5] >= 10) { $feature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]}{'high-coverage'}++;} if ($temp_line[5] == 0) { push (@ {$zerofeature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]}{'zero-coverage'} }, $temp_line[1]+$temp_line[4]-1);} $average_cov=$average_cov+$temp_line[5];$nr_positions++;} } print "$prev_chr\t","$prev_start\t","$prev_end\t",$average_cov/$nr_positions,"\t"; if ($feature{$prev_chr}{$prev_start}{$prev_end}{'high-coverage'}) { print $feature{$prev_chr}{$prev_start}{$prev_end}{'high-coverage'}/$nr_positions, "\t";} else {print "0\t";} if ($zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}) { print scalar(@ {$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'} }), "\t"; for (my $i=0;$i<scalar(@ {$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}});$i++) { print "$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}[$i];";} print "\n"; } else {print "0\t","Na\n";} last;' ?.$outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed_depth_pos "; #InFile using the just created depth_pos file located in the outSampleDirectory
		    print COCCAL "> ".$outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed_depth_pos_collapsed", "\n\n"; #OutFile
		    
		    #Concatenate to 1 file to be able to include info from _coverageBed file
		    print COCCAL "#Concatenate to 1 file to be able to include info from _coverageBed file", "\n";
		    print COCCAL "cat ";
		    print COCCAL $outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed_depth_pos_collapsed "; #InFile
		    print COCCAL $outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed "; #InFile
		    print COCCAL "> ".$outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed_merged_temp", "\n\n"; #OutFile
		    
		    #Sort on chr and then numerically on start position. NOTE makes the two lines from each file end up next to each other. The line order is preserved unless the feature annotation is "-", then it is reversed (exception handled by the next perl one-liner) 
		    print COCCAL "#Sort on chr and then numerically on start position.\n";
		    print COCCAL "sort ";
		    print COCCAL "-k1,1 -k2,2n "; #Sort on chr and then numerically on start position.
		    print COCCAL $outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed_merged_temp "; #Infile
		    print COCCAL "> ".$outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed_merged", "\n\n"; #OutFile
		    
		    #Collpase the two lines from coverageBed_depth_pos_collapsed and _coverageBed
		    print COCCAL "##Collpase the two lines from coverageBed_depth_pos_collapsed and _coverageBed\n";
		    print COCCAL q?perl -nae ' print "#Chr\tStart\tStop\tOverlapping_Reads\tNon-zero_Coverage_Bases\tFeature_Length(Nucleotides)\tAverage_Coverage\tFraction_Non-zero_Coverage_Bases\tFraction_Ten_Coverage_Bases\tNr_Zero_Coverage_Bases\tRelative_position_for_Zero_Coverage_Bases\n";chomp($_);my @prev_line=split("\t",$_);my @temp_line; while (<>) { chomp($_); @temp_line = split("\t",$_); if ($prev_line[0]  && ($prev_line[0] eq $temp_line[0]) && ($prev_line[1] == $temp_line[1]) && ($prev_line[2] == $temp_line[2]) ) { if ($prev_line[3] eq "-") { $temp_line[3]=~s/\./\,/;$prev_line[7]=~s/\./\,/; print $prev_line[0],"\t", $prev_line[1],"\t",$prev_line[2],"\t",$prev_line[4],"\t",$prev_line[5],"\t",$prev_line[6],"\t",$temp_line[3],"\t","$prev_line[7]\t","$temp_line[4]\t","$temp_line[5]\t","$temp_line[6]\n";} else { $prev_line[3]=~s/\./\,/;$temp_line[7]=~s/\./\,/; print $temp_line[0],"\t", $temp_line[1],"\t",$temp_line[2],"\t",$temp_line[4],"\t",$temp_line[5],"\t",$temp_line[6],"\t",$prev_line[3],"\t","$temp_line[7]\t","$prev_line[4]\t","$prev_line[5]\t","$prev_line[6]\n";} } else { @prev_line = @temp_line;} }last;' ?.$outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed_merged "; #InFile using the just created sorted _coverageBed_merged file located in the outSampleDirectory
		    print COCCAL "> ".$outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding."_target_coverage.txt", "\n\n"; #OutFile
		    
		   
		    #Add chr to entry to enable comparison to Gene Db
		    print COCCAL "#Add chr to entry to enable later comparison to Gene Db", "\n";
		    print COCCAL q?perl -i -p -e ' if($_=~/^#/) {} else {s/^(.+)/chr$1/g }' ?.$outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding."_target_coverage.txt", "\n\n"; #Modify in place
		    
		    #Removal of files which the necessary info has been extracted from
		    print COCCAL "#Removal of files which the necessary info has been extracted from\n";
		    print COCCAL "rm ";
		    print COCCAL $outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed", "\n\n";

		    print COCCAL "rm ";
		    print COCCAL $outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed_depth_pos", "\n\n";
		    print COCCAL "rm ";
		    print $outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed_depth_pos_collapsed", "\n\n";
		    print COCCAL "rm ";
		    print COCCAL $outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed_merged_temp", "\n\n";
		    print COCCAL "rm ";
		    print COCCAL $outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed_merged", "\n\n";
		    
		    my $targetCoverageDbFile = $scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/coverageReport/".$sampleID."_lanes_".$mergeLanes;
		    for (my $laneCounter=0;$laneCounter<scalar(@ { $lanes{$sampleID} });$laneCounter++) {
			$targetCoverageDbFile .= $lanes{$sampleID}[$laneCounter];
		    }
		    my $targetCoverageFile =  $targetCoverageDbFile;

		    $targetCoverageDbFile .= $fileEnding."_coverage_target_db_master.txt";		    
		    $targetCoverageFile .= $fileEnding."_target_coverage.txt";    
		    my @targetCoverageDbFiles = ($targetCoverageDbFile); #Db master files	    
		    my @targetCoverageFiles = ($targetCoverageFile); #Target coverage files created above
		
###
#Create target file with EnsembleGeneID (Required for coverage calculation in rank_script)
###
                    #Create db master template 
		    print COCCAL "#Add GeneName to Coverage Report", "\n";
		    for (my $dbFileCounter=0;$dbFileCounter<scalar(@targetCoverageDbFiles);$dbFileCounter++) {
			open (TARCOV, ">".$targetCoverageDbFiles[$dbFileCounter]) or die "Can't write to ".$targetCoverageDbFiles[$dbFileCounter].": $!\n";
			print TARCOV "outcolumns=0_0,0_1,0_2,0_3,0_4,0_5,0_6,0_7,0_8,0_9,0_10,1_4\n"; #Order of columns in outfile
			print TARCOV "outheaders=Chr,Start,Stop,Overlapping_Reads,Non-zero_Coverage_Bases,Feature_Length(Nucleotides),Average_Coverage,Fraction_Non-zero_Coverage_Bases,Fraction_Ten_Coverage_Bases,Nr_Zero_Coverage_Bases,Relative_position_for_Zero_Coverage_Bases,Ensembl_GeneID\n"; #Order and header content in outfile
			print TARCOV $targetCoverageFiles[$dbFileCounter],"\t".'\t'."\t0,1,2\t0\texact\t0,1,2,3,4,5,6,7,8,9,10\tsmall","\n";
			print TARCOV $scriptParameters{'referencesDir'}."/mart_export_Ensembl_GeneID_key_cleaned_noblanks_chr.txt\t".'\t'."\t0,1,2\t0\trange\t4\tsmall", "\n";
			close(TARCOV);
			
                        #Add GeneNameID to Coverage report
			print COCCAL "perl ";
			print COCCAL $scriptParameters{'inScriptDir'}."/intersectCollect.pl ";
			print COCCAL "-o ".$targetCoverageFiles[$dbFileCounter]." "; #OutFile
			print "-db ".$targetCoverageDbFiles[$dbFileCounter]." "; #Db master file (InFile(s))
			print "-prechr 1", "\n\n"; #Use chromosome prefix
		    }   
		}
	    }
	}
    }
    
    elsif ( ($scriptParameters{'pPicardToolsMerge'} > 0) && (scalar( @{ $InfilesLaneNoEnding{$sampleID} }) > 1) ) { #but only if there is more than one mosaikBuild/BWA_Aln file per sample ID (Sanity check)
	if ($scriptParameters{'pPicardToolsMarkduplicates'} == 0) { #Ensure correct file ending
	    $fileEnding = "_sorted_merged";
	}
	else {
	    $fileEnding = "_sorted_merged_pmd";
	}

	print COCCAL "coverageBed "; #Returns the depth and breadth of coverage of features from A
	print COCCAL "-abam ".$inSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding.".bam "; #InFile in BAM format
	print COCCAL "-b ".$scriptParameters{'referencesDir'}."/".$sampleInfo{ $scriptParameters{'familyID'} }{$sampleID}{'exomeTargetBed'}." "; #InFile
	print COCCAL "> ".$outSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed &", "\n\n"; #OutFile

	print COCCAL "coverageBed ";
	print COCCAL "-d "; #Report the depth at each position in each B feature.
	print COCCAL "-abam ".$inSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding.".bam "; #InFile in BAM format
	print COCCAL "-b ".$scriptParameters{'referencesDir'}."/".$sampleInfo{ $scriptParameters{'familyID'} }{$sampleID}{'exomeTargetBed'}." "; #InFile
	print COCCAL "> ".$outSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed_depth_pos &", "\n\n"; #OutFile

	print COCCAL "wait", "\n\n";

	#Generate the average depth per feature and the number of zero-coverage bases as well as their relative to feature position
	print COCCAL "#Generate the average depth per feature and the number of zero-coverage bases as well as their relative to feature position", "\n";
	print COCCAL q?perl -nae'my $prev_chr=0;my $prev_start;my $prev_end;my $average_cov=0;my $nr_positions=0;my %feature;my %zerofeature;my @temp_line; while (<>) { chomp($_);  @temp_line = split("\t",$_); if($feature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]}) { if ($temp_line[5] >= 10) { $feature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]}{'high-coverage'}++;} if ($temp_line[5] == 0) { push (@ {$zerofeature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]}{'zero-coverage'}}, $temp_line[1]+$temp_line[4]-1);} $average_cov=$average_cov+$temp_line[5];$nr_positions++;} else { if ( $prev_chr && ($prev_start != $temp_line[1]) && ($prev_end != $temp_line[2]) ) { print "$prev_chr\t","$prev_start\t","$prev_end\t",$average_cov/$nr_positions,"\t"; if ($feature{$prev_chr}{$prev_start}{$prev_end}{'high-coverage'}) { print $feature{$prev_chr}{$prev_start}{$prev_end}{'high-coverage'}/$nr_positions, "\t";} else {print "0\t";} if ($zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}) { print scalar(@ {$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'} }), "\t"; for (my $i=0;$i<scalar(@ {$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}});$i++) { print "$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}[$i];";} print "\n"; } else {print "0\t","Na\n";} } %feature = ();%zerofeature = ();$average_cov=0;$nr_positions=0;$feature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]} = $_; $prev_chr = $temp_line[0];$prev_start = $temp_line[1];$prev_end = $temp_line[2];if ($temp_line[5] >= 10) { $feature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]}{'high-coverage'}++;} if ($temp_line[5] == 0) { push (@ {$zerofeature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]}{'zero-coverage'} }, $temp_line[1]+$temp_line[4]-1);} $average_cov=$average_cov+$temp_line[5];$nr_positions++;} } print "$prev_chr\t","$prev_start\t","$prev_end\t",$average_cov/$nr_positions,"\t"; if ($feature{$prev_chr}{$prev_start}{$prev_end}{'high-coverage'}) { print $feature{$prev_chr}{$prev_start}{$prev_end}{'high-coverage'}/$nr_positions, "\t";} else {print "0\t";} if ($zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}) { print scalar(@ {$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'} }), "\t"; for (my $i=0;$i<scalar(@ {$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}});$i++) { print "$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}[$i];";} print "\n"; } else {print "0\t","Na\n";} last;' ?.$outSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed_depth_pos "; #InFile using the just created depth_pos file located in the outSampleDirectory
	print COCCAL "> ".$outSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed_depth_pos_collapsed", "\n\n"; #OutFile
	
	#Concatenate to 1 file to be able to include info from _coverageBed file
	print COCCAL "#Concatenate to 1 file to be able to include info from _coverageBed file", "\n";
	print COCCAL "cat ";
	print COCCAL $outSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed_depth_pos_collapsed "; #InFile
	print COCCAL $outSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed "; #InFile
	print COCCAL "> ".$outSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed_merged_temp", "\n\n"; #OutFile
	
	#Sort on chr and then numerically on start position. NOTE makes the two lines from each file end up next to each other. The line order is preserved unless the feature annotation is "-", then it is reversed (exception handled by the next perl one-liner) 
	print COCCAL "#Sort on chr and then numerically on start position.\n";
	print COCCAL "sort ";
	print COCCAL "-k1,1 -k2,2n "; #Sort on chr and then numerically on start position.
	print COCCAL $outSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed_merged_temp "; #Infile
	print COCCAL "> ".$outSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed_merged", "\n\n"; #OutFile
	
	#Collpase the two lines from coverageBed_depth_pos_collapsed and _coverageBed
	print COCCAL "##Collpase the two lines from coverageBed_depth_pos_collapsed and _coverageBed\n";
	print COCCAL q?perl -nae ' print "#Chr\tStart\tStop\tOverlapping_Reads\tNon-zero_Coverage_Bases\tFeature_Length(Nucleotides)\tAverage_Coverage\tFraction_Non-zero_Coverage_Bases\tFraction_Ten_Coverage_Bases\tNr_Zero_Coverage_Bases\tRelative_position_for_Zero_Coverage_Bases\n";chomp($_);my @prev_line=split("\t",$_);my @temp_line; while (<>) { chomp($_); @temp_line = split("\t",$_); if ($prev_line[0]  && ($prev_line[0] eq $temp_line[0]) && ($prev_line[1] == $temp_line[1]) && ($prev_line[2] == $temp_line[2]) ) { if ($prev_line[3] eq "-") { $temp_line[3]=~s/\./\,/;$prev_line[7]=~s/\./\,/; print $prev_line[0],"\t", $prev_line[1],"\t",$prev_line[2],"\t",$prev_line[4],"\t",$prev_line[5],"\t",$prev_line[6],"\t",$temp_line[3],"\t","$prev_line[7]\t","$temp_line[4]\t","$temp_line[5]\t","$temp_line[6]\n";} else { $prev_line[3]=~s/\./\,/;$temp_line[7]=~s/\./\,/; print $temp_line[0],"\t", $temp_line[1],"\t",$temp_line[2],"\t",$temp_line[4],"\t",$temp_line[5],"\t",$temp_line[6],"\t",$prev_line[3],"\t","$temp_line[7]\t","$prev_line[4]\t","$prev_line[5]\t","$prev_line[6]\n";} } else { @prev_line = @temp_line;} }last;' ?.$outSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed_merged "; #InFile using the just created sorted _coverageBed_merged file located in the outSampleDirectory
	print COCCAL "> ".$outSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding."_target_coverage.txt", "\n\n"; #OutFile
	
	#Add chr to entry to enable comparison to Gene Db
	print COCCAL "#Add chr to entry to enable later comparison to Gene Db", "\n";
	print COCCAL q?perl -i -p -e ' if($_=~/^#/) {} else {s/^(.+)/chr$1/g }' ?.$outSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding."_target_coverage.txt", "\n\n"; #Modify in place	   

	#Removal of files which the necessary info has been extracted from
	print COCCAL "#Removal of files which the necessary info has been extracted from\n";
	print COCCAL "rm ";
	print COCCAL $outSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed", "\n\n";

	print COCCAL "rm ";
	print COCCAL $outSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed_depth_pos", "\n\n";
	
	print COCCAL "rm ";
	print COCCAL $outSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed_depth_pos_collapsed", "\n\n";
	
	print COCCAL "rm ";
	print COCCAL $outSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed_merged_temp", "\n\n";
	
	print COCCAL "rm ";
	print COCCAL $outSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,$fileEnding."_coverageBed_merged", "\n\n";

	my $targetCoverageDbFile = $scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/coverageReport/".$sampleID."_lanes_";
	for (my $laneCounter=0;$laneCounter<scalar(@ { $lanes{$sampleID} });$laneCounter++) {
	    $targetCoverageDbFile .= $lanes{$sampleID}[$laneCounter];
	}
	my $targetCoverageFile =  $targetCoverageDbFile;

	$targetCoverageDbFile .= $fileEnding."_coverage_target_db_master.txt";
	$targetCoverageFile .= $fileEnding."_target_coverage.txt";
	my @targetCoverageDbFiles = ($targetCoverageDbFile); #Db master files	    
	my @targetCoverageFiles = ($targetCoverageFile); #Target coverage files created above
###
#Create target file with EnsembleGeneID (Required for coverage calculation in rank_script)
###
	
#Create db master template 
	print COCCAL "#Add GeneName to Coverage Report", "\n";
	for (my $dbFileCounter=0;$dbFileCounter<scalar(@targetCoverageDbFiles);$dbFileCounter++) {
	    
	    open (TARCOV, ">".$targetCoverageDbFiles[$dbFileCounter]) or die "Can't write to ".$targetCoverageDbFiles[$dbFileCounter].": $!\n";
	    print TARCOV "outcolumns=0_0,0_1,0_2,0_3,0_4,0_5,0_6,0_7,0_8,0_9,0_10,1_4\n"; #Order of outcolumns in outfile
	    print TARCOV "outheaders=Chr,Start,Stop,Overlapping_Reads,Non-zero_Coverage_Bases,Feature_Length(Nucleotides),Average_Coverage,Fraction_Non-zero_Coverage_Bases,Fraction_Ten_Coverage_Bases,Nr_Zero_Coverage_Bases,Relative_position_for_Zero_Coverage_Bases,Ensembl_GeneID\n"; #Order and header content in outfile
	    print TARCOV $targetCoverageFiles[$dbFileCounter],"\t".'\t'."\t0,1,2\t0\texact\t0,1,2,3,4,5,6,7,8,9,10\tsmall","\n";
	    print TARCOV $scriptParameters{'referencesDir'}."/mart_export_Ensembl_GeneID_key_cleaned_noblanks_chr.txt\t".'\t'."\t0,1,2\t0\trange\t4\tsmall", "\n";
	    close(TARCOV);
	    
#Add GeneNameID
	    print COCCAL "perl ";
	    print COCCAL $scriptParameters{'inScriptDir'}."/intersectCollect.pl ";
	    print COCCAL "-o ".$targetCoverageFiles[$dbFileCounter]." "; #OutFile
	    print COCCAL "-db ".$targetCoverageDbFiles[$dbFileCounter]. " "; #Db master file (InFile(s))
	    print COCCAL "-prechr 1", "\n\n"; #Chromosome prefix
	}
    }
    else {
	if ($scriptParameters{'pPicardToolsMarkduplicates'} == 0) { #Ensure correct file ending
	    $fileEnding = "_sorted";
	}
	else {
	    $fileEnding = "_sorted_pmd";
	}
	for (my $infileCounter=0;$infileCounter<scalar( @{ $InfilesLaneNoEnding{$sampleID} } );$infileCounter++) { #For all files from MosaikAlign or BWA_Sampe
	    
	    my $tempInfile = $InfilesLaneNoEnding{$sampleID}[$infileCounter];

	    print COCCAL "#Calculate coverage statistics to enable coverage calculation in rank_script", "\n";
	    print COCCAL "coverageBed "; #Returns the depth and breadth of coverage of features from A
	    print COCCAL "-abam ".$inSampleDirectory."/".$tempInfile.$fileEnding.".bam "; #InFiles in BAM format
	    print COCCAL "-b ".$scriptParameters{'referencesDir'}."/".$sampleInfo{ $scriptParameters{'familyID'} }{$sampleID}{'exomeTargetBed'}." "; #Infile
	    print COCCAL "> ".$outSampleDirectory."/".$tempInfile.$fileEnding."_coverageBed &", "\n\n"; #OutFile

	    print COCCAL "coverageBed ";
	    print COCCAL "-d "; #Report the depth at each position in each B feature.
	    print COCCAL "-abam ".$inSampleDirectory."/".$tempInfile.$fileEnding.".bam "; #InFiles in BAM format
	    print COCCAL "-b ".$scriptParameters{'referencesDir'}."/".$sampleInfo{ $scriptParameters{'familyID'} }{$sampleID}{'exomeTargetBed'}." "; #InFile
	    print COCCAL "> ".$outSampleDirectory."/".$tempInfile.$fileEnding."_coverageBed_depth_pos &", "\n\n"; #OutFile

	    print COCCAL "wait", "\n\n";

#Generate the average depth per feature and the number of zero-coverage bases as well as their relative to feature position
	    print COCCAL "#Generate the average depth per feature and the number of zero-coverage bases as well as their relative to feature position", "\n";
	    print COCCAL q?perl -nae'my $prev_chr=0;my $prev_start;my $prev_end;my $average_cov=0;my $nr_positions=0;my %feature;my %zerofeature;my @temp_line; while (<>) { chomp($_);  @temp_line = split("\t",$_); if($feature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]}) { if ($temp_line[5] >= 10) { $feature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]}{'high-coverage'}++;} if ($temp_line[5] == 0) { push (@ {$zerofeature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]}{'zero-coverage'}}, $temp_line[1]+$temp_line[4]-1);} $average_cov=$average_cov+$temp_line[5];$nr_positions++;} else { if ( $prev_chr && ($prev_start != $temp_line[1]) && ($prev_end != $temp_line[2]) ) { print "$prev_chr\t","$prev_start\t","$prev_end\t",$average_cov/$nr_positions,"\t"; if ($feature{$prev_chr}{$prev_start}{$prev_end}{'high-coverage'}) { print $feature{$prev_chr}{$prev_start}{$prev_end}{'high-coverage'}/$nr_positions, "\t";} else {print "0\t";} if ($zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}) { print scalar(@ {$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'} }), "\t"; for (my $i=0;$i<scalar(@ {$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}});$i++) { print "$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}[$i];";} print "\n"; } else {print "0\t","Na\n";} } %feature = ();%zerofeature = ();$average_cov=0;$nr_positions=0;$feature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]} = $_; $prev_chr = $temp_line[0];$prev_start = $temp_line[1];$prev_end = $temp_line[2];if ($temp_line[5] >= 10) { $feature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]}{'high-coverage'}++;} if ($temp_line[5] == 0) { push (@ {$zerofeature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]}{'zero-coverage'} }, $temp_line[1]+$temp_line[4]-1);} $average_cov=$average_cov+$temp_line[5];$nr_positions++;} } print "$prev_chr\t","$prev_start\t","$prev_end\t",$average_cov/$nr_positions,"\t"; if ($feature{$prev_chr}{$prev_start}{$prev_end}{'high-coverage'}) { print $feature{$prev_chr}{$prev_start}{$prev_end}{'high-coverage'}/$nr_positions, "\t";} else {print "0\t";} if ($zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}) { print scalar(@ {$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'} }), "\t"; for (my $i=0;$i<scalar(@ {$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}});$i++) { print "$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}[$i];";} print "\n"; } else {print "0\t","Na\n";} last;' ?.$outSampleDirectory."/".$tempInfile.$fileEnding."_coverageBed_depth_pos "; #InFile using the just created depth_pos file located in the outSampleDirectory
	    print COCCAL "> ".$outSampleDirectory."/".$tempInfile.$fileEnding."_coverageBed_depth_pos_collapsed", "\n\n"; #OutFile
	    
	    #Concatenate to 1 file to be able to include info from _coverageBed file
	    print COCCAL "#Concatenate to 1 file to be able to include info from _coverageBed file", "\n";
	    print COCCAL "cat ";
	    print COCCAL $outSampleDirectory."/".$tempInfile.$fileEnding."_coverageBed_depth_pos_collapsed "; #InFile
	    print COCCAL $outSampleDirectory."/".$tempInfile.$fileEnding."_coverageBed "; #InFile
	    print COCCAL "> ".$outSampleDirectory."/".$tempInfile.$fileEnding."_coverageBed_merged_temp",  "\n\n"; #OutFile
	    
            #Sort on chr and then numerically on start position. NOTE makes the two lines from each file end up next to each other. The line order is preserved unless the feature annotation is "-", then it is reversed (exception handled by the next perl one-liner) 
	    print COCCAL "#Sort on chr and then numerically on start position.\n";
	    print COCCAL "sort ";
	    print COCCAL "-k1,1 -k2,2n "; #Sort on chr and then numerically on start position.
	    print COCCAL $outSampleDirectory."/".$tempInfile.$fileEnding."_coverageBed_merged_temp "; #InFile
	    print COCCAL "> ".$outSampleDirectory."/".$tempInfile.$fileEnding."_coverageBed_merged",  "\n\n"; #OutFile
	    
	    #Collpase the two lines from coverageBed_depth_pos_collapsed and _coverageBed
	    print COCCAL "#Collpase the two lines from coverageBed_depth_pos_collapsed and _coverageBed", "\n";
	    print COCCAL q?perl -nae ' print "#Chr\tStart\tStop\tOverlapping_Reads\tNon-zero_Coverage_Bases\tFeature_Length(Nucleotides)\tAverage_Coverage\tFraction_Non-zero_Coverage_Bases\tFraction_Ten_Coverage_Bases\tNr_Zero_Coverage_Bases\tRelative_position_for_Zero_Coverage_Bases\n";chomp($_);my @prev_line=split("\t",$_);my @temp_line; while (<>) { chomp($_); @temp_line = split("\t",$_); if ($prev_line[0]  && ($prev_line[0] eq $temp_line[0]) && ($prev_line[1] == $temp_line[1]) && ($prev_line[2] == $temp_line[2]) ) { if ($prev_line[3] eq "-") { $temp_line[3]=~s/\./\,/;$prev_line[7]=~s/\./\,/; print $prev_line[0],"\t", $prev_line[1],"\t",$prev_line[2],"\t",$prev_line[4],"\t",$prev_line[5],"\t",$prev_line[6],"\t",$temp_line[3],"\t","$prev_line[7]\t","$temp_line[4]\t","$temp_line[5]\t","$temp_line[6]\n";} else { $prev_line[3]=~s/\./\,/;$temp_line[7]=~s/\./\,/; print $temp_line[0],"\t", $temp_line[1],"\t",$temp_line[2],"\t",$temp_line[4],"\t",$temp_line[5],"\t",$temp_line[6],"\t",$prev_line[3],"\t","$temp_line[7]\t","$prev_line[4]\t","$prev_line[5]\t","$prev_line[6]\n";} } else { @prev_line = @temp_line;} }last;' ?.$outSampleDirectory."/".$tempInfile.$fileEnding."_coverageBed_merged "; #InFile using the just created sorted _coverageBed_merged file located in the outSampleDirectory
	    print COCCAL "> ".$outSampleDirectory."/".$tempInfile.$fileEnding."_target_coverage.txt", "\n\n"; #OutFile
	    
	    #Add chr to entry to enable comparison to Gene Db
	    print COCCAL "#Add chr to entry to enable later comparison to Gene Db", "\n";
	    print COCCAL q?perl -i -p -e ' if($_=~/^#/) {} else {s/^(.+)/chr$1/g }' ?.$outSampleDirectory."/".$tempInfile.$fileEnding."_target_coverage.txt", "\n\n";
            #Removal of files which the necessary info has been extracted from
	    print COCCAL "#Removal of files which the necessary info has been extracted from\n";
	    print COCCAL "rm ";
	    print COCCAL $outSampleDirectory."/".$tempInfile.$fileEnding."_coverageBed", "\n\n";

	    print COCCAL "rm ";
	    print COCCAL $outSampleDirectory."/".$tempInfile.$fileEnding."_coverageBed_depth_pos", "\n\n";

	    print COCCAL "rm ";
	    print COCCAL $outSampleDirectory."/".$tempInfile.$fileEnding."_coverageBed_depth_pos_collapsed", "\n\n";
	    
	    print COCCAL "rm ";
	    print COCCAL $outSampleDirectory."/".$tempInfile.$fileEnding."_coverageBed_merged_temp", "\n\n";
	    
	    print COCCAL "rm ";
	    print COCCAL $outSampleDirectory."/".$tempInfile.$fileEnding."_coverageBed_merged", "\n\n";
	    
###
#Create target file with EnsembleGeneID (Required for coverage calculation in rank_script)
###
	    
	    my @targetCoverageDbFiles = ($scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/coverageReport/".$tempInfile.$fileEnding."_coverage_target_db_master.txt"); #Db master files	    
	    my @targetCoverageFiles = ($scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/coverageReport/".$tempInfile.$fileEnding."_target_coverage.txt"); 

#Create db master template 
	    print COCCAL "#Add GeneName to Coverage Report", "\n";
	    for (my $dbFileCounter=0;$dbFileCounter<scalar(@targetCoverageDbFiles);$dbFileCounter++) { 
		open (TARCOV, ">".$targetCoverageDbFiles[$dbFileCounter]) or die "Can't write to ".$targetCoverageDbFiles[$dbFileCounter].": $!\n";
		print TARCOV "outcolumns=0_0,0_1,0_2,0_3,0_4,0_5,0_6,0_7,0_8,0_9,0_10,1_4\n"; #Order of outcolumns in outfile
		print TARCOV "outheaders=Chr,Start,Stop,Overlapping_Reads,Non-zero_Coverage_Bases,Feature_Length(Nucleotides),Average_Coverage,Fraction_Non-zero_Coverage_Bases,Fraction_Ten_Coverage_Bases,Nr_Zero_Coverage_Bases,Relative_position_for_Zero_Coverage_Bases,Ensembl_GeneID\n"; #Order and header content in outfile
		print TARCOV $targetCoverageFiles[$dbFileCounter],"\t".'\t'."\t0,1,2\t0\texact\t0,1,2,3,4,5,6,7,8,9,10\tsmall","\n";
		print TARCOV $scriptParameters{'referencesDir'}."/mart_export_Ensembl_GeneID_key_cleaned_noblanks_chr.txt\t".'\t'."\t0,1,2\t0\trange\t4\tsmall", "\n";
		close(TARCOV);

#Add GeneNameID to Coverage Report
		print COCCAL "perl ";
		print COCCAL $scriptParameters{'inScriptDir'}."/intersectCollect.pl ";
		print COCCAL "-o ".$targetCoverageFiles[$dbFileCounter]." "; #OutFile
		print COCCAL "-db ".$targetCoverageDbFiles[$dbFileCounter]." "; #Db master file (InFile(s)) 
		print COCCAL "-prechr 1", "\n\n"; #Chromosome prefix
	    }
	}
    }
    close(COCCAL);
    if ( $scriptParameters{'pCalculateCoverage'} == 1) {
	ParallelSampleIDSubmitJob($sampleID,$filename,"all");
    }
    return;
}

sub PicardToolsMarkDuplicates { 
#Mark duplicated reads using PicardTools MarkDuplicates in files generated from alignment (sorted, merged)

    my $sampleID = $_[0];
    my $aligner = $_[1]; 
    
    `mkdir -p $scriptParameters{'outDataDir'}/$sampleID/$aligner/info;`; #Creates the aligner folder and info data file directory
    `mkdir -p $scriptParameters{'outScriptDir'}/$sampleID/$aligner;`; #Creates the aligner script directory
    if ( $scriptParameters{'pPicardToolsMarkduplicates'} == 1 ) {
	$filename = $scriptParameters{'outScriptDir'}."/".$sampleID."/".$aligner."/picardTools_markdup_".$sampleID.".";
    }
    elsif ( $scriptParameters{'pPicardToolsMarkduplicates'} == 2 ) {
	$filename = $scriptParameters{'outScriptDir'}."/".$sampleID."/".$aligner."/dry_run_picardTools_markdup_".$sampleID.".";
	print STDOUT "Dry Run:\n";print MASTERL  "Dry Run:\n";
    }
    Checkfnexists($filename, $fnend);

#Info and Logg
    print STDOUT "Creating sbatch script PicardToolsMarkDuplicates and writing script file(s) to: ".$filename, "\n";print MASTERL "Creating sbatch script PicardToolsMarkDuplicates and writing script file(s) to: ".$filename, "\n";
    print STDOUT "Sbatch script PicardToolsMarkDuplicates data files will be written to: ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner, "\n";print MASTERL "Sbatch script PicardToolsMarkDuplicates data files will be written to: ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner, "\n";

    my $time = ceil(3*scalar( @{ $InfilesBothStrandsNoEnding{$sampleID} })); #One full lane on Hiseq takes approx. 3 h to process, round up to nearest full hour.
    
    open (PMDUP, ">".$filename) or die "Can't write to ".$filename.": $!\n";
    
    print PMDUP "#! /bin/bash -l", "\n";
    print PMDUP "#SBATCH -A ".$scriptParameters{'projectID'}, "\n";
    print PMDUP "#SBATCH -p node -n ".$scriptParameters{'maximumCores'}, "\n";
    print PMDUP "#SBATCH -C thin", "\n";
    
    if ($scriptParameters{'pPicardToolsMerge'} eq 0) { #If No merge has been performed then time requirements goes down
	print PMDUP "#SBATCH -t 3:00:00", "\n";	
    }
    else{
	print PMDUP "#SBATCH -t ".$time.":00:00", "\n";	
    }	
    
    print PMDUP "#SBATCH -J PMDUP_".$sampleID."_".$aligner, "\n";
    if ( $scriptParameters{'pPicardToolsMarkduplicates'} == 1 ) {
	print PMDUP "#SBATCH -e ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/info/picardTools_markdup_".$sampleID.".".$fnt.".stderr.txt", "\n";
	print PMDUP "#SBATCH -o ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/info/picardTools_markdup_".$sampleID.".".$fnt.".stdout.txt", "\n";
    }
    elsif ( $scriptParameters{'pPicardToolsMarkduplicates'} == 2 ) { #Dry run
	print PMDUP "#SBATCH -e ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/info/dry_run_picardTools_markdup_".$sampleID.".".$fnt.".stderr.txt", "\n";
	print PMDUP "#SBATCH -o ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/info/dry_run_picardTools_markdup_".$sampleID.".".$fnt.".stdout.txt", "\n";
    }

    unless ($scriptParameters{'email'} eq 0) {
	print PMDUP "#SBATCH --mail-type=END", "\n";
	print PMDUP "#SBATCH --mail-type=FAIL", "\n";
	print PMDUP "#SBATCH --mail-user=".$scriptParameters{'email'}, "\n\n";
    }
    
    print PMDUP 'echo "Running on: $(hostname)"',"\n\n";
    
    my $inSampleDirectory = $scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner;
    my $outSampleDirectory = $scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner;
    my $mergeLanes; #To pick up merged lanes later if required PicardTools Markduplicates and samTools
    my $coreCounter=1;
###
#PicardToolsMarkDuplicates
###
    for (my $infileCounter=0;$infileCounter<scalar( @{ $InfilesLaneNoEnding{$sampleID} } );$infileCounter++) { #For all files
	
	if ($infileCounter eq $coreCounter*$scriptParameters{'maximumCores'}) { #Using only $scriptParameters{'maximumCores'} cores
	    
	    print PMDUP "wait", "\n\n";
	    $coreCounter=$coreCounter+1;
	}
	my $tempInfile = $InfilesLaneNoEnding{$sampleID}[$infileCounter];
	
	print PMDUP "java -Xmx4g -jar ".$scriptParameters{'picardToolsPath'}."/MarkDuplicates.jar ";
	print PMDUP "ASSUME_SORTED=true ";
	print PMDUP "REMOVE_DUPLICATES=false ";
	print PMDUP "VALIDATION_STRINGENCY=LENIENT ";
	print PMDUP "INPUT=".$inSampleDirectory."/".$tempInfile."_sorted.bam "; #InFile
	print PMDUP "OUTPUT=".$outSampleDirectory."/".$tempInfile."_sorted_pmd.bam "; #OutFile
	print PMDUP "METRICS_FILE=".$outSampleDirectory."/".$tempInfile."_sorted_pmdmetric &","\n\n"; #Metric file  
    }
    print PMDUP "wait", "\n\n";
    
    if ( ($scriptParameters{'pPicardToolsMerge'} > 0) && (scalar( @{ $InfilesLaneNoEnding{$sampleID} }) > 1) ) { #but only if there is more than one mosaikBuild/BWA_Aln file per sample ID (Sanity check)
	
	print PMDUP "java -Xmx4g -jar ".$scriptParameters{'picardToolsPath'}."/MarkDuplicates.jar ";
	print PMDUP "ASSUME_SORTED=true ";
	print PMDUP "REMOVE_DUPLICATES=false ";
	print PMDUP "VALIDATION_STRINGENCY=LENIENT ";
	print PMDUP "INPUT=".$inSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,"_sorted_merged.bam "; #InFile
	print PMDUP "OUTPUT=".$outSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,"_sorted_merged_pmd.bam "; #OutFile
	print PMDUP "METRICS_FILE=".$outSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,"_sorted_merged_pmdmetric &", "\n\n"; #Metric file 
    }
    
    if ( @picardToolsMergePreviousFiles ) { # Coverage report on files merged this round with merged file from previous round
	
	for (my $mergeFileCounter=0;$mergeFileCounter<scalar(@picardToolsMergePreviousFiles);$mergeFileCounter++) {
	    
	    if ($picardToolsMergePreviousFiles[$mergeFileCounter] =~ /$sampleID/) { #Look for sampleID in previously generated file to be merged with current run to be able to merge correct files
		if ($picardToolsMergePreviousFiles[$mergeFileCounter] =~ /lane(\d+)|s_(\d+)/) { #Look for lanes_ or lane\d in previously generated file to be merged with current run to be able to extract previous lanes
		  
		    if($1) {$mergeLanes = $1;} else {$mergeLanes = $2;} #Make sure to always supply lanes from previous regexp  
		    print PMDUP "java -Xmx4g -jar ".$scriptParameters{'picardToolsPath'}."/MarkDuplicates.jar ";
		    print PMDUP "ASSUME_SORTED=true ";
		    print PMDUP "REMOVE_DUPLICATES=false ";
		    print PMDUP "VALIDATION_STRINGENCY=LENIENT ";
		    print PMDUP "INPUT=".$inSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,"_sorted_merged.bam "; #InFile
		    print PMDUP "OUTPUT=".$outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,"_sorted_merged_pmd.bam "; #OutFile
		    print PMDUP "METRICS_FILE=".$outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,"_sorted_merged_pmdmetric &", "\n\n"; #Metric file
		}
	    }
	}
	print PMDUP "wait", "\n\n";
    }
###
#SamTools index on just created _sorted(_merged)_pmd.bam
###
    for (my $infileCounter=0;$infileCounter<scalar( @{ $InfilesLaneNoEnding{$sampleID} } );$infileCounter++) { #For all files from alignment
	
	if ($infileCounter eq $coreCounter*$scriptParameters{'maximumCores'}) { #Using only $scriptParameters{'maximumCores'} cores
	    
	    print PMDUP "wait", "\n\n";
	    $coreCounter=$coreCounter+1;
	}

	my $tempInfile = $InfilesLaneNoEnding{$sampleID}[$infileCounter];

	print PMDUP "samtools index ";
	print PMDUP $outSampleDirectory."/".$tempInfile."_sorted_pmd.bam &","\n\n"; #Just created dedupped inFile located in outSamplesDirectory 
    }
    print PMDUP "wait", "\n\n";
    if ( ($scriptParameters{'pPicardToolsMerge'} > 0) && (scalar( @{ $InfilesLaneNoEnding{$sampleID} }) > 1) ) { #SamTools index on merged file this analysis round
	print PMDUP "samtools index ";
	print PMDUP $outSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,"_sorted_merged_pmd.bam &","\n\n"; #Just created merged and dedupped inFile located in outSamplesDirectory 
    }
    if (@picardToolsMergePreviousFiles ) { # SamTools index on files merged this round with merged file from previous round
	print PMDUP "samtools index ";
	print PMDUP $outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,"_sorted_merged_pmd.bam &","\n\n"; #$mergeLanes should not have changed since if (@picardToolsMergePreviousFiles ) { so no need for regexp
	print PMDUP "wait", "\n\n";    
    }
    print PMDUP "wait", "\n\n";
    close(PMDUP);
    if ( $scriptParameters{'pPicardToolsMarkduplicates'} == 1 ) {
	ParallelSampleIDSubmitJob($sampleID,$filename,"all");
    }
    return;
}

sub PicardToolsMerge { 
#Merges all bam files using PicardTools MergeSamFiles within each sampleid and files generated previously (option if provided with '-picardToolsMergePreviousFiles'). The merged files have to be sorted before attempting to merge.
 
    my $sampleID = $_[0];
    my $aligner = $_[1];   
    
    `mkdir -p $scriptParameters{'outDataDir'}/$sampleID/$aligner/info;`; #Creates the aligner folder and info data file directory
    `mkdir -p $scriptParameters{'outScriptDir'}/$sampleID/$aligner;`; #Creates the aligner folder and info data file directory
    if ( $scriptParameters{'pPicardToolsMerge'} == 1 ) {
	$filename = $scriptParameters{'outScriptDir'}."/".$sampleID."/".$aligner."/picardTools_merge_".$sampleID.".";
    }
    elsif ( $scriptParameters{'pPicardToolsMerge'} == 2) { #Dry run
	$filename = $scriptParameters{'outScriptDir'}."/".$sampleID."/".$aligner."/dry_run_picardTools_merge_".$sampleID.".";
	print STDOUT "Dry Run:\n";print MASTERL  "Dry Run:\n";
    }
    Checkfnexists($filename, $fnend);

#Info and Logg
    print STDOUT "Creating sbatch script PicardToolsMerge and writing script file(s) to: ".$filename, "\n";print MASTERL "Creating sbatch script PicardToolsMerge and writing script file(s) to: ".$filename, "\n";
    print STDOUT "Sbatch script PicardToolsMerge data files will be written to: ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner, "\n";print MASTERL "Sbatch script PicardToolsMerge data files will be written to: ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner, "\n";
    
    open (PMERGE, ">".$filename) or die "Can't write to ".$filename.": $!\n";
    
    print PMERGE "#! /bin/bash -l", "\n";
    print PMERGE "#SBATCH -A ".$scriptParameters{'projectID'}, "\n";
    print PMERGE "#SBATCH -p node -n ".$scriptParameters{'maximumCores'}, "\n";
    print PMERGE "#SBATCH -C thin", "\n";	
    print PMERGE "#SBATCH -t 20:00:00", "\n";
    
    print PMERGE "#SBATCH -J PMerge_".$sampleID."_".$aligner, "\n";
    if ( $scriptParameters{'pPicardToolsMerge'} == 1 ) {
	print PMERGE "#SBATCH -e ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/info/picardTools_merge_".$sampleID.".".$fnt.".stderr.txt", "\n";
	print PMERGE "#SBATCH -o ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/info/picardTools_merge_".$sampleID.".".$fnt.".stdout.txt", "\n";
    }
    elsif ( $scriptParameters{'pPicardToolsMerge'} == 2) { #Dry run
	print PMERGE "#SBATCH -e ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/info/dry_run_picardTools_merge_".$sampleID.".".$fnt.".stderr.txt", "\n";
	print PMERGE "#SBATCH -o ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/info/dry_run_picardTools_merge_".$sampleID.".".$fnt.".stdout.txt", "\n";
    }
    unless ($scriptParameters{'email'} eq 0) {
	print PMERGE "#SBATCH --mail-type=END", "\n";
	print PMERGE "#SBATCH --mail-type=FAIL", "\n";
	print PMERGE "#SBATCH --mail-user=".$scriptParameters{'email'}, "\n\n";
    }
    
    print PMERGE 'echo "Running on: $(hostname)"',"\n\n";
   
    my $inSampleDirectory = $scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner;
    my $outSampleDirectory = $scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner;

    if (scalar( @{ $InfilesLaneNoEnding{$sampleID} } ) > 1) { #Check that we have something to merge and then merge current files before merging with previously merged files
	for (my $infileCounter=0;$infileCounter<scalar( @{ $InfilesLaneNoEnding{$sampleID} } );$infileCounter++) { #For all files from 

	    my $tempInfile = $InfilesLaneNoEnding{$sampleID}[$infileCounter];
	    
	    if ($infileCounter eq 0) {

		print PMERGE "java -Xmx4g -jar ".$scriptParameters{'picardToolsPath'}."/MergeSamFiles.jar ";
		print PMERGE "TMP_DIR=/proj/".$scriptParameters{'projectID'}."/private/nobackup".'/$SLURM_JOB_ID '; #Temp Directory
		print PMERGE "OUTPUT=".$outSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,"_sorted_merged.bam "; #OutFile
	    }
	    
	    print PMERGE "INPUT=".$inSampleDirectory."/".$tempInfile."_sorted.bam "; #InFile
	}
	print PMERGE "\n\n";

	print PMERGE "samtools index ";
	print PMERGE $outSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } ,"_sorted_merged.bam", "\n\n"; #InFile
	print PMERGE "wait", "\n\n";

	print PMERGE "#Remove Temp Directory\n\n";
	print PMERGE "rm ";
	print PMERGE "-rf "; #Remove directory
	print PMERGE "/proj/".$scriptParameters{'projectID'}."/private/nobackup".'/$SLURM_JOB_ID ', "\n\n"; #Remove Temp Directory
    }
    if ( (@picardToolsMergePreviousFiles) && (scalar( @{ $InfilesLaneNoEnding{$sampleID} } ) > 1) ) { #merge previously merged files with merged files generated this run
	
	for (my $mergeFileCounter=0;$mergeFileCounter<scalar(@picardToolsMergePreviousFiles);$mergeFileCounter++) {
	    
	    if ($picardToolsMergePreviousFiles[$mergeFileCounter] =~ /$sampleID/) { #Look for sampleID in previously generated file to be merged with current run to be able to merge correct files
		if ($picardToolsMergePreviousFiles[$mergeFileCounter] =~ /lane(\d+)|s_(\d+)/) { #Look for lanes_ or lane\d in previously generated file to be merged with current run to be able to extract previous lanes
		    
		    my $mergeLanes; if($1) {$mergeLanes = $1;} else {$mergeLanes = $2;} #Make sure to always supply lanes from previous regexp
		    
		    print PMERGE "java -Xmx4g -jar ".$scriptParameters{'picardToolsPath'}."/MergeSamFiles.jar ";
		    print PMERGE "TMP_DIR=/proj/".$scriptParameters{'projectID'}."/private/nobackup/".'$SLURM_JOB_ID '; #Temp directory
		    print PMERGE "OUTPUT=".$outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,"_sorted_merged.bam "; #OutFile
		    print PMERGE "INPUT=".$inSampleDirectory."/".$sampleID."_lanes_", @{ $lanes{$sampleID} } , "_sorted_merged.bam "; #InFile
		    print PMERGE "INPUT=".$picardToolsMergePreviousFiles[$mergeFileCounter], "\n\n"; #$mergeLanes contains lane info on previous merge, $InfilesLaneNoEnding{$sampleID}[0] uses @RG for very first .bam file to include read group for subsequent merges. Complete path. 
		    
		    print PMERGE "samtools index ";
		    print PMERGE $outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,"_sorted_merged.bam ","\n\n"; #InFile

		    print PMERGE "#Remove Temp Directory\n\n";
		    print PMERGE "rm ";
		    print PMERGE "-rf "; #Remove directory
		    print PMERGE "/proj/".$scriptParameters{'projectID'}."/private/nobackup".'/$SLURM_JOB_ID ', "\n\n"; #Remove Temp Directory
		}
	    }
	}
    }
    elsif ( @picardToolsMergePreviousFiles ) { #merge previously merged files with single file generated this run
	
	for (my $mergeFileCounter=0;$mergeFileCounter<scalar(@picardToolsMergePreviousFiles);$mergeFileCounter++) {
	    
	    if ($picardToolsMergePreviousFiles[$mergeFileCounter] =~ /$sampleID/) { #Look for sampleID in previously generated file to be merged with current run to be able to merge correct files
		if ($picardToolsMergePreviousFiles[$mergeFileCounter] =~ /lane(\d+)|s_(\d+)/) { #Look for lanes_ or lane\d in previously generated file to be merged with current run to be able to extract previous lanes
		    
		    my $mergeLanes; if($1) {$mergeLanes = $1;} else {$mergeLanes = $2;} #Make sure to always supply lanes from previous regexp
		    my $tempInfile = $InfilesLaneNoEnding{$sampleID}[0]; #Can only be 1 element in array due to previous if statement
		    
		    print PMERGE "java -Xmx4g -jar ".$scriptParameters{'picardToolsPath'}."/MergeSamFiles.jar ";
		    print PMERGE "TMP_DIR=/proj/".$scriptParameters{'projectID'}."/private/nobackup/".'$SLURM_JOB_ID '; #Temp Directory
		    print PMERGE "OUTPUT=".$outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,"_sorted_merged.bam "; #OutFile
		    print PMERGE "INPUT=".$inSampleDirectory."/".$tempInfile."_sorted.bam "; #InFile
		    print PMERGE "INPUT=".$picardToolsMergePreviousFiles[$mergeFileCounter],"\n\n"; #$mergeLanes contains lane info on previous merge, $InfilesLaneNoEnding{$sampleID}[0] uses @RG for very first .bam file to include read group for subsequent merges. Complete path. 

		    print PMERGE "samtools index ";
		    print PMERGE $outSampleDirectory."/".$sampleID."_lanes_".$mergeLanes, @{ $lanes{$sampleID} } ,"_sorted_merged.bam", "\n\n"; #InFile

		    print PMERGE "#Remove Temp Directory\n\n";
		    print PMERGE "rm ";
		    print PMERGE "-rf "; #Remove directory
		    print PMERGE "/proj/".$scriptParameters{'projectID'}."/private/nobackup".'/$SLURM_JOB_ID ', "\n\n"; #Remove Temp Directory
		}
	    }
	}
    }
    close(PMERGE);
    if ( $scriptParameters{'pPicardToolsMerge'} == 1 ) {
	ParallelSampleIDSubmitJob($sampleID,$filename,"all");
    }
    return;
}

sub SamToolsSortIndex { 
#Sort and indexes bam files using samtools sort and samtools index

    my $sampleID = $_[0];
    my $aligner = $_[1];

    `mkdir -p $scriptParameters{'outDataDir'}/$sampleID/$aligner/info;`; #Creates the aligner folder and info data file directory
    `mkdir -p $scriptParameters{'outScriptDir'}/$sampleID/$aligner;`; #Creates the aligner folder and info data file directory

    my $sbatchScriptTracker=0;
    my $time=0;
    my $infileSize;
    for (my $infileCounter=0;$infileCounter<scalar( @{ $InfilesLaneNoEnding{$sampleID} } );$infileCounter++) { #For all files from

	
	if ( $scriptParameters{'pSamToolsSort'} ==1 ) {
	    $filename = $scriptParameters{'outScriptDir'}."/".$sampleID."/".$aligner."/samToolsSort_index_".$InfilesLaneNoEnding{$sampleID}[$infileCounter].".";
	}
	elsif ( $scriptParameters{'pSamToolsSort'} ==2 ) {
	    $filename = $scriptParameters{'outScriptDir'}."/".$sampleID."/".$aligner."/dry_run_samToolsSort_index_".$InfilesLaneNoEnding{$sampleID}[$infileCounter].".";
	    print STDOUT "Dry Run:\n";print MASTERL  "Dry Run:\n";
	}
	if ($infiles{$sampleID}[$infileCounter] =~/.fastq.gz$/) { #Files are already gz and presently the scalar for compression has not been investigated. Therefore no automatic time allocation can be performed.
	    if ($scriptParameters{'wholeGenomeSequencing'} == 1) {
		$time = 25;  
	    }
	    else {
		$time = 15;
	    }
	}
	else { #Files are in fastq format
	    $infileSize = -s $indirpath{$sampleID}."/".$infiles{$sampleID}[$infileCounter+$sbatchScriptTracker]; # collect .fastq file size to enable estimation of time required for sort & index, +$sbatchScriptTracker for syncing multiple infiles per sampleID. Hence, filesize will be calculated on read1 (should not matter).	   
	    
	    if ($scriptParameters{'pMosaikBuild'} || $scriptParameters{'pMosaikAlign'} || ($scriptParameters{'aligner'} eq "mosaik")) {
		$time = ceil($infileSize/(1700000*60*60)); #1700000 is a constant calculated from the filesize and time needed for procesing in samtools-0.1.12-10 sort and index and 60*60 is to scale to hours.
	    }
	    if ($scriptParameters{'pBwaAln'} || $scriptParameters{'pBwaSampe'} || ($scriptParameters{'aligner'} eq "bwa")) {
		$time = ceil($infileSize/(1700000*60*60)); #1700000 is a constant calculated from the filesize and time needed for procesing in samtools-0.1.12-10 sort and index and 60*60 is to scale to hours.	    
	    }
	}
	
	Checkfnexists($filename, $fnend);
	
#Info and Logg
	print STDOUT "Creating sbatch script SamTools sort & index and writing script file(s) to: ".$filename, "\n";print MASTERL "Creating sbatch script SamTools sort & index and writing script file(s) to: ".$filename, "\n";
	print STDOUT "Sbatch script SamTools sort & index data files will be written to: ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner, "\n";print MASTERL "Sbatch script SamTools sort & index data files will be written to: ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner, "\n";
	
	open (STSI, ">".$filename) or die "Can't write to ".$filename.": $!\n";
	
	print STSI "#! /bin/bash -l", "\n";
	print STSI "#SBATCH -A ".$scriptParameters{'projectID'}, "\n";
	print STSI "#SBATCH -p node -n 1", "\n";
	print STSI "#SBATCH -C thin", "\n";	
	print STSI "#SBATCH -t ".$time.":00:00", "\n";
	
	print STSI "#SBATCH -J STSI_".$sampleID."_".$aligner, "\n";
	if ( $scriptParameters{'pSamToolsSort'} ==1 ) {
	    print STSI "#SBATCH -e ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/info/samToolsSort_index_".$InfilesLaneNoEnding{$sampleID}[$infileCounter].".".$fnt.".stderr.txt", "\n";
	    print STSI "#SBATCH -o ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/info/samToolsSort_index_".$InfilesLaneNoEnding{$sampleID}[$infileCounter].".".$fnt.".stdout.txt", "\n";
	}
	elsif ( $scriptParameters{'pSamToolsSort'} ==1 ) {
	    print STSI "#SBATCH -e ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/info/dry_run_samToolsSort_index_".$InfilesLaneNoEnding{$sampleID}[$infileCounter].".".$fnt.".stderr.txt", "\n";
	    print STSI "#SBATCH -o ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner."/info/dry_run_samToolsSort_index_".$InfilesLaneNoEnding{$sampleID}[$infileCounter].".".$fnt.".stdout.txt", "\n";
	}
    
	unless ($scriptParameters{'email'} eq 0) {
	    print STSI "#SBATCH --mail-type=END", "\n";
	    print STSI "#SBATCH --mail-type=FAIL", "\n";
	    print STSI "#SBATCH --mail-user=".$scriptParameters{'email'}, "\n\n";
	}
	
	print STSI 'echo "Running on: $(hostname)"',"\n\n";
    
###	
#SamTools Sort
###	
	my $inSampleDirectory = $scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner;
	my $outSampleDirectory = $scriptParameters{'outDataDir'}."/".$sampleID."/".$aligner;
	my $tempInfile = $InfilesLaneNoEnding{$sampleID}[$infileCounter];

	print STSI "samtools sort ";
	print STSI $inSampleDirectory."/".$tempInfile.".bam ".$outSampleDirectory."/".$tempInfile."_sorted", "\n\n"; #InFile. SamTools sort adds .bam ending
	
	print STSI "wait", "\n\n";
###	
#SamTools Index
###	
	print STSI "samtools index ";
	print STSI $inSampleDirectory."/".$tempInfile."_sorted.bam", "\n\n"; #InFile
	$sbatchScriptTracker++;
	close(STSI);
	if ( $scriptParameters{'pSamToolsSort'} ==1 ) {
	    ParallelSampleIDSubmitJob($sampleID,$filename,$InfilesLaneNoEnding{$sampleID}[$infileCounter]);
	} 
    }
    return;
}

sub BWA_Sampe {
#Alignments of BWA Aln index reads using BWA sampe
    
    my $sampleID = $_[0];

    `mkdir -p $scriptParameters{'outDataDir'}/$sampleID/bwa/info;`; #Creates the bwa folder and info data file directory
    `mkdir -p $scriptParameters{'outScriptDir'}/$sampleID/bwa;`; #Creates the bwa script directory
    
    my $sbatchScriptTracker=0;
    my $time=0;
    my $infileSize;
    for (my $infileCounter=0;$infileCounter<( scalar( @{ $InfilesLaneNoEnding{$sampleID} }) );$infileCounter++) { #For all files from BWA aln but process in the same command i.e. both reads per align call
	if ($infiles{$sampleID}[$infileCounter] =~/.fastq.gz$/) { #Files are already gz and presently the scalar for compression has not been investigated. Therefore no automatic time allocation can be performed.
	    if ($scriptParameters{'wholeGenomeSequencing'} == 1) {
		$time = 60;  
	    }
	    else {
		$time = 30;
	    }
	}
	else { #Files are in fastq format	
	    $infileSize = -s $indirpath{$sampleID}."/".$infiles{$sampleID}[$infileCounter+$sbatchScriptTracker]; # collect .fastq file size to enable estimation of time required for aligning, +1 for syncing multiple infiles per sampleID. Hence, filesize will be calculated on read2 (should not matter).
	    $time = ceil(($infileSize/238)/(3000*60*60)); #238 is a scalar estimating the number of reads depending on filesize. 3500 is the number of reads/s in Bwa_sampe-0.6.1 plus samtools-0.1.12-10 view sam to bam conversion and 60*60 is to scale to hours. (4600 BWA-0.5.9)
	}
	if ( $scriptParameters{'pBwaSampe'} == 1 ) {
	    $filename = $scriptParameters{'outScriptDir'}."/".$sampleID."/bwa/bwa_sampe_".$InfilesLaneNoEnding{$sampleID}[$infileCounter].".";	
	}
	elsif ( $scriptParameters{'pBwaSampe'} == 2) { #Dry run
	    $filename = $scriptParameters{'outScriptDir'}."/".$sampleID."/bwa/dry_run_bwa_sampe_".$InfilesLaneNoEnding{$sampleID}[$infileCounter].".";
	    print STDOUT "Dry Run:\n";print MASTERL  "Dry Run:\n";
	}
	Checkfnexists($filename, $fnend);

#Info and Logg
	print STDOUT "Creating sbatch script BWA_Sampe and writing script file(s) to: ".$filename, "\n";print MASTERL "Creating sbatch script BWA_Sampe and writing script file(s) to: ".$filename, "\n";
	print STDOUT "Sbatch script BWA_Sampe data files will be written to: ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$scriptParameters{'aligner'}, "\n";print MASTERL "Sbatch script BWA_Sampe data files will be written to: ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$scriptParameters{'aligner'}, "\n";
	
	open (BWAS, ">".$filename) or die "Can't write to ".$filename.": $!\n";
	
	print BWAS "#! /bin/bash -l", "\n";
	print BWAS "#SBATCH -A ".$scriptParameters{'projectID'}, "\n";
	print BWAS "#SBATCH -p node -n ".$scriptParameters{'maximumCores'}, "\n";
	print BWAS "#SBATCH -C thin", "\n";
	print BWAS "#SBATCH -t ".$time.":00:00", "\n";
	print BWAS "#SBATCH -J BWA_Sampe_".$sampleID, "\n";
	if ( $scriptParameters{'pBwaSampe'} == 1 ) {
	    print BWAS "#SBATCH -e ".$scriptParameters{'outDataDir'}."/".$sampleID."/bwa/info/bwa_sampe_".$InfilesLaneNoEnding{$sampleID}[$infileCounter].".".$fnt.".stderr.txt", "\n";
	    print BWAS "#SBATCH -o ".$scriptParameters{'outDataDir'}."/".$sampleID."/bwa/info/bwa_sampe_".$InfilesLaneNoEnding{$sampleID}[$infileCounter].".".$fnt.".stdout.txt", "\n";
	}
	elsif ( $scriptParameters{'pBwaSampe'} == 2 ) {
	    print BWAS "#SBATCH -e ".$scriptParameters{'outDataDir'}."/".$sampleID."/bwa/info/dry_run_bwa_sampe_".$InfilesLaneNoEnding{$sampleID}[$infileCounter].".".$fnt.".stderr.txt", "\n";
	    print BWAS "#SBATCH -o ".$scriptParameters{'outDataDir'}."/".$sampleID."/bwa/info/dry_run_bwa_sampe_".$InfilesLaneNoEnding{$sampleID}[$infileCounter].".".$fnt.".stdout.txt", "\n";
	}
	
	unless ($scriptParameters{'email'} eq 0) {	    
	    print BWAS "#SBATCH --mail-type=END", "\n";
	    print BWAS "#SBATCH --mail-type=FAIL", "\n";
	    print BWAS "#SBATCH --mail-user=".$scriptParameters{'email'}, "\n\n";
	}
	
	print BWAS 'echo "Running on: $(hostname)"',"\n\n";
	
	my $BWAinSampleDirectory = $scriptParameters{'outDataDir'}."/".$sampleID."/bwa";
	my $FASTQinSampleDirectory = $indirpath{$sampleID};
	my $outSampleDirectory = $scriptParameters{'outDataDir'}."/".$sampleID."/bwa";
	my $tempInfile = $infiles{$sampleID}[$infileCounter+$sbatchScriptTracker]; #For required .fastq file
	my $tempInfile2 = $infiles{$sampleID}[ ($infileCounter+$sbatchScriptTracker+1)]; # #For required .fastq file (Paired read)   

#BWA Sampe	
	print BWAS "bwa sampe ";
	print BWAS "-r ".'"@RG\tID:'.$InfilesBothStrandsNoEnding{$sampleID}[$infileCounter+$sbatchScriptTracker].'\tSM:'.$sampleID.'\tPL:ILLUMINA" '.$scriptParameters{'referencesDir'}."/".$scriptParameters{'humanGenomeReference'}." "; #read group header line
	print BWAS $BWAinSampleDirectory."/".$InfilesBothStrandsNoEnding{$sampleID}[$infileCounter+$sbatchScriptTracker].".sai "; #Read 1
	print BWAS $BWAinSampleDirectory."/".$InfilesBothStrandsNoEnding{$sampleID}[ ($infileCounter+$sbatchScriptTracker+1) ].".sai "; #Read 2
	print BWAS $FASTQinSampleDirectory."/".$tempInfile." "; #Fastq read 1
	print BWAS $FASTQinSampleDirectory."/".$tempInfile2." "; #Fastq read 2
	print BWAS "> ".$outSampleDirectory."/".$InfilesLaneNoEnding{$sampleID}[$infileCounter].".sam", "\n\n"; #Outfile (SAM)

#Convert SAM to BAM using samTools view	
	print BWAS "samtools view -bS ".$BWAinSampleDirectory."/".$InfilesLaneNoEnding{$sampleID}[$infileCounter].".sam "; #Infile (SAM)
	print BWAS "> ".$outSampleDirectory."/".$InfilesLaneNoEnding{$sampleID}[$infileCounter].".bam", "\n\n"; #Outfile (BAM)

#Remove SAM file
	print BWAS "Removing temporary SAM-file\n";
	print BWAS "rm ".$BWAinSampleDirectory."/".$InfilesLaneNoEnding{$sampleID}[$infileCounter].".sam";
	$sbatchScriptTracker++;	
	close(BWAS);
	if ( $scriptParameters{'pBwaSampe'} == 1 ) {
	    ParallelSampleIDSubmitJob($sampleID,$filename,$InfilesLaneNoEnding{$sampleID}[$infileCounter]);
	}
    }
    return;
}

sub BWA_Aln {
#Generates BWA aln index on fastq files
    
    my $sampleID = $_[0];

    `mkdir -p $scriptParameters{'outDataDir'}/$sampleID/bwa/info;`; #Creates the bwa folder and info data file directory
    `mkdir -p $scriptParameters{'outScriptDir'}/$sampleID/bwa;`; #Creates the bwa script directory
    if ( $scriptParameters{'pBwaAln'} == 1 ) {
	$filename = $scriptParameters{'outScriptDir'}."/".$sampleID."/bwa/bwa_aln_".$sampleID.".";
    }
    elsif ( $scriptParameters{'pBwaAln'} == 2 ) {
	$filename = $scriptParameters{'outScriptDir'}."/".$sampleID."/bwa/dry_run_bwa_aln_".$sampleID.".";
	print STDOUT "Dry Run:\n";print MASTERL  "Dry Run:\n";
    }
    Checkfnexists($filename, $fnend);

#Info and Logg
    print STDOUT "Creating sbatch script BWA_Aln and writing script file(s) to: ".$filename, "\n";print MASTERL "Creating sbatch script BWA_Aln and writing script file(s) to: ".$filename, "\n";
    print STDOUT "Sbatch script BWA_Aln data files will be written to: ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$scriptParameters{'aligner'}, "\n";print MASTERL "Sbatch script BWA_Aln data files will be written to: ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$scriptParameters{'aligner'}, "\n";
    
    open (BWAA, ">".$filename) or die "Can't write to ".$filename.": $!\n";
    
    my $time = ceil(2.5*scalar( @{ $InfilesLaneNoEnding{$sampleID} })); #One full lane on Hiseq takes approx. 2,5 h for BWA_Aln to process, round up to nearest full hour.
    
    print BWAA "#! /bin/bash -l", "\n";
    print BWAA "#SBATCH -A ".$scriptParameters{'projectID'}, "\n";
    print BWAA "#SBATCH -p node -n ".$scriptParameters{'maximumCores'}, "\n";
    print BWAA "#SBATCH -C thin", "\n";
    print BWAA "#SBATCH -t ".$time.":00:00", "\n";
    print BWAA "#SBATCH -J BWA_Aln_".$sampleID, "\n";
    if ( $scriptParameters{'pBwaAln'} == 1 ) {
	print BWAA "#SBATCH -e ".$scriptParameters{'outDataDir'}."/".$sampleID."/bwa/info/bwa_aln_".$sampleID.".".$fnt.".stderr.txt", "\n";
	print BWAA "#SBATCH -o ".$scriptParameters{'outDataDir'}."/".$sampleID."/bwa/info/bwa_aln_".$sampleID.".".$fnt.".stdout.txt", "\n";
    }
    elsif ( $scriptParameters{'pBwaAln'} == 2 ) { #Dry run
	print BWAA "#SBATCH -e ".$scriptParameters{'outDataDir'}."/".$sampleID."/bwa/info/dry_run_bwa_aln_".$sampleID.".".$fnt.".stderr.txt", "\n";
	print BWAA "#SBATCH -o ".$scriptParameters{'outDataDir'}."/".$sampleID."/bwa/info/dry_run_bwa_aln_".$sampleID.".".$fnt.".stdout.txt", "\n";
    }
    unless ($scriptParameters{'email'} eq 0) {
	print BWAA "#SBATCH --mail-type=END", "\n";
	print BWAA "#SBATCH --mail-type=FAIL", "\n";
	print BWAA "#SBATCH --mail-user=".$scriptParameters{'email'}, "\n\n";
    }
    
    print BWAA 'echo "Running on: $(hostname)"',"\n\n";
    
    my $inSampleDirectory =  $indirpath{$sampleID};
    my $outSampleDirectory = $scriptParameters{'outDataDir'}."/".$sampleID."/bwa";
    my $coreCounter=1;    
    for (my $infileCounter=0;$infileCounter<scalar( @{ $infiles{$sampleID} });$infileCounter++) {
	
	if ($infileCounter eq $coreCounter*$scriptParameters{'maximumCores'}) { #Using only $scriptParameters{'maximumCores'} cores
	    
	    print BWAA "wait", "\n\n";
	    $coreCounter=$coreCounter+1;
	}

	my $tempinfile = $infiles{$sampleID}[$infileCounter];

	print BWAA "bwa aln ";
	print BWAA "-k 1 "; #maximum differences in the seed
	print BWAA "-t 4 "; #number of threads
	print BWAA "-n 3 "; #max #diff (int) or missing prob under 0.02 err rate (float)
	print BWAA "-q ".$scriptParameters{'bwaAlnQualityTrimming'}." "; #Quality trimming
	print BWAA $scriptParameters{'referencesDir'}."/".$scriptParameters{'humanGenomeReference'}." "; #Reference
	print BWAA $inSampleDirectory."/".$tempinfile." "; #InFile
	print BWAA "> ".$outSampleDirectory."/".$InfilesBothStrandsNoEnding{$sampleID}[$infileCounter].".sai &", "\n\n"; #OutFile 
    }
    print BWAA "wait", "\n\n";
    close(BWAA);
    if ( $scriptParameters{'pBwaAln'} == 1 ) {
	SampleIDSubmitJob($sampleID,$filename, 1);   
    }
    return;
}

sub MosaikAlign {
#Aligning reads using MosaikAlign
    
    my $sampleID = $_[0];
    
    `mkdir -p $scriptParameters{'outDataDir'}/$sampleID/mosaik/info;`; #Creates the mosaik folder and info data file directory
    `mkdir -p $scriptParameters{'outScriptDir'}/$sampleID/mosaik;`; #Creates the mosaik script directory
    
    my $sbatchScriptTracker=0;
    my $time=0;
    my $infileSize;
    for (my $infileCounter=0;$infileCounter<scalar( @{ $InfilesLaneNoEnding{$sampleID} } );$infileCounter++) { #For all files from MosaikBuild (platform reads)
	if ($infiles{$sampleID}[$infileCounter] =~/.fastq.gz$/) { #Files are already gz and presently the scalar for compression has not been investigated. Therefore no automatic time allocation can be performed.
	    if ($scriptParameters{'wholeGenomeSequencing'} == 1) {
		$time = 80;  
	    }
	    else {
		$time = 40;
	    }
	}
	else { #Files are in fastq format
	    if (-e $indirpath{$sampleID}."/".$infiles{$sampleID}[$infileCounter+$sbatchScriptTracker]) {
		$infileSize = -s $indirpath{$sampleID}."/".$infiles{$sampleID}[$infileCounter+$sbatchScriptTracker]; # collect .fastq file size to enable estimation of time required for aligning, +$sbatchScriptTracker for syncing multiple infiles per sampleID. Hence, filesize will be calculated on read1 (should not matter).      
		$time = ceil(($infileSize/238)/(650*60*60)); #238 is a scalar estimating the number of reads depending on filesize. 650 is the number of reads/s in MosaikAlign-2.1.52 and 60*60 is to scale to hours.
	    }	    
	} 
	if ( $scriptParameters{'pMosaikAlign'} == 1 ) {
	    $filename = $scriptParameters{'outScriptDir'}."/".$sampleID."/mosaik/mosaikAlign_".$InfilesLaneNoEnding{$sampleID}[$infileCounter].".";
	}
	elsif ( $scriptParameters{'pMosaikAlign'} == 2 ) { #Dry run
	    $filename = $scriptParameters{'outScriptDir'}."/".$sampleID."/mosaik/dry_run_mosaikAlign_".$InfilesLaneNoEnding{$sampleID}[$infileCounter].".";
	    print STDOUT "Dry Run:\n";print MASTERL  "Dry Run:\n";
	}
	Checkfnexists($filename, $fnend);
	
#Info and Logg
	print STDOUT "Creating sbatch script MosaikAlign and writing script file(s) to: ".$filename, "\n";print MASTERL "Creating sbatch script MosaikAlign and writing script file(s) to: ".$filename, "\n";
	print STDOUT "Sbatch script MosaikAlign data files will be written to: ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$scriptParameters{'aligner'}, "\n";print MASTERL "Sbatch script MosaikAlign data files will be written to: ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$scriptParameters{'aligner'}, "\n";
	
	open (MOSA, ">".$filename) or die "Can't write to ".$filename.": $!\n";
	
	print MOSA "#! /bin/bash -l", "\n";
	print MOSA "#SBATCH -A ".$scriptParameters{'projectID'}, "\n";
	print MOSA "#SBATCH -p node -n ".$scriptParameters{'maximumCores'}, "\n";
	print MOSA "#SBATCH -C thin", "\n";
	print MOSA "#SBATCH -t ".$time.":00:00", "\n";
	print MOSA "#SBATCH -J MoA_".$sampleID."_".$sbatchScriptTracker, "\n";
	if ( $scriptParameters{'pMosaikAlign'} == 1 ) {
	    print MOSA "#SBATCH -e ".$scriptParameters{'outDataDir'}."/".$sampleID."/mosaik/info/mosaikAlign_".$InfilesLaneNoEnding{$sampleID}[$infileCounter].".".$fnt.".stderr.txt", "\n";
	    print MOSA "#SBATCH -o ".$scriptParameters{'outDataDir'}."/".$sampleID."/mosaik/info/mosaikAlign_".$InfilesLaneNoEnding{$sampleID}[$infileCounter].".".$fnt.".stdout.txt", "\n";
	}
	elsif ( $scriptParameters{'pMosaikAlign'} == 2 ) { #Dry run
	    print MOSA "#SBATCH -e ".$scriptParameters{'outDataDir'}."/".$sampleID."/mosaik/info/dry_run_mosaikAlign_".$InfilesLaneNoEnding{$sampleID}[$infileCounter].".".$fnt.".stderr.txt", "\n";
	    print MOSA "#SBATCH -o ".$scriptParameters{'outDataDir'}."/".$sampleID."/mosaik/info/dry_run_mosaikAlign_".$InfilesLaneNoEnding{$sampleID}[$infileCounter].".".$fnt.".stdout.txt", "\n";
	}
	unless ($scriptParameters{'email'} eq 0) {
	    print MOSA "#SBATCH --mail-type=END", "\n";
	    print MOSA "#SBATCH --mail-type=FAIL", "\n";
	    print MOSA "#SBATCH --mail-user=".$scriptParameters{'email'}, "\n\n";
	}
	
	print MOSA 'echo "Running on: $(hostname)"',"\n\n";
	print MOSA "mkdir -p /scratch/mosaik_tmp", "\n";
	print MOSA "export MOSAIK_TMP=/scratch/mosaik_tmp", "\n\n";

	my $inSampleDirectory = $scriptParameters{'outDataDir'}."/".$sampleID."/mosaik";
	my $outSampleDirectory = $scriptParameters{'outDataDir'}."/".$sampleID."/mosaik";
	my $tempInfile = $InfilesLaneNoEnding{$sampleID}[$infileCounter];

	print MOSA "MosaikAligner ";
	print MOSA "-in ".$inSampleDirectory."/".$tempInfile.".dat "; #Infile
	print MOSA "-out ".$outSampleDirectory."/".$InfilesLaneNoEnding{$sampleID}[$infileCounter]." "; #OutFile
	print MOSA "-ia ".$scriptParameters{'referencesDir'}."/".$scriptParameters{'mosaikAlignReference'}." "; #Mosaik Reference
	print MOSA "-annpe ".$scriptParameters{'referencesDir'}."/".$scriptParameters{'mosaikAlignNeuralNetworkPeFile'}." "; #NerualNetworkPE
	print MOSA "-annse ".$scriptParameters{'referencesDir'}."/".$scriptParameters{'mosaikAlignNeuralNetworkSeFile'}." "; #NerualNetworkSE
	print MOSA "-hs 15 "; #hash size
	print MOSA "-mm 4 "; #the # of mismatches allowed
	print MOSA "-mhp 100 "; #the maximum # of positions stored per seed
	print MOSA "-ls 100 "; #enable local alignment search for PE reads
	print MOSA "-act 35 "; #the alignment candidate threshold (length)
	print MOSA "-bw 35 "; #specifies the Smith-Waterman bandwidth.
	print MOSA "-j ".$scriptParameters{'referencesDir'}."/".$scriptParameters{'mosaikJumpDbStub'}." "; #JumpDatabase
	print MOSA "-p ".$scriptParameters{'maximumCores'}, "\n\n"; #Nr of cores
	
	$sbatchScriptTracker++; #Tracks nr of sbatch scripts
	close(MOSA);
	if ( $scriptParameters{'pMosaikAlign'} == 1 ) {
	    ParallelSampleIDSubmitJob($sampleID,$filename, $InfilesLaneNoEnding{$sampleID}[$infileCounter]);
	}
    }
    return;
}

sub MosaikBuild {
#Generates Mosaik hash format on reads using MosaikBuild   
    
    my $sampleID = $_[0];
    `mkdir -p $scriptParameters{'outDataDir'}/$sampleID/mosaik/info;`; #Creates the mosaik folder and info data file directory
    `mkdir -p $scriptParameters{'outScriptDir'}/$sampleID/mosaik;`; #Creates the mosaik script directory

    if ($scriptParameters{'pMosaikBuild'} == 1) {	
	$filename = $scriptParameters{'outScriptDir'}."/".$sampleID."/mosaik/mosaikBuild_".$sampleID.".";    
    }
    elsif ( $scriptParameters{'pMosaikBuild'} == 2 ) {
	$filename = $scriptParameters{'outScriptDir'}."/".$sampleID."/mosaik/dry_run_mosaikBuild_".$sampleID.".";
	print STDOUT "Dry Run:\n";print MASTERL  "Dry Run:\n";
    }
    Checkfnexists($filename, $fnend);
    
#Info and Logg
    print STDOUT "Creating sbatch script MosaikBuild and writing script file(s) to: ".$filename, "\n";print MASTERL "Creating sbatch script MosaikBuild and writing script file(s) to: ".$filename, "\n";
    print STDOUT "Sbatch script MosaikBuild data files will be written to: ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$scriptParameters{'aligner'}, "\n";print MASTERL "Sbatch script MosaikBuild data files will be written to: ".$scriptParameters{'outDataDir'}."/".$sampleID."/".$scriptParameters{'aligner'}, "\n";

    my $time = ceil(2.5*scalar( @{ $InfilesLaneNoEnding{$sampleID} })); #One full lane on Hiseq takes approx. 1 h for MosaikBuild to process (compressed format, uncompressed 0.5 h), round up to nearest full hour.
    
    open (MOSB, ">".$filename) or die "Can't write to ".$filename.": $!\n";
    
    print MOSB "#! /bin/bash -l", "\n";
    print MOSB "#SBATCH -A ".$scriptParameters{'projectID'}, "\n";
    print MOSB "#SBATCH -p node -n ".$scriptParameters{'maximumCores'}, "\n";
    print MOSB "#SBATCH -C thin", "\n";
    print MOSB "#SBATCH -t ".$time.":00:00", "\n";
    print MOSB "#SBATCH -J MosB_".$sampleID, "\n";
    if ($scriptParameters{'pMosaikBuild'} == 1) {
	print MOSB "#SBATCH -e ".$scriptParameters{'outDataDir'}."/".$sampleID."/mosaik/info/mosaikBuild_".$sampleID.".".$fnt.".stderr.txt", "\n";
	print MOSB "#SBATCH -o ".$scriptParameters{'outDataDir'}."/".$sampleID."/mosaik/info/mosaikBuild_".$sampleID.".".$fnt.".stdout.txt", "\n";
    }
    elsif ( $scriptParameters{'pMosaikBuild'} == 2 ) {
	print MOSB "#SBATCH -e ".$scriptParameters{'outDataDir'}."/".$sampleID."/mosaik/info/dry_run_mosaikBuild_".$sampleID.".".$fnt.".stderr.txt", "\n";
	print MOSB "#SBATCH -o ".$scriptParameters{'outDataDir'}."/".$sampleID."/mosaik/info/dry_run_mosaikBuild_".$sampleID.".".$fnt.".stdout.txt", "\n";
    }
    unless ($scriptParameters{'email'} eq 0) {
	print MOSB "#SBATCH --mail-type=END", "\n";
	print MOSB "#SBATCH --mail-type=FAIL", "\n";
	print MOSB "#SBATCH --mail-user=".$scriptParameters{'email'}, "\n\n";
    }
    
    print MOSB 'echo "Running on: $(hostname)"',"\n\n";
    
    my $inSampleDirectory = $indirpath{$sampleID};
    my $outSampleDirectory = $scriptParameters{'outDataDir'}."/".$sampleID."/mosaik";
    my $coreCounter=1;
    my $coreTracker=0; #Required to portion out cores and files before wait and to track the MOSB outfiles to correct lane
    
    for (my $infileCounter=0;$infileCounter<( scalar( @{ $infiles{$sampleID} }) -1);$infileCounter++) {
	
	if ($coreTracker eq $coreCounter*$scriptParameters{'maximumCores'}) { #Using only $scriptParameters{'maximumCores'} nr of cores
	    
	    print MOSB "wait", "\n\n";
	    $coreCounter=$coreCounter+1;
	}
	my $tempInfile = $infiles{$sampleID}[$infileCounter];
	my $tempInfile2 = $infiles{$sampleID}[ ($infileCounter+1)]; #Paired read
	$infileCounter = $infileCounter+1; #To correct for reading 2 files at once

	print MOSB "MosaikBuild ";
	print MOSB "-id ".$InfilesBothStrandsNoEnding{$sampleID}[$infileCounter]." "; #Read group ID for BAM Header
	print MOSB "-sam ".$sampleID." "; #Sample name for BAM Header
	print MOSB "-st illumina_long "; #Sequencing technology for BAM Header
	print MOSB "-mfl ".$scriptParameters{'mosaikBuildMedianFragLength'}." "; #Median Fragment Length
	print MOSB "-q ".$inSampleDirectory."/".$tempInfile." "; #Read 1
	print MOSB "-q2 ".$inSampleDirectory."/".$tempInfile2." "; #Read 2
	print MOSB "-out ".$outSampleDirectory."/".$InfilesLaneNoEnding{$sampleID}[($coreTracker)].".dat &", "\n\n"; #OutFile
	$coreTracker++; #Track nr of mosaikBuild calls so that wait can be printed at the correct intervals (dependent on $scriptParameters{'maximumCores'})
    }
    print MOSB "wait", "\n\n";    
    close(MOSB);
    if ($scriptParameters{'pMosaikBuild'} == 1) {
	SampleIDSubmitJob($sampleID,$filename, 1); 
    }
    return;
}   

sub FastQC {
#Raw sequence quality analysis using FASTQC

    my $sampleID = $_[0];
    `mkdir -p $scriptParameters{'outDataDir'}/$sampleID/fastqc/info;`; #Creates the fastqc folder and info data file directory
    `mkdir -p $scriptParameters{'outScriptDir'}/$sampleID/fastqc;`; #Creates the fastqc script directory
    if ($scriptParameters{'pFastQC'} == 1) {
	$filename = $scriptParameters{'outScriptDir'}."/".$sampleID."/fastqc/fastqc_".$sampleID.".";
    }
    elsif ($scriptParameters{'pFastQC'} == 2) { #Dry run
	$filename = $scriptParameters{'outScriptDir'}."/".$sampleID."/fastqc/dry_run_fastqc_".$sampleID.".";
	print STDOUT "Dry Run:\n";print MASTERL  "Dry Run:\n";
    }
    Checkfnexists($filename, $fnend);

#Info and Logg
    print STDOUT "Creating sbatch script FastQC and writing script file(s) to: ".$filename, "\n";print MASTERL "Creating sbatch script Sample check FastQC and writing script file(s) to: ".$filename, "\n";
    print STDOUT "Sbatch script FastQC data files will be written to: ".$scriptParameters{'outDataDir'}."/".$sampleID."/fastqc", "\n";print MASTERL "Sbatch script Sample check FastQC data files will be written to: ".$scriptParameters{'outDataDir'}."/".$sampleID."/fastqc", "\n";

    my $time = ceil(0.5*scalar( @{ $infiles{$sampleID} })); #One full lane on Hiseq takes approx. 0.5 h for FASTQC to process, round up to nearest full hour.
    
    open (FASTQC, ">".$filename) or die "Can't write to ".$filename.": $!\n";
    
    print FASTQC "#! /bin/bash -l", "\n";
    print FASTQC "#SBATCH -A ".$scriptParameters{'projectID'}, "\n";
    print FASTQC "#SBATCH -p node -n ".$scriptParameters{'maximumCores'}, "\n";
    print FASTQC "#SBATCH -C thin", "\n";
    print FASTQC "#SBATCH -t ".$time.":00:00", "\n";
    print FASTQC "#SBATCH -J FQC_".$sampleID, "\n";
    if ($scriptParameters{'pFastQC'} == 1) {
	print FASTQC "#SBATCH -e ".$scriptParameters{'outDataDir'}."/".$sampleID."/fastqc/info/fastqc_".$sampleID.".".$fnt.".stderr.txt", "\n";
	print FASTQC "#SBATCH -o ".$scriptParameters{'outDataDir'}."/".$sampleID."/fastqc/info/fastqc_".$sampleID.".".$fnt.".stdout.txt", "\n";
    }
    elsif ($scriptParameters{'pFastQC'} == 2) { #Dry run
	print FASTQC "#SBATCH -e ".$scriptParameters{'outDataDir'}."/".$sampleID."/fastqc/info/dry_run_fastqc_".$sampleID.".".$fnt.".stderr.txt", "\n";
	print FASTQC "#SBATCH -o ".$scriptParameters{'outDataDir'}."/".$sampleID."/fastqc/info/dry_run_fastqc_".$sampleID.".".$fnt.".stdout.txt", "\n";
    }
    unless ($scriptParameters{'email'} eq 0) {
	print FASTQC "#SBATCH --mail-type=END", "\n";
	print FASTQC "#SBATCH --mail-type=FAIL", "\n";
	print FASTQC "#SBATCH --mail-user=".$scriptParameters{'email'}, "\n\n";
    }
    
    print FASTQC 'echo "Running on: $(hostname)"',"\n\n";
    print FASTQC "cd ".$indirpath{$sampleID}, "\n\n";
    
    my $outSampleDirectory = $scriptParameters{'outDataDir'}."/".$sampleID."/fastqc";
    my $coreCounter=1;
    for (my $infileCounter=0;$infileCounter<scalar( @{ $infiles{$sampleID} });$infileCounter++) {
	
	if ($infileCounter eq $coreCounter*$scriptParameters{'maximumCores'}) { #Using only $scriptParameters{'maximumCores'} cores
	    
	    print FASTQC "wait", "\n\n";
	    $coreCounter=$coreCounter+1;
	}
	my $tempInfile = $infiles{$sampleID}[$infileCounter];
	print FASTQC "fastqc ";
	print FASTQC $tempInfile." "; #InFile
	print FASTQC "-o ".$outSampleDirectory. " &", "\n\n"; #OutFile
    }
    print FASTQC "wait", "\n";    
    
    close(FASTQC);
    if ($scriptParameters{'pFastQC'} == 1) {
	SampleIDSubmitJob($sampleID,$filename, 0);
    }
    return;
}

sub GZipfastq { 
#Automatically gzips fastq files. 
    
    my $sampleID = $_[0];
    `mkdir -p $scriptParameters{'outDataDir'}/$sampleID/gzip/info;`; #Creates the gzip folder and info data file directory
    `mkdir -p $scriptParameters{'outScriptDir'}/$sampleID/gzip;`; #Creates the gzip script folder 
    if ($scriptParameters{'pGZip'} == 1) {
	$filename = $scriptParameters{'outScriptDir'}."/".$sampleID."/gzip/gzipFastq_".$sampleID.".";
    }
    elsif ($scriptParameters{'pGZip'} == 2) { #Dry run  
	$filename = $scriptParameters{'outScriptDir'}."/".$sampleID."/gzip/dry_run_gzipFastq_".$sampleID.".";
	print STDOUT "Dry Run:\n";print MASTERL  "Dry Run:\n";
    }
    Checkfnexists($filename, $fnend);

#Info and Logg
    print STDOUT "Creating sbatch script GzipFastq and writing script file(s) to: ".$filename, "\n";print MASTERL "Creating sbatch script GzipFastq and writing script file(s) to: ".$filename, "\n";
    print STDOUT "Sbatch script GzipFastq data files will be written to: ".$scriptParameters{'outDataDir'}."/".$sampleID."/fastq", "\n";print MASTERL "Sbatch script GzipFastq data files will be written to: ".$scriptParameters{'outDataDir'}."/".$sampleID."/fastq", "\n";
    
    my $time = ceil(1.5*scalar( @{ $infiles{$sampleID} })); #One full lane on Hiseq takes approx. 1.5 h for gzip to process, round up to nearest full hour.
    open (GZFASTQ, ">".$filename) or die "Can't write to ".$filename.": .$!", "\n";
    
    print GZFASTQ "#! /bin/bash -l", "\n";
    print GZFASTQ "#SBATCH -A ".$scriptParameters{'projectID'}, "\n";
    print GZFASTQ "#SBATCH -p node -n ".$scriptParameters{'maximumCores'}, "\n";
    print GZFASTQ "#SBATCH -C thin", "\n";	
    print GZFASTQ "#SBATCH -t ".$time.":00:00", "\n";
    print GZFASTQ "#SBATCH -J GZFQ_".$sampleID, "\n";
    if ($scriptParameters{'pGZip'} == 1) {
	print GZFASTQ "#SBATCH -e ".$scriptParameters{'outDataDir'}."/".$sampleID."/gzip/info/gzipFastq_".$sampleID.".".$fnt.".stderr.txt", "\n";
	print GZFASTQ "#SBATCH -o ".$scriptParameters{'outDataDir'}."/".$sampleID."/gzip/info/gzipFastq_".$sampleID.".".$fnt.".stdout.txt", "\n";
    }
    elsif ($scriptParameters{'pGZip'} == 2) { #Dry run
	print GZFASTQ "#SBATCH -e ".$scriptParameters{'outDataDir'}."/".$sampleID."/gzip/info/dry_run_gzipFastq_".$sampleID.".".$fnt.".stderr.txt", "\n";
	print GZFASTQ "#SBATCH -o ".$scriptParameters{'outDataDir'}."/".$sampleID."/gzip/info/dry_run_gzipFastq_".$sampleID.".".$fnt.".stdout.txt", "\n";
    }
    unless ($scriptParameters{'email'} eq 0) {
	
	print GZFASTQ "#SBATCH --mail-type=END", "\n";
	print GZFASTQ "#SBATCH --mail-type=FAIL", "\n";
	print GZFASTQ "#SBATCH --mail-user=".$scriptParameters{'email'}, "\n\n";
	
    }
    
    print GZFASTQ 'echo "Running on: $(hostname)"',"\n\n";
    print GZFASTQ "cd ".$indirpath{$sampleID}, "\n\n";
    my $inSampleDirectory = $indirpath{$sampleID};
    my $coreCounter=1;
    my $uncompressedFileCounter = 0; #Used to print wait at the right times since infiles cannot be used (can be a mixture of .gz and .fast files)
    for (my $infileCounter=0;$infileCounter<scalar( @{ $infiles{$sampleID} });$infileCounter++) {

	if ($infiles{$sampleID}[$infileCounter] =~/.fastq$/) { #For files ending with .fastq required since there can be a mixture (also .fastq.gz) within the sample dir
	    if ($uncompressedFileCounter eq $coreCounter*$scriptParameters{'maximumCores'}) { #Using only $scriptParameters{'maximumCores'} cores
		
		print GZFASTQ "wait", "\n\n";
		$coreCounter=$coreCounter+1;
	    }
	    my $tempInfile = $infiles{$sampleID}[$infileCounter];
	    print GZFASTQ "gzip ";
	    print GZFASTQ $inSampleDirectory."/".$tempInfile," &", "\n\n"; #InFile
	    $uncompressedFileCounter++;
	    $infiles{$sampleID}[$infileCounter] =~ s/.fastq/.fastq.gz/g; #Replace the .fastq ending with .fastq.gz since this will execute before fastQC screen and mosaikBuild, hence changing the original file name ending from ".fastq" to ".fastq.gz". 

	}
    }
    print GZFASTQ "wait", "\n\n";
    if ($scriptParameters{'pGZip'} == 1) { 
	SampleIDSubmitJob($sampleID,$filename, 1);
    }
    return;
}

sub ReadPedigreeFile {
#Reads famid_pedigree.txt file
#IDN\tSampleID\tMother\tFather\t\Child..n
#$_[0] = filename
#$_[1] = user supplied @sampleIDs info (1 = yes, 0 = no)

    my $filename = $_[0];    
    my $user_sampleid_switch = $_[1];

    open(PEDF, "<$filename") or die "Can't open $filename:$!, \n";    
     
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
	    my @line_info = split("\t",$_);	    #Loads pedigree info
	    if ( $line_info[0] =~ /(\d+)-(\d+|-\d+)-(\d+)(A|U)/) { #Match IDN
		my $local_familyID = $1;
		if ($user_sampleid_switch == 0) {
		    push(@sampleIDs, $line_info[0]); #Save sampleid info
		} 
		if ($3 % 2 == 1) { #Male
#% modulous operator, it gives you an integer representation of the remainder of a division operation. If X / Y divides evenly (that is to say, there's no remainder (or modulus)), the result of X % Y will be zero.
		    $sampleInfo{$local_familyID}{$line_info[0]}{'Sex'} = "M"; #Sex, M=Male
		}
		else { #Female
		   $sampleInfo{$local_familyID}{$line_info[0]}{'Sex'} = "F"; #Sex, F=Female
		}
		if ($4 eq "A") { #Affected
		    $sampleInfo{$local_familyID}{$line_info[0]}{'Disease_status'} = 1; #1=Affected
		}
		else { #Unaffected
		    $sampleInfo{$local_familyID}{$line_info[0]}{'Disease_status'} = 0; #0=Unaffected
		}
		
		if ($line_info[14]) { #Capture kit
		    my @capture_kits = split(";", $line_info[14]);
		    my $capture_kit =  pop(@capture_kits); #Use only the last capture kit since it should be the most interesting
		    
		    for my $supportedCaptureKit (keys %supportedCaptureKits) {
			if ($supportedCaptureKit eq $capture_kit) {
			    if ($exomeTargetBed eq 0) { #No user supplied info on capture kit target BED-file. Add from pedigree file
				$sampleInfo{$local_familyID}{$line_info[0]}{'exomeTargetBed'} = $supportedCaptureKits{$supportedCaptureKit}; #capture kit Bed-file
			    }
			    if ($exomeTargetBedInfileList eq 0) { #No user supplied info on capture kit target BED-file infile list. Add from pedigree file
				$sampleInfo{$local_familyID}{$line_info[0]}{'exomeTargetBedInfileList'} = $supportedCaptureKits{$supportedCaptureKit}.".infile_list"; #capture kit target infile_list
			    }
			    if ($exomeTargetPaddedBedInfileList eq 0) { #No user supplied info on capture kit target BED-file infile list. Add from pedigree file
				$sampleInfo{$local_familyID}{$line_info[0]}{'exomeTargetPaddedBedInfileList'} = $supportedCaptureKits{$supportedCaptureKit}.".pad100.infile_list"; #capture kit padded target infile_list
			    }
			}
		    }
		}	
	    }
	    #push(@{$sampleInfo{$local_familyID}{$line_info[0]}},@line_info[2..4]); #Populate hash of array for Mother/Father/Child. Full hash: hash{FDN}{IDN}[Sex,Affected,Mother/Father/Child]
	}
    } 	
    if ($user_sampleid_switch == 0) {
	@sampleIDs = sort(@sampleIDs); #Lexiographical sort to determine the correct order of ids indata
    }
    print STDOUT "Read pedigree file: $_[0]", "\n\n";
    close(PEDF);
    return;
}

sub ParallelSampleIDSubmitJob {
#Submits parallel jobs within sampleID and infile and includes any previous jobIDs within that sampleID. Use from time of parallelization until time of merge within sampleID. When it is time to merge: the third argument should be all and the subroutine then will only enter the first coding block and hence wait for all jobIDs that previously have been submitted within that sampleID. 
#$_[0] = sampleid
#$_[1] = whole path filename (.sh)
#$_[2] = filename (.sh) or all 
   
    my $parallelsamplejobids=""; #Create string with all previous jobIDs for a sampleID
    my $ret;

    if ( $_[2] eq "all" ) {

	if ( $allsampleIDjobID{$_[0]} ) { #All jobIDs for a sampleID
	    
	    for (my $alljob=0;$alljob<scalar( @{ $allsampleIDjobID{$_[0]} });$alljob++) {	
		
		if ($alljob eq ( scalar( @{ $allsampleIDjobID{$_[0]} }) -1) ) {
		    $parallelsamplejobids .= ":$allsampleIDjobID{$_[0]}[$alljob]"; #last jobID finish without :
		}
		else {
		    $parallelsamplejobids .= ":$allsampleIDjobID{$_[0]}[$alljob]";
		}
	    }
	    $ret = `sbatch --dependency=afterok$parallelsamplejobids $_[1]`; #Supply with dependency of previous within sampleID
	    #$ret = `sbatch --dependency=afterok:$samplejobids $_[1]`; #Supply with dependency of previous within sampleID
	    my ($jobID) = ($ret =~ /Submitted batch job (\d+)/);
	    push ( @{ $paralleljobID{$_[2]} }, $jobID); #Add paralleljobID to hash[sampleID]
	    push ( @{$allsampleIDjobID{$_[0]} }, $jobID); #Add allsamplejobID to hash[sampleID] 
	    print STDOUT "Sbatch script submitted, job id: $jobID\n";
	    print STDOUT "To check status of job, please run \'jobinfo -j $jobID\'\n";
	    print STDOUT "To check status of job, please run \'squeue -j $jobID\'\n";
	}
	else {
	# Iniate chain for sampleID   
	$ret = `sbatch $_[1]`;
	my ($jobID) = ($ret =~ /Submitted batch job (\d+)/);
	push ( @{$jobID{$_[0]} }, $jobID); #Add jobID to hash[sampleID]
	push ( @{$allsampleIDjobID{$_[0]} }, $jobID); #Add allsamplejobID to hash[sampleID]
	print STDOUT "Sbatch script submitted, job id: $jobID\n";
	print STDOUT "To check status of job, please run \'jobinfo -j $jobID\'\n";
	print STDOUT "To check status of job, please run \'squeue -j $jobID\'\n";
	}
    }

    elsif ( ( $jobID{$_[0]} ) || ( $paralleljobID{$_[2]} ) ) { #Any previous jobIDs
	
	if ( $jobID{$_[0]} ) {
	    
	    for (my $job=0;$job<scalar( @{ $jobID{$_[0]} });$job++) {	
		
		if ($job eq ( scalar( @{ $jobID{$_[0]} }) -1) ) {
		    $parallelsamplejobids .= ":$jobID{$_[0]}[$job]"; #last jobID finish without :
		}
		else {
		    $parallelsamplejobids .= ":$jobID{$_[0]}[$job]";
		}
	    }
	}
	if ( $_[2] && $paralleljobID{$_[2]} ) {

	    for (my $parjob=0;$parjob<scalar( @{ $paralleljobID{$_[2]} });$parjob++) {	
		
		if ($parjob eq ( scalar( @{ $paralleljobID{$_[2]} }) -1) ) {
		$parallelsamplejobids .= ":$paralleljobID{$_[2]}[$parjob]"; #last paralleljobID finish without :
		}
		else {
		    $parallelsamplejobids .= ":$paralleljobID{$_[2]}[$parjob]";
		}
	    }
	}
	$ret = `sbatch --dependency=afterok$parallelsamplejobids $_[1]`; #Supply with dependency of previous within sampleID
	#$ret = `sbatch --dependency=afterok:$samplejobids $_[1]`; #Supply with dependency of previous within sampleID
	my ($jobID) = ($ret =~ /Submitted batch job (\d+)/);
	push ( @{ $paralleljobID{$_[2]} }, $jobID); #Add paralleljobID to hash[sampleID] 
	push ( @{$allsampleIDjobID{$_[0]} }, $jobID); #Add allsamplejobID to hash[sampleID]
	print STDOUT "Sbatch script submitted, job id: $jobID\n";
	print STDOUT "To check status of job, please run \'jobinfo -j $jobID\'\n";
	print STDOUT "To check status of job, please run \'squeue -j $jobID\'\n";
    }
    else {
# Iniate chain for paralleljobs within sampleID   
	$ret = `sbatch $_[1]`;
	my ($jobID) = ($ret =~ /Submitted batch job (\d+)/);
	push ( @{$paralleljobID{$_[2]} }, $jobID); #Add paralleljobID to hash[sampleID]
	push ( @{$allsampleIDjobID{$_[0]} }, $jobID); #Add allsamplejobID to hash[sampleID]
	print STDOUT "Sbatch script submitted, job id: $jobID\n";
	print STDOUT "To check status of job, please run \'jobinfo -j $jobID\'\n";
	print STDOUT "To check status of job, please run \'squeue -j $jobID\'\n";
    }
}

sub SampleIDSubmitJob {
#Submits jobs for per sampleID. Use until there is need to parallelize within each sampleID
#$_[0] = sampleid
#$_[1] = filename (.sh)
#$_[2] = add to chain or not (0 = do not add, 1 = add)
    
    my $samplejobids=""; #Create string with all previous jobIDs
    my $ret;
    if ( $jobID{$_[0]} ) {
	
	for (my $job=0;$job<scalar( @{ $jobID{$_[0]} });$job++) {	
	    
	    if ($job eq ( scalar( @{ $jobID{$_[0]} }) -1) ) {
		$samplejobids .= "$jobID{$_[0]}[$job]"; #last jobID finish without :
	    }
	    else {
		$samplejobids .= "$jobID{$_[0]}[$job]:";
	    }
	}
	$ret = `sbatch --dependency=afterok:$samplejobids $_[1]`; #Supply with dependency of previous within sampleID
	my ($jobID) = ($ret =~ /Submitted batch job (\d+)/);
	if ($_[2] == 1) {
	    push ( @{ $jobID{$_[0]} }, $jobID); #Add jobID to hash[sampleID] 
	    push ( @{$allsampleIDjobID{$_[0]} }, $jobID); #Add allsamplejobID to hash[sampleID]
	}
	print STDOUT "Sbatch script submitted, job id: $jobID\n";
	print STDOUT "To check status of job, please run \'jobinfo -j $jobID\'\n";
	print STDOUT "To check status of job, please run \'squeue -j $jobID\'\n";
    }
    else {
# Iniate chain for sampleID   
	$ret = `sbatch $_[1]`;
	my ($jobID) = ($ret =~ /Submitted batch job (\d+)/);
	if ($_[2] == 1) {
	    push ( @{$jobID{$_[0]} }, $jobID); #Add jobID to hash[sampleID]
	    push ( @{$allsampleIDjobID{$_[0]} }, $jobID); #Add allsamplejobID to hash[sampleID]
	}
	print STDOUT "Sbatch script submitted, job id: $jobID\n";
	print STDOUT "To check status of job, please run \'jobinfo -j $jobID\'\n";
	print STDOUT "To check status of job, please run \'squeue -j $jobID\'\n";
    }
}

sub InfilesReFormat {
    
#Code needed to reformat files for mosaik output, which have not yet been created into, correct format so that a sbatch script can be generated with the correct filenames. 

my $uncompressedFileCounter = 0;     

    for my $samplid ( keys %infiles ) { #For every sample id
	
	my $k=1;
	my $itrack=0; #Needed to be able to track when lanes are finished
	for (my $i=0;$i<scalar( @ { $infiles{ $samplid } } );$i++) { #Collects inputfiles for every fastq dir and remakes format
	    if ( $infiles{$samplid}[$i] =~ /\/?([^\.\/]+\.[^\.]+)\.lane(\d+)_([12FfRr])\.fastq.gz/ ) { #Parse fastq.gz 'old' format
		
		push( @ {$lanes{$samplid} }, $2);
		$InfilesLaneNoEnding{ $samplid }[$itrack]= "$1.$2"; #Save new format in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry into hash for both strands and .ending is removed (.fastq).
		$i++; #Skip second direction
		$itrack++; #Track for every lane finished
	    }
	    elsif ( $infiles{$samplid}[$i] =~ /\/?([^\.\/]+\.[^\.]+)\.lane(\d+)_([12FfRr])\.fastq/ ) { #Parse 'old' format

		push( @ {$lanes{$samplid} }, $2);
		$uncompressedFileCounter = 1;
		$InfilesLaneNoEnding{ $samplid }[$itrack]= "$1.$2"; #Save new format in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry into hash for both strands and .ending is removed (.fastq).
		$i++; #Skip second direction
		$itrack++; #Track for every lane finished
	    }
	    elsif ( $infiles{$samplid}[$i] =~ /(\d+)_(\d+)_([^_]+)_([^_]+)_(index[^_]+)_(\d).fastq.gz/ ) { #Parse fastq.gz 'new' format $1=lane, $2=date, $3=Flow-cell, $4=SampleID, $5=index,$6=direction
	
		push( @ {$lanes{$samplid} }, $1);
		$InfilesLaneNoEnding{ $samplid }[$itrack]= "$4.$2_$3_$5."."lane"."$1_$6"; #Save new format (sampleID_date_flow-cell_index_lane_direction) in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry into hash for both strands and .ending is removed (.fastq).
		$i++; #Skip second direction
		$itrack++; #Track for every lane finished
	    }
	    elsif ( $infiles{$samplid}[$i] =~/(\d+)_(\d+)_([^_]+)_([^_]+)_(index[^_]+)_(\d).fastq/ ) { #Parse 'new' format $1=lane, $2=date, $3=Flow-cell, $4=SampleID, $5=index,$6=direction
		
		push( @ {$lanes{$samplid} }, $1);
		$uncompressedFileCounter = 1;
		$InfilesLaneNoEnding{ $samplid }[$itrack]= "$4.$2_$3_$5."."lane"."$1_$6"; #Save new format (sampleID_date_flow-cell_index_lane_direction) in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry into hash for both strands and .ending is removed (.fastq).
		$i++; #Skip second direction
		$itrack++; #Track for every lane finished
	    }
	}
	$k=1;
	for (my $i=0;$i<scalar( @ { $infiles{ $samplid } } );$i++) { #Collects inputfiles for every fastq dir and remakes format
	    if ( $infiles{$samplid}[$i] =~ /\/?([^\.\/]+\.[^\.]+)\.lane(\d+_[12FfRr])\.fastq/ ) { #Parse 'old' format
		
		$InfilesBothStrandsNoEnding{ $samplid }[$i]= "$1.$2"; #Save new format in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry per strand and .ending is removed (.fastq).
		$itrack++; #Track for every lane finished
	    }
	    elsif ( $infiles{$samplid}[$i] =~ /(\d+)_(\d+)_([^_]+)_([^_]+)_(index[^_]+)_(\d).fastq/ ) { #Parse 'new' format
	
		$InfilesBothStrandsNoEnding{ $samplid }[$i]= "$4.$2_$3_$5."."lane"."$1_$6"; #Save new format in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry per strand and .ending is removed (.fastq).
		$itrack++; #Track for every lane finished
	    }
			    
	}
    }
return $uncompressedFileCounter;
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

sub WriteCMDMasterLogg {
    
    open (MASTERL, ">>$masterLoggName") or die "Can't write to $masterLoggName: $!\n"; #Open file run logg
    
    print MASTERL 
	"-p ", $scriptParameters{'projectID'};
    if ( $scriptParameters{'email'} ) {
	print MASTERL " -email ", $scriptParameters{'email'};
    }
    print MASTERL " -familyID ", $scriptParameters{'familyID'};
    if ( $scriptParameters{'pedigreeFile'} ) {
	print MASTERL " -pedigreeFile ", $scriptParameters{'pedigreeFile'};
    }
    print MASTERL
	" -sampleIDs ", $scriptParameters{'sampleIDs'},
	" -inFilesDir ", $scriptParameters{'inFilesDir'},
	" -outDataDir ", $scriptParameters{'outDataDir'},   
	" -inScriptDir ", $scriptParameters{'inScriptDir'},
	" -outScriptDir ", $scriptParameters{'outScriptDir'},
	" -referencesDir ", $scriptParameters{'referencesDir'};
    if ( $scriptParameters{'humanGenomeReference'} ) {
	print MASTERL " -humanGenomeReference ", $scriptParameters{'humanGenomeReference'};
    }
    print MASTERL " -pFastQC ", $scriptParameters{'pFastQC'};
    if ($scriptParameters{'aligner'}) {
	print MASTERL " -aligner ", $scriptParameters{'aligner'};
    }
    print MASTERL
	" -pMosaikBuild ", $scriptParameters{'pMosaikBuild'};
    if ( $scriptParameters{'pMosaikBuild'} == 1  ) {
	print MASTERL " -mosaikBuildMedianFragLength ", $scriptParameters{'mosaikBuildMedianFragLength'};
    }
    print MASTERL
	" -pMosaikAlign ", $scriptParameters{'pMosaikAlign'};
    if ( $scriptParameters{'pMosaikAlign'} == 1  ) {
	print MASTERL 
	    " -mosaikAlignReference ", $scriptParameters{'mosaikAlignReference'},
	    " -mosaikAlignNeuralNetworkPeFile ", $scriptParameters{'mosaikAlignNeuralNetworkPeFile'},
	    " -mosaikAlignNeuralNetworkSeFile ", $scriptParameters{'mosaikAlignNeuralNetworkSeFile'},
	    " -mosaikJumpDbStub ", $scriptParameters{'mosaikJumpDbStub'};
    }
    print MASTERL " -pBwaAln ", $scriptParameters{'pBwaAln'};
    if ( $scriptParameters{'pBwaAln'} == 1  ) {
	print MASTERL " -bwaAlnQualityTrimming ", $scriptParameters{'bwaAlnQualityTrimming'};
    }
    print MASTERL 
	" -pBwaSampe ", $scriptParameters{'pBwaSampe'},
	" -pSamToolsSort ", $scriptParameters{'pSamToolsSort'},
	" -pPicardToolsMerge ", $scriptParameters{'pPicardToolsMerge'};
    if ( scalar(@picardToolsMergePreviousFiles) ) {
	print MASTERL
	    " -picardToolsMergePreviousFiles ", $scriptParameters{'picardToolsMergePreviousFiles'};
    }
    print MASTERL
	" -pPicardToolsMarkduplicates ", $scriptParameters{'pPicardToolsMarkduplicates'};
    if ( $scriptParameters{'picardToolsPath'} ) {
	print MASTERL
	    " -picardToolsPath ", $scriptParameters{'picardToolsPath'};
    }
    print MASTERL
	" -pCalculateCoverage ", $scriptParameters{'pCalculateCoverage'};
    if ( $scriptParameters{'pCalculateCoverage'} ==1) {
	print MASTERL
	    " -pGenomeCoverageBED ", $scriptParameters{'pGenomeCoverageBED'},
	    " -pCoverageBED ", $scriptParameters{'pCoverageBED'},
	    " -pQaCompute ", $scriptParameters{'pQaCompute'},
	    " -pPicardToolsCollectMultipleMetrics ", $scriptParameters{'pPicardToolsCollectMultipleMetrics'},
	    " -pPicardToolsCalculateHSMetrics ", $scriptParameters{'pPicardToolsCalculateHSMetrics'};
	if ( ($scriptParameters{'pGenomeCoverageBED'} ==1) || ($scriptParameters{'pQaCompute'} ==1) ) {
	    print MASTERL
		" -xCoverage ", $scriptParameters{'xCoverage'};
	}    
    }
    if ( $identicalCaptureBedCounter eq scalar(@sampleIDs) ) { #Same capture kit for all sampleIDs
	print MASTERL 
	    " -exomeTargetBed ", $scriptParameters{$sampleIDs[0]}{'exomeTargetBed'};
    }
    elsif ($exomeTargetBed) {
	print MASTERL 
	    " -exomeTargetBed ", $scriptParameters{$sampleIDs[0]}{'exomeTargetBed'};
    }
    if ( $identicalCaptureBedIntervalCounter eq scalar(@sampleIDs) ) { #Same capture kit for all sampleIDs
	print MASTERL 
	    " -exomeTargetBedInfileList ", $scriptParameters{$sampleIDs[0]}{'exomeTargetBedInfileList'},
	    " -exomeTargetPaddedBedInfileList ", $scriptParameters{$sampleIDs[0]}{'exomeTargetPaddedBedInfileList'};
    }
    elsif ($exomeTargetBedInfileList) {
	print MASTERL
	    " -exomeTargetBedInfileList ", $scriptParameters{$sampleIDs[0]}{'exomeTargetBedInfileList'},
	    " -exomeTargetPaddedBedInfileList ", $scriptParameters{$sampleIDs[0]}{'exomeTargetPaddedBedInfileList'};
    }
    print MASTERL
	" -pRCovPlots ", $scriptParameters{'pRCovPlots'},
	" -pGZip ", $scriptParameters{'pGZip'},
	" -pRemovalRedundantFiles ", $scriptParameters{'pRemovalRedundantFiles'},
	" -wholeGenomeSequencing ", $scriptParameters{'wholeGenomeSequencing'},
	" -maximumCores ",$scriptParameters{'maximumCores'};
    if ( $scriptParameters{'configFile'} ne 0) {
	print MASTERL " -configFile ", $scriptParameters{'configFile'};
    }
    if ( $scriptParameters{'environmentUppmax'} == 1) {
	print MASTERL " -environmentUppmax ", $scriptParameters{'environmentUppmax'};
    }
    print MASTERL
	" -writeConfigFile ", $scriptParameters{'writeConfigFile'}, "\n";
    #Note FileHandle MASTERL not closed
    return;
}
