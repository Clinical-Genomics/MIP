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

-pFQC/--pFastQC Flag for running FastQC (defaults to "1" (=yes))

-pMoB/--pMosaikBuild Flag running MosaikBuild (defaults to "1" (=yes))

-mobmfl/--mosaikBuildMeanFragLength Flag for setting the mean fragment length, mfl, (defaults to (=375) bp)

-pMoA/--pMosaikAlign Flag running MosaikAlign (defaults to "1" (=yes))

-moaref/--mosaikAlignReference MosaikAlign reference (defaults to "")

-moaannpe/--mosaikAlignNeuralNetworkPeFile MosaikAlign Neural Network PE File (defaults to "")

-moaannse/--mosaikAlignNeuralNetworkSeFile MosaikAlign Neural Network SE File (defaults to "")

-mojdb/--mosaikJumpDbStub MosaikJump stub (defaults to "")

-pBWA_aln/--pBwaAln Flag running bwa aln (defaults to "0" (=no))

-bwaalnq/--bwaAlnQualityTrimming BWA Aln quality threshold for read trimming (defaults to "20")

-pBWA_sampe/--pBwaSampe Flag running bwa sampe (defaults to "0" (=no))

-pSamT_sort/--pSamToolsSort Flag running samtools sort & index (defaults to "1" (=yes))

-pPicT_merge/--pPicardToolsMerge Flag running picardtools MergeSamFiles (defaults to "0" (=no))

-picT_mergeprev/--picardToolsMergePreviousFiles Flag running picardTools MergeSamFiles on merged current files and previous BAM-file(s) (Supply whole path and name, name must contain sample id, and lanes_Xn info)

-pPicT_MarkDup/--pPicardToolsMarkduplicates Flag running PicardTools MarkDuplicates (defaults to "1" (=yes))

-pic_path/--picardPath  Flag for path to picardtools, must be supplied for picardtools (defaults to "")

-pCC/--pCalculateCoverage Flag running coverage tools: qaCompute, genomeCoverageBED and PicardTools (defaults to "1" (=yes))

-pCC_Bedgc/--pGenomeCoverageBED Flag running genomeCoverageBED under pCC (defaults to "1" (=yes))

-pCC_Bedc/--pCoverageBED Flag running coverageBED under pCC (defaults to "1" (=yes))

-extb/--exomeTargetBed Target BED file of exome capture for coverageBed '-pCC_Bedc'. (defaults to "")

-pCC_Qac/--pQaCompute Flag running qaCompute under pCC (defaults to "1" (=yes))

-xcov/--xCoverage  Flag determining coverage depth genomeCoverageBED, qaCompute (defaults to "30")

-pCC_PicMM/--pPicardToolsCollectMultipleMetrics Flag running picardTools collectMultipleMetrics under pCC (defaults to "1" (=yes))

-pCCE_PicHS/--pPicardToolsCalculateHSMetrics Flag running picardTools calulateHSmetrics under pCC (defaults to "1" (=yes))

-extbl/--exomeTargetBedInfileList Prepared target BED file for picardTools calculateHSMetrics. (defaults to "". File ending should be ".infile_list")
              
-extpbl/--exomeTargetPaddedBedInfileList Prepared padded target BED file for picardTools calculateHSMetrics. (defaults to "". File should be ".padXXX.infile_list")

-pRCP/--pRCovPlots  Flag running rCovPlots (defaults to "1" (=yes))

-pGZ/--pGZip  Flag running gzip for fastq files (defaults to "1" (=yes))

-pREM/--pRemovalRedundantFiles  Flag generating sbatch script of deletion of redundant files (defaults to "1" (=yes);Note: Must be submitted manually)

-al/--aligner  Flag for determining which aligner was used previously if none was specified (defaults to "")

-wgs/--wholeGenomeSequencing Analysis to perform are whole genome sequencing data or not (defaults to "0" (=no))

-mc/--maximumCores The maximum number of cores per node used in the analysis (defaults to "8")

-env_up/--environmentUppmax Sets the environment to UPPMAX. (defaults to "0" (=no))

-c/--configFile YAML config file. (defaults to "")

-wc/--writeConfigFile Write YAML config file. (defaults to "";Supply whole path)

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
	       -pFQC/--pFastQC Flag running FastQC (defaults to "1" (=1))
	       -pMoB/--pMosaikBuild Flag running MosaikBuild (defaults to "1" (=yes))
               -mobmfl/--mosaikBuildMeanFragLength Flag for setting the mean fragment length, mfl, (defaults to (=375) bp)
	       -pMoA/--pMosaikAlign Flag running MosaikAlign (defaults to "1" (=yes))
               -moaref/--mosaikAlignReference MosaikAlign reference (defaults to "")
               -moaannpe/--mosaikAlignNeuralNetworkPeFile MosaikAlign Neural Network PE File (defaults to "")
               -moaannse/--mosaikAlignNeuralNetworkSeFile MosaikAlign Neural Network SE File (defaults to "")
               -mojdb/--mosaikJumpDbStub MosaikJump stub (defaults to "")
               -pBWA_aln/--pBwaAln Flag running bwa aln (defaults to "0" (=no))
               -bwaalnq/--bwaAlnQualityTrimming BWA Aln quality threshold for read trimming (defaults to "20")
               -pBWA_sampe/--pBwaSampe Flag running bwa sampe (defaults to "0" (=no))
               -pSamT_sort/--pSamToolsSort Flag running samtools sort & index (defaults to "1" (=yes))
               -pPicT_merge/--pPicardToolsMerge Flag running picardtools MergeSamFiles (defaults to "0" (=no))
               -picT_mergeprev/--picardToolsMergePreviousFiles Flag running picardTools MergeSamFiles on merged current files and previous BAM-file(s) (Supply whole path and name, name must contain sample id, and lanes_Xn info)
               -pPicT_MarkDup/--pPicardToolsMarkduplicates Flag running PicardTools MarkDuplicates (defaults to "1" (=yes))
               -pic_path/--picardPath  Flag for path to picardtools, must be supplied for picardtools (defaults to "")
               -pCC/--pCalculateCoverage Flag running coverage tools: qaCompute, genomeCoverageBED and PicardTools (defaults to "1" (=yes))
               -pCC_Bedgc/--pGenomeCoverageBED Flag running genomeCoverageBED under pCC (defaults to "1" (=yes))
               -pCC_Bedc/--pCoverageBED Flag running coverageBED under pCC (defaults to "1" (=yes))
               -extb/--exomeTargetBed Target BED file of exome capture for coverageBed '-pCC_Bedc'. (defaults to "")
               -pCC_Qac/--pQaCompute Flag running qaCompute under pCC (defaults to "1" (=yes))
               -xcov/--xCoverage  Flag determining coverage depth genomeCoverageBED, qaCompute (defaults to "30")
               -pCC_PicMM/--pPicardToolsCollectMultipleMetrics Flag running picardTools collectMultipleMetrics under pCC (defaults to "1" (=yes))
               -pCCE_PicHS/--pPicardToolsCalculateHSMetrics Flag running picardTools calulateHSmetrics under pCC (defaults to "1" (=yes))
               -extbl/--exomeTargetBedInfileList Prepared target BED file for picardTools calculateHSMetrics. (defaults to "". File ending should be ".infile_list") 
               -extpbl/--exomeTargetPaddedBedInfileList Prepared padded target BED file for picardTools calculateHSMetrics. (defaults to "". File ending should be ".padXXX.infile_list")
	       -pRCP/--pRCovPlots Flag running rCovPlots (defaults to "1" (=yes))
	       -pGZ/--pGZip  Flag generating gzip sbatch (defaults to "1" (=yes))
               -pREM/--pRemovalRedundantFiles  Flag generating sbatch script of deletion of redundant files (defaults to "1" (=yes);Note: Must be submitted manually)
               -al/--aligner  Flag for determining which aligner was used previously if none was specified (defaults to "")
               -wgs/--wholeGenomeSequencing Analysis to perform are whole genome sequencing data or not (defaults to "0" (=no))
               -mc/--maximumCores The maximum number of cores per node used in the analysis (defaults to "8")
               -env_up/--environmentUppmax Sets the environment to UPPMAX. (defaults to "0" (=no))
               -c/--configFile YAML configuration file. (defaults to "")
               -wc/--writeConfigFile Write YAML configuration file. (defaults to "";Supply whole path)
	   };
}

###
#Program parameters
###

#Project specific
my ($projectID,$em, $inScriptDir, $referencesDir, $outDataDir, $outScriptDir, $familyID, $pedigreeFile, $configFile, $writeConfigFile) = (0,0,0,0,0,0,0,0,0,0);
my (@inFilesDir,@sampleIDs); #Arrays for input file directorys,sampleIDs

#GZip
my ($pGZ) = (-1);
#FastQC
my ($pFQC) = (-1);

#Mosaik
my ($pMoB, $pMoA) = (-1,-1);
my ($mosaikBuildMeanFragLength, $mosaikAlignReference, $mosaikAlignNeuralNetworkPeFile, $mosaikAlignNeuralNetworkSeFile, $mosaikJumpDbStub) = (375,0,0,0,0);

#BWA
my ($pBWA_aln, $pBWA_sampe) = (-1,-1);
my ($bwaAlnQualityTrimming) = (-1);

#SamTools
my ($pSamT_sort) = (-1);

#PicardTools
my ($pPicT_merge, $pPicT_MarkDup) = (-1,-1);
my (@picardToolsMergePreviousFiles);

#Coverage
my ($pCC, $pCC_Bedgc, $pCC_Bedc, $pCC_Qac, $pCC_PicMM, $pCCE_PicHS, $pGdb, $pRCP) = (-1,-1,-1,-1,-1,-1,-1,-1);
my ($exomeTargetBed, $exomeTargetBedInfileList, $exomeTargetPaddedBedInfileList, $picardPath, $xcoverage, $identicalCaptureBedCounter, $identicalCaptureBedIntervalCounter, $identicalCaptureBedPaddedIntervalCounter) = (0,0,0,0,30,0,0,0);

#Pipe
my ($pREM) = (-1);
my ($humanGenomeReference, $humanGenomeReferenceSource, $humanGenomeRefereceChromosomePrefix, $humanGenomeReferenceVersion, $fnend, $aligner, $wholeGenomeSequencing, $maximumCores, $environmentUppmax, $filename, $fnt, $fnt2, $help) = (0, 0, 0, 0, ".sh", 0, 0,8,0); #Arguments for project
my (@chromosomes);
my (%infiles, %indirpath, %Infiles_lane_noending, %lanes, %Infiles_bothstrands_noending, %jobID, %paralleljobID, %allsampleIDjobID, %sample_info, %script_parameters); 
#%infiles=from platform (Illumina), %indirpath for the path to infiles, %Infiles_lane_noending for MosaikBuild (one entry for both strands), %lanes for sample lanes, Infiles_bothstrands_noending for bwa_aln (one entry per strand)

###
#Staging Area
###

#Capture kits supported from pedigree file.
my %supported_capture_kits = (
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
	   'em|email:s'  => \$em,
	   'odd|outDataDir:s'  => \$outDataDir, #One dir above sample id, must supply whole path i.e. /proj/...
	   'osd|outScriptDir:s'  => \$outScriptDir,  #One dir above sample id, must supply whole path i.e. /proj/...
	   'f|familyID:s' => \$familyID, #Family group ID (Merged to same vcf file after GATK Base Recalibration)
	   'pedigree|pedigreeFile:s' => \$pedigreeFile, #Pedigree file
	   'huref|humanGenomeReference:s' => \$humanGenomeReference, #Human genome reference
	   'pFQC|pFastQC:n' => \$pFQC,
	   'pMoB|pMosaikBuild:n' => \$pMoB,
	   'mobmfl|mosaikBuildMeanFragLength:n' => \$mosaikBuildMeanFragLength, #for fragment length estimation and local search
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
	   'pic_path|picardPath:s' => \$picardPath, #Path to picardtools
	   'pCC|pCalculateCoverage:n' => \$pCC,
	   'pCC_Bedgc|pGenomeCoverageBED:n' => \$pCC_Bedgc,
	   'pCC_Bedc|pCoverageBED:n' => \$pCC_Bedc,
	   'extb|exomeTargetBed:s' => \$exomeTargetBed, #target file for coverageBed
	   'pCC_Qac|pQaCompute:n' => \$pCC_Qac,
	   'xcov|xCoverage:n' => \$xcoverage, #Sets max depth to calculate coverage
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
    
    if ($environmentUppmax == 1) {
	use lib '/bubo/home/h11/henriks/lib/'; #YAML not installed at @UPPMAX
    }
    else {
	use lib '/bubo/home/h12/henriks/lib/'; #Remove later 
    }
    use YAML;
    open (YAML, '<', $configFile) or die "can't open $configFile: $!";
    %script_parameters = YAML::LoadFile($configFile);
    close(YAML);
}

if ($projectID eq 0) { #Command line has precedence.   
    unless ( defined($script_parameters{'projectID'}) ) { #No input from config file   
	print STDERR "\n";
	print STDERR "Must supply a project ID. Use '-projectID'", "\n\n";
	die $USAGE;
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $script_parameters{'projectID'} = $projectID;   
}

if ( $familyID eq 0 ) {
     if ( defined($script_parameters{'familyID'}) ) { #Input from config file - Do nothing
     }
     else {
	 print STDERR "\n";
	 print STDERR "Must supply a family id. Use '-familyID'. -If not applicable supply the same familyID as the sampleid ", "\n\n";
	 die $USAGE;
     }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $script_parameters{'familyID'} = $familyID;   
}

if ($outDataDir eq 0) {
    
    if ( defined($script_parameters{'outDataDir'}) ) { #Input from config file - do nothing
	
    }
    elsif ($environmentUppmax == 1) { #Use Uppmax default
	print STDOUT "\n";
	if ($wholeGenomeSequencing == 1) {
	    $outDataDir = "/proj/$script_parameters{'projectID'}/private/genomes";
	}
	else {
	    $outDataDir = "/proj/$script_parameters{'projectID'}/private/exomes";
	}
	print STDOUT "Setting MIP output data dir to: $outDataDir", "\n\n";
	$script_parameters{'outDataDir'} = $outDataDir;
    }
    else {
	print STDERR "\n";
	print STDERR "Must supply a MIP output data dir", "\n\n";
	die $USAGE;
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $script_parameters{'outDataDir'} = $outDataDir;   
}

if ($pedigreeFile eq 0) {
    
    if ( defined($script_parameters{'pedigreeFile'}) ) { #Input from config file - ReadPedigreefile
	if ( scalar(@sampleIDs) == 0 ) { #No user supplied sample info
	    if ( defined($script_parameters{'sampleIDs'}) ) { #sampleIDs info in config file
		ReadPedigreeFile($script_parameters{'pedigreeFile'}, 1);  # scalar(@sampleIDs) = 0:No user supplied sample info, but present in config file do NOT overwrite using info from pedigree file
	    }
	    else { #No sampleIDs info in config file
		ReadPedigreeFile($script_parameters{'pedigreeFile'}, scalar(@sampleIDs) );  # scalar(@sampleIDs) = 0:No user supplied sample info, not defined $script_parameters{'sampleIDs'} in config file, add it from pedigree file
	    }
	}
	else {
	    ReadPedigreeFile($script_parameters{'pedigreeFile'}, scalar(@sampleIDs) );  # User supplied sample info, do NOT overwrite using info from pedigree file
	}
    }
    elsif ($environmentUppmax == 1) {
	print STDOUT "\n";
	$pedigreeFile = "$script_parameters{'outDataDir'}/$script_parameters{'familyID'}/$script_parameters{'familyID'}"."_pedigree.txt";
	print STDOUT "Assuming location of pedigree file to be: $pedigreeFile", "\n\n";
	if (-e $pedigreeFile) { #if file exists 
	    print STDOUT "Found pedigree file at: $pedigreeFile", "\n\n";
	    $script_parameters{'pedigreeFile'} = $pedigreeFile; #Add to enable recreation of cmd line later
	    ReadPedigreeFile($pedigreeFile, scalar(@sampleIDs) ); #  scalar(@sampleIDs)= 0:No user supplied sample info, add it from pedigree file
	}
	else { 
	    print STDERR "Could not find pedigree file at: $pedigreeFile \n";
	    die $USAGE;
	} 
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    ReadPedigreeFile($pedigreeFile, scalar(@sampleIDs) );  # scalar(@sampleIDs) = 0:No user supplied sample info, add it from pedigree file
    $script_parameters{'pedigreeFile'} = $pedigreeFile;   
}

if (scalar (@inFilesDir) == 0) {
    if ( defined($script_parameters{'inFilesDir'}) ) { #Input from config file - do nothing
    }
    elsif ($environmentUppmax == 1 && $script_parameters{'pedigreeFile'} ne 0) { #Locate infile directory. Order dictated by @sampleIDs (lexiographical)
	for (my $indir_Count=0;$indir_Count<scalar(@sampleIDs);$indir_Count++) {
	    if ($wholeGenomeSequencing == 0) {
		push(@inFilesDir, "/proj/$script_parameters{'projectID'}/private/exomes/$sampleIDs[$indir_Count]/fastq");
	    }
	    else {
		push(@inFilesDir, "/proj/$script_parameters{'projectID'}/private/genomes/$sampleIDs[$indir_Count]/fastq");
	    }
	}
	$script_parameters{'inFilesDir'} = join(',',@inFilesDir); #Add to enable recreation of cmd line later
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
    $script_parameters{'inFilesDir'} = join(',',@inFilesDir);   
}

if ( scalar(@sampleIDs) == 0) {
    if ( defined($script_parameters{'sampleIDs'}) ) { #Input from config file - do nothing
	@sampleIDs = split(/,/, $script_parameters{'sampleIDs'}); #Transfer to array
    }
    else {
	print STDERR "\n";
	print STDERR "Must supply a sample ID as a comma separated list", "\n\n";
	die $USAGE;
    }
}
else {
    $script_parameters{'sampleIDs'} = join(',',@sampleIDs); #Add to enable recreation of cmd line later
    @sampleIDs = split(/,/,join(',',@sampleIDs)); #Enables comma separated list of sample IDs from user supplied cmd info
}

if ( $inScriptDir eq 0) {
    if ( defined($script_parameters{'inScriptDir'}) ) { #Input from config file - do nothing
	
    }
    elsif ($environmentUppmax == 1) {
	print STDOUT "\n";
	$inScriptDir = "/proj/$script_parameters{'projectID'}/private/mip_scripts_master";
	print STDOUT "Setting the MIP scripts dir to: $inScriptDir", "\n\n";
	$script_parameters{'inScriptDir'} = $inScriptDir; #Add to enable recreation of cmd line later
    }
    else {
	print STDERR "\n";
	print STDERR "Must supply the MIP scripts dir", "\n\n";
	die $USAGE;
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $script_parameters{'inScriptDir'} = $inScriptDir;   
}

if ( $referencesDir eq 0) {
    if ( defined($script_parameters{'referencesDir'}) ) { #Input from config file - do nothing
	
    }
    elsif ($environmentUppmax == 1) {
	print STDOUT "\n";
	$referencesDir = "/proj/$script_parameters{'projectID'}/private/mip_references";
	print STDOUT "Setting MIP reference dir to: $referencesDir", "\n\n";
	$script_parameters{'referencesDir'} = $referencesDir; #Add to enable recreation of cmd line later
    }
    else {
	print STDERR "\n";
	print STDERR "Must supply a MIP reference dir", "\n\n";
	die $USAGE;
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $script_parameters{'referencesDir'} = $referencesDir; #Add to enable recreation of cmd line later
}

if ($outScriptDir eq 0) {
    if ( defined($script_parameters{'outScriptDir'}) ) { #Input from config file - do nothing	
    }
    elsif ($environmentUppmax == 1) {
	print STDOUT "\n";
	if ($wholeGenomeSequencing == 1) {
	    $outScriptDir = "/proj/$script_parameters{'projectID'}/private/genomes_scripts";
	}
	else {
	    $outScriptDir = "/proj/$script_parameters{'projectID'}/private/exomes_scripts";
	}
	print STDOUT "Setting MIP output scripts dir to: $outScriptDir", "\n\n";
	$script_parameters{'outScriptDir'} = $outScriptDir; #Add to enable recreation of cmd line later
    }
    else {
	print STDERR "\n";
	print STDERR "Must supply a MIP output script dir", "\n\n";
	die $USAGE;
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $script_parameters{'outScriptDir'} = $outScriptDir; #Add to enable recreation of cmd line later
}

#
#GZip
#
if ( $pGZ == -1) { #Not set by cmd 
    if ( defined($script_parameters{'pGZip'}) ) { #Input from config file - do nothing	
    }
    else {
	$script_parameters{'pGZip'} = 1; #Default
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $script_parameters{'pGZip'} = $pGZ;
}

#
#FASTQC
#
if ( $pFQC == -1) { #Not set by cmd 
    if ( defined($script_parameters{'pFastQC'}) ) { #Input from config file - do nothing	
    }
    else {
	$script_parameters{'pFastQC'} = 1; #Default
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $script_parameters{'pFastQC'} = $pFQC;
}

#
#MosaikBuild
#
if ( $pMoB == -1) { #Not set by cmd 
    if ( defined($script_parameters{'pMosaikBuild'}) ) { #Input from config file - do nothing	
    }
    else {
	$script_parameters{'pMosaikBuild'} = 1; #Default
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $script_parameters{'pMosaikBuild'} = $pMoB;
}

#
#MosaikAlign
#
if ( $pMoA == -1) { #Not set by cmd
    if ( defined($script_parameters{'pMosaikAlign'}) ) { #Input from config file - do nothing	
    }
    else {
	$script_parameters{'pMosaikAlign'} = 1; #Default
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $script_parameters{'pMosaikAlign'} = $pMoA;
}
if ( $script_parameters{'pMosaikAlign'} eq 1 ) { #MosaikAlign is to be used - check prerequisets
    if ($mosaikAlignReference eq 0) { #No input from cmd
	if ( defined($script_parameters{'mosaikAlignReference'}) ) { #Input from config file - do nothing	
	}
	elsif ($environmentUppmax == 1) { #Use default 
	    $mosaikAlignReference = "Homo_sapiens.GRCh37.70_nochr.dat";
	    $script_parameters{'mosaikAlignReference'} = $mosaikAlignReference;
	}
	else {
	    print STDERR "\nSupply '-mosaikAlignReference' if you want to run MosaikAlign\n\n";
	    die $USAGE;
	}
    }
    if ($mosaikAlignNeuralNetworkPeFile eq 0) { #No input from cmd
	if ( defined($script_parameters{'mosaikAlignNeuralNetworkPeFile'}) ) { #Input from config file - do nothing	
	}
	elsif ($environmentUppmax == 1) {
	    $mosaikAlignNeuralNetworkPeFile = "2.1.78.pe.ann";
	    $script_parameters{'mosaikAlignNeuralNetworkPeFile'} = $mosaikAlignNeuralNetworkPeFile;
	}
	else {
	    print STDERR "\nSupply '-mosaikAlignNeuralNetworkPeFile' if you want to run MosaikAlign\n\n";
	    die $USAGE;
	}
    }
    if ($mosaikAlignNeuralNetworkSeFile eq 0) { #No input from cmd
	if ( defined($script_parameters{'mosaikAlignNeuralNetworkSeFile'}) ) { #Input from config file - do nothing	
	}
	elsif ($environmentUppmax == 1) {
	    $mosaikAlignNeuralNetworkSeFile = "2.1.78.se.ann";
	    $script_parameters{'mosaikAlignNeuralNetworkSeFile'} = $mosaikAlignNeuralNetworkSeFile;
	}
	else {
	    print STDERR "\nSupply '-mosaikAlignNeuralNetworkSeFile' if you want to run MosaikAlign\n\n";
	    die $USAGE;
	}
    }
    if ($mosaikJumpDbStub eq 0) { #No input from cmd
	if ( defined($script_parameters{'mosaikJumpDbStub'}) ) { #Input from config file - do nothing	
	}
	elsif ($environmentUppmax == 1) {
	    $mosaikJumpDbStub = "Homo_sapiens.GRCh37.70_nochr_jdb_15";
	    $script_parameters{'mosaikJumpDbStub'} = $mosaikJumpDbStub;
	}
	else {
	    print STDERR "\nSupply '-mosaikJumpDbStub' if you want to run MosaikAlign\n\n";
	    die $USAGE;
	}
    }
#File existence checks
    unless (-e "$script_parameters{'referencesDir'}/$script_parameters{'mosaikAlignReference'}" ) { #Check existence of mosaik hashed reference in supplied reference dir
	print STDERR "\nCould not find intended Mosaik Reference file: $script_parameters{'referencesDir'}/$script_parameters{'mosaikAlignReference'}\n\n";
	die $USAGE;		
    }
    unless (-e "$script_parameters{'referencesDir'}/$script_parameters{'mosaikAlignNeuralNetworkPeFile'}" ) { #Check existence of mosaik neural network file in supplied reference dir
	print STDERR "\nCould not find intended Mosaik Neural Network file: $script_parameters{'referencesDir'}/$script_parameters{'mosaikAlignNeuralNetworkPeFile'}\n\n";
	die $USAGE;		
    }
    unless (-e "$script_parameters{'referencesDir'}/$script_parameters{'mosaikAlignNeuralNetworkSeFile'}" ) { #Check existence of mosaik neural network file in supplied reference dir
	print STDERR "\nCould not find intended Mosaik Neural Network file: $script_parameters{'referencesDir'}/$script_parameters{'mosaikAlignNeuralNetworkSeFile'}\n\n";
	die $USAGE;		
    }
    unless (-e "$script_parameters{'referencesDir'}/$script_parameters{'mosaikJumpDbStub'}"."_keys.jmp" ) { #Check existence of Mosaik jumpDb in supplied reference dir
	print STDERR "\nCould not find intended Mosaik JumpDb: $script_parameters{'referencesDir'}/$script_parameters{'mosaikJumpDbStub'}"."_keys.jmp\n\n";
	die $USAGE;		
    }
    unless (-e "$script_parameters{'referencesDir'}/$script_parameters{'mosaikJumpDbStub'}"."_meta.jmp" ) { #Check existence of Mosaik jumpDb in supplied reference dir
	print STDERR "\nCould not find intended Mosaik JumpDb: $script_parameters{'referencesDir'}/$script_parameters{'mosaikJumpDbStub'}"."_meta.jmp\n\n";
	die $USAGE;		
    }
    unless (-e "$script_parameters{'referencesDir'}/$script_parameters{'mosaikJumpDbStub'}"."_positions.jmp" ) { #Check existence of Mosaik jumpDb in supplied reference dir
	print STDERR "\nCould not find intended Mosaik JumpDb: $script_parameters{'referencesDir'}/$script_parameters{'mosaikJumpDbStub'}"."_positions.jmp\n\n";
	die $USAGE;		
    }
}
#
#BWA_aln
#
if ( $pBWA_aln == -1) { #No input from cmd
    if ( defined($script_parameters{'pBwaAln'}) ) { #Input from config file - do nothing	
    }
    else {
	$script_parameters{'pBwaAln'} = 0; #Default
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $script_parameters{'pBwaAln'} = $pBWA_aln;
}
if ( $script_parameters{'pBwaAln'} ==1 ) {
    if ( $bwaAlnQualityTrimming == -1) { #No input from cmd
	if ( defined($script_parameters{'bwaAlnQualityTrimming'}) ) { #Input from config file - do nothing	
	}
	else {
	    $script_parameters{'bwaAlnQualityTrimming'} = 20; #Default
	}
    }
    else {
	$script_parameters{'bwaAlnQualityTrimming'} = $bwaAlnQualityTrimming;
    }
}

#
#BWA_sampe
#
if ( $pBWA_sampe == -1) { #No input from cmd
    if ( defined($script_parameters{'pBwaSampe'}) ) { #Input from config file - do nothing	
    }
    else {
	$script_parameters{'pBwaSampe'} = 0; #Default
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $script_parameters{'pBwaSampe'} = $pBWA_sampe;
}

if ( $aligner eq 0 ) { #No input from cmd
    unless ( defined($script_parameters{'aligner'}) ) { #No input from config file
	if ( ($script_parameters{'pMosaikBuild'} ==1) || ($script_parameters{'pMosaikAlign'} ==1) ) { #Mosaik track
	    if ( ($script_parameters{'pBwaAln'} ==0) && ($script_parameters{'pBwaSampe'} ==0) ) {
		$script_parameters{'aligner'} = "mosaik";
	    }
	    else {
		print STDERR "\n";
		print STDERR "You have to choose either mosaik or bwa to perform alignments or specify which aligner (-aligner 'mosaik' or 'bwa') was used if you want to only run programs after alignment.", "\n\n";
		die $USAGE;
	    }
	}
	if ( ($script_parameters{'pBwaAln'} ==1) || ($script_parameters{'pBwaSampe'} ==1) ) { #BWA track
	    if ( ($script_parameters{'pMosaikBuild'} ==0) && ($script_parameters{'pMosaikAlign'} ==0) ) {
		$script_parameters{'aligner'} = "bwa";
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
    $script_parameters{'aligner'} = $aligner;
}

#
#SamTools Sort/Index
#
if ( $pSamT_sort == -1) { #No input from cmd
    if ( defined($script_parameters{'pSamToolsSort'}) ) { #Input from config file - do nothing	
    }
    else {
	$script_parameters{'pSamToolsSort'} = 1; #Default
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $script_parameters{'pSamToolsSort'} = $pSamT_sort;
}

#
#PicardToolsMerge
#
if ( $pPicT_merge == -1) { #No input from cmd
    if ( defined($script_parameters{'pPicardToolsMerge'}) ) { #Input from config file - do nothing	
    }
    else {
	$script_parameters{'pPicardToolsMerge'} = 0; #Default
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $script_parameters{'pPicardToolsMerge'} = $pPicT_merge;
}
if ( $script_parameters{'pPicardToolsMerge'} ==1) {
    if (scalar(@picardToolsMergePreviousFiles) == 0 ) {
	if ( defined($script_parameters{'picardToolsMergePreviousFiles'}) ) { #Input from config file - transfer to array
	    @picardToolsMergePreviousFiles = split(/,/, $script_parameters{'picardToolsMergePreviousFiles'}); #Transfer to array
	}
    }
    else {
	$script_parameters{'picardToolsMergePreviousFiles'} = join(',',@picardToolsMergePreviousFiles);
    }
}

#
#PicardToolsMarkduplicates
#
if ( $pPicT_MarkDup == -1) { #No input from cmd
    if ( defined($script_parameters{'pPicardToolsMarkduplicates'}) ) { #Input from config file - do nothing	
    }
    else {
	$script_parameters{'pPicardToolsMarkduplicates'} = 1; #Default
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $script_parameters{'pPicardToolsMarkduplicates'} = $pPicT_MarkDup;
}

#
#Coverage Calculations
#
if ( $pCC == -1) { #No input from cmd
    if ( defined($script_parameters{'pCalculateCoverage'}) ) { #Input from config file - do nothing	
    }
    else {
	$script_parameters{'pCalculateCoverage'} = 1; #Default
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $script_parameters{'pCalculateCoverage'} = $pCC;
}
if ( $script_parameters{'pCalculateCoverage'} == 1) {
    if ( $pCC_Bedgc == -1) { #No input from cmd
	if ( defined($script_parameters{'pGenomeCoverageBED'}) ) { #Input from config file - do nothing	
	}
	else {
	    $script_parameters{'pGenomeCoverageBED'} = 1; #Default
	}
    }
    else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
	$script_parameters{'pGenomeCoverageBED'} = $pCC_Bedgc;
    }
    if ( $pCC_Bedc == -1) { #No input from cmd
	if ( defined($script_parameters{'pCoverageBED'}) ) { #Input from config file - do nothing	
	}
	else {
	    $script_parameters{'pCoverageBED'} = 1; #Default
	}
    }
    else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
	$script_parameters{'pCoverageBED'} = $pCC_Bedc;
    }
    if ( $pCC_Qac == -1) { #No input from cmd
	if ( defined($script_parameters{'pQaCompute'}) ) { #Input from config file - do nothing	
	}
	else {
	    $script_parameters{'pQaCompute'} = 1; #Default
	}
    }
    else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
	$script_parameters{'pQaCompute'} = $pCC_Qac;
    }
    if ( $pCC_PicMM == -1) { #No input from cmd
	if ( defined($script_parameters{'pPicardToolsCollectMultipleMetrics'}) ) { #Input from config file - do nothing	
	}
	else {
	    $script_parameters{'pPicardToolsCollectMultipleMetrics'} = 1; #Default
	}
    }
    else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
	$script_parameters{'pPicardToolsCollectMultipleMetrics'} = $pCC_PicMM;
    }
    if ( $pCCE_PicHS == -1) { #No input from cmd
	if ( defined($script_parameters{'pPicardToolsCalculateHSMetrics'}) ) { #Input from config file - do nothing	
	}
	else {
	    $script_parameters{'pPicardToolsCalculateHSMetrics'} = 1; #Default
	}
    }
    else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
	$script_parameters{'pPicardToolsCalculateHSMetrics'} = $pCCE_PicHS;
    }
}

#
#RCovPlots
#
if ( $pRCP == -1) { #No input from cmd
    if ( defined($script_parameters{'pRCovPlots'}) ) { #Input from config file - do nothing	
    }
    else {
	$script_parameters{'pRCovPlots'} = 1; #Default
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $script_parameters{'pRCovPlots'} = $pRCP;
}

#
#RemovalRedundantFiles
#
if ( $pREM == -1) { #No input from cmd
    if ( defined($script_parameters{'pRemovalRedundantFiles'}) ) { #Input from config file - do nothing	
    }
    else {
	$script_parameters{'pRemovalRedundantFiles'} = 1; #Default
    }
}
else { #Add to enable or overwrite info gathered from config and use in recreation of cmd line later
    $script_parameters{'pRemovalRedundantFiles'} = $pREM;
}

if ($humanGenomeReference eq 0) { #No input from cmd

    if ( ($script_parameters{'pPicardToolsCalculateHSMetrics'} ==1) || ($script_parameters{'pPicardToolsCollectMultipleMetrics'} ==1) || ($script_parameters{'pBwaAln'} ==1)|| ($script_parameters{'pBwaSampe'} ==1) ) { #Programs that uses $humanGenomeReference
	if ( defined($script_parameters{'humanGenomeReference'}) ) { #Input from config file - do nothing
	    if ( $script_parameters{'humanGenomeReference'} =~/^Homo_sapiens.GRCh(\d+\.\d+)/) { #Used to change capture kit genome reference version later
		$humanGenomeReferenceVersion = $1;
		$humanGenomeReferenceSource = "GRCh"; #Ensembl
		$humanGenomeRefereceChromosomePrefix = "nochr";
	    }
	    elsif ( $script_parameters{'humanGenomeReference'} =~/^Homo_sapiens.hg(\d+)/) { #Used to change capture kit genome reference version later
		$humanGenomeReferenceVersion = $1;
		$humanGenomeReferenceSource = "hg"; #Refseq
		$humanGenomeRefereceChromosomePrefix = "chr";
	    }	
	}
	elsif ($environmentUppmax == 1) {
	    $humanGenomeReference = "Homo_sapiens.GRCh37.70_nochr.fasta";
	    if ( $humanGenomeReference =~/Homo_sapiens.GRCh(\d+\.\d+)/) {
		$humanGenomeReferenceVersion = $1;
		$humanGenomeReferenceSource = "GRCh"; #Ensembl
		$humanGenomeRefereceChromosomePrefix = "nochr";
	    }   
	    $script_parameters{'humanGenomeReference'} = $humanGenomeReference; #Add to enable recreation of cmd line later
	}
	else {
	    print STDERR "\nSupply human genome reference fasta file to run 'pCCE_PicHS' and/or 'pCC_PicMM'\n\n";
	    die $USAGE
	}
    }
    else {
	$script_parameters{'humanGenomeReference'} = $humanGenomeReference; #Add to enable recreation of cmd line later
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
    $script_parameters{'humanGenomeReference'} = $humanGenomeReference;
}
unless (-e "$script_parameters{'referencesDir'}/$script_parameters{'humanGenomeReference'}" ) { #Check for human genome reference in supplied reference dir
    print STDERR "\nCould not find human genome reference fasta file: $script_parameters{'referencesDir'}/$script_parameters{'humanGenomeReference'}\n\n";
    die $USAGE;
}

if ($picardPath eq 0 ) { #No input from cmd
    if ( ($script_parameters{'pPicardToolsMerge'} ==1) || ($script_parameters{'pPicardToolsMarkduplicates'} ==1) || ($script_parameters{'pPicardToolsCalculateHSMetrics'} ==1) || ($script_parameters{'pPicardToolsCollectMultipleMetrics'} ==1) ) {
	if ( defined($script_parameters{'picardPath'}) ) { #Input from config file - do nothing	
	}
	elsif ($environmentUppmax == 1) {
	    $picardPath = "/bubo/home/h12/henriks/programs/picard-tools-1.59";
	    print STDOUT "\nSet picard path to: $picardPath\n";
	    $script_parameters{'picardPath'} = $picardPath; #Add to enable recreation of cmd line later
	}
	else {
	    print STDERR "\nMust supply picardTools path to use PicardTools", "\n";
	    die $USAGE;
	}
    }
}
else { #Add to enable recreation of cmd line later
    $script_parameters{'picardPath'} = $picardPath;
}

#
#exomeTargetBed
#

if ( $exomeTargetBed eq 0 ) { #No input from cmd
    if ( $script_parameters{'pCoverageBED'} ==1) {
	my $uncorrect_capture_Counter = 0; #Track no entries or wrong format entry in pedigree file
	for (my $sampleid_Counter=0;$sampleid_Counter<scalar(@sampleIDs);$sampleid_Counter++) { #Check all samples
	    if ( defined($script_parameters{ $sampleIDs[$sampleid_Counter] }{'exomeTargetBed'}) ) { #Input from config file - transfer to sample_info
		$sample_info{$script_parameters{'familyID'}}{$sampleIDs[$sampleid_Counter]}{'exomeTargetBed'} = $script_parameters{ $sampleIDs[$sampleid_Counter] }{'exomeTargetBed'};
	    }	
	    elsif ($environmentUppmax == 1) {
		
		if ( defined( $sample_info{ $script_parameters{'familyID'} }{$sampleIDs[$sampleid_Counter]}{'exomeTargetBed'} ) ) { #Capture kit check
		    $sample_info{ $script_parameters{'familyID'} }{ $sampleIDs[$sampleid_Counter] }{'exomeTargetBed'} =~ s/GenomeReferenceSource/$humanGenomeReferenceSource/; #Replace with Refseq genome or Ensembl genome
		    $sample_info{ $script_parameters{'familyID'} }{ $sampleIDs[$sampleid_Counter] }{'exomeTargetBed'} =~ s/Version/$humanGenomeReferenceVersion/; #Replace with actual version 
		    $sample_info{ $script_parameters{'familyID'} }{ $sampleIDs[$sampleid_Counter] }{'exomeTargetBed'} =~ s/ChromosomePrefix/$humanGenomeRefereceChromosomePrefix/; #Replace with chromosome prefix
		    $script_parameters{ $sampleIDs[$sampleid_Counter] }{'exomeTargetBed'} = $sample_info{ $script_parameters{'familyID'} }{ $sampleIDs[$sampleid_Counter] }{'exomeTargetBed'}; #Add to enable recreation of cmd line later
		}
		else {
		    print STDERR "\nCould not find a capture kit entry for sample: $sampleIDs[$sampleid_Counter]\n";
		    print STDERR "\nYou must supply a capture kit bed file when running 'pCC_Bedc'\n";
		    $uncorrect_capture_Counter++;
		}
	    }
	    else { #No capture kit information   
		print STDERR "\nYou must supply a capture kit bed file when running 'pCC_Bedc'\n";
		print STDERR "\n";
		die $USAGE;
	    }
	}
	if ($uncorrect_capture_Counter > 0) { #If lacking or not supported
	    print STDERR "\nChange/add capture kit record in pedigree file: $script_parameters{'pedigreeFile'}\n";
	    print STDERR "List of pedigree supported capture kits records:\n\n";
	    print STDERR "Pedigree record", "\t", "Capture kit BED-file\n";
	    for my $supported_capture_kit (keys %supported_capture_kits) {
		print STDERR $supported_capture_kit, "\t", $supported_capture_kits{$supported_capture_kit}, "\n";
	    }	    
	    print STDERR "\n";
	    die $USAGE;
	}
	for (my $sampleid_Counter=0;$sampleid_Counter<scalar(@sampleIDs);$sampleid_Counter++) { #Check that all samples in pedigree have a capture kit file
	    unless (-e "$script_parameters{'referencesDir'}/$sample_info{$script_parameters{'familyID'}}{$sampleIDs[$sampleid_Counter]}{'exomeTargetBed'}" ) { #Check existence of capture file in supplied reference dir
		print STDERR "\nCould not find capture kit target BED-file: $script_parameters{'referencesDir'}/$sample_info{$script_parameters{'familyID'}}{$sampleIDs[$sampleid_Counter]}{'exomeTargetBed'}\n\n";
		die $USAGE;
	    }	
	    if ($sample_info{$script_parameters{'familyID'}}{$sampleIDs[0]}{'exomeTargetBed'} eq $sample_info{$script_parameters{'familyID'}}{$sampleIDs[$sampleid_Counter]}{'exomeTargetBed'}) { #Enable later print of cmd to be executable if all IDN have been sequenced using the same capture kit
		$identicalCaptureBedCounter++;
	    }
	} 	
    }
}
else {
    for (my $sampleid_Counter=0;$sampleid_Counter<scalar(@sampleIDs);$sampleid_Counter++) { #Add capture kit to all samples
	$sample_info{ $script_parameters{'familyID'} }{ $sampleIDs[$sampleid_Counter] }{'exomeTargetBed'} = $exomeTargetBed; #Add capture target BED-file to sample_info info to enable individal adjusted capture calculation for each family member
#Check for file existence
	unless (-e "$script_parameters{'referencesDir'}/$sample_info{ $script_parameters{'familyID'} }{ $sampleIDs[$sampleid_Counter] }{'exomeTargetBed'}" ) { #Check for capture file in supplied reference dir
	    print STDERR "\nCould not find capture kit target BED-file: $script_parameters{'referencesDir'}/$sample_info{ $script_parameters{'familyID'} }{ $sampleIDs[$sampleid_Counter] }{'exomeTargetBed'}\n\n";
	    die $USAGE;		
	}
	$script_parameters{ $sampleIDs[$sampleid_Counter] }{'exomeTargetBed'} = $exomeTargetBed; #Add to enable recreation of cmd line later
    }
}

#
#exomeTargetBedInfileList
#
if ( $exomeTargetBedInfileList eq 0 ) {
    
    if ( $script_parameters{'pPicardToolsCalculateHSMetrics'} ==1 ) {
	my $uncorrect_capture_Counter = 0; #Track no entries or wrong format entry in pedigree file
	for (my $sampleid_Counter=0;$sampleid_Counter<scalar(@sampleIDs);$sampleid_Counter++) { #Check all samples
	    if ( defined($script_parameters{ $sampleIDs[$sampleid_Counter] }{'exomeTargetBedInfileList'}) ) { #Input from config file - transfer to sample_info
		$sample_info{$script_parameters{'familyID'}}{$sampleIDs[$sampleid_Counter]}{'exomeTargetBedInfileList'} = $script_parameters{ $sampleIDs[$sampleid_Counter] }{'exomeTargetBedInfileList'};
	    }	
	    elsif ($environmentUppmax == 1) {		
		
		if ( defined( $sample_info{$script_parameters{'familyID'} }{$sampleIDs[$sampleid_Counter]}{'exomeTargetBedInfileList'} ) ) { #Capture kit check
		    #Add actual run info to the capture kit file names
		    $sample_info{ $script_parameters{'familyID'} }{ $sampleIDs[$sampleid_Counter] }{'exomeTargetBedInfileList'} =~ s/GenomeReferenceSource/$humanGenomeReferenceSource/; #Replace with Refseq genome or Ensembl genome
		    $sample_info{ $script_parameters{'familyID'} }{ $sampleIDs[$sampleid_Counter] }{'exomeTargetBedInfileList'} =~ s/Version/$humanGenomeReferenceVersion/; #Replace with actual version
		    $sample_info{ $script_parameters{'familyID'} }{ $sampleIDs[$sampleid_Counter] }{'exomeTargetBedInfileList'} =~ s/ChromosomePrefix/$humanGenomeRefereceChromosomePrefix/; #Replace with chromosome prefix
		    $script_parameters{ $sampleIDs[$sampleid_Counter] }{'exomeTargetBedInfileList'} = $sample_info{ $script_parameters{'familyID'} }{ $sampleIDs[$sampleid_Counter] }{'exomeTargetBedInfileList'}; #Add to enable recreation of cmd line later
		}
		else {
		    print STDERR "\nCould not find a capture kit entry for sample: $sampleIDs[$sampleid_Counter]\n";
		    print STDERR "\nYou must supply a capture kit target bed infile_list when running 'pCCE_PicHS'\n";
		    $uncorrect_capture_Counter++;
		}
	    }
	    else { #Lacking capture kit information   
		print STDERR "\nYou must supply a capture kit target bed infile_list when running 'pCCE_PicHS'\n";
		print STDERR "\n";
		die $USAGE;
	    }
	}
	if ($uncorrect_capture_Counter > 0) { #If lacking or not supported
	    print STDERR "\nChange/add capture kit record in pedigree file: $script_parameters{'pedigreeFile'}\n";
	    print STDERR "List of pedigree supported capture kits records:\n\n";
	    print STDERR "Pedigree record", "\t", "Capture kit BED-file\n";
	    for my $supported_capture_kit (keys %supported_capture_kits) {
		print STDERR $supported_capture_kit, "\t", $supported_capture_kits{$supported_capture_kit}, "\n";
	    }	 
	    print "\n";
	    die $USAGE;
	}
	for (my $sampleid_Counter=0;$sampleid_Counter<scalar(@sampleIDs);$sampleid_Counter++) { #Check all samples
	    unless (-e "$script_parameters{'referencesDir'}/$sample_info{ $script_parameters{'familyID'} }{ $sampleIDs[$sampleid_Counter] }{'exomeTargetBedInfileList'}" ) { #Check existence of capture file in supplied reference dir
		print STDERR "\nCould not find capture kit target bed infile_list file: $script_parameters{'referencesDir'}/$sample_info{ $script_parameters{'familyID'} }{ $sampleIDs[$sampleid_Counter] }{'exomeTargetBedInfileList'}\n\n";
		die $USAGE;		
	    }
	    if ($sample_info{$script_parameters{'familyID'}}{$sampleIDs[0]}{'exomeTargetBedInfileList'} eq $sample_info{ $script_parameters{'familyID'} }{ $sampleIDs[$sampleid_Counter] }{'exomeTargetBedInfileList'}) { #Enable later print of cmd to be executable if all IDN have been sequenced using the same capture kit
		$identicalCaptureBedIntervalCounter++;
	    }
	}
    }
}   
else {
    
    for (my $sampleid_Counter=0;$sampleid_Counter<scalar(@sampleIDs);$sampleid_Counter++) { #Add capture kit to all samples using user supplied info
	$sample_info{ $script_parameters{'familyID'} }{ $sampleIDs[$sampleid_Counter] }{'exomeTargetBedInfileList'} = $exomeTargetBedInfileList; #Add capture target file to sample_info  to enable individal adjusted capture calculation for each family member
	
#Check for file existence 
	unless (-e "$script_parameters{'referencesDir'}/$sample_info{ $script_parameters{'familyID'} }{ $sampleIDs[$sampleid_Counter] }{'exomeTargetBedInfileList'}" ) { #Check for capture file in supplied reference dir
	    print STDERR "\nCould not find capture kit target file: $script_parameters{'referencesDir'}/$sample_info{ $script_parameters{'familyID'} }{ $sampleIDs[$sampleid_Counter] }{'exomeTargetBedInfileList'}\n\n";
	    die $USAGE;		
	}
	$script_parameters{$sampleIDs[$sampleid_Counter]}{'exomeTargetBedInfileList'} = $sample_info{ $script_parameters{'familyID'} }{ $sampleIDs[$sampleid_Counter] }{'exomeTargetBedInfileList'}; #Add to enable recreation of cmd line later
    }
}

#
#exomeTargetPaddedBedInfileList
#
if ( $exomeTargetPaddedBedInfileList eq 0 ) {
    if ( $script_parameters{'pPicardToolsCalculateHSMetrics'} ==1 ) {
	my $uncorrect_capture_Counter = 0; #Track no entries or wrong format entry in pedigree file
	for (my $sampleid_Counter=0;$sampleid_Counter<scalar(@sampleIDs);$sampleid_Counter++) { #Check all samples
	    if ( defined($script_parameters{ $sampleIDs[$sampleid_Counter] }{'exomeTargetPaddedBedInfileList'}) ) { #Input from config file - transfer to sample_info
		$sample_info{$script_parameters{'familyID'}}{$sampleIDs[$sampleid_Counter]}{'exomeTargetPaddedBedInfileList'} = $script_parameters{ $sampleIDs[$sampleid_Counter] }{'exomeTargetPaddedBedInfileList'};
	    }	
	    elsif ($environmentUppmax == 1) {		
		if ( defined( $sample_info{$script_parameters{'familyID'} }{$sampleIDs[$sampleid_Counter]}{'exomeTargetPaddedBedInfileList'} ) ) { #Capture kit check
		    #Add actual run info to the capture kit file names
		    $sample_info{ $script_parameters{'familyID'} }{ $sampleIDs[$sampleid_Counter] }{'exomeTargetPaddedBedInfileList'} =~ s/GenomeReferenceSource/$humanGenomeReferenceSource/; #Replace with Refseq genome or Ensembl genome
		    $sample_info{ $script_parameters{'familyID'} }{ $sampleIDs[$sampleid_Counter] }{'exomeTargetPaddedBedInfileList'} =~ s/Version/$humanGenomeReferenceVersion/; #Replace with actual version
		    $sample_info{ $script_parameters{'familyID'} }{ $sampleIDs[$sampleid_Counter] }{'exomeTargetPaddedBedInfileList'} =~ s/ChromosomePrefix/$humanGenomeRefereceChromosomePrefix/; #Replace with chromosome prefix
		    $script_parameters{ $sampleIDs[$sampleid_Counter] }{'exomeTargetPaddedBedInfileList'} = $sample_info{ $script_parameters{'familyID'} }{ $sampleIDs[$sampleid_Counter] }{'exomeTargetPaddedBedInfileList'}; #Add to enable recreation of cmd line later
		}
		else {
		    print STDERR "\nCould not find a capture kit target bed padded infile_list entry for sample: $sampleIDs[$sampleid_Counter]\n";
		    print STDERR "\nYou must supply a capture kit target bed padded infile_list when running 'pCCE_PicHS'\n";
		    $uncorrect_capture_Counter++;
		}
	    }
	    else { #Lacking capture kit information   
		print STDERR "\nYou must supply a capture kit target bed padded infile_list when running 'pCCE_PicHS'\n";
		print STDERR "\n";
		die $USAGE;
	    }
	}
	if ($uncorrect_capture_Counter > 0) { #If lacking or not supported
	    print STDERR "\nChange/add capture kit record in pedigree file: $script_parameters{'pedigreeFile'}\n";
	    print STDERR "List of pedigree supported capture kits records:\n\n";
	    print STDERR "Pedigree record", "\t", "Capture kit BED-file\n";
	    for my $supported_capture_kit (keys %supported_capture_kits) {
		print STDERR $supported_capture_kit, "\t", $supported_capture_kits{$supported_capture_kit}, "\n";
	    }	 
	    print "\n";
	    die $USAGE;
	}
	for (my $sampleid_Counter=0;$sampleid_Counter<scalar(@sampleIDs);$sampleid_Counter++) { #Check all samples
	    unless (-e "$script_parameters{'referencesDir'}/$sample_info{ $script_parameters{'familyID'} }{ $sampleIDs[$sampleid_Counter] }{'exomeTargetPaddedBedInfileList'}" ) { #Check existence of capture file in supplied reference dir
		print STDERR "\nCould not find padded capture kit target file: $script_parameters{'referencesDir'}/$sample_info{ $script_parameters{'familyID'} }{ $sampleIDs[$sampleid_Counter] }{'exomeTargetPaddedBedInfileList'}\n\n";
		die $USAGE;		
	    }
	    if ($sample_info{$script_parameters{'familyID'}}{$sampleIDs[0]}{'exomeTargetPaddedBedInfileList'} eq $sample_info{ $script_parameters{'familyID'} }{ $sampleIDs[$sampleid_Counter] }{'exomeTargetPaddedBedInfileList'}) { #Enable later print of cmd to be executable if all IDN have been sequenced using the same capture kit
		$identicalCaptureBedPaddedIntervalCounter++;
	    }
	}
    }
}   
else {
    
    for (my $sampleid_Counter=0;$sampleid_Counter<scalar(@sampleIDs);$sampleid_Counter++) { #Add capture kit to all samples using user supplied info
	$sample_info{ $script_parameters{'familyID'} }{ $sampleIDs[$sampleid_Counter] }{'exomeTargetPaddedBedInfileList'} = $exomeTargetPaddedBedInfileList; #Add padded capture target file to sample_info  to enable individal adjusted capture calculation for each family member
	
#Check for file existence 
	unless (-e "$script_parameters{'referencesDir'}/$sample_info{ $script_parameters{'familyID'} }{ $sampleIDs[$sampleid_Counter] }{'exomeTargetPaddedBedInfileList'}" ) { #Check for capture file in supplied reference dir
	    print STDERR "\nCould not find padded capture kit target file: $script_parameters{'referencesDir'}/$sample_info{ $script_parameters{'familyID'} }{ $sampleIDs[$sampleid_Counter] }{'exomeTargetPaddedBedInfileList'}\n\n";
	    die $USAGE;		
	}
	$script_parameters{$sampleIDs[$sampleid_Counter]}{'exomeTargetPaddedBedInfileList'} = $sample_info{ $script_parameters{'familyID'} }{ $sampleIDs[$sampleid_Counter] }{'exomeTargetPaddedBedInfileList'}; #Add to enable recreation of cmd line later
    }
}

if( $writeConfigFile eq 0) { #No input from cmd
    if ( defined($script_parameters{'writeConfigFile'}) ) { #Input from config file - do nothing
    }
}
else {
    $script_parameters{'writeConfigFile'} = $writeConfigFile;
}
if ( $script_parameters{'writeConfigFile'} ne 0) { #Write config file for family
    open (YAML, '>', $script_parameters{'writeConfigFile'}) or die "can't open $script_parameters{'writeConfigFile'}: $!";
    print YAML Dump(%script_parameters), "\n";
    close (YAML);
}
###
#Creates master_logg for the master script 
###
`mkdir -p $script_parameters{'outDataDir'}/$script_parameters{'familyID'}/master_logg;`; #Creates the master_logg dir
my ($base,$script) = (`date +%Y%m%d`,`basename $0`); #Catches current date and script name
chomp($base,$script); #Remove \n;
my $master_logg_name="$script_parameters{'outDataDir'}/$script_parameters{'familyID'}/master_logg/$script"."_"."$base.txt"; #concatenates master_logg filename
open (MASTERL, ">>$master_logg_name") or die "Can't write to $master_logg_name: $!\n"; #Open file run logg
#Add parameters
print MASTERL "\n$script "; #Adds script name to recontruct command line
WriteCMDMasterLogg();

print STDOUT "\nScript parameters and info from $script are saved in file: $master_logg_name", "\n";

#@inFilesDir = split(/,/,join(',',@inFilesDir)); #Enables comma separated indir(s)
@inFilesDir = split(/,/,$script_parameters{'inFilesDir'}); #Enables comma separated indir(s)
@picardToolsMergePreviousFiles = split(/,/,join(',',@picardToolsMergePreviousFiles)); #Enables comma separated list of previously generated _sorted.bam files

for (my $inputdir=0;$inputdir<scalar(@inFilesDir);$inputdir++) { #Collects inputfiles
    
    my @infiles = `cd $inFilesDir[ $inputdir ];ls *.fastq*;`; #cd to input dir and collect fastq files and fastq.gz files
   
    print STDOUT "\nReads from Platform", "\n";print MASTERL "\nReads from Platform", "\n";
    print STDOUT "\nSample ID", "\t", $sampleIDs[$inputdir],"\n";print MASTERL "\nSample ID", "\t", $sampleIDs[$inputdir],"\n";
    print STDOUT "Inputfiles", "\n", @ { $infiles{ $sampleIDs[$inputdir] }  =[@infiles] }, "\n"; #hash with sample id as key and inputfiles in dir as array 
    print MASTERL "Inputfiles", "\n", @ { $infiles{ $sampleIDs[$inputdir] }  =[@infiles] }, "\n";
    $indirpath{$sampleIDs[$inputdir]} = $inFilesDir[ $inputdir ];  #Catch inputdir path
    chomp(@infiles);    #Remove newline from every entry in array
    $infiles{ $sampleIDs[$inputdir] }  =[@infiles]; #Reload files into hash (kept above newline just for print STDOUT)
}
close(MASTERL);

#Set chr prefix and chromosome names depending on reference used
if ($humanGenomeReference=~/hg\d+/) { #Refseq - prefix and M
    @chromosomes = ("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM"); #Chr for filtering of bam file
}
elsif ($humanGenomeReference=~/GRCh\d+/) { #Ensembl - no prefix and MT
    @chromosomes = ("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT"); #Chr for filtering of bam file
}

my $uncompressed_file_switch = InfilesReFormat(); #Required to format infiles correctly for subsequent input into aligners

#########################
###Run program part######
#########################

open (MASTERL, ">>$master_logg_name") or die "Can't write to $master_logg_name: $!\n"; #Open file run logg

if ( ($script_parameters{'pGZip'} eq 1) && ($uncompressed_file_switch eq 1) ) { #GZip of fastq files

    print STDOUT "\nGZip for fastq files", "\n";print MASTERL "\nGZip for fastq files", "\n";

    for (my $sampleid_Counter=0;$sampleid_Counter<scalar(@sampleIDs);$sampleid_Counter++) {  

	for (my $infile=0;$infile<scalar( @{ $infiles{$sampleIDs[$sampleid_Counter]} });$infile++) { #To determine which sampleID had the uncompressed files
	    
	    if ($infiles{$sampleIDs[$sampleid_Counter]}[$infile] =~/.fastq$/) {
	
		GZipfastq($sampleIDs[$sampleid_Counter]);
		last; #Return to sampleID loop i.e. only call sunroutine GZipfastq once per sampleID
	    }
	}
    }
}

if ($script_parameters{'pFastQC'} eq 1) { #FASTQC
    
    print STDOUT "\nFastQC", "\n";print MASTERL "\nFastQC", "\n";
    
    for (my $sampleid_Counter=0;$sampleid_Counter<scalar(@sampleIDs);$sampleid_Counter++) {  
	
	FastQC($sampleIDs[$sampleid_Counter]);	
    }
}

if ($script_parameters{'pMosaikBuild'} eq 1) { #Run MosaikBuild
    
    print STDOUT "\nMosaikBuild", "\n";print MASTERL "\nMosaikBuild", "\n";

    for (my $sampleid_Counter=0;$sampleid_Counter<scalar(@sampleIDs);$sampleid_Counter++) {  
	
	MosaikBuild($sampleIDs[$sampleid_Counter]);	
    }
}


if ($script_parameters{'pMosaikAlign'} eq 1) { #Run MosaikAlign
    
    print STDOUT "\nMosaikAlign", "\n"; print MASTERL "\nMosaikAlign", "\n";

    for (my $sampleid_Counter=0;$sampleid_Counter<scalar(@sampleIDs);$sampleid_Counter++) {  
	
	MosaikAlign($sampleIDs[$sampleid_Counter]);	
    }
}

if ($script_parameters{'pBwaAln'} eq 1) { #Run bwa aln
    
    print STDOUT "\nBWA aln", "\n";print MASTERL "\nBWA aln", "\n";
    
    for (my $sampleid_Counter=0;$sampleid_Counter<scalar(@sampleIDs);$sampleid_Counter++) {  
	
	BWA_aln($sampleIDs[$sampleid_Counter]);	
    }    
}

if ($script_parameters{'pBwaSampe'} eq 1) { #Run bwa sampe
    
    print STDOUT "\nBWA sampe", "\n";print MASTERL "\nBWA sampe", "\n";
    
    for (my $sampleid_Counter=0;$sampleid_Counter<scalar(@sampleIDs);$sampleid_Counter++) {  
	
	BWA_sampe($sampleIDs[$sampleid_Counter]);
    }
}

if ($script_parameters{'pSamToolsSort'} eq 1) { #Run samtools index and sort

    print STDOUT "\nSamtools sort & index", "\n";

    for (my $sampleid_Counter=0;$sampleid_Counter<scalar(@sampleIDs);$sampleid_Counter++) {  
    
	SamtoolsSortIndex($sampleIDs[$sampleid_Counter], $script_parameters{'aligner'});
	
    }
}

if ($script_parameters{'pPicardToolsMerge'} eq 1) { #Run picardtools merge (Requires sorted bam files)

    print STDOUT "\nPicardtools MergeSamFiles", "\n";

    for (my $sampleid_Counter=0;$sampleid_Counter<scalar(@sampleIDs);$sampleid_Counter++) {  

	PicardMerge($sampleIDs[$sampleid_Counter], $script_parameters{'aligner'});	
    }
}

if ($script_parameters{'pPicardToolsMarkduplicates'} eq 1) { #PicardTools MarkDuplicates

    print STDOUT "\nPicard MarkDuplicates", "\n";print MASTERL "\nPicard MarkDuplicates", "\n";

    for (my $sampleid_Counter=0;$sampleid_Counter<scalar(@sampleIDs);$sampleid_Counter++) {  
    
	PicardMarkDup($sampleIDs[$sampleid_Counter], $script_parameters{'aligner'});	
    }
}

if ($script_parameters{'pCalculateCoverage'} eq 1) { #Run GenomeCoverageBED, qaCompute (Paul Costea), Picard (CollectAlignmentSummaryMetrics, CalculateHsMetrics)
    
    print STDOUT "\nCalculate Coverage", "\n";print MASTERL "\nCalculate Coverage", "\n";    
    
    for (my $sampleid_Counter=0;$sampleid_Counter<scalar(@sampleIDs);$sampleid_Counter++) {  

	Cal_Coverage($sampleIDs[$sampleid_Counter], $script_parameters{'aligner'}, $script_parameters{'familyID'}); #FamilyID for reporting coverage to CMMS	
    }
}

if ( ($script_parameters{'pRCovPlots'} eq 1) && ( $script_parameters{'pCalculateCoverage'} eq 1) ) { #Run Rcovplot scripts after calculateCoverage with option genomeCoverageBED Rscript: covplots_genome.R & coverageBED Rscript:covplots_exome_all.R      
    print STDOUT "\nRCovPlots", "\n";print MASTERL "\nRCovPlots", "\n";	

    for (my $sampleid_Counter=0;$sampleid_Counter<scalar(@sampleIDs);$sampleid_Counter++) {  
	
	RcoveragePlots($sampleIDs[$sampleid_Counter], $script_parameters{'aligner'});	
    }
}

if ($script_parameters{'pRemovalRedundantFiles'} eq 1) { #Sbatch generation of removal of alignment files
    
    print STDOUT "\nRemoval of alignment files", "\n"; print MASTERL "\nRemoval of alignment files", "\n";
    
    for (my $sampleid_Counter=0;$sampleid_Counter<scalar(@sampleIDs);$sampleid_Counter++) {  
	
	Rem($sampleIDs[$sampleid_Counter], $script_parameters{'aligner'});	
    }
}

close(MASTERL); #Close Master_logg file

######################
###Sub Routines#######
######################

sub Rem {
#Generates a sbatch script, which removes alignment files after they have been sorted and indexed.
#$_[0] = sampleid
#$_[1] = aligner

    `mkdir -p $script_parameters{'outDataDir'}/$_[0]/$_[1]/info;`; #Creates the aligner and info data file directory
    `mkdir -p $script_parameters{'outScriptDir'}/$_[0]/$_[1];`; #Creates the aligner script directory
    $filename = "$script_parameters{'outScriptDir'}/$_[0]/$_[1]/rem_$_[0].";

    Checkfnexists($filename, $fnend);

#Info and Logg
    print STDOUT "Creating sbatch script Remove Aligned Files and writing script file(s) to: ", $filename, "\n";print MASTERL "Creating sbatch script Remove Aligned Files and writing script file(s) to: ", $filename, "\n";
    print STDOUT "Sbatch script Remove Aligned data files will be removed in: ", $script_parameters{'outDataDir'},"/$_[0]/$_[1]", "\n";print MASTERL "Sbatch script Remove Aligned data files will be removed in: ", $script_parameters{'outDataDir'},"/$_[0]/$_[1]", "\n";

    open (REM, ">$filename") or die "Can't write to $filename: $!\n";
    
    print REM "#! /bin/bash -l", "\n";
    print REM "#SBATCH -A ", $script_parameters{'projectID'}, "\n";
    print REM "#SBATCH -n 1", "\n";
    print REM "#SBATCH -C thin", "\n";
    print REM "#SBATCH -t 00:15:00", "\n";
    print REM "#SBATCH -J RMMoB_", $_[0], "\n";

    print REM "#SBATCH -e $script_parameters{'outDataDir'}/$_[0]/$_[1]/info/rem_$_[0].", $fnt ,".stderr.txt", "\n";
    print REM "#SBATCH -o $script_parameters{'outDataDir'}/$_[0]/$_[1]/info/rem_$_[0].", $fnt ,".stdout.txt", "\n";

    unless ($em eq 0) {
	
	print REM "#SBATCH --mail-type=END", "\n";
	print REM "#SBATCH --mail-type=FAIL", "\n";
	print REM "#SBATCH --mail-user=$em", "\n\n";
	
    }
    print REM 'echo "Running on: $(hostname)"',"\n\n";
    print REM "cd $script_parameters{'outDataDir'}/$_[0]/$_[1]", "\n\n";
    print REM "#Samples", "\n\n";
    print REM 'inSampleDir="',"$script_parameters{'outDataDir'}/$_[0]/$_[1]", '"', "\n\n";
###
#Remove Mosaik files
###
    if ($script_parameters{'pMosaikBuild'} || $script_parameters{'pMosaikAlign'} || ($script_parameters{'aligner'} eq "mosaik")) {
	for (my $infile=0;$infile < scalar( @{ $Infiles_lane_noending{$_[0]} } );$infile++) { #MosaikBuild takes both reads at once
	    
	    my $tempinfile = $Infiles_lane_noending{$_[0]}[$infile]; 
	    print REM "rm  ", '${inSampleDir}', "/$tempinfile", ".dat", "\n\n"; #MosaikBuild
	    print REM "rm  ", '${inSampleDir}', "/$tempinfile", ".stat", "\n\n"; #MosaikAlign
	    print REM "rm  ", '${inSampleDir}', "/$tempinfile", ".bam", "\n\n"; #MosaikAlign
	    #print REM "rm  ", '${inSampleDir}', "/$tempinfile", ".multiple.bam", "\n\n"; #MosaikAlign
	    #print REM "rm  ", '${inSampleDir}', "/$tempinfile", "_sorted.bam", "\n\n"; #MosaikAlign/samtools
	    #print REM "rm  ", '${inSampleDir}', "/$tempinfile", "_sorted.bam.bai", "\n\n"; #MosaikAlign/samtools
	    print REM "rm  ", '${inSampleDir}', "/coverageReport/$tempinfile", "_sorted_pmd_coverageBed", "\n\n"; #Coverage of features in BED-file
	    print REM "rm  ", '${inSampleDir}', "/coverageReport/$tempinfile", "_sorted_pmd_coverageBed_hist", "\n\n"; #bedtools histogram of BED-file
	    print REM "rm  ", '${inSampleDir}', "/coverageReport/$tempinfile", "_sorted_pmd_rmdup_coverageBed_hist", "\n\n"; #bedtools histogram of BED-file
	    #print REM "rm  ", '${inSampleDir}', "/coverageReport/$tempinfile", "_sorted_pmd_coverageBed_depth_pos", "\n\n"; #Coverage of features BED-file per position
	    if ( ($script_parameters{'pPicardToolsMerge'} eq 1) && ( $infile == scalar( @{ $Infiles_lane_noending{$_[0]} } )-1 ) && ($infile >= 1) ) { #If merged, last file and that there exits at least 2 files
		print REM "rm  ", '${inSampleDir}', "/coverageReport/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, "_sorted_merged_pmd_coverageBed", "\n\n"; #Coverage of features in BED-file
		print REM "rm  ", '${inSampleDir}', "/coverageReport/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, "_sorted_merged_pmd_coverageBed_hist", "\n\n"; #bedtools histogram of BED-file
		print REM "rm  ", '${inSampleDir}', "/coverageReport/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, "_sorted_merged_pmd_rmdup_coverageBed_hist", "\n\n"; #bedtools histogram of BED-file
	#	print REM "rm  ", '${inSampleDir}', "/coverageReport/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, "_sorted_merged_pmd_coverageBed_depth_pos", "\n\n"; #Coverage of features BED-file per position
	    }
	}    
	print REM "rm -rf ", '${inSampleDir}', "/per_chr", "\n\n"; #samtools/GATK (real/recal)
    }
###
#Remove BWA files
###
    if ($script_parameters{'pBwaAln'} || $script_parameters{'pBwaSampe'} || ($script_parameters{'aligner'} eq "bwa")) {

	for (my $infile=0;$infile < scalar( @{ $Infiles_bothstrands_noending{$_[0]} } );$infile++) { #BWA_aln takes 1 read at a time 
	    
	    my $tempinfile = $Infiles_bothstrands_noending{$_[0]}[$infile]; 
	    print REM "rm  ", '${inSampleDir}', "/$tempinfile", ".sai", "\n\n"; #BWA_aln
	}
	for (my $infile=0;$infile < scalar( @{ $Infiles_lane_noending{$_[0]} } );$infile++) { #BWA_sampe 
	    
	    my $tempinfile = $Infiles_lane_noending{$_[0]}[$infile]; 
	    print REM "rm  ", '${inSampleDir}', "/$tempinfile", ".bam", "\n\n"; #BWA_sampe
	}    
    }
    close(REM);
    return;
}

sub GZip {
#Generates a sbatch script and runs gZip on aligned reads
#$_[0] = sampleid
#$_[1] = aligner

    `mkdir -p $script_parameters{'outDataDir'}/$_[0]/gzip/info;`; #Creates the gzip folder and info data file directory
    `mkdir -p $script_parameters{'outScriptDir'}/$_[0]/gzip;`; #Creates the gzip script directory

    if ($script_parameters{'pMosaikBuild'} || $script_parameters{'pMosaikAlign'} || ($script_parameters{'aligner'} eq "mosaik")) {
	$filename = "$script_parameters{'outScriptDir'}/$_[0]/gzip/gzipMoA_$_[0].";
    }
    if ($script_parameters{'pBwaAln'} || $script_parameters{'pBwaSampe'} || ($script_parameters{'aligner'} eq "bwa")) {
	$filename = "$script_parameters{'outScriptDir'}/$_[0]/gzip/gzipbwa_sampe_$_[0].";
    }
    Checkfnexists($filename, $fnend);

#Info and Logg
    print STDOUT "Creating sbatch script GzipAligned and writing script file(s) to: ", $filename, "\n";print MASTERL "Creating sbatch script GzipAligned and writing script file(s) to: ", $filename, "\n";
    print STDOUT "Sbatch script GzipAligned data files will be written to: ", $script_parameters{'outDataDir'},"/$_[0]/$_[1]", "\n";print MASTERL "Sbatch script GzipAligned data files will be written to: ", $script_parameters{'outDataDir'},"/$_[0]/$_[1]", "\n";

    my $t = ceil(3*scalar( @{ $Infiles_lane_noending{$_[0]} })); #One full aligned lane takes approx. 3 h for gzip to process, round up to nearest full hour.

    open (GZMOA, ">$filename") or die "Can't write to $filename: $!\n";
    
    print GZMOA "#! /bin/bash -l", "\n";
    print GZMOA "#SBATCH -A ", $script_parameters{'projectID'}, "\n";
    print GZMOA "#SBATCH -p node -n $maximumCores", "\n";
    print GZMOA "#SBATCH -C thin", "\n";
    print GZMOA "#SBATCH -t $t:00:00", "\n";
    print GZMOA "#SBATCH -J GZMoA_", $_[0], "\n";
    
    if ($script_parameters{'pMosaikBuild'} || $script_parameters{'pMosaikAlign'} || ($script_parameters{'aligner'} eq "mosaik")) {
	print GZMOA "#SBATCH -e $script_parameters{'outDataDir'}/$_[0]/gzip/info/gzipMoA_$_[0].", $fnt ,".stderr.txt", "\n";
	print GZMOA "#SBATCH -o $script_parameters{'outDataDir'}/$_[0]/gzip/info/gzipMoA_$_[0].", $fnt ,".stdout.txt", "\n";
    }
    if ($script_parameters{'pBwaAln'} || $script_parameters{'pBwaSampe'} || ($script_parameters{'aligner'} eq "bwa")) {
	print GZMOA "#SBATCH -e $script_parameters{'outDataDir'}/$_[0]/gzip/info/gzipMoA_$_[0].", $fnt ,".stderr.txt", "\n";
	print GZMOA "#SBATCH -o $script_parameters{'outDataDir'}/$_[0]/gzip/info/gzipMoA_$_[0].", $fnt ,".stdout.txt", "\n";
    }

    unless ($em eq 0) {
	
	print GZMOA "#SBATCH --mail-type=END", "\n";
	print GZMOA "#SBATCH --mail-type=FAIL", "\n";
	print GZMOA "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print GZMOA 'echo "Running on: $(hostname)"',"\n\n";
    print GZMOA "cd $script_parameters{'outDataDir'}/$_[0]/$_[1]", "\n\n";
    print GZMOA "#Samples", "\n\n";
    print GZMOA 'inSampleDir="',"$script_parameters{'outDataDir'}/$_[0]/$_[1]", '"', "\n\n";

    my $core_Counter=1;
    for (my $infile=0;$infile < scalar( @{ $Infiles_lane_noending{$_[0]} } );$infile++) { #All alignemnt files
	
	if ($infile eq $core_Counter*$maximumCores) { #Using only $maximumCores cores
	    
	    print GZMOA "wait", "\n\n";
	    $core_Counter=$core_Counter+1;
	}
	my $tempinfile = $Infiles_lane_noending{$_[0]}[$infile]; #Read 1 and read 2
	print GZMOA "gzip  ", '${inSampleDir}', "/$tempinfile", ".bam &", "\n\n";
    }
    print GZMOA "wait", "\n";
    close(GZMOA);
    return;
}

sub RcoveragePlots { 
#Generates sbatch scripts for R scripts:
#1. covplots_genome.R 
# on files generated from calculateCoverage genomeCoverageBED
#$_[0] = sampleid
#$_[1] = aligner
    
    `mkdir -p $script_parameters{'outDataDir'}/$_[0]/$_[1]/info;`; #Creates the aligner folder and info data file directory
    `mkdir -p $script_parameters{'outDataDir'}/$_[0]/$_[1]/coverageReport;`; #Creates the aligner and coverageReport folder
    `mkdir -p $script_parameters{'outScriptDir'}/$_[0]/$_[1];`; #Creates the aligner script directory
    $filename = "$script_parameters{'outScriptDir'}/$_[0]/$_[1]/rCovPlots_$_[0].";
      
    Checkfnexists($filename, $fnend);

#Info and Logg
    print STDOUT "Creating sbatch script RcoveragePlots and writing script file(s) to: ", $filename, "\n";print MASTERL "Creating sbatch script RcoveragePlots and writing script file(s) to: ", $filename, "\n";
    print STDOUT "Sbatch script RcoveragePlots data files will be written to: ", $script_parameters{'outDataDir'},"/$_[0]/$_[1]/coverageReport", "\n";print MASTERL "Sbatch script RcoveragePlots data files will be written to: ", $script_parameters{'outDataDir'},"/$_[0]/$_[1]/coverageReport", "\n";
    
    open (RCovP, ">$filename") or die "Can't write to $filename: $!\n";
    
    print RCovP "#! /bin/bash -l", "\n";
    print RCovP "#SBATCH -A ", $script_parameters{'projectID'}, "\n";
    print RCovP "#SBATCH -n 1 ", "\n";
    print RCovP "#SBATCH -C thin", "\n";	
    print RCovP "#SBATCH -t 01:00:00", "\n"; 

    print RCovP "#SBATCH -J RCovPlots_$_[1]", $_[0], "\n";
    print RCovP "#SBATCH -e $script_parameters{'outDataDir'}/$_[0]/$_[1]/info/rCovPlots_$_[0].", $fnt ,".stderr.txt", "\n";
    print RCovP "#SBATCH -o $script_parameters{'outDataDir'}/$_[0]/$_[1]/info/rCovPlots_$_[0].", $fnt ,".stdout.txt", "\n";
  
    unless ($em eq 0) {
	
	print RCovP "#SBATCH --mail-type=END", "\n";
	print RCovP "#SBATCH --mail-type=FAIL", "\n";
	print RCovP "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print RCovP 'echo "Running on: $(hostname)"',"\n\n";
    print RCovP "module load bioinfo-tools", "\n\n"; 
    print RCovP "module load R/2.12.2", "\n\n";
    print RCovP "#Samples", "\n";
    
    print RCovP 'inSampleDir="',"$script_parameters{'outDataDir'}/$_[0]/$_[1]/coverageReport", '"', "\n";
    print RCovP 'outSampleDir="', "$script_parameters{'outDataDir'}/$_[0]/$_[1]/coverageReport", '"', "\n\n"; 
 
    my $fileending;
    if ($script_parameters{'pPicardToolsMarkduplicates'} eq 0) { #Ensure correct file ending
	$fileending = "_sorted";
    }
    else {
	$fileending = "_sorted_pmd";
    }	
    for (my $infile=0;$infile<scalar( @{ $Infiles_lane_noending{$_[0]} } );$infile++) { #For all files from MosaikBuild
	
	my $tempinfile = $Infiles_lane_noending{$_[0]}[$infile];
	
	if ($script_parameters{'pGenomeCoverageBED'} eq 1) {
	    print RCovP "Rscript $script_parameters{'inScriptDir'}/covplots_genome.R  ", '${inSampleDir}', "/$tempinfile", $fileending, "_genomeCoverageBed $tempinfile $xcoverage ", '${outSampleDir}', "\n\n";
	}
	if ($script_parameters{'pCoverageBED'} eq 1) {
	    #my $file = '${inSampleDir}', "/$tempinfile", $fileending, "_coverageBed_hist";
	    print RCovP "grep ^all ",'${inSampleDir}', "/$tempinfile", $fileending, "_coverageBed_hist > ",'${inSampleDir}', "/$tempinfile", $fileending, "_coverageBed_all_hist", "\n\n"; #Prepp indata file to contain only all features 
	    print RCovP "Rscript $script_parameters{'inScriptDir'}/covplots_exome_all.R  ", '${inSampleDir}', "/$tempinfile", $fileending, "_coverageBed_all_hist $tempinfile ", '${outSampleDir}', "\n\n";
#Duplicates removed
	    print RCovP "grep ^all ",'${inSampleDir}', "/$tempinfile", $fileending, "_rmdup_coverageBed_hist > ",'${inSampleDir}', "/$tempinfile", $fileending, "_rmdup_coverageBed_all_hist", "\n\n"; #Prepp indata file to contain only all features 
	    print RCovP "Rscript $script_parameters{'inScriptDir'}/covplots_exome_all.R  ", '${inSampleDir}', "/$tempinfile", $fileending, "_rmdup_coverageBed_all_hist $tempinfile","_rmdup ", '${outSampleDir}', "\n\n";	    
	}
    }
    
    if ($script_parameters{'pPicardToolsMerge'} eq 1 && scalar( @{ $Infiles_lane_noending{$_[0]} }) > 1) { #but only if there is more than one mosaikBuild file
	if ($script_parameters{'pPicardToolsMarkduplicates'} eq 0) { #Ensure correct file ending
	    $fileending = "_sorted_merged";
	}
	else {
	    $fileending = "_sorted_merged_pmd";
	}
	if ($script_parameters{'pGenomeCoverageBED'} eq 1) {
	print RCovP "Rscript $script_parameters{'inScriptDir'}/covplots_genome.R  ", '${inSampleDir}',"/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, $fileending,"_genomeCoverageBed $_[0]_lanes_", @{ $lanes{$_[0]} }," $xcoverage ", '${outSampleDir}', "\n\n";
	}
	if ($script_parameters{'pCoverageBED'} eq 1) {
	    #my $file = '${inSampleDir}',"/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, $fileending, "_coverageBed_hist";	    
	    print RCovP "grep ^all ", '${inSampleDir}',"/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, $fileending, "_coverageBed_hist > ",'${inSampleDir}',"/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, $fileending, "_coverageBed_all_hist", "\n\n"; #Prepp indata file to contain only all features 
	    print RCovP "Rscript $script_parameters{'inScriptDir'}/covplots_exome_all.R  ", '${inSampleDir}',"/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, $fileending, "_coverageBed_all_hist $_[0]_lanes_", @{ $lanes{$_[0]} }, '${outSampleDir}', "\n\n";
	    #Duplicates removed
	    print RCovP "grep ^all ", '${inSampleDir}',"/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, $fileending, "_rmdup_coverageBed_hist > ",'${inSampleDir}',"/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, $fileending, "_rmdup_coverageBed_all_hist", "\n\n"; #Prepp indata file to contain only all features 
	    print RCovP "Rscript $script_parameters{'inScriptDir'}/covplots_exome_all.R  ", '${inSampleDir}',"/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, $fileending, "_rmdup_coverageBed_all_hist $_[0]_lanes_", @{ $lanes{$_[0]} }, "_rmdup ", '${outSampleDir}', "\n\n";    
	}
    }
    if (@picardToolsMergePreviousFiles ) { # Coverage report R plots on files merged this round with merged file from previous round
	if ( ($script_parameters{'pPicardToolsMerge'} == 0) && ($script_parameters{'pPicardToolsMarkduplicates'} == 0) ) { #Ensure correct file ending
	    $fileending = "_sorted";
	}
	if ( ($script_parameters{'pPicardToolsMerge'} == 1) && ($script_parameters{'pPicardToolsMarkduplicates'} eq 0 ) ) {
	    $fileending = "_sorted_merged";
	}
	if ( ($script_parameters{'pPicardToolsMerge'} == 0) && ($script_parameters{'pPicardToolsMarkduplicates'} == 1) ) {
	    $fileending = "_sorted_pmd";
	}
	if ( ($script_parameters{'pPicardToolsMerge'} == 1) && ($script_parameters{'pPicardToolsMarkduplicates'} == 1) ) {
	    $fileending = "_sorted_merged_pmd";
	}
	for (my $mergefile=0;$mergefile<scalar(@picardToolsMergePreviousFiles);$mergefile++) {
	    
	    if ($picardToolsMergePreviousFiles[$mergefile] =~ /$_[0]/) { #Look for sampleID in previously generated file to be merged with current run to be able to merge correct files
		if ($picardToolsMergePreviousFiles[$mergefile] =~ /lane(\d+)|s_(\d+)/) { #Look for lanes_ or lane\d in previously generated file to be merged with current run to be able to extract previous lanes

		    my $mergelanes; if($1) {$mergelanes = $1;} else {$mergelanes = $2;} #Make sure to always supply lanes from previous regexp

		    if ($script_parameters{'pGenomeCoverageBED'} eq 1) {
			print RCovP "Rscript $script_parameters{'inScriptDir'}/covplots_genome.R ", '${inSampleDir}',"/$_[0]_", "lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending, "_genomeCoverageBed $_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }," $xcoverage ", '${outSampleDir}', "\n\n";
		    }
		    if ($script_parameters{'pCoverageBED'} eq 1) {
			print RCovP "grep ^all ", '${inSampleDir}',"/$_[0]_", "lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending, "_coverageBed_hist > ", '${inSampleDir}',"/$_[0]_", "lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending, "_coverageBed_all_hist", "\n\n"; #Prepp indata file to contain only all features
			print RCovP "Rscript $script_parameters{'inScriptDir'}/covplots_exome_all.R  ", '${inSampleDir}',"/$_[0]_", "lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending, "_coverageBed_all_hist $_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, " ", '${outSampleDir}', "\n\n";
			#Duplicates removed
			print RCovP "grep ^all ", '${inSampleDir}',"/$_[0]_", "lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending, "_rmdup_coverageBed_hist > ", '${inSampleDir}',"/$_[0]_", "lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending, "_rmdup_coverageBed_all_hist", "\n\n"; #Prepp indata file to contain only all features
			print RCovP "Rscript $script_parameters{'inScriptDir'}/covplots_exome_all.R  ", '${inSampleDir}',"/$_[0]_", "lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending, "rmdup_coverageBed_all_hist $_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, "_rmdup ", '${outSampleDir}', "\n\n";
		    }
		}
	    }
	}
	print RCovP "wait", "\n\n";
    }
    close(RCovP);
    ParallelSampleIDSubmitJob($_[0],$filename,"all");
    return;
}

sub Cal_Coverage { 
#Generates sbatch scripts and runs calculates coverage on files generated from MosaikAlign or BWA (sorted). 
#NOTE:Collect_info.pl collects key metric reference file from .alignment_summary_metrics. If not processed genome build will be missing in key metric file.
#$_[0] = sampleid
#$_[1] = aligner
    
    `mkdir -p $script_parameters{'outDataDir'}/$_[0]/$_[1]/info;`; #Creates the aligner folder and info data file directory
    `mkdir -p $script_parameters{'outDataDir'}/$_[0]/$_[1]/coverageReport;`; #Creates the aligner and coverageReport folder
    `mkdir -p $script_parameters{'outDataDir'}/$_[2]/$_[1]/GATK/candidates/coverageReport;`; #Creates the aligner and coverageReport folder. NOTE: FamilyID
    `mkdir -p $script_parameters{'outScriptDir'}/$_[0]/$_[1];`; #Creates the aligner script directory
    $filename = "$script_parameters{'outScriptDir'}/$_[0]/$_[1]/cal_coverage_$_[0].";

    Checkfnexists($filename, $fnend);

#Info and Logg
    print STDOUT "Creating sbatch script Calculate Coverage and writing script file(s) to: ", $filename, "\n";print MASTERL "Creating sbatch script Calculate Coverage and writing script file(s) to: ", $filename, "\n";
    print STDOUT "Sbatch script Calculate Coverage data files will be written to: ", $script_parameters{'outDataDir'},"/$_[0]/$_[1]/coverageReport", "\n";print MASTERL "Sbatch script Calculate Coverage data files will be written to: ", $script_parameters{'outDataDir'},"/$_[0]/$_[1]/coverageReport", "\n";

    my $t = ceil(3*scalar( @{ $Infiles_bothstrands_noending{$_[0]} })); #One full lane on Hiseq takes approx. 2 h to process, round up to nearest full hour.
    
    open (CRG, ">$filename") or die "Can't write to $filename: $!\n";
    
    print CRG "#! /bin/bash -l", "\n";
    print CRG "#SBATCH -A ", $script_parameters{'projectID'}, "\n";
    print CRG "#SBATCH -p node -n $maximumCores ", "\n";
    print CRG "#SBATCH -C thin", "\n";
    if ($script_parameters{'pPicardToolsMerge'} eq 0) {	
	print CRG "#SBATCH -t 4:00:00", "\n";	
    }
    else{
	print CRG "#SBATCH -t $t:00:00", "\n";	
    }	

    print CRG "#SBATCH -J CRG_$_[1]", $_[0], "\n";
    print CRG "#SBATCH -e $script_parameters{'outDataDir'}/$_[0]/$_[1]/info/cal_coverage_$_[0].", $fnt ,".stderr.txt", "\n";
    print CRG "#SBATCH -o $script_parameters{'outDataDir'}/$_[0]/$_[1]/info/cal_coverage_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print CRG "#SBATCH --mail-type=END", "\n";
	print CRG "#SBATCH --mail-type=FAIL", "\n";
	print CRG "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print CRG 'echo "Running on: $(hostname)"',"\n\n";
    print CRG "module load bioinfo-tools", "\n\n"; 
    print CRG "#Samples", "\n";
    
    print CRG 'inSampleDir="',"$script_parameters{'outDataDir'}/$_[0]/$_[1]", '"', "\n";
    print CRG 'outSampleDir="', "$script_parameters{'outDataDir'}/$_[0]/$_[1]/coverageReport", '"', "\n\n"; 
    print CRG 'outFamilyDir="', "$script_parameters{'outDataDir'}/$_[2]/$_[1]/GATK/candidates/coverageReport", '"', "\n\n";
   
    my $core_Counter=1;
    my $fileending;
    if ($script_parameters{'pPicardToolsMarkduplicates'} eq 0) { #Ensure correct file ending
	$fileending = "_sorted";
    }
    else {
	$fileending = "_sorted_pmd";
    }
    for (my $infile=0;$infile<scalar( @{ $Infiles_lane_noending{$_[0]} } );$infile++) { #For all files from MosaikAlign or BWA_sampe
	
	my $tempinfile = $Infiles_lane_noending{$_[0]}[$infile];
	if ($script_parameters{'pGenomeCoverageBED'} eq 1) {
	    print CRG "genomeCoverageBed -max $xcoverage -ibam ", '${inSampleDir}', "/$tempinfile",$fileending,".bam > ", '${outSampleDir}',"/$tempinfile", $fileending,"_genomeCoverageBed &", "\n\n";
	}
	if ($script_parameters{'pQaCompute'} eq 1) {
	    print CRG "qaCompute -m -d -i -c $xcoverage ", '${inSampleDir}', "/$tempinfile",$fileending,".bam ", '${outSampleDir}',"/$tempinfile", $fileending,"_qaCompute &", "\n\n";
	}
	if ($script_parameters{'pPicardToolsCollectMultipleMetrics'} eq 1) {
	    print CRG "java -Xmx4g -jar $picardPath/CollectMultipleMetrics.jar I=", '${inSampleDir}', "/$tempinfile",$fileending,".bam O=", '${outSampleDir}',"/$tempinfile", $fileending," R=$script_parameters{'referencesDir'}/$humanGenomeReference &", "\n\n";
	}
	if ($script_parameters{'pPicardToolsCalculateHSMetrics'} eq 1) { #Run CalculateHsMetrics (exome)
	    print CRG "java -Xmx4g -jar $picardPath/CalculateHsMetrics.jar INPUT=", '${inSampleDir}', "/$tempinfile",$fileending,".bam OUTPUT=", '${outSampleDir}',"/$tempinfile", $fileending,"_CalculateHsMetrics REFERENCE_SEQUENCE=$script_parameters{'referencesDir'}/$humanGenomeReference BAIT_INTERVALS=$script_parameters{'referencesDir'}/",$sample_info{$script_parameters{'familyID'}}{$_[0]}{'exomeTargetPaddedBedInfileList'}," TARGET_INTERVALS=$script_parameters{'referencesDir'}/",$sample_info{$script_parameters{'familyID'}}{$_[0]}{'exomeTargetBedInfileList'}," &", "\n\n";
	}
	if ($script_parameters{'pCoverageBED'} == 1) { #Run coverageBed (exome)
	    print CRG "coverageBed -hist -abam ", '${inSampleDir}', "/$tempinfile",$fileending,".bam -b $script_parameters{'referencesDir'}/",$sample_info{$script_parameters{'familyID'}}{$_[0]}{'exomeTargetBed'}," > ", '${outSampleDir}',"/$tempinfile", $fileending,"_coverageBed_hist &", "\n\n";
	    #Remove PCR and Optical duplicates
	    print CRG "samtools view -F 0x400 -b ", '${inSampleDir}', "/$tempinfile",$fileending,".bam | coverageBed -hist -abam stdin -b $script_parameters{'referencesDir'}/",$sample_info{$script_parameters{'familyID'}}{$_[0]}{'exomeTargetBed'}," > ", '${outSampleDir}',"/$tempinfile", $fileending,"_rmdup_coverageBed_hist &", "\n\n";
	    
	}
	print CRG "wait", "\n\n";
    }

    if ( $script_parameters{'pPicardToolsMerge'} eq 1 && scalar( @{ $Infiles_lane_noending{$_[0]} }) > 1 ) { #but only if there is more than one mosaikBuild/BWA_aln file per sample ID (Sanity check)
	if ($script_parameters{'pPicardToolsMarkduplicates'} eq 0) { #Ensure correct file ending
	    $fileending = "_sorted_merged";
	}
	else {
	    $fileending = "_sorted_merged_pmd";
	}	
	if ($script_parameters{'pGenomeCoverageBED'} eq 1) {
	    print CRG "genomeCoverageBed -max $xcoverage -ibam ", '${inSampleDir}',"/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, $fileending,".bam > ", '${outSampleDir}', "/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, $fileending,"_genomeCoverageBed &", "\n\n";
	}
	if ($script_parameters{'pQaCompute'} eq 1) {
	    print CRG "qaCompute -m -d -i -c $xcoverage ", '${inSampleDir}',"/$_[0]_lanes_", @{ $lanes{$_[0]} }, $fileending,".bam ", '${outSampleDir}', "/$_[0]_lanes_", @{ $lanes{$_[0]} }, $fileending,"_qaCompute &", "\n\n";
	}
	if ($script_parameters{'pPicardToolsCollectMultipleMetrics'} eq 1) {
	    print CRG "java -Xmx4g -jar $picardPath/CollectMultipleMetrics.jar I=", '${inSampleDir}',"/$_[0]_lanes_", @{ $lanes{$_[0]} }, $fileending,".bam O=", '${outSampleDir}', "/$_[0]_lanes_", @{ $lanes{$_[0]} }, $fileending," R=$script_parameters{'referencesDir'}/$humanGenomeReference &", "\n\n";
	}
	if ($script_parameters{'pPicardToolsCalculateHSMetrics'} eq 1) { #Run CalculateHsMetrics
	    print CRG "java -Xmx4g -jar $picardPath/CalculateHsMetrics.jar INPUT=", '${inSampleDir}',"/$_[0]_lanes_", @{ $lanes{$_[0]} }, $fileending,".bam OUTPUT=", '${outSampleDir}', "/$_[0]_lanes_", @{ $lanes{$_[0]} }, $fileending,"_CalculateHsMetrics REFERENCE_SEQUENCE=$script_parameters{'referencesDir'}/$humanGenomeReference BAIT_INTERVALS=$script_parameters{'referencesDir'}/",$sample_info{$script_parameters{'familyID'}}{$_[0]}{'exomeTargetPaddedBedInfileList'}," TARGET_INTERVALS=$script_parameters{'referencesDir'}/",$sample_info{$script_parameters{'familyID'}}{$_[0]}{'exomeTargetBedInfileList'}," &", "\n\n";
	}
	if ($script_parameters{'pCoverageBED'} == 1) { #Run coverageBed (exome)
	    
	    print CRG "coverageBed -hist -abam ", '${inSampleDir}',"/$_[0]_lanes_", @{ $lanes{$_[0]} }, $fileending,".bam -b $script_parameters{'referencesDir'}/",$sample_info{$script_parameters{'familyID'}}{$_[0]}{'exomeTargetBed'}," > ", '${outSampleDir}', "/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, $fileending,"_coverageBed_hist &", "\n\n";
	    #Remove PCR and Optical duplicates
	    print CRG "samtools view -F 0x400 -b ", '${inSampleDir}',"/$_[0]_lanes_", @{ $lanes{$_[0]} }, $fileending,".bam | coverageBed -hist -abam stdin -b $script_parameters{'referencesDir'}/",$sample_info{$script_parameters{'familyID'}}{$_[0]}{'exomeTargetBed'}," > ", '${outSampleDir}', "/$_[0]_lanes_", @{ $lanes{$_[0]} }, $fileending,"_rmdup_coverageBed_hist &", "\n\n";
	}
	print CRG "wait", "\n\n"; 
    }
    
    if (@picardToolsMergePreviousFiles ) { # Coverage report on files merged this round with merged file from previous round
	if ( ($script_parameters{'pPicardToolsMerge'} == 0) && ($script_parameters{'pPicardToolsMarkduplicates'} == 0) ) { #Ensure correct file ending
	    $fileending = "_sorted";
	}
	if ( ($script_parameters{'pPicardToolsMerge'} == 1) && ($script_parameters{'pPicardToolsMarkduplicates'} eq 0 ) ) {
	    $fileending = "_sorted_merged";
	}
	if ( ($script_parameters{'pPicardToolsMerge'} == 0) && ($script_parameters{'pPicardToolsMarkduplicates'} == 1) ) {
	    $fileending = "_sorted_pmd";
	}
	if ( ($script_parameters{'pPicardToolsMerge'} == 1) && ($script_parameters{'pPicardToolsMarkduplicates'} == 1) ) {
	    $fileending = "_sorted_merged_pmd";
	}
	print CRG "wait", "\n\n"; #For previous runs
	for (my $mergefile=0;$mergefile<scalar(@picardToolsMergePreviousFiles);$mergefile++) {
	    
	    if ($picardToolsMergePreviousFiles[$mergefile] =~ /$_[0]/) { #Look for sampleID in previously generated file to be merged with current run to be able to merge correct files
		if ($picardToolsMergePreviousFiles[$mergefile] =~ /lane(\d+)|s_(\d+)/) { #Look for lanes_ or lane\d in previously generated file to be merged with current run to be able to extract previous lanes

		    my $mergelanes; if($1) {$mergelanes = $1;} else {$mergelanes = $2;} #Make sure to always supply lanes from previous regexp

		    if ($script_parameters{'pGenomeCoverageBED'} eq 1) {
			print CRG "genomeCoverageBed -max $xcoverage -ibam ", '${inSampleDir}',"/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending,".bam > ", '${outSampleDir}', "/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending,"_genomeCoverageBed &", "\n\n";
		    }
		    if ($script_parameters{'pQaCompute'} eq 1) {
			print CRG "qaCompute -m -d -i -c $xcoverage ", '${inSampleDir}',"/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending,".bam ", '${outSampleDir}', "/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending,"_qaCompute &", "\n\n";
		    }
		    if ($script_parameters{'pPicardToolsCollectMultipleMetrics'} eq 1) {
			print CRG "java -Xmx4g -jar $picardPath/CollectMultipleMetrics.jar I=", '${inSampleDir}',"/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending,".bam O=", '${outSampleDir}', "/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending," R=$script_parameters{'referencesDir'}/$humanGenomeReference &", "\n\n";
		    }
		    if ($script_parameters{'pPicardToolsCalculateHSMetrics'} eq 1) { #Run CalculateHsMetrics
			print CRG "java -Xmx4g -jar $picardPath/CalculateHsMetrics.jar INPUT=", '${inSampleDir}',"/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending,".bam OUTPUT=", '${outSampleDir}', "/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending,"_CalculateHsMetrics REFERENCE_SEQUENCE=$script_parameters{'referencesDir'}/$humanGenomeReference BAIT_INTERVALS=$script_parameters{'referencesDir'}/",$sample_info{$script_parameters{'familyID'}}{$_[0]}{'exomeTargetPaddedBedInfileList'}," TARGET_INTERVALS=$script_parameters{'referencesDir'}/",$sample_info{$script_parameters{'familyID'}}{$_[0]}{'exomeTargetBedInfileList'}," &", "\n\n";
		    }
		    if ($script_parameters{'pCoverageBED'} == 1) { #Run coverageBed (exome)
			
			print CRG "coverageBed -hist -abam ", '${inSampleDir}',"/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending,".bam -b $script_parameters{'referencesDir'}/",$sample_info{$script_parameters{'familyID'}}{$_[0]}{'exomeTargetBed'}," > ", '${outSampleDir}', "/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending,"_coverageBed_hist &", "\n\n";
			#Remove PCR and Optical duplicates
			print CRG "samtools view -F 0x400 -b ", '${inSampleDir}',"/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending,".bam | coverageBed -hist -abam stdin -b $script_parameters{'referencesDir'}/",$sample_info{$script_parameters{'familyID'}}{$_[0]}{'exomeTargetBed'}," > ", '${outSampleDir}', "/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending,"_rmdup_coverageBed_hist &", "\n\n";
		    }
		}
	    }
	}
    }
    print CRG "wait", "\n\n";

###
#Coverage Report
###    
#Only 1 set of files should be generated with the BAM infile that contains the most information.
    if (@picardToolsMergePreviousFiles ) { # Coverage report on files merged this round with merged file from previous round
	if ( ($script_parameters{'pPicardToolsMerge'} == 0) && ($script_parameters{'pPicardToolsMarkduplicates'} == 0) ) { #Ensure correct file ending
	    $fileending = "_sorted";
	}
	if ( ($script_parameters{'pPicardToolsMerge'} == 1) && ($script_parameters{'pPicardToolsMarkduplicates'} eq 0 ) ) {
	    $fileending = "_sorted_merged";
	}
	if ( ($script_parameters{'pPicardToolsMerge'} == 0) && ($script_parameters{'pPicardToolsMarkduplicates'} == 1) ) {
	    $fileending = "_sorted_pmd";
	}
	if ( ($script_parameters{'pPicardToolsMerge'} == 1) && ($script_parameters{'pPicardToolsMarkduplicates'} == 1) ) {
	    $fileending = "_sorted_merged_pmd";
	}
	for (my $mergefile=0;$mergefile<scalar(@picardToolsMergePreviousFiles);$mergefile++) {
	    
	    if ($picardToolsMergePreviousFiles[$mergefile] =~ /$_[0]/) { #Look for sampleID in previously generated file to be merged with current run to be able to merge correct files
		if ($picardToolsMergePreviousFiles[$mergefile] =~ /lane(\d+)|s_(\d+)/) { #Look for lanes_ or lane\d in previously generated file to be merged with current run to be able to extract previous lanes
		    
		    my $mergelanes; if($1) {$mergelanes = $1;} else {$mergelanes = $2;} #Make sure to always supply lanes from previous regexp
		    
		    print CRG "coverageBed -abam ", '${inSampleDir}',"/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending,".bam -b $script_parameters{'referencesDir'}/",$sample_info{$script_parameters{'familyID'}}{$_[0]}{'exomeTargetBed'}," > ", '${outSampleDir}', "/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending,"_coverageBed &", "\n\n";
		    print CRG "coverageBed -d -abam ", '${inSampleDir}',"/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending,".bam -b $script_parameters{'referencesDir'}/",$sample_info{$script_parameters{'familyID'}}{$_[0]}{'exomeTargetBed'}," > ", '${outSampleDir}', "/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending,"_coverageBed_depth_pos &", "\n\n";
		    
		    print CRG "wait", "\n\n";
		    
                    #Generate the average depth per feature and the number of zero-coverage bases as well as their relative to feature position
		    print CRG "#Generate the average depth per feature and the number of zero-coverage bases as well as their relative to feature position", "\n";
		    print CRG q?perl -nae'my $prev_chr=0;my $prev_start;my $prev_end;my $average_cov=0;my $nr_positions=0;my %feature;my %zerofeature;my @temp_line; while (<>) { chomp($_);  @temp_line = split("\t",$_); if($feature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]}) { if ($temp_line[5] >= 10) { $feature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]}{'high-coverage'}++;} if ($temp_line[5] == 0) { push (@ {$zerofeature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]}{'zero-coverage'}}, $temp_line[1]+$temp_line[4]-1);} $average_cov=$average_cov+$temp_line[5];$nr_positions++;} else { if ( $prev_chr && ($prev_start != $temp_line[1]) && ($prev_end != $temp_line[2]) ) { print "$prev_chr\t","$prev_start\t","$prev_end\t",$average_cov/$nr_positions,"\t"; if ($feature{$prev_chr}{$prev_start}{$prev_end}{'high-coverage'}) { print $feature{$prev_chr}{$prev_start}{$prev_end}{'high-coverage'}/$nr_positions, "\t";} else {print "0\t";} if ($zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}) { print scalar(@ {$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'} }), "\t"; for (my $i=0;$i<scalar(@ {$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}});$i++) { print "$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}[$i];";} print "\n"; } else {print "0\t","Na\n";} } %feature = ();%zerofeature = ();$average_cov=0;$nr_positions=0;$feature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]} = $_; $prev_chr = $temp_line[0];$prev_start = $temp_line[1];$prev_end = $temp_line[2];if ($temp_line[5] >= 10) { $feature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]}{'high-coverage'}++;} if ($temp_line[5] == 0) { push (@ {$zerofeature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]}{'zero-coverage'} }, $temp_line[1]+$temp_line[4]-1);} $average_cov=$average_cov+$temp_line[5];$nr_positions++;} } print "$prev_chr\t","$prev_start\t","$prev_end\t",$average_cov/$nr_positions,"\t"; if ($feature{$prev_chr}{$prev_start}{$prev_end}{'high-coverage'}) { print $feature{$prev_chr}{$prev_start}{$prev_end}{'high-coverage'}/$nr_positions, "\t";} else {print "0\t";} if ($zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}) { print scalar(@ {$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'} }), "\t"; for (my $i=0;$i<scalar(@ {$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}});$i++) { print "$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}[$i];";} print "\n"; } else {print "0\t","Na\n";} last;' ?, '${outSampleDir}', "/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending,"_coverageBed_depth_pos > ", '${outSampleDir}', "/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending,"_coverageBed_depth_pos_collapsed", "\n\n";
		    
		    #Concatenate to 1 file to be able to include info from _coverageBed file
		    print CRG "#Concatenate to 1 file to be able to include info from _coverageBed file", "\n";
		    print CRG "cat ", '${outSampleDir}', "/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending,"_coverageBed_depth_pos_collapsed ", '${outSampleDir}', "/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending,"_coverageBed > ", '${outSampleDir}', "/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending,"_coverageBed_merged_temp", "\n\n";
		    
		    #Sort on chr and then numerically on start position. NOTE makes the two lines from each file end up next to each other. The line order is preserved unless the feature annotation is "-", then it is reversed (exception handled by the next perl one-liner) 
		    print CRG "sort -k1,1 -k2,2n ", '${outSampleDir}', "/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending,"_coverageBed_merged_temp > ", '${outSampleDir}', "/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending,"_coverageBed_merged", "\n\n";
		    
		    #Collpase the two lines from coverageBed_depth_pos_collapsed and _coverageBed
		    print CRG q?perl -nae ' print "#Chr\tStart\tStop\tOverlapping_Reads\tNon-zero_Coverage_Bases\tFeature_Length(Nucleotides)\tAverage_Coverage\tFraction_Non-zero_Coverage_Bases\tFraction_Ten_Coverage_Bases\tNr_Zero_Coverage_Bases\tRelative_position_for_Zero_Coverage_Bases\n";chomp($_);my @prev_line=split("\t",$_);my @temp_line; while (<>) { chomp($_); @temp_line = split("\t",$_); if ($prev_line[0]  && ($prev_line[0] eq $temp_line[0]) && ($prev_line[1] == $temp_line[1]) && ($prev_line[2] == $temp_line[2]) ) { if ($prev_line[3] eq "-") { $temp_line[3]=~s/\./\,/;$prev_line[7]=~s/\./\,/; print $prev_line[0],"\t", $prev_line[1],"\t",$prev_line[2],"\t",$prev_line[4],"\t",$prev_line[5],"\t",$prev_line[6],"\t",$temp_line[3],"\t","$prev_line[7]\t","$temp_line[4]\t","$temp_line[5]\t","$temp_line[6]\n";} else { $prev_line[3]=~s/\./\,/;$temp_line[7]=~s/\./\,/; print $temp_line[0],"\t", $temp_line[1],"\t",$temp_line[2],"\t",$temp_line[4],"\t",$temp_line[5],"\t",$temp_line[6],"\t",$prev_line[3],"\t","$temp_line[7]\t","$prev_line[4]\t","$prev_line[5]\t","$prev_line[6]\n";} } else { @prev_line = @temp_line;} }last;' ?, '${outSampleDir}', "/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending,"_coverageBed_merged > ", '${outSampleDir}',"/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending,"_target_coverage.txt", "\n\n";
		    
		   
		    #Add chr to entry to enable comparison to Im_Db
		    print CRG "#Add chr to entry to enable later comparison to Im_Db", "\n";
		    print CRG q?perl -i -p -e ' if($_=~/^#/) {} else {s/^(.+)/chr$1/g }' ?, '${outSampleDir}',"/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending,"_target_coverage.txt", "\n\n";
		    
		    #Removal of files which the necessary info has been extracted from
		    print CRG "rm ", '${outSampleDir}', "/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending,"_coverageBed", "\n\n";
		    print CRG "rm ", '${outSampleDir}', "/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending,"_coverageBed_depth_pos", "\n\n";
		    print CRG "rm ", '${outSampleDir}', "/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending,"_coverageBed_depth_pos_collapsed", "\n\n";
		    print CRG "rm ", '${outSampleDir}', "/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending,"_coverageBed_merged_temp", "\n\n";
		    print CRG "rm ", '${outSampleDir}', "/$_[0]_lanes_",$mergelanes, @{ $lanes{$_[0]} }, $fileending,"_coverageBed_merged", "\n\n";
		    
		    my $target_coverage_db_file = "$script_parameters{'outDataDir'}/$_[0]/$_[1]/coverageReport/$_[0]_lanes_".$mergelanes;
		    for (my $i=0;$i<scalar(@ { $lanes{$_[0]} });$i++) {
			$target_coverage_db_file .= $lanes{$_[0]}[$i];
		    }
		    my $target_coverage_file =  $target_coverage_db_file;

		    $target_coverage_db_file .= $fileending."_coverage_target_db_master.txt";		    
		    $target_coverage_file .= $fileending."_target_coverage.txt";    
		    my @target_coverage_db_files = ($target_coverage_db_file); #Db master files	    
		    my @target_coverage_files = ($target_coverage_file); #Target coverage files created above
		
###
#Create target file with EnsembleGeneID (Required for coverage calculation in rank_script)
###
                    #Create db master template 
		    print CRG "#Add GeneName to coverage files", "\n";
		    for (my $db_fileCounter=0;$db_fileCounter<scalar(@target_coverage_db_files);$db_fileCounter++) {
			open (TARCOV, ">$target_coverage_db_files[$db_fileCounter]") or die "Can't write to $target_coverage_db_files[$db_fileCounter]: $!\n";
			print TARCOV "outcolumns=0_0,0_1,0_2,0_3,0_4,0_5,0_6,0_7,0_8,0_9,0_10,1_4\n";
			print TARCOV "outheaders=Chr,Start,Stop,Overlapping_Reads,Non-zero_Coverage_Bases,Feature_Length(Nucleotides),Average_Coverage,Fraction_Non-zero_Coverage_Bases,Fraction_Ten_Coverage_Bases,Nr_Zero_Coverage_Bases,Relative_position_for_Zero_Coverage_Bases,Ensembl_GeneID\n";
			print TARCOV $target_coverage_files[$db_fileCounter],"\t",'\t',"\t0,1,2\t0\texact\t0,1,2,3,4,5,6,7,8,9,10\tsmall","\n";
			print TARCOV "$script_parameters{'referencesDir'}/mart_export_Ensembl_GeneID_key_cleaned_noblanks_chr.txt\t",'\t',"\t0,1,2\t0\trange\t4\tsmall", "\n";
			close(TARCOV);
			
                        #Add GeneNameID
			print CRG "perl $script_parameters{'inScriptDir'}/intersectCollect.pl -o $target_coverage_files[$db_fileCounter] -db $target_coverage_db_files[$db_fileCounter] -prechr 1", "\n\n";
		    }
		    
		}
	    }
	}
    }
    
    elsif ( $script_parameters{'pPicardToolsMerge'} eq 1 && scalar( @{ $Infiles_lane_noending{$_[0]} }) > 1 ) { #but only if there is more than one mosaikBuild/BWA_aln file per sample ID (Sanity check)
	if ($script_parameters{'pPicardToolsMarkduplicates'} eq 0) { #Ensure correct file ending
	    $fileending = "_sorted_merged";
	}
	else {
	    $fileending = "_sorted_merged_pmd";
	}

	print CRG "coverageBed -abam ", '${inSampleDir}',"/$_[0]_lanes_", @{ $lanes{$_[0]} }, $fileending,".bam -b $script_parameters{'referencesDir'}/",$sample_info{$script_parameters{'familyID'}}{$_[0]}{'exomeTargetBed'}," > ", '${outSampleDir}', "/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, $fileending,"_coverageBed &", "\n\n";
	print CRG "coverageBed -d -abam ", '${inSampleDir}',"/$_[0]_lanes_", @{ $lanes{$_[0]} }, $fileending,".bam -b $script_parameters{'referencesDir'}/",$sample_info{$script_parameters{'familyID'}}{$_[0]}{'exomeTargetBed'}," > ", '${outSampleDir}', "/$_[0]_lanes_", @{ $lanes{$_[0]} }, $fileending,"_coverageBed_depth_pos &", "\n\n";

	print CRG "wait", "\n\n";

	#Generate the average depth per feature and the number of zero-coverage bases as well as their relative to feature position
	print CRG "#Generate the average depth per feature and the number of zero-coverage bases as well as their relative to feature position", "\n";
	print CRG q?perl -nae'my $prev_chr=0;my $prev_start;my $prev_end;my $average_cov=0;my $nr_positions=0;my %feature;my %zerofeature;my @temp_line; while (<>) { chomp($_);  @temp_line = split("\t",$_); if($feature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]}) { if ($temp_line[5] >= 10) { $feature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]}{'high-coverage'}++;} if ($temp_line[5] == 0) { push (@ {$zerofeature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]}{'zero-coverage'}}, $temp_line[1]+$temp_line[4]-1);} $average_cov=$average_cov+$temp_line[5];$nr_positions++;} else { if ( $prev_chr && ($prev_start != $temp_line[1]) && ($prev_end != $temp_line[2]) ) { print "$prev_chr\t","$prev_start\t","$prev_end\t",$average_cov/$nr_positions,"\t"; if ($feature{$prev_chr}{$prev_start}{$prev_end}{'high-coverage'}) { print $feature{$prev_chr}{$prev_start}{$prev_end}{'high-coverage'}/$nr_positions, "\t";} else {print "0\t";} if ($zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}) { print scalar(@ {$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'} }), "\t"; for (my $i=0;$i<scalar(@ {$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}});$i++) { print "$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}[$i];";} print "\n"; } else {print "0\t","Na\n";} } %feature = ();%zerofeature = ();$average_cov=0;$nr_positions=0;$feature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]} = $_; $prev_chr = $temp_line[0];$prev_start = $temp_line[1];$prev_end = $temp_line[2];if ($temp_line[5] >= 10) { $feature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]}{'high-coverage'}++;} if ($temp_line[5] == 0) { push (@ {$zerofeature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]}{'zero-coverage'} }, $temp_line[1]+$temp_line[4]-1);} $average_cov=$average_cov+$temp_line[5];$nr_positions++;} } print "$prev_chr\t","$prev_start\t","$prev_end\t",$average_cov/$nr_positions,"\t"; if ($feature{$prev_chr}{$prev_start}{$prev_end}{'high-coverage'}) { print $feature{$prev_chr}{$prev_start}{$prev_end}{'high-coverage'}/$nr_positions, "\t";} else {print "0\t";} if ($zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}) { print scalar(@ {$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'} }), "\t"; for (my $i=0;$i<scalar(@ {$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}});$i++) { print "$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}[$i];";} print "\n"; } else {print "0\t","Na\n";} last;' ?, '${outSampleDir}', "/$_[0]_lanes_", @{ $lanes{$_[0]} }, $fileending,"_coverageBed_depth_pos > ", '${outSampleDir}', "/$_[0]_lanes_", @{ $lanes{$_[0]} }, $fileending,"_coverageBed_depth_pos_collapsed", "\n\n";
	
	#Concatenate to 1 file to be able to include info from _coverageBed file
	print CRG "#Concatenate to 1 file to be able to include info from _coverageBed file", "\n";
	print CRG "cat ", '${outSampleDir}', "/$_[0]_lanes_", @{ $lanes{$_[0]} }, $fileending,"_coverageBed_depth_pos_collapsed ", '${outSampleDir}', "/$_[0]_lanes_", @{ $lanes{$_[0]} }, $fileending,"_coverageBed > ", '${outSampleDir}', "/$_[0]_lanes_", @{ $lanes{$_[0]} }, $fileending,"_coverageBed_merged_temp", "\n\n";
	
	#Sort on chr and then numerically on start position. NOTE makes the two lines from each file end up next to each other. The line order is preserved unless the feature annotation is "-", then it is reversed (exception handled by the next perl one-liner) 
	print CRG "sort -k1,1 -k2,2n ", '${outSampleDir}', "/$_[0]_lanes_", @{ $lanes{$_[0]} }, $fileending,"_coverageBed_merged_temp > ", '${outSampleDir}', "/$_[0]_lanes_", @{ $lanes{$_[0]} }, $fileending,"_coverageBed_merged", "\n\n";
	
	#Collpase the two lines from coverageBed_depth_pos_collapsed and _coverageBed
	print CRG q?perl -nae ' print "#Chr\tStart\tStop\tOverlapping_Reads\tNon-zero_Coverage_Bases\tFeature_Length(Nucleotides)\tAverage_Coverage\tFraction_Non-zero_Coverage_Bases\tFraction_Ten_Coverage_Bases\tNr_Zero_Coverage_Bases\tRelative_position_for_Zero_Coverage_Bases\n";chomp($_);my @prev_line=split("\t",$_);my @temp_line; while (<>) { chomp($_); @temp_line = split("\t",$_); if ($prev_line[0]  && ($prev_line[0] eq $temp_line[0]) && ($prev_line[1] == $temp_line[1]) && ($prev_line[2] == $temp_line[2]) ) { if ($prev_line[3] eq "-") { $temp_line[3]=~s/\./\,/;$prev_line[7]=~s/\./\,/; print $prev_line[0],"\t", $prev_line[1],"\t",$prev_line[2],"\t",$prev_line[4],"\t",$prev_line[5],"\t",$prev_line[6],"\t",$temp_line[3],"\t","$prev_line[7]\t","$temp_line[4]\t","$temp_line[5]\t","$temp_line[6]\n";} else { $prev_line[3]=~s/\./\,/;$temp_line[7]=~s/\./\,/; print $temp_line[0],"\t", $temp_line[1],"\t",$temp_line[2],"\t",$temp_line[4],"\t",$temp_line[5],"\t",$temp_line[6],"\t",$prev_line[3],"\t","$temp_line[7]\t","$prev_line[4]\t","$prev_line[5]\t","$prev_line[6]\n";} } else { @prev_line = @temp_line;} }last;' ?, '${outSampleDir}', "/$_[0]_lanes_", @{ $lanes{$_[0]} }, $fileending,"_coverageBed_merged > ", '${outSampleDir}',"/$_[0]_lanes_", @{ $lanes{$_[0]} }, $fileending,"_target_coverage.txt", "\n\n";
	
	#Add chr to entry to enable comparison to Im_Db
	print CRG "#Add chr to entry to enable later comparison to Im_Db", "\n";
	print CRG q?perl -i -p -e ' if($_=~/^#/) {} else {s/^(.+)/chr$1/g }' ?, '${outSampleDir}',"/$_[0]_lanes_", @{ $lanes{$_[0]} }, $fileending,"_target_coverage.txt", "\n\n";	   

	#Removal of files which the necessary info has been extracted from
	print CRG "rm ", '${outSampleDir}', "/$_[0]_lanes_", @{ $lanes{$_[0]} }, $fileending,"_coverageBed", "\n\n";
	print CRG "rm ", '${outSampleDir}', "/$_[0]_lanes_", @{ $lanes{$_[0]} }, $fileending,"_coverageBed_depth_pos", "\n\n";
	print CRG "rm ", '${outSampleDir}', "/$_[0]_lanes_", @{ $lanes{$_[0]} }, $fileending,"_coverageBed_depth_pos_collapsed", "\n\n";
	print CRG "rm ", '${outSampleDir}', "/$_[0]_lanes_", @{ $lanes{$_[0]} }, $fileending,"_coverageBed_merged_temp", "\n\n";
	print CRG "rm ", '${outSampleDir}', "/$_[0]_lanes_", @{ $lanes{$_[0]} }, $fileending,"_coverageBed_merged", "\n\n";

	my $target_coverage_db_file = "$script_parameters{'outDataDir'}/$_[0]/$_[1]/coverageReport/$_[0]_lanes_";
	for (my $i=0;$i<scalar(@ { $lanes{$_[0]} });$i++) {
	    $target_coverage_db_file .= $lanes{$_[0]}[$i];
	}
	my $target_coverage_file =  $target_coverage_db_file;

	$target_coverage_db_file .= $fileending."_coverage_target_db_master.txt";
	$target_coverage_file .= $fileending."_target_coverage.txt";
	my @target_coverage_db_files = ($target_coverage_db_file); #Db master files	    
	my @target_coverage_files = ($target_coverage_file); #Target coverage files created above
###
#Create target file with EnsembleGeneID (Required for coverage calculation in rank_script)
###
	
#Create db master template 
	print CRG "#Add GeneName to coverage files", "\n";
	for (my $db_fileCounter=0;$db_fileCounter<scalar(@target_coverage_db_files);$db_fileCounter++) {
	    
	    open (TARCOV, ">$target_coverage_db_files[$db_fileCounter]") or die "Can't write to $target_coverage_db_files[$db_fileCounter]: $!\n";
	    print TARCOV "outcolumns=0_0,0_1,0_2,0_3,0_4,0_5,0_6,0_7,0_8,0_9,0_10,1_4\n";
	    print TARCOV "outheaders=Chr,Start,Stop,Overlapping_Reads,Non-zero_Coverage_Bases,Feature_Length(Nucleotides),Average_Coverage,Fraction_Non-zero_Coverage_Bases,Fraction_Ten_Coverage_Bases,Nr_Zero_Coverage_Bases,Relative_position_for_Zero_Coverage_Bases,Ensembl_GeneID\n";
	    print TARCOV $target_coverage_files[$db_fileCounter],"\t",'\t',"\t0,1,2\t0\texact\t0,1,2,3,4,5,6,7,8,9,10\tsmall","\n";
	    print TARCOV "$script_parameters{'referencesDir'}/mart_export_Ensembl_GeneID_key_cleaned_noblanks_chr.txt\t",'\t',"\t0,1,2\t0\trange\t4\tsmall", "\n";
	    close(TARCOV);
	    
#Add GeneNameID
	    print CRG "perl $script_parameters{'inScriptDir'}/intersectCollect.pl -o $target_coverage_files[$db_fileCounter] -db $target_coverage_db_files[$db_fileCounter] -prechr 1", "\n\n";
	}
    }
    else {
	if ($script_parameters{'pPicardToolsMarkduplicates'} eq 0) { #Ensure correct file ending
	    $fileending = "_sorted";
	}
	else {
	    $fileending = "_sorted_pmd";
	}
	for (my $infile=0;$infile<scalar( @{ $Infiles_lane_noending{$_[0]} } );$infile++) { #For all files from MosaikAlign or BWA_sampe
	    
	    my $tempinfile = $Infiles_lane_noending{$_[0]}[$infile];

	    print CRG "#Calculate coverage statistics to enable coverage calculation in rank_script", "\n";
	    print CRG "coverageBed -abam ", '${inSampleDir}', "/$tempinfile",$fileending,".bam -b $script_parameters{'referencesDir'}/",$sample_info{$script_parameters{'familyID'}}{$_[0]}{'exomeTargetBed'}," > ", '${outSampleDir}',"/$tempinfile", $fileending,"_coverageBed &", "\n\n";
	    print CRG "coverageBed -d -abam ", '${inSampleDir}', "/$tempinfile",$fileending,".bam -b $script_parameters{'referencesDir'}/",$sample_info{$script_parameters{'familyID'}}{$_[0]}{'exomeTargetBed'}," > ", '${outSampleDir}',"/$tempinfile", $fileending,"_coverageBed_depth_pos &", "\n\n";

	    print CRG "wait", "\n\n";

#Generate the average depth per feature and the number of zero-coverage bases as well as their relative to feature position
	    print CRG "#Generate the average depth per feature and the number of zero-coverage bases as well as their relative to feature position", "\n";
	    print CRG q?perl -nae'my $prev_chr=0;my $prev_start;my $prev_end;my $average_cov=0;my $nr_positions=0;my %feature;my %zerofeature;my @temp_line; while (<>) { chomp($_);  @temp_line = split("\t",$_); if($feature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]}) { if ($temp_line[5] >= 10) { $feature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]}{'high-coverage'}++;} if ($temp_line[5] == 0) { push (@ {$zerofeature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]}{'zero-coverage'}}, $temp_line[1]+$temp_line[4]-1);} $average_cov=$average_cov+$temp_line[5];$nr_positions++;} else { if ( $prev_chr && ($prev_start != $temp_line[1]) && ($prev_end != $temp_line[2]) ) { print "$prev_chr\t","$prev_start\t","$prev_end\t",$average_cov/$nr_positions,"\t"; if ($feature{$prev_chr}{$prev_start}{$prev_end}{'high-coverage'}) { print $feature{$prev_chr}{$prev_start}{$prev_end}{'high-coverage'}/$nr_positions, "\t";} else {print "0\t";} if ($zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}) { print scalar(@ {$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'} }), "\t"; for (my $i=0;$i<scalar(@ {$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}});$i++) { print "$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}[$i];";} print "\n"; } else {print "0\t","Na\n";} } %feature = ();%zerofeature = ();$average_cov=0;$nr_positions=0;$feature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]} = $_; $prev_chr = $temp_line[0];$prev_start = $temp_line[1];$prev_end = $temp_line[2];if ($temp_line[5] >= 10) { $feature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]}{'high-coverage'}++;} if ($temp_line[5] == 0) { push (@ {$zerofeature{$temp_line[0]}{$temp_line[1]}{$temp_line[2]}{'zero-coverage'} }, $temp_line[1]+$temp_line[4]-1);} $average_cov=$average_cov+$temp_line[5];$nr_positions++;} } print "$prev_chr\t","$prev_start\t","$prev_end\t",$average_cov/$nr_positions,"\t"; if ($feature{$prev_chr}{$prev_start}{$prev_end}{'high-coverage'}) { print $feature{$prev_chr}{$prev_start}{$prev_end}{'high-coverage'}/$nr_positions, "\t";} else {print "0\t";} if ($zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}) { print scalar(@ {$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'} }), "\t"; for (my $i=0;$i<scalar(@ {$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}});$i++) { print "$zerofeature{$prev_chr}{$prev_start}{$prev_end}{'zero-coverage'}[$i];";} print "\n"; } else {print "0\t","Na\n";} last;' ${outSampleDir}/?.$tempinfile.$fileending.q?_coverageBed_depth_pos > ${outSampleDir}/?.$tempinfile.$fileending.q?_coverageBed_depth_pos_collapsed?, "\n\n";
	    
	    #Concatenate to 1 file to be able to include info from _coverageBed file
	    print CRG "#Concatenate to 1 file to be able to include info from _coverageBed file", "\n";
	    print CRG q?cat ${outSampleDir}/?.$tempinfile.$fileending.q?_coverageBed_depth_pos_collapsed ${outSampleDir}/?.$tempinfile.$fileending.q?_coverageBed > ${outSampleDir}/?.$tempinfile.$fileending.q?_coverageBed_merged_temp?,  "\n\n";
	    
            #Sort on chr and then numerically on start position. NOTE makes the two lines from each file end up next to each other. The line order is preserved unless the feature annotation is "-", then it is reversed (exception handled by the next perl one-liner) 
	    print CRG q?sort -k1,1 -k2,2n ${outSampleDir}/?.$tempinfile.$fileending.q?_coverageBed_merged_temp > ${outSampleDir}/?.$tempinfile.$fileending.q?_coverageBed_merged?,  "\n\n";
	    
	    #Collpase the two lines from coverageBed_depth_pos_collapsed and _coverageBed
	    print CRG "#Collpase the two lines from coverageBed_depth_pos_collapsed and _coverageBed", "\n";
	    print CRG q?perl -nae ' print "#Chr\tStart\tStop\tOverlapping_Reads\tNon-zero_Coverage_Bases\tFeature_Length(Nucleotides)\tAverage_Coverage\tFraction_Non-zero_Coverage_Bases\tFraction_Ten_Coverage_Bases\tNr_Zero_Coverage_Bases\tRelative_position_for_Zero_Coverage_Bases\n";chomp($_);my @prev_line=split("\t",$_);my @temp_line; while (<>) { chomp($_); @temp_line = split("\t",$_); if ($prev_line[0]  && ($prev_line[0] eq $temp_line[0]) && ($prev_line[1] == $temp_line[1]) && ($prev_line[2] == $temp_line[2]) ) { if ($prev_line[3] eq "-") { $temp_line[3]=~s/\./\,/;$prev_line[7]=~s/\./\,/; print $prev_line[0],"\t", $prev_line[1],"\t",$prev_line[2],"\t",$prev_line[4],"\t",$prev_line[5],"\t",$prev_line[6],"\t",$temp_line[3],"\t","$prev_line[7]\t","$temp_line[4]\t","$temp_line[5]\t","$temp_line[6]\n";} else { $prev_line[3]=~s/\./\,/;$temp_line[7]=~s/\./\,/; print $temp_line[0],"\t", $temp_line[1],"\t",$temp_line[2],"\t",$temp_line[4],"\t",$temp_line[5],"\t",$temp_line[6],"\t",$prev_line[3],"\t","$temp_line[7]\t","$prev_line[4]\t","$prev_line[5]\t","$prev_line[6]\n";} } else { @prev_line = @temp_line;} }last;' ${outSampleDir}/?.$tempinfile.$fileending.q?_coverageBed_merged > ${outSampleDir}/?.$tempinfile.$fileending.q?_target_coverage.txt?, "\n\n";
	    
	    #Add chr to entry to enable comparison to Im_Db
	    print CRG "#Add chr to entry to enable later comparison to Im_Db", "\n";
	    print CRG q?perl -i -p -e ' if($_=~/^#/) {} else {s/^(.+)/chr$1/g }' ${outSampleDir}/?.$tempinfile.$fileending.q?_target_coverage.txt?, "\n\n";
            #Removal of files which the necessary info has been extracted from
	    print CRG q?rm ${outSampleDir}/?.$tempinfile.$fileending.q?_coverageBed?, "\n\n";
	    print CRG q?rm ${outSampleDir}/?.$tempinfile.$fileending.q?_coverageBed_depth_pos?, "\n\n";
	    print CRG q?rm ${outSampleDir}/?.$tempinfile.$fileending.q?_coverageBed_depth_pos_collapsed?, "\n\n";
	    print CRG q?rm ${outSampleDir}/?.$tempinfile.$fileending.q?_coverageBed_merged_temp?, "\n\n";
	    print CRG q?rm ${outSampleDir}/?.$tempinfile.$fileending.q?_coverageBed_merged?, "\n\n";
	    
###
#Create target file with EnsembleGeneID (Required for coverage calculation in rank_script)
###
	    
	    my @target_coverage_db_files = ("$script_parameters{'outDataDir'}/$_[0]/$_[1]/coverageReport/$tempinfile$fileending"."_coverage_target_db_master.txt"); #Db master files	    
	    my @target_coverage_files = ("$script_parameters{'outDataDir'}/$_[0]/$_[1]/coverageReport/$tempinfile"."$fileending"."_target_coverage.txt"); 

#Create db master template 
	    print CRG "#Add GeneName to coverage files", "\n";
	    for (my $db_fileCounter=0;$db_fileCounter<scalar(@target_coverage_db_files);$db_fileCounter++) { 
		open (TARCOV, ">$target_coverage_db_files[$db_fileCounter]") or die "Can't write to $target_coverage_db_files[$db_fileCounter]: $!\n";
		print TARCOV "outcolumns=0_0,0_1,0_2,0_3,0_4,0_5,0_6,0_7,0_8,0_9,0_10,1_4\n";
		print TARCOV "outheaders=Chr,Start,Stop,Overlapping_Reads,Non-zero_Coverage_Bases,Feature_Length(Nucleotides),Average_Coverage,Fraction_Non-zero_Coverage_Bases,Fraction_Ten_Coverage_Bases,Nr_Zero_Coverage_Bases,Relative_position_for_Zero_Coverage_Bases,Ensembl_GeneID\n";
		print TARCOV $target_coverage_files[$db_fileCounter],"\t",'\t',"\t0,1,2\t0\texact\t0,1,2,3,4,5,6,7,8,9,10\tsmall","\n";
		print TARCOV "$script_parameters{'referencesDir'}/mart_export_Ensembl_GeneID_key_cleaned_noblanks_chr.txt\t",'\t',"\t0,1,2\t0\trange\t4\tsmall", "\n";
		close(TARCOV);

#Add GeneNameID
		print CRG "perl $script_parameters{'inScriptDir'}/intersectCollect.pl -o $target_coverage_files[$db_fileCounter] -db $target_coverage_db_files[$db_fileCounter] -prechr 1", "\n\n";
	    }
	}
    }
    close(CRG);
    ParallelSampleIDSubmitJob($_[0],$filename,"all");
    return;
}

sub PicardMarkDup { 
#Generates sbatch scripts and runs PicardTools MarkDuplicates on files generated from MosaikAlign or BWA (sorted, merged)
#$_[0] = sampleid
#$_[1] = aligner
    
    `mkdir -p $script_parameters{'outDataDir'}/$_[0]/$_[1]/info;`; #Creates the aligner folder and info data file directory
    `mkdir -p $script_parameters{'outScriptDir'}/$_[0]/$_[1];`; #Creates the aligner script directory
    $filename = "$script_parameters{'outScriptDir'}/$_[0]/$_[1]/picard_markdup_$_[0].";
   
    Checkfnexists($filename, $fnend);

#Info and Logg
    print STDOUT "Creating sbatch script PicardMarkDuplicates and writing script file(s) to: ", $filename, "\n";print MASTERL "Creating sbatch script PicardMarkDuplicates and writing script file(s) to: ", $filename, "\n";
    print STDOUT "Sbatch script PicardMarkDuplicates data files will be written to: ", $script_parameters{'outDataDir'},"/$_[0]/$_[1]", "\n";print MASTERL "Sbatch script PicardMarkDuplicates data files will be written to: ", $script_parameters{'outDataDir'},"/$_[0]/$_[1]", "\n";

    my $t = ceil(3*scalar( @{ $Infiles_bothstrands_noending{$_[0]} })); #One full lane on Hiseq takes approx. 3 h to process, round up to nearest full hour.
    
    open (PMDUP, ">$filename") or die "Can't write to $filename: $!\n";
    
    print PMDUP "#! /bin/bash -l", "\n";
    print PMDUP "#SBATCH -A ", $script_parameters{'projectID'}, "\n";
    print PMDUP "#SBATCH -p node -n $maximumCores ", "\n";
    print PMDUP "#SBATCH -C thin", "\n";
    if ($script_parameters{'pPicardToolsMerge'} eq 0) {	
	print PMDUP "#SBATCH -t 3:00:00", "\n";	
    }
    else{
	print PMDUP "#SBATCH -t $t:00:00", "\n";	
    }	
    
    print PMDUP "#SBATCH -J PMDUP_$_[1]", $_[0], "\n";
    print PMDUP "#SBATCH -e $script_parameters{'outDataDir'}/$_[0]/$_[1]/info/picard_markdup_$_[0].", $fnt ,".stderr.txt", "\n";
    print PMDUP "#SBATCH -o $script_parameters{'outDataDir'}/$_[0]/$_[1]/info/picard_markdup_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print PMDUP "#SBATCH --mail-type=END", "\n";
	print PMDUP "#SBATCH --mail-type=FAIL", "\n";
	print PMDUP "#SBATCH --mail-user=$em", "\n\n";
    }
    
    print PMDUP 'echo "Running on: $(hostname)"',"\n\n";
    print PMDUP "#Samples", "\n";
    print PMDUP 'inSampleDir="',"$script_parameters{'outDataDir'}/$_[0]/$_[1]", '"', "\n";
    print PMDUP 'outSampleDir="', "$script_parameters{'outDataDir'}/$_[0]/$_[1]", '"', "\n\n"; 
    
    my $mergelanes; #To pick up mergedlanes later if required Markduolicates and samtools
    my $core_Counter=1;
    for (my $infile=0;$infile<scalar( @{ $Infiles_lane_noending{$_[0]} } );$infile++) { #For all files from MosaikAlign or BWA_sampe
	
	if ($infile eq $core_Counter*$maximumCores) { #Using only $maximumCores cores
	    
	    print PMDUP "wait", "\n\n";
	    $core_Counter=$core_Counter+1;
	}
	my $tempinfile = $Infiles_lane_noending{$_[0]}[$infile];
	
	print PMDUP "java -Xmx4g -jar $picardPath/MarkDuplicates.jar ASSUME_SORTED=true REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=LENIENT INPUT=", '${inSampleDir}',"/$tempinfile","_sorted.bam OUTPUT=", '${outSampleDir}', "/$tempinfile","_sorted_pmd.bam METRICS_FILE=", '${outSampleDir}', "/$tempinfile","_sorted_pmdmetric &","\n\n"; 
    }
    print PMDUP "wait", "\n\n";
    
    if ( $script_parameters{'pPicardToolsMerge'} eq 1 && scalar( @{ $Infiles_lane_noending{$_[0]} }) > 1 ) { #but only if there is more than one mosaikBuild/BWA_aln file per sample ID (Sanity check)
	
	print PMDUP "java -Xmx4g -jar $picardPath/MarkDuplicates.jar ASSUME_SORTED=true REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=LENIENT INPUT=",, '${inSampleDir}',"/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, "_sorted_merged.bam OUTPUT=", '${outSampleDir}', "/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, "_sorted_merged_pmd.bam METRICS_FILE=", '${outSampleDir}', "/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, "_sorted_merged_pmdmetric &", "\n\n"; 
    }
    
    if (@picardToolsMergePreviousFiles ) { # Coverage report on files merged this round with merged file from previous round
	
	for (my $mergefile=0;$mergefile<scalar(@picardToolsMergePreviousFiles);$mergefile++) {
	    
	    if ($picardToolsMergePreviousFiles[$mergefile] =~ /$_[0]/) { #Look for sampleID in previously generated file to be merged with current run to be able to merge correct files
		if ($picardToolsMergePreviousFiles[$mergefile] =~ /lane(\d+)|s_(\d+)/) { #Look for lanes_ or lane\d in previously generated file to be merged with current run to be able to extract previous lanes
		  
		    if($1) {$mergelanes = $1;} else {$mergelanes = $2;} #Make sure to always supply lanes from previous regexp  
		    print PMDUP "java -Xmx4g -jar $picardPath/MarkDuplicates.jar ASSUME_SORTED=true REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=LENIENT INPUT=", '${inSampleDir}',"/$_[0]_", "lanes_",$mergelanes, @{ $lanes{$_[0]} }, "_sorted_merged.bam OUTPUT=", '${outSampleDir}', "/$_[0]_", "lanes_",$mergelanes, @{ $lanes{$_[0]} }, "_sorted_merged_pmd.bam METRICS_FILE=", '${outSampleDir}', "/$_[0]_", "lanes_",$mergelanes, @{ $lanes{$_[0]} }, "_sorted_merged_pmdmetric &", "\n\n";
		}
	    }
	}
	print PMDUP "wait", "\n\n";
    }
###
#SamTools index on just created _sorted(_merged)_pmd.bam
###
    for (my $infile=0;$infile<scalar( @{ $Infiles_lane_noending{$_[0]} } );$infile++) { #For all files from MosaikAlign or BWA_sampe
	
	if ($infile eq $core_Counter*$maximumCores) { #Using only $maximumCores cores
	    
	    print PMDUP "wait", "\n\n";
	    $core_Counter=$core_Counter+1;
	}
	my $tempinfile = $Infiles_lane_noending{$_[0]}[$infile];
	print PMDUP "samtools index ", '${outSampleDir}', "/$tempinfile","_sorted_pmd.bam &","\n\n";  
    }
    print PMDUP "wait", "\n\n";
    if ( $script_parameters{'pPicardToolsMerge'} eq 1 && scalar( @{ $Infiles_lane_noending{$_[0]} }) > 1 ) { #SamTools index on merged file this round
	print PMDUP "samtools index ", '${outSampleDir}', "/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, "_sorted_merged_pmd.bam &","\n\n";
    }
    if (@picardToolsMergePreviousFiles ) { # SamTools index on files merged this round with merged file from previous round
	print PMDUP "samtools index ", '${outSampleDir}', "/$_[0]_", "lanes_",$mergelanes, @{ $lanes{$_[0]} }, "_sorted_merged_pmd.bam &","\n\n"; #$mergelanes should not have changed since if (@picardToolsMergePreviousFiles ) { so no need for regexp
	print PMDUP "wait", "\n\n";    
    }
    print PMDUP "wait", "\n\n";
    close(PMDUP);
    ParallelSampleIDSubmitJob($_[0],$filename,"all");
    return;
}

sub PicardMerge { 
#Merges all bam files using PicradToolds MergeSamFiles within each sampleid and optinally files generated previously. The merged files have to be sorted before attemting to merge.
#$_[0]= $sampleid
#$_[1] = aligner 
   
    `mkdir -p $script_parameters{'outDataDir'}/$_[0]/$_[1]/info;`; #Creates the aligner folder and info data file directory
    `mkdir -p $script_parameters{'outScriptDir'}/$_[0]/$_[1];`; #Creates the aligner folder and info data file directory
    $filename = "$script_parameters{'outScriptDir'}/$_[0]/$_[1]/picard_merge_$_[0].";
    
    Checkfnexists($filename, $fnend);

#Info and Logg
    print STDOUT "Creating sbatch script PicardMerge and writing script file(s) to: ", $filename, "\n";print MASTERL "Creating sbatch script PicardMerge and writing script file(s) to: ", $filename, "\n";
    print STDOUT "Sbatch script PicardMerge data files will be written to: ", $script_parameters{'outDataDir'},"/$_[0]/$_[1]", "\n";print MASTERL "Sbatch script PicardMerge data files will be written to: ", $script_parameters{'outDataDir'},"/$_[0]/$_[1]", "\n";
    
    open (PM, ">$filename") or die "Can't write to $filename: $!\n";
    
    print PM "#! /bin/bash -l", "\n";
    print PM "#SBATCH -A ", $script_parameters{'projectID'}, "\n";
    print PM "#SBATCH -p node -n $maximumCores", "\n";
    print PM "#SBATCH -C thin", "\n";	
    print PM "#SBATCH -t 20:00:00", "\n";
    
    print PM "#SBATCH -J PMerge_$_[1]", $_[0], "\n";
    print PM "#SBATCH -e $script_parameters{'outDataDir'}/$_[0]/$_[1]/info/picard_merge_$_[0].", $fnt ,".stderr.txt", "\n";
    print PM "#SBATCH -o $script_parameters{'outDataDir'}/$_[0]/$_[1]/info/picard_merge_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print PM "#SBATCH --mail-type=END", "\n";
	print PM "#SBATCH --mail-type=FAIL", "\n";
	print PM "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print PM 'echo "Running on: $(hostname)"',"\n\n";
    print PM "#Samples", "\n";
	
   
    print PM 'inSampleDir="',"$script_parameters{'outDataDir'}/$_[0]/$_[1]", '"', "\n";
    print PM 'outSampleDir="', "$script_parameters{'outDataDir'}/$_[0]/$_[1]", '"', "\n\n";

    if (scalar( @{ $Infiles_lane_noending{$_[0]} } ) > 1) { #Check that we have something to merge and then merge before merging with previously merged files
	
	for (my $infile=0;$infile<scalar( @{ $Infiles_lane_noending{$_[0]} } );$infile++) { #For all files from 
	    my $tempinfile = $Infiles_lane_noending{$_[0]}[$infile];
	    
	    if ($infile eq 0) {

		print PM "java -Xmx4g -jar $picardPath/MergeSamFiles.jar TMP_DIR=/proj/$script_parameters{'projectID'}/private/nobackup",'/$SLURM_JOB_ID ', "OUTPUT=", '${outSampleDir}', "/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, "_sorted_merged.bam ";
	    }
	    
	    print PM "INPUT=", '${inSampleDir}', "/$tempinfile","_sorted.bam ";   
	}
	print PM "\n\nsamtools index ", '${outSampleDir}', "/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, "_sorted_merged.bam", "\n\n";
	print PM "wait", "\n\n";
    }
    if ( (@picardToolsMergePreviousFiles) && (scalar( @{ $Infiles_lane_noending{$_[0]} } ) > 1) ) { #merge previously merged files with merged files generated this run
	
	for (my $mergefile=0;$mergefile<scalar(@picardToolsMergePreviousFiles);$mergefile++) {
	    
	    if ($picardToolsMergePreviousFiles[$mergefile] =~ /$_[0]/) { #Look for sampleID in previously generated file to be merged with current run to be able to merge correct files
		if ($picardToolsMergePreviousFiles[$mergefile] =~ /lane(\d+)|s_(\d+)/) { #Look for lanes_ or lane\d in previously generated file to be merged with current run to be able to extract previous lanes
		    
		    my $mergelanes; if($1) {$mergelanes = $1;} else {$mergelanes = $2;} #Make sure to always supply lanes from previous regexp
		    
		    print PM "java -Xmx4g -jar $picardPath/MergeSamFiles.jar TMP_DIR=/proj/$script_parameters{'projectID'}/private/nobackup",'/$SLURM_JOB_ID ', "OUTPUT=", '${outSampleDir}', "/$_[0]_", "lanes_",$mergelanes, @{ $lanes{$_[0]} }, "_sorted_merged.bam INPUT=", '${outSampleDir}', "/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, "_sorted_merged.bam INPUT=", $picardToolsMergePreviousFiles[$mergefile], "\n\n"; #$mergelanes contains lane info on previous merge, $Infiles_lane_noending{$_[0]}[0] uses @RG for very first .bam file to include read group for subsequent merges. 
		    print PM "samtools index ", '${outSampleDir}', "/$_[0]_", "lanes_",$mergelanes, @{ $lanes{$_[0]} }, "_sorted_merged.bam ","\n\n";
		}
	    }
	}
    }
    elsif ( @picardToolsMergePreviousFiles ) { #merge previously merged files with single file generated this run
	
	for (my $mergefile=0;$mergefile<scalar(@picardToolsMergePreviousFiles);$mergefile++) {
	    
	    if ($picardToolsMergePreviousFiles[$mergefile] =~ /$_[0]/) { #Look for sampleID in previously generated file to be merged with current run to be able to merge correct files
		if ($picardToolsMergePreviousFiles[$mergefile] =~ /lane(\d+)|s_(\d+)/) { #Look for lanes_ or lane\d in previously generated file to be merged with current run to be able to extract previous lanes
		    
		    my $mergelanes; if($1) {$mergelanes = $1;} else {$mergelanes = $2;} #Make sure to always supply lanes from previous regexp
		    my $tempinfile = $Infiles_lane_noending{$_[0]}[0]; #Can only be 1 element in array due to previous if statement
		    
		    print PM q?java -Xmx4g -jar ?.$picardPath.q?/MergeSamFiles.jar TMP_DIR=/proj/?.$script_parameters{'projectID'}.q?/private/nobackup/$SLURM_JOB_ID OUTPUT=${outSampleDir}/?.$_[0].q?_lanes_?.$mergelanes.q??.qq?@{ $lanes{$_[0]} }?.q?_sorted_merged.bam INPUT=${inSampleDir}/?.$tempinfile.q?_sorted.bam INPUT=?.$picardToolsMergePreviousFiles[$mergefile],"\n\n"; #$mergelanes contains lane info on previous merge, $Infiles_lane_noending{$_[0]}[0] uses @RG for very first .bam file to include read group for subsequent merges. 
		    print PM q?samtools index ${outSampleDir}/?.$_[0].q?_lanes_?.$mergelanes.q??.qq?@{ $lanes{$_[0]} }?.q?_sorted_merged.bam?, "\n\n";
		}
	    }
	}
    }
    close(PM);
    ParallelSampleIDSubmitJob($_[0],$filename,"all");
    return;
}

sub SamtoolsSortIndex { 
#Sort and indexes bam files using samtools sort and samtools index
#$_[0]= $sampleid
#$_[1] = aligner

    my $t=0;
    my $k=0;
    my $infilesize;
    for (my $infile=0;$infile<scalar( @{ $Infiles_lane_noending{$_[0]} } );$infile++) { #For all files from

	`mkdir -p $script_parameters{'outDataDir'}/$_[0]/$_[1]/info;`; #Creates the aligner folder and info data file directory
	`mkdir -p $script_parameters{'outScriptDir'}/$_[0]/$_[1];`; #Creates the aligner folder and info data file directory
	$filename = "$script_parameters{'outScriptDir'}/$_[0]/$_[1]/samToolsSort_index_$Infiles_lane_noending{$_[0]}[$infile].";
	if ($infiles{$_[0]}[$infile] =~/.fastq.gz$/) { #Files are already gz and presently the scalar for compression has not been investigated. Therefore no automatic time allocation can be performed.
	    if ($wholeGenomeSequencing == 1) {
		$t = 25;  
	    }
	    else {
		$t = 15;
	    }
	}
	else { #Files are in fastq format
	    $infilesize = -s "$indirpath{$_[0]}/$infiles{$_[0]}[$infile+$k]"; # collect .fastq file size to enable estimation of time required for sort & index, +$k for syncing multiple infiles per sampleID. Hence, filesize will be calculated on read1 (should not matter).	   
	    
	    if ($script_parameters{'pMosaikBuild'} || $script_parameters{'pMosaikAlign'} || ($script_parameters{'aligner'} eq "mosaik")) {
		$t = ceil($infilesize/(1700000*60*60)); #1700000 is a constant calculated from the filesize and time needed for procesing in samtools-0.1.12-10 sort and index and 60*60 is to scale to hours.
	    }
	    if ($script_parameters{'pBwaAln'} || $script_parameters{'pBwaSampe'} || ($script_parameters{'aligner'} eq "bwa")) {
		$t = ceil($infilesize/(1700000*60*60)); #1700000 is a constant calculated from the filesize and time needed for procesing in samtools-0.1.12-10 sort and index and 60*60 is to scale to hours.	    
	    }
	}
	
	Checkfnexists($filename, $fnend);
	
#Info and Logg
	print STDOUT "Creating sbatch script Samtools sort & index and writing script file(s) to: ", $filename, "\n";print MASTERL "Creating sbatch script Samtools sort & index and writing script file(s) to: ", $filename, "\n";
	print STDOUT "Sbatch script Samtools sort & index data files will be written to: ", $script_parameters{'outDataDir'},"/$_[0]/$_[1]", "\n";print MASTERL "Sbatch script Samtools sort & index data files will be written to: ", $script_parameters{'outDataDir'},"/$_[0]/$_[1]", "\n";
	
	open (STSI, ">$filename") or die "Can't write to $filename: $!\n";
	
	print STSI "#! /bin/bash -l", "\n";
	print STSI "#SBATCH -A ", $script_parameters{'projectID'}, "\n";
	print STSI "#SBATCH -p node -n 1", "\n";
	print STSI "#SBATCH -C thin", "\n";	
	print STSI "#SBATCH -t $t:00:00", "\n";
	
	print STSI "#SBATCH -J STSI_$_[1]", $_[0], "\n";
	print STSI "#SBATCH -e $script_parameters{'outDataDir'}/$_[0]/$_[1]/info/samToolsSort_index_$Infiles_lane_noending{$_[0]}[$infile].", $fnt ,".stderr.txt", "\n";
	print STSI "#SBATCH -o $script_parameters{'outDataDir'}/$_[0]/$_[1]/info/samToolsSort_index_$Infiles_lane_noending{$_[0]}[$infile].", $fnt ,".stdout.txt", "\n";
    
	unless ($em eq 0) {
	    
	    print STSI "#SBATCH --mail-type=END", "\n";
	    print STSI "#SBATCH --mail-type=FAIL", "\n";
	    print STSI "#SBATCH --mail-user=$em", "\n\n";
	    
	}
	
	print STSI 'echo "Running on: $(hostname)"',"\n\n";
	print STSI "#Samples", "\n";
	print STSI 'inSampleDir="',"$script_parameters{'outDataDir'}/$_[0]/$_[1]", '"', "\n";
	print STSI 'outSampleDir="', "$script_parameters{'outDataDir'}/$_[0]/$_[1]", '"', "\n\n";
    
###	
#SamTools Sort
###	
	my $tempinfile = $Infiles_lane_noending{$_[0]}[$infile];
	print STSI "samtools sort ", '${inSampleDir}', "/$tempinfile",'.', "bam ", '${outSampleDir}', "/$tempinfile","_sorted", "\n\n"; #Samtools sort adds .bam ending
	
	print STSI "wait", "\n\n";
###	
#SamTools Index
###	
	print STSI "samtools index ", '${inSampleDir}', "/$tempinfile","_sorted.bam", "\n\n";
	$k++;
	close(STSI);
	ParallelSampleIDSubmitJob($_[0],$filename,$Infiles_lane_noending{$_[0]}[$infile]);
    } 
    return;
}

sub BWA_sampe {
#Generates sbatch scripts for BWA sampe on files generated by BWA aln
    
    `mkdir -p $script_parameters{'outDataDir'}/$_[0]/bwa/info;`; #Creates the bwa folder and info data file directory
    `mkdir -p $script_parameters{'outScriptDir'}/$_[0]/bwa;`; #Creates the bwa script directory

    my $k=0;
    my $t=0;
    my $infilesize;
    for (my $infile=0;$infile<( scalar( @{ $Infiles_lane_noending{$_[0]} }) );$infile++) { #For all files from BWA aln but process in the same command i.e. both reads per align call
	if ($infiles{$_[0]}[$infile] =~/.fastq.gz$/) { #Files are already gz and presently the scalar for compression has not been investigated. Therefore no automatic time allocation can be performed.
	    if ($wholeGenomeSequencing == 1) {
		$t = 60;  
	    }
	    else {
		$t = 30;
	    }
	}
	else { #Files are in fastq format	
	    $infilesize = -s "$indirpath{$_[0]}/$infiles{$_[0]}[$infile+$k]"; # collect .fastq file size to enable estimation of time required for aligning, +1 for syncing multiple infiles per sampleID. Hence, filesize will be calculated on read2 (should not matter).
	    $t = ceil(($infilesize/238)/(3000*60*60)); #238 is a scalar estimating the number of reads depending on filesize. 3500 is the number of reads/s in Bwa_sampe-0.6.1 plus samtools-0.1.12-10 view sam to bam conversion and 60*60 is to scale to hours. (4600 BWA-0.5.9)
	}
	$filename = "$script_parameters{'outScriptDir'}/$_[0]/bwa/bwa_sampe_$Infiles_lane_noending{$_[0]}[$infile].";
	Checkfnexists($filename, $fnend);

#Info and Logg
	print STDOUT "Creating sbatch script BWA_Sampe and writing script file(s) to: ", $filename, "\n";print MASTERL "Creating sbatch script BWA_Sampe and writing script file(s) to: ", $filename, "\n";
	print STDOUT "Sbatch script BWA_Sampe data files will be written to: ", $script_parameters{'outDataDir'},"/$_[0]/$script_parameters{'aligner'}", "\n";print MASTERL "Sbatch script BWA_Sampe data files will be written to: ", $script_parameters{'outDataDir'},"/$_[0]/$script_parameters{'aligner'}", "\n";
	
	open (BWAS, ">$filename") or die "Can't write to $filename: $!\n";
	
	print BWAS "#! /bin/bash -l", "\n";
	print BWAS "#SBATCH -A ", $script_parameters{'projectID'}, "\n";
	print BWAS "#SBATCH -p node -n $maximumCores ", "\n";
	print BWAS "#SBATCH -C thin", "\n";
	print BWAS "#SBATCH -t $t:00:00", "\n";
	print BWAS "#SBATCH -J BWA_sampe", "$_[0]_", "\n";
	print BWAS "#SBATCH -e $script_parameters{'outDataDir'}/$_[0]/bwa/info/bwa_sampe_$Infiles_lane_noending{$_[0]}[$infile].", $fnt, ".stderr.txt", "\n";
	print BWAS "#SBATCH -o $script_parameters{'outDataDir'}/$_[0]/bwa/info/bwa_sampe_$Infiles_lane_noending{$_[0]}[$infile].", $fnt, ".stdout.txt", "\n";
	
	unless ($em eq 0) {
	    
	    print BWAS "#SBATCH --mail-type=END", "\n";
	    print BWAS "#SBATCH --mail-type=FAIL", "\n";
	    print BWAS "#SBATCH --mail-user=$em", "\n\n";
	    
	}
	
	print BWAS 'echo "Running on: $(hostname)"',"\n\n";
	print BWAS "module load bioinfo-tools", "\n\n"; 
	print BWAS "#Samples", "\n";
	print BWAS 'inSampleDir="', "$script_parameters{'outDataDir'}/$_[0]/bwa", '"', "\n\n";
	print BWAS 'inSampleDir2="', $indirpath{$_[0]}, '"', "\n\n"; #Fastq path
	print BWAS 'outSampleDir="', "$script_parameters{'outDataDir'}/$_[0]/bwa", '"', "\n\n";
	print BWAS "#Reference archive", "\n";
	print BWAS 'referenceArchive="', "$script_parameters{'referencesDir'}", '"', "\n\n";
	
	my $tempinfile = $infiles{$_[0]}[$infile+$k]; #For required .fastq file
	my $tempinfile2 = $infiles{$_[0]}[ ($infile+$k+1)]; # #For required .fastq file (Paired read)   
	
	print BWAS "bwa sampe -r ", '"@RG\tID:',$Infiles_bothstrands_noending{$_[0]}[$infile+$k],'\tSM:',"$_[0]",'\tPL:ILLUMINA" ','${referenceArchive}/', $humanGenomeReference,' ${inSampleDir}', "/$Infiles_bothstrands_noending{$_[0]}[$infile+$k].sai ", '${inSampleDir}', "/$Infiles_bothstrands_noending{$_[0]}[ ($infile+$k+1) ].sai ", '${inSampleDir2}', "/$tempinfile ", '${inSampleDir2}', "/$tempinfile2 > ", '${outSampleDir}', "/$Infiles_lane_noending{$_[0]}[$infile].sam", "\n\n";
	
	print BWAS "samtools view -bS ",'${inSampleDir}', "/$Infiles_lane_noending{$_[0]}[$infile].sam > ", '${outSampleDir}', "/$Infiles_lane_noending{$_[0]}[$infile].bam", "\n\n"; 
	print BWAS "rm ", '${inSampleDir}', "/$Infiles_lane_noending{$_[0]}[$infile].sam";
	$k++;	
	close(BWAS);
	ParallelSampleIDSubmitJob($_[0],$filename,$Infiles_lane_noending{$_[0]}[$infile]);
    }
    return;
}

sub BWA_aln {
    
#Generates sbatch scripts for BWA aln on files generated by platform
    
    `mkdir -p $script_parameters{'outDataDir'}/$_[0]/bwa/info;`; #Creates the bwa folder and info data file directory
    `mkdir -p $script_parameters{'outScriptDir'}/$_[0]/bwa;`; #Creates the bwa script directory
    
    $filename = "$script_parameters{'outScriptDir'}/$_[0]/bwa/bwa_aln_$_[0].";
    Checkfnexists($filename, $fnend);

#Info and Logg
    print STDOUT "Creating sbatch script BWA_Aln and writing script file(s) to: ", $filename, "\n";print MASTERL "Creating sbatch script BWA_Aln and writing script file(s) to: ", $filename, "\n";
    print STDOUT "Sbatch script BWA_Aln data files will be written to: ", $script_parameters{'outDataDir'},"/$_[0]/$script_parameters{'aligner'}", "\n";print MASTERL "Sbatch script BWA_Aln data files will be written to: ", $script_parameters{'outDataDir'},"/$_[0]/$script_parameters{'aligner'}", "\n";
    
    open (BWAA, ">$filename") or die "Can't write to $filename: $!\n";
    
    my $t = ceil(2.5*scalar( @{ $Infiles_lane_noending{$_[0]} })); #One full lane on Hiseq takes approx. 2,5 h for BWA_aln to process, round up to nearest full hour.
    
    print BWAA "#! /bin/bash -l", "\n";
    print BWAA "#SBATCH -A ", $script_parameters{'projectID'}, "\n";
    print BWAA "#SBATCH -p node -n $maximumCores ", "\n";
    print BWAA "#SBATCH -C thin", "\n";
    print BWAA "#SBATCH -t $t:00:00", "\n";
    print BWAA "#SBATCH -J BWA_aln", "$_[0]_", "\n";
    print BWAA "#SBATCH -e $script_parameters{'outDataDir'}/$_[0]/bwa/info/bwa_aln_$_[0].", $fnt, ".stderr.txt", "\n";
    print BWAA "#SBATCH -o $script_parameters{'outDataDir'}/$_[0]/bwa/info/bwa_aln_$_[0].", $fnt, ".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print BWAA "#SBATCH --mail-type=END", "\n";
	print BWAA "#SBATCH --mail-type=FAIL", "\n";
	print BWAA "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print BWAA 'echo "Running on: $(hostname)"',"\n\n";
    #print BWAA "module load bioinfo-tools", "\n\n"; 
    #print BWAA "module load bwa/0.5.9", "\n\n";
    print BWAA "#Samples", "\n";
    print BWAA 'inSampleDir="', $indirpath{$_[0]}, '"', "\n\n";
    print BWAA 'outSampleDir="', "$script_parameters{'outDataDir'}/$_[0]/bwa", '"', "\n\n";
    print BWAA "#Reference archive", "\n";
    print BWAA 'referenceArchive="', "$script_parameters{'referencesDir'}", '"', "\n\n";
    
    my $core_Counter=1;    
    for (my $infile=0;$infile<scalar( @{ $infiles{$_[0]} });$infile++) {
	
	if ($infile eq $core_Counter*$maximumCores) { #Using only $maximumCores cores
	    
	    print BWAA "wait", "\n\n";
	    $core_Counter=$core_Counter+1;
	}
	my $tempinfile = $infiles{$_[0]}[$infile];

	print BWAA "bwa aln -k 1 -t 4 -n 3 -I -q $bwaAlnQualityTrimming ", '${referenceArchive}/', $humanGenomeReference,' ${inSampleDir}', "/$tempinfile > ",'${outSampleDir}', "/$Infiles_bothstrands_noending{$_[0]}[$infile].sai &", "\n\n"; 
    }
    print BWAA "wait", "\n\n";
    close(BWAA);
    SampleIDSubmitJob($_[0],$filename, 1);   
    return;
}

sub MosaikAlign {
#Generates sbatch scripts for MosaikAlign on files generated from MosaikBuild
#$_[0] = SampleID
    
    `mkdir -p $script_parameters{'outDataDir'}/$_[0]/mosaik/info;`; #Creates the mosaik folder and info data file directory
    `mkdir -p $script_parameters{'outScriptDir'}/$_[0]/mosaik;`; #Creates the mosaik script directory
    
    my $k=0;
    my $t=0;
    my $infilesize;
    for (my $infile=0;$infile<scalar( @{ $Infiles_lane_noending{$_[0]} } );$infile++) { #For all files from MosaikBuild (platform reads)
	if ($infiles{$_[0]}[$infile] =~/.fastq.gz$/) { #Files are already gz and presently the scalar for compression has not been investigated. Therefore no automatic time allocation can be performed.
	    if ($wholeGenomeSequencing == 1) {
		$t = 80;  
	    }
	    else {
		$t = 40;
	    }
	}
	else { #Files are in fastq format
	    if (-e "$indirpath{$_[0]}/$infiles{$_[0]}[$infile+$k]") {
		$infilesize = -s "$indirpath{$_[0]}/$infiles{$_[0]}[$infile+$k]"; # collect .fastq file size to enable estimation of time required for aligning, +$k for syncing multiple infiles per sampleID. Hence, filesize will be calculated on read1 (should not matter).      
		$t = ceil(($infilesize/238)/(650*60*60)); #238 is a scalar estimating the number of reads depending on filesize. 650 is the number of reads/s in MosaikAlign-2.1.52 and 60*60 is to scale to hours.
	    }	    
	} 
	$filename = "$script_parameters{'outScriptDir'}/$_[0]/mosaik/mosaikAlign_$Infiles_lane_noending{$_[0]}[$infile].";
	Checkfnexists($filename, $fnend);
	
#Info and Logg
	print STDOUT "Creating sbatch script MosaikAlign and writing script file(s) to: ", $filename, "\n";print MASTERL "Creating sbatch script MosaikAlign and writing script file(s) to: ", $filename, "\n";
	print STDOUT "Sbatch script MosaikAlign data files will be written to: ", $script_parameters{'outDataDir'},"/$_[0]/$script_parameters{'aligner'}", "\n";print MASTERL "Sbatch script MosaikAlign data files will be written to: ", $script_parameters{'outDataDir'},"/$_[0]/$script_parameters{'aligner'}", "\n";
	
	open (MosA, ">$filename") or die "Can't write to $filename: $!\n";
	
	print MosA "#! /bin/bash -l", "\n";
	print MosA "#SBATCH -A ", $script_parameters{'projectID'}, "\n";
	print MosA "#SBATCH -p node -n $maximumCores ", "\n";
	print MosA "#SBATCH -C thin", "\n";
	print MosA "#SBATCH -t $t:00:00", "\n";
	print MosA "#SBATCH -J MoA", "$_[0]_",$k, "\n";
	print MosA "#SBATCH -e $script_parameters{'outDataDir'}/$_[0]/mosaik/info/mosaikAlign_$Infiles_lane_noending{$_[0]}[$infile].$fnt.stderr.txt", "\n";
	print MosA "#SBATCH -o $script_parameters{'outDataDir'}/$_[0]/mosaik/info/mosaikAlign_$Infiles_lane_noending{$_[0]}[$infile].$fnt.stdout.txt", "\n";
	
	unless ($em eq 0) {
	    
	    print MosA "#SBATCH --mail-type=END", "\n";
	    print MosA "#SBATCH --mail-type=FAIL", "\n";
	    print MosA "#SBATCH --mail-user=$em", "\n\n";
	    
	}
	
	print MosA 'echo "Running on: $(hostname)"',"\n\n";
	print MosA "mkdir -p /scratch/mosaik_tmp", "\n";
	print MosA "export MOSAIK_TMP=/scratch/mosaik_tmp", "\n\n";
	print MosA "#Samples", "\n";
	print MosA 'inSampleDir="',"$script_parameters{'outDataDir'}/$_[0]/mosaik", '"', "\n";
	print MosA 'outSampleDir="', "$script_parameters{'outDataDir'}/$_[0]/mosaik", '"', "\n\n";
	print MosA "#Reference archive", "\n";
	print MosA 'referenceArchive="', "$script_parameters{'referencesDir'}", '"', "\n\n";
	
	my $tempinfile = $Infiles_lane_noending{$_[0]}[$infile];
	print MosA "MosaikAligner -in ", '${inSampleDir}', "/$tempinfile",'.', "dat -out ",'${outSampleDir}', "/$Infiles_lane_noending{$_[0]}[$infile] -ia ", '${referenceArchive}', "/$mosaikAlignReference -annpe ",'${referenceArchive}', "/$mosaikAlignNeuralNetworkPeFile -annse ",'${referenceArchive}', "/$mosaikAlignNeuralNetworkSeFile -hs 15 -mm 4 -mhp 100 -ls 100 -act 35 -bw 35 -j ", '${referenceArchive}', "/$mosaikJumpDbStub -p $maximumCores", "\n\n";
	
	$k++; #Tracks nr of sbatch scripts
	close(MosA);
	ParallelSampleIDSubmitJob($_[0],$filename, $Infiles_lane_noending{$_[0]}[$infile]);
    }
    return;
}

sub MosaikBuild {
    
#Generates a sbatch script and runs MosaikBuild on reads from Illumina platform or fastq_filtered reads
#Starts mosaikAlign from sbatch when mosaikBuild is finished for all files and submits mosaikSort with dependency afterrok of mosaikAlign  
#$_[0] = sampleid
    
    `mkdir -p $script_parameters{'outDataDir'}/$_[0]/mosaik/info;`; #Creates the mosaik folder and info data file directory
    `mkdir -p $script_parameters{'outScriptDir'}/$_[0]/mosaik;`; #Creates the mosaik script directory
    
    $filename = "$script_parameters{'outScriptDir'}/$_[0]/mosaik/mosaikBuild_$_[0].";
    Checkfnexists($filename, $fnend);
    
#Info and Logg
    print STDOUT "Creating sbatch script MosaikBuild and writing script file(s) to: ", $filename, "\n";print MASTERL "Creating sbatch script MosaikBuild and writing script file(s) to: ", $filename, "\n";
    print STDOUT "Sbatch script MosaikBuild data files will be written to: ", $script_parameters{'outDataDir'},"/$_[0]/$script_parameters{'aligner'}", "\n";print MASTERL "Sbatch script MosaikBuild data files will be written to: ", $script_parameters{'outDataDir'},"/$_[0]/$script_parameters{'aligner'}", "\n";

    my $t = ceil(2.5*scalar( @{ $Infiles_lane_noending{$_[0]} })); #One full lane on Hiseq takes approx. 1 h for MosaikBuild to process (compressed format, uncompressed 0.5 h), round up to nearest full hour.
    
    open (MosB, ">$filename") or die "Can't write to $filename: $!\n";
    
    print MosB "#! /bin/bash -l", "\n";
    print MosB "#SBATCH -A ", $script_parameters{'projectID'}, "\n";
    print MosB "#SBATCH -p node -n $maximumCores ", "\n";
    print MosB "#SBATCH -C thin", "\n";
    print MosB "#SBATCH -t $t:00:00", "\n";
    print MosB "#SBATCH -J MoB", $_[0], "\n";
    print MosB "#SBATCH -e $script_parameters{'outDataDir'}/$_[0]/mosaik/info/mosaikBuild_$_[0].", $fnt ,".stderr.txt", "\n";
    print MosB "#SBATCH -o $script_parameters{'outDataDir'}/$_[0]/mosaik/info/mosaikBuild_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print MosB "#SBATCH --mail-type=END", "\n";
	print MosB "#SBATCH --mail-type=FAIL", "\n";
	print MosB "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print MosB 'echo "Running on: $(hostname)"',"\n\n";
    print MosB 'outSampleDir="', "$script_parameters{'outDataDir'}/$_[0]/mosaik", '"', "\n\n";
    
    my $core_Counter=1;
    my $allp=0; #Required to portion out $maximumCores files before wait and to track the MosB outfiles to correct lane
    
    print MosB 'inSampleDir="', $indirpath{$_[0]}, '"', "\n\n";
    
    for (my $infile=0;$infile<( scalar( @{ $infiles{$_[0]} }) -1);$infile++) {
	
	if ($allp eq $core_Counter*$maximumCores) { #Using only $maximumCores cores
	    
	    print MosB "wait", "\n\n";
	    $core_Counter=$core_Counter+1;
	}
	my $tempinfile = $infiles{$_[0]}[$infile];
	my $tempinfile2 = $infiles{$_[0]}[ ($infile+1)]; #Paired read
	
	$infile = $infile+1; #To correct for reading 2 files at once
	print MosB "MosaikBuild -id $Infiles_bothstrands_noending{$_[0]}[$infile] -sam $_[0] -st illumina_long -mfl $mosaikBuildMeanFragLength -q ", '${inSampleDir}', "/$tempinfile -q2 ", '${inSampleDir}', "/$tempinfile2 -out ", 
	'${outSampleDir}', "/$Infiles_lane_noending{$_[0]}[ ($allp)]",'.',"dat &", "\n\n";
	$allp++; #Track nr of printed so that wait can be printed between hashes
    }
    print MosB "wait", "\n\n";    
    close(MosB);
    SampleIDSubmitJob($_[0],$filename, 1); 
    return;
}   

sub FastQC {
    
#Generates a sbatch script and runs FASTQC
#$_[0] = sampleid

    `mkdir -p $script_parameters{'outDataDir'}/$_[0]/fastqc/info;`; #Creates the fastqc folder and info data file directory
    `mkdir -p $script_parameters{'outScriptDir'}/$_[0]/fastqc;`; #Creates the fastqc script directory
    
    $filename = "$script_parameters{'outScriptDir'}/$_[0]/fastqc/fastqc_$_[0].";
    Checkfnexists($filename, $fnend);

#Info and Logg
    print STDOUT "Creating sbatch script FastQC and writing script file(s) to: ", $filename, "\n";print MASTERL "Creating sbatch script Sample check FastQC and writing script file(s) to: ", $filename, "\n";
    print STDOUT "Sbatch script FastQC data files will be written to: ", $script_parameters{'outDataDir'},"/$_[0]/fastqc", "\n";print MASTERL "Sbatch script Sample check FastQC data files will be written to: ", $script_parameters{'outDataDir'},"/$_[0]/fastqc", "\n";


    my $t = ceil(0.5*scalar( @{ $infiles{$_[0]} })); #One full lane on Hiseq takes approx. 0.5 h for FASTQC to process, round up to nearest full hour.
    
    open (FASTQC, ">$filename") or die "Can't write to $filename: $!\n";
    
    print FASTQC "#! /bin/bash -l", "\n";
    print FASTQC "#SBATCH -A ", $script_parameters{'projectID'}, "\n";
    print FASTQC "#SBATCH -p node -n $maximumCores ", "\n";
    print FASTQC "#SBATCH -C thin", "\n";
    print FASTQC "#SBATCH -t $t:00:00", "\n";
    print FASTQC "#SBATCH -J FQC", $_[0], "\n";
    print FASTQC "#SBATCH -e $script_parameters{'outDataDir'}/$_[0]/fastqc/info/fastqc_$_[0].", $fnt ,".stderr.txt", "\n";
    print FASTQC "#SBATCH -o $script_parameters{'outDataDir'}/$_[0]/fastqc/info/fastqc_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print FASTQC "#SBATCH --mail-type=END", "\n";
	print FASTQC "#SBATCH --mail-type=FAIL", "\n";
	print FASTQC "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print FASTQC 'echo "Running on: $(hostname)"',"\n\n";
    print FASTQC "cd $indirpath{$_[0]}", "\n\n";
    print FASTQC "#Samples", "\n";
    print FASTQC 'outSampleDir="', "$script_parameters{'outDataDir'}/$_[0]/fastqc", '"', "\n\n";
    
    my $core_Counter=1;
    
    for (my $infile=0;$infile<scalar( @{ $infiles{$_[0]} });$infile++) {
	
	if ($infile eq $core_Counter*$maximumCores) { #Using only $maximumCores cores
	    
	    print FASTQC "wait", "\n\n";
	    $core_Counter=$core_Counter+1;
	}
	my $tempinfile = $infiles{$_[0]}[$infile];
	print FASTQC "fastqc ", $tempinfile, ' -o ${outSampleDir}';
	print FASTQC " &", "\n\n";

    }
    
    print FASTQC "wait", "\n";    
    
    close(FASTQC);
    SampleIDSubmitJob($_[0],$filename, 0);
    return;
}

sub GZipfastq { 
#Generates sbatch scripts for gziping, which will not be started but can be used whenever needed.

    `mkdir -p $script_parameters{'outDataDir'}/$_[0]/gzip/info;`; #Creates the gzip folder and info data file directory
    `mkdir -p $script_parameters{'outScriptDir'}/$_[0]/gzip;`; #Creates the gzip script folder    
    $filename = "$script_parameters{'outScriptDir'}/$_[0]/gzip/gzipFastq_$_[0].";
    Checkfnexists($filename, $fnend);

#Info and Logg
    print STDOUT "Creating sbatch script GzipFastq and writing script file(s) to: ", $filename, "\n";print MASTERL "Creating sbatch script GzipFastq and writing script file(s) to: ", $filename, "\n";
    print STDOUT "Sbatch script GzipFastq data files will be written to: ", $script_parameters{'outDataDir'},"/$_[0]/fastq", "\n";print MASTERL "Sbatch script GzipFastq data files will be written to: ", $script_parameters{'outDataDir'},"/$_[0]/fastq", "\n";
    
    my $t = ceil(1.5*scalar( @{ $infiles{$_[0]} })); #One full lane on Hiseq takes approx. 1.5 h for gzip to process, round up to nearest full hour.
    open (GZFASTQ, ">$filename") or die "Can't write to $filename: $!\n";
    
    print GZFASTQ "#! /bin/bash -l", "\n";
    print GZFASTQ "#SBATCH -A ", $script_parameters{'projectID'}, "\n";
    print GZFASTQ "#SBATCH -p node -n $maximumCores", "\n";
    print GZFASTQ "#SBATCH -C thin", "\n";	
    print GZFASTQ "#SBATCH -t $t:00:00", "\n";
    print GZFASTQ "#SBATCH -J GZFQ", $_[0], "\n";
    print GZFASTQ "#SBATCH -e $script_parameters{'outDataDir'}/$_[0]/gzip/info/gzipFastq_$_[0].", $fnt ,".stderr.txt", "\n";
    print GZFASTQ "#SBATCH -o $script_parameters{'outDataDir'}/$_[0]/gzip/info/gzipFastq_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print GZFASTQ "#SBATCH --mail-type=END", "\n";
	print GZFASTQ "#SBATCH --mail-type=FAIL", "\n";
	print GZFASTQ "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print GZFASTQ 'echo "Running on: $(hostname)"',"\n\n";
    print GZFASTQ "#Samples", "\n";
    print GZFASTQ "cd $indirpath{$_[0]}", "\n\n";
    print GZFASTQ 'inSampleDir="',"$indirpath{$_[0]}", '"', "\n\n";

    my $core_Counter=1;
    my $uncompressed_file_Counter = 0; #Used to print wait at the right times since infiles cannot be used (can be a mixture of .gz and .fast files)
    for (my $infile=0;$infile<scalar( @{ $infiles{$_[0]} });$infile++) {

	if ($infiles{$_[0]}[$infile] =~/.fastq$/) { #For files ending with .fastq required since there can be a mixture (also .fastq.gz) within the sample dir
	    if ($uncompressed_file_Counter eq $core_Counter*$maximumCores) { #Using only $maximumCores cores
		
		print GZFASTQ "wait", "\n\n";
		$core_Counter=$core_Counter+1;
	    }
	    my $tempinfile = $infiles{$_[0]}[$infile];
	    print GZFASTQ "gzip ", '${inSampleDir}', "/$tempinfile"," &", "\n\n";
	    $uncompressed_file_Counter++;
	    $infiles{$_[0]}[$infile] =~ s/.fastq/.fastq.gz/g; #Replace the .fastq ending with .fastq.gz since this will execute before fastQC screen and mosaikBuild, hence changing the original file name ending from ".fastq" to ".fastq.gz". 

	}
    }
    print GZFASTQ "wait", "\n\n";
    SampleIDSubmitJob($_[0],$filename, 1);
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
		    $sample_info{$local_familyID}{$line_info[0]}{'Sex'} = "M"; #Sex, M=Male
		}
		else { #Female
		   $sample_info{$local_familyID}{$line_info[0]}{'Sex'} = "F"; #Sex, F=Female
		}
		if ($4 eq "A") { #Affected
		    $sample_info{$local_familyID}{$line_info[0]}{'Disease_status'} = 1; #1=Affected
		}
		else { #Unaffected
		    $sample_info{$local_familyID}{$line_info[0]}{'Disease_status'} = 0; #0=Unaffected
		}
		
		if ($line_info[14]) { #Capture kit
		    my @capture_kits = split(";", $line_info[14]);
		    my $capture_kit =  pop(@capture_kits); #Use only the last capture kit since it should be the most interesting
		    
		    for my $supported_capture_kit (keys %supported_capture_kits) {
			if ($supported_capture_kit eq $capture_kit) {
			    if ($exomeTargetBed eq 0) { #No user supplied info on capture kit target BED-file. Add from pedigree file
				$sample_info{$local_familyID}{$line_info[0]}{'exomeTargetBed'} = $supported_capture_kits{$supported_capture_kit}; #capture kit Bed-file
			    }
			    if ($exomeTargetBedInfileList eq 0) { #No user supplied info on capture kit target BED-file infile list. Add from pedigree file
				$sample_info{$local_familyID}{$line_info[0]}{'exomeTargetBedInfileList'} = $supported_capture_kits{$supported_capture_kit}.".infile_list"; #capture kit target infile_list
			    }
			    if ($exomeTargetPaddedBedInfileList eq 0) { #No user supplied info on capture kit target BED-file infile list. Add from pedigree file
				$sample_info{$local_familyID}{$line_info[0]}{'exomeTargetPaddedBedInfileList'} = $supported_capture_kits{$supported_capture_kit}.".pad100.infile_list"; #capture kit padded target infile_list
			    }
			}
		    }
		}	
	    }
	    #push(@{$sample_info{$local_familyID}{$line_info[0]}},@line_info[2..4]); #Populate hash of array for Mother/Father/Child. Full hash: hash{FDN}{IDN}[Sex,Affected,Mother/Father/Child]
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

my $uncompressed_file_counter = 0;     

    for my $samplid ( keys %infiles ) { #For every sample id
	
	my $k=1;
	my $itrack=0; #Needed to be able to track when lanes are finished
	for (my $i=0;$i<scalar( @ { $infiles{ $samplid } } );$i++) { #Collects inputfiles for every fastq dir and remakes format
	    if ( $infiles{$samplid}[$i] =~ /\/?([^\.\/]+\.[^\.]+)\.lane(\d+)_([12FfRr])\.fastq.gz/ ) { #Parse fastq.gz 'old' format
		
		push( @ {$lanes{$samplid} }, $2);
		$Infiles_lane_noending{ $samplid }[$itrack]= "$1.$2"; #Save new format in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry into hash for both strands and .ending is removed (.fastq).
		$i++; #Skip second direction
		$itrack++; #Track for every lane finished
	    }
	    elsif ( $infiles{$samplid}[$i] =~ /\/?([^\.\/]+\.[^\.]+)\.lane(\d+)_([12FfRr])\.fastq/ ) { #Parse 'old' format

		push( @ {$lanes{$samplid} }, $2);
		$uncompressed_file_counter = 1;
		$Infiles_lane_noending{ $samplid }[$itrack]= "$1.$2"; #Save new format in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry into hash for both strands and .ending is removed (.fastq).
		$i++; #Skip second direction
		$itrack++; #Track for every lane finished
	    }
	    elsif ( $infiles{$samplid}[$i] =~ /(\d+)_(\d+)_([^_]+)_([^_]+)_(index[^_]+)_(\d).fastq.gz/ ) { #Parse fastq.gz 'new' format $1=lane, $2=date, $3=Flow-cell, $4=SampleID, $5=index,$6=direction
	
		push( @ {$lanes{$samplid} }, $1);
		$Infiles_lane_noending{ $samplid }[$itrack]= "$4.$2_$3_$5."."lane"."$1_$6"; #Save new format (sampleID_date_flow-cell_index_lane_direction) in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry into hash for both strands and .ending is removed (.fastq).
		$i++; #Skip second direction
		$itrack++; #Track for every lane finished
	    }
	    elsif ( $infiles{$samplid}[$i] =~/(\d+)_(\d+)_([^_]+)_([^_]+)_(index[^_]+)_(\d).fastq/ ) { #Parse 'new' format $1=lane, $2=date, $3=Flow-cell, $4=SampleID, $5=index,$6=direction
		
		push( @ {$lanes{$samplid} }, $1);
		$uncompressed_file_counter = 1;
		$Infiles_lane_noending{ $samplid }[$itrack]= "$4.$2_$3_$5."."lane"."$1_$6"; #Save new format (sampleID_date_flow-cell_index_lane_direction) in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry into hash for both strands and .ending is removed (.fastq).
		$i++; #Skip second direction
		$itrack++; #Track for every lane finished
	    }
	}
	$k=1;
	for (my $i=0;$i<scalar( @ { $infiles{ $samplid } } );$i++) { #Collects inputfiles for every fastq dir and remakes format
	    if ( $infiles{$samplid}[$i] =~ /\/?([^\.\/]+\.[^\.]+)\.lane(\d+_[12FfRr])\.fastq/ ) { #Parse 'old' format
		
		$Infiles_bothstrands_noending{ $samplid }[$i]= "$1.$2"; #Save new format in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry per strand and .ending is removed (.fastq).
		$itrack++; #Track for every lane finished
	    }
	    elsif ( $infiles{$samplid}[$i] =~ /(\d+)_(\d+)_([^_]+)_([^_]+)_(index[^_]+)_(\d).fastq/ ) { #Parse 'new' format
	
		$Infiles_bothstrands_noending{ $samplid }[$i]= "$4.$2_$3_$5."."lane"."$1_$6"; #Save new format in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry per strand and .ending is removed (.fastq).
		$itrack++; #Track for every lane finished
	    }
			    
	}
    }
return $uncompressed_file_counter;
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
    
    open (MASTERL, ">>$master_logg_name") or die "Can't write to $master_logg_name: $!\n"; #Open file run logg
    
    print MASTERL 
	"-p ", $script_parameters{'projectID'},
	" -familyID ", $script_parameters{'familyID'};
    if ( $script_parameters{'pedigreeFile'} ) {
	print MASTERL
	    " -pedigreeFile ", $script_parameters{'pedigreeFile'};
    }
    print MASTERL
	" -sampleIDs ", $script_parameters{'sampleIDs'},
	" -inFilesDir ", $script_parameters{'inFilesDir'},
	" -outDataDir ", $script_parameters{'outDataDir'},   
	" -inScriptDir ", $script_parameters{'inScriptDir'},
	" -outScriptDir ", $script_parameters{'outScriptDir'},
	" -referencesDir ", $script_parameters{'referencesDir'};
    if ( $script_parameters{'humanGenomeReference'} ) {
	print MASTERL
	    " -humanGenomeReference ", $script_parameters{'humanGenomeReference'};
    }
    print MASTERL
	" -pFastQC ", $script_parameters{'pFastQC'};
    if ($script_parameters{'aligner'}) {
	print MASTERL
	    " -aligner ", $script_parameters{'aligner'};
    }
    print MASTERL
	" -pMosaikBuild ", $script_parameters{'pMosaikBuild'},
	" -pMosaikAlign ", $script_parameters{'pMosaikAlign'},
	" -pBwaAln ", $script_parameters{'pBwaAln'},
	" -pBwaSampe ", $script_parameters{'pBwaSampe'},
	" -pSamToolsSort ", $script_parameters{'pSamToolsSort'},
	" -pPicardToolsMerge ", $script_parameters{'pPicardToolsMerge'};
    if ( scalar(@picardToolsMergePreviousFiles) ) {
	print MASTERL
	    " -picardToolsMergePreviousFiles ", $script_parameters{'picardToolsMergePreviousFiles'};
    }
    print MASTERL
	" -pPicardToolsMarkduplicates ", $script_parameters{'pPicardToolsMarkduplicates'};
    if ( $picardPath ) {
	print MASTERL
	    " -picardPath ", $script_parameters{'picardPath'};
    }
    print MASTERL
	" -pCalculateCoverage ", $script_parameters{'pCalculateCoverage'},
	" -pGenomeCoverageBED ", $script_parameters{'pGenomeCoverageBED'},
	" -pCoverageBED ", $script_parameters{'pCoverageBED'},
	" -pQaCompute ", $script_parameters{'pQaCompute'},
	" -pPicardToolsCollectMultipleMetrics ", $script_parameters{'pPicardToolsCollectMultipleMetrics'},
	" -pPicardToolsCalculateHSMetrics ", $script_parameters{'pPicardToolsCalculateHSMetrics'};
    if ( $identicalCaptureBedCounter eq scalar(@sampleIDs) ) { #Same capture kit for all sampleIDs
	print MASTERL 
	    " -exomeTargetBed ", $script_parameters{$sampleIDs[0]}{'exomeTargetBed'};
    }
    elsif ($exomeTargetBed) {
	print MASTERL 
	    " -exomeTargetBed ", $script_parameters{$sampleIDs[0]}{'exomeTargetBed'};
    }
    if ( $identicalCaptureBedIntervalCounter eq scalar(@sampleIDs) ) { #Same capture kit for all sampleIDs
	print MASTERL 
	    " -exomeTargetBedInfileList ", $script_parameters{$sampleIDs[0]}{'exomeTargetBedInfileList'},
	    " -exomeTargetPaddedBedInfileList ", $script_parameters{$sampleIDs[0]}{'exomeTargetPaddedBedInfileList'};
    }
    elsif ($exomeTargetBedInfileList) {
	print MASTERL
	    " -exomeTargetBedInfileList ", $script_parameters{$sampleIDs[0]}{'exomeTargetBedInfileList'},
	    " -exomeTargetPaddedBedInfileList ", $script_parameters{$sampleIDs[0]}{'exomeTargetPaddedBedInfileList'};
    }
    print MASTERL
	" -pRCovPlots ", $script_parameters{'pRCovPlots'},
	" -pGZip ", $script_parameters{'pGZip'},
	" -pRemovalRedundantFiles ", $script_parameters{'pRemovalRedundantFiles'}, "\n";
    #Note FileHandle MASTERL not closed
    return;
}
