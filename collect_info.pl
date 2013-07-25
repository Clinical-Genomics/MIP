#!/usr/bin/perl - w

use strict;
use warnings;

=for comment
Collects sample summary info and writes to a tab-sep file.
Works on data structure created by wgs_align and wgs_var_call pipelines.
Will look for merged files (when suitable) and if none is found within sampleID,
the highest ordered STDERR/STDOUT file or is used for collecting info.
=cut

# Copyright 2012 Henrik Stranneheim

=head1 SYNOPSIS

collect_info.pl -a [project ID] -p [projdir] -f [family ID...n] -o [out file]

=head2 COMMANDS AND OPTIONS

-ifd/--inFilesDir Infile directory

-a/--projectid The project ID (Mandatory)

-p/--projectdir Project dir (For instance: exomes)

-pedigree/--pedigreeFile (Supply whole path, defaults to "")

-f/--familyid Group id of samples to be compared, comma separated or just all for all familyIDs (defaults to "")

-o/--outfile The data file output (Supply whole path, defaults to summary_info.txt)

=head3 I/O

Input format ( familyID/sampleID )

Output format

1. tab-sep

=head4 Dependencies

Data structure created by wgs pipelines.

=cut
    
use Pod::Usage;
use Pod::Text;
use Getopt::Long;
use POSIX;

use vars qw($USAGE);

BEGIN {
    $USAGE =
	qq{collect_info.pl -a [project ID] -p [projdir] -s [sample ID...n] -f [family ID...n] -o [out file]
               -ifd/--inFilesDir Infile directory
	       -a/--projectid The project ID  (Mandatory)
               -p/--projectdir Project dir (For instance: exomes)  
               -pedigree/--pedigreeFile (Supply whole path, defaults to "")
               -f/--familyid Group id of samples to be compared, comma separated or just all for all familyIDs (defaults to "")
	       -o/--outfile The data file output (Supply whole path, defaults to summary_info.txt)};
    
}

my ($inFilesDir, $aid,$pid, $o, $projpath, $pedigreeFile, $help) = ("", 0,0,"summary_info.txt","");
my (@inid,@fid);
my (%sampleData, %familyMembers, %infiles, %indirpath, %Infiles_lane_noending, %lanes, %Infiles_bothstrands_noending, %sampleFlowCells);
#%infiles=from platform (Illumina), %indirpath for the path to infiles, %Infiles_lane_noending for MosaikBuild (one entry for both strands), %lanes for sample lanes, Infiles_bothstrands_noending for bwa_aln (one entry per strand)

GetOptions('ifd|inFilesDir:s' => \$inFilesDir,
           'a|projectid:s'  => \$aid,
	   'p|projectdir:s'  => \$pid,
	   'pedigree|pedigreeFile:s' => \$pedigreeFile, 
	   'f|familyid:s'  => \@fid, #Comma separated list
	   'o|outfile:s'  => \$o, #The data file output (Supply whole path, defaults to summary_info.txt)
	   'h|help' => \$help,
    );

die $USAGE if( $help );

if ($aid eq 0) {
    
    print STDERR "\n";
    print STDERR "Must supply a project ID", "\n\n";
    die $USAGE;
}
if ( scalar($pid) eq 0) {
    
    print STDERR "\n";
    print STDERR "Must supply a project dir", "\n\n";
    die $USAGE;
}
if ( scalar(@fid) eq 0) {
    
    print STDERR "\n";
    print STDERR "Must supply a family ID as a comma separated list", "\n\n";
    die $USAGE;
}
if ( $fid[0] =~/all/i ) {
@fid = `cd /bubo/proj/$aid/private/$pid;ls --ignore=*-*;`; #cd to proj dir and collect all familyIDs (ls -d) only directories
chomp(@fid);
}
else {
    @fid = split(/,/,join(',',@fid)); #Enables comma separeted list of family IDs
}

#Set projectpath 
$projpath = "/bubo/proj/$aid/private/$pid";

###
#MAIN
###

for (my $familyid=0;$familyid<scalar(@fid);$familyid++) {
    #ReadPedigreeFile($projpath, $fid[$familyid]); #Collects sampleID, pedigree and disease status etc
    ReadPedigreeFile($pedigreeFile, $fid[$familyid]); #Collects sampleID, pedigree and disease status etc
}
#Collect input files for all family and sampleIDs
for my $familyid ( keys %familyMembers ) { #For every family id

    print STDOUT "\nFamily ID: ", $familyid,"\n";

    for my $sampleid ( keys %{ $familyMembers{$familyid} } ) { #For every sample id
	#my @infiles = `cd /bubo/proj/$aid/private/$pid/$sampleid/fastq;ls *.fastq*;`; #cd to sample fastq dir and collect *.fastq* files
	my @infiles = `cd $inFilesDir/$sampleid/fastq;ls *.fastq*;`; #cd to sample fastq dir and collect *.fastq* files
	
	if (scalar(@infiles) > 0) {
	    print STDOUT "\nReads from Platform", "\n";
	    print STDOUT "\nSample ID: ", $sampleid,"\n";
	    print STDERR "Inputfiles", "\n", @ { $infiles{$familyid}{$sampleid}  =[@infiles] }, "\n"; #hash with {familyid}{sampleid} as keys and inputfiles as array 
	    chomp(@infiles);    #Remove newline from every entry in array
	    $infiles{$familyid}{$sampleid}  =[@infiles]; #Reload files into hash (kept above newline just for print STDERR)
	}
	else {
	    print STDERR "Could not find any infiles for Family ID: $familyid in /bubo/proj/$aid/private/$pid/$sampleid/fastq\n";
	}
    }
}

InfilesReFormat(); #Required to format infiles correctly for subsequent write to summary_info file

ReadFastqc($projpath);
ReadMosaik($projpath);
MarkDupMetrics($projpath);
AlignSumMetrics($projpath);
CalHSMetrics($projpath);
GenderCheck($projpath);
ReadGATK($projpath);
GATKVareValAll($projpath);
GATKVareValExome($projpath);
PedigreeCheck($projpath);
InbreedCheck($projpath);
WriteSummary($o); #Writes a summary report for all familyID and samplID supplied

###
#Sub Routines
###

sub InfilesReFormat {
    
#Code needed to reformat files for summary_info output. Extracts name and lane from infiles
    
    for my $familyid ( keys %infiles ) { #For every family id
	
	for my $sampleid ( keys %{ $infiles{$familyid} } ) { #For every sample id
	    
	    my $k=1;
	    my $itrack=0; #Needed to be able to track when lanes are finished
	    for (my $i=0;$i<scalar( @ { $infiles{$familyid}{$sampleid} } );$i++) { #Collects inputfiles for every fastq dir and remakes format
		if ( $infiles{$familyid}{$sampleid}[$i] =~ /\/?[^\.\/]+\.([^\.]+)\.lane(\d+)_([12FfRr])\.\S+/ ) { #Parse 'old' format
		    chomp($2);
		    $sampleData{$familyid}{$sampleid}{'Flow-cell'} .= "$1;";
		    $sampleData{$familyid}{$sampleid}{'Lane'} .= "$2;";
		    $sampleData{$familyid}{$sampleid}{'Date'} .= "Na;";
		    $sampleData{$familyid}{$sampleid}{'Index'} .= "Na;";
		    push( @ {$lanes{$sampleid} }, $2);
		    push(@ {$sampleFlowCells{$sampleid}}, $1);
		}
		elsif ( $infiles{$familyid}{$sampleid}[$i] =~ /(\d+)_(\d+)_([^_]+)_([^_]+)_(index[^_]+)_(\d).fastq/ ) { #Parse 'new' format
		    chomp($1);
		    $sampleData{$familyid}{$sampleid}{'Lane'} .= "$1;";
		    $sampleData{$familyid}{$sampleid}{'Date'} .= "$2;";
		    $sampleData{$familyid}{$sampleid}{'Flow-cell'} .= "$3;";
		    $sampleData{$familyid}{$sampleid}{'Index'} .= "$5;";  
		    $sampleData{$familyid}{$sampleid}{'Unmerged_post_MosA'} = "$4.$2_$3_$5."."lane"."$1_$6"; #For unmerged files after mosaikAlign
		    push( @ {$lanes{$sampleid} }, $1);
		    push(@ {$sampleFlowCells{$sampleid}}, $3);
		}
		if ( $infiles{$familyid}{$sampleid}[$i] =~ /\/?([^\.\/]+\.[^\.]+\.lane(\d+)_[12FfRr])\.\S+/ ) { #Parse 'old' format
		    chomp($1);
		    $Infiles_lane_noending{$familyid}{$sampleid}[$itrack] = $1; #Save new format in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry into hash for both strands and .ending is removed (.fastq).
		    #print $Infiles_lane_noending{$familyid}{$sampleid}[$itrack], "\n";
		    $i++; #Skip second direction
		    $itrack++; #Track for every lane finished
		}
		elsif ( $infiles{$familyid}{$sampleid}[$i] =~ /(\d+_\d+_[^_]+_[^_]+_(index[^_]+)_\d).fastq/ ) { #Parse 'new' format
		    $Infiles_lane_noending{$familyid}{$sampleid}[$itrack] = $1; #Save new format in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry into hash for both strands and .ending is removed (.fastq).
		    $i++; #Skip second direction
		    $itrack++; #Track for every lane finished
		}
	    }
	    $k=1;
	    for (my $i=0;$i<scalar( @ { $infiles{$familyid}{$sampleid} } );$i++) { #Collects inputfiles for every fastq dir and remakes format
	       if ( $infiles{$familyid}{$sampleid}[$i] =~ /\/?([^\.\/]+\.[^\.]+\.lane\d+_[12FfRr])\.fastq/ ) { #Parse 'old' format
	    
	    
		   $Infiles_bothstrands_noending{$familyid}{$sampleid}[$i]= $1; #Save new format in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry per strand and .ending is removed (.fastq).
	       }
	       elsif ( $infiles{$familyid}{$sampleid}[$i] =~ /(\d+_\d+_[^_]+_[^_]+_(index[^_]+)_\d).fastq/ ) { #Parse 'new' format
	    
	    
		   $Infiles_bothstrands_noending{$familyid}{$sampleid}[$i]= $1; #Save new format in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry per strand and .ending is removed (.fastq).
	       }
	    
	    }
	}
    }
    return;
}

sub ReadPedigreeFile {
#Reads famid_pedigree.txt file
#Famid\tSampleID\tSex\tAffected\tMother\tFather\t\Child
#$_[0] = projpath
#$_[1] = familyid
    
    open(PEDF, "<$_[0]/$_[1]/$_[1]_pedigree.txt") or die "Can't open $_[0]/$_[1]/$_[1]_pedigree.txt:$!, \n";    
    
    while (<PEDF>) {
	chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
            next;
        }
	if (m/^\#/) {		# Avoid #
            next;
        }		
	if ( ($_ =~/\S+/) ) {	
	    chomp($_);
	    my @temp = split("\t",$_);	    #Loads pedigree info
	    #print "Temp 1. $temp[1]", "\n";
	    if ( ($temp[0]) && ($temp[1]) ) {
	    #Parse IDN
		if ( $temp[0] =~ /(\d+)-(\d+|-\d+)-(\d+)(A|U)/) { #Match IDN
		    #print "$1-$2-$3$4", "\n"; 
		    $familyMembers{$1}{$temp[0]}=""; #Hash for determining which sampleIDs belong to a certain family
		    #Test for odd and even numbers to determine gender Odd = Male
#% modulous operator, it gives you an integer representation of the remainder of a division operation. If X / Y divides evenly (that is to say, there's no remainder (or modulus)), the result of X % Y will be zero.
		    if ($3 % 2 == 1) { #Male
			$sampleData{$1}{$temp[0]}{'Sex'} = "M";
		    }
		    else { #Female
			$sampleData{$1}{$temp[0]}{'Sex'} = "F";
		    }
		    if ($4 eq "A") { #Affected
			$sampleData{$1}{$temp[0]}{'Affected'} = 1;
		    }
		    elsif ($4 eq "U") { #Unaffected
			$sampleData{$1}{$temp[0]}{'Affected'} = 0;
		    }
		    else {
			$sampleData{$1}{$temp[0]}{'Affected'} = "Na";
		    }
		    if ($temp[1]) {
			$sampleData{$1}{$temp[0]}{'SampleID'} = $temp[1];
		    }
		    else {
			$sampleData{$1}{$temp[0]}{'Affected'} = "Na";
		    }
		    if ($temp[2] && $temp[2] == 1) {
			$sampleData{$1}{$temp[0]}{'Member'} = "Mother";
		    }
		    elsif ($temp[3] && $temp[3] == 1) {
			$sampleData{$1}{$temp[0]}{'Member'} = "Father";
		    }
		    elsif ($temp[4] && $temp[4] == 1) {
			$sampleData{$1}{$temp[0]}{'Member'} = "Child";
		    }
		    else {
			$sampleData{$1}{$temp[0]}{'Member'} = "Na";
		    }
		    if ($temp[5]) {
			$sampleData{$1}{$temp[0]}{'Origin'} = $temp[5];
		    }
		    else {
			$sampleData{$1}{$temp[0]}{'Origin'} = "Na";
		    }
		    if ($temp[6]) {
			$sampleData{$1}{$temp[0]}{'Isolation_kit'} = $temp[6];
		    }
		    else {
			$sampleData{$1}{$temp[0]}{'Isolation_kit'} = "Na";
		    }
		    if ($temp[7]) {
			$sampleData{$1}{$temp[0]}{'Date_of_Isolation'} = $temp[7];
		    }
		    else {
			$sampleData{$1}{$temp[0]}{'Date_of_Isolation'} = "Na";
		    }
		    if ($temp[8]) {
			$sampleData{$1}{$temp[0]}{'Personnel'} = $temp[8];
		    }
		    else {
			$sampleData{$1}{$temp[0]}{'Personnel'} = "Na";
		    }
		    if ($temp[9]) {
			$sampleData{$1}{$temp[0]}{'MD'} = $temp[9];
		    }
		    else {
			$sampleData{$1}{$temp[0]}{'MD'} = "Na";
		    }
		    if ($temp[10]) {
			$sampleData{$1}{$temp[0]}{'Inheritance_Model'} = $temp[10];
		    }
		    else {
			$sampleData{$1}{$temp[0]}{'Inheritance_Model'} = "Na";
		    }
		    if ($temp[11]) {
			$sampleData{$1}{$temp[0]}{'Phenotype_terms'} = $temp[11];
		    }
		    else {
			$sampleData{$1}{$temp[0]}{'Phenotype_terms'} = "Na";
		    }
		    if ($temp[12]) {
			$sampleData{$1}{$temp[0]}{'CMMS_SeqID'} = $temp[12];
		    }
		    else {
			$sampleData{$1}{$temp[0]}{'CMMS_SeqID'} = "Na";
		    }
		    if ($temp[13]) {
			$sampleData{$1}{$temp[0]}{'SciLifeID'} = $temp[13];
		    }
		    else {
			$sampleData{$1}{$temp[0]}{'SciLifeID'} = "Na";
		    }
		    if ($temp[14]) {
			$sampleData{$1}{$temp[0]}{'Capture_kit'} = $temp[14];
		    }
		    else {
			$sampleData{$1}{$temp[0]}{'Capture_kit'} = "Na";
		    }
		    if ($temp[15]) {
			$sampleData{$1}{$temp[0]}{'Capture_date'} = $temp[15];
		    }
		    else {
			$sampleData{$1}{$temp[0]}{'Capture_date'} = "Na";
		    }
		    if ($temp[16]) {
			$sampleData{$1}{$temp[0]}{'Capture_personnel'} = $temp[16];
		    }
		    else {
			$sampleData{$1}{$temp[0]}{'Capture_personnel'} = "Na";
		    }
		    if ($temp[17]) {
			$sampleData{$1}{$temp[0]}{'Clustering_date'} = $temp[17];
		    }
		    else {
			$sampleData{$1}{$temp[0]}{'Clustering_date'} = "Na";
		    }
		    if ($temp[18]) {
			$sampleData{$1}{$temp[0]}{'Sequencing_kit'} = $temp[18];
		    }
		    else {
			$sampleData{$1}{$temp[0]}{'Sequencing_kit'} = "Na";
		    }
		}
	    }
	}
    }
    close(PEDF);
    return;
}

sub ReadFastqc {
#Collects sequence length, Total number of reads, and %GC
#$_[0] = projpath

    my @temp;
    my ($temp_Enc) = "";
    my $pqEnc = q?perl -F'\t' -nae' if ($_=~/Encoding\s+(\S+\s\S+\s\S+\s\S+|\S+\s\S+)/) {print $1;}' ?; #Collect Encoding
    my ($temp_SeqL) = "";
    my $pqSeqL = q?perl -F'\t' -nae' if ($_=~/Sequence length\s(\d+)/) {print $1;}' ?; #Collect Sequence length
    my ($temp_TS,$total_TS) = (0,0);
    my $pqTS = q?perl -F'\t' -nae' if ($_=~/Total Sequences\s(\d+)/) {print $1;}' ?; #Collect Total sequences
    my ($temp_GC,$total_GC,$total_GCc) = (0,0,0);
    my $pqGC = q?perl -F'\t' -nae' if ($_=~/%GC\s(\d+)/) {print $1;}' ?; #Collect GC content
    my ($SeqDup) = "";
    my $pqSeqDup = q?perl -F'\t' -nae' if ($_=~/#Total Duplicate Percentage\s+(\d+.\d)/) {print $1;}' ?; #Collect Sequence duplication level
#Modules
    my $Basic_stat = "";
    my $pqBasic_stat = q?perl -F'\t' -nae' if ($_=~/>>Basic Statistics\s+(\S+)/) {print $1;}' ?; #Collect Basic Statistics
    my $Per_base_seq_qu = "";
    my $pqPer_base_seq_qu = q?perl -F'\t' -nae' if ($_=~/>>Per base sequence quality\s+(\S+)/) {print $1;}' ?; #Collect Per base sequence quality
    my $Per_seq_qu_sco = "";
    my $pqPer_seq_qu_sco = q?perl -F'\t' -nae' if ($_=~/>>Per sequence quality scores\s+(\S+)/) {print $1;}' ?; #Collect Per sequence quality scores
    my $Per_base_seq_co = "";
    my $pqPer_base_seq_co = q?perl -F'\t' -nae' if ($_=~/>>Per base sequence content\s+(\S+)/) {print $1;}' ?; #Collect Per base sequence content
    my $Per_base_gc_co = "";
    my $pqPer_base_gc_co = q?perl -F'\t' -nae' if ($_=~/>>Per base GC content\s+(\S+)/) {print $1;}' ?; #Collect Per base GC content
    my $Per_seq_gc_co = "";
    my $pqPer_seq_gc_co = q?perl -F'\t' -nae' if ($_=~/>>Per sequence GC content\s+(\S+)/) {print $1;}' ?; #Collect Per sequence GC content
    my $Per_base_N_co = "";
    my $pqPer_base_N_co = q?perl -F'\t' -nae' if ($_=~/>>Per base N content\s+(\S+)/) {print $1;}' ?; #Collect Per base N content
    my $Seq_dupl = "";
    my $pqSeq_dupl = q?perl -F'\t' -nae' if ($_=~/>>Sequence Duplication Levels\s+(\S+)/) {print $1;}' ?; #Collect Sequence Duplication Levels
    my $Overrep_seq = "";
    my $pqOverrep_seq = q?perl -F'\t' -nae' if ($_=~/>>Overrepresented sequences\s+(\S+)/) {print $1;}' ?; #Collect Overrepresented sequences
    my $Kmer_co ="";
    my $pqKmer_co = q?perl -F'\t' -nae' if ($_=~/>>Kmer Content\s+(\S+)/) {print $1;}' ?; #Collect Kmer Content

    for my $familyid ( keys %familyMembers ) { #For every family id
	
	for my $sampleid ( keys %{ $familyMembers{$familyid} } ) { #For every sample id
	    
	    if ( $Infiles_bothstrands_noending{$familyid}{$sampleid} ) { #To catch families with no infiles present (solved cases usually)
		for (my $i=0;$i<scalar( @ { $Infiles_bothstrands_noending{$familyid}{$sampleid} } );$i++) {#All infiles
		    
		    #Encoding 
		    $temp_Enc = `$pqEnc $_[0]/$sampleid/fastqc/$Infiles_bothstrands_noending{$familyid}{$sampleid}[$i]_fastqc/fastqc_data.txt;`;
		    #Sequence length
		    $temp_SeqL = `$pqSeqL $_[0]/$sampleid/fastqc/$Infiles_bothstrands_noending{$familyid}{$sampleid}[$i]_fastqc/fastqc_data.txt;`; 
		    #Total Sequences
		    $temp_TS = `$pqTS $_[0]/$sampleid/fastqc/$Infiles_bothstrands_noending{$familyid}{$sampleid}[$i]_fastqc/fastqc_data.txt;`; #Total sequences for one direction
		    #$total_TS = $total_TS+$temp_TS; #Collapse all infiles to one merged file total read count
#Sequence length
		    $SeqDup = `$pqSeqDup $_[0]/$sampleid/fastqc/$Infiles_bothstrands_noending{$familyid}{$sampleid}[$i]_fastqc/fastqc_data.txt;`;
		    #Modules
		    $Basic_stat = `$pqBasic_stat $_[0]/$sampleid/fastqc/$Infiles_bothstrands_noending{$familyid}{$sampleid}[$i]_fastqc/fastqc_data.txt;`;
		    $Per_base_seq_qu = `$pqPer_base_seq_qu $_[0]/$sampleid/fastqc/$Infiles_bothstrands_noending{$familyid}{$sampleid}[$i]_fastqc/fastqc_data.txt;`;
		    $Per_seq_qu_sco = `$pqPer_seq_qu_sco $_[0]/$sampleid/fastqc/$Infiles_bothstrands_noending{$familyid}{$sampleid}[$i]_fastqc/fastqc_data.txt;`;
		    $Per_base_seq_co = `$pqPer_base_seq_co $_[0]/$sampleid/fastqc/$Infiles_bothstrands_noending{$familyid}{$sampleid}[$i]_fastqc/fastqc_data.txt;`;
		    $Per_base_gc_co = `$pqPer_base_gc_co $_[0]/$sampleid/fastqc/$Infiles_bothstrands_noending{$familyid}{$sampleid}[$i]_fastqc/fastqc_data.txt;`;
		    $Per_seq_gc_co = `$pqPer_seq_gc_co $_[0]/$sampleid/fastqc/$Infiles_bothstrands_noending{$familyid}{$sampleid}[$i]_fastqc/fastqc_data.txt;`;
		    $Per_base_N_co = `$pqPer_base_N_co $_[0]/$sampleid/fastqc/$Infiles_bothstrands_noending{$familyid}{$sampleid}[$i]_fastqc/fastqc_data.txt;`;
		    $Seq_dupl = `$pqSeq_dupl $_[0]/$sampleid/fastqc/$Infiles_bothstrands_noending{$familyid}{$sampleid}[$i]_fastqc/fastqc_data.txt;`;
		    $Overrep_seq = `$pqOverrep_seq $_[0]/$sampleid/fastqc/$Infiles_bothstrands_noending{$familyid}{$sampleid}[$i]_fastqc/fastqc_data.txt;`;
		    $Kmer_co = `$pqKmer_co $_[0]/$sampleid/fastqc/$Infiles_bothstrands_noending{$familyid}{$sampleid}[$i]_fastqc/fastqc_data.txt;`;
		    $temp_GC = `$pqGC $_[0]/$sampleid/fastqc/$Infiles_bothstrands_noending{$familyid}{$sampleid}[$i]_fastqc/fastqc_data.txt;`; #Total sequences for one direction
		    if ($temp_Enc) {
			$sampleData{$familyid}{$sampleid}{'Encoding'} .= "$i"."_$temp_Enc;";
		    }
		    else {
			$sampleData{$familyid}{$sampleid}{'Encoding'} .= "$i"."_Na;";
		    }
		    if ($temp_SeqL) {
			$sampleData{$familyid}{$sampleid}{'Sequence_run'} .= "$i"."_2x"."$temp_SeqL;";
		    }
		    else {
			$sampleData{$familyid}{$sampleid}{'Sequence_run'} .= "$i"."_Na;";
		    }
		    if ($temp_TS) {
			$sampleData{$familyid}{$sampleid}{'Nr_of_reads(Total)'} .= "$i"."_$temp_TS;";
		    }
		    else {
			$sampleData{$familyid}{$sampleid}{'Nr_of_reads(Total)'} .= "$i"."_Na";
		    }
		    if ($SeqDup) {
			$sampleData{$familyid}{$sampleid}{'FASTQC_Seq_Dup(%)'} .= "$i"."_$SeqDup;";
			$sampleData{$familyid}{$sampleid}{'FASTQC_Seq_Dup(%)'} =~s/\./\,/; #Substitute . for ,
		    }
		    else {
			$sampleData{$familyid}{$sampleid}{'FASTQC_Seq_Dup(%)'} .= "$i"."_Na;";
		    }
		    if ($Basic_stat) {
			$sampleData{$familyid}{$sampleid}{'Basic_statistics'} .= "$i"."_$Basic_stat;";
		    }
		    else {
			$sampleData{$familyid}{$sampleid}{'Basic_statistics'} .= "$i"."_Na;";
		    }
		    if ($Per_base_seq_qu) {
			$sampleData{$familyid}{$sampleid}{'Per_base_sequence_quality'} .= "$i"."_$Per_base_seq_qu;";
		    }
		    else {
			$sampleData{$familyid}{$sampleid}{'Per_base_sequence_quality'} .= "$i"."_Na;";
		    }
		    if ($Per_seq_qu_sco) {
			$sampleData{$familyid}{$sampleid}{'Per_sequence_quality_scores'} .= "$i"."_$Per_seq_qu_sco;";
		    }
		    else {
			$sampleData{$familyid}{$sampleid}{'Per_sequence_quality_scores'} .= "$i"."_Na;";
		    }
		    if ($Per_base_seq_co) {
			$sampleData{$familyid}{$sampleid}{'Per_base_sequence_content'} .= "$i"."_$Per_base_seq_co;";
		    }
		    else {
			$sampleData{$familyid}{$sampleid}{'Per_base_sequence_content'} .= "$i"."_Na;";
		    }
		    if ($Per_base_gc_co) {
			$sampleData{$familyid}{$sampleid}{'Per_base_GC_content'} .= "$i"."_$Per_base_gc_co;";
		    }
		    else {
			$sampleData{$familyid}{$sampleid}{'Per_base_GC_content'} .= "$i"."_Na;";
		    }
		    if ($Per_seq_gc_co) {
			$sampleData{$familyid}{$sampleid}{'Per_sequence_GC_content'} .= "$i"."_$Per_seq_gc_co;";
		    }
		    else {
			$sampleData{$familyid}{$sampleid}{'Per_sequence_GC_content'} .= "$i"."_Na;";
		    }
		    if ($Per_base_N_co) {
			$sampleData{$familyid}{$sampleid}{'Per_base_N_content'} .= "$i"."_$Per_base_N_co;";
		    }
		    else {
			$sampleData{$familyid}{$sampleid}{'Per_base_N_content'} .= "$i"."_Na;";
		    }
		    if ($Seq_dupl) {
			$sampleData{$familyid}{$sampleid}{'Sequence_Duplication_Levels'} .= "$i"."_$Seq_dupl;";
		    }
		    else {
			$sampleData{$familyid}{$sampleid}{'Sequence_Duplication_Levels'} .= "$i"."_Na;";
		    }
		    if ($Overrep_seq) {
			$sampleData{$familyid}{$sampleid}{'Overrepresented_sequences'} .= "$i"."_$Overrep_seq;";
		    }
		    else {
			$sampleData{$familyid}{$sampleid}{'Overrepresented_sequences'} .= "$i"."_Na;";
		    }
		    if ($Kmer_co) {
			$sampleData{$familyid}{$sampleid}{'Kmer_content'} .= "$i"."_$Kmer_co;";
		    }
		    else {
			$sampleData{$familyid}{$sampleid}{'Kmer_content'} .= "$i"."_Na;";
		    }
		    if($temp_GC) {
			$sampleData{$familyid}{$sampleid}{'%GC'} .= "$i"."_$temp_GC;";
			$sampleData{$familyid}{$sampleid}{'%GC'} =~s/\./\,/; #Substitute . for ,
		    }
		    else {
			$sampleData{$familyid}{$sampleid}{'%GC'} .= "$i"."_Na";
		    }	
		}
	    }
	}
    }
}

sub ReadMosaik {
#Collects alignment statistics
#$_[0] = projpath
    
    my ($temp_MOSVER) = ""; #Mosaik Version
    my $pqMOSVER = q?perl -nae' if ($_=~/(\d+\.\d+\.\d+)\s/) {print $1;}' ?; #Collect Mosaik Version
    my ($temp_UNAM,$total_UNAM) = ("",0);
    my $pqUNAM = q?perl -F'\t' -nae' if ($_=~/# unaligned mates\S+\s+(\d+)\s\(\s+(\d+\.\d+)/) {print $2;}' ?; #Collect Nr of unaligned mates
    my ($temp_FTO,$total_FTO) = ("",0);
    my $pqFTO = q?perl -F'\t' -nae' if ($_=~/# filtered out\S+\s+(\d+)\s\(\s+(\d+\.\d+)/) {print $2;}' ?; #Collect Nr of filtered out reads
    my ($temp_UAL, $total_UAL) = ("",0);
    my $pqUAL = q?perl -F'\t' -nae' if ($_=~/# uniquely aligned mates\S+\s+(\d+)\s\(\s+(\d+\.\d+)/) {print $2;}' ?; #Collect Uniquely aligned mates
    my ($temp_MAL, $total_MAL) = ("",0);
    my $pqMAL = q?perl -F'\t' -nae' if ($_=~/# multiply aligned mates\S+\s+(\d+)\s\(\s+(\d+\.\d+)/) {print $2;}' ?; #Collect Multiply aligned mates
    my ($temp_TS, $total_TS) = ("",0);
    my $pqTS = q?perl -F'\t' -nae' if ($_=~/total aligned:\s+\S+\s+(\S+)\s\(\S+\s(\d+.\d+)/ ) {print $2;} elsif ($_=~/total aligned:\s+(\S+)\s\(\S+\s(\d+.\d+)/ ) { print $2}' ?; #Collect total aligned sequences
    my $fileC = 0;

#Find correct file for every individual within each family
    for my $familyid ( keys %familyMembers ) { #For every family id
	
	for my $sampleid ( keys %{ $familyMembers{$familyid} } ) { #For every sample id
	    my @infiles;
	    if ( $Infiles_lane_noending{$familyid}{$sampleid} ) { #To catch families with no infiles present (solved cases usually)
		for (my $i=0;$i<scalar( @ { $Infiles_lane_noending{$familyid}{$sampleid} } );$i++) {#All infiles
		    
		    my $ret = `ls  $_[0]/$sampleid/mosaik/info/;`; #To look for mosaikAlign files
		    if ($ret =~ /mosaikAlign_$sampleid.$sampleFlowCells{$sampleid}[$i]/) { #Old format
			@infiles = `ls  $_[0]/$sampleid/mosaik/info/mosaikAlign_$sampleid.$sampleFlowCells{$sampleid}[$i].*.stdout.txt;`; #Find all files mosaikAlign files with sampleid and flow-cell.		
		    }
		    unless (@infiles) { #New format
			
			@infiles = `ls  $_[0]/$sampleid/mosaik/info/mosaikAlign_$sampleData{$familyid}{$sampleid}{'Unmerged_post_MosA'}.*.stdout.txt;`; #Find all files mosaikAlign files with sampleid
		    }
		    if (@infiles) {
			chomp(@infiles); #Remove newline
			my $temp = pop(@infiles); #Only last file (should be the most interesting)
			#Mosaik Version
			$temp_MOSVER = `$pqMOSVER $temp;`; #last record i.e. file with highest version (.Y.stdout.txt) which should be the most interesting (however, beaware that older version like .1.0. will not a appear last in @infiles)
			#Unaligned mates		
			$temp_UNAM = `$pqUNAM $temp;`; #last record i.e. file with highest version (.Y.stdout.txt) which should be the most interesting (however, beaware that older version like .1.0. will not a appear last in @infiles)
			if ($temp_UNAM) {
			    $total_UNAM = $total_UNAM+$temp_UNAM; #Sum for average later
			}
			#Nr of sequences filtered out
			$temp_FTO = `$pqFTO $temp;`;
			if ($temp_FTO) {
			    $total_FTO = $total_FTO+$temp_FTO; #Sum for average later
			}
			#Uniquely aligned
			$temp_UAL = `$pqUAL $temp;`;
			if ($temp_UAL) {
			    $total_UAL = $total_UAL+$temp_UAL; #Sum for average later
			}
			#Multiply aligned
			$temp_MAL = `$pqMAL $temp;`; 
			if ($temp_MAL) {
			    $total_MAL = $total_MAL+$temp_MAL; #Sum for average later
			}
			#Total aligned sequences
			$temp_TS = `$pqTS $temp;`; 
			if ($temp_TS) {
			    $total_TS = $total_TS+$temp_TS; #Sum for average later
			}
			$fileC++;
		    }
		}
		if ($temp_MOSVER) {
		$sampleData{$familyid}{$sampleid}{'Mosaik_version'} = $temp_MOSVER;
		}
		else {
		    $sampleData{$familyid}{$sampleid}{'Mosaik_version'} = "Na";
		}
		if ($total_UNAM) {
		$sampleData{$familyid}{$sampleid}{'Unaligned_mates(%)'} = $total_UNAM/$fileC;
	    	$sampleData{$familyid}{$sampleid}{'Unaligned_mates(%)'} =~s/\./\,/; #Substitute . for ,
		}
		else {
		    $sampleData{$familyid}{$sampleid}{'Unaligned_mates(%)'} = "Na";
		}
		if ($total_FTO) {
		    $sampleData{$familyid}{$sampleid}{'Filtered out(%)'} = $total_FTO/$fileC;
		    $sampleData{$familyid}{$sampleid}{'Filtered out(%)'} =~s/\./\,/; #Substitute . for , 
		}
		else {
		    $sampleData{$familyid}{$sampleid}{'Filtered out(%)'} = "Na";
		}
		if ($total_UAL) {
		    $sampleData{$familyid}{$sampleid}{'Uniquely_aligned(%)'} = $total_UAL/$fileC;
		    $sampleData{$familyid}{$sampleid}{'Uniquely_aligned(%)'} =~s/\./\,/; #Substitute . for ,
		}
		else {
		    $sampleData{$familyid}{$sampleid}{'Uniquely_aligned(%)'} = "Na";
		}
		if ($total_MAL) {
		    $sampleData{$familyid}{$sampleid}{'Multiply_aligned(%)'} = $total_MAL/$fileC; 	
		    $sampleData{$familyid}{$sampleid}{'Multiply_aligned(%)'} =~s/\./\,/; #Substitute . for ,
		}
		else {
		    $sampleData{$familyid}{$sampleid}{'Multiply_aligned(%)'} = "Na";
		}
		if ($total_TS) {
		    $sampleData{$familyid}{$sampleid}{'Total_aligned(%)'} = $total_TS/$fileC; 	
		    $sampleData{$familyid}{$sampleid}{'Total_aligned(%)'} =~s/\./\,/; #Substitute . for ,
		}
		else {
		    $sampleData{$familyid}{$sampleid}{'Total_aligned(%)'} = "Na";
		}
#Clear when sampleid infiles are finished
		$total_UNAM= 0; 
		$total_FTO = 0;
		$total_UAL= 0;
		$total_MAL= 0;
		$total_TS= 0;
		$fileC = 0;
	    }
	}
    }
    return;
}

sub MarkDupMetrics {
#Set up which MarlDuplicates metric file to be read and then calls sub routine ReadMDup to collect info. NOTE: Assumes that if no merged file was found that there is only 1 single infile and will report metrics from that. If you have to files you need to merge and run Markduplicate on that bam-file to get correct metrics. 
#LIBRARY	UNPAIRED_READS_EXAMINED	READ_PAIRS_EXAMINED	UNMAPPED_READS	UNPAIRED_READ_DUPLICATES	READ_PAIR_DUPLICATES	READ_PAIR_OPTICAL_DUPLICATES	PERCENT_DUPLICATION	ESTIMATED_LIBRARY_SIZE
#$_[0] = projpath    

    for my $familyid ( keys %familyMembers ) { #For every family id
	
	for my $sampleid ( keys %{ $familyMembers{$familyid} } ) { #For every sample id
	   
	    if ( ($lanes{$sampleid}) || ($Infiles_lane_noending{$familyid}{$sampleid}) ) { 
#Check for any merged files, first
		my $merged_metricfile = "$_[0]/$sampleid/mosaik/$sampleid"."_lanes_";
		for (my $i=0;$i<scalar(@ { $lanes{$sampleid} });$i++) {
		    $merged_metricfile .= $lanes{$sampleid}[$i];
		}
		$merged_metricfile .= "_sorted_merged_pmdmetric";
		if (-e $merged_metricfile) {
		    ReadMDup($merged_metricfile,$familyid,$sampleid );
		}
		else { # Unless merged - perform on singel/individual files
		    
		    for (my $i=0;$i<scalar( @ { $Infiles_lane_noending{$familyid}{$sampleid} } );$i++) {#All infiles  
			my $metricfile = "$_[0]/$sampleid/mosaik/$sampleid.$sampleFlowCells{$sampleid}[$i].@{ $lanes{$sampleid} }_sorted_pmdmetric"; #Old format
			if (-e $metricfile) {
			    ReadMDup($metricfile, $familyid,$sampleid);
			}
			else { #New format	
			    $metricfile = "$_[0]/$sampleid/mosaik/$sampleData{$familyid}{$sampleid}{'Unmerged_post_MosA'}_sorted_pmdmetric";
			    if (-e $metricfile) {
				ReadMDup($metricfile, $familyid,$sampleid);
			    }
			}
		    }
		}
	    }
	}
    }
    return;
}

sub ReadMDup {
#Reads Markduplicate metrics file. Info is on line 8 in file. 
#$_[0] = metricfile
#$_[1] = familyid
#$_[2] = sampleid
    
    open(MDUP, "<$_[0]") or die "Can't open $_[0]:$!, \n";
    
    while (<MDUP>) {
	chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
	    next;
	}
	if (m/^\#/) {		# Avoid #
	    next;
	}		
	if ( ($. ==8) && ($_ =~/(\S+)/) ) { #Only process line 8 and if there is non whitespace there	
	    my @temp = split("\t",$_);	    #Loads MarkDuplicates info
	    
	    if ($temp[4]) {
		$sampleData{$_[1]}{$_[2]}{'UNPAIRED_READ_DUPLICATES'} = $temp[4];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'UNPAIRED_READ_DUPLICATES'} = "Na";
	    }
	    if ($temp[5]) {
		$sampleData{$_[1]}{$_[2]}{'READ_PAIR_DUPLICATES'} = $temp[5];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'READ_PAIR_DUPLICATES'} = "Na";
	    }
	    if ($temp[6]) {
		$sampleData{$_[1]}{$_[2]}{'READ_PAIR_OPTICAL_DUPLICATES'} = $temp[6];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'READ_PAIR_OPTICAL_DUPLICATES'} = "Na";
	    }
	    if ($temp[7]) {
		$sampleData{$_[1]}{$_[2]}{'PERCENT_DUPLICATION'} = $temp[7];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'PERCENT_DUPLICATION'} = "Na";
	    }
	    if ($temp[8]) {
		$sampleData{$_[1]}{$_[2]}{'ESTIMATED_LIBRARY_SIZE'} = $temp[8];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'ESTIMATED_LIBRARY_SIZE'} = "Na";
	    }
	}
    }
    close(MDUP);
    return;	
}	    	    

sub CalHSMetrics {
#Set up which CalculateHS metric file to be read and then calls sub routine ReadCalHS to collect info. NOTE: Assumes that if no merged file was found that there is only 1 single infile and will report metrics from that. If you have to files you need to merge and run CalculateHS on that bam-file to get correct metrics. 
#$_[0] = projpath    

    for my $familyid ( keys %familyMembers ) { #For every family id
	
	for my $sampleid ( keys %{ $familyMembers{$familyid} } ) { #For every sample id

	    if ( ($lanes{$sampleid}) || ($Infiles_lane_noending{$familyid}{$sampleid}) ) { 	    
#Check for any merged files, first
		my $merged_metricfile = "$_[0]/$sampleid/mosaik/coverageReport/$sampleid"."_lanes_";
		for (my $i=0;$i<scalar(@ { $lanes{$sampleid} });$i++) {
		    $merged_metricfile .= $lanes{$sampleid}[$i];
		}
		$merged_metricfile .="_sorted_merged_pmd_CalculateHsMetrics";
		if (-e $merged_metricfile) {
		    ReadCalHS($merged_metricfile,$familyid,$sampleid );
		}
		else { # Unless merged - perform on singel/individual files
		    
		    for (my $i=0;$i<scalar( @ { $Infiles_lane_noending{$familyid}{$sampleid} } );$i++) {#All infiles
			my $metricfile = "$_[0]/$sampleid/mosaik/coverageReport/$sampleid.$sampleFlowCells{$sampleid}[$i].@{ $lanes{$sampleid} }_sorted_pmd_CalculateHsMetrics";	    
#my $metricfile = "$_[0]/$sampleid/mosaik/coverageReport/$sampleid.$sampleData{$familyid}{$sampleid}{'Flow-cell'}.@{ $lanes{$sampleid} }_sorted_pmd_CalculateHsMetrics";
			if (-e $metricfile) {
			    ReadCalHS($metricfile, $familyid,$sampleid);
			}
			else { #New format
			    
			    $metricfile = "$_[0]/$sampleid/mosaik/coverageReport/$sampleData{$familyid}{$sampleid}{'Unmerged_post_MosA'}_sorted_pmd_CalculateHsMetrics";
			    if (-e $metricfile) {
				ReadCalHS($metricfile, $familyid,$sampleid);
			    }
			}
		    }
		}
	    }
	}
    }
    return;
}

sub ReadCalHS {
#Reads CalculateHS metrics file. Info is on line 8 in file. 
#$_[0] = metricfile
#$_[1] = familyid
#$_[2] = sampleid
    
    open(CHS, "<$_[0]") or die "Can't open $_[0]:$!, \n";
    
    while (<CHS>) {
	chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
	    next;
	}
	if (m/^\#/) {		# Avoid #
	    next;
	}		
	if ( ($. ==8) && ($_ =~/(\S+)/) ) { #Only process line 8 and if there is non whitespace there	
	    my @temp = split("\t",$_);	    #Loads CalHs info

	    if ($temp[0]) {
		$sampleData{$_[1]}{$_[2]}{'BAIT_SET'} = $temp[0];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'BAIT_SET'} = "Na";
	    }
	    if ($temp[2]) {
		$sampleData{$_[1]}{$_[2]}{'BAIT_TERRITORY'} = $temp[2];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'BAIT_TERRITORY'} = "Na";
	    }
	    if ($temp[3]) {
		$sampleData{$_[1]}{$_[2]}{'TARGET_TERRITORY'} = $temp[3];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'TARGET_TERRITORY'} = "Na";
	    }
	    if ($temp[7]) {
		$sampleData{$_[1]}{$_[2]}{'PF_UNIQUE_READS'} = $temp[7];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'PF_UNIQUE_READS'} = "Na";
	    }
	    if ($temp[9]) {
		$sampleData{$_[1]}{$_[2]}{'PCT_PF_UQ_READS'} = $temp[9];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'PCT_PF_UQ_READS'} = "Na";
	    }
	    if ($temp[10]) {
		$sampleData{$_[1]}{$_[2]}{'PF_UQ_READS_ALIGNED'} = $temp[10];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'PF_UQ_READS_ALIGNED'} = "Na";
	    }
	    if ($temp[11]) {
		$sampleData{$_[1]}{$_[2]}{'PCT_PF_UQ_READS_ALIGNED'} = $temp[11];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'PCT_PF_UQ_READS_ALIGNED'} = "Na";
	    }
	    if ($temp[12]) {
		$sampleData{$_[1]}{$_[2]}{'PF_UQ_BASES_ALIGNED'} = $temp[12];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'PF_UQ_BASES_ALIGNED'} = "Na";
	    }
	    if ($temp[13]) {
		$sampleData{$_[1]}{$_[2]}{'ON_BAIT_BASES'} = $temp[13];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'ON_BAIT_BASES'} = "Na";
	    }
	    if ($temp[15]) {
		$sampleData{$_[1]}{$_[2]}{'OFF_BAIT_BASES'} = $temp[15];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'OFF_BAIT_BASES'} = "Na";
	    }
	    if ($temp[16]) {
		$sampleData{$_[1]}{$_[2]}{'ON_TARGET_BASES'} = $temp[16];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'ON_TARGET_BASES'} = "Na";
	    }
	    if ($temp[17]) {
		$sampleData{$_[1]}{$_[2]}{'PCT_SELECTED_BASES'} = $temp[17];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'PCT_SELECTED_BASES'} = "Na";
	    }
	    if ($temp[18]) {
		$sampleData{$_[1]}{$_[2]}{'PCT_OFF_BAIT'} = $temp[18];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'PCT_OFF_BAIT'} = "Na";
	    }
	    if ($temp[21]) {
		$sampleData{$_[1]}{$_[2]}{'MEAN_TARGET_COVERAGE'} = $temp[21];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'MEAN_TARGET_COVERAGE'} = "Na";
	    }
	    
	    if ($temp[22]) {
		$sampleData{$_[1]}{$_[2]}{'PCT_USABLE_BASES_ON_BAIT'} = $temp[22];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'PCT_USABLE_BASES_ON_BAIT'} = "Na";
	    }
	    if ($temp[23]) {
		$sampleData{$_[1]}{$_[2]}{'PCT_USABLE_BASES_ON_TARGET'} = $temp[23];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'PCT_USABLE_BASES_ON_TARGET'} = "Na";
	    }
	    if ($temp[24]) {
		$sampleData{$_[1]}{$_[2]}{'FOLD_ENRICHMENT'} = $temp[24];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'FOLD_ENRICHMENT'} = "Na";
	    }
	    if ($temp[26]) {
		$sampleData{$_[1]}{$_[2]}{'FOLD_80_BASE_PENALTY'} = $temp[26];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'FOLD_80_BASE_PENALTY'} = "Na";
	    }
	    if ($temp[28]) {
		$sampleData{$_[1]}{$_[2]}{'PCT_TARGET_BASES_10X'} = $temp[28];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'PCT_TARGET_BASES_10X'} = "Na";
	    }
	    if ($temp[30]) {
		$sampleData{$_[1]}{$_[2]}{'PCT_TARGET_BASES_30X'} = $temp[30];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'PCT_TARGET_BASES_30X'} = "Na";
	    }
	    if ($temp[31]) {
		$sampleData{$_[1]}{$_[2]}{'HS_LIBRARY_SIZE'} = $temp[31];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'HS_LIBRARY_SIZE'} = "Na";
	    }
	    if ($temp[32]) {
		$sampleData{$_[1]}{$_[2]}{'HS_PENALTY_10X'} = $temp[32];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'HS_PENALTY_10X'} = "Na";
	    }
	    if ( ($temp[3]) && ($temp[12]) && ($temp[23]) ) { #TARGET_TERRITORY,PF_UQ_BASES_ALIGNED, PCT_USABLE_BASES_ON_TARGET
		$temp[23] =~s/\,/\./; #Substitute , for ., Perl numeric calculation
		$sampleData{$_[1]}{$_[2]}{'USABLE_MEAN_TARGET_COVERAGE'} = ($temp[12]*$temp[23])/$temp[3];
		$sampleData{$_[1]}{$_[2]}{'USABLE_MEAN_TARGET_COVERAGE'} =~s/\./\,/; #Substitute . for ,
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'USABLE_MEAN_TARGET_COVERAGE'} = "Na";
	    }
	}
    }
    close(CHS);
    return;	
}

sub AlignSumMetrics {
#Set up which CollectAlignmentSummaryMetrics file to be read and then calls sub routine ReadAlSum to collect info. NOTE: Assumes that if no merged file was found that there is only 1 single infile and will report metrics from that. If you have to files you need to merge and run CollectAlignmentSummaryMetrics on that bam-file to get correct metrics. 
#$_[0] = projpath    

    for my $familyid ( keys %familyMembers ) { #For every family id
	
	for my $sampleid ( keys %{ $familyMembers{$familyid} } ) { #For every sample id

	    if ( ($lanes{$sampleid}) || ($Infiles_lane_noending{$familyid}{$sampleid}) ) { 	    
#Check for any merged files, first
		my $merged_metricfile = "$_[0]/$sampleid/mosaik/coverageReport/$sampleid"."_lanes_";
		for (my $i=0;$i<scalar(@ { $lanes{$sampleid} });$i++) {
		    $merged_metricfile .= $lanes{$sampleid}[$i];
		}
		$merged_metricfile .="_sorted_merged_pmd.alignment_summary_metrics";
		if (-e $merged_metricfile) {
		    ReadAlSum($merged_metricfile,$familyid,$sampleid );
		}
		else { # Unless merged - perform on singel/individual files
		    
		    for (my $i=0;$i<scalar( @ { $Infiles_lane_noending{$familyid}{$sampleid} } );$i++) {#All infiles
			my $metricfile = "$_[0]/$sampleid/mosaik/coverageReport/$sampleid.$sampleFlowCells{$sampleid}[$i].@{ $lanes{$sampleid} }_sorted_pmd.alignment_summary_metrics";
			#my $metricfile = "$_[0]/$sampleid/mosaik/coverageReport/$sampleid.$sampleData{$familyid}{$sampleid}{'Flow-cell'}.@{ $lanes{$sampleid} }_sorted_pmd.alignment_summary_metrics";
			if (-e $metricfile) {
			    ReadAlSum($metricfile, $familyid,$sampleid);
			}
			else { #New format	
			    $metricfile = "$_[0]/$sampleid/mosaik/coverageReport/$sampleData{$familyid}{$sampleid}{'Unmerged_post_MosA'}_sorted_pmd.alignment_summary_metrics";
			    if (-e $metricfile) {
				ReadAlSum($metricfile, $familyid,$sampleid);
			    }
			}
		    }
		}
	    }
	}
    }
    return;
}

sub ReadAlSum {
#Reads CollectAlignmentSummaryMetrics file. PAIR info is on line 10 in file. 
#$_[0] = metricfile
#$_[1] = familyid
#$_[2] = sampleid
    
    open(ALS, "<$_[0]") or die "Can't open $_[0]:$!, \n";

    while (<ALS>) {
	chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
	    next;
	}
	if (m/^\#/) {		# Collect reference
	    if ( $_=~/REFERENCE_SEQUENCE=\/(\S+)/ ) { #Whole path
		if( $1=~/([^\/]+$)/ ) { #Just last string i.e. /../(reference.fa)
		    $sampleData{$_[1]}{$_[2]}{'Reference'} = $1;
		}
		else {
		    $sampleData{$_[1]}{$_[2]}{'Reference'} = "Na";
		}
	    }
	    next; #Avoid rest of line
	}		
	if ( ($. ==10) && ($_ =~/(\S+)/) ) { #Only process line 10 (PAIR) and if there is non whitespace there	
	    my @temp = split("\t",$_);	    #Loads CalHs info

	    if ($temp[8]) {
		$sampleData{$_[1]}{$_[2]}{'PF_HQ_ALIGNED_READS'} = $temp[8];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'PF_HQ_ALIGNED_READS'} = "Na";
	    }
	    if ($temp[9]) {
		$sampleData{$_[1]}{$_[2]}{'PF_HQ_ALIGNED_BASES'} = $temp[9];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'PF_HQ_ALIGNED_BASES'} = "Na";
	    }
	    if ($temp[10]) {
		$sampleData{$_[1]}{$_[2]}{'PF_HQ_ALIGNED_Q20_BASES'} = $temp[10];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'PF_HQ_ALIGNED_Q20_BASES'} = "Na";
	    }
	    if ($temp[19]) {
		$sampleData{$_[1]}{$_[2]}{'STRAND_BALANCE'} = $temp[19];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'STRAND_BALANCE'} = "Na";
	    }
	    if ($temp[20]) {
		$sampleData{$_[1]}{$_[2]}{'PCT_CHIMERAS'} = $temp[20];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'PCT_CHIMERAS'} = "Na";
	    }
	    if ($temp[21]) {
		$sampleData{$_[1]}{$_[2]}{'PCT_ADAPTER'} = $temp[21];
	    }
	    else {
		$sampleData{$_[1]}{$_[2]}{'PCT_ADAPTER'} = "Na";
	    }
	}
    }
    close(ALS);
    return;	
}

sub ReadGATK {
#Collects GATK version
#$_[0] = projpath

    my ($temp_GATKVER) = ""; #GATK Version
    my $pqGATKVER = q?perl -nae' if ($_=~/The Genome Analysis Toolkit \(GATK\) (v[^,]+)/) {print $1;last; }' ?; #Collect GATK Version, break if found once in file
     my @infiles;

#Find correct file for every individual within each family
    for my $familyid ( keys %familyMembers ) { #For every family id
	
	for my $sampleid ( keys %{ $familyMembers{$familyid} } ) { #For every sample id

	    if ( $Infiles_lane_noending{$familyid}{$sampleid} ) {
 	    
		for (my $i=0;$i<scalar( @ { $Infiles_lane_noending{$familyid}{$sampleid} } );$i++) {#All infiles
		    
		    my $ret = `ls  $_[0]/$sampleid/mosaik/info/;`; #To look for GATK realign files (Assumes it has been run and same GATK version throughout analyses)
		    if ($ret =~ /gatk_realign_$sampleid/) { #Old format
			@infiles = `ls  $_[0]/$sampleid/mosaik/info/gatk_realign_$sampleid.*.stdout.txt;`; #Find all GATK realign files with sampleid.		
		    }
		    if (@infiles) {
			chomp(@infiles); #Remove newline
			my $temp = pop(@infiles); #Only last file (should be the most interesting)
			#GATK Version		
			$temp_GATKVER = `$pqGATKVER $temp;`; #last record i.e. file with highest version (.Y.stdout.txt) which should be the most interesting
		    }
		}
		if ($temp_GATKVER) {
		    $sampleData{$familyid}{$sampleid}{'GATK_version'} = $temp_GATKVER;
		}
		else {
		    $sampleData{$familyid}{$sampleid}{'GATK_version'} = "Na";
		}
	    }
	}
    }
}

sub GATKVareValExome {
#For exonic variants! Sets up which GATK Variantevaluator file to be read and then calls sub routine ReadVareVal to collect info. NOTE: Assumes that if no merged file was found that there is only 1 single infile and will report metrics from that. If you have to files you need to merge and run GATK variantevaluator on that corresponding vcf-file to get correct metrics. 
#$_[0] = projpath    

    for my $familyid ( keys %familyMembers ) { #For every family id
	
	for my $sampleid ( keys %{ $familyMembers{$familyid} } ) { #For every sample id

	    if ( ($lanes{$sampleid}) || ($Infiles_lane_noending{$familyid}{$sampleid}) ) {	    
#Check for any merged files, first
		my $merged_metricfile = "$_[0]/$sampleid/mosaik/GATK/varianteval/$sampleid"."_lanes_";
		for (my $i=0;$i<scalar(@ { $lanes{$sampleid} });$i++) {
		    $merged_metricfile .= $lanes{$sampleid}[$i];
		}
		#my $merged_BOTH_metricfile ="$merged_metricfile"."_sorted_merged_pmd_allchr_real_recal_resrt_varrecal_BOTH_filt.vcf.varianteval_exonic"; #Test if BOTh INDEL and SNV were analyzed together
		my $merged_BOTH_metricfile ="$merged_metricfile"."_sorted_pmd_rreal_brecal_vrecal_BOTH_exome.vcf.varianteval"; #Test if BOTh INDEL and SNV were analyzed together
		if (-e $merged_BOTH_metricfile) {
		    ReadVareValExome($merged_BOTH_metricfile,$familyid,$sampleid );
		    next;
		}
		$merged_metricfile .="_sorted_merged_pmd_allchr_real_recal_resrt_varrecal_SNV_filt.vcf.varianteval_exonic"; #SNV were analyzed seperatly from INDELs
		if (-e $merged_metricfile) {
		    ReadVareValExome($merged_metricfile,$familyid,$sampleid );
		}
		else { # Unless merged - perform on singel/individual files
		    
		    for (my $i=0;$i<scalar( @ { $Infiles_lane_noending{$familyid}{$sampleid} } );$i++) {#All infiles 
			#my $BOTH_metricfile = "$_[0]/$sampleid/mosaik/GATK/varianteval/$sampleid.$sampleFlowCells{$sampleid}[$i].@{ $lanes{$sampleid} }_sorted_pmd_allchr_real_recal_resrt_varrecal_BOTH_filt.vcf.varianteval_exonic";
			my $BOTH_metricfile = "$_[0]/$sampleid/mosaik/GATK/varianteval/$sampleid.$sampleFlowCells{$sampleid}[$i].@{ $lanes{$sampleid} }_sorted_pmd_rreal_brecal_vrecal_BOTH.vcf.varianteval_exonic";
			if (-e $BOTH_metricfile) {
			    ReadVareValExome($BOTH_metricfile,$familyid,$sampleid );
			    next;
			}
			my $metricfile = "$_[0]/$sampleid/mosaik/GATK/varianteval/$sampleid.$sampleFlowCells{$sampleid}[$i].@{ $lanes{$sampleid} }_sorted_pmd_allchr_real_recal_resrt_varrecal_SNV_filt.vcf.varianteval_exonic";
			if ( -e $metricfile ) {
			    ReadVareValExome($metricfile, $familyid,$sampleid);
			}
			else { #New format
			    if (-e "$_[0]/$sampleid/mosaik/GATK/varianteval/") {	
				$metricfile = "$_[0]/$sampleid/mosaik/GATK/varianteval/$sampleData{$familyid}{$sampleid}{'Unmerged_post_MosA'}_sorted_pmd_allchr_real_recal_resrt_varrecal_BOTH_filt.vcf.varianteval_exonic"; #BOTH
				if (-e $metricfile) {
				    ReadVareValExome($metricfile, $familyid,$sampleid);
				    next;
				}	
				$metricfile = "$_[0]/$sampleid/mosaik/GATK/varianteval/$sampleData{$familyid}{$sampleid}{'Unmerged_post_MosA'}_sorted_pmd_allchr_real_recal_resrt_varrecal_SNV_filt.vcf.varianteval_exonic"; #SNV
				if (-e $metricfile) {
				    ReadVareValExome($metricfile, $familyid,$sampleid);
				}
			    }
			}
		    }
		}
	    }
	}
    }
    return;
}

sub ReadVareValExome {
#Reads GATK Variantevaluation Report file. 
#$_[0] = metricfile
#$_[1] = familyid
#$_[2] = sampleid
    
    open(VV, "<$_[0]") or die "Can't open $_[0]:$!, \n";
    
    while (<VV>) {
	chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
	    next;
	}
	if (m/^\#/) {		# Avoid #
	    next;
	}		
	if ( $_ =~/CompOverlap/ ) { # CompOverlap module	
	    my @temp = split("\t",$_);	    #Loads all info in temp[0] (GATK Variation evaluation file is not really tab-sep)

	    if ($temp[0] =~/\S+\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+.\d+)\s+(\d+)\s+(\d+.\d+)/ ) {
		if ($2 eq "all") {
		    $sampleData{$_[1]}{$_[2]}{'CompRod(Exome)'} = $1;
		    $sampleData{$_[1]}{$_[2]}{'nEvalVariants(Exome)'} = $3;
		    $sampleData{$_[1]}{$_[2]}{'novelSites(Exome)'} = $4;
		    $sampleData{$_[1]}{$_[2]}{'nVariantsAtComp(Exome)'} = $5;
		    $sampleData{$_[1]}{$_[2]}{'compRate(Exome)'} = $6;
		    $sampleData{$_[1]}{$_[2]}{'nConcordant(Exome)'} = $7;
		    $sampleData{$_[1]}{$_[2]}{'concordantRate(Exome)'} = $8;
		    $sampleData{$_[1]}{$_[2]}{'compRate(Exome)'}  =~s/\./\,/; #Substitute . for ,
		    $sampleData{$_[1]}{$_[2]}{'concordantRate(Exome)'} =~s/\./\,/; #Substitute . for ,

		}		
	    }
	}
	if ( $_ =~/CountVariants/ ) { # CountVariants module	
	    my @temp = split("\t",$_);	    #Loads all info in temp[0] (GATK Variation evaluation file is not really tab-sep)
	 
	    if ($temp[0] =~/\S+\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+(\S+)/ ) {
		if ($2 eq "all") {
		    $sampleData{$_[1]}{$_[2]}{'nSNVs(Exome)'} = $3;
		    $sampleData{$_[1]}{$_[2]}{'nInsertions(Exome)'} = $4;
		    $sampleData{$_[1]}{$_[2]}{'nDeletions(Exome)'} = $5;
		    $sampleData{$_[1]}{$_[2]}{'nHets(Exome)'} = $6;
		    $sampleData{$_[1]}{$_[2]}{'nHomVar(Exome)'} = $7;
		    $sampleData{$_[1]}{$_[2]}{'heterozygosity(Exome)'} = $8;
		    $sampleData{$_[1]}{$_[2]}{'heterozygosityPerBp(Exome)'} = $9;
		    $sampleData{$_[1]}{$_[2]}{'hetHomRatio(Exome)'} = $10;	
		    $sampleData{$_[1]}{$_[2]}{'insertionDeletionRatio(Exome)'} = $11;	    
		    $sampleData{$_[1]}{$_[2]}{'heterozygosity(Exome)'} =~s/\./\,/; #Substitute . for ,
		    $sampleData{$_[1]}{$_[2]}{'heterozygosityPerBp(Exome)'} =~s/\./\,/; #Substitute . for ,
		    $sampleData{$_[1]}{$_[2]}{'hetHomRatio(Exome)'} =~s/\./\,/; #Substitute . for ,
		}
		if ($2 eq "novel") {
		    $sampleData{$_[1]}{$_[2]}{'nSNVs(Novel;Exome)'} = $3;
		    $sampleData{$_[1]}{$_[2]}{'nInsertions(Novel;Exome)'} = $4;
		    $sampleData{$_[1]}{$_[2]}{'nDeletions(Novel;Exome)'} = $5;
		    $sampleData{$_[1]}{$_[2]}{'nHets(Novel;Exome)'} = $6;
		    $sampleData{$_[1]}{$_[2]}{'nHomVar(Novel;Exome)'} = $7;
		    $sampleData{$_[1]}{$_[2]}{'heterozygosity(Novel;Exome)'} = $8;
		    $sampleData{$_[1]}{$_[2]}{'heterozygosityPerBp(Novel;Exome)'} = $9;
		    $sampleData{$_[1]}{$_[2]}{'hetHomRatio(Novel;Exome)'} = $10;		    
		    $sampleData{$_[1]}{$_[2]}{'insertionDeletionRatio(Novel;Exome)'} = $11;
		    $sampleData{$_[1]}{$_[2]}{'heterozygosity(Novel;Exome)'} =~s/\./\,/; #Substitute . for ,
		    $sampleData{$_[1]}{$_[2]}{'heterozygosityPerBp(Exome)'} =~s/\./\,/; #Substitute . for ,
		    $sampleData{$_[1]}{$_[2]}{'hetHomRatio(Novel;Exome)'} =~s/\./\,/; #Substitute . for ,
		}		
	    } 
	}  
	if ( $_ =~/TiTvVariantEvaluator/ ) { # TiTvVariantEvaluator module	
	    my @temp = split("\t",$_);	    #Loads all info in temp[0] (GATK Variation evaluation file is not really tab-sep)
	    
	    if ($temp[0] =~/\S+\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+/ ) {
		if ($2 eq "all") {
		    $sampleData{$_[1]}{$_[2]}{'tiTvRatio(Exome)'} = $3;	
		    $sampleData{$_[1]}{$_[2]}{'tiTvRatio(Exome)'}  =~s/\./\,/; #Substitute . for ,	    
		}
		if ($2 eq "novel") {
		    $sampleData{$_[1]}{$_[2]}{'tiTvRatio(Novel;Exome)'} = $3;
		    $sampleData{$_[1]}{$_[2]}{'tiTvRatio(Novel;Exome)'} =~s/\./\,/; #Substitute . for ,		    
		}		
	    }
	}
	if ( $_ =~/IndelSummary/ ) { # IndelSummary module	
	    my @temp = split("\t",$_);	    #Loads all info in temp[0] (GATK Variation evaluation file is not really tab-sep)
	    
	    if ($temp[0] =~/\S+\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)/ ) {
		if ($2 eq "all") {
		    $sampleData{$_[1]}{$_[2]}{'n_indels(Exome)'} = $3;
		    $sampleData{$_[1]}{$_[2]}{'n_singleton_indels(Exome)'} = $4;
		    $sampleData{$_[1]}{$_[2]}{'n_indels_matching_gold_standard(Exome)'} = $5;
		    $sampleData{$_[1]}{$_[2]}{'gold_standard_matching_rate(Exome)'} = $6;
		    $sampleData{$_[1]}{$_[2]}{'SNV_to_indel_ratio(Exome)'} = $7;
		    $sampleData{$_[1]}{$_[2]}{'gold_standard_matching_rate(Exome)'} =~s/\./\,/; #Substitute . for ,
		    $sampleData{$_[1]}{$_[2]}{'SNV_to_indel_ratio(Exome)'} =~s/\./\,/; #Substitute . for ,
		    
		}
		if ($2 eq "novel") {
		    $sampleData{$_[1]}{$_[2]}{'n_indels(Novel;Exome)'} = $3;
		    $sampleData{$_[1]}{$_[2]}{'n_singleton_indels(Novel;Exome)'} = $4;
		    $sampleData{$_[1]}{$_[2]}{'n_indels_matching_gold_standard(Novel;Exome)'} = $5;
		    $sampleData{$_[1]}{$_[2]}{'gold_standard_matching_rate(Novel;Exome)'} = $6;
		    $sampleData{$_[1]}{$_[2]}{'SNV_to_indel_ratio(Novel;Exome)'} = $7;
		    $sampleData{$_[1]}{$_[2]}{'gold_standard_matching_rate(Novel;Exome)'} =~s/\./\,/; #Substitute . for ,
		    $sampleData{$_[1]}{$_[2]}{'SNV_to_indel_ratio(Novel;Exome)'} =~s/\./\,/; #Substitute . for ,	    
		}		
	    }
	}
    }
    close(VV);
    return;	
}

sub GATKVareValAll {
#For all variants! Sets up which GATK Variantevaluator file to be read and then calls sub routine ReadVareVal to collect info. NOTE: Assumes that if no merged file was found that there is only 1 single infile and will report metrics from that. If you have to files you need to merge and run GATK variantevaluator on that corresponding vcf-file to get correct metrics. 
#$_[0] = projpath    

    for my $familyid ( keys %familyMembers ) { #For every family id
	
	for my $sampleid ( keys %{ $familyMembers{$familyid} } ) { #For every sample id

	    if ( ($lanes{$sampleid}) || ($Infiles_lane_noending{$familyid}{$sampleid}) ) {
#Check for any merged files, first
		my $merged_metricfile = "$_[0]/$sampleid/mosaik/GATK/varianteval/$sampleid"."_lanes_";
		for (my $i=0;$i<scalar(@ { $lanes{$sampleid} });$i++) {
		    $merged_metricfile .= $lanes{$sampleid}[$i];
		}
		#my $merged_BOTH_metricfile ="$merged_metricfile"."_sorted_merged_pmd_allchr_real_recal_resrt_varrecal_BOTH_filt.vcf.varianteval"; #Test if BOTH INDEL and SNV were analyzed together
		my $merged_BOTH_metricfile ="$merged_metricfile"."_sorted_pmd_rreal_brecal_vrecal_BOTH.vcf.varianteval"; #Test if BOTH INDEL and SNV were analyzed together
		if (-e $merged_BOTH_metricfile) {
		    ReadVareValAll($merged_BOTH_metricfile,$familyid,$sampleid );
		    next;
		}
		$merged_metricfile .="_sorted_merged_pmd_allchr_real_recal_resrt_varrecal_SNV_filt.vcf.varianteval";	  
		if (-e $merged_metricfile ) {
		    ReadVareValAll($merged_metricfile,$familyid,$sampleid );
		}
		else { # Unless merged - perform on singel/individual files
		    
		    for (my $i=0;$i<scalar( @ { $Infiles_lane_noending{$familyid}{$sampleid} } );$i++) {#All infiles 
			my $BOTH_metricfile = "$_[0]/$sampleid/mosaik/GATK/varianteval/$sampleid.$sampleFlowCells{$sampleid}[$i].@{ $lanes{$sampleid} }_sorted_pmd_allchr_real_recal_resrt_varrecal_BOTH_filt.vcf.varianteval";
			if (-e $BOTH_metricfile) {
			    ReadVareValAll($BOTH_metricfile,$familyid,$sampleid );
			    next;
			}
			my $metricfile = "$_[0]/$sampleid/mosaik/GATK/varianteval/$sampleid.$sampleFlowCells{$sampleid}[$i].@{ $lanes{$sampleid} }_sorted_pmd_allchr_real_recal_resrt_varrecal_SNV_filt.vcf.varianteval";	      
			if ( -e $metricfile ) {
			    ReadVareValAll($metricfile, $familyid,$sampleid);
			}
			else { #New format
			    if (-e "$_[0]/$sampleid/mosaik/GATK/varianteval/") {	
				$metricfile = "$_[0]/$sampleid/mosaik/GATK/varianteval/$sampleData{$familyid}{$sampleid}{'Unmerged_post_MosA'}_sorted_pmd_allchr_real_recal_resrt_varrecal_BOTH_filt.vcf.varianteval"; #BOTH
				if (-e $metricfile) {
				    ReadVareValAll($metricfile, $familyid,$sampleid);
				    next;
				}	
				$metricfile = "$_[0]/$sampleid/mosaik/GATK/varianteval/$sampleData{$familyid}{$sampleid}{'Unmerged_post_MosA'}_sorted_pmd_allchr_real_recal_resrt_varrecal_SNV_filt.vcf.varianteval"; #SNV
				if (-e $metricfile) {
				    ReadVareValAll($metricfile, $familyid,$sampleid);
				}
			    }
			}
		    }
		}
	    }
	}
    }
    return;
}

sub ReadVareValAll {
#Reads GATK Variantevaluation Report file. 
#$_[0] = metricfile
#$_[1] = familyid
#$_[2] = sampleid
    
    open(VV, "<$_[0]") or die "Can't open $_[0]:$!, \n";
    
    while (<VV>) {
	chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
	    next;
	}
	if (m/^\#/) {		# Avoid #
	    next;
	}	
	if ( $_ =~/CompOverlap/ ) { # CompOverlap module	
	    my @temp = split("\t",$_);	    #Loads all info in temp[0] (GATK Variation evaluation file is not really tab-sep)

	    if ($temp[0] =~/\S+\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+.\d+)\s+(\S+)\s+(\d+.\d+)/ ) {
		if ($2 eq "all") {
		    $sampleData{$_[1]}{$_[2]}{'CompRod(All)'} = $1;
		    $sampleData{$_[1]}{$_[2]}{'nEvalVariants(All)'} = $3;
		    $sampleData{$_[1]}{$_[2]}{'novelSites(All)'} = $4;
		    $sampleData{$_[1]}{$_[2]}{'nVariantsAtComp(All)'} = $5;
		    $sampleData{$_[1]}{$_[2]}{'compRate(All)'} = $6;
		    $sampleData{$_[1]}{$_[2]}{'nConcordant(All)'} = $7;
		    $sampleData{$_[1]}{$_[2]}{'concordantRate(All)'} = $8;
		    $sampleData{$_[1]}{$_[2]}{'compRate(All)'} =~s/\./\,/; #Substitute . for ,
		    $sampleData{$_[1]}{$_[2]}{'concordantRate(All)'} =~s/\./\,/; #Substitute . for ,

		}		
	    }
	}
	if ( $_ =~/CountVariants/ ) { # CountVariants module	
	    my @temp = split("\t",$_);	    #Loads all info in temp[0] (GATK Variation evaluation file is not really tab-sep)
	    
	    if ($temp[0] =~/\S+\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+(\S+)/ ) {
		if ($2 eq "all") {
		    $sampleData{$_[1]}{$_[2]}{'nSNVs(All)'} = $3;
		    $sampleData{$_[1]}{$_[2]}{'nInsertions(All)'} = $4;
		    $sampleData{$_[1]}{$_[2]}{'nDeletions(All)'} = $5;
		    $sampleData{$_[1]}{$_[2]}{'nHets(All)'} = $6;
		    $sampleData{$_[1]}{$_[2]}{'nHomVar(All)'} = $7;
		    $sampleData{$_[1]}{$_[2]}{'heterozygosity(All)'} = $8;
		    $sampleData{$_[1]}{$_[2]}{'heterozygosityPerBp(All)'} = $9;
		    $sampleData{$_[1]}{$_[2]}{'hetHomRatio(All)'} = $10;	
		    $sampleData{$_[1]}{$_[2]}{'insertionDeletionRatio(All)'} = $11;	    
		    $sampleData{$_[1]}{$_[2]}{'heterozygosity(All)'} =~s/\./\,/; #Substitute . for ,
		    $sampleData{$_[1]}{$_[2]}{'heterozygosityPerBp(All)'} =~s/\./\,/; #Substitute . for ,
		    $sampleData{$_[1]}{$_[2]}{'hetHomRatio(All)'} =~s/\./\,/; #Substitute . for ,
		}
		if ($2 eq "novel") {
		    $sampleData{$_[1]}{$_[2]}{'nSNVs(Novel;All)'} = $3;
		    $sampleData{$_[1]}{$_[2]}{'nInsertions(Novel;All)'} = $4;
		    $sampleData{$_[1]}{$_[2]}{'nDeletions(Novel;All)'} = $5;
		    $sampleData{$_[1]}{$_[2]}{'nHets(Novel;All)'} = $6;
		    $sampleData{$_[1]}{$_[2]}{'nHomVar(Novel;All)'} = $7;
		    $sampleData{$_[1]}{$_[2]}{'heterozygosity(Novel;All)'} = $8;
		    $sampleData{$_[1]}{$_[2]}{'heterozygosityPerBp(Novel;All)'} = $9;
		    $sampleData{$_[1]}{$_[2]}{'hetHomRatio(Novel;All)'} = $10;		    
		    $sampleData{$_[1]}{$_[2]}{'insertionDeletionRatio(Novel;All)'} = $11;
		    $sampleData{$_[1]}{$_[2]}{'heterozygosity(Novel;All)'} =~s/\./\,/; #Substitute . for ,
		    $sampleData{$_[1]}{$_[2]}{'heterozygosityPerBp(All)'} =~s/\./\,/; #Substitute . for ,
		    $sampleData{$_[1]}{$_[2]}{'hetHomRatio(Novel;All)'} =~s/\./\,/; #Substitute . for ,
		}		
	    }
	}
	if ( $_ =~/TiTvVariantEvaluator/ ) { # TiTvVariantEvaluator module	
	    my @temp = split("\t",$_);	    #Loads all info in temp[0] (GATK Variation evaluation file is not really tab-sep)
	    
	    if ($temp[0] =~/\S+\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+/ ) {
		if ($2 eq "all") {
		    $sampleData{$_[1]}{$_[2]}{'tiTvRatio(All)'} = $3;
		    $sampleData{$_[1]}{$_[2]}{'tiTvRatio(All)'}	=~s/\./\,/; #Substitute . for ,	    
		}
		if ($2 eq "novel") {
		    $sampleData{$_[1]}{$_[2]}{'tiTvRatio(Novel;All)'} = $3;
		    $sampleData{$_[1]}{$_[2]}{'tiTvRatio(Novel;All)'} =~s/\./\,/; #Substitute . for ,		    
		}		
	    }
	}
	if ( $_ =~/IndelSummary/ ) { # IndelSummary module	
	    my @temp = split("\t",$_);	    #Loads all info in temp[0] (GATK Variation evaluation file is not really tab-sep)
	    
	    if ($temp[0] =~/\S+\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)/ ) {
		if ($2 eq "all") {
		    $sampleData{$_[1]}{$_[2]}{'n_indels(All)'} = $3;
		    $sampleData{$_[1]}{$_[2]}{'n_singleton_indels(All)'} = $4;
		    $sampleData{$_[1]}{$_[2]}{'n_indels_matching_gold_standard(All)'} = $5;
		    $sampleData{$_[1]}{$_[2]}{'gold_standard_matching_rate(All)'} = $6;
		    $sampleData{$_[1]}{$_[2]}{'SNV_to_indel_ratio(All)'} = $7;
		    $sampleData{$_[1]}{$_[2]}{'gold_standard_matching_rate(All)'} =~s/\./\,/; #Substitute . for ,
		    $sampleData{$_[1]}{$_[2]}{'SNV_to_indel_ratio(All)'} =~s/\./\,/; #Substitute . for ,
		       
		}
		if ($2 eq "novel") {
		    $sampleData{$_[1]}{$_[2]}{'n_indels(Novel;All)'} = $3;
		    $sampleData{$_[1]}{$_[2]}{'n_singleton_indels(Novel;All)'} = $4;
		    $sampleData{$_[1]}{$_[2]}{'n_indels_matching_gold_standard(Novel;All)'} = $5;
		    $sampleData{$_[1]}{$_[2]}{'gold_standard_matching_rate(Novel;All)'} = $6;
		    $sampleData{$_[1]}{$_[2]}{'SNV_to_indel_ratio(Novel;All)'} = $7;
		    $sampleData{$_[1]}{$_[2]}{'gold_standard_matching_rate(Novel;All)'} =~s/\./\,/; #Substitute . for ,
		    $sampleData{$_[1]}{$_[2]}{'SNV_to_indel_ratio(Novel;All)'} =~s/\./\,/; #Substitute . for ,	    
		}		
	    }
	}	
    }
    close(VV);
    return;	
}	    

sub GenderCheck {
#Uses the coverage to check that the sample sequenced has the expected gender
#$_[0] = projpath    

    my (@temp_cov) = ""; #ChrX coverage
    my $pqcov = q?perl -nae' if ($F[0]=~/^X/ || $F[0]=~/^Y/ ) {print "$F[2],";}' ?; #Collect X and Y coverage. "," required for later split
    my $correct_cov=0;
    my $incorrect_cov=0;
    
    for my $familyid ( keys %familyMembers ) { #For every family id
	
	for my $sampleid ( keys %{ $familyMembers{$familyid} } ) { #For every sample id
	    
	    my @qaCompute_infiles = `ls  $_[0]/$sampleid/mosaik/coverageReport/$sampleid*qaCompute;`; #Collect all qaCompute infiles	    
	    for (my $files=0;$files<scalar(@qaCompute_infiles);$files++) {
		@temp_cov = split(/\,/, `$pqcov $qaCompute_infiles[$files]`);
		print $qaCompute_infiles[$files], "\n";
		if ( ($temp_cov[0]/$temp_cov[1] >= 8) && ($sampleData{$familyid}{$sampleid}{'Sex'} eq "F") ) {
		    $correct_cov++;
		    print $sampleid, "\t", $temp_cov[0]/$temp_cov[1], "\n";
		}
		elsif ( ($temp_cov[0]/$temp_cov[1] < 8) && ($sampleData{$familyid}{$sampleid}{'Sex'} eq "M") ) {
		    $correct_cov++;
		    print  $sampleid, "\t", $temp_cov[0]/$temp_cov[1], "\n";
		}
		else {
		    print $sampleid, "\t", $temp_cov[0]/$temp_cov[1], "\n";
		    $incorrect_cov++;
		}
	    }
	    if ( $correct_cov == scalar(@qaCompute_infiles) ) {
		$sampleData{$familyid}{$sampleid}{'GenderCheck'} = "Pass/($correct_cov/".scalar(@qaCompute_infiles).")";
	    }
	    else {
		$sampleData{$familyid}{$sampleid}{'GenderCheck'} = "Fail/($incorrect_cov/".scalar(@qaCompute_infiles).")";
	    }
	    $correct_cov=0; #Reset for next sampleid
	    $incorrect_cov=0; #Reset for next sampleid
	}	
    }
    return;
}

sub PedigreeCheck {
#Uses the pedigree to check that the sample sequenced has the expected relatives
#$_[0] = projpath    

    my (@sample_order) = ""; #Sample order
    my $psample_order = q?perl -nae 'if ($_=~/^#CHROM/ && $_=~/(\d+)-(\d+|-\d+)-(\d+)(A|U)/) {chomp $_; my @line = split(/\t/,$_); for (my $sample=9;$sample<scalar(@line);$sample++) { print $line[$sample], "\t";}last;}' ?; #Collect sample order from vcf file used to create ".ped", ".map" and hence ".mibs". 
    
    for my $familyid ( keys %familyMembers ) { #For every family id
	my $incorrect_rel=0;
	my $order_file;
	if (-e "$_[0]/$familyid/mosaik/GATK/".$familyid."_sorted_pmd_rreal_brecal_vrecal_BOTH.vcf") {	
	    $order_file = "$_[0]/$familyid/mosaik/GATK/".$familyid."_sorted_pmd_rreal_brecal_vrecal_BOTH.vcf";
	}
	else {
	    $order_file = "$_[0]/$familyid/mosaik/GATK/".$familyid."_allchr_real_recal_resrt_varrecal_SNV_filt.vcf";
	}
	@sample_order = split(/\t/,`$psample_order $order_file`);
	my %family; #Stores family relations and pairwise comparisons family{$IDN}{mother/father/child}{mother/father/child}["Nr"] -> [pairwise]

	if (-e "$_[0]/$familyid/mosaik/samplecheck/$familyid.mibs") {
	    open (PCHECK, "<$_[0]/$familyid/mosaik/samplecheck/$familyid.mibs") or die "Can't open $_[0]/$familyid/mosaik/samplecheck/$familyid.mibs: $!\n"; #Collect prepared family mibs matrix file before pedigree check
	    while (<PCHECK>) {
		chomp $_;
		if (m/^\s+$/) {		# Avoid blank lines
		    next;
		}
		if (m/^\#/) {		# Avoid #
		    next;
		}	
		if ( $_ =~/\S+/ ) { #Pairwise comparisons 
		    
		    my @temp_pairwise = split(/\s/,$_); 
		    for (my $column=0;$column<scalar(@sample_order);$column++) { #All columns in .mibs file
			push ( @{ $family{$sample_order[$.-1]}{ $sampleData{$familyid}{$sample_order[$.-1]}{'Member'} }{ $sampleData{$familyid}{$sample_order[$column]}{'Member'} } }, $temp_pairwise[$column]); #Store sampleID, individual status (Mother/Father/Child) and other family members status (Mother/Father/Child) and each pairwise comparison. Uses array for to accomodate sibling info.
			#print $familyid, "\t", $sample_order[$.-1], "\t", $sampleData{$familyid}{$sample_order[$.-1]}{'Member'}, "\t", $sampleData{$familyid}{$sample_order[$column]}{'Member'}, "\t", $family{$sample_order[$.-1]}{ $sampleData{$familyid}{$sample_order[$.-1]}{'Member'} }{ $sampleData{$familyid}{$sample_order[$column]}{'Member'} }[0], "\n";
		    }
		}
	    }
	    close(PCHECK);
	    for my $IDN ( keys %family ) { #For all IDNs
		
		for my $sampleid ( keys %{$family{$IDN} } ) { #For every relation within family (mother/father/child)
		    
		    for my $members ( keys %{$family{$IDN}{$sampleid}} ) { #For all members including self
			
			for (my $i=0;$i<scalar( @{$family{$IDN}{$sampleid}{$members} } );$i++) { #@ Necessary for siblings
			    
			    if ($family{$IDN}{$sampleid}{$members}[$i] == 1 ) { #Should only hit self
				if ( $sampleid eq  $members) {
				    #print $IDN, "\t", $sampleid,"\t", $members, "\t", $family{$IDN}{$sampleid}{$members}[$i], "\n";
				}
				else {
				    $incorrect_rel++;
				    $sampleData{$familyid}{$IDN}{'PedigreeCheck'} .= "FAIL;Duplicated sample?";
				    #print $IDN, "\t", $sampleid,"\t", $members, "\t", $family{$IDN}{$sampleid}{$members}[$i], "\n";
				}
			    }
			    elsif ($family{$IDN}{$sampleid}{$members}[$i] > 0.74 ) { #Should include parent to child and child to siblings unless inbreed parents
				if ($sampleid eq "Child" || $members eq "Child") { #Correct
				    #print $IDN, "\t", $sampleid,"\t", $members, "\t", $family{$IDN}{$sampleid}{$members}[$i], "\n";
				}
				else {
				    $incorrect_rel++;
				    $sampleData{$familyid}{$IDN}{'PedigreeCheck'} .= "FAIL;Parents related?";
				    #print $IDN, "\t", $sampleid,"\t", $members, "\t", $family{$IDN}{$sampleid}{$members}[$i], "\n";
				}
			    }
			    elsif ($family{$IDN}{$sampleid}{$members}[$i] < 0.70 ) { #Parents unless inbreed
				if ( ($sampleid eq "Father") && ($members eq "Mother") ) {
				    #print $IDN, "\t", $sampleid,"\t", $members, "\t", $family{$IDN}{$sampleid}{$members}[$i], "\n";
				}
				elsif ( ($sampleid eq "Mother") && ($members eq "Father") ) {
				    #print $IDN, "\t", $sampleid,"\t", $members, "\t", $family{$IDN}{$sampleid}{$members}[$i], "\n";
				}
				else {
				    $incorrect_rel++;
				    $sampleData{$familyid}{$IDN}{'PedigreeCheck'} .= "FAIL;$IDN not related to $members;";
				    #print $IDN, "\t", $sampleid,"\t", $members, "\t", $family{$IDN}{$sampleid}{$members}[$i], "\n";
				}
			    }
			}
		    }	
		}
		if ($incorrect_rel == 0) {
		    $sampleData{$familyid}{$IDN}{'PedigreeCheck'} = "PASS";    
		}
	    }
	}
    }
    return;
}

sub InbreedCheck {
#Collects the inbreeding factor calculated by vcfTools --het. 
#$_[0] = projpath     
    
    for my $familyid ( keys %familyMembers ) { #For every family id

	if (-e "$_[0]/$familyid/mosaik/samplecheck/$familyid.het") {
	    open (IBCHECK, "<$_[0]/$familyid/mosaik/samplecheck/$familyid.het") or die "Can't open $_[0]/$familyid/mosaik/samplecheck/$familyid.het: $!\n"; #Collect prepared family mibs matrix file before pedigree check
	    while (<IBCHECK>) {
		chomp $_;
		if (m/^\s+$/) {		# Avoid blank lines
		    next;
		}
		if (m/^\#/) {		# Avoid #
		    next;
		}	
		if ( $_ =~/\S+/ ) { #Pairwise comparisons 
		    
		    my @temp = split(/\s/,$_);
		    if ($temp[4]) { #Column with inbreeding coefficient F
			$sampleData{$familyid}{$temp[0]}{'Inbreeding_F'} = $temp[4];
			$sampleData{$familyid}{$temp[0]}{'Inbreeding_F'} =~s/\./\,/; #Substitute . for ,
		    }
		    else {
			$sampleData{$familyid}{$temp[0]}{'Inbreeding_F'} = "Na";
		    }
		}
	    }
	    close(IBCHECK);
	}
    }
    return;
}


sub WriteSummary {
#Write summary_info output.

    open (WS, ">$_[0]") or die "Can't write to $_[0]: $!\n";
    
#Print headers
    print WS "IDN\t","SampleID\t","Sex\t","Affected\t","Member\t","MD\t","Inheritance_Model\t","Phenotype_terms\t","CMMS_SeqID\t","SciLifeID\t","Capture_kit\t","Capture_date\t","Capture_personnel\t","Clustering_date\t","Sequencing_kit\t",
    "Nr_of_runs\t","Date(s)\t","Flow-Cell(s)\t","Lane(s)\t","Index(s)\t","Reference\t",
    
    "GenderCheck(Pass|Fail(Runs/Total runs per SampleID))\t",

    "PedigreeCheck(Pass|Fail;Reason)\t",

    "Encoding\t","Direction_Sequence_run;\t","Nr_of_reads/run(Direction_reads;)\t","%GC/run(Direction_Seq;)\t","Basic_statistics\t","Per_base_sequence_quality\t","Per_sequence_quality_scores\t","Per_base_sequence_content\t","Per_base_GC_content\t","Per_sequence_GC_content\t","Per_base_N_content\t","Sequence_Duplication_Levels\t","FASTQC_Seq_Dup(%)\t","Overrepresented_sequences\t","Kmer_content\t",
 
    "Mosaik_version\t", "Unaligned_mates(%)(Average/SampleID)\t", "Filtered out(%)(Average/SampleID)\t", "Uniquely_aligned(%)(Average/SampleID)\t","Multiply_aligned(%)(Average/SampleID)\t", "Total_aligned(%;Average/SampleID)\t", 
    "UNPAIRED_READ_DUPLICATES\t", "READ_PAIR_DUPLICATES\t", "READ_PAIR_OPTICAL_DUPLICATES\t", "PERCENT_DUPLICATION\t", "ESTIMATED_LIBRARY_SIZE\t",
    "PF_HQ_ALIGNED_READS\t","PF_HQ_ALIGNED_BASES\t","PF_HQ_ALIGNED_Q20_BASES\t","STRAND_BALANCE\t","PCT_CHIMERAS\t","PCT_ADAPTER\t",  
    "BAIT_SET\t","BAIT_TERRITORY\t", "TARGET_TERRITORY\t","PF_UNIQUE_READS\t","PCT_PF_UQ_READS\t","PF_UQ_READS_ALIGNED\t","ON_BAIT_BASES\t","OFF_BAIT_BASES\t","ON_TARGET_BASES\t","PCT_SELECTED_BASES\t","PCT_OFF_BAIT\t","MEAN_TARGET_COVERAGE\t","USABLE_MEAN_TARGET_COVERAGE\t","PCT_USABLE_BASES_ON_BAIT\t","PCT_USABLE_BASES_ON_TARGET\t","FOLD_ENRICHMENT\t","FOLD_80_BASE_PENALTY\t","PCT_TARGET_BASES_10X\t","PCT_TARGET_BASES_30X\t","HS_LIBRARY_SIZE\t","HS_PENALTY_10X\t",
    
    "GATK_version\t","CompRod(Exome)\t","nEvalVariants(Exome)\t","novelSites(Exome)\t","nVariantsAtComp(Exome)\t","compRate(Exome)\t","nConcordant(Exome)\t","concordantRate(Exome)\t","nSNVs(Exome)\t","nSNVs(Novel;Exome)\t","nInsertions(Exome)\t","nInsertions(Novel;Exome)\t","nDeletions(Exome)\t","nDeletions(Novel;Exome)\t","nHets(Exome)\t","nHets(Novel;Exome)\t","nHomVar(Exome)\t","nHomVar(Novel;Exome)\t","heterozygosity(Exome)\t","heterozygosity(Novel;Exome)\t","heterozygosityPerBp(Exome)\t","heterozygosityPerBp(Novel;Exome)\t","hetHomRatio(Exome)\t","hetHomRatio(Novel;Exome)\t","insertionDeletionRatio(Exome)\t","insertionDeletionRatio(Novel;Exome)\t","tiTvRatio(Exome)\t","tiTvRatio(Novel;Exome)\t","n_indels(Exome)\t","n_indels(Novel;Exome)\t","n_singleton_indels(Exome),\t","n_singleton_indels(Novel;Exome)\t","n_indels_matching_gold_standard(Exome)\t","n_indels_matching_gold_standard(Novel;Exome)\t","gold_standard_matching_rate(Exome)\t","gold_standard_matching_rate(Novel;Exome)\t","SNV_to_indel_ratio(Exome)\t","SNV_to_indel_ratio(Novel;Exome)\t",
    
    "CompRod(All)\t","nEvalVariants(All)\t","novelSites(All)\t","nVariantsAtComp(All)\t","compRate(All)\t","nConcordant(All)\t","concordantRate(All)\t","nSNVs(All)\t","nSNVs(Novel;All)\t","nInsertions(All)\t","nInsertions(Novel;All)\t","nDeletions(All)\t","nDeletions(Novel;All)\t","nHets(All)\t","nHets(Novel;All)\t","nHomVar(All)\t","nHomVar(Novel;All)\t","heterozygosity(All)\t","heterozygosity(Novel;All)\t","heterozygosityPerBp(All)\t","heterozygosityPerBp(Novel;All)\t","hetHomRatio(All)\t","hetHomRatio(Novel;All)\t","insertionDeletionRatio(All)\t","insertionDeletionRatio(Novel;All)\t","tiTvRatio(All)\t","tiTvRatio(Novel;All)\t","n_indels(All)\t","n_indels(Novel;All)\t","n_singleton_indels(All),\t","n_singleton_indels(Novel;All)\t","n_indels_matching_gold_standard(All)\t","n_indels_matching_gold_standard(Novel;All)\t","gold_standard_matching_rate(All)\t","gold_standard_matching_rate(Novel;All)\t","SNV_to_indel_ratio(All)\t","SNV_to_indel_ratio(Novel;All)\t",

    "Inbreeding_Coeff_F\n";
 
    for my $familyid ( keys %familyMembers ) { #For every family id
	
	for my $sampleid ( keys %{ $familyMembers{$familyid} } ) { #For every sample id
	    
	    print WS 
		$sampleid, "\t", #IDN
		$sampleData{$familyid}{$sampleid}{'SampleID'}, "\t", #CMMS SampleID
		$sampleData{$familyid}{$sampleid}{'Sex'}, "\t",
		$sampleData{$familyid}{$sampleid}{'Affected'}, "\t",
		$sampleData{$familyid}{$sampleid}{'Member'}, "\t",
		$sampleData{$familyid}{$sampleid}{'MD'}, "\t",
		$sampleData{$familyid}{$sampleid}{'Inheritance_Model'}, "\t",
		$sampleData{$familyid}{$sampleid}{'Phenotype_terms'}, "\t",
		$sampleData{$familyid}{$sampleid}{'CMMS_SeqID'}, "\t",
		$sampleData{$familyid}{$sampleid}{'SciLifeID'}, "\t",
		$sampleData{$familyid}{$sampleid}{'Capture_kit'}, "\t",
		$sampleData{$familyid}{$sampleid}{'Capture_date'}, "\t",
		$sampleData{$familyid}{$sampleid}{'Capture_personnel'}, "\t",
		$sampleData{$familyid}{$sampleid}{'Clustering_date'}, "\t",
		$sampleData{$familyid}{$sampleid}{'Sequencing_kit'}, "\t";
	    if ( $Infiles_lane_noending{$familyid}{ $sampleid } ) { 
		print WS scalar( @ { $Infiles_lane_noending{$familyid}{ $sampleid } } ), "\t"; #Nr of runs
	    }
	    else {
		print WS "Na\t";
	    }
	    print WS 
		$sampleData{$familyid}{$sampleid}{'Date'}, "\t",
		$sampleData{$familyid}{$sampleid}{'Flow-cell'}, "\t",
		$sampleData{$familyid}{$sampleid}{'Lane'}, "\t",
		$sampleData{$familyid}{$sampleid}{'Index'}, "\t",
		$sampleData{$familyid}{$sampleid}{'Reference'}, "\t",

		$sampleData{$familyid}{$sampleid}{'GenderCheck'}, "\t",
		$sampleData{$familyid}{$sampleid}{'PedigreeCheck'}, "\t",

		$sampleData{$familyid}{$sampleid}{'Encoding'}, "\t",
		$sampleData{$familyid}{$sampleid}{'Sequence_run'}, "\t",
		$sampleData{$familyid}{$sampleid}{'Nr_of_reads(Total)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'%GC'}, "\t",
		$sampleData{$familyid}{$sampleid}{'Basic_statistics'}, "\t",
		$sampleData{$familyid}{$sampleid}{'Per_base_sequence_quality'}, "\t",
		$sampleData{$familyid}{$sampleid}{'Per_sequence_quality_scores'}, "\t",
		$sampleData{$familyid}{$sampleid}{'Per_base_sequence_content'}, "\t",
		$sampleData{$familyid}{$sampleid}{'Per_base_GC_content'}, "\t",
		$sampleData{$familyid}{$sampleid}{'Per_sequence_GC_content'}, "\t",
		$sampleData{$familyid}{$sampleid}{'Per_base_N_content'}, "\t",
		$sampleData{$familyid}{$sampleid}{'Sequence_Duplication_Levels'}, "\t",
		$sampleData{$familyid}{$sampleid}{'FASTQC_Seq_Dup(%)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'Overrepresented_sequences'}, "\t",
		$sampleData{$familyid}{$sampleid}{'Kmer_content'}, "\t",
		
		$sampleData{$familyid}{$sampleid}{'Mosaik_version'}, "\t",
		$sampleData{$familyid}{$sampleid}{'Unaligned_mates(%)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'Filtered out(%)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'Uniquely_aligned(%)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'Multiply_aligned(%)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'Total_aligned(%)'}, "\t",

		$sampleData{$familyid}{$sampleid}{'UNPAIRED_READ_DUPLICATES'}, "\t",
		$sampleData{$familyid}{$sampleid}{'READ_PAIR_DUPLICATES'}, "\t",
		$sampleData{$familyid}{$sampleid}{'READ_PAIR_OPTICAL_DUPLICATES'}, "\t",
		$sampleData{$familyid}{$sampleid}{'PERCENT_DUPLICATION'}, "\t",
		$sampleData{$familyid}{$sampleid}{'ESTIMATED_LIBRARY_SIZE'}, "\t",

		$sampleData{$familyid}{$sampleid}{'PF_HQ_ALIGNED_READS'}, "\t",
		$sampleData{$familyid}{$sampleid}{'PF_HQ_ALIGNED_BASES'}, "\t",
		$sampleData{$familyid}{$sampleid}{'PF_HQ_ALIGNED_Q20_BASES'}, "\t",
		$sampleData{$familyid}{$sampleid}{'STRAND_BALANCE'}, "\t",
		$sampleData{$familyid}{$sampleid}{'PCT_CHIMERAS'}, "\t",
		$sampleData{$familyid}{$sampleid}{'PCT_ADAPTER'}, "\t",

		$sampleData{$familyid}{$sampleid}{'BAIT_SET'}, "\t",
		$sampleData{$familyid}{$sampleid}{'BAIT_TERRITORY'}, "\t",
		$sampleData{$familyid}{$sampleid}{'TARGET_TERRITORY'}, "\t",
		$sampleData{$familyid}{$sampleid}{'PF_UNIQUE_READS'}, "\t",
		$sampleData{$familyid}{$sampleid}{'PCT_PF_UQ_READS'}, "\t",
		$sampleData{$familyid}{$sampleid}{'PF_UQ_READS_ALIGNED'}, "\t",
		$sampleData{$familyid}{$sampleid}{'ON_BAIT_BASES'}, "\t",
		$sampleData{$familyid}{$sampleid}{'OFF_BAIT_BASES'}, "\t",
		$sampleData{$familyid}{$sampleid}{'ON_TARGET_BASES'}, "\t",
		$sampleData{$familyid}{$sampleid}{'PCT_SELECTED_BASES'}, "\t",
		$sampleData{$familyid}{$sampleid}{'PCT_OFF_BAIT'}, "\t",
		$sampleData{$familyid}{$sampleid}{'MEAN_TARGET_COVERAGE'}, "\t",
		$sampleData{$familyid}{$sampleid}{'USABLE_MEAN_TARGET_COVERAGE'}, "\t",
		$sampleData{$familyid}{$sampleid}{'PCT_USABLE_BASES_ON_BAIT'}, "\t",
		$sampleData{$familyid}{$sampleid}{'PCT_USABLE_BASES_ON_TARGET'}, "\t",
		$sampleData{$familyid}{$sampleid}{'FOLD_ENRICHMENT'}, "\t",
		$sampleData{$familyid}{$sampleid}{'FOLD_80_BASE_PENALTY'}, "\t",
		$sampleData{$familyid}{$sampleid}{'PCT_TARGET_BASES_10X'}, "\t",
		$sampleData{$familyid}{$sampleid}{'PCT_TARGET_BASES_30X'}, "\t",
		$sampleData{$familyid}{$sampleid}{'HS_LIBRARY_SIZE'}, "\t",
		$sampleData{$familyid}{$sampleid}{'HS_PENALTY_10X'}, "\t",

		$sampleData{$familyid}{$sampleid}{'GATK_version'}, "\t",
		$sampleData{$familyid}{$sampleid}{'CompRod(Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'nEvalVariants(Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'novelSites(Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'nVariantsAtComp(Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'compRate(Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'nConcordant(Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'concordantRate(Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'nSNVs(Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'nSNVs(Novel;Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'nInsertions(Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'nInsertions(Novel;Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'nDeletions(Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'nDeletions(Novel;Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'nHets(Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'nHets(Novel;Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'nHomVar(Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'nHomVar(Novel;Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'heterozygosity(Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'heterozygosity(Novel;Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'heterozygosityPerBp(Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'heterozygosityPerBp(Novel;Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'hetHomRatio(Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'hetHomRatio(Novel;Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'insertionDeletionRatio(Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'insertionDeletionRatio(Novel;Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'tiTvRatio(Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'tiTvRatio(Novel;Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'n_indels(Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'n_indels(Novel;Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'n_singleton_indels(Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'n_singleton_indels(Novel;Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'n_indels_matching_gold_standard(Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'n_indels_matching_gold_standard(Novel;Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'gold_standard_matching_rate(Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'gold_standard_matching_rate(Novel;Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'SNV_to_indel_ratio(Exome)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'SNV_to_indel_ratio(Novel;Exome)'}, "\t",
		
		$sampleData{$familyid}{$sampleid}{'CompRod(All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'nEvalVariants(All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'novelSites(All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'nVariantsAtComp(All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'compRate(All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'nConcordant(All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'concordantRate(All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'nSNVs(All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'nSNVs(Novel;All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'nInsertions(All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'nInsertions(Novel;All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'nDeletions(All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'nDeletions(Novel;All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'nHets(All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'nHets(Novel;All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'nHomVar(All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'nHomVar(Novel;All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'heterozygosity(All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'heterozygosity(Novel;All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'heterozygosityPerBp(All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'heterozygosityPerBp(Novel;All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'hetHomRatio(All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'hetHomRatio(Novel;All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'insertionDeletionRatio(All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'insertionDeletionRatio(Novel;All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'tiTvRatio(All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'tiTvRatio(Novel;All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'n_indels(All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'n_indels(Novel;All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'n_singleton_indels(All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'n_singleton_indels(Novel;All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'n_indels_matching_gold_standard(All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'n_indels_matching_gold_standard(Novel;All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'gold_standard_matching_rate(All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'gold_standard_matching_rate(Novel;All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'SNV_to_indel_ratio(All)'}, "\t",
		$sampleData{$familyid}{$sampleid}{'SNV_to_indel_ratio(Novel;All)'}, "\t",

		$sampleData{$familyid}{$sampleid}{'Inbreeding_F'}, "\t", "\n";
	}
    }
    close(WS);
    return;
}
