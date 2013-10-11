# ---------------------------------------------------------------------- 
#  File   :  chanjo.pl
#  History:  27-sep-2013 (robinandeer) Created file
# ----------------------------------------------------------------------
#
#  Sets up a SBATCH script for running the Chanjo coverage analysis
#  tool that calculates coverage for each genetic element in the CCDS
#  database.
# 
# ----------------------------------------------------------------------

sub chanjo {
  # Generates a SBATCH script for running Chanjo.
  # Required paramters:
  # ----------------------------------
  # <store>:   A location to find/store a SQLite database with genetic elements
  # <bam>:     The location of the BAM-alignment file to read coverage from
  # --cutoff:  The limit for adequate covearage for completeness [default: 10]
  # --sample:  The sample ID of the individual
  # --group:   The family ID for the family
  # --ccds:    Only required if building a new SQLite datastore on the fly
  # --json:    Location to temporarily store the exon annotations (JSON-
  #            formatted, parallelizable)

  # Get input parameters and setup other stuff (dependencies)
  my @sampleIDs = $_[0];
  my $familyID = $_[1];
  my $aligner = $_[2];
  my $outDataDir = $_[3];
  my $storePath = $_[4];
  my $cutoff = $_[5];
  my $runMode = $_[6];
  my $dryRunAll = $_[7];
  my %sampleInfo = $_[8];

  # Estimating runtime (hours)
  # ---------------------------
  # For Chanjo this should be pretty much the same for all BAM-files when all
  # genes are included. The same number of bases are always parsed.
  my $runtimeEst = 1.5 * scalar(@sampleIDs);

  # -------------------------------------------------------
  #  Writing SBATCH headers
  # -------------------------------------------------------
  ProgramPreRequisites($familyID, "chanjo", "$aligner/coverageReport", 0, *CHANJO, 1, $runtimeEst);

  # -------------------------------------------------------
  #  Writing body of the SBATCH script
  # -------------------------------------------------------
  print CHANJO "
  # ------------------------------------------------------------
  #  \033[93mChanjo\033[0m - Coverage analysis tool
  #  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  Creates a temporary JSON file with exon coverage
  #  annotations.
  # ------------------------------------------------------------\n";
  # One central location to store coverage for the family
  my $storePath = "$outDataDir/$familyID/$aligner/coverage.sqlite";

  # Build a new empty database
  print CHANJO "chanjo build $storePath using $ccdsPath \n\n";

  # We run all samples in parallel using the ``&`` and ``wait``-commands
  for my $sampleID (@sampleIDs) {

    # -------------------------------------------------------
    #  Figuring out in- and out-files
    # -------------------------------------------------------
    my $baseDir = "$outDataDir/$sampleID/$aligner";
    my $inDir = $baseDir;
    my $outDir = "$baseDir/coverageReport";

    my $infileEnding = $sampleInfo{ $familyID }{ $sampleID }{'pPicardToolsMarkduplicates'}{'fileEnding'};
    my $outfile = "$sampleID.coverage.json";

    # Files might have been merged from previous analyses
    my ($infile, $mergeSwitch) = CheckIfMergedFiles($sampleID);

    if ($mergeSwitch == 1) {
      # Files were merged previously
      my $bamPath = "$inDir/$infile"."$infileEnding.bam"
    } else {
      # Files haven't been merged previously
      # TODO: Don't know how this case differs from the one above...
      # Q: Is ``$infilesLaneNoEnding`` a global var? Where from?
      for my $infile ($infilesLaneNoEnding{ $sampleID } }) {
        #For all infiles per lane
        #print CHANJO "$inDir/"."$infile.$infileEnding.bam";
      }
    }

    print CHANJO "chanjo annotate $storePath using $bamPath ";
    print CHANJO "--cutoff $cutoff ";
    print CHANJO "--sample $sampleID ";
    print CHANJO "--group $familyID ";
    print CHANJO "--json $outDir/$outfile ";
    print CHANJO "&", "\n\n";

  } # END (for each sample)

  # Now wait for all the sample BAM-files to be processed
  print CHANJO "wait", "\n\n";

  # And sequencially load all the JSON files into the database
  my $jsonPaths = "";
  for my $sampleID (@sampleIDs) {
    $jsonPaths .= "$baseDir/coverageReport$sampleID.coverage.json ";
  }

  print CHANJO "chanjo $storePath import $jsonPaths";
  print CHANJO "wait", "\n\n";

  # Drop the file handle
  close(CHANJO);

  # If the program was set to run and dry run is disabled
  if ( ($runMode == 1) && ($dryRunAll == 0) ) {
    # Chanjo is a terminally branching job: linear dependencies/no follow up
    FIDSubmitJob(0, $familyID, 2, $parameter{'pChanjo'}{'chain'}, $filename, 0);
  }

  return 1;
}

# Return true
1;
__END__

=head1 NAME

Coverage - SBATCH script generator for Chanjo.

=head1 SYNOPSIS

USE CASE: CREATE SCRIPTS FOR EACH FAMILY MEMBER IN A TRIO

  # Subroutines *must* be specifically imported
  use Coverage qw(&chanjo);

  my familyID = 11;
  @sampleIDs = ("11-1-1A", "11-2-1U", "11-2-2U");
  foreach (@sampleIDs) {
    chanjo($_, $familyID, "mosaik");
  }

=head1 DESCRIPTION

This module implements a SBATCH script generator for Chanjo.

=head1 AUTHOR - Your Name

Robin Andeer, robin.andeer@scilifelab.se

=head1 FILES

Chanjo itself will read from a SQLite database with genetic
element definitions. This file is specified through ...

Chanjo also writes a JSON-file to the "coverageReport" folder
of the  "outDataDir" for each individual.

=cut
