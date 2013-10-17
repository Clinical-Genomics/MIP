# ---------------------------------------------------------------------- 
#  File   :  <filename>.pl
#  History:  <date> (<username>) Created file
# ----------------------------------------------------------------------
#
#  <Short overview/description of what the subroutine does>
# 
# ----------------------------------------------------------------------

sub <sub-name> {
  # Required paramters:
  # ----------------------------------
  # <Explain input parameters>

  # Get input parameters and setup other stuff (dependencies)
  # At least these three are usually always required
  my $sampleID = $_[0];
  my $familyID = $_[1];
  my $aligner = $_[2];
  my $runMode = $_[3];

  # Estimating runtime (hours)
  # ---------------------------
  my $runtimeEst = <time>;

  # -------------------------------------------------------
  #  Writing SBATCH headers
  # -------------------------------------------------------
  ProgramPreRequisites($sampleID, <program-name>, <output-path>, 0, *<HANDLE>, 1, $runtimeEst);

  # -------------------------------------------------------
  #  Figuring out in- and out-files
  # -------------------------------------------------------
  my $baseDir = "<output-path>";  # e.g. /$sampleID/$aligner
  my $inDir = "<where the input files are stored>";
  my $outDir = "<where the output files will be put>";

  my $infileEnding = $sampleInfo{ $familyID }{ $sampleID }{'<upstream program>'}{'fileEnding'};
  my $outfile = "<your outfile name>";  # This is valid for processes without downstream dependencies, otherwise you should also specify file-ending separately

  # Files might have been merged from previous analyses
  my ($infile, $mergeSwitch) = CheckIfMergedFiles($sampleID);

  if ($mergeSwitch == 1) {
    # Files were merged previously
    my $bamPath = "$inDir/$infile"."$infileEnding.<file extension>"
  } else {
    # Files haven't been merged previously
    # You need to loop over and process each of the files!
    for my $infile ($infilesLaneNoEnding{ $sampleID } }) {
      #For all infiles per lane
      #print CHANJO "$inDir/"."$infile.$infileEnding.bam";
    }
  }

  # -------------------------------------------------------
  #  Writing body of the SBATCH script
  # -------------------------------------------------------
  print <handle> "
  # ------------------------------------------------------------
  #  <program-name>
  #  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  <short description of what will happen>
  # ------------------------------------------------------------\n";

  # Drop the file handle
  close(<handle>);

  # If the program was set to run and dry run is disabled
  if ( ($runMode == 1) && ($dryRunAll == 0) ) {
    # Chanjo is a terminally branching job: linear dependencies/no follow up
    FIDSubmitJob($sampleID, $familyID, 2, $parameter{'<program-name>'}{'chain'}, $filename, 0);
  }

  return 1;
}

# Return true
1;
__END__
