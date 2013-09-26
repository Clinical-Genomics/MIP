package Coverage;

# Followed: http://www.perlmonks.org/?node_id=102347
use Exporter;
use strict;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 1.00;
@ISA         = qw(Exporter);
@EXPORT      = ();             # Auto export (not recommended)
@EXPORT_OK   = qw(chanjo);     # Export on request

sub chanjo {
  # Generates a SBATCH script for running Chanjo.
  # Required paramters:
  # --------------------
  # <store>: A location to find/store a SQLite database with genetic elements
  # <bam>: The location of the BAM-alignment file to read coverage from
  # --cutoff: The limit for adequate covearage for completeness [default: 10]
  # --sample: The sample ID of the individual
  # --group: The family ID for the family
  # --ccds: Only required if building a new SQLite datastore on the fly
  # --json: Location to temporarily store the exon annotations (JSON-formatted)

  # Getting the function parameters
  # First make some local variables to avoid polluting global scope 
  local($sampleID, $famliyID, $aligner, $cutoff, $jsonPath, $storePath,
        $runtimeEstimate);

  # Get input parameters and setup other stuff
  $sampleID = $_[0];
  $familyID = $_[1];
  $aligner = $_[2];
  $cutoff = $scriptParameter{'chanjoCutoff'};

  $jsonPath = "$scriptParameter{'outDataDir'}/$sampleID/$aligner/"
              . "coverageReport/$sampleID.coverage.json";

  # This is the SQLite database used as a reference for all samples
  $storePath = "$scriptParameter{'chanjoSQL'}";

  # Estimating runtime (min)
  # -------------------------
  # For Chanjo this should be pretty much the same for all BAM-files when all
  # genes are included. The same number of bases are always parsed.
  $runtimeEstimate = 120;

  # -------------------------------------------------------
  #  Writing SBATCH headers
  #  ~~~~~~~~~~~~~~~~~~~~~~~~
  #  Hard to remember how this worked...
  # -------------------------------------------------------
  ProgramPreRequisites($sampleID, "Chanjo", $aligner."/coverageReport", 0, *RCOVP, 1, 1);

  # -------------------------------------------------------
  #  Writing body of the SBATCH script
  # -------------------------------------------------------
  print CHANJO "
  # ============================================================
  #  Create a temp JSON file with exon coverage annotations
  # ------------------------------------------------------------\n";
  print CHANJO "chanjo annotate $storePath using $bamFile";
  print CHANJO "--cutoff $cutoff";
  print CHANJO "--sample $sampleID";
  print CHANJO "--group $familyID";
  print CHANJO "--json $jsonPath";
}

# Return true
1;