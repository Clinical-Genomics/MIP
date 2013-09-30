# ---------------------------------------------------------------------- 
#  File   :  Coverage.pm
#  History:  27-sep-2013 (robinandeer) Created module
# ----------------------------------------------------------------------
#
#  Sets up a SBATCH script for running the Chanjo coverage analysis
#  tool that calculates coverage for each genetic element in the CCDS
#  database.
# 
# ----------------------------------------------------------------------
use strict;
use warnings;

package Coverage;

# Followed: http://www.perlmonks.org/?node_id=102347
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK);

$VERSION     = 1.00;
@ISA         = qw(Exporter);
@EXPORT      = ();             # Auto export (not recommended)
@EXPORT_OK   = qw(chanjo);     # Export on request

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

  # Get input parameters and setup other stuff
  my $sampleID = $_[0];
  my $familyID = $_[1];
  my $aligner = $_[2];
  my $cutoff = $scriptParameter{'chanjoCutoff'};

  my $jsonPath = "$scriptParameter{'outDataDir'}/$sampleID/$aligner/"
                 . "coverageReport/$sampleID.coverage.json";

  # This is the SQLite database used as a reference for all samples
  my $storePath = "$scriptParameter{'chanjoSQL'}";

  # Estimating runtime (min)
  # -------------------------
  # For Chanjo this should be pretty much the same for all BAM-files when all
  # genes are included. The same number of bases are always parsed.
  my $runtimeEst = 120;

  # -------------------------------------------------------
  #  Writing SBATCH headers
  #  ~~~~~~~~~~~~~~~~~~~~~~~~
  #  Hard to remember how this worked...
  # -------------------------------------------------------
  ProgramPreRequisites($sampleID, "chanjo", "$aligner/coverageReport", 0, *CHANJO, 1, 1);

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
