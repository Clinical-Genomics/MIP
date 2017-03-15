package Program::Alignment::Gatk;

use strict;
use warnings;
use warnings qw( FATAL utf8 );
use utf8;  #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );

BEGIN {
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.00;

    # Inherit from Exporter to export functions and variables
    our @ISA = qw(Exporter);

    # Functions and variables which are exported by default
    our @EXPORT = qw();

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(realignertargetcreator indelrealigner baserecalibrator printreads);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub realignertargetcreator {

##realignertargetcreator

##Function : Perl wrapper for writing GATK realignertargetcreator recipe to $FILEHANDLE. Based on GATK 3.7.0.
##Returns  : "@commands"
##Arguments: $known_alleles_ref, $intervals_ref, $infile_path, $outfile_path, $referencefile_path, $stderrfile_path, $FILEHANDLE, $downsample_to_coverage, $gatk_disable_auto_index_and_file_lock, $logging_level
##         : $known_alleles_ref                     => Input VCF file(s) with known indels {REF}
##         : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##         : $infile_path                           => Infile paths
##         : $outfile_path                          => Outfile path
##         : $referencefile_path                    => Reference sequence file
##         : $stderrfile_path                       => Stderrfile path
##         : $FILEHANDLE                            => Sbatch filehandle to write to
##         : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##         : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##         : $logging_level                         => Set the minimum level of logging

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $known_alleles_ref;
    my $intervals_ref;
    my $infile_path;
    my $outfile_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $downsample_to_coverage;
    my $gatk_disable_auto_index_and_file_lock;
    my $logging_level;
    
    my $tmpl = {
	known_alleles_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$known_alleles_ref},
	intervals_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$intervals_ref},
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	outfile_path => { required => 1, defined => 1, strict_type => 1, store => \$outfile_path },
	referencefile_path => { required => 1, defined => 1, strict_type => 1, store => \$referencefile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	downsample_to_coverage => { allow => qr/^\d+$/,
				    strict_type => 1, store => \$downsample_to_coverage },
	gatk_disable_auto_index_and_file_lock => { default => 0,
						   allow => [0, 1],
						   strict_type => 1, store => \$gatk_disable_auto_index_and_file_lock },
	logging_level => { default => "INFO",
			   allow => ["INFO", "ERROR", "FATAL"],
			   strict_type => 1, store => \$logging_level },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Gatk realignertargetcreator
    my @commands = qw(--analysis_type RealignerTargetCreator);  #Stores commands depending on input parameters

    ## Options
    if($logging_level) {
	
	push(@commands, "--logging_level ".$logging_level);
    }
    if ($downsample_to_coverage) {

	push(@commands, "--downsample_to_coverage ".$downsample_to_coverage);
    }
    if ($gatk_disable_auto_index_and_file_lock) {

	push(@commands, "--disable_auto_index_creation_and_locking_when_reading_rods");
    }
    if ($intervals_ref) {

	push(@commands, "--intervals ".join(" --intervals ", @$intervals_ref));
    }
    if ($referencefile_path) {

	push(@commands, "--reference_sequence ".$referencefile_path);  #Reference sequence file
    }

    ##Known alleles reference
    push(@commands, "--known ".join(" --known ", @$known_alleles_ref));

    ## Infile
    push(@commands, "--input_file ".$infile_path);

    ## Output
    push(@commands, "--out ".$outfile_path);  #Specify output filename

    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


sub indelrealigner {

##indelrealigner

##Function : Perl wrapper for writing GATK indelrealigner recipe to $FILEHANDLE. Based on GATK 3.7.0.
##Returns  : "@commands"
##Arguments: $known_alleles_ref, $intervals_ref, $infile_path, $outfile_path, $target_intervals_file, $referencefile_path, $stderrfile_path, $FILEHANDLE, $downsample_to_coverage, $gatk_disable_auto_index_and_file_lock, $logging_level, $consensus_determination_model
##         : $known_alleles_ref                     => Input VCF file(s) with known indels {REF}
##         : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##         : $infile_path                           => Infile paths
##         : $outfile_path                          => Outfile path
##         : $target_intervals_file                 => Intervals file output from RealignerTargetCreator
##         : $referencefile_path                    => Reference sequence file
##         : $stderrfile_path                       => Stderrfile path
##         : $FILEHANDLE                            => Sbatch filehandle to write to
##         : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##         : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##         : $logging_level                         => Set the minimum level of logging
##         : $consensus_determination_model         => Determines how to compute the possible alternate consenses

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $known_alleles_ref;
    my $intervals_ref;
    my $infile_path;
    my $outfile_path;
    my $target_intervals_file;
    my $referencefile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $downsample_to_coverage;
    my $gatk_disable_auto_index_and_file_lock;
    my $logging_level;
    my $consensus_determination_model;
    
    my $tmpl = {
	known_alleles_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$known_alleles_ref},
	intervals_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$intervals_ref},
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	outfile_path => { required => 1, defined => 1, strict_type => 1, store => \$outfile_path },
	target_intervals_file => { required => 1, defined => 1, strict_type => 1, store => \$target_intervals_file},
	referencefile_path => { required => 1, defined => 1, strict_type => 1, store => \$referencefile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	downsample_to_coverage => { allow => qr/^\d+$/,
				    strict_type => 1, store => \$downsample_to_coverage },
	gatk_disable_auto_index_and_file_lock => { default => 0,
						   allow => [0, 1],
						   strict_type => 1, store => \$gatk_disable_auto_index_and_file_lock },
	logging_level => { default => "INFO",
			   allow => ["INFO", "ERROR", "FATAL"],
			   strict_type => 1, store => \$logging_level },
	consensus_determination_model => { default => "USE_READS",
					   allow => ["KNOWNS_ONLY", "USE_READS", "USE_SW"],
					   strict_type => 1, store => \$consensus_determination_model },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Gatk indelrealigner
    my @commands = qw(--analysis_type IndelRealigner);  #Stores commands depending on input parameters

    ## Options
    if($logging_level) {
	
	push(@commands, "--logging_level ".$logging_level);
    }
    if ($downsample_to_coverage) {

	push(@commands, "--downsample_to_coverage ".$downsample_to_coverage);
    }
    if ($gatk_disable_auto_index_and_file_lock) {

	push(@commands, "--disable_auto_index_creation_and_locking_when_reading_rods");
    }
    if ($intervals_ref) {

	push(@commands, "--intervals ".join(" --intervals ", @$intervals_ref));
    }
    if ($referencefile_path) {

	push(@commands, "--reference_sequence ".$referencefile_path);  #Reference sequence file
    }
    if ($consensus_determination_model) {
	
	push(@commands, "--consensusDeterminationModel ".$consensus_determination_model);
    }
    if ($target_intervals_file) {
	
	push(@commands, "--targetIntervals ".$target_intervals_file);
    }
    
    ##Known alleles reference
    push(@commands, "--knownAlleles ".join(" --knownAlleles ", @$known_alleles_ref));

    ## Infile
    push(@commands, "--input_file ".$infile_path);

    ## Output
    push(@commands, "--out ".$outfile_path);  #Specify output filename

    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


sub baserecalibrator {

##baserecalibrator

##Function : Perl wrapper for writing GATK baserecalibrator recipe to $FILEHANDLE. Based on GATK 3.7.0.
##Returns  : "@commands"
##Arguments: $known_alleles_ref, $covariates_ref, $intervals_ref, $infile_path, $outfile_path, $referencefile_path, $stderrfile_path, $FILEHANDLE, $num_cpu_threads_per_data_thread, $downsample_to_coverage, $gatk_disable_auto_index_and_file_lock, $logging_level
##         : $known_alleles_ref                     => Input VCF file(s) with known indels {REF}
##         : $covariates_ref                        => One or more covariate to be used in the recalibration. Can be specified multiple times {REF}
##         : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##         : $infile_path                           => Infile paths
##         : $outfile_path                          => Outfile path
##         : $referencefile_path                    => Reference sequence file
##         : $stderrfile_path                       => Stderrfile path
##         : $FILEHANDLE                            => Sbatch filehandle to write to
##         : $num_cpu_threads_per_data_thread       => Number of CPU threads to allocate per data thread
##         : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##         : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##         : $logging_level                         => Set the minimum level of logging

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $known_alleles_ref;
    my $covariates_ref;
    my $intervals_ref;
    my $infile_path;
    my $outfile_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $num_cpu_threads_per_data_thread;
    my $downsample_to_coverage;
    my $gatk_disable_auto_index_and_file_lock;
    my $logging_level;
    
    my $tmpl = {
	known_alleles_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$known_alleles_ref},
	covariates_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$covariates_ref},
	intervals_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$intervals_ref},
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	outfile_path => { required => 1, defined => 1, strict_type => 1, store => \$outfile_path },
	referencefile_path => { required => 1, defined => 1, strict_type => 1, store => \$referencefile_path },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	num_cpu_threads_per_data_thread => { allow => qr/^\d+$/,
				    strict_type => 1, store => \$num_cpu_threads_per_data_thread },
	downsample_to_coverage => { allow => qr/^\d+$/,
				    strict_type => 1, store => \$downsample_to_coverage },
	gatk_disable_auto_index_and_file_lock => { default => 0,
						   allow => [0, 1],
						   strict_type => 1, store => \$gatk_disable_auto_index_and_file_lock },
	logging_level => { default => "INFO",
			   allow => ["INFO", "ERROR", "FATAL"],
			   strict_type => 1, store => \$logging_level },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Gatk baserecalibrator
    my @commands = qw(--analysis_type BaseRecalibrator);  #Stores commands depending on input parameters

    ## Options
    if($logging_level) {
	
	push(@commands, "--logging_level ".$logging_level);
    }
    if ($downsample_to_coverage) {

	push(@commands, "--downsample_to_coverage ".$downsample_to_coverage);
    }
    if ($gatk_disable_auto_index_and_file_lock) {

	push(@commands, "--disable_auto_index_creation_and_locking_when_reading_rods");
    }
    if ($num_cpu_threads_per_data_thread) {

	push(@commands, "--num_cpu_threads_per_data_thread ".$num_cpu_threads_per_data_thread);
    }
    if ($intervals_ref) {

	push(@commands, "--intervals ".join(" --intervals ", @$intervals_ref));
    }
    if ($referencefile_path) {

	push(@commands, "--reference_sequence ".$referencefile_path);  #Reference sequence file
    }

    ##Known alleles reference
    push(@commands, "--knownSites ".join(" --knownSites ", @$known_alleles_ref));

    ##Covariates
    push(@commands, "--covariate ".join(" --covariate ", @$covariates_ref));

    ## Infile
    push(@commands, "--input_file ".$infile_path);

    ## Output
    push(@commands, "--out ".$outfile_path);  #Specify output filename

    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


sub printreads {

##printreads

##Function : Perl wrapper for writing GATK printreads recipe to $FILEHANDLE. Based on GATK 3.7.0.
##Returns  : "@commands"
##Arguments: $read_filters_ref, $static_quantized_quals_ref, $intervals_ref, $infile_path, $outfile_path, $referencefile_path, $base_quality_score_recalibration_file, $stderrfile_path, $FILEHANDLE, $num_cpu_threads_per_data_thread, $downsample_to_coverage, $gatk_disable_auto_index_and_file_lock, $disable_indel_qual, $logging_level
##         : $read_filters_ref                      => Filters to apply to reads before analysis {REF}
##         : $static_quantized_quals_ref            => Use static quantized quality scores to a given number of levels (with -BQSR) {REF}
##         : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##         : $infile_path                           => Infile paths
##         : $outfile_path                          => Outfile path
##         : $referencefile_path                    => Reference sequence file
##         : $base_quality_score_recalibration_file => Input covariates table file for on-the-fly base quality score recalibration
##         : $stderrfile_path                       => Stderrfile path
##         : $FILEHANDLE                            => Sbatch filehandle to write to
##         : $num_cpu_threads_per_data_thread       => Number of CPU threads to allocate per data thread
##         : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##         : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##         : $disable_indel_qual                    => Disable printing of base insertion and deletion tags (with -BQSR)
##         : $logging_level                         => Set the minimum level of logging

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $read_filters_ref;
    my $static_quantized_quals_ref;
    my $intervals_ref;
    my $infile_path;
    my $outfile_path;
    my $referencefile_path;
    my $base_quality_score_recalibration_file;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $num_cpu_threads_per_data_thread;
    my $downsample_to_coverage;
    my $gatk_disable_auto_index_and_file_lock;
    my $disable_indel_qual;
    my $logging_level;
    
    my $tmpl = {
	read_filters_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$read_filters_ref},
	static_quantized_quals_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$static_quantized_quals_ref},
	intervals_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$intervals_ref},
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	outfile_path => { required => 1, defined => 1, strict_type => 1, store => \$outfile_path },
	referencefile_path => { required => 1, defined => 1, strict_type => 1, store => \$referencefile_path },
	base_quality_score_recalibration_file => { strict_type => 1, store => \$base_quality_score_recalibration_file },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	num_cpu_threads_per_data_thread => { allow => qr/^\d+$/,
				    strict_type => 1, store => \$num_cpu_threads_per_data_thread },
	downsample_to_coverage => { allow => qr/^\d+$/,
				    strict_type => 1, store => \$downsample_to_coverage },
	gatk_disable_auto_index_and_file_lock => { default => 0,
						   allow => [0, 1],
						   strict_type => 1, store => \$gatk_disable_auto_index_and_file_lock },
	disable_indel_qual => { default => 0,
				    allow => [0, 1],
				    strict_type => 1, store => \$disable_indel_qual },
	logging_level => { default => "INFO",
			   allow => ["INFO", "ERROR", "FATAL"],
			   strict_type => 1, store => \$logging_level },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ## Gatk printreads
    my @commands = qw(--analysis_type PrintReads);  #Stores commands depending on input parameters

    ## Options
    if($logging_level) {
	
	push(@commands, "--logging_level ".$logging_level);
    }
    if ($downsample_to_coverage) {

	push(@commands, "--downsample_to_coverage ".$downsample_to_coverage);
    }
    if ($gatk_disable_auto_index_and_file_lock) {

	push(@commands, "--disable_auto_index_creation_and_locking_when_reading_rods");
    }
    if ($num_cpu_threads_per_data_thread) {

	push(@commands, "--num_cpu_threads_per_data_thread ".$num_cpu_threads_per_data_thread);
    }
    if ($read_filters_ref) {

	push(@commands, "--read_filter ".join(" --read_filter ", @$read_filters_ref));
    }
    if ($intervals_ref) {

	push(@commands, "--intervals ".join(" --intervals ", @$intervals_ref));
    }
    if ($referencefile_path) {

	push(@commands, "--reference_sequence ".$referencefile_path);  #Reference sequence file
    }
    if ($base_quality_score_recalibration_file) {

	push(@commands, "--BQSR ".$base_quality_score_recalibration_file);
    }
    if ($disable_indel_qual) {

	push(@commands, "--disable_indel_quals");
    }
    if ($static_quantized_quals_ref) {
	
	push(@commands, "--static_quantized_quals ".join(" --static_quantized_quals ", @$static_quantized_quals_ref));
    }

    ## Infile
    push(@commands, "--input_file ".$infile_path);

    ## Output
    push(@commands, "--out ".$outfile_path);  #Specify output filename

    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


1;
