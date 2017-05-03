package Program::Variantcalling::Gatk;

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
    our @EXPORT_OK = qw(catvariants genotypegvcfs selectvariants variantrecalibrator applyrecalibration calculategenotypeposteriors combinevariants);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case


sub genotypegvcfs {

##genotypegvcfs

##Function : Perl wrapper for writing GATK genotypegvcfs recipe to $FILEHANDLE. Based on GATK 3.7.0.
##Returns  : "@commands"
##Arguments: $gatk_path, $intervals_ref, $infile_paths_ref, $outfile_path, $referencefile_path, $stderrfile_path, $FILEHANDLE, $pedigree, $dbsnp, $downsample_to_coverage, $gatk_disable_auto_index_and_file_lock, $logging_level, $pedigree_validation_type, $include_nonvariant_sites
##         : $gatk_path                             => Path to java jar and analysis to run
##         : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##         : $infile_paths_ref                      => Infile paths {REF}
##         : $outfile_path                          => Outfile path
##         : $referencefile_path                    => Reference sequence file
##         : $stderrfile_path                       => Stderrfile path
##         : $FILEHANDLE                            => Sbatch filehandle to write to
##         : $pedigree                              => Pedigree files for samples
##         : $dbsnp                                 => DbSNP file
##         : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##         : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##         : $logging_level                         => Set the minimum level of logging
##         : $pedigree_validation_type              => Validation strictness for pedigree
##         : $include_nonvariant_sites              => Include loci found to be non-variant after genotyping


    my ($arg_href) = @_;

    ## Default(s)
    my $gatk_disable_auto_index_and_file_lock;
    my $logging_level;
    my $pedigree_validation_type;
    my $include_nonvariant_sites;

    ## Flatten argument(s)
    my $intervals_ref;
    my $infile_paths_ref;
    my $outfile_path;
    my $referencefile_path;
    my $gatk_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $pedigree;
    my $dbsnp;
    my $downsample_to_coverage;
    
    my $tmpl = {
	intervals_ref => { default => [], strict_type => 1, store => \$intervals_ref},
	infile_paths_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$infile_paths_ref},
	outfile_path => { required => 1, defined => 1, strict_type => 1, store => \$outfile_path },
	referencefile_path => { required => 1, defined => 1, strict_type => 1, store => \$referencefile_path },
	gatk_path => { strict_type => 1, store => \$gatk_path},
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	pedigree => { strict_type => 1, store => \$pedigree },
	dbsnp => { strict_type => 1, store => \$dbsnp },
	downsample_to_coverage => { allow => qr/^\d+$/,
				    strict_type => 1, store => \$downsample_to_coverage },
	gatk_disable_auto_index_and_file_lock => { default => 0,
						   allow => [0, 1],
						   strict_type => 1, store => \$gatk_disable_auto_index_and_file_lock },
	logging_level => { default => "INFO",
			   allow => ["INFO", "ERROR", "FATAL"],
			   strict_type => 1, store => \$logging_level },
	pedigree_validation_type => { default => "SILENT",
				      allow => ["STRICT", "SILENT"],
				      strict_type => 1, store => \$pedigree_validation_type },
	include_nonvariant_sites  => { default => 0,
				       allow => [0, 1],
				       strict_type => 1, store => \$include_nonvariant_sites },

    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ### Gatk genotypegvcfs

    ## Common gatk options
    my @commands = qw(--analysis_type GenotypeGVCFs);  #stores commands depending on input parameters

    if ($gatk_path) {

	push(@commands, $gatk_path);
    }

    push(@commands, "--logging_level ".$logging_level);

    push(@commands, "--pedigreeValidationType ".$pedigree_validation_type);

    if ($pedigree) {

	push(@commands, "--pedigree ".$pedigree);
    }
    if ($downsample_to_coverage) {

	push(@commands, "--downsample_to_coverage ".$downsample_to_coverage);
    }
    if ($gatk_disable_auto_index_and_file_lock) {

	push(@commands, "--disable_auto_index_creation_and_locking_when_reading_rods");
    }
    if (@$intervals_ref) {

	push(@commands, "--intervals ".join(" --intervals ", @$intervals_ref));
    }
    if ($referencefile_path) {

	push(@commands, "--reference_sequence ".$referencefile_path);  #Reference sequence file
    }

    ## Tool specific options
    if ($include_nonvariant_sites) {
	
	push(@commands, "--includeNonVariantSites");
    }
    if ($dbsnp) {
	
	push(@commands, "--dbsnp ".$dbsnp);
    }

    ## Infile
    if (@$infile_paths_ref) {

	push(@commands, "--variant ".join(" --variant ", @$infile_paths_ref));
    }

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


sub selectvariants {

##selectvariants

##Function : Perl wrapper for writing GATK selectvariants recipe to $FILEHANDLE. Based on GATK 3.7.0.
##Returns  : "@commands"
##Arguments: $intervals_ref, $sample_names_ref, $infile_path, $outfile_path, $referencefile_path, $stderrfile_path, $FILEHANDLE, $gatk_path, $pedigree, $downsample_to_coverage, $gatk_disable_auto_index_and_file_lock, $logging_level, $pedigree_validation_type, $exclude_nonvariants
##         : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##         : $sample_names_ref                      => Include genotypes from this sample {REF}
##         : $infile_path                           => Infile paths
##         : $outfile_path                          => Outfile path
##         : $referencefile_path                    => Reference sequence file
##         : $stderrfile_path                       => Stderrfile path
##         : $FILEHANDLE                            => Sbatch filehandle to write to
##         : $gatk_path                             => Path to java jar and analysis to run
##         : $pedigree                              => Pedigree files for samples
##         : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##         : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##         : $logging_level                         => Set the minimum level of logging
##         : $pedigree_validation_type              => Validation strictness for pedigree
##         : $exclude_nonvariants                   => Exclude non-variant sites


    my ($arg_href) = @_;

    ## Default(s)
    my $gatk_disable_auto_index_and_file_lock;
    my $logging_level;
    my $pedigree_validation_type;
    my $exclude_nonvariants;

    ## Flatten argument(s)
    my $intervals_ref;
    my $sample_names_ref;
    my $infile_path;
    my $outfile_path;
    my $referencefile_path;
    my $gatk_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $pedigree;
    my $downsample_to_coverage;
    
    my $tmpl = {
	intervals_ref => { default => [], strict_type => 1, store => \$intervals_ref},
	sample_names_ref  => { required => 1, defined => 1, default => [], strict_type => 1, store => \$sample_names_ref },
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	outfile_path => { required => 1, defined => 1, strict_type => 1, store => \$outfile_path },
	referencefile_path => { required => 1, defined => 1, strict_type => 1, store => \$referencefile_path },
	gatk_path => { strict_type => 1, store => \$gatk_path},
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	pedigree => { strict_type => 1, store => \$pedigree },
	downsample_to_coverage => { allow => qr/^\d+$/,
				    strict_type => 1, store => \$downsample_to_coverage },
	gatk_disable_auto_index_and_file_lock => { default => 0,
						   allow => [0, 1],
						   strict_type => 1, store => \$gatk_disable_auto_index_and_file_lock },
	logging_level => { default => "INFO",
			   allow => ["INFO", "ERROR", "FATAL"],
			   strict_type => 1, store => \$logging_level },
	pedigree_validation_type => { default => "SILENT",
				      allow => ["STRICT", "SILENT"],
				      strict_type => 1, store => \$pedigree_validation_type },
	exclude_nonvariants => { default => 0,
				 allow => [0, 1],
				 strict_type => 1, store => \$exclude_nonvariants },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ### Gatk selectvariants

    ## Common gatk options
    my @commands = qw(--analysis_type SelectVariants);  #stores commands depending on input parameters

    if ($gatk_path) {

	push(@commands, $gatk_path);
    }

    push(@commands, "--logging_level ".$logging_level);

    push(@commands, "--pedigreeValidationType ".$pedigree_validation_type);

    if ($pedigree) {

	push(@commands, "--pedigree ".$pedigree);
    }
    if ($downsample_to_coverage) {

	push(@commands, "--downsample_to_coverage ".$downsample_to_coverage);
    }
    if ($gatk_disable_auto_index_and_file_lock) {

	push(@commands, "--disable_auto_index_creation_and_locking_when_reading_rods");
    }
    if (@$intervals_ref) {

	push(@commands, "--intervals ".join(" --intervals ", @$intervals_ref));
    }
    if ($referencefile_path) {

	push(@commands, "--reference_sequence ".$referencefile_path);  #Reference sequence file
    }

    ## Tool specific options
    if ($exclude_nonvariants) {

	push(@commands, "--excludeNonVariants");
    }
    if (@$sample_names_ref) {
	
	push(@commands, "--sample_name ".join(" --sample_name ", @$sample_names_ref));
    }

    ## Infile
    if ($infile_path) {

	push(@commands, "--variant ".$infile_path);
    }

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


sub catvariants {

##catvariants

##Function : Perl wrapper for writing GATK catvariants recipe to $FILEHANDLE. Based on GATK 3.7.0.
##Returns  : "@commands"
##Arguments: $gatk_path, $intervals_ref, $infile_paths_ref, $outfile_path, $referencefile_path, $stderrfile_path, $FILEHANDLE, $downsample_to_coverage, $gatk_disable_auto_index_and_file_lock, $logging_level, $assume_sorted
##         : $gatk_path                             => Path to java jar and analysis to run
##         : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##         : $infile_paths_ref                      => Infile paths {REF}
##         : $outfile_path                          => Outfile path
##         : $referencefile_path                    => Reference sequence file
##         : $stderrfile_path                       => Stderrfile path
##         : $FILEHANDLE                            => Sbatch filehandle to write to
##         : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##         : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##         : $logging_level                         => Set the minimum level of logging
##         : $assume_sorted                         => Assume_sorted should be true if the input files are already sorted

    my ($arg_href) = @_;

    ## Default(s)
    my $gatk_disable_auto_index_and_file_lock;
    my $logging_level;
    my $assume_sorted;

    ## Flatten argument(s)
    my $gatk_path;
    my $intervals_ref;
    my $infile_paths_ref;
    my $outfile_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $downsample_to_coverage;
    
    my $tmpl = {
	gatk_path => { required => 1, defined => 1, strict_type => 1, store => \$gatk_path},
	intervals_ref => { default => [], strict_type => 1, store => \$intervals_ref},
	infile_paths_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$infile_paths_ref},
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
	assume_sorted => { default => 1,
			   allow => [undef, 0, 1],
			   strict_type => 1, store => \$assume_sorted },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ### Gatk catvariants
    my @commands = qw(-cp);  #Stores commands depending on input parameters

    push(@commands, $gatk_path);

    ## Common gatk options
    push(@commands, "--logging_level ".$logging_level);

    if ($downsample_to_coverage) {

	push(@commands, "--downsample_to_coverage ".$downsample_to_coverage);
    }
    if ($gatk_disable_auto_index_and_file_lock) {

	push(@commands, "--disable_auto_index_creation_and_locking_when_reading_rods");
    }
    if (@$intervals_ref) {

	push(@commands, "--intervals ".join(" --intervals ", @$intervals_ref));
    }
    if ($referencefile_path) {

	push(@commands, "--reference_sequence ".$referencefile_path);  #Reference sequence file
    }

    ## Tool specific options
    if ($assume_sorted) {

	push(@commands, "--assumeSorted");  #AssumeSorted should be true if the input files are already sorted
    }

    ## Infile
    if (@$infile_paths_ref) {

	push(@commands, "--variant: ".join(" --variant: ", @$infile_paths_ref));
    }

    ## Output
    push(@commands, "--outputFile ".$outfile_path);  #Specify output filename

    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


sub variantrecalibrator {

##variantrecalibrator

##Function : Perl wrapper for writing GATK variantrecalibrator recipe to $FILEHANDLE. Based on GATK 3.7.0.
##Returns  : "@commands"
##Arguments: $infile_paths_ref, $resources_ref, $intervals_ref, $read_filters_ref, $static_quantized_quals_ref, $annotations_ref, $referencefile_path, $base_quality_score_recalibration_file, $stderrfile_path, $FILEHANDLE, $pedigree, $rscript_file_path, $recal_file_path, $tranches_file_path, $num_cpu_threads_per_data_thread, $downsample_to_coverage, $max_gaussian_level, $gatk_disable_auto_index_and_file_lock, $disable_indel_qual, $logging_level, $pedigree_validation_type, $mode
##         : $infile_paths_ref                      => Infile paths
##         : $resources_ref                         => A list of sites for which to apply a prior probability of being correct but which aren't used by the algorithm (training and truth setsare required to run)
##         : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##         : $read_filters_ref                      => Filters to apply to reads before analysis {REF}
##         : $static_quantized_quals_ref            => Use static quantized quality scores to a given number of levels (with -BQSR) {REF}
##         : $annotations_ref                       => One or more specific annotations to apply to variant calls
##         : $referencefile_path                    => Reference sequence file
##         : $mode                                  => Mode for emitting reference confidence scores
##         : $base_quality_score_recalibration_file => Input covariates table file for on-the-fly base quality score recalibration
##         : $stderrfile_path                       => Stderrfile path
##         : $FILEHANDLE                            => Sbatch filehandle to write to
##         : $pedigree                              => Pedigree files for samples
##         : $rscript_file_path                     => Rscript file path
##         : $recal_file_path                       => The output recal file used by ApplyRecalibration
##         : $tranches_file_path                    => The output tranches file used by ApplyRecalibration
##         : $num_cpu_threads_per_data_thread       => Number of CPU threads to allocate per data thread
##         : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##         : $max_gaussian_level                    => Max number of Gaussians for the positive model
##         : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##         : $disable_indel_qual                    => Disable printing of base insertion and deletion tags (with -BQSR)
##         : $logging_level                         => Set the minimum level of logging
##         : $pedigree_validation_type              => Validation strictness for pedigree

    my ($arg_href) = @_;

    ## Default(s)
    my $gatk_disable_auto_index_and_file_lock;
    my $disable_indel_qual;
    my $logging_level;
    my $pedigree_validation_type;

    ## Flatten argument(s)
    my $infile_paths_ref;
    my $resources_ref;
    my $intervals_ref;
    my $read_filters_ref;
    my $static_quantized_quals_ref;
    my $annotations_ref;
    my $referencefile_path;
    my $mode;
    my $base_quality_score_recalibration_file;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $pedigree;
    my $rscript_file_path;
    my $recal_file_path;
    my $tranches_file_path;
    my $num_cpu_threads_per_data_thread;
    my $downsample_to_coverage;
    my $max_gaussian_level;
    
    my $tmpl = {
	infile_paths_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$infile_paths_ref},
	resources_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$resources_ref },
	intervals_ref => { default => [], strict_type => 1, store => \$intervals_ref},
	read_filters_ref => { default => [], strict_type => 1, store => \$read_filters_ref},
	static_quantized_quals_ref => { default => [], strict_type => 1, store => \$static_quantized_quals_ref},
	annotations_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$annotations_ref},
	referencefile_path => { required => 1, defined => 1, strict_type => 1, store => \$referencefile_path },
	mode => { required => 1, defined => 1,
		  allow => ["SNP", "INDEL", "BOTH"],
		  strict_type => 1, store => \$mode },
	base_quality_score_recalibration_file => { strict_type => 1, store => \$base_quality_score_recalibration_file },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	pedigree => { strict_type => 1, store => \$pedigree },
	rscript_file_path => { strict_type => 1, store => \$rscript_file_path },
	recal_file_path => { strict_type => 1, store => \$recal_file_path },
	tranches_file_path => { strict_type => 1, store => \$tranches_file_path },
	num_cpu_threads_per_data_thread => { allow => qr/^\d+$/,
				    strict_type => 1, store => \$num_cpu_threads_per_data_thread },
	downsample_to_coverage => { allow => qr/^\d+$/,
				    strict_type => 1, store => \$downsample_to_coverage },
	max_gaussian_level => { allow => [undef, qr/^\d+$/],
				strict_type => 1, store => \$max_gaussian_level },
	gatk_disable_auto_index_and_file_lock => { default => 0,
						   allow => [0, 1],
						   strict_type => 1, store => \$gatk_disable_auto_index_and_file_lock },
	disable_indel_qual => { default => 0,
				allow => [0, 1],
				strict_type => 1, store => \$disable_indel_qual },
	logging_level => { default => "INFO",
			   allow => ["INFO", "ERROR", "FATAL"],
			   strict_type => 1, store => \$logging_level },
	pedigree_validation_type => { default => "SILENT",
				      allow => ["STRICT", "SILENT"],
				      strict_type => 1, store => \$pedigree_validation_type },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ### Gatk variantrecalibrator

    ## Common gatk options
    my @commands = qw(--analysis_type VariantRecalibrator);  #Stores commands depending on input parameters

    push(@commands, "--logging_level ".$logging_level);

    push(@commands, "--pedigreeValidationType ".$pedigree_validation_type);

    if ($pedigree) {

	push(@commands, "--pedigree ".$pedigree);
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
    if (@$read_filters_ref) {

	push(@commands, "--read_filter ".join(" --read_filter ", @$read_filters_ref));
    }
    if (@$intervals_ref) {

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
    if (@$static_quantized_quals_ref) {
	
	push(@commands, "--static_quantized_quals ".join(" --static_quantized_quals ", @$static_quantized_quals_ref));
    }

    ## Tool specific options 
    if ($rscript_file_path) {
	
	push(@commands, "--rscript_file ".$rscript_file_path);
    }
    if ($recal_file_path) {
	
	push(@commands, "--recal_file ".$recal_file_path);
    }
    if ($tranches_file_path) {
	
	push(@commands, "--tranches_file ".$tranches_file_path);
    }
    if ($max_gaussian_level) {
	
	push(@commands, "--maxGaussians ".$max_gaussian_level);
    }	
    if (@$annotations_ref) {

	push(@commands, "--use_annotation ".join(" --use_annotation ", @$annotations_ref));
    }
    if (@$resources_ref) {
	
	push(@commands, "--resource:".join(" --resource:", @$resources_ref));
    }
    
    push(@commands, "--mode ".$mode);

    ## Infile
    push(@commands, "--input ".join("--input_file ", @$infile_paths_ref));

    if ($stderrfile_path) {

	push(@commands, "2> ".$stderrfile_path);  #Redirect stderr output to program specific stderr file
    }
    if($FILEHANDLE) {
	
	print $FILEHANDLE join(" ", @commands)." ";
    }
    return @commands;
}


sub applyrecalibration {

##applyrecalibration

##Function : Perl wrapper for writing GATK applyrecalibration recipe to $FILEHANDLE. Based on GATK 3.7.0.
##Returns  : "@commands"
##Arguments: $infile_path, $resources_ref, $intervals_ref, $read_filters_ref, $static_quantized_quals_ref, $annotations_ref, $outfile_path, $referencefile_path, $base_quality_score_recalibration_file, $stderrfile_path, $FILEHANDLE, $pedigree, $ts_filter_level, $recal_file_path, $tranches_file_path, $num_cpu_threads_per_data_thread, $downsample_to_coverage, $max_gaussian_level, $gatk_disable_auto_index_and_file_lock, $disable_indel_qual, $logging_level, $pedigree_validation_type, $mode
##         : $infile_path                      => Infile paths
##         : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##         : $read_filters_ref                      => Filters to apply to reads before analysis {REF}
##         : $static_quantized_quals_ref            => Use static quantized quality scores to a given number of levels (with -BQSR) {REF}
##         : $outfile_path                          => Outfile path
##         : $referencefile_path                    => Reference sequence file
##         : $mode                                  => Mode for emitting reference confidence scores
##         : $base_quality_score_recalibration_file => Input covariates table file for on-the-fly base quality score recalibration
##         : $stderrfile_path                       => Stderrfile path
##         : $FILEHANDLE                            => Sbatch filehandle to write to
##         : $pedigree                              => Pedigree files for samples
##         : $recal_file_path                       => The output recal file used by ApplyRecalibration
##         : $tranches_file_path                    => The output tranches file used by ApplyRecalibration
##         : $num_cpu_threads_per_data_thread       => Number of CPU threads to allocate per data thread
##         : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##         : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##         : $disable_indel_qual                    => Disable printing of base insertion and deletion tags (with -BQSR)
##         : $logging_level                         => Set the minimum level of logging
##         : $pedigree_validation_type              => Validation strictness for pedigree
##         : $ts_filter_level                     => Rscript file path

    my ($arg_href) = @_;

    ## Default(s)
    my $gatk_disable_auto_index_and_file_lock;
    my $disable_indel_qual;
    my $logging_level;
    my $pedigree_validation_type;

    ## Flatten argument(s)
    my $infile_path;
    my $intervals_ref;
    my $read_filters_ref;
    my $static_quantized_quals_ref;
    my $outfile_path;
    my $referencefile_path;
    my $mode;
    my $base_quality_score_recalibration_file;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $pedigree;
    my $recal_file_path;
    my $tranches_file_path;
    my $num_cpu_threads_per_data_thread;
    my $downsample_to_coverage;
    my $ts_filter_level;
    
    my $tmpl = {
	intervals_ref => { default => [], strict_type => 1, store => \$intervals_ref },
	read_filters_ref => { default => [], strict_type => 1, store => \$read_filters_ref },
	static_quantized_quals_ref => { default => [], strict_type => 1, store => \$static_quantized_quals_ref },
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path },
	outfile_path => { required => 1, defined => 1, strict_type => 1, store => \$outfile_path },
	referencefile_path => { required => 1, defined => 1, strict_type => 1, store => \$referencefile_path },
	mode => { required => 1, defined => 1,
		  allow => ["SNP", "INDEL", "BOTH"],
		  strict_type => 1, store => \$mode },
	base_quality_score_recalibration_file => { strict_type => 1, store => \$base_quality_score_recalibration_file },
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	pedigree => { strict_type => 1, store => \$pedigree },
	tranches_file_path => { strict_type => 1, store => \$tranches_file_path },
	recal_file_path => { strict_type => 1, store => \$recal_file_path },
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
	pedigree_validation_type => { default => "SILENT",
				      allow => ["STRICT", "SILENT"],
				      strict_type => 1, store => \$pedigree_validation_type },
	ts_filter_level => { default => 99.0,
			     allow => qr/^\d+$|^\d+.\d+$/,
			     strict_type => 1, store => \$ts_filter_level },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ### Gatk applyrecalibration

    ## Common gatk options
    my @commands = qw(--analysis_type ApplyRecalibration);  #Stores commands depending on input parameters

    push(@commands, "--logging_level ".$logging_level);

    push(@commands, "--pedigreeValidationType ".$pedigree_validation_type);

    if ($pedigree) {

	push(@commands, "--pedigree ".$pedigree);
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
    if (@$read_filters_ref) {

	push(@commands, "--read_filter ".join(" --read_filter ", @$read_filters_ref));
    }
    if (@$intervals_ref) {

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
    if (@$static_quantized_quals_ref) {
	
	push(@commands, "--static_quantized_quals ".join(" --static_quantized_quals ", @$static_quantized_quals_ref));
    }

    ## Tool specific options 
    if ($ts_filter_level) {
	
	push(@commands, "--ts_filter_level ".$ts_filter_level);
    }
    if ($recal_file_path) {
	
	push(@commands, "--recal_file ".$recal_file_path);
    }
    if ($tranches_file_path) {
	
	push(@commands, "--tranches_file ".$tranches_file_path);
    }	
    
    push(@commands, "--mode ".$mode);

    ## Infile
    push(@commands, "--input ".$infile_path);

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


sub calculategenotypeposteriors {

##calculategenotypeposteriors

##Function : Perl wrapper for writing GATK calculategenotypeposteriors recipe to $FILEHANDLE. Based on GATK 3.7.0.
##Returns  : "@commands"
##Arguments: $intervals_ref, $infile_path, $outfile_path, $referencefile_path, $stderrfile_path, $FILEHANDLE, $gatk_path, $pedigree, $supporting_callset_file_path, $downsample_to_coverage, $gatk_disable_auto_index_and_file_lock, $logging_level, $pedigree_validation_type
##         : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##         : $infile_path                           => Infile paths
##         : $outfile_path                          => Outfile path
##         : $referencefile_path                    => Reference sequence file
##         : $stderrfile_path                       => Stderrfile path
##         : $FILEHANDLE                            => Sbatch filehandle to write to
##         : $gatk_path                             => Path to java jar and analysis to run
##         : $pedigree                              => Pedigree files for samples
##         : $supporting_callset_file_path          => Other callsets to use in generating genotype posteriors
##         : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##         : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##         : $logging_level                         => Set the minimum level of logging
##         : $pedigree_validation_type              => Validation strictness for pedigree


    my ($arg_href) = @_;

    ## Default(s)
    my $gatk_disable_auto_index_and_file_lock;
    my $logging_level;
    my $pedigree_validation_type;

    ## Flatten argument(s)
    my $intervals_ref;
    my $infile_path;
    my $outfile_path;
    my $referencefile_path;
    my $gatk_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $pedigree;
    my $supporting_callset_file_path;
    my $downsample_to_coverage;
    
    my $tmpl = {
	intervals_ref => { default => [], strict_type => 1, store => \$intervals_ref},
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	outfile_path => { required => 1, defined => 1, strict_type => 1, store => \$outfile_path },
	referencefile_path => { required => 1, defined => 1, strict_type => 1, store => \$referencefile_path },
	gatk_path => { strict_type => 1, store => \$gatk_path},
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	pedigree => { strict_type => 1, store => \$pedigree },
	supporting_callset_file_path => { strict_type => 1, store => \$supporting_callset_file_path },
	downsample_to_coverage => { allow => qr/^\d+$/,
				    strict_type => 1, store => \$downsample_to_coverage },
	gatk_disable_auto_index_and_file_lock => { default => 0,
						   allow => [0, 1],
						   strict_type => 1, store => \$gatk_disable_auto_index_and_file_lock },
	logging_level => { default => "INFO",
			   allow => ["INFO", "ERROR", "FATAL"],
			   strict_type => 1, store => \$logging_level },
	pedigree_validation_type => { default => "SILENT",
				      allow => ["STRICT", "SILENT"],
				      strict_type => 1, store => \$pedigree_validation_type },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ### Gatk calculategenotypeposteriors

    ## Common gatk options
    my @commands = qw(--analysis_type CalculateGenotypePosteriors);  #stores commands depending on input parameters

    if ($gatk_path) {

	push(@commands, $gatk_path);
    }

    push(@commands, "--logging_level ".$logging_level);

    push(@commands, "--pedigreeValidationType ".$pedigree_validation_type);

    if ($pedigree) {

	push(@commands, "--pedigree ".$pedigree);
    }
    if ($downsample_to_coverage) {

	push(@commands, "--downsample_to_coverage ".$downsample_to_coverage);
    }
    if ($gatk_disable_auto_index_and_file_lock) {

	push(@commands, "--disable_auto_index_creation_and_locking_when_reading_rods");
    }
    if (@$intervals_ref) {

	push(@commands, "--intervals ".join(" --intervals ", @$intervals_ref));
    }
    if ($referencefile_path) {

	push(@commands, "--reference_sequence ".$referencefile_path);  #Reference sequence file
    }

    ## Tool specific options
    if ($supporting_callset_file_path) {

	push(@commands, "--supporting ".$supporting_callset_file_path);
    }

    ## Infile
    if ($infile_path) {

	push(@commands, "--variant ".$infile_path);
    }

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


sub combinevariants {

##combinevariants

##Function : Perl wrapper for writing GATK combinevariants recipe to $FILEHANDLE. Based on GATK 3.7.0.
##Returns  : "@commands"
##Arguments: $gatk_path, $intervals_ref, $infile_paths_ref, $outfile_path, $referencefile_path, $stderrfile_path, $FILEHANDLE, $pedigree, $prioritize_caller, $downsample_to_coverage, $gatk_disable_auto_index_and_file_lock, $logging_level, $pedigree_validation_type, $exclude_nonvariants, $genotype_merge_option,
##         : $gatk_path                             => Path to java jar and analysis to run
##         : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##         : $infile_paths_ref                      => Infile paths {REF}
##         : $outfile_path                          => Outfile path
##         : $referencefile_path                    => Reference sequence file
##         : $stderrfile_path                       => Stderrfile path
##         : $FILEHANDLE                            => Sbatch filehandle to write to
##         : $pedigree                              => Pedigree files for samples
##         : $prioritize_caller                     => Comma seperated string specifying priority for merging
##         : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##         : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##         : $logging_level                         => Set the minimum level of logging
##         : $pedigree_validation_type              => Validation strictness for pedigree
##         : $exclude_nonvariants                   => Exclude non-variant sites
##         : $genotype_merge_option                 => Determines how we should merge genotype records for samples shared across the ROD files


    my ($arg_href) = @_;

    ## Default(s)
    my $gatk_disable_auto_index_and_file_lock;
    my $logging_level;
    my $pedigree_validation_type;
    my $exclude_nonvariants;

    ## Flatten argument(s)
    my $intervals_ref;
    my $infile_paths_ref;
    my $outfile_path;
    my $referencefile_path;
    my $gatk_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $pedigree;
    my $prioritize_caller;
    my $genotype_merge_option;
    my $downsample_to_coverage;
    
    my $tmpl = {
	intervals_ref => { default => [], strict_type => 1, store => \$intervals_ref},
	infile_paths_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$infile_paths_ref},
	outfile_path => { required => 1, defined => 1, strict_type => 1, store => \$outfile_path },
	referencefile_path => { required => 1, defined => 1, strict_type => 1, store => \$referencefile_path },
	gatk_path => { strict_type => 1, store => \$gatk_path},
	stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
	FILEHANDLE => { store => \$FILEHANDLE },
	pedigree => { strict_type => 1, store => \$pedigree },
	prioritize_caller => { strict_type => 1, store => \$prioritize_caller },
	downsample_to_coverage => { allow => qr/^\d+$/,
				    strict_type => 1, store => \$downsample_to_coverage },
	gatk_disable_auto_index_and_file_lock => { default => 0,
						   allow => [0, 1],
						   strict_type => 1, store => \$gatk_disable_auto_index_and_file_lock },
	logging_level => { default => "INFO",
			   allow => ["INFO", "ERROR", "FATAL"],
			   strict_type => 1, store => \$logging_level },
	pedigree_validation_type => { default => "SILENT",
				      allow => ["STRICT", "SILENT"],
				      strict_type => 1, store => \$pedigree_validation_type },
	exclude_nonvariants  => { default => 0,
				  allow => [0, 1],
				  strict_type => 1, store => \$exclude_nonvariants },
	genotype_merge_option => { default => "PRIORITIZE",
				   allow => [undef, "UNIQUIFY", "PRIORITIZE", "UNSORTED", "REQUIRE_UNIQUE"],
				   strict_type => 1, store => \$genotype_merge_option },
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    ### Gatk combinevariants

    ## Common gatk options
    my @commands = qw(--analysis_type CombineVariants);  #stores commands depending on input parameters

    if ($gatk_path) {

	push(@commands, $gatk_path);
    }

    push(@commands, "--logging_level ".$logging_level);

    push(@commands, "--pedigreeValidationType ".$pedigree_validation_type);

    if ($pedigree) {

	push(@commands, "--pedigree ".$pedigree);
    }
    if ($downsample_to_coverage) {

	push(@commands, "--downsample_to_coverage ".$downsample_to_coverage);
    }
    if ($gatk_disable_auto_index_and_file_lock) {

	push(@commands, "--disable_auto_index_creation_and_locking_when_reading_rods");
    }
    if (@$intervals_ref) {

	push(@commands, "--intervals ".join(" --intervals ", @$intervals_ref));
    }
    if ($referencefile_path) {

	push(@commands, "--reference_sequence ".$referencefile_path);  #Reference sequence file
    }

    ## Tool specific options
    if ($exclude_nonvariants) {
	
	push(@commands, "--excludeNonVariants");
    }
    if ($genotype_merge_option) {
	
	push(@commands, "--genotypemergeoption ".$genotype_merge_option);
    }
    if ($prioritize_caller) {
	
	push(@commands, "--rod_priority_list ".$prioritize_caller);
    }

    ## Infile
    if (@$infile_paths_ref) {

	push(@commands, "--variant:".join(" --variant:", @$infile_paths_ref));
    }

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
