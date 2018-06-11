package MIP::Program::Alignment::Gatk;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };

## CPANM
use Readonly;

## MIPs lib/
use MIP::Language::Java qw{ java_core };
use MIP::Program::Base::Gatk qw{ gatk_base };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ gatk_baserecalibrator gatk_haplotypecaller gatk_indelrealigner gatk_printreads gatk_realignertargetcreator gatk_splitncigarreads };
}

## Constants
Readonly my $SPACE => q{ };

sub gatk_baserecalibrator {

## Function : Perl wrapper for writing GATK baserecalibrator recipe to $FILEHANDLE. Based on GATK 3.7.0.
## Returns  : @commands
## Arguments: $covariates_ref                        => One or more covariate to be used in the recalibration. Can be specified multiple times {REF}
##          : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##          : $FILEHANDLE                            => Sbatch filehandle to write to
##          : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##          : $infile_path                           => Infile paths
##          : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##          : $java_jar                              => Java jar
##          : $java_use_large_pages                  => Use java large pages
##          : $known_alleles_ref                     => Input VCF file(s) with known indels {REF}
##          : $logging_level                         => Set the minimum level of logging
##          : $memory_allocation                     => Memory allocation to run Gatk
##          : $num_cpu_threads_per_data_thread       => Number of CPU threads to allocate per data thread
##          : $outfile_path                          => Outfile path
##          : $pedigree                              => Pedigree files
##          : $pedigree_validation_type              => Validation strictness for pedigree
##          : $referencefile_path                    => Reference sequence file
##          : $stderrfile_path                       => Stderrfile path
##          : $temp_directory                        => Redirect tmp files to java temp

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $covariates_ref;
    my $downsample_to_coverage;
    my $FILEHANDLE;
    my $infile_path;
    my $intervals_ref;
    my $java_jar;
    my $java_use_large_pages;
    my $known_alleles_ref;
    my $memory_allocation;
    my $num_cpu_threads_per_data_thread;
    my $outfile_path;
    my $pedigree;
    my $pedigree_validation_type;
    my $referencefile_path;
    my $stderrfile_path;
    my $temp_directory;

    ## Default(s)
    my $gatk_disable_auto_index_and_file_lock;
    my $logging_level;

    my $tmpl = {
        covariates_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$covariates_ref,
            strict_type => 1,
        },
        downsample_to_coverage => {
            allow       => qr/ ^\d+$ /sxm,
            store       => \$downsample_to_coverage,
            strict_type => 1,
        },

        FILEHANDLE                            => { store => \$FILEHANDLE },
        gatk_disable_auto_index_and_file_lock => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$gatk_disable_auto_index_and_file_lock,
            strict_type => 1,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        intervals_ref => {
            default     => [],
            defined     => 1,
            store       => \$intervals_ref,
            strict_type => 1,
        },
        java_jar             => { store => \$java_jar, strict_type => 1, },
        java_use_large_pages => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$java_use_large_pages,
            strict_type => 1,
        },
        known_alleles_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$known_alleles_ref,
            strict_type => 1,
        },
        logging_level => {
            allow       => [qw{ INFO ERROR FATAL }],
            default     => q{INFO},
            store       => \$logging_level,
            strict_type => 1,
        },
        memory_allocation =>
          { store => \$memory_allocation, strict_type => 1, },
        num_cpu_threads_per_data_thread => {
            allow       => qr/ ^\d+$ /sxm,
            default     => undef,
            store       => \$num_cpu_threads_per_data_thread,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        pedigree                 => { store => \$pedigree, strict_type => 1, },
        pedigree_validation_type => {
            allow       => [qw{ STRICT SILENT}],
            default     => q{SILENT},
            store       => \$pedigree_validation_type,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
            required    => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        temp_directory  => { store => \$temp_directory,  strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Gatk baserecalibrator
    # Stores commands depending on input parameters
    my @commands;

    ## Write java core commands to filehandle.
    if ($java_jar) {

        @commands = java_core(
            {
                java_jar             => $java_jar,
                java_use_large_pages => $java_use_large_pages,
                memory_allocation    => $memory_allocation,
                temp_directory       => $temp_directory,
            }
        );
    }

    ### Gatk base args
    @commands = gatk_base(
        {
            analysis_type          => q{BaseRecalibrator},
            commands_ref           => \@commands,
            downsample_to_coverage => $downsample_to_coverage,
            gatk_disable_auto_index_and_file_lock =>
              $gatk_disable_auto_index_and_file_lock,
            intervals_ref                   => $intervals_ref,
            logging_level                   => $logging_level,
            num_cpu_threads_per_data_thread => $num_cpu_threads_per_data_thread,
            pedigree                        => $pedigree,
            pedigree_validation_type        => $pedigree_validation_type,
            referencefile_path              => $referencefile_path,
        }
    );

    ## Known alleles reference
    push @commands,
      q{--knownSites} . $SPACE . join $SPACE . q{--knownSites} . $SPACE,
      @{$known_alleles_ref};

    ## Covariates
    push @commands,
      q{--covariate} . $SPACE . join $SPACE . q{--covariate} . $SPACE,
      @{$covariates_ref};

    ## Infile
    push @commands, q{--input_file} . $SPACE . $infile_path;

    ## Output
    push @commands, q{--out} . $SPACE . $outfile_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gatk_haplotypecaller {

## Function : Perl wrapper for writing GATK haplotypecaller recipe to $FILEHANDLE. Based on GATK 3.7.0.
## Returns  : @commands
## Arguments: $annotations_ref                               => One or more specific annotations to apply to variant calls
##          : $base_quality_score_recalibration_file         => Input covariates table file for on-the-fly base quality score recalibration
##          : $dbsnp                                         => DbSNP file
##          : $disable_indel_qual                            => Disable printing of base insertion and deletion tags (with -BQSR)
##          : $dont_use_soft_clipped_bases                   => Do not analyze soft clipped bases in the reads
##          : $downsample_to_coverage                        => Target coverage threshold for downsampling to coverage
##          : $emit_ref_confidence                           => Mode for emitting reference confidence scores
##          : $FILEHANDLE                                    => Sbatch filehandle to write to
##          : $gatk_disable_auto_index_and_file_lock         => Disable both auto-generation of index files and index file locking
##          : $infile_path                                   => Infile paths
##          : $intervals_ref                                 => One or more genomic intervals over which to operate {REF}
##          : $java_jar                                      => Java jar
##          : $java_use_large_pages                          => Use java large pages
##          : $logging_level                                 => Set the minimum level of logging
##          : $memory_allocation                             => Memory allocation to run Gatk
##          : $num_cpu_threads_per_data_thread               => Number of CPU threads to allocate per data thread
##          : $outfile_path                                  => Outfile path
##          : $pcr_indel_model                               => The PCR indel model to use
##          : $pedigree                                      => Pedigree files for samples
##          : $pedigree_validation_type                      => Validation strictness for pedigree
##          : $read_filters_ref                              => Filters to apply to reads before analysis {REF}
##          : $referencefile_path                            => Reference sequence file
##          : $standard_min_confidence_threshold_for_calling => The minimum phred-scaled confidence threshold at which variants should be called
##          : $static_quantized_quals_ref                    => Use static quantized quality scores to a given number of levels (with -BQSR) {REF}
##          : $stderrfile_path                               => Stderrfile path
##          : $temp_directory                                => Redirect tmp files to java temp
##          : $variant_index_parameter                       => Parameter to pass to the VCF/BCF IndexCreator
##          : $variant_index_type                            => Type of IndexCreator to use for VCF/BCF indices

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $annotations_ref;
    my $base_quality_score_recalibration_file;
    my $dbsnp;
    my $downsample_to_coverage;
    my $FILEHANDLE;
    my $infile_path;
    my $intervals_ref;
    my $java_jar;
    my $java_use_large_pages;
    my $memory_allocation;
    my $num_cpu_threads_per_data_thread;
    my $outfile_path;
    my $pcr_indel_model;
    my $pedigree;
    my $read_filters_ref;
    my $referencefile_path;
    my $standard_min_confidence_threshold_for_calling;
    my $static_quantized_quals_ref;
    my $stderrfile_path;
    my $temp_directory;
    my $variant_index_parameter;

    ## Default(s)
    my $disable_indel_qual;
    my $dont_use_soft_clipped_bases;
    my $emit_ref_confidence;
    my $gatk_disable_auto_index_and_file_lock;
    my $logging_level;
    my $pedigree_validation_type;
    my $variant_index_type;

    my $tmpl = {
        annotations_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$annotations_ref,
            strict_type => 1,
        },
        base_quality_score_recalibration_file => {
            store       => \$base_quality_score_recalibration_file,
            strict_type => 1,
        },
        dbsnp              => { store => \$dbsnp, strict_type => 1, },
        disable_indel_qual => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$disable_indel_qual,
            strict_type => 1,
        },
        dont_use_soft_clipped_bases => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$dont_use_soft_clipped_bases,
            strict_type => 1,
        },
        downsample_to_coverage => {
            allow       => qr/ ^\d+$ /sxm,
            store       => \$downsample_to_coverage,
            strict_type => 1,
        },
        emit_ref_confidence => {
            allow       => [qw{ NONE BP_RESOLUTION GVCF }],
            default     => q{GVCF},
            store       => \$emit_ref_confidence,
            strict_type => 1,
        },
        FILEHANDLE                            => { store => \$FILEHANDLE },
        gatk_disable_auto_index_and_file_lock => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$gatk_disable_auto_index_and_file_lock,
            strict_type => 1,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        intervals_ref => {
            default     => [],
            defined     => 1,
            store       => \$intervals_ref,
            strict_type => 1,
        },
        java_jar             => { store => \$java_jar, strict_type => 1, },
        java_use_large_pages => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$java_use_large_pages,
            strict_type => 1,
        },
        logging_level => {
            allow       => [qw{ INFO ERROR FATAL }],
            default     => q{INFO},
            store       => \$logging_level,
            strict_type => 1,
        },
        memory_allocation =>
          { store => \$memory_allocation, strict_type => 1, },
        num_cpu_threads_per_data_thread => {
            allow       => qr/ ^\d+$ /sxm,
            store       => \$num_cpu_threads_per_data_thread,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        pedigree                 => { store => \$pedigree, strict_type => 1, },
        pedigree_validation_type => {
            allow       => [qw{ STRICT SILENT }],
            default     => q{SILENT},
            store       => \$pedigree_validation_type,
            strict_type => 1,
        },
        pcr_indel_model => {
            allow => [ undef, qw{ NONE HOSTILE AGGRESSIVE CONSERVATIVE } ],
            store       => \$pcr_indel_model,
            strict_type => 1,
        },
        read_filters_ref => {
            default     => [],
            store       => \$read_filters_ref,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
            required    => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        standard_min_confidence_threshold_for_calling => {
            allow       => qr/ ^\d+$ /sxm,
            store       => \$standard_min_confidence_threshold_for_calling,
            strict_type => 1,
        },
        static_quantized_quals_ref => {
            default     => [],
            store       => \$static_quantized_quals_ref,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        temp_directory  => { store => \$temp_directory,  strict_type => 1, },
        variant_index_parameter => {
            allow       => qr/ ^\d+$ /sxm,
            store       => \$variant_index_parameter,
            strict_type => 1,
        },
        variant_index_type => {
            allow       => [qw{ 0 DYNAMIC_SEEK DYNAMIC_SIZE LINEAR INTERVAL }],
            default     => q{LINEAR},
            store       => \$variant_index_type,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Gatk haplotypecaller
    # Stores commands depending on input parameters
    my @commands;

    ## Write java core commands to filehandle.
    if ($java_jar) {

        @commands = java_core(
            {
                java_jar             => $java_jar,
                java_use_large_pages => $java_use_large_pages,
                memory_allocation    => $memory_allocation,
                temp_directory       => $temp_directory,
            }
        );
    }

    ### Gatk base args
    @commands = gatk_base(
        {
            analysis_type => q{HaplotypeCaller},
            base_quality_score_recalibration_file =>
              $base_quality_score_recalibration_file,
            commands_ref           => \@commands,
            disable_indel_qual     => $disable_indel_qual,
            downsample_to_coverage => $downsample_to_coverage,
            gatk_disable_auto_index_and_file_lock =>
              $gatk_disable_auto_index_and_file_lock,
            intervals_ref                   => $intervals_ref,
            logging_level                   => $logging_level,
            num_cpu_threads_per_data_thread => $num_cpu_threads_per_data_thread,
            pedigree                        => $pedigree,
            pedigree_validation_type        => $pedigree_validation_type,
            read_filters_ref                => $read_filters_ref,
            referencefile_path              => $referencefile_path,
            static_quantized_quals_ref      => $static_quantized_quals_ref,
        }
    );

    if ($dbsnp) {

        push @commands, q{--dbsnp} . $SPACE . $dbsnp;
    }
    if ($standard_min_confidence_threshold_for_calling) {

        push @commands,
            q{--standard_min_confidence_threshold_for_calling}
          . $SPACE
          . $standard_min_confidence_threshold_for_calling;
    }
    if ($dont_use_soft_clipped_bases) {

        push @commands, q{--dontUseSoftClippedBases};
    }
    if ($pcr_indel_model) {

        push @commands, q{--pcr_indel_model} . $SPACE . $pcr_indel_model;
    }
    if ( @{$annotations_ref} ) {

        push @commands,
          q{--annotation} . $SPACE . join $SPACE . q{--annotation} . $SPACE,
          @{$annotations_ref};
    }
    if ($variant_index_parameter) {

        push @commands,
          q{--variant_index_parameter} . $SPACE . $variant_index_parameter;
    }

    if ($variant_index_type) {
        push @commands, q{--variant_index_type} . $SPACE . $variant_index_type;
    }

    push @commands, q{--emitRefConfidence} . $SPACE . $emit_ref_confidence;

    ## Infile
    push @commands, q{--input_file} . $SPACE . $infile_path;

    ## Output
    push @commands, q{--out} . $SPACE . $outfile_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gatk_indelrealigner {

## Function : Perl wrapper for writing GATK indelrealigner recipe to $FILEHANDLE. Based on GATK 3.7.0.
## Returns  : @commands
## Arguments: $consensus_determination_model         => Determines how to compute the possible alternate consenses
##          : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##          : $FILEHANDLE                            => Sbatch filehandle to write to
##          : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##          : $infile_path                           => Infile paths
##          : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##          : $java_jar                              => Java jar
##          : $java_use_large_pages                  => Use java large pages
##          : $known_alleles_ref                     => Input VCF file(s) with known indels {REF}
##          : $logging_level                         => Set the minimum level of logging
##          : $memory_allocation                     => Memory allocation to run Gatk
##          : $outfile_path                          => Outfile path
##          : $pedigree                              => Pedigree files
##          : $pedigree_validation_type              => Validation strictness for pedigree
##          : $referencefile_path                    => Reference sequence file
##          : $stderrfile_path                       => Stderrfile path
##          : $target_intervals_file                 => Intervals file output from RealignerTargetCreator
##          : $temp_directory                        => Redirect tmp files to java temp

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $downsample_to_coverage;
    my $FILEHANDLE;
    my $infile_path;
    my $intervals_ref;
    my $java_jar;
    my $java_use_large_pages;
    my $known_alleles_ref;
    my $memory_allocation;
    my $outfile_path;
    my $pedigree;
    my $pedigree_validation_type;
    my $referencefile_path;
    my $stderrfile_path;
    my $target_intervals_file;
    my $temp_directory;

    ## Default(s)
    my $consensus_determination_model;
    my $gatk_disable_auto_index_and_file_lock;
    my $logging_level;

    my $tmpl = {
        consensus_determination_model => {
            allow       => [qw{ KNOWNS_ONLY USE_READS USE_SW }],
            default     => q{USE_READS},
            store       => \$consensus_determination_model,
            strict_type => 1,
        },
        downsample_to_coverage => {
            allow       => qr/ ^\d+$ /sxm,
            store       => \$downsample_to_coverage,
            strict_type => 1,
        },
        FILEHANDLE                            => { store => \$FILEHANDLE },
        gatk_disable_auto_index_and_file_lock => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$gatk_disable_auto_index_and_file_lock,
            strict_type => 1,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        intervals_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$intervals_ref,
            strict_type => 1,
        },
        java_jar             => { store => \$java_jar, strict_type => 1, },
        java_use_large_pages => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$java_use_large_pages,
            strict_type => 1,
        },
        known_alleles_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$known_alleles_ref,
            strict_type => 1,
        },
        logging_level => {
            default     => q{INFO},
            allow       => [qw{ INFO ERROR FATAL }],
            strict_type => 1,
            store       => \$logging_level
        },
        memory_allocation => { strict_type => 1, store => \$memory_allocation },
        outfile_path      => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        pedigree                 => { store => \$pedigree, strict_type => 1, },
        pedigree_validation_type => {
            allow       => [qw{ STRICT SILENT}],
            default     => q{SILENT},
            store       => \$pedigree_validation_type,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
            required    => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        target_intervals_file => {
            defined     => 1,
            required    => 1,
            store       => \$target_intervals_file,
            strict_type => 1,
        },
        temp_directory => { store => \$temp_directory, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Gatk indelrealigne
    # Stores commands depending on input parameters
    my @commands;

    ## Write java core commands to filehandle.
    if ($java_jar) {

        @commands = java_core(
            {
                java_jar             => $java_jar,
                java_use_large_pages => $java_use_large_pages,
                memory_allocation    => $memory_allocation,
                temp_directory       => $temp_directory,
            }
        );
    }

    ### Gatk base args
    @commands = gatk_base(
        {
            analysis_type          => q{IndelRealigner},
            commands_ref           => \@commands,
            downsample_to_coverage => $downsample_to_coverage,
            gatk_disable_auto_index_and_file_lock =>
              $gatk_disable_auto_index_and_file_lock,
            intervals_ref            => $intervals_ref,
            logging_level            => $logging_level,
            pedigree                 => $pedigree,
            pedigree_validation_type => $pedigree_validation_type,
            referencefile_path       => $referencefile_path,
        }
    );

    if ($consensus_determination_model) {

        push @commands,
            q{--consensusDeterminationModel}
          . $SPACE
          . $consensus_determination_model;
    }

    push @commands, q{--targetIntervals} . $SPACE . $target_intervals_file;

    ## Known alleles reference
    push @commands,
      q{--knownAlleles} . $SPACE . join $SPACE . q{--knownAlleles} . $SPACE,
      @{$known_alleles_ref};

    ## Infile
    push @commands, q{--input_file} . $SPACE . $infile_path;

    ## Output
    push @commands, q{--out} . $SPACE . $outfile_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gatk_printreads {

## Function : Perl wrapper for writing GATK printreads recipe to $FILEHANDLE. Based on GATK 3.7.0.
## Returns  : @commands
## Arguments: $base_quality_score_recalibration_file => Input covariates table file for on-the-fly base quality score recalibration
##          : $disable_indel_qual                    => Disable printing of base insertion and deletion tags (with -BQSR)
##          : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##          : $FILEHANDLE                            => Sbatch filehandle to write to
##          : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##          : $infile_path                           => Infile paths
##          : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##          : $java_jar                              => Java jar
##          : $java_use_large_pages                  => Use java large pages
##          : $logging_level                         => Set the minimum level of logging
##          : $memory_allocation                     => Memory allocation to run Gatk
##          : $num_cpu_threads_per_data_thread       => Number of CPU threads to allocate per data thread
##          : $outfile_path                          => Outfile path
##          : $pedigree                              => Pedigree files
##          : $pedigree_validation_type              => Validation strictness for pedigree
##          : $read_filters_ref                      => Filters to apply to reads before analysis {REF}
##          : $referencefile_path                    => Reference sequence file
##          : $static_quantized_quals_ref            => Use static quantized quality scores to a given number of levels (with -BQSR) {REF}
##          : $stderrfile_path                       => Stderrfile path
##          : $temp_directory                        => Redirect tmp files to java temp

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $base_quality_score_recalibration_file;
    my $downsample_to_coverage;
    my $FILEHANDLE;
    my $infile_path;
    my $intervals_ref;
    my $java_jar;
    my $java_use_large_pages;
    my $memory_allocation;
    my $num_cpu_threads_per_data_thread;
    my $outfile_path;
    my $pedigree;
    my $pedigree_validation_type;
    my $read_filters_ref;
    my $referencefile_path;
    my $static_quantized_quals_ref;
    my $stderrfile_path;
    my $temp_directory;

    ## Default(s)
    my $disable_indel_qual;
    my $gatk_disable_auto_index_and_file_lock;
    my $logging_level;

    my $tmpl = {
        base_quality_score_recalibration_file => {
            store       => \$base_quality_score_recalibration_file,
            strict_type => 1,
        },
        disable_indel_qual => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$disable_indel_qual,
            strict_type => 1,
        },
        downsample_to_coverage => {
            allow       => qr/ ^\d+$ /sxm,
            store       => \$downsample_to_coverage,
            strict_type => 1,
        },
        FILEHANDLE                            => { store => \$FILEHANDLE },
        gatk_disable_auto_index_and_file_lock => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$gatk_disable_auto_index_and_file_lock,
            strict_type => 1,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        intervals_ref => {
            default     => [],
            defined     => 1,
            store       => \$intervals_ref,
            strict_type => 1,
        },
        java_jar             => { store => \$java_jar, strict_type => 1, },
        java_use_large_pages => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$java_use_large_pages,
            strict_type => 1,
        },
        logging_level => {
            allow       => [qw{ INFO ERROR FATAL }],
            default     => q{INFO},
            store       => \$logging_level,
            strict_type => 1,
        },
        memory_allocation =>
          { store => \$memory_allocation, strict_type => 1, },
        num_cpu_threads_per_data_thread => {
            allow       => qr/ ^\d+$ /sxm,
            store       => \$num_cpu_threads_per_data_thread,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        pedigree                 => { store => \$pedigree, strict_type => 1, },
        pedigree_validation_type => {
            allow       => [qw{ STRICT SILENT}],
            default     => q{SILENT},
            store       => \$pedigree_validation_type,
            strict_type => 1,
        },
        read_filters_ref => {
            default     => [],
            store       => \$read_filters_ref,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
            required    => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        static_quantized_quals_ref => {
            default     => [],
            store       => \$static_quantized_quals_ref,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        temp_directory  => { store => \$temp_directory,  strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Gatk printreads
    # Stores commands depending on input parameters
    my @commands;

    ## Write java core commands to filehandle.
    if ($java_jar) {

        @commands = java_core(
            {
                java_jar             => $java_jar,
                java_use_large_pages => $java_use_large_pages,
                memory_allocation    => $memory_allocation,
                temp_directory       => $temp_directory,
            }
        );
    }

    ### Gatk base args
    @commands = gatk_base(
        {
            analysis_type => q{PrintReads},
            base_quality_score_recalibration_file =>
              $base_quality_score_recalibration_file,
            commands_ref           => \@commands,
            disable_indel_qual     => $disable_indel_qual,
            downsample_to_coverage => $downsample_to_coverage,
            gatk_disable_auto_index_and_file_lock =>
              $gatk_disable_auto_index_and_file_lock,
            intervals_ref                   => $intervals_ref,
            logging_level                   => $logging_level,
            num_cpu_threads_per_data_thread => $num_cpu_threads_per_data_thread,
            pedigree                        => $pedigree,
            pedigree_validation_type        => $pedigree_validation_type,
            read_filters_ref                => $read_filters_ref,
            referencefile_path              => $referencefile_path,
            static_quantized_quals_ref      => $static_quantized_quals_ref,
        }
    );

    ## Infile
    push @commands, q{--input_file} . $SPACE . $infile_path;

    ## Output
    push @commands, q{--out} . $SPACE . $outfile_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gatk_realignertargetcreator {

## Function : Perl wrapper for writing GATK realignertargetcreator recipe to $FILEHANDLE. Based on GATK 3.7.0.
## Returns  : @commands
## Arguments: $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##          : $FILEHANDLE                            => Sbatch filehandle to write to
##          : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##          : $infile_path                           => Infile paths
##          : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##          : $java_jar                              => Java jar
##          : $java_use_large_pages                  => Use java large pages
##          : $known_alleles_ref                     => Input VCF file(s) with known indels {REF}
##          : $logging_level                         => Set the minimum level of logging
##          : $memory_allocation                     => Memory allocation to run Gatk
##          : $outfile_path                          => Outfile path
##          : $pedigree                              => Pedigree files
##          : $pedigree_validation_type              => Validation strictness for pedigree
##          : $referencefile_path                    => Reference sequence file
##          : $stderrfile_path                       => Stderrfile path
##          : $temp_directory                        => Redirect tmp files to java temp

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $downsample_to_coverage;
    my $FILEHANDLE;
    my $infile_path;
    my $intervals_ref;
    my $java_jar;
    my $java_use_large_pages;
    my $known_alleles_ref;
    my $memory_allocation;
    my $outfile_path;
    my $pedigree;
    my $pedigree_validation_type;
    my $referencefile_path;
    my $stderrfile_path;
    my $temp_directory;

    ## Default(s)
    my $gatk_disable_auto_index_and_file_lock;
    my $logging_level;

    my $tmpl = {
        downsample_to_coverage => {
            allow       => qr/ ^\d+$ /sxm,
            store       => \$downsample_to_coverage,
            strict_type => 1,
        },
        FILEHANDLE                            => { store => \$FILEHANDLE },
        gatk_disable_auto_index_and_file_lock => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$gatk_disable_auto_index_and_file_lock,
            strict_type => 1,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        intervals_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$intervals_ref,
            strict_type => 1,
        },
        java_jar             => { store => \$java_jar, strict_type => 1, },
        java_use_large_pages => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$java_use_large_pages,
            strict_type => 1,
        },
        known_alleles_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$known_alleles_ref,
            strict_type => 1,
        },
        logging_level => {
            allow       => [qw{ INFO ERROR FATAL }],
            default     => q{INFO},
            store       => \$logging_level,
            strict_type => 1,
        },
        memory_allocation =>
          { store => \$memory_allocation, strict_type => 1, },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        pedigree                 => { store => \$pedigree, strict_type => 1, },
        pedigree_validation_type => {
            allow       => [qw{ STRICT SILENT}],
            default     => q{SILENT},
            store       => \$pedigree_validation_type,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
            required    => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        temp_directory  => { store => \$temp_directory,  strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Gatk realignertargetcreator
    # Stores commands depending on input parameters
    my @commands;

    ## Write java core commands to filehandle.
    if ($java_jar) {

        @commands = java_core(
            {
                java_jar             => $java_jar,
                java_use_large_pages => $java_use_large_pages,
                memory_allocation    => $memory_allocation,
                temp_directory       => $temp_directory,
            }
        );
    }

    ## Gatk base args
    @commands = gatk_base(
        {
            analysis_type          => q{RealignerTargetCreator},
            commands_ref           => \@commands,
            downsample_to_coverage => $downsample_to_coverage,
            gatk_disable_auto_index_and_file_lock =>
              $gatk_disable_auto_index_and_file_lock,
            intervals_ref            => $intervals_ref,
            logging_level            => $logging_level,
            pedigree                 => $pedigree,
            pedigree_validation_type => $pedigree_validation_type,
            referencefile_path       => $referencefile_path,
        }
    );

    ## Known alleles reference
    push @commands, q{--known} . $SPACE . join $SPACE . q{--known} . $SPACE,
      @{$known_alleles_ref};

    ## Infile
    push @commands, q{--input_file} . $SPACE . $infile_path;

    ## Output
    push @commands, q{--out} . $SPACE . $outfile_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gatk_splitncigarreads {

## Function : Perl wrapper for writing GATK splitNCigarReads recipe to $FILEHANDLE. Based on GATK 3.7.0.
## Returns  : @commands
## Arguments: $FILEHANDLE           => Sbatch filehandle to write to
##          : $infile_path          => Infile paths
##          : $java_jar             => Java jar
##          : $java_use_large_pages => Use java large pages
##          : $logging_level        => Set the minimum level of logging
##          : $maping_quality_from  => Input mapping quality values
##          : $maping_quality_to    => Output mapping quality values
##          : $memory_allocation    => Memory allocation to run Gatk
##          : $outfile_path         => Outfile path
##          : $operation            => The read operation that is to be applied on the bam file
##          : $readfilter           => Read-filter settings
##          : $referencefile_path   => Reference sequence file
##          : $stderrfile_path      => Stderrfile path
##          : $temp_directory       => Redirect tmp files to java temp

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $java_jar;
    my $java_use_large_pages;
    my $memory_allocation;
    my $outfile_path;
    my $stderrfile_path;
    my $referencefile_path;
    my $temp_directory;

    ## Default(s)
    my $logging_level;
    my $maping_quality_from;
    my $maping_quality_to;
    my $operation;
    my $readfilter;

    ## Constants
    Readonly my $MAPPING_QUALITY_FROM => 255;
    Readonly my $MAPPING_QUALITY_TO   => 60;

    my $tmpl = {
        FILEHANDLE  => { store => \$FILEHANDLE, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        java_jar             => { store => \$java_jar, strict_type => 1, },
        java_use_large_pages => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$java_use_large_pages,
            strict_type => 1,
        },
        logging_level => {
            allow       => [qw{ INFO ERROR FATAL }],
            default     => q{INFO},
            store       => \$logging_level,
            strict_type => 1,
        },
        maping_quality_from => {
            default     => $MAPPING_QUALITY_FROM,
            defined     => 1,
            store       => \$maping_quality_from,
            strict_type => 1,
        },
        maping_quality_to => {
            default     => $MAPPING_QUALITY_TO,
            defined     => 1,
            store       => \$maping_quality_to,
            strict_type => 1,
        },
        memory_allocation =>
          { store => \$memory_allocation, strict_type => 1, },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        operation => {
            default     => q{ALLOW_N_CIGAR_READS},
            defined     => 1,
            store       => \$operation,
            strict_type => 1,
        },
        readfilter => {
            default     => q{ReassignOneMappingQuality},
            defined     => 1,
            store       => \$readfilter,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
            required    => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        temp_directory  => { store => \$temp_directory,  strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands;

    ## Write java core commands to filehandle.
    if ($java_jar) {

        @commands = java_core(
            {
                java_jar             => $java_jar,
                java_use_large_pages => $java_use_large_pages,
                memory_allocation    => $memory_allocation,
                temp_directory       => $temp_directory,

            }
        );
    }

    ### Gatk base args
    @commands = gatk_base(
        {
            analysis_type      => q{SplitNCigarReads},
            commands_ref       => \@commands,
            logging_level      => $logging_level,
            referencefile_path => $referencefile_path,
        }
    );

    push @commands, q{-I} . $SPACE . $infile_path;

    push @commands, q{-o} . $SPACE . $outfile_path;

    push @commands, q{-rf} . $SPACE . $readfilter;

    push @commands, q{-RMQF} . $SPACE . $maping_quality_from;

    push @commands, q{-RMQT} . $SPACE . $maping_quality_to;

    push @commands, q{-U} . $SPACE . $operation;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );
    return @commands;
}

1;
