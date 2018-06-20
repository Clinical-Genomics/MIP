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
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };
use MIP::Program::Base::Gatk qw{ gatk_base };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ gatk_realignertargetcreator gatk_indelrealigner gatk_baserecalibrator gatk_printreads gatk_haplotypecaller };
}

## Constants
Readonly my $SPACE => q{ };

sub gatk_realignertargetcreator {

## Function : Perl wrapper for writing GATK realignertargetcreator recipe to $FILEHANDLE. Based on GATK 3.7.0.
## Returns  : @commands
## Arguments: $known_alleles_ref                     => Input VCF file(s) with known indels {REF}
##          : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##          : $memory_allocation                     => Memory allocation to run Gatk
##          : $java_use_large_pages                  => Use java large pages
##          : $temp_directory                        => Redirect tmp files to java temp
##          : $java_jar                              => Java jar
##          : $infile_path                           => Infile paths
##          : $outfile_path                          => Outfile path
##          : $referencefile_path                    => Reference sequence file
##          : $pedigree                              => Pedigree files
##          : $pedigree_validation_type              => Validation strictness for pedigree
##          : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##          : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##          : $logging_level                         => Set the minimum level of logging
##          : $stderrfile_path                       => Stderrfile path
##          : $FILEHANDLE                            => Sbatch filehandle to write to

    my ($arg_href) = @_;

    ## Default(s)
    my $gatk_disable_auto_index_and_file_lock;
    my $logging_level;

    ## Flatten argument(s)
    my $known_alleles_ref;
    my $intervals_ref;
    my $memory_allocation;
    my $java_use_large_pages;
    my $temp_directory;
    my $java_jar;
    my $infile_path;
    my $outfile_path;
    my $referencefile_path;
    my $pedigree;
    my $pedigree_validation_type;
    my $downsample_to_coverage;
    my $stderrfile_path;
    my $FILEHANDLE;

    my $tmpl = {
        known_alleles_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$known_alleles_ref
        },
        intervals_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$intervals_ref
        },
        memory_allocation => { strict_type => 1, store => \$memory_allocation },
        java_use_large_pages => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$java_use_large_pages
        },
        temp_directory => { strict_type => 1, store => \$temp_directory },
        java_jar       => { strict_type => 1, store => \$java_jar },
        infile_path    => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        outfile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfile_path
        },
        referencefile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$referencefile_path
        },
        pedigree                 => { strict_type => 1, store => \$pedigree },
        pedigree_validation_type => {
            default     => q{SILENT},
            allow       => [qw{ STRICT SILENT}],
            strict_type => 1,
            store       => \$pedigree_validation_type
        },
        downsample_to_coverage => {
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$downsample_to_coverage
        },
        gatk_disable_auto_index_and_file_lock => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$gatk_disable_auto_index_and_file_lock
        },
        logging_level => {
            default     => q{INFO},
            allow       => [qw{ INFO ERROR FATAL }],
            strict_type => 1,
            store       => \$logging_level
        },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        FILEHANDLE => { store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Gatk realignertargetcreator
    # Stores commands depending on input parameters
    my @commands;

    ## Write java core commands to filehandle.
    if ($java_jar) {

        @commands = java_core(
            {
                memory_allocation    => $memory_allocation,
                java_use_large_pages => $java_use_large_pages,
                temp_directory       => $temp_directory,
                java_jar             => $java_jar,
            }
        );
    }

    ## Gatk base args
    @commands = gatk_base(
        {
            commands_ref             => \@commands,
            analysis_type            => q{RealignerTargetCreator},
            logging_level            => $logging_level,
            intervals_ref            => $intervals_ref,
            referencefile_path       => $referencefile_path,
            pedigree                 => $pedigree,
            pedigree_validation_type => $pedigree_validation_type,
            downsample_to_coverage   => $downsample_to_coverage,
            gatk_disable_auto_index_and_file_lock =>
              $gatk_disable_auto_index_and_file_lock,
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    return @commands;
}

sub gatk_indelrealigner {

## Function : Perl wrapper for writing GATK indelrealigner recipe to $FILEHANDLE. Based on GATK 3.7.0.
## Returns  : @commands
## Arguments: $known_alleles_ref                     => Input VCF file(s) with known indels {REF}
##          : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##          : $memory_allocation                     => Memory allocation to run Gatk
##          : $java_use_large_pages                  => Use java large pages
##          : $temp_directory                        => Redirect tmp files to java temp
##          : $java_jar                              => Java jar
##          : $infile_path                           => Infile paths
##          : $outfile_path                          => Outfile path
##          : $target_intervals_file                 => Intervals file output from RealignerTargetCreator
##          : $referencefile_path                    => Reference sequence file
##          : $pedigree                              => Pedigree files
##          : $pedigree_validation_type              => Validation strictness for pedigree
##          : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##          : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##          : $logging_level                         => Set the minimum level of logging
##          : $consensus_determination_model         => Determines how to compute the possible alternate consenses
##          : $stderrfile_path                       => Stderrfile path
##          : $FILEHANDLE                            => Sbatch filehandle to write to

    my ($arg_href) = @_;

    ## Default(s)
    my $gatk_disable_auto_index_and_file_lock;
    my $logging_level;
    my $consensus_determination_model;

    ## Flatten argument(s)
    my $known_alleles_ref;
    my $intervals_ref;
    my $memory_allocation;
    my $java_use_large_pages;
    my $temp_directory;
    my $java_jar;
    my $infile_path;
    my $outfile_path;
    my $target_intervals_file;
    my $referencefile_path;
    my $pedigree;
    my $pedigree_validation_type;
    my $downsample_to_coverage;
    my $stderrfile_path;
    my $FILEHANDLE;

    my $tmpl = {
        known_alleles_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$known_alleles_ref
        },
        intervals_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$intervals_ref
        },
        memory_allocation => { strict_type => 1, store => \$memory_allocation },
        java_use_large_pages => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$java_use_large_pages
        },
        temp_directory => { strict_type => 1, store => \$temp_directory },
        java_jar       => { strict_type => 1, store => \$java_jar },
        infile_path    => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        outfile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfile_path
        },
        target_intervals_file => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$target_intervals_file
        },
        referencefile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$referencefile_path
        },
        pedigree                 => { strict_type => 1, store => \$pedigree },
        pedigree_validation_type => {
            default     => q{SILENT},
            allow       => [qw{ STRICT SILENT}],
            strict_type => 1,
            store       => \$pedigree_validation_type
        },
        downsample_to_coverage => {
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$downsample_to_coverage
        },
        gatk_disable_auto_index_and_file_lock => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$gatk_disable_auto_index_and_file_lock
        },
        logging_level => {
            default     => q{INFO},
            allow       => [qw{ INFO ERROR FATAL }],
            strict_type => 1,
            store       => \$logging_level
        },
        consensus_determination_model => {
            default     => q{USE_READS},
            allow       => [qw{ KNOWNS_ONLY USE_READS USE_SW }],
            strict_type => 1,
            store       => \$consensus_determination_model
        },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        FILEHANDLE => { store => \$FILEHANDLE },

    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Gatk indelrealigne
    # Stores commands depending on input parameters
    my @commands;

    ## Write java core commands to filehandle.
    if ($java_jar) {

        @commands = java_core(
            {
                memory_allocation    => $memory_allocation,
                java_use_large_pages => $java_use_large_pages,
                temp_directory       => $temp_directory,
                java_jar             => $java_jar,
            }
        );
    }

    ### Gatk base args
    @commands = gatk_base(
        {
            commands_ref             => \@commands,
            analysis_type            => q{IndelRealigner},
            logging_level            => $logging_level,
            intervals_ref            => $intervals_ref,
            referencefile_path       => $referencefile_path,
            pedigree                 => $pedigree,
            pedigree_validation_type => $pedigree_validation_type,
            downsample_to_coverage   => $downsample_to_coverage,
            gatk_disable_auto_index_and_file_lock =>
              $gatk_disable_auto_index_and_file_lock,
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    return @commands;
}

sub gatk_baserecalibrator {

## Function : Perl wrapper for writing GATK baserecalibrator recipe to $FILEHANDLE. Based on GATK 3.7.0.
## Returns  : @commands
## Arguments: $known_alleles_ref                     => Input VCF file(s) with known indels {REF}
##          : $covariates_ref                        => One or more covariate to be used in the recalibration. Can be specified multiple times {REF}
##          : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##          : $memory_allocation                     => Memory allocation to run Gatk
##          : $java_use_large_pages                  => Use java large pages
##          : $temp_directory                        => Redirect tmp files to java temp
##          : $java_jar                              => Java jar
##          : $infile_path                           => Infile paths
##          : $outfile_path                          => Outfile path
##          : $referencefile_path                    => Reference sequence file
##          : $pedigree                              => Pedigree files
##          : $pedigree_validation_type              => Validation strictness for pedigree
##          : $num_cpu_threads_per_data_thread       => Number of CPU threads to allocate per data thread
##          : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##          : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##          : $logging_level                         => Set the minimum level of logging
##          : $stderrfile_path                       => Stderrfile path
##          : $FILEHANDLE                            => Sbatch filehandle to write to

    my ($arg_href) = @_;

    ## Default(s)
    my $gatk_disable_auto_index_and_file_lock;
    my $logging_level;

    ## Flatten argument(s)
    my $known_alleles_ref;
    my $covariates_ref;
    my $intervals_ref;
    my $memory_allocation;
    my $java_use_large_pages;
    my $temp_directory;
    my $java_jar;
    my $infile_path;
    my $outfile_path;
    my $referencefile_path;
    my $pedigree;
    my $pedigree_validation_type;
    my $num_cpu_threads_per_data_thread;
    my $downsample_to_coverage;
    my $stderrfile_path;
    my $FILEHANDLE;

    my $tmpl = {
        known_alleles_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$known_alleles_ref
        },
        covariates_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$covariates_ref
        },
        intervals_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$intervals_ref
        },
        memory_allocation => { strict_type => 1, store => \$memory_allocation },
        java_use_large_pages => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$java_use_large_pages
        },
        temp_directory => { strict_type => 1, store => \$temp_directory },
        java_jar       => { strict_type => 1, store => \$java_jar },
        infile_path    => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        outfile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfile_path
        },
        referencefile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$referencefile_path
        },
        pedigree                 => { strict_type => 1, store => \$pedigree },
        pedigree_validation_type => {
            default     => q{SILENT},
            allow       => [qw{ STRICT SILENT}],
            strict_type => 1,
            store       => \$pedigree_validation_type
        },
        num_cpu_threads_per_data_thread => {
            default     => undef,
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$num_cpu_threads_per_data_thread
        },
        downsample_to_coverage => {
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$downsample_to_coverage
        },
        gatk_disable_auto_index_and_file_lock => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$gatk_disable_auto_index_and_file_lock
        },
        logging_level => {
            default     => q{INFO},
            allow       => [qw{ INFO ERROR FATAL }],
            strict_type => 1,
            store       => \$logging_level
        },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        FILEHANDLE => { store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Gatk baserecalibrator
    # Stores commands depending on input parameters
    my @commands;

    ## Write java core commands to filehandle.
    if ($java_jar) {

        @commands = java_core(
            {
                memory_allocation    => $memory_allocation,
                java_use_large_pages => $java_use_large_pages,
                temp_directory       => $temp_directory,
                java_jar             => $java_jar,
            }
        );
    }

    ### Gatk base args
    @commands = gatk_base(
        {
            commands_ref                    => \@commands,
            analysis_type                   => q{BaseRecalibrator},
            logging_level                   => $logging_level,
            intervals_ref                   => $intervals_ref,
            num_cpu_threads_per_data_thread => $num_cpu_threads_per_data_thread,
            referencefile_path              => $referencefile_path,
            pedigree                        => $pedigree,
            pedigree_validation_type        => $pedigree_validation_type,
            downsample_to_coverage          => $downsample_to_coverage,
            gatk_disable_auto_index_and_file_lock =>
              $gatk_disable_auto_index_and_file_lock,
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    return @commands;
}

sub gatk_printreads {

## Function : Perl wrapper for writing GATK printreads recipe to $FILEHANDLE. Based on GATK 3.7.0.
## Returns  : @commands
## Arguments: $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##          : $read_filters_ref                      => Filters to apply to reads before analysis {REF}
##          : $static_quantized_quals_ref            => Use static quantized quality scores to a given number of levels (with -BQSR) {REF}
##          : $memory_allocation                     => Memory allocation to run Gatk
##          : $java_use_large_pages                  => Use java large pages
##          : $temp_directory                        => Redirect tmp files to java temp
##          : $java_jar                              => Java jar
##          : $infile_path                           => Infile paths
##          : $outfile_path                          => Outfile path
##          : $referencefile_path                    => Reference sequence file
##          : $pedigree                              => Pedigree files
##          : $pedigree_validation_type              => Validation strictness for pedigree
##          : $base_quality_score_recalibration_file => Input covariates table file for on-the-fly base quality score recalibration
##          : $num_cpu_threads_per_data_thread       => Number of CPU threads to allocate per data thread
##          : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##          : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##          : $disable_indel_qual                    => Disable printing of base insertion and deletion tags (with -BQSR)
##          : $logging_level                         => Set the minimum level of logging
##          : $stderrfile_path                       => Stderrfile path
##          : $FILEHANDLE                            => Sbatch filehandle to write to

    my ($arg_href) = @_;

    ## Default(s)
    my $gatk_disable_auto_index_and_file_lock;
    my $disable_indel_qual;
    my $logging_level;

    ## Flatten argument(s)
    my $intervals_ref;
    my $read_filters_ref;
    my $static_quantized_quals_ref;
    my $memory_allocation;
    my $java_use_large_pages;
    my $temp_directory;
    my $java_jar;
    my $infile_path;
    my $outfile_path;
    my $referencefile_path;
    my $pedigree;
    my $pedigree_validation_type;
    my $base_quality_score_recalibration_file;
    my $num_cpu_threads_per_data_thread;
    my $downsample_to_coverage;
    my $stderrfile_path;
    my $FILEHANDLE;

    my $tmpl = {
        intervals_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$intervals_ref
        },
        read_filters_ref =>
          { default => [], strict_type => 1, store => \$read_filters_ref },
        static_quantized_quals_ref => {
            default     => [],
            strict_type => 1,
            store       => \$static_quantized_quals_ref
        },
        memory_allocation => { strict_type => 1, store => \$memory_allocation },
        java_use_large_pages => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$java_use_large_pages
        },
        temp_directory => { strict_type => 1, store => \$temp_directory },
        java_jar       => { strict_type => 1, store => \$java_jar },
        infile_path    => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        outfile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfile_path
        },
        referencefile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$referencefile_path
        },
        base_quality_score_recalibration_file => {
            strict_type => 1,
            store       => \$base_quality_score_recalibration_file
        },
        pedigree                 => { strict_type => 1, store => \$pedigree },
        pedigree_validation_type => {
            default     => q{SILENT},
            allow       => [qw{ STRICT SILENT}],
            strict_type => 1,
            store       => \$pedigree_validation_type
        },
        num_cpu_threads_per_data_thread => {
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$num_cpu_threads_per_data_thread
        },
        downsample_to_coverage => {
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$downsample_to_coverage
        },
        gatk_disable_auto_index_and_file_lock => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$gatk_disable_auto_index_and_file_lock
        },
        disable_indel_qual => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$disable_indel_qual
        },
        logging_level => {
            default     => q{INFO},
            allow       => [qw{ INFO ERROR FATAL }],
            strict_type => 1,
            store       => \$logging_level
        },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        FILEHANDLE => { store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Gatk printreads
    # Stores commands depending on input parameters
    my @commands;

    ## Write java core commands to filehandle.
    if ($java_jar) {

        @commands = java_core(
            {
                memory_allocation    => $memory_allocation,
                java_use_large_pages => $java_use_large_pages,
                temp_directory       => $temp_directory,
                java_jar             => $java_jar,
            }
        );
    }

    ### Gatk base args
    @commands = gatk_base(
        {
            commands_ref                    => \@commands,
            analysis_type                   => q{PrintReads},
            logging_level                   => $logging_level,
            intervals_ref                   => $intervals_ref,
            read_filters_ref                => $read_filters_ref,
            static_quantized_quals_ref      => $static_quantized_quals_ref,
            num_cpu_threads_per_data_thread => $num_cpu_threads_per_data_thread,
            referencefile_path              => $referencefile_path,
            pedigree                        => $pedigree,
            pedigree_validation_type        => $pedigree_validation_type,
            downsample_to_coverage          => $downsample_to_coverage,
            gatk_disable_auto_index_and_file_lock =>
              $gatk_disable_auto_index_and_file_lock,
            base_quality_score_recalibration_file =>
              $base_quality_score_recalibration_file,
            disable_indel_qual => $disable_indel_qual,
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    return @commands;
}

sub gatk_haplotypecaller {

## Function : Perl wrapper for writing GATK haplotypecaller recipe to $FILEHANDLE. Based on GATK 3.7.0.
## Returns  : @commands
## Arguments: $intervals_ref                                 => One or more genomic intervals over which to operate {REF}
##          : $read_filters_ref                              => Filters to apply to reads before analysis {REF}
##          : $static_quantized_quals_ref                    => Use static quantized quality scores to a given number of levels (with -BQSR) {REF}
##          : $annotations_ref                               => One or more specific annotations to apply to variant calls
##          : $memory_allocation                             => Memory allocation to run Gatk
##          : $java_use_large_pages                          => Use java large pages
##          : $temp_directory                                => Redirect tmp files to java temp
##          : $java_jar                                      => Java jar
##          : $infile_path                                   => Infile paths
##          : $outfile_path                                  => Outfile path
##          : $referencefile_path                            => Reference sequence file
##          : $base_quality_score_recalibration_file         => Input covariates table file for on-the-fly base quality score recalibration
##          : $pedigree                                      => Pedigree files for samples
##          : $dbsnp                                         => DbSNP file
##          : $standard_min_confidence_threshold_for_calling => The minimum phred-scaled confidence threshold at which variants should be called
##          : $num_cpu_threads_per_data_thread               => Number of CPU threads to allocate per data thread
##          : $downsample_to_coverage                        => Target coverage threshold for downsampling to coverage
##          : $pcr_indel_model                               => The PCR indel model to use
##          : $variant_index_parameter                       => Parameter to pass to the VCF/BCF IndexCreator
##          : $gatk_disable_auto_index_and_file_lock         => Disable both auto-generation of index files and index file locking
##          : $disable_indel_qual                            => Disable printing of base insertion and deletion tags (with -BQSR)
##          : $logging_level                                 => Set the minimum level of logging
##          : $pedigree_validation_type                      => Validation strictness for pedigree
##          : $dont_use_soft_clipped_bases                   => Do not analyze soft clipped bases in the reads
##          : $emit_ref_confidence                           => Mode for emitting reference confidence scores
##          : $variant_index_type                            => Type of IndexCreator to use for VCF/BCF indices
##          : $stderrfile_path                               => Stderrfile path
##          : $FILEHANDLE                                    => Sbatch filehandle to write to
##          : $sample_ploidy                                 => Ploidy per sample

    my ($arg_href) = @_;

    ## Default(s)
    my $gatk_disable_auto_index_and_file_lock;
    my $disable_indel_qual;
    my $logging_level;
    my $pedigree_validation_type;
    my $dont_use_soft_clipped_bases;
    my $emit_ref_confidence;
    my $variant_index_type;

    ## Flatten argument(s)
    my $intervals_ref;
    my $read_filters_ref;
    my $static_quantized_quals_ref;
    my $annotations_ref;
    my $memory_allocation;
    my $java_use_large_pages;
    my $temp_directory;
    my $java_jar;
    my $infile_path;
    my $outfile_path;
    my $referencefile_path;
    my $base_quality_score_recalibration_file;
    my $pedigree;
    my $dbsnp;
    my $standard_min_confidence_threshold_for_calling;
    my $num_cpu_threads_per_data_thread;
    my $downsample_to_coverage;
    my $pcr_indel_model;
    my $variant_index_parameter;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $sample_ploidy;

    my $tmpl = {
        intervals_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$intervals_ref
        },
        read_filters_ref =>
          { default => [], strict_type => 1, store => \$read_filters_ref },
        static_quantized_quals_ref => {
            default     => [],
            strict_type => 1,
            store       => \$static_quantized_quals_ref
        },
        annotations_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$annotations_ref
        },
        memory_allocation => { strict_type => 1, store => \$memory_allocation },
        java_use_large_pages => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$java_use_large_pages
        },
        temp_directory => { strict_type => 1, store => \$temp_directory },
        java_jar       => { strict_type => 1, store => \$java_jar },
        infile_path    => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        outfile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfile_path
        },
        referencefile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$referencefile_path
        },
        base_quality_score_recalibration_file => {
            strict_type => 1,
            store       => \$base_quality_score_recalibration_file
        },
        pedigree => { strict_type => 1, store => \$pedigree },
        dbsnp    => { strict_type => 1, store => \$dbsnp },
        standard_min_confidence_threshold_for_calling => {
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$standard_min_confidence_threshold_for_calling
        },
        num_cpu_threads_per_data_thread => {
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$num_cpu_threads_per_data_thread
        },
        downsample_to_coverage => {
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$downsample_to_coverage
        },
        pcr_indel_model => {
            allow => [ undef, qw{ NONE HOSTILE AGGRESSIVE CONSERVATIVE } ],
            strict_type => 1,
            store       => \$pcr_indel_model
        },
        variant_index_parameter => {
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$variant_index_parameter
        },
        sample_ploidy => {
            allow       => [ undef, qr/ ^\d+$ /sxm ],
            strict_type => 1,
            store       => \$sample_ploidy,
        },
        gatk_disable_auto_index_and_file_lock => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$gatk_disable_auto_index_and_file_lock
        },
        disable_indel_qual => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$disable_indel_qual
        },
        logging_level => {
            default     => q{INFO},
            allow       => [qw{ INFO ERROR FATAL }],
            strict_type => 1,
            store       => \$logging_level
        },
        pedigree_validation_type => {
            default     => q{SILENT},
            allow       => [qw{ STRICT SILENT }],
            strict_type => 1,
            store       => \$pedigree_validation_type
        },
        dont_use_soft_clipped_bases => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$dont_use_soft_clipped_bases
        },
        emit_ref_confidence => {
            default     => q{GVCF},
            allow       => [qw{ NONE BP_RESOLUTION GVCF }],
            strict_type => 1,
            store       => \$emit_ref_confidence
        },
        variant_index_type => {
            default     => q{LINEAR},
            allow       => [qw{ DYNAMIC_SEEK DYNAMIC_SIZE LINEAR INTERVAL }],
            strict_type => 1,
            store       => \$variant_index_type
        },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        FILEHANDLE => { store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Gatk haplotypecaller
    # Stores commands depending on input parameters
    my @commands;

    ## Write java core commands to filehandle.
    if ($java_jar) {

        @commands = java_core(
            {
                memory_allocation    => $memory_allocation,
                java_use_large_pages => $java_use_large_pages,
                temp_directory       => $temp_directory,
                java_jar             => $java_jar,
            }
        );
    }

    ### Gatk base args
    @commands = gatk_base(
        {
            commands_ref                    => \@commands,
            analysis_type                   => q{HaplotypeCaller},
            logging_level                   => $logging_level,
            intervals_ref                   => $intervals_ref,
            read_filters_ref                => $read_filters_ref,
            static_quantized_quals_ref      => $static_quantized_quals_ref,
            referencefile_path              => $referencefile_path,
            num_cpu_threads_per_data_thread => $num_cpu_threads_per_data_thread,
            pedigree                        => $pedigree,
            pedigree_validation_type        => $pedigree_validation_type,
            downsample_to_coverage          => $downsample_to_coverage,
            gatk_disable_auto_index_and_file_lock =>
              $gatk_disable_auto_index_and_file_lock,
            base_quality_score_recalibration_file =>
              $base_quality_score_recalibration_file,
            disable_indel_qual => $disable_indel_qual,
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
    if ($sample_ploidy) {

        push @commands, q{--sample_ploidy} . $SPACE . $sample_ploidy;
    }

    push @commands, q{--emitRefConfidence} . $SPACE . $emit_ref_confidence;

    push @commands, q{--variant_index_type} . $SPACE . $variant_index_type;

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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    return @commands;
}

1;
