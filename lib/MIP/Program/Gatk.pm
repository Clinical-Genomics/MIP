package MIP::Program::Gatk;

use 5.026;
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
use MIP::Constants
  qw{ $ASTERISK $AMPERSAND $COLON $DOT $DOUBLE_QUOTE $EMPTY_STR $ESCAPE $NEWLINE $SPACE $UNDERSCORE };
use MIP::Language::Java qw{ java_core };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.20;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      gatk_applybqsr
      gatk_applyvqsr
      gatk_asereadcounter
      gatk_base
      gatk_baserecalibrator
      gatk_calculategenotypeposteriors
      gatk_cnnscorevariants
      gatk_combinevariants
      gatk_common_options
      gatk_concatenate_variants
      gatk_gatherbqsrreports
      gatk_gathervcfscloud
      gatk_genomicsdbimport
      gatk_genotypegvcfs
      gatk_haplotypecaller
      gatk_indexfeaturefile
      gatk_java_options
      gatk_leftalignandtrimvariants
      gatk_printreads
      gatk_selectvariants
      gatk_splitncigarreads
      gatk_varianteval
      gatk_variantfiltration
      gatk_variantrecalibrator
    };
}

sub gatk_applybqsr {

## Function : Perl wrapper for writing GATK ApplyBQSR recipe to $filehandle. Based on GATK 4.0.8
## Returns  : @commands
## Arguments: $base_quality_score_recalibration_file => Input recalibration table for BQSR
##          : $filehandle                            => Sbatch filehandle to write to
##          : $infile_path                           => Infile paths
##          : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##          : $java_use_large_pages                  => Use java large pages
##          : $memory_allocation                     => Memory allocation to run Gatk
##          : $outfile_path                          => Outfile path
##          : $read_filters_ref                      => Filters to apply on reads {REF}
##          : $referencefile_path                    => Reference sequence file
##          : $static_quantized_quals_ref            => Use static quantized quality scores to a given number of levels {REF}
##          : $stderrfile_path                       => Stderrfile path
##          : $temp_directory                        => Redirect tmp files to java temp
##          : $verbosity                             => Set the minimum level of logging
##          : $xargs_mode                            => Set if the program will be executed via xargs

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $base_quality_score_recalibration_file;
    my $filehandle;
    my $infile_path;
    my $intervals_ref;
    my $memory_allocation;
    my $outfile_path;
    my $read_filters_ref;
    my $referencefile_path;
    my $static_quantized_quals_ref;
    my $stderrfile_path;
    my $temp_directory;

    ## Default(s)
    my $java_use_large_pages;
    my $verbosity;
    my $xargs_mode;

    my $tmpl = {
        base_quality_score_recalibration_file => {
            defined     => 1,
            required    => 1,
            store       => \$base_quality_score_recalibration_file,
            strict_type => 1,
        },
        filehandle  => { store => \$filehandle, },
        infile_path => {
            allow       => qr/ (?: bam | sam | cram )$ /xms,
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
        java_use_large_pages => {
            allow       => [ undef, 0, 1 ],
            store       => \$java_use_large_pages,
            strict_type => 1,
        },
        memory_allocation => {
            store       => \$memory_allocation,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        read_filters_ref => {
            default     => [],
            store       => \$read_filters_ref,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
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
        verbosity       => {
            allow       => [qw{ INFO ERROR FATAL }],
            default     => q{INFO},
            store       => \$verbosity,
            strict_type => 1,
        },
        xargs_mode => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$xargs_mode,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ gatk };

    ## Add java options
    gatk_java_options(
        {
            commands_ref         => \@commands,
            java_use_large_pages => $java_use_large_pages,
            memory_allocation    => $memory_allocation,
            xargs_mode           => $xargs_mode,
        }
    );

    ## Add tool command
    push @commands, q{ApplyBQSR};

    ## Add infile
    push @commands, q{--input} . $SPACE . $infile_path;

    ## Add common options
    gatk_common_options(
        {
            commands_ref       => \@commands,
            intervals_ref      => $intervals_ref,
            read_filters_ref   => $read_filters_ref,
            referencefile_path => $referencefile_path,
            temp_directory     => $temp_directory,
            verbosity          => $verbosity,
        }
    );

    ## Add static_quantized_quals
    if ( @{$static_quantized_quals_ref} ) {
        push
          @commands,
          q{--static-quantized-quals}
          . $SPACE
          . join $SPACE
          . q{--static-quantized-quals}
          . $SPACE, @{$static_quantized_quals_ref};
    }

    ## Add BQSR table
    push @commands,
      q{--bqsr-recal-file} . $SPACE . $base_quality_score_recalibration_file;

    ## Output
    push @commands, q{--output} . $SPACE . $outfile_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gatk_applyvqsr {

## Function : Perl wrapper for writing GATK ApplyVQSR recipe to $filehandle. Based on GATK 4.0.8.
## Returns  : @commands
## Arguments: $filehandle           => Sbatch filehandle to write to
##          : $infile_path          => Infile paths
##          : $intervals_ref        => One or more genomic intervals over which to operate {REF}
##          : $java_use_large_pages => Use java large pages
##          : $memory_allocation    => Memory allocation to run Gatk
##          : $mode                 => Mode for emitting reference confidence scores
##          : $outfile_path         => Outfile path
##          : $read_filters_ref     => Filters to apply to reads before analysis {REF}
##          : $recal_file_path      => Output recal file used by ApplyVQSR
##          : $referencefile_path   => Reference sequence file
##          : $stderrfile_path      => Stderrfile path
##          : $temp_directory       => Redirect tmp files to java temp
##          : $tranches_file_path   => Output tranches file used by ApplyRecalibration
##          : $ts_filter_level      => Ts filter level
##          : $verbosity	          => Set the minimum level of logging
##          : $xargs_mode   		    => Set if the program will be executed via xargs

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $intervals_ref;
    my $java_use_large_pages;
    my $memory_allocation;
    my $mode;
    my $outfile_path;
    my $read_filters_ref;
    my $recal_file_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $temp_directory;
    my $tranches_file_path;
    my $ts_filter_level;

    ## Default(s)
    my $pedigree_validation_type;
    my $verbosity;
    my $xargs_mode;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        intervals_ref => {
            default     => [],
            store       => \$intervals_ref,
            strict_type => 1,
        },
        java_use_large_pages => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$java_use_large_pages,
            strict_type => 1,
        },
        memory_allocation => {
            store       => \$memory_allocation,
            strict_type => 1,
        },
        mode => {
            allow       => [qw{ SNP INDEL BOTH }],
            store       => \$mode,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        read_filters_ref => {
            default     => [],
            store       => \$read_filters_ref,
            strict_type => 1,
        },
        recal_file_path => {
            required    => 1,
            store       => \$recal_file_path,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
        tranches_file_path => {
            required    => 1,
            store       => \$tranches_file_path,
            strict_type => 1,
        },
        ts_filter_level => {
            allow       => qr/ \A \d+ \z | \A \d+.\d+ \z /sxm,
            default     => 99.0,
            store       => \$ts_filter_level,
            strict_type => 1,
        },
        verbosity => {
            allow       => [qw{ INFO ERROR FATAL }],
            default     => q{INFO},
            store       => \$verbosity,
            strict_type => 1,
        },
        xargs_mode => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$xargs_mode,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ gatk };

    ## Add java options
    gatk_java_options(
        {
            commands_ref         => \@commands,
            java_use_large_pages => $java_use_large_pages,
            memory_allocation    => $memory_allocation,
            xargs_mode           => $xargs_mode,
        }
    );

    ## Add tool command
    push @commands, q{ApplyVQSR};

    push @commands, q{--variant} . $SPACE . $infile_path;

    ## Add common options
    gatk_common_options(
        {
            commands_ref       => \@commands,
            intervals_ref      => $intervals_ref,
            read_filters_ref   => $read_filters_ref,
            referencefile_path => $referencefile_path,
            temp_directory     => $temp_directory,
            verbosity          => $verbosity,
        }
    );

    if ($mode) {

        push @commands, q{--mode} . $SPACE . $mode;
    }

    push @commands, q{--recal-file} . $SPACE . $recal_file_path;

    push @commands, q{--tranches-file} . $SPACE . $tranches_file_path;

    if ($ts_filter_level) {

        push @commands, q{--truth-sensitivity-filter-level} . $SPACE . $ts_filter_level;
    }

    push @commands, q{--output} . $SPACE . $outfile_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gatk_asereadcounter {

## Function : Perl wrapper for writing GATK ASEReadCounter recipe to $filehandle. Based on GATK 4.0.8.
## Returns  : @commands
## Arguments: $filehandle           => Sbatch filehandle to write to
##          : $infile_path          => Infile path
##          : $intervals_ref        => One or more genomic intervals over which to operate {REF}
##          : $java_use_large_pages => Use java large pages
##          : $memory_allocation    => Memory allocation to run Gatk
##          : $outfile_path         => Outfile path
##          : $read_filters_ref     => Filters to apply on reads {REF}
##          : $referencefile_path   => Reference sequence file
##          : $stderrfile_path      => Stderrfile path
##          : $temp_directory       => Redirect tmp files to java temp
##          : $variant_infile_path  => Sites for which to count allele specific expression
##          : $verbosity            => Set the minimum level of logging
##          : $xargs_mode           => Set if the program will be executed via xargs

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $intervals_ref;
    my $memory_allocation;
    my $outfile_path;
    my $read_filters_ref;
    my $referencefile_path;
    my $stderrfile_path;
    my $temp_directory;
    my $variant_infile_path;

    ## Default(s)
    my $java_use_large_pages;
    my $verbosity;
    my $xargs_mode;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        infile_path => {
            allow       => qr/ (?: bam | sam | cram )$ /xms,
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        intervals_ref => {
            default     => [],
            store       => \$intervals_ref,
            strict_type => 1,
        },
        java_use_large_pages => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$java_use_large_pages,
            strict_type => 1,
        },
        memory_allocation => {
            store       => \$memory_allocation,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        read_filters_ref => {
            default     => [],
            store       => \$read_filters_ref,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
        variant_infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$variant_infile_path,
            strict_type => 1,
        },
        verbosity => {
            allow       => [qw{ INFO ERROR FATAL }],
            default     => q{INFO},
            store       => \$verbosity,
            strict_type => 1,
        },
        xargs_mode => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$xargs_mode,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ gatk };

    ## Add java options
    gatk_java_options(
        {
            commands_ref         => \@commands,
            java_use_large_pages => $java_use_large_pages,
            memory_allocation    => $memory_allocation,
            xargs_mode           => $xargs_mode,
        }
    );

    ## Add tool command
    push @commands, q{ASEReadCounter};

    ## Add infile
    push @commands, q{--input} . $SPACE . $infile_path;

    ## Add variant infile
    push @commands, q{--variant} . $SPACE . $variant_infile_path;

    ## Add common options
    gatk_common_options(
        {
            commands_ref       => \@commands,
            intervals_ref      => $intervals_ref,
            read_filters_ref   => $read_filters_ref,
            referencefile_path => $referencefile_path,
            temp_directory     => $temp_directory,
            verbosity          => $verbosity,
        }
    );

    ## Add output
    if ($outfile_path) {

        push @commands, q{--output} . $SPACE . $outfile_path;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gatk_base {

## Function : Perl wrapper for Gatk base parameters. Based on Gatk 3.7
## Returns  : @commands
## Arguments: $analysis_type                         => Analysis type
##          : $base_quality_score_recalibration_file => Base quality score recalibration file
##          : $commands_ref                          => List of commands added earlier
##          : $disable_indel_qual                    => Disable indel quality
##          : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##          : $filehandle                            => Filehandle to write to
##          : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##          : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##          : $logging_level                         => Logging level
##          : $num_cpu_threads_per_data_thread       => Number of CPU threads to allocate per data thread
##          : $pedigree                              => Pedigree files
##          : $pedigree_validation_type              => Validation strictness for pedigree
##          : $read_filters_ref                      => Filters to apply to reads before analysis {REF}
##          : $referencefile_path                    => Reference sequence file
##          : $static_quantized_quals_ref            => Use static quantized quality scores [ref]

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $analysis_type;
    my $base_quality_score_recalibration_file;
    my $commands_ref;
    my $disable_indel_qual;
    my $downsample_to_coverage;
    my $filehandle;
    my $intervals_ref;
    my $pedigree;
    my $pedigree_validation_type;
    my $read_filters_ref;
    my $referencefile_path;
    my $static_quantized_quals_ref;

    ## Default(s)
    my $gatk_disable_auto_index_and_file_lock;
    my $num_cpu_threads_per_data_thread;
    my $logging_level;

    my $tmpl = {
        analysis_type => {
            defined     => 1,
            required    => 1,
            store       => \$analysis_type,
            strict_type => 1,
        },
        base_quality_score_recalibration_file => {
            store       => \$base_quality_score_recalibration_file,
            strict_type => 1,
        },
        commands_ref => { default => [], store => \$commands_ref, strict_type => 1, },
        disable_indel_qual => {
            allow       => [ 0, 1 ],
            store       => \$disable_indel_qual,
            strict_type => 1,
        },
        downsample_to_coverage => {
            store       => \$downsample_to_coverage,
            strict_type => 1,
        },
        filehandle                            => { store => \$filehandle, },
        gatk_disable_auto_index_and_file_lock => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$gatk_disable_auto_index_and_file_lock,
            strict_type => 1,
        },
        intervals_ref => { default => [], store => \$intervals_ref, strict_type => 1, },
        num_cpu_threads_per_data_thread => {
            allow       => [ undef, qr/ \A \d+ \z /sxm ],
            default     => 0,
            store       => \$num_cpu_threads_per_data_thread,
            strict_type => 1,
        },
        logging_level => {
            allow       => [qw{ INFO ERROR FATAL }],
            default     => q{INFO},
            store       => \$logging_level,
            strict_type => 1,
        },
        pedigree                 => { store => \$pedigree, strict_type => 1, },
        pedigree_validation_type => {
            allow       => [qw{ STRICT SILENT}],
            default     => q{SILENT},
            store       => \$pedigree_validation_type,
            strict_type => 1,
        },
        read_filters_ref =>
          { default => [], store => \$read_filters_ref, strict_type => 1, },
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = @{$commands_ref};

    push @commands, ( qw{ --analysis_type }, $analysis_type, );

    push @commands, q{--logging_level} . $SPACE . $logging_level;

    if ($pedigree_validation_type) {

        push @commands, q{--pedigreeValidationType} . $SPACE . $pedigree_validation_type;
    }

    if ($pedigree) {

        push @commands, q{--pedigree} . $SPACE . $pedigree;
    }

    if ($num_cpu_threads_per_data_thread) {

        push @commands,
            q{--num_cpu_threads_per_data_thread}
          . $SPACE
          . $num_cpu_threads_per_data_thread;
    }

    if ($downsample_to_coverage) {

        push @commands, q{--downsample_to_coverage} . $SPACE . $downsample_to_coverage;
    }

    if ($gatk_disable_auto_index_and_file_lock) {

        push @commands, q{--disable_auto_index_creation_and_locking_when_reading_rods};
    }

    if ( @{$intervals_ref} ) {

        push @commands,
          q{--intervals} . $SPACE . join $SPACE . q{--intervals} . $SPACE,
          @{$intervals_ref};
    }

    if ( @{$read_filters_ref} ) {

        push @commands,
          q{--read_filter} . $SPACE . join $SPACE . q{--read_filter} . $SPACE,
          @{$read_filters_ref};
    }

    if ($referencefile_path) {

        push @commands, q{--reference_sequence} . $SPACE . $referencefile_path;
    }

    if ($base_quality_score_recalibration_file) {

        push @commands, q{--BQSR} . $SPACE . $base_quality_score_recalibration_file;
    }

    if ($disable_indel_qual) {

        push @commands, q{--disable_indel_quals};
    }

    if ( @{$static_quantized_quals_ref} ) {

        push
          @commands,
          q{--static_quantized_quals}
          . $SPACE
          . join $SPACE
          . q{--static_quantized_quals}
          . $SPACE, @{$static_quantized_quals_ref};
    }
    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gatk_baserecalibrator {

## Function : Perl wrapper for writing GATK BaseRecalibrator recipe to $filehandle. Based on GATK 4.0.8
## Returns  : @commands
## Arguments: $filehandle           => Sbatch filehandle to write to
##          : $infile_path          => Infile paths
##          : $intervals_ref        => One or more genomic intervals over which to operate {REF}
##          : $java_use_large_pages => Use java large pages
##          : $known_alleles_ref    => Input VCF file(s) with known indels {REF}
##          : $memory_allocation    => Memory allocation to run Gatk
##          : $outfile_path         => Outfile path
##          : $pedigree             => Pedigree files
##          : $read_filters_ref     => Filters to apply on reads {REF}
##          : $referencefile_path   => Reference sequence file
##          : $stderrfile_path      => Stderrfile path
##          : $temp_directory       => Redirect tmp files to java temp
##          : $verbosity            => Set the minimum level of logging
##          : $xargs_mode           => Set if the program will be executed via xargs

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $intervals_ref;
    my $known_sites_ref;
    my $memory_allocation;
    my $outfile_path;
    my $read_filters_ref;
    my $referencefile_path;
    my $static_quantized_quals_ref;
    my $stderrfile_path;
    my $temp_directory;

    ## Default(s)
    my $java_use_large_pages;
    my $verbosity;
    my $xargs_mode;

    my $tmpl = {
        filehandle  => { store => \$filehandle, },
        infile_path => {
            allow       => qr/ (?: bam | sam | cram )$ /xms,
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
        java_use_large_pages => {
            allow       => [ undef, 0, 1 ],
            store       => \$java_use_large_pages,
            strict_type => 1,
        },
        known_sites_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$known_sites_ref,
            strict_type => 1,
        },
        memory_allocation => {
            store       => \$memory_allocation,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        read_filters_ref => {
            default     => [],
            store       => \$read_filters_ref,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
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
        verbosity       => {
            allow       => [qw{ INFO ERROR FATAL }],
            default     => q{INFO},
            store       => \$verbosity,
            strict_type => 1,
        },
        xargs_mode => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$xargs_mode,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ gatk };

    ## Add java options
    gatk_java_options(
        {
            commands_ref         => \@commands,
            java_use_large_pages => $java_use_large_pages,
            memory_allocation    => $memory_allocation,
            xargs_mode           => $xargs_mode,
        }
    );

    ## Add tool command
    push @commands, q{BaseRecalibrator};

    ## Add infile
    push @commands, q{--input} . $SPACE . $infile_path;

    ## Add common options
    gatk_common_options(
        {
            commands_ref       => \@commands,
            intervals_ref      => $intervals_ref,
            read_filters_ref   => $read_filters_ref,
            referencefile_path => $referencefile_path,
            temp_directory     => $temp_directory,
            verbosity          => $verbosity,
        }
    );

    ## Add known sites reference
    push @commands,
      q{--known-sites} . $SPACE . join $SPACE . q{--known-sites} . $SPACE,
      @{$known_sites_ref};

    ## Output
    push @commands, q{--output} . $SPACE . $outfile_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gatk_calculategenotypeposteriors {

## Function : Perl wrapper for writing GATK CalculateGenotypePosteriors recipe to $filehandle. Based on GATK 4.1.0.
## Returns  : @commands
##          : $filehandle                   => Sbatch filehandle to write to
##          : $infile_path                  => Infile paths
##          : $intervals_ref                => One or more genomic intervals over which to operate {REF}
##          : $java_use_large_pages         => Use java large pages
##          : $memory_allocation            => Memory allocation to run Gatk
##          : $num_ref_samples_if_no_call   => Number of hom-ref genotypes to infer at sites not present in a panel
##          : $outfile_path                 => Outfile path
##          : $pedigree                     => Pedigree files for samples
##          : $read_filters_ref             => Filters to apply on reads {REF}
##          : $referencefile_path           => Reference sequence file
##          : $stderrfile_path              => Stderrfile path
##          : $supporting_callset_file_path => Other callsets to use in generating genotype posteriors
##          : $temp_directory               => Redirect tmp files to java temp
##          : $verbosity                    => Set the minimum level of logging
##          : $xargs_mode                   => Set if the program will be executed via xargs

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $intervals_ref;
    my $memory_allocation;
    my $num_ref_samples_if_no_call;
    my $outfile_path;
    my $pedigree;
    my $read_filters_ref;
    my $referencefile_path;
    my $stderrfile_path;
    my $supporting_callset_file_path;
    my $temp_directory;

    ## Default(s)
    my $java_use_large_pages;
    my $verbosity;
    my $xargs_mode;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        intervals_ref => {
            default     => [],
            store       => \$intervals_ref,
            strict_type => 1,
        },
        java_use_large_pages => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$java_use_large_pages,
            strict_type => 1,
        },
        memory_allocation => {
            store       => \$memory_allocation,
            strict_type => 1,
        },
        num_ref_samples_if_no_call => {
            allow       => [ undef, qr/ \A \d+ \z /sxm ],
            store       => \$num_ref_samples_if_no_call,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        pedigree => {
            store       => \$pedigree,
            strict_type => 1,
        },
        read_filters_ref => {
            default     => [],
            store       => \$read_filters_ref,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        supporting_callset_file_path => {
            store       => \$supporting_callset_file_path,
            strict_type => 1,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
        verbosity => {
            allow       => [qw{ INFO ERROR FATAL }],
            default     => q{INFO},
            store       => \$verbosity,
            strict_type => 1,
        },
        xargs_mode => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$xargs_mode,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ gatk };

    ## Add java options
    gatk_java_options(
        {
            commands_ref         => \@commands,
            java_use_large_pages => $java_use_large_pages,
            memory_allocation    => $memory_allocation,
            xargs_mode           => $xargs_mode,
        }
    );

    push @commands, q{CalculateGenotypePosteriors};

    push @commands, q{--variant} . $SPACE . $infile_path;

    ## Add common options
    gatk_common_options(
        {
            commands_ref       => \@commands,
            intervals_ref      => $intervals_ref,
            pedigree           => $pedigree,
            read_filters_ref   => $read_filters_ref,
            referencefile_path => $referencefile_path,
            temp_directory     => $temp_directory,
            verbosity          => $verbosity,
        }
    );

    if ($num_ref_samples_if_no_call) {

        push @commands,
          q{--num-reference-samples-if-no-call} . $SPACE . $num_ref_samples_if_no_call;
    }

    if ($supporting_callset_file_path) {

        push @commands, q{--supporting-callsets} . $SPACE . $supporting_callset_file_path;
    }

    push @commands, q{--output} . $SPACE . $outfile_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gatk_cnnscorevariants {

## Function : Perl wrapper for writing GATK CNNScoreVariants recipe to $filehandle. Based on GATK 4.1.0
## Returns  : @commands
## Arguments: $alignment_infile_paths_ref => BAM/SAM/CRAM file containing reads {REF}
##          : $filehandle                 => Sbatch filehandle to write to
##          : $infile_path                => Infile paths
##          : $intervals_ref              => One or more genomic intervals over which to operate {REF}
##          : $java_use_large_pages       => Use java large pages
##          : $memory_allocation          => Memory allocation to run Gatk
##          : $outfile_path               => Outfile path
##          : $pedigree                   => Pedigree files for samples
##          : $referencefile_path         => Reference sequence file
##          : $stderrfile_path            => Stderrfile path
##          : $temp_directory             => Redirect tmp files to java temp
##          : $verbosity                  => Set the minimum level of logging
##          : $xargs_mode                 => Set if the program will be executed via xargs

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $alignment_infile_paths_ref;
    my $filehandle;
    my $infile_path;
    my $intervals_ref;
    my $memory_allocation;
    my $outfile_path;
    my $pedigree;
    my $referencefile_path;
    my $stderrfile_path;
    my $temp_directory;
    my $xargs_mode;

    ## Default(s)
    my $java_use_large_pages;
    my $verbosity;

    my $tmpl = {
        alignment_infile_paths_ref => {
            default     => [],
            store       => \$alignment_infile_paths_ref,
            strict_type => 1,
        },
        filehandle  => { store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        intervals_ref => {
            default     => [],
            store       => \$intervals_ref,
            strict_type => 1,
        },
        java_use_large_pages => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$java_use_large_pages,
            strict_type => 1,
        },
        memory_allocation => {
            store       => \$memory_allocation,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        pedigree => {
            store       => \$pedigree,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
        verbosity => {
            allow       => [qw{ INFO ERROR FATAL }],
            default     => q{INFO},
            store       => \$verbosity,
            strict_type => 1,
        },
        xargs_mode => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$xargs_mode,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ gatk };

    ## Add java options
    gatk_java_options(
        {
            commands_ref         => \@commands,
            java_use_large_pages => $java_use_large_pages,
            memory_allocation    => $memory_allocation,
            xargs_mode           => $xargs_mode,
        }
    );

    ## Add tool command
    push @commands, q{CNNScoreVariants};

    push @commands, q{--variant} . $SPACE . $infile_path;

    ## Add common options
    gatk_common_options(
        {
            commands_ref       => \@commands,
            intervals_ref      => $intervals_ref,
            pedigree           => $pedigree,
            referencefile_path => $referencefile_path,
            temp_directory     => $temp_directory,
            verbosity          => $verbosity,
        }
    );

    if ( @{$alignment_infile_paths_ref} ) {

        push @commands, q{--input} . $SPACE . join $SPACE . q{--input} . $SPACE,
          @{$alignment_infile_paths_ref};

        ## Modify tensor-type accordingly
        push @commands, q{--tensor-type} . $SPACE . q{read_tensor};
    }

    push @commands, q{--output} . $SPACE . $outfile_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gatk_combinevariants {

## Function : Perl wrapper for writing GATK combinevariants recipe to $filehandle. Based on GATK 3.7.0.
## Returns  : @commands
## Arguments: $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##          : $exclude_nonvariants                   => Exclude non-variant sites
##          : $filehandle                            => Sbatch filehandle to write to
##          : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##          : $genotype_merge_option                 => Determines how we should merge genotype records for samples shared across the ROD files
##          : $infile_paths_ref                      => Infile paths {REF}
##          : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##          : $java_jar                              => Java jar
##          : $java_use_large_pages                  => Use java large pages
##          : $logging_level                         => Set the minimum level of logging
##          : $memory_allocation                     => Memory allocation to run Gatk
##          : $outfile_path                          => Outfile path
##          : $pedigree                              => Pedigree files for samples
##          : $pedigree_validation_type              => Validation strictness for pedigree
##          : $prioritize_caller                     => Comma seperated string specifying priority for merging
##          : $referencefile_path                    => Reference sequence file
##          : $stderrfile_path                       => Stderrfile path
##          : $temp_directory                        => Redirect tmp files to java temp

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $downsample_to_coverage;
    my $filehandle;
    my $genotype_merge_option;
    my $infile_paths_ref;
    my $intervals_ref;
    my $java_jar;
    my $java_use_large_pages;
    my $memory_allocation;
    my $outfile_path;
    my $pedigree;
    my $prioritize_caller;
    my $referencefile_path;
    my $stderrfile_path;
    my $temp_directory;

    ## Default(s)
    my $exclude_nonvariants;
    my $gatk_disable_auto_index_and_file_lock;
    my $logging_level;
    my $pedigree_validation_type;

    my $tmpl = {
        downsample_to_coverage => {
            allow       => qr/ \A \d+ \z /sxm,
            store       => \$downsample_to_coverage,
            strict_type => 1,
        },
        exclude_nonvariants => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$exclude_nonvariants,
            strict_type => 1,
        },
        filehandle                            => { store => \$filehandle, },
        gatk_disable_auto_index_and_file_lock => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$gatk_disable_auto_index_and_file_lock,
            strict_type => 1,
        },
        genotype_merge_option => {
            allow       => [ undef, qw{ UNIQUIFY PRIORITIZE UNSORTED REQUIRE_UNIQUE } ],
            default     => q{PRIORITIZE},
            store       => \$genotype_merge_option,
            strict_type => 1,
        },
        infile_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infile_paths_ref,
            strict_type => 1,
        },
        intervals_ref => { default => [], store => \$intervals_ref, strict_type => 1, },
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
        memory_allocation => { store => \$memory_allocation, strict_type => 1, },
        outfile_path      => {
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
        prioritize_caller  => { store => \$prioritize_caller, strict_type => 1, },
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

    @commands = gatk_base(
        {
            analysis_type          => q{CombineVariants},
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

    ## Add binary to beginning
    unshift @commands, q{gatk3};

    if ($exclude_nonvariants) {

        push @commands, q{--excludeNonVariants};
    }

    if ($genotype_merge_option) {

        push @commands, q{--genotypemergeoption} . $SPACE . $genotype_merge_option;
    }

    if ($prioritize_caller) {

        push @commands, q{--rod_priority_list} . $SPACE . $prioritize_caller;
    }

    if ( @{$infile_paths_ref} ) {

        push @commands, q{--variant:} . join $SPACE . q{--variant:}, @{$infile_paths_ref};

    }

    if ($outfile_path) {

        push @commands, q{--out} . $SPACE . $outfile_path;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gatk_common_options {

## Function : Perl wrapper for adding common GATK options to commands_ref. Based on GATK 4.0.
## Returns  : $commands_ref
## Arguments: $commands_ref       => List of commands added earlier {REF}
##          : $intervals_ref      => One or more genomic intervals over which to operate {REF}
##          : $pedigree           => Pedigree file
##          : $read_filters_ref   => Filters to apply to reads before analysis {REF}
##          : $referencefile_path => Path to reference sequence file
##          : $temp_dir_path      => Path to temporary directory to use
##          : $verbosity          => Verbosity level

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $commands_ref;
    my $intervals_ref;
    my $pedigree;
    my $read_filters_ref;
    my $referencefile_path;
    my $temp_directory;
    my $verbosity;

    my $tmpl = {
        commands_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$commands_ref,
            strict_type => 1,
        },
        intervals_ref => {
            default     => [],
            defined     => 1,
            store       => \$intervals_ref,
            strict_type => 1,
        },
        pedigree => {
            store       => \$pedigree,
            strict_type => 1,
        },
        read_filters_ref => {
            default     => [],
            store       => \$read_filters_ref,
            strict_type => 1,
        },
        referencefile_path => {
            store       => \$referencefile_path,
            strict_type => 1,
        },
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
        verbosity => {
            allow       => [ undef, qw{ ERROR INFO WARNING DEBUG } ],
            store       => \$verbosity,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    if ( @{$intervals_ref} ) {

        push @{$commands_ref},
          q{--intervals} . $SPACE . join $SPACE . q{--intervals} . $SPACE,
          @{$intervals_ref};
    }

    if ($pedigree) {

        push @{$commands_ref}, q{--pedigree} . $SPACE . $pedigree;
    }

    if ( @{$read_filters_ref} ) {

        push @{$commands_ref},
          q{--read-filter} . $SPACE . join $SPACE . q{--read-filter} . $SPACE,
          @{$read_filters_ref};
    }

    if ($referencefile_path) {

        push @{$commands_ref}, q{--reference} . $SPACE . $referencefile_path;
    }

    if ($temp_directory) {

        push @{$commands_ref}, q{--tmp-dir} . $SPACE . $temp_directory;
    }

    if ($verbosity) {

        push @{$commands_ref}, q{--verbosity} . $SPACE . $verbosity;
    }
    return @{$commands_ref};
}

sub gatk_concatenate_variants {

## Function : Writes sbatch code to supplied filehandle to concatenate variants in vcf format. Each array element is combined with the infile prefix and postfix.
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $continue              => Adds an ampersand to the end of the command
##          : $filehandle            => SBATCH script filehandle to print to
##          : $elements_ref          => Holding the number and part of file names to be combined
##          : $infile_prefix         => Will be combined with the each array element
##          : $infile_postfix        => Will be combined with the each array element
##          : $outfile_path_prefix   => Combined outfile path prefix
##          : $outfile_suffix        => Combined outfile suffix

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $continue;
    my $elements_ref;
    my $filehandle;
    my $infile_prefix;
    my $infile_postfix;
    my $outfile_path_prefix;
    my $outfile_suffix;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        continue => {
            allow => [ 0, 1 ],
            store => \$continue,
        },
        elements_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$elements_ref,
            strict_type => 1,
        },
        filehandle => {
            defined  => 1,
            required => 1,
            store    => \$filehandle,
        },
        infile_prefix => {
            defined     => 1,
            required    => 1,
            store       => \$infile_prefix,
            strict_type => 1,
        },
        infile_postfix => {
            store       => \$infile_postfix,
            strict_type => 1,
        },
        outfile_path_prefix => {
            store       => \$outfile_path_prefix,
            strict_type => 1,
        },
        outfile_suffix => {
            allow       => [qw{ .all.vcf .selected.vcf .vcf }],
            default     => q{.vcf},
            store       => \$outfile_suffix,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Gatk qw(gatk_gathervcfscloud);

    ## Outfile path to be built
    my $outfile_path;

    ## No postfix
    if ( not defined $infile_postfix ) {
        $infile_postfix = $EMPTY_STR;
    }

    ## Build $outfile_path
    if ( defined $outfile_path_prefix ) {
        $outfile_path = $outfile_path_prefix . $outfile_suffix;
    }
    else {
        $outfile_path = $infile_prefix . $outfile_suffix;
    }

    say {$filehandle} q{## GATK GatherVCFs};

    ## Assemble infile paths
    my @infile_paths =
      map { $infile_prefix . $DOT . $_ . $infile_postfix } @{$elements_ref};

    gatk_gathervcfscloud(
        {
            filehandle           => $filehandle,
            infile_paths_ref     => \@infile_paths,
            java_use_large_pages => $active_parameter_href->{java_use_large_pages},
            memory_allocation    => q{Xmx4g},
            outfile_path         => $outfile_path,
            temp_directory       => $active_parameter_href->{temp_directory},
            verbosity            => $active_parameter_href->{gatk_logging_level},
        }
    );

    ## Launch multiple processes
    if ($continue) {

        print {$filehandle} $AMPERSAND;
    }
    say {$filehandle} $NEWLINE;
    return;
}

sub gatk_gatherbqsrreports {

## Function : Perl wrapper for writing GATK GatherBQSRReports recipe to $filehandle. Based on GATK 4.0.11
## Returns  : @commands
## Arguments: $base_quality_score_recalibration_files_ref => Input recalibration table for BQSR {REF}
##          : $filehandle                                 => Sbatch filehandle to write to
##          : $infile_path                                => Infile paths
##          : $intervals_ref                              => One or more genomic intervals over which to operate {REF}
##          : $java_use_large_pages                       => Use java large pages
##          : $memory_allocation                          => Memory allocation to run Gatk
##          : $outfile_path                               => Outfile path
##          : $read_filters_ref                           => Filters to apply on reads {REF}
##          : $referencefile_path                         => Reference sequence file
##          : $static_quantized_quals_ref                 => Use static quantized quality scores to a given number of levels {REF}
##          : $stderrfile_path                            => Stderrfile path
##          : $temp_directory                             => Redirect tmp files to java temp
##          : $verbosity                                  => Set the minimum level of logging
##          : $xargs_mode                                 => Set if the program will be executed via xargs

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $base_quality_score_recalibration_files_ref;
    my $filehandle;
    my $infile_path;
    my $intervals_ref;
    my $memory_allocation;
    my $outfile_path;
    my $read_filters_ref;
    my $referencefile_path;
    my $stderrfile_path;
    my $temp_directory;

    ## Default(s)
    my $java_use_large_pages;
    my $verbosity;
    my $xargs_mode;

    my $tmpl = {
        base_quality_score_recalibration_files_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$base_quality_score_recalibration_files_ref,
            strict_type => 1,
        },
        filehandle    => { store => \$filehandle, },
        intervals_ref => {
            default     => [],
            defined     => 1,
            store       => \$intervals_ref,
            strict_type => 1,
        },
        java_use_large_pages => {
            allow       => [ undef, 0, 1 ],
            store       => \$java_use_large_pages,
            strict_type => 1,
        },
        memory_allocation => {
            store       => \$memory_allocation,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        read_filters_ref => {
            default     => [],
            store       => \$read_filters_ref,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        temp_directory  => { store => \$temp_directory,  strict_type => 1, },
        verbosity       => {
            allow       => [qw{ INFO ERROR FATAL }],
            default     => q{INFO},
            store       => \$verbosity,
            strict_type => 1,
        },
        xargs_mode => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$xargs_mode,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ gatk };

    ## Add java options
    gatk_java_options(
        {
            commands_ref         => \@commands,
            java_use_large_pages => $java_use_large_pages,
            memory_allocation    => $memory_allocation,
            xargs_mode           => $xargs_mode,
        }
    );

    ## Add tool command
    push @commands, q{GatherBQSRReports};

    ## Add infile
    push @commands, q{--input} . $SPACE . join $SPACE . q{--input} . $SPACE,
      @{$base_quality_score_recalibration_files_ref};

    ## Add common options
    gatk_common_options(
        {
            commands_ref       => \@commands,
            intervals_ref      => $intervals_ref,
            read_filters_ref   => $read_filters_ref,
            referencefile_path => $referencefile_path,
            temp_directory     => $temp_directory,
            verbosity          => $verbosity,
        }
    );

    ## Output
    push @commands, q{--output} . $SPACE . $outfile_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gatk_gathervcfscloud {

## Function : Perl wrapper for writing GATK GatherVcfsCloud recipe to $filehandle. Based on GATK 4.0.8.
## Returns  : @commands
## Arguments: $filehandle           => Sbatch filehandle to write to
##          : $ignore_safety_checks => Disable sanity checks to improve performance
##          : $infile_paths_ref     => VCF files to gather {REF}
##          : $java_use_large_pages => Use java large pages
##          : $memory_allocation    => Memory allocation to run Gatk
##          : $outfile_path         => Outfile path
##          : $stderrfile_path      => Stderrfile path
##          : $temp_directory       => Redirect tmp files to java temp
##          : $verbosity            => Set the minimum level of logging

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $ignore_safety_checks;
    my $infile_paths_ref;
    my $memory_allocation;
    my $outfile_path;
    my $stderrfile_path;
    my $temp_directory;

    ## Default(s)
    my $java_use_large_pages;
    my $verbosity;

    my $tmpl = {
        filehandle           => { store => \$filehandle, },
        ignore_safety_checks => {
            allow       => [ undef, 0, 1 ],
            store       => \$ignore_safety_checks,
            strict_type => 1,
        },
        infile_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infile_paths_ref,
            strict_type => 1,
        },
        java_use_large_pages => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$java_use_large_pages,
            strict_type => 1,
        },
        memory_allocation => {
            store       => \$memory_allocation,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
        verbosity => {
            allow       => [qw{ INFO ERROR FATAL }],
            default     => q{INFO},
            store       => \$verbosity,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ gatk };

    ## Add java options
    gatk_java_options(
        {
            commands_ref         => \@commands,
            java_use_large_pages => $java_use_large_pages,
            memory_allocation    => $memory_allocation,
        }
    );

    push @commands, q{GatherVcfsCloud};

    push @commands,
      q{--input} . $SPACE . join $SPACE . q{--input} . $SPACE,
      @{$infile_paths_ref};

    ## Add common options
    gatk_common_options(
        {
            commands_ref   => \@commands,
            temp_directory => $temp_directory,
            verbosity      => $verbosity,
        }
    );

    if ($ignore_safety_checks) {

        push @commands, q{--ignore-safety-checks};
    }

    push @commands, q{--output} . $SPACE . $outfile_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gatk_genomicsdbimport {

## Function : Perl wrapper for writing GATK GenomicsDBImport recipe to $filehandle. Based on GATK 4.0.8
## Returns  : @commands
## Arguments: $filehandle                => Sbatch filehandle to write to
##          : $genomicsdb_workspace_path => Workspace for GenomicsDB
##          : $infile_paths_ref          => GVCF files to be imported to GenomicsDB {REF}
##          : $intervals_ref             => One or more genomic intervals over which to operate {REF}
##          : $java_use_large_pages      => Use java large pages
##          : $memory_allocation         => Memory allocation to run Gatk
##          : $pedigree                  => Pedigree files
##          : $read_filters_ref          => Filters to apply on reads {REF}
##          : $referencefile_path        => Reference sequence file
##          : $sample_name_map_path      => Sample name map for merged references (Format: sample_id\tfile.vcf )
##          : $stderrfile_path           => Stderrfile path
##          : $temp_directory            => Redirect tmp files to java temp
##          : $verbosity                 => Set the minimum level of logging
##          : $xargs_mode                => Set if the program will be executed via xargs

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $genomicsdb_workspace_path;
    my $infile_paths_ref;
    my $intervals_ref;
    my $java_use_large_pages;
    my $memory_allocation;
    my $read_filters_ref;
    my $referencefile_path;
    my $sample_name_map_path;
    my $stderrfile_path;
    my $temp_directory;

    ## Default(s)
    my $verbosity;
    my $xargs_mode;

    my $tmpl = {
        filehandle                => { store => \$filehandle, },
        genomicsdb_workspace_path => {
            defined     => 1,
            required    => 1,
            store       => \$genomicsdb_workspace_path,
            strict_type => 1,
        },
        infile_paths_ref => {
            default     => [],
            store       => \$infile_paths_ref,
            strict_type => 1,
        },
        intervals_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$intervals_ref,
            strict_type => 1,
        },
        java_use_large_pages => {
            allow       => [ undef, 0, 1 ],
            store       => \$java_use_large_pages,
            strict_type => 1,
        },
        memory_allocation => {
            store       => \$memory_allocation,
            strict_type => 1,
        },
        read_filters_ref => {
            default     => [],
            store       => \$read_filters_ref,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        sample_name_map_path => {
            store       => \$sample_name_map_path,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        temp_directory  => { store => \$temp_directory,  strict_type => 1, },
        verbosity       => {
            allow       => [qw{ INFO ERROR FATAL }],
            default     => q{INFO},
            store       => \$verbosity,
            strict_type => 1,
        },
        xargs_mode => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$xargs_mode,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ gatk };

    ## Add java options
    gatk_java_options(
        {
            commands_ref         => \@commands,
            java_use_large_pages => $java_use_large_pages,
            memory_allocation    => $memory_allocation,
            xargs_mode           => $xargs_mode,
        }
    );

    push @commands, q{GenomicsDBImport};

    if ( scalar @{$infile_paths_ref} ) {
        push @commands,
          q{--variant} . $SPACE . join $SPACE . q{--variant} . $SPACE,
          @{$infile_paths_ref};
    }

    if ($sample_name_map_path) {

        push @commands, q{--sample-name-map} . $SPACE . $sample_name_map_path;
    }

    ## Add common options
    gatk_common_options(
        {
            commands_ref       => \@commands,
            intervals_ref      => $intervals_ref,
            referencefile_path => $referencefile_path,
            temp_directory     => $temp_directory,
            verbosity          => $verbosity,
        }
    );

    push @commands, q{--genomicsdb-workspace-path} . $SPACE . $genomicsdb_workspace_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gatk_genotypegvcfs {

## Function : Perl wrapper for writing GATK GenoTypeGVCFs recipe to $filehandle. Based on GATK 4.1.0.
## Returns  : @commands
## Arguments: $dbsnp_path               => Path to DbSNP file
##          : $filehandle               => Sbatch filehandle to write to
##          : $infile_path              => Path to variant input
##          : $include_nonvariant_sites => Include loci found to be non-variant after genotyping
##          : $intervals_ref            => One or more genomic intervals over which to operate {REF}
##          : $java_use_large_pages     => Use java large pages
##          : $memory_allocation        => Memory allocation to run Gatk
##          : $outfile_path             => Outfile path
##          : $pedigree                 => Pedigree files for samples
##          : $referencefile_path       => Reference sequence file
##          : $stderrfile_path          => Stderrfile path
##          : $temp_directory           => Redirect tmp files to java temp
##          : $verbosity                => Set the minimum level of logging
##          : $xargs_mode               => Set if the program will be executed via xargs

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $dbsnp_path;
    my $filehandle;
    my $infile_path;
    my $intervals_ref;
    my $memory_allocation;
    my $outfile_path;
    my $pedigree;
    my $referencefile_path;
    my $stderrfile_path;
    my $temp_directory;

    ## Default(s)
    my $include_nonvariant_sites;
    my $java_use_large_pages;
    my $verbosity;
    my $xargs_mode;

    my $tmpl = {
        dbsnp_path => {
            store       => \$dbsnp_path,
            strict_type => 1,
        },
        filehandle               => { store => \$filehandle, },
        include_nonvariant_sites => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$include_nonvariant_sites,
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
            store       => \$intervals_ref,
            strict_type => 1,
        },
        java_use_large_pages => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$java_use_large_pages,
            strict_type => 1,
        },
        memory_allocation => {
            store       => \$memory_allocation,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        pedigree => {
            store       => \$pedigree,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
            required    => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
        verbosity => {
            allow       => [qw{ INFO ERROR FATAL }],
            default     => q{INFO},
            store       => \$verbosity,
            strict_type => 1,
        },
        xargs_mode => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$xargs_mode,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ gatk };

    ## Add java options
    gatk_java_options(
        {
            commands_ref         => \@commands,
            java_use_large_pages => $java_use_large_pages,
            memory_allocation    => $memory_allocation,
            xargs_mode           => $xargs_mode,
        }
    );

    ## Add tool command
    push @commands, q{GenotypeGVCFs};

    push @commands, q{--variant} . $SPACE . $infile_path;

    ## Add common options
    gatk_common_options(
        {
            commands_ref       => \@commands,
            intervals_ref      => $intervals_ref,
            pedigree           => $pedigree,
            referencefile_path => $referencefile_path,
            temp_directory     => $temp_directory,
            verbosity          => $verbosity,
        }
    );

    if ($include_nonvariant_sites) {

        push @commands, q{--include-non-variant-sites};
    }

    if ($dbsnp_path) {
        push @commands, q{--dbsnp} . $SPACE . $dbsnp_path;
    }

    push @commands, q{--output} . $SPACE . $outfile_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gatk_haplotypecaller {

## Function : Perl wrapper for writing GATK haplotypecaller recipe to $filehandle. Based on GATK 4.1.8.1.
## Returns  : @commands
## Arguments: $annotations_ref                               => One or more specific annotations to apply to variant calls
##          : $dbsnp_path                                    => Path to DbSNP file
##          : $dont_use_soft_clipped_bases                   => Do not analyze soft clipped bases in the reads
##          : $emit_ref_confidence                           => Mode for emitting reference confidence scores
##          : $filehandle                                    => Sbatch filehandle to write to
##          : $infile_path                                   => Infile paths
##          : $intervals_ref                                 => One or more genomic intervals over which to operate {REF}
##          : $java_use_large_pages                          => Use java large pages
###         : $linked_de_bruijn_graph                        => Use linked de bruijn graph assembly mode
##          : $memory_allocation                             => Memory allocation to run Gatk
##          : $outfile_path                                  => Outfile path
##          : $pcr_indel_model                               => The PCR indel model to use
##          : $pedigree                                      => Pedigree files for samples
##          : $read_filters_ref                              => Filters to apply to reads before analysis {REF}
##          : $referencefile_path                            => Reference sequence file
##          : $sample_ploidy                                 => Ploidy per sample
##          : $standard_min_confidence_threshold_for_calling => The minimum phred-scaled confidence threshold at which variants should be called
##          : $stderrfile_path                               => Stderrfile path
##          : $temp_directory                                => Redirect tmp files to java temp
##          : $verbosity                                     => Set the minimum level of logging
##          : $xargs_mode                                    => Set if the program will be executed via xargs

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $annotations_ref;
    my $dbsnp_path;
    my $dont_use_soft_clipped_bases;
    my $filehandle;
    my $infile_path;
    my $intervals_ref;
    my $memory_allocation;
    my $outfile_path;
    my $pcr_indel_model;
    my $pedigree;
    my $read_filters_ref;
    my $referencefile_path;
    my $sample_ploidy;
    my $standard_min_confidence_threshold_for_calling;
    my $stderrfile_path;
    my $temp_directory;

    ## Default(s)
    my $emit_ref_confidence;
    my $java_use_large_pages;
    my $linked_de_bruijn_graph;
    my $verbosity;
    my $xargs_mode;

    my $tmpl = {
        annotations_ref => {
            default     => [],
            defined     => 1,
            store       => \$annotations_ref,
            strict_type => 1,
        },
        dbsnp_path => {
            store       => \$dbsnp_path,
            strict_type => 1,
        },
        dont_use_soft_clipped_bases => {
            allow       => [ undef, 0, 1 ],
            store       => \$dont_use_soft_clipped_bases,
            strict_type => 1,
        },
        emit_ref_confidence => {
            allow       => [qw{ NONE BP_RESOLUTION GVCF }],
            default     => q{GVCF},
            store       => \$emit_ref_confidence,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        infile_path => {
            allow       => qr/ (?: bam | sam | cram )$ /xms,
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
        java_use_large_pages => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$java_use_large_pages,
            strict_type => 1,
        },
        linked_de_bruijn_graph => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$linked_de_bruijn_graph,
            strict_type => 1,
        },
        memory_allocation => {
            store       => \$memory_allocation,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        pedigree => {
            store       => \$pedigree,
            strict_type => 1,
        },
        pcr_indel_model => {
            allow       => [ undef, qw{ NONE HOSTILE AGGRESSIVE CONSERVATIVE } ],
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
            store       => \$referencefile_path,
            strict_type => 1,
        },
        sample_ploidy => {
            allow       => [ undef, qr/ ^\d+$ /sxm ],
            store       => \$sample_ploidy,
            strict_type => 1,
        },
        standard_min_confidence_threshold_for_calling => {
            allow       => qr/ ^\d+$ /sxm,
            store       => \$standard_min_confidence_threshold_for_calling,
            strict_type => 1,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
        verbosity => {
            allow       => [qw{ INFO ERROR FATAL }],
            default     => q{INFO},
            store       => \$verbosity,
            strict_type => 1,
        },
        xargs_mode => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$xargs_mode,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ gatk };

    ## Add java options
    gatk_java_options(
        {
            commands_ref         => \@commands,
            java_use_large_pages => $java_use_large_pages,
            memory_allocation    => $memory_allocation,
            xargs_mode           => $xargs_mode,
        }
    );

    ## Add tool command
    push @commands, q{HaplotypeCaller};

    ## Add infile
    push @commands, q{--input} . $SPACE . $infile_path;

    ## Add common options
    gatk_common_options(
        {
            commands_ref       => \@commands,
            intervals_ref      => $intervals_ref,
            pedigree           => $pedigree,
            read_filters_ref   => $read_filters_ref,
            referencefile_path => $referencefile_path,
            temp_directory     => $temp_directory,
            verbosity          => $verbosity,
        }
    );

    ## Add annotaions
    if ( @{$annotations_ref} ) {
        push @commands,
          q{--annotation} . $SPACE . join $SPACE . q{--annotation} . $SPACE,
          @{$annotations_ref};
    }

    ## Add dbsnp
    if ($dbsnp_path) {

        push @commands, q{--dbsnp} . $SPACE . $dbsnp_path;
    }

    if ($linked_de_bruijn_graph) {
        push @commands, q{--linked-de-bruijn-graph};
    }

    ## No soft clipped bases
    if ($dont_use_soft_clipped_bases) {

        push @commands, q{--dont-use-soft-clipped-bases};
    }

    ## Set output mode
    push @commands, q{--emit-ref-confidence} . $SPACE . $emit_ref_confidence;

    ## Add PCR indel model
    if ($pcr_indel_model) {

        push @commands, q{--pcr-indel-model} . $SPACE . $pcr_indel_model;
    }

    ## Add sample ploidy
    if ($sample_ploidy) {

        push @commands, q{--sample-ploidy} . $SPACE . $sample_ploidy;
    }

    ## Add minimum phred-scaled confidence threshold
    if ($standard_min_confidence_threshold_for_calling) {
        push @commands,
            q{--standard-min-confidence-threshold-for-calling}
          . $SPACE
          . $standard_min_confidence_threshold_for_calling;
    }

    ## Output
    push @commands, q{--output} . $SPACE . $outfile_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gatk_indexfeaturefile {

## Function : Perl wrapper for writing GATK IndexFeatureFile recipe to $filehandle. Based on GATK 4.1.6.
## Returns  : @commands
## Arguments: $filehandle           => Filehandle to write to
##          : $infile_path          => Path to feature file
##          : $java_use_large_pages => Use java large pages
##          : $memory_allocation    => Memory allocation to run Gatk
##          : $outfile_path         => Path to index
##          : $stderrfile_path      => Stderrfile path
##          : $temp_directory       => Redirect tmp files to java temp
##          : $verbosity            => Set the minimum level of logging

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $memory_allocation;
    my $outfile_path;
    my $stderrfile_path;
    my $temp_directory;

    ## Default(s)
    my $java_use_large_pages;
    my $verbosity;

    my $tmpl = {
        filehandle  => { store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        java_use_large_pages => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$java_use_large_pages,
            strict_type => 1,
        },
        memory_allocation => {
            store       => \$memory_allocation,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
        verbosity => {
            allow       => [qw{ INFO ERROR FATAL }],
            default     => q{INFO},
            store       => \$verbosity,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ gatk };

    ## Add java options
    gatk_java_options(
        {
            commands_ref         => \@commands,
            java_use_large_pages => $java_use_large_pages,
            memory_allocation    => $memory_allocation,
        }
    );

    push @commands, q{IndexFeatureFile};

    push @commands, q{--feature-file} . $SPACE . $infile_path;

    ## Add common options
    gatk_common_options(
        {
            commands_ref   => \@commands,
            temp_directory => $temp_directory,
            verbosity      => $verbosity,
        }
    );

    if ($outfile_path) {

        push @commands, q{--output} . $SPACE . $outfile_path;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gatk_java_options {

## Function : Push java commands for GATK to $commands_ref. Based on GATK 4.0
## Returns  : $commands_ref
## Arguments: $commands_ref         => List of commands added earlier {REF}
##          : $java_use_large_pages => Use java large pages
##          : $memory_allocation    => Memmory allocation for java
##          : $xargs_mode           => Escape quotation marks when running in xargs mode

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $commands_ref;
    my $java_use_large_pages;
    my $memory_allocation;
    my $xargs_mode;

    my $tmpl = {
        commands_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$commands_ref,
            strict_type => 1,
        },
        java_use_large_pages => {
            allow       => [ undef, 0, 1 ],
            store       => \$java_use_large_pages,
            strict_type => 1,
        },
        memory_allocation => {
            allow       => [ undef, qr/ ^Xm[sx]\d+[MG]$ /xmsi ],
            store       => \$memory_allocation,
            strict_type => 1,
        },
        xargs_mode => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$xargs_mode,
            strict_type => 1
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Return command array unchanged if no java options were given
    return $commands_ref
      if not $java_use_large_pages and not $memory_allocation;

    my @java_commands;

    # UseLargePages for requiring large memory pages (cross-platform flag)
    if ($java_use_large_pages) {

        push @java_commands, q{-XX:-UseLargePages};
    }

    if ($memory_allocation) {

        push @java_commands, q{-} . $memory_allocation;
    }

    my $java_command;

    ## Escape quotation marks if the program will be executed via xargs
    if ($xargs_mode) {
        $java_command =
            q{--java-options}
          . $SPACE
          . $ESCAPE
          . $DOUBLE_QUOTE
          . join( $SPACE, @java_commands )
          . $ESCAPE
          . $DOUBLE_QUOTE;
    }
    else {
        $java_command =
            q{--java-options}
          . $SPACE
          . $DOUBLE_QUOTE
          . join( $SPACE, @java_commands )
          . $DOUBLE_QUOTE;
    }

    push @{$commands_ref}, $java_command;

    return $commands_ref;
}

sub gatk_leftalignandtrimvariants {

## Function : Perl wrapper for writing GATK LeftAlignAndTrimVariants recipe to $filehandle. Based on GATK 4.0.0.
## Returns  : @commands
##          : $filehandle           => Sbatch filehandle to write to
##          : $infile_path          => Infile path
##          : $intervals_ref        => One or more genomic intervals over which to operate {REF}
##          : $java_use_large_pages => Use java large pages
##          : $memory_allocation    => Memory allocation to run Gatk
##          : $outfile_path         => Outfile path
##          : $referencefile_path   => Reference sequence file
##          : $split_multiallelics  => Split multiallelic records and left-align individual alleles
##          : $stderrfile_path      => Stderrfile path
##          : $temp_directory       => Redirect tmp files to java temp
##          : $verbosity            => Set the minimum level of logging
##          : $xargs_mode           => Set if the program will be executed via xargs

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $intervals_ref;
    my $memory_allocation;
    my $outfile_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $temp_directory;

    ## Default(s)
    my $java_use_large_pages;
    my $split_multiallelics;
    my $verbosity;
    my $xargs_mode;

    my $tmpl = {
        filehandle  => { store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        java_use_large_pages => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$java_use_large_pages,
            strict_type => 1,
        },
        memory_allocation => {
            store       => \$memory_allocation,
            strict_type => 1,
        },
        temp_directory => {
            strict_type => 1,
            store       => \$temp_directory,
        },
        intervals_ref => {
            default     => [],
            store       => \$intervals_ref,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
            required    => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        split_multiallelics => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$split_multiallelics,
            strict_type => 1,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        verbosity => {
            allow       => [qw{ INFO ERROR FATAL }],
            default     => q{INFO},
            store       => \$verbosity,
            strict_type => 1,
        },
        xargs_mode => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$xargs_mode,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ gatk };

    ## Add java options
    gatk_java_options(
        {
            commands_ref         => \@commands,
            java_use_large_pages => $java_use_large_pages,
            memory_allocation    => $memory_allocation,
            xargs_mode           => $xargs_mode,
        }
    );

    push @commands, q{LeftAlignAndTrimVariants};

    push @commands, q{--variant} . $SPACE . $infile_path;

    ## Add common options
    gatk_common_options(
        {
            commands_ref       => \@commands,
            intervals_ref      => $intervals_ref,
            referencefile_path => $referencefile_path,
            temp_directory     => $temp_directory,
            verbosity          => $verbosity,
        }
    );

    if ($split_multiallelics) {

        push @commands, q{--split-multi-allelics};
    }

    if ($outfile_path) {

        push @commands, q{--output} . $SPACE . $outfile_path;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gatk_printreads {

## Function : Perl wrapper for writing GATK PrintReads recipe to $filehandle. Based on GATK 4.0.8.
## Returns  : @commands
##          : $filehandle           => Sbatch filehandle to write to
##          : $infile_path          => Infile paths
##          : $intervals_ref        => One or more genomic intervals over which to operate {REF}
##          : $java_use_large_pages => Use java large pages
##          : $memory_allocation    => Memory allocation to run Gatk
##          : $outfile_path         => Outfile path
##          : $read_filters_ref     => Filters to apply to reads before analysis {REF}
##          : $referencefile_path   => Reference sequence file
##          : $stderrfile_path      => Stderrfile path
##          : $temp_directory       => Redirect tmp files to java temp
##          : $verbosity            => Set the minimum level of logging

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $intervals_ref;
    my $memory_allocation;
    my $outfile_path;
    my $read_filters_ref;
    my $referencefile_path;
    my $stderrfile_path;
    my $temp_directory;
    my $xargs_mode;

    ## Default(s)
    my $java_use_large_pages;
    my $verbosity;

    my $tmpl = {
        infile_path => {
            allow       => qr/ (?: bam | sam | cram )$ /xms,
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        filehandle    => { store => \$filehandle, },
        intervals_ref => {
            default     => [],
            defined     => 1,
            store       => \$intervals_ref,
            strict_type => 1,
        },
        java_use_large_pages => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$java_use_large_pages,
            strict_type => 1,
        },
        memory_allocation => {
            store       => \$memory_allocation,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        read_filters_ref => {
            default     => [],
            store       => \$read_filters_ref,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        temp_directory  => { store => \$temp_directory,  strict_type => 1, },
        verbosity       => {
            allow       => [qw{ INFO ERROR FATAL }],
            default     => q{INFO},
            store       => \$verbosity,
            strict_type => 1,
        },
        xargs_mode => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$xargs_mode,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ gatk };

    ## Add java options
    gatk_java_options(
        {
            commands_ref         => \@commands,
            java_use_large_pages => $java_use_large_pages,
            memory_allocation    => $memory_allocation,
            xargs_mode           => $xargs_mode,
        }
    );

    ## Add tool command
    push @commands, q{PrintReads};

    ## Add infile
    push @commands, q{--input} . $SPACE . $infile_path;

    ## Add common options
    gatk_common_options(
        {
            commands_ref       => \@commands,
            intervals_ref      => $intervals_ref,
            read_filters_ref   => $read_filters_ref,
            referencefile_path => $referencefile_path,
            temp_directory     => $temp_directory,
            verbosity          => $verbosity,
        }
    );

    ## Add Output
    push @commands, q{--output} . $SPACE . $outfile_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gatk_selectvariants {

## Function : Perl wrapper for writing GATK SelectVariants recipe to $filehandle. Based on GATK 4.0.8
## Returns  : @commands
## Arguments: $exclude_non_variants       => Exclude non-variant sites
##          : $filehandle                 => Sbatch filehandle to write to
##          : $infile_path                => Infile paths
##          : $intervals_ref              => One or more genomic intervals over which to operate {REF}
##          : $java_use_large_pages       => Use java large pages
##          : $memory_allocation          => Memory allocation to run Gatk
##          : $outfile_path               => Outfile path
##          : $pedigree                   => Pedigree files for samples
##          : $referencefile_path         => Reference sequence file
##          : $restrict_alleles_to        => Restrict allele output
##          : $sample_names_ref           => Include genotypes from this sample {REF}
##          : $select_type_to_include_ref => Select which variants to include {REF}
##          : $stderrfile_path            => Stderrfile path
##          : $temp_directory             => Redirect tmp files to java temp
##          : $verbosity                  => Set the minimum level of logging
##          : $xargs_mode                 => Set if the program will be executed via xargs

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $intervals_ref;
    my $memory_allocation;
    my $outfile_path;
    my $pedigree;
    my $referencefile_path;
    my $restrict_alleles_to;
    my $sample_names_ref;
    my $select_type_to_include_ref;
    my $stderrfile_path;
    my $temp_directory;
    my $xargs_mode;

    ## Default(s)
    my $exclude_non_variants;
    my $java_use_large_pages;
    my $verbosity;

    my $tmpl = {
        exclude_non_variants => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$exclude_non_variants,
            strict_type => 1,
        },
        filehandle  => { store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        intervals_ref => {
            default     => [],
            store       => \$intervals_ref,
            strict_type => 1,
        },
        java_use_large_pages => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$java_use_large_pages,
            strict_type => 1,
        },
        memory_allocation => {
            store       => \$memory_allocation,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        pedigree => {
            store       => \$pedigree,
            strict_type => 1,
        },
        restrict_alleles_to => {
            allow       => [qw{ ALL BIALLELIC MULTIALLELIC }],
            store       => \$restrict_alleles_to,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        sample_names_ref => {
            default     => [],
            defined     => 1,
            store       => \$sample_names_ref,
            strict_type => 1,
        },
        select_type_to_include_ref => {
            default     => [],
            defined     => 1,
            store       => \$select_type_to_include_ref,
            strict_type => 1,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
        verbosity => {
            allow       => [qw{ INFO ERROR FATAL }],
            default     => q{INFO},
            store       => \$verbosity,
            strict_type => 1,
        },
        xargs_mode => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$xargs_mode,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ gatk };

    ## Add java options
    gatk_java_options(
        {
            commands_ref         => \@commands,
            java_use_large_pages => $java_use_large_pages,
            memory_allocation    => $memory_allocation,
            xargs_mode           => $xargs_mode,
        }
    );

    ## Add tool command
    push @commands, q{SelectVariants};

    push @commands, q{--variant} . $SPACE . $infile_path;

    ## Add common options
    gatk_common_options(
        {
            commands_ref       => \@commands,
            intervals_ref      => $intervals_ref,
            pedigree           => $pedigree,
            referencefile_path => $referencefile_path,
            temp_directory     => $temp_directory,
            verbosity          => $verbosity,
        }
    );

    if ($exclude_non_variants) {

        push @commands, q{--exclude-non-variants};
    }

    if ($restrict_alleles_to) {

        push @commands, q{--restrict-alleles-to} . $SPACE . $restrict_alleles_to;
    }

    if ( @{$sample_names_ref} ) {

        push @commands,
          q{--sample-name} . $SPACE . join $SPACE . q{--sample-name} . $SPACE,
          @{$sample_names_ref};
    }

    if ( @{$select_type_to_include_ref} ) {

        push @commands,
            q{--select-type-to-include}
          . $SPACE
          . join $SPACE
          . q{--select-type-to-include}
          . $SPACE, @{$select_type_to_include_ref};
    }

    push @commands, q{--output} . $SPACE . $outfile_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gatk_splitncigarreads {

## Function : Perl wrapper for writing GATK splitNCigarReads recipe to $filehandle. Based on GATK 4.0.8.
## Returns  : @commands
## Arguments: $filehandle           => Sbatch filehandle to write to
##          : $infile_path          => Infile paths
##          : $java_use_large_pages => Use java large pages
##          : $maping_quality_from  => Input mapping quality values
##          : $maping_quality_to    => Output mapping quality values
##          : $memory_allocation    => Memory allocation to run Gatk
##          : $outfile_path         => Outfile path
##          : $operation            => The read operation that is to be applied on the bam file
##          : $read_filters_ref     => Filters to apply to reads before analysis {REF}
##          : $referencefile_path   => Reference sequence file
##          : $stderrfile_path      => Stderrfile path
##          : $temp_directory       => Redirect tmp files to java temp
##          : $verbosity            => Set the minimum level of logging
##          : $xargs_mode           => Set if the program will be executed via xargs

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $memory_allocation;
    my $outfile_path;
    my $read_filters_ref;
    my $referencefile_path;
    my $stderrfile_path;
    my $temp_directory;

    ## Default(s)
    my $java_use_large_pages;
    my $verbosity;
    my $xargs_mode;

    my $tmpl = {
        filehandle  => { store => \$filehandle, },
        infile_path => {
            allow       => qr/ (?: bam | sam | cram )$ /xms,
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        java_use_large_pages => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$java_use_large_pages,
            strict_type => 1,
        },
        memory_allocation => {
            store       => \$memory_allocation,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
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
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
        verbosity => {
            allow       => [qw{ INFO ERROR FATAL }],
            default     => q{INFO},
            store       => \$verbosity,
            strict_type => 1,
        },
        xargs_mode => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$xargs_mode,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ gatk };

    ## Add java options
    gatk_java_options(
        {
            commands_ref         => \@commands,
            java_use_large_pages => $java_use_large_pages,
            memory_allocation    => $memory_allocation,
            xargs_mode           => $xargs_mode,
        }
    );

    ## Add tool command
    push @commands, q{SplitNCigarReads};

    ## Add input file
    push @commands, q{--input} . $SPACE . $infile_path;

    ## Add common options
    gatk_common_options(
        {
            commands_ref       => \@commands,
            read_filters_ref   => $read_filters_ref,
            referencefile_path => $referencefile_path,
            temp_directory     => $temp_directory,
            verbosity          => $verbosity,
        }
    );

    ## Add output file
    push @commands, q{--output} . $SPACE . $outfile_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gatk_varianteval {

## Function : Perl wrapper for writing GATK varianteval recipe to $filehandle. Based on GATK 4.1.0.
## Returns  : @commands
## Arguments: $dbsnp_file_path               => DbSNP file path
##          : $filehandle                    => Sbatch filehandle to write to
##          : $indel_gold_standard_file_path => Evaluations that count calls at sites of true variation (e.g., indel calls) will use this argument as their gold standard for comparison
##          : $infile_paths_ref              => Infile paths
##          : $intervals_ref                 => One or more genomic intervals over which to operate {REF}
##          : $java_use_large_pages          => Use java large pages
##          : $memory_allocation             => Memory allocation to run Gatk
##          : $outfile_path                  => Outfile path
##          : $pedigree                      => Pedigree files for samples
##          : $referencefile_path            => Reference sequence file
##          : $stderrfile_path               => Stderrfile path
##          : $temp_directory                => Redirect tmp files to java temp
##          : $verbosity                     => Set the minimum level of logging
##          : $xargs_mode                    => Set if the program will be executed via xargs

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $dbsnp_file_path;
    my $filehandle;
    my $indel_gold_standard_file_path;
    my $infile_paths_ref;
    my $intervals_ref;
    my $java_use_large_pages;
    my $memory_allocation;
    my $outfile_path;
    my $pedigree;
    my $referencefile_path;
    my $stderrfile_path;
    my $temp_directory;

    ## Default(s)
    my $verbosity;
    my $xargs_mode;

    my $tmpl = {
        dbsnp_file_path => { store => \$dbsnp_file_path, strict_type => 1, },
        filehandle      => { store => \$filehandle, },
        indel_gold_standard_file_path =>
          { store => \$indel_gold_standard_file_path, strict_type => 1, },
        intervals_ref => { default => [], store => \$intervals_ref, strict_type => 1, },
        infile_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infile_paths_ref,
            strict_type => 1,
        },
        java_use_large_pages => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$java_use_large_pages,
            strict_type => 1,
        },
        memory_allocation => { store => \$memory_allocation, strict_type => 1, },
        outfile_path      => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        pedigree           => { store => \$pedigree, strict_type => 1, },
        referencefile_path => {
            defined     => 1,
            required    => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        temp_directory  => { store => \$temp_directory,  strict_type => 1, },
        verbosity       => {
            allow       => [qw{ INFO ERROR FATAL }],
            default     => q{INFO},
            store       => \$verbosity,
            strict_type => 1,
        },
        xargs_mode => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$xargs_mode,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ gatk };

    ## Add java options
    gatk_java_options(
        {
            commands_ref         => \@commands,
            java_use_large_pages => $java_use_large_pages,
            memory_allocation    => $memory_allocation,
            xargs_mode           => $xargs_mode,
        }
    );

    push @commands, q{VariantEval};

    ## Add common options
    gatk_common_options(
        {
            commands_ref       => \@commands,
            intervals_ref      => $intervals_ref,
            pedigree           => $pedigree,
            referencefile_path => $referencefile_path,
            temp_directory     => $temp_directory,
            verbosity          => $verbosity,
        }
    );

    if ($dbsnp_file_path) {

        push @commands, q{--dbsnp} . $SPACE . $dbsnp_file_path;
    }

    if ($indel_gold_standard_file_path) {

        push @commands, q{--gold-standard} . $SPACE . $indel_gold_standard_file_path;
    }

    if ( @{$infile_paths_ref} ) {

        push @commands, q{--eval} . $SPACE . join $SPACE . q{--eval} . $SPACE,
          @{$infile_paths_ref};
    }

    if ($outfile_path) {
        push @commands, q{--output} . $SPACE . $outfile_path;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gatk_variantfiltration {

## Function : Perl wrapper for writing GATK VariantFiltration recipe to $filehandle. Based on GATK 4.1.0.
## Returns  : @commands
##          : $cluster_size         => Number of SNPs which make up a cluster
##          : $cluster_window_size  => Window size (in bases) in which to evaluate clustered SNPs
##          : $filehandle           => Sbatch filehandle to write to
##          : $filter_href          => Hash with the name of the filter as key and the filter expression as value {REF}
##          : $infile_path          => Infile paths
##          : $intervals_ref        => One or more genomic intervals over which to operate {REF}
##          : $java_use_large_pages => Use java large pages
##          : $memory_allocation    => Memory allocation to run Gatk
##          : $outfile_path         => Outfile path
##          : $read_filters_ref     => Filters to apply on reads {REF}
##          : $referencefile_path   => Reference sequence file
##          : $stderrfile_path      => Stderrfile path
##          : $temp_directory       => Redirect tmp files to java temp
##          : $verbosity            => Set the minimum level of logging
##          : $xargs_mode           => Set if the program will be executed via xargs

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $cluster_size;
    my $cluster_window_size;
    my $filehandle;
    my $filter_href;
    my $infile_path;
    my $intervals_ref;
    my $memory_allocation;
    my $outfile_path;
    my $read_filters_ref;
    my $referencefile_path;
    my $stderrfile_path;
    my $temp_directory;

    ## Default(s)
    my $java_use_large_pages;
    my $verbosity;
    my $xargs_mode;

    my $tmpl = {
        cluster_size => {
            allow       => qr/ \A \d+ \z /sxm,
            store       => \$cluster_size,
            strict_type => 1,
        },
        cluster_window_size => {
            allow       => qr/ \A \d+ \z /xms,
            store       => \$cluster_window_size,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        filter_href => {
            default     => {},
            store       => \$filter_href,
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
            store       => \$intervals_ref,
            strict_type => 1,
        },
        java_use_large_pages => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$java_use_large_pages,
            strict_type => 1,
        },
        memory_allocation => {
            store       => \$memory_allocation,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        read_filters_ref => {
            default     => [],
            store       => \$read_filters_ref,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
        verbosity => {
            allow       => [qw{INFO ERROR FATAL }],
            default     => q{INFO},
            store       => \$verbosity,
            strict_type => 1,
        },
        xargs_mode => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$xargs_mode,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ gatk };

    ## Add java options
    gatk_java_options(
        {
            commands_ref         => \@commands,
            java_use_large_pages => $java_use_large_pages,
            memory_allocation    => $memory_allocation,
            xargs_mode           => $xargs_mode,
        }
    );

    push @commands, q{VariantFiltration};

    push @commands, q{--variant} . $SPACE . $infile_path;

    ## Add common options
    gatk_common_options(
        {
            commands_ref       => \@commands,
            intervals_ref      => $intervals_ref,
            read_filters_ref   => $read_filters_ref,
            referencefile_path => $referencefile_path,
            temp_directory     => $temp_directory,
            verbosity          => $verbosity,
        }
    );

    if ($cluster_size) {

        push @commands, q{--cluster-size} . $SPACE . $cluster_size;
    }

    if ($cluster_window_size) {

        push @commands, q{--cluster-window-size} . $SPACE . $cluster_window_size;
    }

    if ($filter_href) {

      FILTERNAME:
        foreach my $filtername ( keys %{$filter_href} ) {
            push @commands,
                q{--filter-name}
              . $SPACE
              . $filtername
              . $SPACE
              . q{--filter-expression}
              . $SPACE
              . $DOUBLE_QUOTE
              . $filter_href->{$filtername}
              . $DOUBLE_QUOTE;
        }
    }

    push @commands, q{--output} . $SPACE . $outfile_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gatk_variantrecalibrator {

## Function : Perl wrapper for writing GATK variantrecalibrator recipe to $filehandle. Based on GATK 4.1.0.
## Returns  : @commands
##          : $annotations_ref       => One or more specific annotations to apply to variant calls
##          : $filehandle            => Sbatch filehandle to write to
##          : $infile_path           => Infile path
##          : $intervals_ref         => One or more genomic intervals over which to operate {REF}
##          : $java_use_large_pages  => Use java large pages
##          : $max_attempts          => Number of attempts to build the model
##          : $max_gaussian_level    => Max number of Gaussians for the positive model
##          : $memory_allocation     => Memory allocation to run Gatk
##          : $mode                  => Mode for emitting reference confidence scores
##          : $outfile_path          => The output recal file used by ApplyRecalibration
##          : $read_filters_ref      => Filters to apply to reads before analysis {REF}
##          : $referencefile_path    => Reference sequence file
##          : $resources_ref         => A list of sites for which to apply a prior probability of being correct but which aren't used by the algorithm (training and truth setsare required to run)
##          : $rscript_file_path     => Rscript file path
##          : $stderrfile_path       => Stderrfile path
##          : $temp_directory        => Redirect tmp files to java temp
##          : $tranches_file_path    => The output tranches file used by ApplyRecalibration
##          : $trust_all_polymorphic => Trust that all the input training sets' unfiltered records contain only polymorphic sites to speed up computation
##          : $ts_tranches_ref       => Levels of truth sensitivity at which to slice the data. (in percent, that is 1.0 for 1 percent)
##          : $verbosity	           => Set the minimum level of logging
##          : $xargs_mode            => Set if the program will be executed via xargs

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $annotations_ref;
    my $filehandle;
    my $infile_path;
    my $intervals_ref;
    my $max_gaussian_level;
    my $memory_allocation;
    my $mode;
    my $outfile_path;
    my $read_filters_ref;
    my $referencefile_path;
    my $resources_ref;
    my $rscript_file_path;
    my $stderrfile_path;
    my $temp_directory;
    my $tranches_file_path;
    my $trust_all_polymorphic;
    my $ts_tranches_ref;

    ## Default(s)
    my $java_use_large_pages;
    my $max_attempts;
    my $verbosity;
    my $xargs_mode;

    my $tmpl = {
        annotations_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$annotations_ref,
            strict_type => 1,
        },
        filehandle  => { store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        intervals_ref => {
            default     => [],
            store       => \$intervals_ref,
            strict_type => 1,
        },
        java_use_large_pages => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$java_use_large_pages,
            strict_type => 1,
        },
        max_attempts => {
            allow       => [qr/ \A \d+ \z /sxm],
            default     => 3,
            store       => \$max_attempts,
            strict_type => 1,
        },
        max_gaussian_level => {
            allow       => [ undef, qr/ ^\d+$ /sxm ],
            store       => \$max_gaussian_level,
            strict_type => 1,
        },
        memory_allocation => {
            store       => \$memory_allocation,
            strict_type => 1,
        },
        mode => {
            allow       => [ undef, qw{ SNP INDEL BOTH } ],
            defined     => 1,
            store       => \$mode,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        read_filters_ref => {
            default     => [],
            store       => \$read_filters_ref,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        resources_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$resources_ref,
            strict_type => 1,
        },
        rscript_file_path => {
            store       => \$rscript_file_path,
            strict_type => 1,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
        tranches_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$tranches_file_path,
            strict_type => 1,
        },
        trust_all_polymorphic => {
            allow       => [ undef, 0, 1 ],
            store       => \$trust_all_polymorphic,
            strict_type => 1,
        },
        ts_tranches_ref => {
            default     => [],
            store       => \$ts_tranches_ref,
            strict_type => 1,
        },
        verbosity => {
            allow       => [qw{ INFO ERROR FATAL }],
            default     => q{INFO},
            store       => \$verbosity,
            strict_type => 1,
        },
        xargs_mode => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$xargs_mode,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ gatk };

    ## Add java options
    gatk_java_options(
        {
            commands_ref         => \@commands,
            java_use_large_pages => $java_use_large_pages,
            memory_allocation    => $memory_allocation,
            xargs_mode           => $xargs_mode,
        }
    );

    ## Add tool command
    push @commands, q{VariantRecalibrator};

    push @commands, q{--variant} . $SPACE . $infile_path;

    ## Add common options
    gatk_common_options(
        {
            commands_ref       => \@commands,
            intervals_ref      => $intervals_ref,
            read_filters_ref   => $read_filters_ref,
            referencefile_path => $referencefile_path,
            temp_directory     => $temp_directory,
            verbosity          => $verbosity,
        }
    );

    push @commands,
      q{--use-annotation} . $SPACE . join $SPACE . q{--use-annotation} . $SPACE,
      @{$annotations_ref};

    if ($max_attempts) {

        push @commands, q{--max-attempts} . $SPACE . $max_attempts;
    }

    if ($max_gaussian_level) {

        push @commands, q{--max-gaussians} . $SPACE . $max_gaussian_level;
    }

    if ($mode) {

        push @commands, q{--mode} . $SPACE . $mode;
    }

    if ($rscript_file_path) {

        push @commands, q{--rscript-file} . $SPACE . $rscript_file_path;
    }

    push @commands,
      q{--resource} . $COLON . join $SPACE . q{--resource} . $COLON,
      @{$resources_ref};

    if ($ts_tranches_ref) {

        push @commands, q{-tranche} . $SPACE . join $SPACE . q{-tranche} . $SPACE,
          @{$ts_tranches_ref};
    }

    push @commands, q{--tranches-file} . $SPACE . $tranches_file_path;

    if ($trust_all_polymorphic) {

        push @commands, q{--trust-all-polymorphic};
    }

    push @commands, q{--output} . $SPACE . $outfile_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

1;
