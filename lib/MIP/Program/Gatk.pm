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
use MIP::Constants qw{ $SPACE $DOUBLE_QUOTE };
use MIP::Program::Base::Gatk qw{ gatk_common_options gatk_java_options };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.13;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      gatk_applybqsr
      gatk_asereadcounter
      gatk_baserecalibrator
      gatk_gatherbqsrreports
      gatk_haplotypecaller
      gatk_printreads
      gatk_splitncigarreads
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

sub gatk_haplotypecaller {

## Function : Perl wrapper for writing GATK haplotypecaller recipe to $filehandle. Based on GATK 4.1.0.
## Returns  : @commands
## Arguments: $annotations_ref                               => One or more specific annotations to apply to variant calls
##          : $dbsnp_path                                    => Path to DbSNP file
##          : $dont_use_soft_clipped_bases                   => Do not analyze soft clipped bases in the reads
##          : $emit_ref_confidence                           => Mode for emitting reference confidence scores
##          : $filehandle                                    => Sbatch filehandle to write to
##          : $infile_path                                   => Infile paths
##          : $intervals_ref                                 => One or more genomic intervals over which to operate {REF}
##          : $java_use_large_pages                          => Use java large pages
##          : $memory_allocation                             => Memory allocation to run Gatk
##          : $num_ref_samples_if_no_call                    => Number of hom-ref genotypes to infer at sites not present in a panel
##          : $outfile_path                                  => Outfile path
##          : $pcr_indel_model                               => The PCR indel model to use
##          : $pedigree                                      => Pedigree files for samples
##          : $population_callset                            => Callset to use in calculating genotype priors
##          : $read_filters_ref                              => Filters to apply to reads before analysis {REF}
##          : $referencefile_path                            => Reference sequence file
##          : $sample_ploidy                                 => Ploidy per sample
##          : $standard_min_confidence_threshold_for_calling => The minimum phred-scaled confidence threshold at which variants should be called
##          : $stderrfile_path                               => Stderrfile path
##          : $temp_directory                                => Redirect tmp files to java temp
##          : $use_new_qual_calculator                       => Use the new AF model instead of the so-called exact model
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
    my $num_ref_samples_if_no_call;
    my $outfile_path;
    my $pcr_indel_model;
    my $pedigree;
    my $population_callset;
    my $read_filters_ref;
    my $referencefile_path;
    my $sample_ploidy;
    my $standard_min_confidence_threshold_for_calling;
    my $stderrfile_path;
    my $temp_directory;

    ## Default(s)
    my $emit_ref_confidence;
    my $java_use_large_pages;
    my $use_new_qual_calculator;
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
        memory_allocation => {
            store       => \$memory_allocation,
            strict_type => 1,
        },
        num_ref_samples_if_no_call => {
            allow       => [ undef, qr/ ^\d+$ /sxm ],
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
        pcr_indel_model => {
            allow       => [ undef, qw{ NONE HOSTILE AGGRESSIVE CONSERVATIVE } ],
            store       => \$pcr_indel_model,
            strict_type => 1,
        },
        population_callset => {
            store       => \$population_callset,
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
        use_new_qual_calculator => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$use_new_qual_calculator,
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

    if ($population_callset) {

        push @commands, q{--population-callset} . $SPACE . $population_callset;
    }

    if ($num_ref_samples_if_no_call) {
        push @commands,
          q{--num-reference-samples-if-no-call} . $SPACE . $num_ref_samples_if_no_call;
    }

    if ($use_new_qual_calculator) {

        push @commands, q{--use-new-qual-calculator};
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

1;
