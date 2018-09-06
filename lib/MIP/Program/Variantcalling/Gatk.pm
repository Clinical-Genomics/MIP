package MIP::Program::Variantcalling::Gatk;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use FindBin qw{ $Bin };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };

## CPANM
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Language::Java qw{ java_core };
use MIP::Program::Base::Gatk
  qw{ gatk_base gatk_common_options gatk_java_options };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.06;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      gatk_applyrecalibration
      gatk_asereadcounter
      gatk_calculategenotypeposteriors
      gatk_catvariants gatk_combinevariants
      gatk_concatenate_variants
      gatk_genomicsdbimport
      gatk_genotypegvcfs
      gatk_leftalignandtrimvariants
      gatk_selectvariants
      gatk_varianteval
      gatk_variantfiltration
      gatk_variantrecalibrator
    };

}

## Constants
Readonly my $AMPERSAND    => q{&};
Readonly my $DOUBLE_QOUTE => q{"};
Readonly my $EMPTY_STR    => q{};
Readonly my $NEWLINE      => qq{\n};
Readonly my $SPACE        => q{ };

sub gatk_genotypegvcfs {

## Function : Perl wrapper for writing GATK genotypegvcfs recipe to $FILEHANDLE. Based on GATK 3.7.0.
## Returns  : @commands
## Arguments: $memory_allocation                     => Memory allocation to run Gatk
##          : $java_use_large_pages                  => Use java large pages
##          : $temp_directory                        => Redirect tmp files to java temp
##          : $java_jar                              => Java jar
##          : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##          : $infile_paths_ref                      => Infile paths {REF}
##          : $outfile_path                          => Outfile path
##          : $referencefile_path                    => Reference sequence file
##          : $stderrfile_path                       => Stderrfile path
##          : $FILEHANDLE                            => Sbatch filehandle to write to
##          : $pedigree                              => Pedigree files for samples
##          : $dbsnp_file_path                       => DbSNP file
##          : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##          : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##          : $logging_level                         => Set the minimum level of logging
##          : $pedigree_validation_type              => Validation strictness for pedigree
##          : $include_nonvariant_sites              => Include loci found to be non-variant after genotyping

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $memory_allocation;
    my $java_use_large_pages;
    my $temp_directory;
    my $java_jar;
    my $intervals_ref;
    my $infile_paths_ref;
    my $outfile_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $pedigree;
    my $dbsnp_file_path;
    my $downsample_to_coverage;

    ## Default(s)
    my $gatk_disable_auto_index_and_file_lock;
    my $logging_level;
    my $pedigree_validation_type;
    my $include_nonvariant_sites;

    my $tmpl = {
        memory_allocation => { strict_type => 1, store => \$memory_allocation },
        java_use_large_pages => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$java_use_large_pages
        },
        temp_directory => { strict_type => 1, store => \$temp_directory },
        java_jar       => { strict_type => 1, store => \$java_jar },
        intervals_ref =>
          { default => [], strict_type => 1, store => \$intervals_ref },
        infile_paths_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$infile_paths_ref
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
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        FILEHANDLE      => { store       => \$FILEHANDLE },
        pedigree        => { strict_type => 1, store => \$pedigree },
        dbsnp_file_path => { strict_type => 1, store => \$dbsnp_file_path },
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
        pedigree_validation_type => {
            default     => q{SILENT},
            allow       => [qw{ STRICT SILENT }],
            strict_type => 1,
            store       => \$pedigree_validation_type
        },
        include_nonvariant_sites => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$include_nonvariant_sites
        },

    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands;

    if ($java_jar) {    # Write java core commands to filehandle.
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
            analysis_type            => q{GenotypeGVCFs},
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

    # Tool-specific options

    if ($include_nonvariant_sites) {

        push @commands, q{--includeNonVariantSites};
    }

    if ($dbsnp_file_path) {

        push @commands, q{--dbsnp} . $SPACE . $dbsnp_file_path;
    }

    # Infile
    if ( @{$infile_paths_ref} ) {

        push @commands,
          q{--variant} . $SPACE . join $SPACE . q{--variant} . $SPACE,
          @{$infile_paths_ref};
    }

    ## Output
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;

}

sub gatk_selectvariants {

## Function : Perl wrapper for writing GATK selectvariants recipe to $FILEHANDLE. Based on GATK 3.7.0.
## Returns  : @commands
## Arguments: $memory_allocation                     => Memory allocation to run Gatk
##          : $java_use_large_pages                  => Use java large pages
##          : $temp_directory                        => Redirect tmp files to java temp
##          : $java_jar                              => Java jar
##          : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##          : $sample_names_ref                      => Include genotypes from this sample {REF}
##          : $infile_path                           => Infile paths
##          : $outfile_path                          => Outfile path
##          : $referencefile_path                    => Reference sequence file
##          : $stderrfile_path                       => Stderrfile path
##          : $FILEHANDLE                            => Sbatch filehandle to write to
##          : $pedigree                              => Pedigree files for samples
##          : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##          : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##          : $logging_level                         => Set the minimum level of logging
##          : $pedigree_validation_type              => Validation strictness for pedigree
##          : $exclude_nonvariants                   => Exclude non-variant sites

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $memory_allocation;
    my $java_use_large_pages;
    my $temp_directory;
    my $java_jar;
    my $intervals_ref;
    my $sample_names_ref;
    my $infile_path;
    my $outfile_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $pedigree;
    my $downsample_to_coverage;

    ## Default(s)
    my $gatk_disable_auto_index_and_file_lock;
    my $logging_level;
    my $pedigree_validation_type;
    my $exclude_nonvariants;

    my $tmpl = {
        memory_allocation => { strict_type => 1, store => \$memory_allocation },
        java_use_large_pages => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$java_use_large_pages
        },
        temp_directory => { strict_type => 1, store => \$temp_directory },
        java_jar       => { strict_type => 1, store => \$java_jar },
        intervals_ref =>
          { default => [], strict_type => 1, store => \$intervals_ref },
        sample_names_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$sample_names_ref
        },
        infile_path => {
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
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        FILEHANDLE => { store       => \$FILEHANDLE },
        pedigree   => { strict_type => 1, store => \$pedigree },
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
        pedigree_validation_type => {
            default     => q{SILENT},
            allow       => [qw{ STRICT SILENT }],
            strict_type => 1,
            store       => \$pedigree_validation_type
        },
        exclude_nonvariants => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$exclude_nonvariants
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands;

    if ($java_jar) {    # Write java core commands to filehandle.
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
            analysis_type            => q{SelectVariants},
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

    ## Tool-specific options
    if ($exclude_nonvariants) {

        push @commands, q{--excludeNonVariants};
    }

    if ( @{$sample_names_ref} ) {

        push @commands,
          q{--sample_name} . $SPACE . join $SPACE . q{--sample_name} . $SPACE,
          @{$sample_names_ref};
    }

    ## Infile
    if ($infile_path) {

        push @commands, q{--variant} . $SPACE . $infile_path;
    }

    ## Output
    if ($outfile_path) {

        push @commands, q{--out} . $SPACE . $outfile_path;
    }

    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
        }
      );

    return @commands;

}

sub gatk_catvariants {

## Function : Perl wrapper for writing GATK catvariants recipe to $FILEHANDLE. Based on GATK 3.7.0.
## Returns  : @commands
## Arguments: $memory_allocation                     => Memory allocation to run Gatk
##          : $java_use_large_pages                  => Use java large pages
##          : $temp_directory                        => Redirect tmp files to java temp
##          : $java_jar                              => Java jar
##          : $gatk_path                             => Path to java jar and analysis to run
##          : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##          : $infile_paths_ref                      => Infile paths {REF}
##          : $outfile_path                          => Outfile path
##          : $referencefile_path                    => Reference sequence file
##          : $stderrfile_path                       => Stderrfile path
##          : $FILEHANDLE                            => Sbatch filehandle to write to
##          : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##          : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##          : $logging_level                         => Set the minimum level of logging
##          : $assume_sorted                         => Assume_sorted should be true if the input files are already sorted

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $memory_allocation;
    my $java_use_large_pages;
    my $temp_directory;
    my $gatk_path;
    my $intervals_ref;
    my $infile_paths_ref;
    my $outfile_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $downsample_to_coverage;

    ## Default(s)
    my $gatk_disable_auto_index_and_file_lock;
    my $logging_level;
    my $assume_sorted;

    my $tmpl = {
        memory_allocation => { strict_type => 1, store => \$memory_allocation },
        java_use_large_pages => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$java_use_large_pages
        },
        temp_directory => { strict_type => 1, store => \$temp_directory },
        gatk_path      => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$gatk_path
        },
        intervals_ref =>
          { default => [], strict_type => 1, store => \$intervals_ref },
        infile_paths_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$infile_paths_ref
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
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        FILEHANDLE             => { store => \$FILEHANDLE },
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
        assume_sorted => {
            default     => 1,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$assume_sorted
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands;

    ## Write java core commands to filehandle
    if ($memory_allocation) {

        @commands = java_core(
            {
                memory_allocation    => $memory_allocation,
                java_use_large_pages => $java_use_large_pages,
                temp_directory       => $temp_directory,
            }
        );
    }
    ## Base command
    push @commands, q{-cp};

    ## Gatk path
    push @commands, $gatk_path;

    push @commands, q{--logging_level} . $SPACE . $logging_level;

    if ( @{$intervals_ref} ) {

        push @commands,
          q{--intervals} . $SPACE . join $SPACE . q{--intervals} . $SPACE,
          @{$intervals_ref};
    }

    if ($referencefile_path) {

        push @commands, q{--reference} . $SPACE . $referencefile_path;
    }

    if ($downsample_to_coverage) {

        push @commands,
          q{--downsample_to_coverage} . $SPACE . $downsample_to_coverage;
    }

    if ($gatk_disable_auto_index_and_file_lock) {

        push @commands,
          q{--disable_auto_index_creation_and_locking_when_reading_rods};
    }

    ## Tool-specific options
    if ($assume_sorted) {

        push @commands, q{--assumeSorted};
    }

    ## Infile
    if ( @{$infile_paths_ref} ) {

        push @commands,
          q{--variant} . $SPACE . join $SPACE . q{--variant} . $SPACE,
          @{$infile_paths_ref};

    }

    ## Output
    if ($outfile_path) {

        push @commands, q{--outputFile} . $SPACE . $outfile_path;
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;

}

sub gatk_variantrecalibrator {

## Function : Perl wrapper for writing GATK variantrecalibrator recipe to $FILEHANDLE. Based on GATK 3.7.0.
## Returns  : @commands
## Arguments: $memory_allocation                     => Memory allocation to run Gatk
##          : $java_use_large_pages                  => Use java large pages
##          : $temp_directory                        => Redirect tmp files to java temp
##          : $java_jar                              => Java jar
##          : $infile_paths_ref                      => Infile paths
##          : $resources_ref                         => A list of sites for which to apply a prior probability of being correct but which aren't used by the algorithm (training and truth setsare required to run)
##          : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##          : $read_filters_ref                      => Filters to apply to reads before analysis {REF}
##          : $static_quantized_quals_ref            => Use static quantized quality scores to a given number of levels (with -BQSR) {REF}
##          : $annotations_ref                       => One or more specific annotations to apply to variant calls
##          : $referencefile_path                    => Reference sequence file
##          : $mode                                  => Mode for emitting reference confidence scores
##          : $base_quality_score_recalibration_file => Input covariates table file for on-the-fly base quality score recalibration
##          : $stderrfile_path                       => Stderrfile path
##          : $FILEHANDLE                            => Sbatch filehandle to write to
##          : $pedigree                              => Pedigree files for samples
##          : $rscript_file_path                     => Rscript file path
##          : $recal_file_path                       => The output recal file used by ApplyRecalibration
##          : $tranches_file_path                    => The output tranches file used by ApplyRecalibration
##          : $num_cpu_threads_per_data_thread       => Number of CPU threads to allocate per data thread
##          : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##          : $max_gaussian_level                    => Max number of Gaussians for the positive model
##          : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##          : $disable_indel_qual                    => Disable printing of base insertion and deletion tags (with -BQSR)
##          : $logging_level                         => Set the minimum level of logging
##          : $pedigree_validation_type              => Validation strictness for pedigree

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $memory_allocation;
    my $java_use_large_pages;
    my $temp_directory;
    my $java_jar;
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

    ## Default(s)
    my $gatk_disable_auto_index_and_file_lock;
    my $disable_indel_qual;
    my $logging_level;
    my $pedigree_validation_type;

    my $tmpl = {
        memory_allocation => { strict_type => 1, store => \$memory_allocation },
        java_use_large_pages => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$java_use_large_pages
        },
        temp_directory   => { strict_type => 1, store => \$temp_directory },
        java_jar         => { strict_type => 1, store => \$java_jar },
        infile_paths_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$infile_paths_ref
        },
        resources_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$resources_ref
        },
        intervals_ref =>
          { default => [], strict_type => 1, store => \$intervals_ref },
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
        referencefile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$referencefile_path
        },
        mode => {
            required    => 1,
            defined     => 1,
            allow       => [qw{ SNP INDEL BOTH }],
            strict_type => 1,
            store       => \$mode
        },
        base_quality_score_recalibration_file => {
            strict_type => 1,
            store       => \$base_quality_score_recalibration_file
        },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        FILEHANDLE        => { store       => \$FILEHANDLE },
        pedigree          => { strict_type => 1, store => \$pedigree },
        rscript_file_path => { strict_type => 1, store => \$rscript_file_path },
        recal_file_path   => { strict_type => 1, store => \$recal_file_path },
        tranches_file_path =>
          { strict_type => 1, store => \$tranches_file_path },
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
        max_gaussian_level => {
            allow       => [ undef, qr/ ^\d+$ /sxm ],
            strict_type => 1,
            store       => \$max_gaussian_level
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands;

    if ($java_jar) {    # Write java core commands to filehandle.
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
            analysis_type            => q{VariantRecalibrator},
            logging_level            => $logging_level,
            intervals_ref            => $intervals_ref,
            referencefile_path       => $referencefile_path,
            pedigree                 => $pedigree,
            pedigree_validation_type => $pedigree_validation_type,
            downsample_to_coverage   => $downsample_to_coverage,
            gatk_disable_auto_index_and_file_lock =>
              $gatk_disable_auto_index_and_file_lock,
            base_quality_score_recalibration_file =>
              $base_quality_score_recalibration_file,
            disable_indel_qual         => $disable_indel_qual,
            static_quantized_quals_ref => $static_quantized_quals_ref,
        }
    );

    if ($num_cpu_threads_per_data_thread) {

        push @commands,
            q{--num_cpu_threads_per_data_thread}
          . $SPACE
          . $num_cpu_threads_per_data_thread;
    }

    if ( @{$read_filters_ref} ) {

        push @commands,
          q{--read_filter} . $SPACE . join $SPACE . q{--read_filter} . $SPACE,
          @{$read_filters_ref};
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

    if ($rscript_file_path) {

        push @commands, q{--rscript_file} . $SPACE . $rscript_file_path;
    }

    if ($recal_file_path) {

        push @commands, q{--recal_file} . $SPACE . $recal_file_path;
    }

    if ($tranches_file_path) {

        push @commands, q{--tranches_file} . $SPACE . $tranches_file_path;
    }

    if ($max_gaussian_level) {

        push @commands, q{--maxGaussians} . $SPACE . $max_gaussian_level;
    }

    if ( @{$annotations_ref} ) {

        push @commands,
            q{--use_annotation}
          . $SPACE
          . join $SPACE
          . q{--use_annotation}
          . $SPACE, @{$annotations_ref};
    }

    if ( @{$resources_ref} ) {

        push @commands,
          q{--resource:} . join $SPACE . q{--resource:}, @{$resources_ref};
    }

    if ($mode) {

        push @commands, q{--mode} . $SPACE . $mode;
    }

    ## Infile
    if ( @{$infile_paths_ref} ) {

        push @commands, q{--input} . $SPACE . join q{--input_file} . $SPACE,
          @{$infile_paths_ref};
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;

}

sub gatk_applyrecalibration {

## Function : Perl wrapper for writing GATK applyrecalibration recipe to $FILEHANDLE. Based on GATK 3.7.0.
## Returns  : @commands
## Arguments: $memory_allocation                     => Memory allocation to run Gatk
##          : $java_use_large_pages                  => Use java large pages
##          : $temp_directory                        => Redirect tmp files to java temp
##          : $java_jar                              => Java jar
##          : $infile_path                           => Infile paths
##          : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##          : $read_filters_ref                      => Filters to apply to reads before analysis {REF}
##          : $static_quantized_quals_ref            => Use static quantized quality scores to a given number of levels (with -BQSR) {REF}
##          : $outfile_path                          => Outfile path
##          : $referencefile_path                    => Reference sequence file
##          : $mode                                  => Mode for emitting reference confidence scores
##          : $base_quality_score_recalibration_file => Input covariates table file for on-the-fly base quality score recalibration
##          : $stderrfile_path                       => Stderrfile path
##          : $FILEHANDLE                            => Sbatch filehandle to write to
##          : $pedigree                              => Pedigree files for samples
##          : $recal_file_path                       => The output recal file used by ApplyRecalibration
##          : $tranches_file_path                    => The output tranches file used by ApplyRecalibration
##          : $num_cpu_threads_per_data_thread       => Number of CPU threads to allocate per data thread
##          : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##          : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##          : $disable_indel_qual                    => Disable printing of base insertion and deletion tags (with -BQSR)
##          : $logging_level                         => Set the minimum level of logging
##          : $pedigree_validation_type              => Validation strictness for pedigree
##          : $ts_filter_level                       => Ts filter level

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $memory_allocation;
    my $java_use_large_pages;
    my $temp_directory;
    my $java_jar;
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

    ## Default(s)
    my $gatk_disable_auto_index_and_file_lock;
    my $disable_indel_qual;
    my $logging_level;
    my $pedigree_validation_type;

    my $tmpl = {
        memory_allocation => { strict_type => 1, store => \$memory_allocation },
        java_use_large_pages => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$java_use_large_pages
        },
        temp_directory => { strict_type => 1, store => \$temp_directory },
        java_jar       => { strict_type => 1, store => \$java_jar },
        intervals_ref =>
          { default => [], strict_type => 1, store => \$intervals_ref },
        read_filters_ref =>
          { default => [], strict_type => 1, store => \$read_filters_ref },
        static_quantized_quals_ref => {
            default     => [],
            strict_type => 1,
            store       => \$static_quantized_quals_ref
        },
        infile_path => {
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
        mode => {
            required    => 1,
            defined     => 1,
            allow       => [qw{ SNP INDEL BOTH }],
            strict_type => 1,
            store       => \$mode
        },
        base_quality_score_recalibration_file => {
            strict_type => 1,
            store       => \$base_quality_score_recalibration_file
        },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        FILEHANDLE => { store       => \$FILEHANDLE },
        pedigree   => { strict_type => 1, store => \$pedigree },
        tranches_file_path =>
          { strict_type => 1, store => \$tranches_file_path },
        recal_file_path => { strict_type => 1, store => \$recal_file_path },
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
        pedigree_validation_type => {
            default     => q{SILENT},
            allow       => [qw{ STRICT SILENT}],
            strict_type => 1,
            store       => \$pedigree_validation_type
        },
        ts_filter_level => {
            default     => 99.0,
            allow       => qr/ ^\d+$ | ^\d+.\d+$ /sxm,
            strict_type => 1,
            store       => \$ts_filter_level
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands;

    if ($java_jar) {    # Write java core commands to filehandle.
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
            analysis_type            => q{ApplyRecalibration},
            logging_level            => $logging_level,
            intervals_ref            => $intervals_ref,
            referencefile_path       => $referencefile_path,
            pedigree                 => $pedigree,
            pedigree_validation_type => $pedigree_validation_type,
            downsample_to_coverage   => $downsample_to_coverage,
            gatk_disable_auto_index_and_file_lock =>
              $gatk_disable_auto_index_and_file_lock,
            base_quality_score_recalibration_file =>
              $base_quality_score_recalibration_file,
            disable_indel_qual         => $disable_indel_qual,
            static_quantized_quals_ref => $static_quantized_quals_ref,
        }
    );

    if ($num_cpu_threads_per_data_thread) {

        push @commands,
            q{--num_cpu_threads_per_data_thread}
          . $SPACE
          . $num_cpu_threads_per_data_thread;
    }

    if ( @{$read_filters_ref} ) {

        push @commands,
          q{--read_filter} . $SPACE . join $SPACE . q{--read_filter} . $SPACE,
          @{$read_filters_ref};
    }

    if ($ts_filter_level) {

        push @commands, q{--ts_filter_level} . $SPACE . $ts_filter_level;
    }

    if ($recal_file_path) {

        push @commands, q{--recal_file} . $SPACE . $recal_file_path;
    }

    if ($tranches_file_path) {

        push @commands, q{--tranches_file} . $SPACE . $tranches_file_path;
    }

    if ($mode) {

        push @commands, q{--mode} . $SPACE . $mode;
    }

    ## Infile
    if ($infile_path) {

        push @commands, q{--input} . $SPACE . $infile_path;
    }

    ## Output
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;

}

sub gatk_calculategenotypeposteriors {

## Function : Perl wrapper for writing GATK calculategenotypeposteriors recipe to $FILEHANDLE. Based on GATK 3.7.0.
## Returns  : @commands
## Arguments: $memory_allocation                     => Memory allocation to run Gatk
##          : $java_use_large_pages                  => Use java large pages
##          : $temp_directory                        => Redirect tmp files to java temp
##          : $java_jar                              => Java jar
##          : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##          : $infile_path                           => Infile paths
##          : $outfile_path                          => Outfile path
##          : $referencefile_path                    => Reference sequence file
##          : $stderrfile_path                       => Stderrfile path
##          : $FILEHANDLE                            => Sbatch filehandle to write to
##          : $pedigree                              => Pedigree files for samples
##          : $supporting_callset_file_path          => Other callsets to use in generating genotype posteriors
##          : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##          : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##          : $logging_level                         => Set the minimum level of logging
##          : $pedigree_validation_type              => Validation strictness for pedigree

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $memory_allocation;
    my $java_use_large_pages;
    my $temp_directory;
    my $java_jar;
    my $intervals_ref;
    my $infile_path;
    my $outfile_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $pedigree;
    my $supporting_callset_file_path;
    my $downsample_to_coverage;

    ## Default(s)
    my $gatk_disable_auto_index_and_file_lock;
    my $logging_level;
    my $pedigree_validation_type;

    my $tmpl = {
        memory_allocation => { strict_type => 1, store => \$memory_allocation },
        java_use_large_pages => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$java_use_large_pages
        },
        temp_directory => { strict_type => 1, store => \$temp_directory },
        java_jar       => { strict_type => 1, store => \$java_jar },
        intervals_ref =>
          { default => [], strict_type => 1, store => \$intervals_ref },
        infile_path => {
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
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        FILEHANDLE => { store       => \$FILEHANDLE },
        pedigree   => { strict_type => 1, store => \$pedigree },
        supporting_callset_file_path =>
          { strict_type => 1, store => \$supporting_callset_file_path },
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
        pedigree_validation_type => {
            default     => q{SILENT},
            allow       => [qw{ STRICT SILENT }],
            strict_type => 1,
            store       => \$pedigree_validation_type
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands;

    if ($java_jar) {    # Write java core commands to filehandle.
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
            analysis_type            => q{CalculateGenotypePosteriors},
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

    ## Tool-specific options
    if ($supporting_callset_file_path) {

        push @commands,
          q{--supporting} . $SPACE . $supporting_callset_file_path;
    }

    ## Infile
    if ($infile_path) {

        push @commands, q{--variant} . $SPACE . $infile_path;
    }

    ## Output
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;

}

sub gatk_combinevariants {

## Function : Perl wrapper for writing GATK combinevariants recipe to $FILEHANDLE. Based on GATK 3.7.0.
## Returns  : @commands
## Arguments: $memory_allocation                     => Memory allocation to run Gatk
##          : $java_use_large_pages                  => Use java large pages
##          : $temp_directory                        => Redirect tmp files to java temp
##          : $java_jar                              => Java jar
##          : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##          : $infile_paths_ref                      => Infile paths {REF}
##          : $outfile_path                          => Outfile path
##          : $referencefile_path                    => Reference sequence file
##          : $stderrfile_path                       => Stderrfile path
##          : $FILEHANDLE                            => Sbatch filehandle to write to
##          : $pedigree                              => Pedigree files for samples
##          : $prioritize_caller                     => Comma seperated string specifying priority for merging
##          : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##          : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##          : $logging_level                         => Set the minimum level of logging
##          : $pedigree_validation_type              => Validation strictness for pedigree
##          : $exclude_nonvariants                   => Exclude non-variant sites
##          : $genotype_merge_option                 => Determines how we should merge genotype records for samples shared across the ROD files

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $memory_allocation;
    my $java_use_large_pages;
    my $temp_directory;
    my $java_jar;
    my $intervals_ref;
    my $infile_paths_ref;
    my $outfile_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $pedigree;
    my $prioritize_caller;
    my $genotype_merge_option;
    my $downsample_to_coverage;

    ## Default(s)
    my $gatk_disable_auto_index_and_file_lock;
    my $logging_level;
    my $pedigree_validation_type;
    my $exclude_nonvariants;

    my $tmpl = {
        memory_allocation => { strict_type => 1, store => \$memory_allocation },
        java_use_large_pages => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$java_use_large_pages
        },
        temp_directory => { strict_type => 1, store => \$temp_directory },
        java_jar       => { strict_type => 1, store => \$java_jar },
        intervals_ref =>
          { default => [], strict_type => 1, store => \$intervals_ref },
        infile_paths_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$infile_paths_ref
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
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        FILEHANDLE        => { store       => \$FILEHANDLE },
        pedigree          => { strict_type => 1, store => \$pedigree },
        prioritize_caller => { strict_type => 1, store => \$prioritize_caller },
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
        pedigree_validation_type => {
            default     => q{SILENT},
            allow       => [qw{ STRICT SILENT }],
            strict_type => 1,
            store       => \$pedigree_validation_type
        },
        exclude_nonvariants => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$exclude_nonvariants
        },
        genotype_merge_option => {
            default => q{PRIORITIZE},
            allow =>
              [ undef, qw{ UNIQUIFY PRIORITIZE UNSORTED REQUIRE_UNIQUE } ],
            strict_type => 1,
            store       => \$genotype_merge_option
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands;

    if ($java_jar) {    # Write java core commands to filehandle.
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
            analysis_type            => q{CombineVariants},
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

    ## Tool-specific options
    if ($exclude_nonvariants) {

        push @commands, q{--excludeNonVariants};
    }

    if ($genotype_merge_option) {

        push @commands,
          q{--genotypemergeoption} . $SPACE . $genotype_merge_option;
    }

    if ($prioritize_caller) {

        push @commands, q{--rod_priority_list} . $SPACE . $prioritize_caller;
    }

    ## Infile
    if ( @{$infile_paths_ref} ) {

        push @commands,
          q{--variant:} . join $SPACE . q{--variant:},
          @{$infile_paths_ref};

    }

    ## Output
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;

}

sub gatk_varianteval {

## Function : Perl wrapper for writing GATK varianteval recipe to $FILEHANDLE. Based on GATK 3.7.0.
## Returns  : @commands
## Arguments: $memory_allocation                     => Memory allocation to run Gatk
##          : $java_use_large_pages                  => Use java large pages
##          : $temp_directory                        => Redirect tmp files to java temp
##          : $java_jar                              => Java jar
##          : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##          : $infile_paths_ref                      => Infile paths
##          : $outfile_path                          => Outfile path
##          : $referencefile_path                    => Reference sequence file
##          : $stderrfile_path                       => Stderrfile path
##          : $FILEHANDLE                            => Sbatch filehandle to write to
##          : $pedigree                              => Pedigree files for samples
##          : $dbsnp_file_path                       => DbSNP file path
##          : $indel_gold_standard_file_path         => Evaluations that count calls at sites of true variation (e.g., indel calls) will use this argument as their gold standard for comparison
##          : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##          : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##          : $logging_level                         => Set the minimum level of logging
##          : $pedigree_validation_type              => Validation strictness for pedigree

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $memory_allocation;
    my $java_use_large_pages;
    my $temp_directory;
    my $java_jar;
    my $intervals_ref;
    my $infile_paths_ref;
    my $outfile_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $pedigree;
    my $dbsnp_file_path;
    my $indel_gold_standard_file_path;
    my $downsample_to_coverage;

    ## Default(s)
    my $gatk_disable_auto_index_and_file_lock;
    my $logging_level;
    my $pedigree_validation_type;

    my $tmpl = {
        memory_allocation => { strict_type => 1, store => \$memory_allocation },
        java_use_large_pages => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$java_use_large_pages
        },
        temp_directory => { strict_type => 1, store => \$temp_directory },
        java_jar       => { strict_type => 1, store => \$java_jar },
        intervals_ref =>
          { default => [], strict_type => 1, store => \$intervals_ref },
        infile_paths_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$infile_paths_ref
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
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        FILEHANDLE      => { store       => \$FILEHANDLE },
        pedigree        => { strict_type => 1, store => \$pedigree },
        dbsnp_file_path => { strict_type => 1, store => \$dbsnp_file_path },
        indel_gold_standard_file_path =>
          { strict_type => 1, store => \$indel_gold_standard_file_path },
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
            allow       => [qw{INFO ERROR FATAL }],
            strict_type => 1,
            store       => \$logging_level
        },
        pedigree_validation_type => {
            default     => q{SILENT},
            allow       => [qw{ STRICT SILENT }],
            strict_type => 1,
            store       => \$pedigree_validation_type
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands;

    if ($java_jar) {    # Write java core commands to filehandle.
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
            analysis_type            => q{VariantEval},
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

    ## Tool specific options
    if ($dbsnp_file_path) {

        push @commands, q{--dbsnp} . $SPACE . $dbsnp_file_path;
    }

    if ($indel_gold_standard_file_path) {

        push @commands,
          q{--goldStandard} . $SPACE . $indel_gold_standard_file_path;
    }

    ## Infile
    if ( @{$infile_paths_ref} ) {

        push @commands, q{--eval} . $SPACE . join q{--eval} . $SPACE,
          @{$infile_paths_ref};
    }

    ## Output
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;

}

sub gatk_leftalignandtrimvariants {

## Function : Perl wrapper for writing GATK leftalignandtrimvariants recipe to $FILEHANDLE. Based on GATK 3.7.0.
## Returns  : @commands
## Arguments: $memory_allocation                     => Memory allocation to run Gatk
##          : $java_use_large_pages                  => Use java large pages
##          : $temp_directory                        => Redirect tmp files to java temp
##          : $java_jar                              => Java jar
##          : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##          : $infile_path                           => Infile path
##          : $outfile_path                          => Outfile path
##          : $referencefile_path                    => Reference sequence file
##          : $stderrfile_path                       => Stderrfile path
##          : $FILEHANDLE                            => Sbatch filehandle to write to
##          : $pedigree                              => Pedigree files for samples
##          : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##          : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##          : $logging_level                         => Set the minimum level of logging
##          : $pedigree_validation_type              => Validation strictness for pedigree
##          : $split_multiallelics                   => Split multiallelic records and left-align individual alleles

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $memory_allocation;
    my $java_use_large_pages;
    my $temp_directory;
    my $java_jar;
    my $intervals_ref;
    my $infile_path;
    my $outfile_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $pedigree;
    my $downsample_to_coverage;

    ## Default(s)
    my $gatk_disable_auto_index_and_file_lock;
    my $logging_level;
    my $pedigree_validation_type;
    my $split_multiallelics;

    my $tmpl = {
        memory_allocation => { strict_type => 1, store => \$memory_allocation },
        java_use_large_pages => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$java_use_large_pages
        },
        temp_directory => { strict_type => 1, store => \$temp_directory },
        java_jar       => { strict_type => 1, store => \$java_jar },
        intervals_ref =>
          { default => [], strict_type => 1, store => \$intervals_ref },
        infile_path => {
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
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        FILEHANDLE => { store       => \$FILEHANDLE },
        pedigree   => { strict_type => 1, store => \$pedigree },
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
        pedigree_validation_type => {
            default     => q{SILENT},
            allow       => [qw{ STRICT SILENT }],
            strict_type => 1,
            store       => \$pedigree_validation_type
        },
        split_multiallelics => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$split_multiallelics
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands;

    if ($java_jar) {    # Write java core commands to filehandle.
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
            analysis_type            => q{LeftAlignAndTrimVariants},
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

    ## Tool-specific options
    if ($split_multiallelics) {

        push @commands, q{--splitMultiallelics};
    }

    ## Infile
    if ($infile_path) {

        push @commands, q{--variant} . $SPACE . $infile_path;

    }

    ## Output
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;

}

sub gatk_concatenate_variants {

## Function : Writes sbatch code to supplied filehandle to concatenate variants in vcf format. Each array element is combined with the infile prefix and postfix.
## Returns  :
## Arguments: $active_parameter_href  => Active parameters for this analysis hash {REF}
##          : $continue               => Adds an ampersand to the end of the command
##          : $FILEHANDLE             => SBATCH script FILEHANDLE to print to
##          : $elements_ref           => Holding the number and part of file names to be combined
##          : $infile_prefix          => Will be combined with the each array element
##          : $infile_postfix         => Will be combined with the each array element
##          : $outfile_path_prefix    => Combined outfile path prefix
##          : $outfile_suffix         => Combined outfile suffix
##          : $human_genome_reference => Human genome reference {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $continue;
    my $elements_ref;
    my $FILEHANDLE;
    my $infile_prefix;
    my $infile_postfix;
    my $outfile_path_prefix;
    my $outfile_suffix;

    ## Default(s)
    my $human_genome_reference;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        continue => {
            allow => [ 0, 1 ],
            store => \$continue,
        },
        elements_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$elements_ref
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE, },
        infile_prefix => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_prefix
        },
        infile_postfix => { strict_type => 1, store => \$infile_postfix },
        outfile_path_prefix =>
          { strict_type => 1, store => \$outfile_path_prefix, },
        outfile_suffix => {
            default     => q{.vcf},
            allow       => [qw{ .vcf }],
            strict_type => 1,
            store       => \$outfile_suffix,
        },
        human_genome_reference => {
            default =>
              $arg_href->{active_parameter_href}{human_genome_reference},
            strict_type => 1,
            store       => \$human_genome_reference
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Variantcalling::Gatk qw(gatk_catvariants);

    ## Outfile path be built
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

    say {$FILEHANDLE} q{## GATK CatVariants};

    ## Assemble infile paths
    my @infile_paths =
      map { $infile_prefix . $_ . $infile_postfix } @{$elements_ref};

    my $gatk_jar =
        catfile( $active_parameter_href->{gatk_path}, q{GenomeAnalysisTK.jar} )
      . $SPACE
      . q{org.broadinstitute.gatk.tools.CatVariants};
    gatk_catvariants(
        {
            memory_allocation => q{Xmx4g},
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            temp_directory     => $active_parameter_href->{temp_directory},
            gatk_path          => $gatk_jar,
            infile_paths_ref   => \@infile_paths,
            outfile_path       => $outfile_path,
            referencefile_path => $human_genome_reference,
            logging_level      => $active_parameter_href->{gatk_logging_level},
            FILEHANDLE         => $FILEHANDLE,
        }
    );

    ## Launch multiple processes
    if ($continue) {
        print {$FILEHANDLE} $AMPERSAND;
    }

    say {$FILEHANDLE} $NEWLINE;
    return;
}

sub gatk_asereadcounter {

## Function : Perl wrapper for writing GATK ASEReadCounter recipe to $FILEHANDLE. Based on GATK 3.8.0.
## Returns  : @commands
## Arguments: $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##          : $FILEHANDLE                            => Sbatch filehandle to write to
##          : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##          : $gatk_sites_vcffile                    => Sites for which to count allele specific expression
##          : $infile_path                           => Infile path
##          : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##          : $java_jar                              => Java jar
##          : $java_use_large_pages                  => Use java large pages
##          : $logging_level                         => Set the minimum level of logging
##          : $memory_allocation                     => Memory allocation to run Gatk
##          : $num_cpu_threads_per_data_thread       => Number of CPU threads to allocate per data thread
##          : $outfile_path                          => Outfile path
##          : $pedigree                              => Pedigree files for samples
##          : $pedigree_validation_type              => Validation strictness for pedigree
##          : $referencefile_path                    => Reference sequence file
##          : $stderrfile_path                       => Stderrfile path
##          : $temp_directory                        => Redirect tmp files to java temp

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $downsample_to_coverage;
    my $FILEHANDLE;
    my $gatk_sites_vcffile;
    my $infile_path;
    my $intervals_ref;
    my $java_jar;
    my $memory_allocation;
    my $num_cpu_threads_per_data_thread;
    my $outfile_path;
    my $pedigree;
    my $referencefile_path;
    my $stderrfile_path;
    my $temp_directory;

    ## Default(s)
    my $gatk_disable_auto_index_and_file_lock;
    my $java_use_large_pages;
    my $logging_level;
    my $pedigree_validation_type;

    my $tmpl = {
        downsample_to_coverage => {
            allow       => qr/ ^\d+$ /sxm,
            store       => \$downsample_to_coverage,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        gatk_disable_auto_index_and_file_lock => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$gatk_disable_auto_index_and_file_lock,
            strict_type => 1,
        },
        gatk_sites_vcffile => {
            defined     => 1,
            required    => 1,
            strict_type => 1,
            store       => \$gatk_sites_vcffile,
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
        java_jar => {
            store       => \$java_jar,
            strict_type => 1,
        },
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
        memory_allocation => {
            store       => \$memory_allocation,
            strict_type => 1,
        },
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
        pedigree => {
            store       => \$pedigree,
            strict_type => 1,
        },
        pedigree_validation_type => {
            allow       => [qw{ STRICT SILENT }],
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
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands;

    ## Write java core commands to filehandle
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
            analysis_type          => q{ASEReadCounter},
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

    ## Tool-specific options
    push @commands, q{--sitesVCFFile} . $SPACE . $gatk_sites_vcffile;

    ## Infile
    if ($infile_path) {
        push @commands, q{--input_file} . $SPACE . $infile_path;
    }

    ## Output
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
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );

    return @commands;
}

sub gatk_variantfiltration {

## Function : Perl wrapper for writing GATK VariantFiltration recipe to $FILEHANDLE. Based on GATK 3.8.0.
## Returns  : @commands
## Arguments: $cluster_size                          => Number of SNPs which make up a cluster
##          : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##          : $FILEHANDLE                            => Sbatch filehandle to write to
##          : $filter_href                           => Hash with the name of the filter as key and the filter expression as value {REF}
##          : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##          : $infile_path                           => Infile paths
##          : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##          : $java_jar                              => Java jar
##          : $java_use_large_pages                  => Use java large pages
##          : $logging_level                         => Set the minimum level of logging
##          : $memory_allocation                     => Memory allocation to run Gatk
##          : $outfile_path                          => Outfile path
##          : $pedigree                              => Pedigree files for samples
##          : $pedigree_validation_type              => Validation strictness for pedigree
##          : $referencefile_path                    => Reference sequence file
##          : $stderrfile_path                       => Stderrfile path
##          : $temp_directory                        => Redirect tmp files to java temp
##          : $cluster_window_size                   => Window size (in bases) in which to evaluate clustered SNPs

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $cluster_size;
    my $downsample_to_coverage;
    my $FILEHANDLE;
    my $filter_href;
    my $infile_path;
    my $intervals_ref;
    my $java_jar;
    my $java_use_large_pages;
    my $memory_allocation;
    my $outfile_path;
    my $pedigree;
    my $referencefile_path;
    my $stderrfile_path;
    my $temp_directory;
    my $cluster_window_size;

    ## Default(s)
    my $gatk_disable_auto_index_and_file_lock;
    my $logging_level;
    my $pedigree_validation_type;

    my $tmpl = {
        cluster_size => {
            allow       => qr/ ^\d+$ /sxm,
            store       => \$cluster_size,
            strict_type => 1,
        },
        downsample_to_coverage => {
            allow       => qr/ ^\d+$ /sxm,
            store       => \$downsample_to_coverage,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        filter_href => {
            default     => {},
            store       => \$filter_href,
            strict_type => 1,
        },
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
            store       => \$intervals_ref,
            strict_type => 1,
        },
        logging_level => {
            allow       => [qw{INFO ERROR FATAL }],
            default     => q{INFO},
            store       => \$logging_level,
            strict_type => 1,
        },
        java_jar => {
            store       => \$java_jar,
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
        pedigree => {
            store       => \$pedigree,
            strict_type => 1,
        },
        pedigree_validation_type => {
            allow       => [qw{ STRICT SILENT }],
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
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
        cluster_window_size => {
            allow       => qr/ ^\d+$ /xms,
            store       => \$cluster_window_size,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands;

    # Write java core commands to filehandle.
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
            analysis_type          => q{VariantFiltration},
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

    ## Tool specific options

    ## Infile
    push @commands, q{--variant} . $SPACE . $infile_path;

    ## Output
    push @commands, q{--out} . $SPACE . $outfile_path;

    if ($cluster_size) {
        push @commands, q{--clusterSize} . $SPACE . $cluster_size;
    }
    if ($filter_href) {
      FILTERNAME:
        foreach my $filtername ( keys %{$filter_href} ) {
            push @commands,
                q{--filterName}
              . $SPACE
              . $filtername
              . $SPACE
              . q{--filterExpression}
              . $SPACE
              . $DOUBLE_QOUTE
              . $filter_href->{$filtername}
              . $DOUBLE_QOUTE;
        }
    }
    if ($cluster_window_size) {
        push @commands, q{--clusterWindowSize} . $SPACE . $cluster_window_size;
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;
}

sub gatk_genomicsdbimport {

## Function : Perl wrapper for writing GATK GenomicsDBImport recipe to $FILEHANDLE. Based on GATK 4.0.8
## Returns  : @commands
## Arguments: $FILEHANDLE                => Sbatch filehandle to write to
##          : $genomicsdb_workspace_path => Workspace for GenomicsDB
##          : $intervals_ref             => One or more genomic intervals over which to operate {REF}
##          : $java_use_large_pages      => Use java large pages
##          : $memory_allocation         => Memory allocation to run Gatk
##          : $pedigree                  => Pedigree files
##          : $read_filters_ref          => Filters to apply on reads {REF}
##          : $referencefile_path        => Reference sequence file
##          : $stderrfile_path           => Stderrfile path
##          : $temp_directory            => Redirect tmp files to java temp
##          : $variantfile_paths_ref     => GVCF files to be imported to GenomicsDB {REF}
##          : $verbosity                 => Set the minimum level of logging
##          : $xargs_mode                => Set if the program will be executed via xargs

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $genomicsdb_workspace_path;
    my $intervals_ref;
    my $java_use_large_pages;
    my $memory_allocation;
    my $read_filters_ref;
    my $referencefile_path;
    my $stderrfile_path;
    my $temp_directory;
    my $variantfile_paths_ref;

    ## Default(s)
    my $verbosity;
    my $xargs_mode;

    my $tmpl = {
        FILEHANDLE                => { store => \$FILEHANDLE },
        genomicsdb_workspace_path => {
            defined     => 1,
            required    => 1,
            store       => \$genomicsdb_workspace_path,
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
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        temp_directory  => { store => \$temp_directory,  strict_type => 1, },
        variantfile_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$variantfile_paths_ref,
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

    ## GATK GenomicsDBImport

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
    push @commands, q{GenomicsDBImport};

    ## Add GVCF files
    push @commands,
      q{--variant} . $SPACE . join $SPACE . q{--variant} . $SPACE,
      @{$variantfile_paths_ref};

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
    push @commands,
      q{--genomicsdb-workspace-path} . $SPACE . $genomicsdb_workspace_path;

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
