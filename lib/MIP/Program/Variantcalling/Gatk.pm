package MIP::Program::Variantcalling::Gatk;

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
use FindBin qw{ $Bin };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };

## CPANM
use Readonly;

## MIPs lib/
use MIP::Constants
  qw{ $ASTERISK $AMPERSAND $COLON $DOT $DOUBLE_QUOTE $EMPTY_STR $NEWLINE $SPACE $UNDERSCORE };
use MIP::Language::Java qw{ java_core };
use MIP::Program::Base::Gatk qw{ gatk_base gatk_common_options gatk_java_options };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.18;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      gatk_applyvqsr
      gatk_calculategenotypeposteriors
      gatk_cnnscorevariants
      gatk_combinevariants
      gatk_concatenate_variants
      gatk_gathervcfscloud
      gatk_genomicsdbimport
      gatk_genotypegvcfs
      gatk_indexfeaturefile
      gatk_leftalignandtrimvariants
      gatk_selectvariants
      gatk_varianteval
      gatk_variantfiltration
      gatk_variantrecalibrator
    };

}

sub gatk_genotypegvcfs {

## Function : Perl wrapper for writing GATK GenoTypeGVCFs recipe to $FILEHANDLE. Based on GATK 4.1.0.
## Returns  : @commands
## Arguments: $dbsnp_path               => Path to DbSNP file
##          : $FILEHANDLE               => Sbatch filehandle to write to
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
##          : $use_new_qual_calculator  => Use the new AF model instead of the so-called exact model
##          : $verbosity                => Set the minimum level of logging
##          : $xargs_mode               => Set if the program will be executed via xargs

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $dbsnp_path;
    my $FILEHANDLE;
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
    my $use_new_qual_calculator;
    my $verbosity;
    my $xargs_mode;

    my $tmpl = {
        dbsnp_path => {
            store       => \$dbsnp_path,
            strict_type => 1,
        },
        FILEHANDLE               => { store => \$FILEHANDLE, },
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

    ## GATK GenotyeGVCFs

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
    push @commands, q{GenotypeGVCFs};

    ## Add infile
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

    if ($use_new_qual_calculator) {

        push @commands, q{--use-new-qual-calculator};
    }

    ## Add dbsnp
    if ($dbsnp_path) {
        push @commands, q{--dbsnp} . $SPACE . $dbsnp_path;
    }

    ## Add output
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;
}

sub gatk_selectvariants {

## Function : Perl wrapper for writing GATK SelectVariants recipe to $FILEHANDLE. Based on GATK 4.0.8
## Returns  : @commands
## Arguments: $exclude_non_variants       => Exclude non-variant sites
##          : $FILEHANDLE                 => Sbatch filehandle to write to
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
    my $FILEHANDLE;
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
        FILEHANDLE  => { store => \$FILEHANDLE, },
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

    ## GATK SelectVariants

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
    push @commands, q{SelectVariants};

    ## Add infile
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

    ## Exclude non variants
    if ($exclude_non_variants) {
        push @commands, q{--exclude-non-variants};
    }

    ## Restrict allele output
    if ($restrict_alleles_to) {
        push @commands, q{--restrict-alleles-to} . $SPACE . $restrict_alleles_to;
    }

    ## Sample names to include
    if ( @{$sample_names_ref} ) {
        push @commands,
          q{--sample-name} . $SPACE . join $SPACE . q{--sample-name} . $SPACE,
          @{$sample_names_ref};
    }

    ## Select variant types to include
    if ( @{$select_type_to_include_ref} ) {
        push @commands,
            q{--select-type-to-include}
          . $SPACE
          . join $SPACE
          . q{--select-type-to-include}
          . $SPACE, @{$select_type_to_include_ref};
    }

    ## Output
    push @commands, q{--output} . $SPACE . $outfile_path;

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

sub gatk_variantrecalibrator {

## Function : Perl wrapper for writing GATK variantrecalibrator recipe to $FILEHANDLE. Based on GATK 4.1.0.
## Returns  : @commands
##          : $annotations_ref       => One or more specific annotations to apply to variant calls
##          : $FILEHANDLE            => Sbatch filehandle to write to
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
##          : $ts_tranches_ref            => Levels of truth sensitivity at which to slice the data. (in percent, that is 1.0 for 1 percent)
##          : $verbosity	     => Set the minimum level of logging
##          : $xargs_mode            => Set if the program will be executed via xargs

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $annotations_ref;
    my $FILEHANDLE;
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
        FILEHANDLE  => { store => \$FILEHANDLE, },
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
            allow       => [qr/ ^\d+$ /sxm],
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

    ## GATK VariantRecalibrator

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
    push @commands, q{VariantRecalibrator};

    ## Add infile
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

    ## Add annotations
    push @commands,
      q{--use-annotation} . $SPACE . join $SPACE . q{--use-annotation} . $SPACE,
      @{$annotations_ref};

    ## Add max attempts to build model
    if ($max_attempts) {
        push @commands, q{--max-attempts} . $SPACE . $max_attempts;
    }

    ## Add max gaussians for positive model
    if ($max_gaussian_level) {
        push @commands, q{--max-gaussians} . $SPACE . $max_gaussian_level;
    }

    ## Add mode
    if ($mode) {
        push @commands, q{--mode} . $SPACE . $mode;
    }

    ## Add path to r-script
    if ($rscript_file_path) {
        push @commands, q{--rscript-file} . $SPACE . $rscript_file_path;
    }

    ## Add list of resources
    push @commands,
      q{--resource} . $COLON . join $SPACE . q{--resource} . $COLON,
      @{$resources_ref};

    if ($ts_tranches_ref) {

        push @commands, q{-tranche} . $SPACE . join $SPACE . q{-tranche} . $SPACE,
          @{$ts_tranches_ref};
    }

    ## Add path to tranches file
    push @commands, q{--tranches-file} . $SPACE . $tranches_file_path;

    ## Trust all training sites to be polymorphic
    if ($trust_all_polymorphic) {
        push @commands, q{--trust-all-polymorphic};
    }

    ## Add path to output recal file
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;
}

sub gatk_applyvqsr {

## Function : Perl wrapper for writing GATK ApplyVQSR recipe to $FILEHANDLE. Based on GATK 4.0.8.
## Returns  : @commands
## Arguments: $FILEHANDLE                            => Sbatch filehandle to write to
##          : $infile_path                           => Infile paths
##          : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##          : $java_use_large_pages                  => Use java large pages
##          : $memory_allocation                     => Memory allocation to run Gatk
##          : $mode                                  => Mode for emitting reference confidence scores
##          : $outfile_path                          => Outfile path
##          : $read_filters_ref                      => Filters to apply to reads before analysis {REF}
##          : $recal_file_path                       => The output recal file used by ApplyVQSR
##          : $referencefile_path                    => Reference sequence file
##          : $stderrfile_path                       => Stderrfile path
##          : $temp_directory                        => Redirect tmp files to java temp
##          : $tranches_file_path                    => The output tranches file used by ApplyRecalibration
##          : $ts_filter_level                       => Ts filter level
##          : $verbosity	                         => Set the minimum level of logging
##          : $xargs_mode   		                 => Set if the program will be executed via xargs

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
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
        FILEHANDLE => {
            store => \$FILEHANDLE,
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
            allow       => qr/ ^\d+$ | ^\d+.\d+$ /sxm,
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

    ## GATK ApplyVQSR

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
    push @commands, q{ApplyVQSR};

    ## Add infile
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

    ## Set run mode
    if ($mode) {
        push @commands, q{--mode} . $SPACE . $mode;
    }

    ## Add recalibration file
    push @commands, q{--recal-file} . $SPACE . $recal_file_path;

    ## Add tranches file
    push @commands, q{--tranches-file} . $SPACE . $tranches_file_path;

    ## Add truth sensitivity level
    if ($ts_filter_level) {
        push @commands, q{--truth-sensitivity-filter-level} . $SPACE . $ts_filter_level;
    }

    ## Add outfile
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;
}

sub gatk_calculategenotypeposteriors {

## Function : Perl wrapper for writing GATK CalculateGenotypePosteriors recipe to $FILEHANDLE. Based on GATK 4.1.0.
## Returns  : @commands
##          : $FILEHANDLE                            => Sbatch filehandle to write to
##          : $infile_path                           => Infile paths
##          : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##          : $java_use_large_pages                  => Use java large pages
##          : $memory_allocation                     => Memory allocation to run Gatk
##          : $num_ref_samples_if_no_call            => Number of hom-ref genotypes to infer at sites not present in a panel
##          : $outfile_path                          => Outfile path
##          : $pedigree                              => Pedigree files for samples
##          : $read_filters_ref                      => Filters to apply on reads {REF}
##          : $referencefile_path                    => Reference sequence file
##          : $stderrfile_path                       => Stderrfile path
##          : $supporting_callset_file_path          => Other callsets to use in generating genotype posteriors
##          : $temp_directory                        => Redirect tmp files to java temp
##          : $verbosity                             => Set the minimum level of logging
##          : $xargs_mode                            => Set if the program will be executed via xargs

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
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
        FILEHANDLE => {
            store => \$FILEHANDLE,
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

    ## GATK CalculateGenotypePosteriors

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
    push @commands, q{CalculateGenotypePosteriors};

    ## Add infile
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

    ## Add supporting data
    if ($supporting_callset_file_path) {
        push @commands, q{--supporting-callsets} . $SPACE . $supporting_callset_file_path;
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;
}

sub gatk_cnnscorevariants {

## Function : Perl wrapper for writing GATK CNNScoreVariants recipe to $FILEHANDLE. Based on GATK 4.1.0
## Returns  : @commands
## Arguments: $alignment_infile_paths_ref => BAM/SAM/CRAM file containing reads {REF}
##          : $FILEHANDLE                 => Sbatch filehandle to write to
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
    my $FILEHANDLE;
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
        FILEHANDLE  => { store => \$FILEHANDLE, },
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

    ## GATK SelectVariants

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
    push @commands, q{CNNScoreVariants};

    ## Add infile
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

    ## Output
    push @commands, q{--output} . $SPACE . $outfile_path;

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
        memory_allocation    => { strict_type => 1, store => \$memory_allocation },
        java_use_large_pages => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$java_use_large_pages
        },
        temp_directory => { strict_type => 1, store => \$temp_directory },
        java_jar       => { strict_type => 1, store => \$java_jar },
        intervals_ref    => { default => [], strict_type => 1, store => \$intervals_ref },
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
        stderrfile_path        => { strict_type => 1, store => \$stderrfile_path },
        FILEHANDLE             => { store       => \$FILEHANDLE },
        pedigree               => { strict_type => 1, store => \$pedigree },
        prioritize_caller      => { strict_type => 1, store => \$prioritize_caller },
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
            default     => q{PRIORITIZE},
            allow       => [ undef, qw{ UNIQUIFY PRIORITIZE UNSORTED REQUIRE_UNIQUE } ],
            strict_type => 1,
            store       => \$genotype_merge_option
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ gatk3 };

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

        push @commands, q{--genotypemergeoption} . $SPACE . $genotype_merge_option;
    }

    if ($prioritize_caller) {

        push @commands, q{--rod_priority_list} . $SPACE . $prioritize_caller;
    }

    ## Infile
    if ( @{$infile_paths_ref} ) {

        push @commands, q{--variant:} . join $SPACE . q{--variant:}, @{$infile_paths_ref};

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

## Function : Perl wrapper for writing GATK varianteval recipe to $FILEHANDLE. Based on GATK 4.1.0.
## Returns  : @commands
## Arguments: $dbsnp_file_path                       => DbSNP file path
##          : $FILEHANDLE                            => Sbatch filehandle to write to
##          : $indel_gold_standard_file_path         => Evaluations that count calls at sites of true variation (e.g., indel calls) will use this argument as their gold standard for comparison
##          : $infile_paths_ref                      => Infile paths
##          : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##          : $java_use_large_pages                  => Use java large pages
##          : $memory_allocation                     => Memory allocation to run Gatk
##          : $outfile_path                          => Outfile path
##          : $pedigree                              => Pedigree files for samples
##          : $referencefile_path                    => Reference sequence file
##          : $stderrfile_path                       => Stderrfile path
##          : $temp_directory                        => Redirect tmp files to java temp
   #          : $verbosity                             => Set the minimum level of logging
##          : $xargs_mode                            => Set if the program will be executed via xargs

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $dbsnp_file_path;
    my $FILEHANDLE;
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
        FILEHANDLE      => { store => \$FILEHANDLE, },
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

    ## GATK VariantEval

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

    ## Tool specific options
    if ($dbsnp_file_path) {

        push @commands, q{--dbsnp} . $SPACE . $dbsnp_file_path;
    }

    if ($indel_gold_standard_file_path) {

        push @commands, q{--gold-standard} . $SPACE . $indel_gold_standard_file_path;
    }

    ## Infile
    if ( @{$infile_paths_ref} ) {

        push @commands, q{--eval} . $SPACE . join $SPACE . q{--eval} . $SPACE,
          @{$infile_paths_ref};
    }

    ## Output
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
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );

    return @commands;
}

sub gatk_leftalignandtrimvariants {

## Function : Perl wrapper for writing GATK LeftAlignAndTrimVariants recipe to $FILEHANDLE. Based on GATK 4.0.0.
## Returns  : @commands
##          : $FILEHANDLE                            => Sbatch filehandle to write to
##          : $infile_path                           => Infile path
##          : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##          : $java_use_large_pages                  => Use java large pages
##          : $memory_allocation                     => Memory allocation to run Gatk
##          : $outfile_path                          => Outfile path
##          : $referencefile_path                    => Reference sequence file
##          : $split_multiallelics                   => Split multiallelic records and left-align individual alleles
##          : $stderrfile_path                       => Stderrfile path
##          : $temp_directory                        => Redirect tmp files to java temp
##          : $verbosity                             => Set the minimum level of logging
##          : $xargs_mode                            => Set if the program will be executed via xargs

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
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
        FILEHANDLE  => { store => \$FILEHANDLE },
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

    ## GATK LeftAlignAndTrimVariants

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
    push @commands, q{LeftAlignAndTrimVariants};

    ## Add infile
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

    ## Split Multiallelic
    if ($split_multiallelics) {
        push @commands, q{--split-multi-allelics};
    }

    ## Output
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;

}

sub gatk_concatenate_variants {

## Function : Writes sbatch code to supplied filehandle to concatenate variants in vcf format. Each array element is combined with the infile prefix and postfix.
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $continue              => Adds an ampersand to the end of the command
##          : $FILEHANDLE            => SBATCH script FILEHANDLE to print to
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
    my $FILEHANDLE;
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
        FILEHANDLE => {
            defined  => 1,
            required => 1,
            store    => \$FILEHANDLE,
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
            allow       => [qw{ .vcf .selected.vcf }],
            default     => q{.vcf},
            store       => \$outfile_suffix,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Variantcalling::Gatk qw(gatk_gathervcfscloud);

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

    say {$FILEHANDLE} q{## GATK GatherVCFs};

    ## Assemble infile paths
    my @infile_paths =
      map { $infile_prefix . $DOT . $_ . $infile_postfix } @{$elements_ref};

    gatk_gathervcfscloud(
        {
            FILEHANDLE           => $FILEHANDLE,
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

        print {$FILEHANDLE} $AMPERSAND;
    }
    say {$FILEHANDLE} $NEWLINE;
    return;
}

sub gatk_variantfiltration {

## Function : Perl wrapper for writing GATK VariantFiltration recipe to $FILEHANDLE. Based on GATK 4.1.0.
## Returns  : @commands
##          : $cluster_size         => Number of SNPs which make up a cluster
##          : $cluster_window_size  => Window size (in bases) in which to evaluate clustered SNPs
##          : $FILEHANDLE           => Sbatch filehandle to write to
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
    my $FILEHANDLE;
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
            allow       => qr/ ^\d+$ /sxm,
            store       => \$cluster_size,
            strict_type => 1,
        },
        cluster_window_size => {
            allow       => qr/ ^\d+$ /xms,
            store       => \$cluster_window_size,
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

    ## GATK VariantFiltration

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
    push @commands, q{VariantFiltration};

    ## Add infile
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

    ## Add number of cluster
    if ($cluster_size) {
        push @commands, q{--cluster-size} . $SPACE . $cluster_size;
    }

    ## Add window size
    if ($cluster_window_size) {
        push @commands, q{--cluster-window-size} . $SPACE . $cluster_window_size;
    }

    ## Add filters
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
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );

    return @commands;
}

sub gatk_genomicsdbimport {

## Function : Perl wrapper for writing GATK GenomicsDBImport recipe to $FILEHANDLE. Based on GATK 4.0.8
## Returns  : @commands
## Arguments: $FILEHANDLE                => Sbatch filehandle to write to
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
    my $FILEHANDLE;
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
        FILEHANDLE                => { store => \$FILEHANDLE },
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
    if ( scalar @{$infile_paths_ref} ) {
        push @commands,
          q{--variant} . $SPACE . join $SPACE . q{--variant} . $SPACE,
          @{$infile_paths_ref};
    }

    ## Add merged reference files
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

    ## Output
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
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gatk_gathervcfscloud {

## Function : Perl wrapper for writing GATK GatherVcfsCloud recipe to $FILEHANDLE. Based on GATK 4.0.8.
## Returns  : @commands
## Arguments: $FILEHANDLE                            => Sbatch filehandle to write to
##          : $ignore_safety_checks                  => Disable sanity checks to improve performance
##          : $infile_paths_ref                      => VCF files to gather {REF}
##          : $java_use_large_pages                  => Use java large pages
##          : $memory_allocation                     => Memory allocation to run Gatk
##          : $outfile_path                          => Outfile path
##          : $stderrfile_path                       => Stderrfile path
##          : $temp_directory                        => Redirect tmp files to java temp
##          : $verbosity                             => Set the minimum level of logging

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
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
        FILEHANDLE           => { store => \$FILEHANDLE, },
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

    ## GATK GatherVcfsCloud

    # Stores commands depending on input parameters
    my @commands = qw{ gatk };

    ## Add java options
    gatk_java_options(
        {
            commands_ref         => \@commands,
            java_use_large_pages => $java_use_large_pages,
            memory_allocation    => $memory_allocation,
        }
    );

    ## Add tool command
    push @commands, q{GatherVcfsCloud};

    ## Add infile
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

    ## Disable sanity check
    if ($ignore_safety_checks) {
        push @commands, q{--ignore-safety-checks};
    }

    ## Add output
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;
}

sub gatk_indexfeaturefile {

## Function : Perl wrapper for writing GATK IndexFeatureFile recipe to $FILEHANDLE. Based on GATK 4.0.10.
## Returns  : @commands
## Arguments: $FILEHANDLE                            => Sbatch filehandle to write to
##          : $infile_path                           => Path to feature file
##          : $java_use_large_pages                  => Use java large pages
##          : $memory_allocation                     => Memory allocation to run Gatk
##          : $outfile_path                          => Path to index
##          : $stderrfile_path                       => Stderrfile path
##          : $temp_directory                        => Redirect tmp files to java temp
##          : $verbosity                             => Set the minimum level of logging

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $memory_allocation;
    my $outfile_path;
    my $stderrfile_path;
    my $temp_directory;

    ## Default(s)
    my $java_use_large_pages;
    my $verbosity;

    my $tmpl = {
        FILEHANDLE  => { store => \$FILEHANDLE, },
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

    ## GATK GatherVcfsCloud

    # Stores commands depending on input parameters
    my @commands = qw{ gatk };

    ## Add java options
    gatk_java_options(
        {
            commands_ref         => \@commands,
            java_use_large_pages => $java_use_large_pages,
            memory_allocation    => $memory_allocation,
        }
    );

    ## Add tool command
    push @commands, q{IndexFeatureFile};

    ## Add infile
    push @commands, q{--feature-file} . $SPACE . $infile_path;

    ## Add common options
    gatk_common_options(
        {
            commands_ref   => \@commands,
            temp_directory => $temp_directory,
            verbosity      => $verbosity,
        }
    );

    ## Add output path
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;
}

1;
