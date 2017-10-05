package MIP::Program::Base::Gatk;

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
use File::Spec::Functions qw{ catdir };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ gatk_base };
}

## Constants
Readonly my $SPACE => q{ };

sub gatk_base {

## pgatk_base

## Function : Perl wrapper for Gatk base. Based on Gatk v?????
## Returns  : @commands
## Arguments: $commands_ref, $FILEHANDLE, $analysis_type, $gatk_path, $referencefile_path, $infile_path, $infile_paths_ref, $sample_names_ref, $outfile_path, $num_cpu_threads_per_data_thread, $read_filters_ref, $intervals_ref, $resources_ref, annotations_ref, $pedigree, $downsample_to_coverage, $max_gaussian_level, $static_quantized_quals_ref, $mode, $base_quality_score_recalibration_file, $rscript_file_path, $recal_file_path, $tranches_file_path, $ts_filter_level, $supporting_callset_file_path, $prioritize_caller, $genotype_merge_option, $dbsnp_file_path, $indel_gold_standard_file_path, $logging_level, $gatk_disable_auto_index_and_file_lock, $pedigree_validation_type, $include_nonvariant_sites, $exclude_nonvariants, $assume_sorted, $disable_indel_qual, $split_multiallelics
##          : $commands_ref                          => List of commands added earlier
##          : $FILEHANDLE                            => Filehandle to write to
##          : $analysis_type                         => Analysis type
##          : $gatk_path                             => Path to Gatk
##          : $referencefile_path                    => Reference sequence file
##          : $infile_path                           => Infile path
##          : $infile_paths_ref                      => Infile paths {REF}
##          : $sample_names_ref                      => Include genotypes from this sample {REF}
##          : $outfile_path                          => Outfile path
##          : $num_cpu_threads_per_data_thread       => Number of CPU threads to allocate per data thread
##          : $read_filters_ref                      => Filters to apply to reads before analysis {REF}
##          : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##          : $resources_ref                         => A list of sites for which to apply a prior probability of being correct but which aren't used by the algorithm (training and truth setsare required to run)
##          : $annotations_ref                       => One or more specific annotations to apply to variant calls
##          : $pedigree                              => Pedigree files
##          : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##          : $max_gaussian_level                    => Max number of Gaussians for the positive mode
##          : $static_quantized_quals_ref            => Use static quantized quality scores to a given number of levels (with -BQSR) {REF}
##          : $mode                                  => Mode for emitting reference confidence scores
##          : $base_quality_score_recalibration_file => Input covariates table file for on-the-fly base quality score recalibration
##          : $rscript_file_path                     => Rscript file path
##          : $recal_file_path                       => The output recal file used by ApplyRecalibration
##          : $tranches_file_path                    => The output tranches file used by ApplyRecalibration
##          : $supporting_callset_file_path          => Other callsets to use in generating genotype posteriors
##          : $prioritize_caller                     => Comma seperated string specifying priority for merging
##          : $genotype_merge_option                 => Determines how we should merge genotype records for samples shared across the ROD files
##          : $dbsnp_file_path                       => DbSNP file path
##          : $indel_gold_standard_file_path         => Evaluations that count calls at sites of true variation (e.g., indel calls) will use this argument as their gold standard for comparison
##          : $logging_level                         => Logging level
##          : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##          : $pedigree_validation_type              => Validation strictness for pedigree
##          : $include_nonvariant_sites              => Include loci found to be non-variant after genotyping
##          : $exclude_nonvariants                   => Exclude non-variant sites
##          : $assume_sorted                         => Assume_sorted should be true if the input files are already sorted
##          : $disable_indel_qual                    => Disable printing of base insertion and deletion tags (with -BQSR)
##          : $split_multiallelics                   => Split multiallelic records and left-align individual alleles

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $commands_ref;
    my $FILEHANDLE;
    my $analysis_type;
    my $gatk_path;
    my $referencefile_path;
    my $infile_path;
    my $infile_paths_ref;
    my $sample_names_ref;
    my $outfile_path;
    my $num_cpu_threads_per_data_thread;
    my $read_filters_ref;
    my $intervals_ref;
    my $resources_ref;
    my $annotations_ref;
    my $pedigree;
    my $downsample_to_coverage;
    my $max_gaussian_level;
    my $static_quantized_quals_ref;
    my $mode;
    my $base_quality_score_recalibration_file;
    my $rscript_file_path;
    my $recal_file_path;
    my $tranches_file_path;
    my $ts_filter_level;
    my $supporting_callset_file_path;
    my $prioritize_caller;
    my $genotype_merge_option;
    my $dbsnp_file_path;
    my $indel_gold_standard_file_path;

    ## Default(s)
    my $logging_level;
    my $gatk_disable_auto_index_and_file_lock;
    my $pedigree_validation_type;
    my $include_nonvariant_sites;
    my $exclude_nonvariants;
    my $assume_sorted;
    my $disable_indel_qual;
    my $split_multiallelics;

    my $tmpl = {
        commands_ref =>
          { default => [], strict_type => 1, store => \$commands_ref },
        FILEHANDLE    => { store => \$FILEHANDLE },
        analysis_type => {
            defined     => 1,
            strict_type => 1,
            store       => \$analysis_type
        },
        gatk_path => {
            defined     => 1,
            strict_type => 1,
            store       => \$gatk_path
        },
        referencefile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$referencefile_path
        },
        infile_path =>
          { defined => 1, strict_type => 1, store => \$infile_path },
        infile_paths_ref => {
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$infile_paths_ref
        },
        sample_names_ref => {
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$sample_names_ref
        },
        outfile_path =>
          { defined => 1, strict_type => 1, store => \$outfile_path },
        num_cpu_threads_per_data_thread => {
            default     => undef,
            strict_type => 1,
            store       => \$num_cpu_threads_per_data_thread
        },
        read_filters_ref =>
          { default => [], strict_type => 1, store => \$read_filters_ref },
        intervals_ref =>
          { default => [], strict_type => 1, store => \$intervals_ref },
        resources_ref => {
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$resources_ref
        },
        annotations_ref => {
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$annotations_ref
        },
        pedigree                 => { strict_type => 1, store => \$pedigree },
        pedigree_validation_type => {
            default     => q{SILENT},
            allow       => [qw{ STRICT SILENT}],
            strict_type => 1,
            store       => \$pedigree_validation_type
        },
        downsample_to_coverage => {
            strict_type => 1,
            store       => \$downsample_to_coverage
        },
        max_gaussian_level => {
            allow       => [ undef, qr/^\d+$/ ],
            strict_type => 1,
            store       => \$max_gaussian_level
        },
        static_quantized_quals_ref => {
            default     => [],
            strict_type => 1,
            store       => \$static_quantized_quals_ref
        },
        mode => {
            defined     => 1,
            allow       => [qw{ SNP INDEL BOTH }],
            strict_type => 1,
            store       => \$mode
        },
        base_quality_score_recalibration_file => {
            strict_type => 1,
            store       => \$base_quality_score_recalibration_file
        },
        rscript_file_path => { strict_type => 1, store => \$rscript_file_path },
        recal_file_path   => { strict_type => 1, store => \$recal_file_path },
        tranches_file_path =>
          { strict_type => 1, store => \$tranches_file_path },
        ts_filter_level => {
            default     => undef,
            allow       => qr/^\d+$|^\d+.\d+$/,
            strict_type => 1,
            store       => \$ts_filter_level
        },
        supporting_callset_file_path =>
          { strict_type => 1, store => \$supporting_callset_file_path },
        prioritize_caller => { strict_type => 1, store => \$prioritize_caller },
        genotype_merge_option => {
            default => q{PRIORITIZE},
            allow =>
              [ undef, qw{ UNIQUIFY PRIORITIZE UNSORTED REQUIRE_UNIQUE } ],
            strict_type => 1,
            store       => \$genotype_merge_option
        },
        dbsnp_file_path => { strict_type => 1, store => \$dbsnp_file_path },
        indel_gold_standard_file_path =>
          { strict_type => 1, store => \$indel_gold_standard_file_path },
        logging_level => {
            default     => q{INFO},
            allow       => [qw{ INFO ERROR FATAL }],
            strict_type => 1,
            store       => \$logging_level
        },
        gatk_disable_auto_index_and_file_lock => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$gatk_disable_auto_index_and_file_lock
        },
        include_nonvariant_sites => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$include_nonvariant_sites
        },
        exclude_nonvariants => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$exclude_nonvariants
        },
        assume_sorted => {
            default     => undef,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$assume_sorted
        },
        disable_indel_qual => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$disable_indel_qual
        },
        split_multiallelics => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$split_multiallelics
        },

    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = @{$commands_ref};

    # catvariants does not need --analysis_type to be written in the commands
    if ( $analysis_type ne q{catvariants} ) {

        push @commands, q{--analysis_type} . $SPACE . $analysis_type;
    }

  # shared by all Gatk modules except VariantRecalibrator and applyrecalibration
    if ($gatk_path) {

        push @commands, $gatk_path;
    }

    ## Parameters shared by all Gatk modules:
    push @commands, q{--logging_level} . $SPACE . $logging_level;

    if ($gatk_disable_auto_index_and_file_lock) {

        push @commands,
          q{--disable_auto_index_creation_and_locking_when_reading_rods};
    }

    if ($downsample_to_coverage) {

        push @commands,
          q{--downsample_to_coverage} . $SPACE . $downsample_to_coverage;
    }

    push @commands, q{--reference_sequence} . $SPACE . $referencefile_path;

    if ( @{$intervals_ref} ) {

        push @commands,
          q{--intervals} . $SPACE . join $SPACE . q{--intervals} . $SPACE,
          @{$intervals_ref};
    }

    ## Parameters specific to different Gatk modules:
    # specific for all except gatk_catvariants
    if ( $analysis_type ne q{catvariants} ) {

        if ($pedigree) {

            push @commands, q{--pedigree} . $SPACE . $pedigree;
        }

        push @commands,
          q{--pedigreeValidationType} . $SPACE . $pedigree_validation_type;
    }

# applies to gatk_selectvariants, gatk_applyrecalibration, gatk_calculategenotypeposteriors, gatk_leftalignandtrimvariants
    if ($infile_path) {

        if ( $analysis_type eq q{ApplyRecalibration} ) {

            push @commands, q{--input} . $SPACE . $infile_path;
        }
        else {

            push @commands, q{--variant} . $SPACE . $infile_path;
        }
    }

# applies to all except gatk_selectvariants, gatk_applyrecalibration, gatk_calculategenotypeposteriors, gatk_leftalignandtrimvariants
    if ( @{$infile_paths_ref} ) {

        if ( $analysis_type eq q{VariantRecalibrator} ) {

            push @commands, q{--input} . $SPACE . join q{--input_file} . $SPACE,
              @{$infile_paths_ref};
        }
        elsif ( $analysis_type eq q{VariantEval} ) {

            push @commands, q{--eval} . $SPACE . join q{--eval} . $SPACE,
              @{$infile_paths_ref};
        }
        else {

            push @commands,
              q{--variant} . $SPACE . join $SPACE . q{--variant} . $SPACE,
              @{$infile_paths_ref};
        }
    }

    # applies to all except gatk_variantrecalibrator
    if ($outfile_path) {

        if ( $analysis_type eq q{-cp} ) {

            push @commands, q{--outputFile} . $SPACE . $outfile_path;
        }
        else {

            push @commands, q{--out} . $SPACE . $outfile_path;
        }
    }

    #applies to gatk_genotypegvcfs and gatk_varianteval
    if ($dbsnp_file_path) {

        push @commands, q{--dbsnp} . $SPACE . $dbsnp_file_path;
    }

    # applies to gatk_variantrecalibrator, gatk_applyrecalibration
    if ($num_cpu_threads_per_data_thread) {

        push @commands,
            q{--num_cpu_threads_per_data_thread}
          . $SPACE
          . $num_cpu_threads_per_data_thread;
    }

    # applies to gatk_variantrecalibrator, gatk_applyrecalibration
    if ( @{$read_filters_ref} ) {

        push @commands,
          q{--read_filter} . $SPACE . join $SPACE . q{--read_filter} . $SPACE,
          @{$read_filters_ref};
    }

    # applies to gatk_variantrecalibrator, gatk_applyrecalibration
    if ( @{$static_quantized_quals_ref} ) {

        push
          @commands,
          q{--static_quantized_quals}
          . $SPACE
          . join $SPACE
          . q{--static_quantized_quals}
          . $SPACE, @{$static_quantized_quals_ref};
    }

    # applies to gatk_variantrecalibrator, gatk_applyrecalibration
    if ($base_quality_score_recalibration_file) {

        push @commands,
          q{--BQSR} . $SPACE . $base_quality_score_recalibration_file;
    }

    # applies to gatk_variantrecalibrator, gatk_applyrecalibration
    if ($recal_file_path) {

        push @commands, q{--recal_file} . $SPACE . $recal_file_path;
    }

    # applies to gatk_variantrecalibrator, gatk_applyrecalibration
    if ($tranches_file_path) {

        push @commands, q{--tranches_file} . $SPACE . $tranches_file_path;
    }

    # applies to gatk_variantrecalibrator, gatk_applyrecalibration
    if ($mode) {

        push @commands, q{--mode} . $SPACE . $mode;
    }

    # applies to gatk_variantrecalibrator, gatk_applyrecalibration
    if ($disable_indel_qual) {

        push @commands, q{--disable_indel_quals};
    }

    # specific to gatk_selectvariants, gatk_combinevariants
    if ($exclude_nonvariants) {

        push @commands, q{--excludeNonVariants};
    }

    # specific to gatk_genotypegvcfs
    if ($include_nonvariant_sites) {

        push @commands, q{--includeNonVariantSites};
    }

    # specific to gatk_selectvariants
    if ( @{$sample_names_ref} ) {

        push @commands,
          q{--sample_name} . $SPACE . join $SPACE . q{--sample_name} . $SPACE,
          @{$sample_names_ref};
    }

    # specific to gatk_catvariants
    if ($assume_sorted) {

        push @commands, q{--assumeSorted};
    }

    # specific for gatk_variantrecalibrator
    if ( @{$resources_ref} ) {

        push @commands,
          q{--resource:} . join $SPACE . q{--resource:}, @{$resources_ref};
    }

    # specific for gatk_variantrecalibrator'
    if ( @{$annotations_ref} ) {

        push @commands,
            q{--use_annotation}
          . $SPACE
          . join $SPACE
          . q{--use_annotation}
          . $SPACE, @{$annotations_ref};
    }

    # specific for gatk_variantrecalibrator
    if ($max_gaussian_level) {

        push @commands, q{--maxGaussians} . $SPACE . $max_gaussian_level;
    }

    # specific for gatk_variantrecalibrator
    if ($rscript_file_path) {

        push @commands, q{--rscript_file} . $SPACE . $rscript_file_path;
    }

    # specific for gatk_applyrecalibration
    if ($ts_filter_level) {

        push @commands, q{--ts_filter_level} . $SPACE . $ts_filter_level;
    }

    # specific for gatk_calculategenotypeposteriors
    if ($supporting_callset_file_path) {

        push @commands,
          q{--supporting} . $SPACE . $supporting_callset_file_path;
    }

    # specific to gatk_combinevariants
    if ( $analysis_type eq q{CombineVariants} ) {

        if ($prioritize_caller) {

            push @commands,
              q{--rod_priority_list} . $SPACE . $prioritize_caller;
        }
        push @commands,
          q{-genotypemergeoption} . $SPACE . $genotype_merge_option;

    }

    # specific for gatk_varianteval
    if ($indel_gold_standard_file_path) {

        push @commands,
          q{--goldStandard} . $SPACE . $indel_gold_standard_file_path;
    }

    # specific to gatk_leftalignandtrimvariants
    if ($split_multiallelics) {

        push @commands, q{--splitMultiallelics};
    }

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
