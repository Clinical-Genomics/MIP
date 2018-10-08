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
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ gatk_base gatk_common_options gatk_java_options };
}

## Constants
Readonly my $DOUBLE_QUOTE => q{"};
Readonly my $ESCAPE       => q{\\};
Readonly my $SPACE        => q{ };

sub gatk_base {

## Function : Perl wrapper for Gatk base parameters. Based on Gatk 3.7
## Returns  : @commands
## Arguments: $commands_ref                          => List of commands added earlier
##          : $intervals_ref                         => One or more genomic intervals over which to operate {REF}
##          : $read_filters_ref                      => Filters to apply to reads before analysis {REF}
##          : $static_quantized_quals_ref            => Use static quantized quality scores [ref]
##          : $analysis_type                         => Analysis type
##          : $referencefile_path                    => Reference sequence file
##          : $pedigree                              => Pedigree files
##          : $pedigree_validation_type              => Validation strictness for pedigree
##          : $downsample_to_coverage                => Target coverage threshold for downsampling to coverage
##          : $gatk_disable_auto_index_and_file_lock => Disable both auto-generation of index files and index file locking
##          : $base_quality_score_recalibration_file => Base quality score recalibration file
##          : $disable_indel_qual                    => Disable indel quality
##          : $num_cpu_threads_per_data_thread       => Number of CPU threads to allocate per data thread
##          : $logging_level                         => Logging level
##          : $FILEHANDLE                            => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $commands_ref;
    my $analysis_type;
    my $intervals_ref;
    my $read_filters_ref;
    my $static_quantized_quals_ref;
    my $referencefile_path;
    my $pedigree;
    my $pedigree_validation_type;
    my $downsample_to_coverage;
    my $base_quality_score_recalibration_file;
    my $disable_indel_qual;
    my $FILEHANDLE;

    ## Default(s)
    my $num_cpu_threads_per_data_thread;
    my $logging_level;
    my $gatk_disable_auto_index_and_file_lock;

    my $tmpl = {
        commands_ref =>
          { default => [], strict_type => 1, store => \$commands_ref },
        analysis_type => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$analysis_type
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
            strict_type => 1,
            store       => \$downsample_to_coverage
        },
        base_quality_score_recalibration_file => {
            strict_type => 1,
            store       => \$base_quality_score_recalibration_file
        },
        disable_indel_qual => {
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$disable_indel_qual
        },
        num_cpu_threads_per_data_thread => {
            default     => 0,
            allow       => [ undef, qr/ ^\d+$ /sxm ],
            strict_type => 1,
            store       => \$num_cpu_threads_per_data_thread
        },
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
        FILEHANDLE => { store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = @{$commands_ref};

    push @commands, q{--analysis_type} . $SPACE . $analysis_type;

    push @commands, q{--logging_level} . $SPACE . $logging_level;

    if ($pedigree_validation_type) {

        push @commands,
          q{--pedigreeValidationType} . $SPACE . $pedigree_validation_type;
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

        push @commands,
          q{--downsample_to_coverage} . $SPACE . $downsample_to_coverage;
    }

    if ($gatk_disable_auto_index_and_file_lock) {

        push @commands,
          q{--disable_auto_index_creation_and_locking_when_reading_rods};
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

        push @commands,
          q{--BQSR} . $SPACE . $base_quality_score_recalibration_file;
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
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

    ## Add to command array
    push @{$commands_ref}, $java_command;

    return $commands_ref;
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
            strict_type => 1,
            store       => \$pedigree
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

    ## Add intervals
    if ( @{$intervals_ref} ) {
        push @{$commands_ref},
          q{--intervals} . $SPACE . join $SPACE . q{--intervals} . $SPACE,
          @{$intervals_ref};
    }

    ## Add Pedigree
    if ($pedigree) {
        push @{$commands_ref}, q{--pedigree} . $SPACE . $pedigree;
    }

    ## Add read filters
    if ( @{$read_filters_ref} ) {
        push @{$commands_ref},
          q{--read-filter} . $SPACE . join $SPACE . q{--read-filter} . $SPACE,
          @{$read_filters_ref};
    }

    ## Add path to reference
    if ($referencefile_path) {
        push @{$commands_ref}, q{--reference} . $SPACE . $referencefile_path;
    }

    ## Add path to temporary directory
    if ($temp_directory) {
        push @{$commands_ref}, q{--tmp-dir} . $SPACE . $temp_directory;
    }

    ## Add verbosity level
    if ($verbosity) {
        push @{$commands_ref}, q{--verbosity} . $SPACE . $verbosity;
    }

    return @{$commands_ref};
}
1;
