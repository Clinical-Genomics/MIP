package MIP::Program::Alignment::Picardtools;

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

## CPANM
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };
use MIP::Language::Java qw{java_core};
use MIP::Program::Base::Picardtools qw{ picardtools_base};

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ picardtools_addorreplacereadgroups picardtools_mergesamfiles  picardtools_markduplicates picardtools_gatherbamfiles picardtools_collectmultiplemetrics picardtools_collecthsmetrics };
}

## Constants
Readonly my $SPACE => q{ };

sub picardtools_addorreplacereadgroups {
## Function : Perl wrapper for writing picardtools addorreplacereadgroups recipe to $FILEHANDLE. Based on picardtools 2.5.0.
## Returns  : @commands
## Arguments: $FILEHANDLE                     => Sbatch filehandle to write to
##          : $infile_path                    => Infile paths
##          : $java_jar                       => Java jar
##          : $java_use_large_pages           => Use java large pages
##          : $memory_allocation              => Memory allocation for java
##          : $outfile_path                   => Outfile path
##          : $readgroup_id                   => Readgroup id value
##          : $readgroup_library              => Readgroup library
##          : $readgroup_platform             => Reagroup platform
##          : $readgroup_platform_unit        => ID of the sequencing unit
##          : $readgroup_sample               => Sample id
##          : $stderrfile_path                => Stderrfile path
##          : $temp_directory                 => Redirect tmp files to java temp

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $java_jar;
    my $memory_allocation;
    my $outfile_path;
    my $readgroup_sample;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $temp_directory;

    ## Default(s)
    my $java_use_large_pages;
    my $readgroup_id;
    my $readgroup_library;
    my $readgroup_platform;
    my $readgroup_platform_unit;

    my $tmpl = {
        FILEHANDLE  => { store => \$FILEHANDLE },
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
        memory_allocation => { strict_type => 1, store => \$memory_allocation },
        readgroup_id      => {
            default     => 1,
            defined     => 1,
            required    => 1,
            strict_type => 1,
            store       => \$readgroup_id,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            strict_type => 1,
            store       => \$outfile_path,
        },
        readgroup_library => {
            default     => q{undefined},
            defined     => 1,
            required    => 1,
            strict_type => 1,
            store       => \$readgroup_library,
        },
        readgroup_platform => {
            default     => q{Illumina},
            defined     => 1,
            required    => 1,
            strict_type => 1,
            store       => \$readgroup_platform,
        },
        readgroup_platform_unit => {
            default     => q{undefined},
            defined     => 1,
            required    => 1,
            strict_type => 1,
            store       => \$readgroup_platform_unit,
        },
        readgroup_sample => {
            defined     => 1,
            required    => 1,
            strict_type => 1,
            store       => \$readgroup_sample,
        },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },

        temp_directory => { strict_type => 1, store => \$temp_directory },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands;

    ## Return java core commands
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

    ## Picardtools addorreaplacereadgroups
    push @commands, q{AddOrReplaceReadGroups};

    ## Infile
    push @commands, q{INPUT=} . $infile_path;

    ## Output
    push @commands, q{OUTPUT=} . $outfile_path;

    ## Readgroup id
    push @commands, q{RGID=} . $readgroup_id;

    ## Readgroup library
    push @commands, q{RGLB=} . $readgroup_library;

    ## Readgroup platform
    push @commands, q{RGPL=} . $readgroup_platform;

    ## Readgroup platform unit
    push @commands, q{RGPU=} . $readgroup_platform_unit;

    ## Readgroup sample
    push @commands, q{RGSM=} . $readgroup_sample;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
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

sub picardtools_mergesamfiles {

## Function : Perl wrapper for writing picardtools mergesamfiles recipe to $FILEHANDLE. Based on picardtools 2.5.0.
## Returns  : @commands
## Arguments: $infile_paths_ref       => Infile paths {REF}
##          : $outfile_path           => Outfile path
##          : $referencefile_path     => Genome reference file
##          : $regionsfile_path       => The regions to process {REF}
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $FILEHANDLE             => Sbatch filehandle to write to
##          : $memory_allocation      => Memory allocation for java
##          : $temp_directory         => Redirect tmp files to java temp
##          : $java_use_large_pages   => Use java large pages
##          : $java_jar               => Java jar
##          : $create_index           => Create index
##          : $threading              => Create a background thread to encode, compress and write to disk the output file.

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_paths_ref;
    my $outfile_path;
    my $referencefile_path;
    my $regionsfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;
    my $memory_allocation;
    my $temp_directory;
    my $java_jar;
    my $threading;

    ## Default(s)
    my $java_use_large_pages;
    my $create_index;

    my $tmpl = {
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
        regionsfile_path => { strict_type => 1, store => \$regionsfile_path },
        stderrfile_path  => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        FILEHANDLE        => { store       => \$FILEHANDLE },
        memory_allocation => { strict_type => 1, store => \$memory_allocation },
        temp_directory    => { strict_type => 1, store => \$temp_directory },
        java_jar          => { strict_type => 1, store => \$java_jar },
        java_use_large_pages => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$java_use_large_pages
        },
        create_index => {
            default     => q{false},
            allow       => [qw{ true false }],
            strict_type => 1,
            store       => \$create_index
        },
        threading => {
            allow       => [qw{ true false }],
            strict_type => 1,
            store       => \$threading
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands;

    ## Return java core commands
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

    ## Picardtools mergesamfiles
    push @commands, q{MergeSamFiles};

    ## Picardtools base args
    @commands = picardtools_base(
        {
            commands_ref       => \@commands,
            referencefile_path => $referencefile_path,
            create_index       => $create_index,
        }
    );

    if ($threading) {

        ## Create a background thread to encode, compress and write to disk the output file
        push @commands, q{USE_THREADING=} . $threading;
    }

    ## Infile
    push @commands, q{INPUT=} . join $SPACE . q{INPUT=}, @{$infile_paths_ref};

    ## Output
    push @commands, q{OUTPUT=} . $outfile_path;

    # Limit output to regions
    if ($regionsfile_path) {

        push @commands, q{INTERVALS=} . $regionsfile_path;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
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

sub picardtools_markduplicates {

## Function : Perl wrapper for writing picardtools markduplicates recipe to $FILEHANDLE. Based on picardtools 2.5.0.
## Returns  : @commands
## Arguments: $infile_paths_ref       => Infile paths {REF}
##          : $outfile_path           => Outfile path
##          : $referencefile_path     => Genome reference file
##          : $metrics_file           => File to write duplication metrics to
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $FILEHANDLE             => Sbatch filehandle to write to
##          : $memory_allocation      => Memory allocation for java
##          : $temp_directory         => Redirect tmp files to java temp
##          : $java_use_large_pages   => Use java large pages
##          : $java_jar               => Java jar
##          : $create_index           => Create index

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_paths_ref;
    my $outfile_path;
    my $referencefile_path;
    my $metrics_file;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;
    my $memory_allocation;
    my $temp_directory;
    my $java_jar;

    ## Default(s)
    my $java_use_large_pages;
    my $create_index;

    my $tmpl = {
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
        metrics_file => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$metrics_file
        },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        FILEHANDLE        => { store       => \$FILEHANDLE },
        memory_allocation => { strict_type => 1, store => \$memory_allocation },
        temp_directory    => { strict_type => 1, store => \$temp_directory },
        java_jar          => { strict_type => 1, store => \$java_jar },
        java_use_large_pages => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$java_use_large_pages
        },
        create_index => {
            default     => q{false},
            allow       => [qw{ true false }],
            strict_type => 1,
            store       => \$create_index
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands;

    ## Return java core commands
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

    ## Picardtools markduplicates
    push @commands, q{MarkDuplicates};

    ## Picardtools base args
    @commands = picardtools_base(
        {
            commands_ref       => \@commands,
            referencefile_path => $referencefile_path,
            create_index       => $create_index,
        }
    );

    if ($metrics_file) {

        # File to write duplication metrics to
        push @commands, q{METRICS_FILE=} . $metrics_file;
    }

    ## Infile
    push @commands, q{INPUT=} . join $SPACE . q{INPUT=}, @{$infile_paths_ref};

    ## Output
    push @commands, q{OUTPUT=} . $outfile_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
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

sub picardtools_gatherbamfiles {

## Function : Perl wrapper for writing picardtools gatherbamfiles recipe to $FILEHANDLE. Based on picardtools 2.5.0.
## Returns  : @commands
## Arguments: $infile_paths_ref       => Infile paths {REF}
##          : $outfile_path           => Outfile path
##          : $referencefile_path     => Genome reference file
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $FILEHANDLE             => Sbatch filehandle to write to
##          : $memory_allocation      => Memory allocation for java
##          : $temp_directory         => Redirect tmp files to java temp
##          : $java_use_large_pages   => Use java large pages
##          : $java_jar               => Java jar
##          : $create_index           => Create index

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_paths_ref;
    my $outfile_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;
    my $memory_allocation;
    my $temp_directory;
    my $java_jar;

    ## Default(s)
    my $java_use_large_pages;
    my $create_index;

    my $tmpl = {
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
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        FILEHANDLE        => { store       => \$FILEHANDLE },
        memory_allocation => { strict_type => 1, store => \$memory_allocation },
        temp_directory    => { strict_type => 1, store => \$temp_directory },
        java_jar          => { strict_type => 1, store => \$java_jar },
        java_use_large_pages => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$java_use_large_pages
        },
        create_index => {
            default     => q{false},
            allow       => [qw{ true false }],
            strict_type => 1,
            store       => \$create_index
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands;

    ## Return java core commands
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

    ## Picardtools gatherbamfiles
    push @commands, q{GatherBamFiles};

    ## Picardtools base args
    @commands = picardtools_base(
        {
            commands_ref       => \@commands,
            referencefile_path => $referencefile_path,
            create_index       => $create_index,
        }
    );

    ## Infile
    push @commands, q{INPUT=} . join $SPACE . q{INPUT=}, @{$infile_paths_ref};

    ## Output
    push @commands, q{OUTPUT=} . $outfile_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
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

sub picardtools_collectmultiplemetrics {

## Function : Perl wrapper for writing picardtools collectmultiplemetrics recipe to $FILEHANDLE. Based on picardtools 2.5.0.
## Returns  : @commands
## Arguments: $infile_path            => Infile paths
##          : $outfile_path           => Outfile path
##          : $referencefile_path     => Genome reference file
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $FILEHANDLE             => Sbatch filehandle to write to
##          : $memory_allocation      => Memory allocation for java
##          : $temp_directory         => Redirect tmp files to java temp
##          : $java_use_large_pages   => Use java large pages
##          : $java_jar               => Java jar
##          : $create_index           => Create index

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $outfile_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;
    my $memory_allocation;
    my $temp_directory;
    my $java_jar;

    ## Default(s)
    my $java_use_large_pages;
    my $create_index;

    my $tmpl = {
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
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        FILEHANDLE        => { store       => \$FILEHANDLE },
        memory_allocation => { strict_type => 1, store => \$memory_allocation },
        temp_directory    => { strict_type => 1, store => \$temp_directory },
        java_jar          => { strict_type => 1, store => \$java_jar },
        java_use_large_pages => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$java_use_large_pages
        },
        create_index => {
            default     => q{false},
            allow       => [qw{ true false }],
            strict_type => 1,
            store       => \$create_index
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands;

    ## Return java core commands
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

    ## Picardtools collectmultiplemetrics
    push @commands, q{CollectMultipleMetrics};

    ## Picardtools base args
    @commands = picardtools_base(
        {
            commands_ref       => \@commands,
            referencefile_path => $referencefile_path,
            create_index       => $create_index,
        }
    );

    ## Infile
    push @commands, q{INPUT=} . $infile_path;

    ## Output
    push @commands, q{OUTPUT=} . $outfile_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
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

sub picardtools_collecthsmetrics {

## Function : Perl wrapper for writing picardtools collecthsmetrics recipe to $FILEHANDLE. Based on picardtools 2.5.0.
## Returns  : @commands
## Arguments: $bait_interval_file_paths_ref   => Interval list file(s) that contains the locations of the baits used {REF}
##          : $target_interval_file_paths_ref => Interval list file(s) that contains the locations of the targets
##          : $infile_path                    => Infile paths
##          : $outfile_path                   => Outfile path
##          : $referencefile_path             => Genome reference file
##          : $stderrfile_path                => Stderrfile path
##          : $FILEHANDLE                     => Sbatch filehandle to write to
##          : $memory_allocation              => Memory allocation for java
##          : $temp_directory                 => Redirect tmp files to java temp
##          : $java_use_large_pages           => Use java large pages
##          : $java_jar                       => Java jar
##          : $create_index                   => Create index

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bait_interval_file_paths_ref;
    my $target_interval_file_paths_ref;
    my $infile_path;
    my $outfile_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;
    my $memory_allocation;
    my $temp_directory;
    my $java_jar;

    ## Default(s)
    my $java_use_large_pages;
    my $create_index;

    my $tmpl = {
        bait_interval_file_paths_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$bait_interval_file_paths_ref
        },
        target_interval_file_paths_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$target_interval_file_paths_ref
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
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        FILEHANDLE        => { store       => \$FILEHANDLE },
        memory_allocation => { strict_type => 1, store => \$memory_allocation },
        temp_directory    => { strict_type => 1, store => \$temp_directory },
        java_jar          => { strict_type => 1, store => \$java_jar },
        java_use_large_pages => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$java_use_large_pages
        },
        create_index => {
            default     => q{false},
            allow       => [qw{ true false }],
            strict_type => 1,
            store       => \$create_index
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands;

    ## Return java core commands
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

    ## Picardtools collecthsmetrics
    push @commands, q{CollectHsMetrics};

    ## Picardtools base args
    @commands = picardtools_base(
        {
            commands_ref       => \@commands,
            referencefile_path => $referencefile_path,
            create_index       => $create_index,
        }
    );

    ## Baits
    push @commands,
      q{BAIT_INTERVALS=} . join $SPACE . q{BAIT_INTERVALS=},
      @{$bait_interval_file_paths_ref};

    ## Targets
    push @commands,
      q{TARGET_INTERVALS=} . join $SPACE . q{TARGET_INTERVALS=},
      @{$target_interval_file_paths_ref};

    ## Infile
    push @commands, q{INPUT=} . $infile_path;

    ## Output
    push @commands, q{OUTPUT=} . $outfile_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
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
