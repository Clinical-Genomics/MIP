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
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ picardtools_addorreplacereadgroups picardtools_mergesamfiles picardtools_markduplicates picardtools_gatherbamfiles picardtools_collectmultiplemetrics picardtools_collecthsmetrics };
}

## Constants
Readonly my $SPACE => q{ };

sub picardtools_addorreplacereadgroups {
## Function : Perl wrapper for writing picardtools addorreplacereadgroups recipe to $FILEHANDLE. Based on picardtools 2.5.0.
## Returns  : @commands
## Arguments: $create_index            => Create index
##          : $FILEHANDLE              => Sbatch filehandle to write to
##          : $infile_path             => Infile paths
##          : $java_jar                => Java jar
##          : $java_use_large_pages    => Use java large pages
##          : $memory_allocation       => Memory allocation for java
##          : $outfile_path            => Outfile path
##          : $readgroup_id            => Readgroup id value
##          : $readgroup_library       => Readgroup library
##          : $readgroup_platform      => Reagroup platform
##          : $readgroup_platform_unit => ID of the sequencing unit
##          : $readgroup_sample        => Sample id
##          : $stderrfile_path         => Stderrfile path
##          : $temp_directory          => Redirect tmp files to java temp

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
    my $create_index;
    my $java_use_large_pages;
    my $readgroup_id;
    my $readgroup_library;
    my $readgroup_platform;
    my $readgroup_platform_unit;

    my $tmpl = {
        create_index => {
            allow       => [qw{ true false }],
            default     => q{false},
            store       => \$create_index,
            strict_type => 1,
        },
        FILEHANDLE  => { store => \$FILEHANDLE },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        java_jar             => { strict_type => 1, store => \$java_jar, },
        java_use_large_pages => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$java_use_large_pages,
            strict_type => 1,
        },
        memory_allocation => { store => \$memory_allocation, strict_type => 1, },
        readgroup_id      => {
            default     => 1,
            defined     => 1,
            required    => 1,
            store       => \$readgroup_id,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        readgroup_library => {
            default     => q{undefined},
            defined     => 1,
            required    => 1,
            store       => \$readgroup_library,
            strict_type => 1,
        },
        readgroup_platform => {
            default     => q{Illumina},
            defined     => 1,
            required    => 1,
            store       => \$readgroup_platform,
            strict_type => 1,
        },
        readgroup_platform_unit => {
            default     => q{undefined},
            defined     => 1,
            required    => 1,
            store       => \$readgroup_platform_unit,
            strict_type => 1,
        },
        readgroup_sample => {
            defined     => 1,
            required    => 1,
            store       => \$readgroup_sample,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },

        temp_directory => { store => \$temp_directory, strict_type => 1, },
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

    ## Picardtools base args
    @commands = picardtools_base(
        {
            commands_ref => \@commands,
            create_index => $create_index,
        }
    );

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
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );
    return @commands;

}

sub picardtools_mergesamfiles {

## Function : Perl wrapper for writing picardtools mergesamfiles recipe to $FILEHANDLE. Based on picardtools 2.5.0.
## Returns  : @commands
## Arguments: $create_index           => Create index
##          : $FILEHANDLE             => Sbatch filehandle to write to
##          : $infile_paths_ref       => Infile paths {REF}
##          : $java_use_large_pages   => Use java large pages
##          : $java_jar               => Java jar
##          : $memory_allocation      => Memory allocation for java
##          : $outfile_path           => Outfile path
##          : $referencefile_path     => Genome reference file
##          : $regionsfile_path       => The regions to process {REF}
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $temp_directory         => Redirect tmp files to java temp
##          : $threading              => Create a background thread to encode, compress and write to disk the output file.

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_paths_ref;
    my $java_jar;
    my $memory_allocation;
    my $outfile_path;
    my $referencefile_path;
    my $regionsfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $temp_directory;
    my $threading;

    ## Default(s)
    my $create_index;
    my $java_use_large_pages;

    my $tmpl = {
        create_index => {
            allow       => [qw{ true false }],
            default     => q{false},
            store       => \$create_index,
            strict_type => 1,
        },
        FILEHANDLE       => { store => \$FILEHANDLE },
        infile_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infile_paths_ref,
            strict_type => 1,
        },
        java_jar             => { store => \$java_jar, strict_type => 1, },
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
        referencefile_path => {
            defined     => 1,
            required    => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        regionsfile_path => { store => \$regionsfile_path, strict_type => 1, },
        stderrfile_path  => { store => \$stderrfile_path,  strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        temp_directory => { store => \$temp_directory, strict_type => 1, },
        threading      => {
            allow       => [qw{ true false }],
            store       => \$threading,
            strict_type => 1,
        },
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

    ## Picardtools mergesamfiles
    push @commands, q{MergeSamFiles};

    ## Picardtools base args
    @commands = picardtools_base(
        {
            commands_ref       => \@commands,
            create_index       => $create_index,
            referencefile_path => $referencefile_path,
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
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub picardtools_markduplicates {

## Function : Perl wrapper for writing picardtools markduplicates recipe to $FILEHANDLE. Based on picardtools 2.5.0.
## Returns  : @commands
## Arguments: $create_index           => Create index
##          : $FILEHANDLE             => Sbatch filehandle to write to
##          : $infile_paths_ref       => Infile paths {REF}
##          : $java_jar               => Java jar
##          : $java_use_large_pages   => Use java large pages
##          : $memory_allocation      => Memory allocation for java
##          : $metrics_file           => File to write duplication metrics to
##          : $outfile_path           => Outfile path
##          : $referencefile_path     => Genome reference file
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $temp_directory         => Redirect tmp files to java temp

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $memory_allocation;
    my $metrics_file;
    my $infile_paths_ref;
    my $java_jar;
    my $outfile_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $temp_directory;

    ## Default(s)
    my $create_index;
    my $java_use_large_pages;

    my $tmpl = {
        create_index => {
            allow       => [qw{ true false }],
            default     => q{false},
            store       => \$create_index,
            strict_type => 1,
        },
        FILEHANDLE       => { store => \$FILEHANDLE },
        infile_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infile_paths_ref,
            strict_type => 1,
        },
        java_jar             => { store => \$java_jar, strict_type => 1, },
        java_use_large_pages => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$java_use_large_pages,
            strict_type => 1,
        },
        memory_allocation => { store => \$memory_allocation, strict_type => 1, },
        metrics_file      => {
            defined     => 1,
            required    => 1,
            store       => \$metrics_file,
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
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        temp_directory => { store => \$temp_directory, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands;

    ## Return java core commands
    if ($java_jar) {

        @commands = java_core(
            {
                java_use_large_pages => $java_use_large_pages,
                java_jar             => $java_jar,
                memory_allocation    => $memory_allocation,
                temp_directory       => $temp_directory,
            }
        );
    }

    ## Picardtools markduplicates
    push @commands, q{MarkDuplicates};

    ## Picardtools base args
    @commands = picardtools_base(
        {
            commands_ref       => \@commands,
            create_index       => $create_index,
            referencefile_path => $referencefile_path,
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
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub picardtools_gatherbamfiles {

## Function : Perl wrapper for writing picardtools gatherbamfiles recipe to $FILEHANDLE. Based on picardtools 2.5.0.
## Returns  : @commands
## Arguments: $create_index           => Create index
##          : $FILEHANDLE             => Sbatch filehandle to write to
##          : $infile_paths_ref       => Infile paths {REF}
##          : $java_jar               => Java jar
##          : $java_use_large_pages   => Use java large pages
##          : $memory_allocation      => Memory allocation for java
##          : $outfile_path           => Outfile path
##          : $referencefile_path     => Genome reference file
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $temp_directory         => Redirect tmp files to java temp

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_paths_ref;
    my $java_jar;
    my $memory_allocation;
    my $outfile_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $temp_directory;

    ## Default(s)
    my $create_index;
    my $java_use_large_pages;

    my $tmpl = {
        create_index => {
            allow       => [qw{ true false }],
            default     => q{false},
            store       => \$create_index,
            strict_type => 1,
        },
        FILEHANDLE       => { store => \$FILEHANDLE },
        infile_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infile_paths_ref,
            strict_type => 1,
        },
        java_jar             => { store => \$java_jar, strict_type => 1, },
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
        referencefile_path => {
            defined     => 1,
            required    => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        temp_directory => { store => \$temp_directory, strict_type => 1, },
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

    ## Picardtools gatherbamfiles
    push @commands, q{GatherBamFiles};

    ## Picardtools base args
    @commands = picardtools_base(
        {
            commands_ref       => \@commands,
            create_index       => $create_index,
            referencefile_path => $referencefile_path,
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
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub picardtools_collectmultiplemetrics {

## Function : Perl wrapper for writing picardtools collectmultiplemetrics recipe to $FILEHANDLE. Based on picardtools 2.5.0.
## Returns  : @commands
## Arguments: $create_index           => Create index
##          : $FILEHANDLE             => Sbatch filehandle to write to
##          : $infile_path            => Infile paths
##          : $java_jar               => Java jar
##          : $java_use_large_pages   => Use java large pages
##          : $memory_allocation      => Memory allocation for java
##          : $outfile_path           => Outfile path
##          : $referencefile_path     => Genome reference file
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $temp_directory         => Redirect tmp files to java temp

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $java_jar;
    my $memory_allocation;
    my $outfile_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $temp_directory;

    ## Default(s)
    my $create_index;
    my $java_use_large_pages;

    my $tmpl = {
        create_index => {
            allow       => [qw{ true false }],
            default     => q{false},
            store       => \$create_index,
            strict_type => 1,
        },
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
        memory_allocation => { store => \$memory_allocation, strict_type => 1, },
        outfile_path      => {
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
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        temp_directory => { store => \$temp_directory, strict_type => 1, },
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

    ## Picardtools collectmultiplemetrics
    push @commands, q{CollectMultipleMetrics};

    ## Picardtools base args
    @commands = picardtools_base(
        {
            commands_ref       => \@commands,
            create_index       => $create_index,
            referencefile_path => $referencefile_path,
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
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub picardtools_collecthsmetrics {

## Function : Perl wrapper for writing picardtools collecthsmetrics recipe to $FILEHANDLE. Based on picardtools 2.5.0.
## Returns  : @commands
## Arguments: $bait_interval_file_paths_ref   => Interval list file(s) that contains the locations of the baits used {REF}
##          : $create_index                   => Create index
##          : $FILEHANDLE                     => Sbatch filehandle to write to
##          : $infile_path                    => Infile paths
##          : $java_jar                       => Java jar
##          : $java_use_large_pages           => Use java large pages
##          : $memory_allocation              => Memory allocation for java
##          : $outfile_path                   => Outfile path
##          : $referencefile_path             => Genome reference file
##          : $stderrfile_path                => Stderrfile path
##          : $temp_directory                 => Redirect tmp files to java temp
##          : $target_interval_file_paths_ref => Interval list file(s) that contains the locations of the targets

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bait_interval_file_paths_ref;
    my $FILEHANDLE;
    my $infile_path;
    my $java_jar;
    my $memory_allocation;
    my $outfile_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $target_interval_file_paths_ref;
    my $temp_directory;

    ## Default(s)
    my $create_index;
    my $java_use_large_pages;

    my $tmpl = {
        bait_interval_file_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$bait_interval_file_paths_ref,
            strict_type => 1,
        },
        create_index => {
            allow       => [qw{ true false }],
            default     => q{false},
            store       => \$create_index,
            strict_type => 1,
        },
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
        memory_allocation => { store => \$memory_allocation, strict_type => 1, },
        outfile_path      => {
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
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        target_interval_file_paths_ref => {
            default     => [],
            defined     => 1,
            store       => \$target_interval_file_paths_ref,
            strict_type => 1,
        },
        temp_directory => { store => \$temp_directory, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands;

    ## Return java core commands
    if ($java_jar) {

        @commands = java_core(
            {
                java_use_large_pages => $java_use_large_pages,
                java_jar             => $java_jar,
                memory_allocation    => $memory_allocation,
                temp_directory       => $temp_directory,
            }
        );
    }

    ## Picardtools collecthsmetrics
    push @commands, q{CollectHsMetrics};

    ## Picardtools base args
    @commands = picardtools_base(
        {
            commands_ref       => \@commands,
            create_index       => $create_index,
            referencefile_path => $referencefile_path,
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
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );
    return @commands;
}

1;
