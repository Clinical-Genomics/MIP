package MIP::Program::Picardtools;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $SPACE };
use MIP::Language::Java qw{java_core};
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ picardtools_addorreplacereadgroups
      picardtools_base
      picardtools_createsequencedictionary
      picardtools_collecthsmetrics
      picardtools_collectmultiplemetrics
      picardtools_gatherbamfiles
      picardtools_intervallisttools
      picardtools_markduplicates
      picardtools_mergesamfiles
      picardtools_sortvcf
      sort_vcf
    };
}

sub picardtools_addorreplacereadgroups {

## Function : Perl wrapper for writing picardtools addorreplacereadgroups recipe to $FILEHANDLE. Based on picardtools 2.20.7.
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

sub picardtools_base {

## Function : Perl wrapper for picardtools base. Based on Picardtools v2.9.2-SNAPSHOT
## Returns  : @commands
## Arguments: $commands_ref       => List of commands added earlier
##          : $FILEHANDLE         => Filehandle to write to
##          : $create_index       => Create a BAM index when writing a coordinate-sorted BAM file
##          : $referencefile_path => Genome reference file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $commands_ref;
    my $referencefile_path;
    my $FILEHANDLE;

    ## Default(s)
    my $create_index;

    my $tmpl = {
        commands_ref => {
            default     => [],
            store       => \$commands_ref,
            strict_type => 1,
        },
        create_index => {
            allow       => [qw{ true false }],
            default     => q{false},
            store       => \$create_index,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        FILEHANDLE => { store => \$FILEHANDLE, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = @{$commands_ref};

    if ( $create_index ne q{false} ) {

        push @commands, q{CREATE_INDEX=} . $create_index;
    }
    if ($referencefile_path) {

        push @commands, q{R=} . $referencefile_path;
    }
    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub picardtools_createsequencedictionary {

## Function : Perl wrapper for writing picardtools createsequencedictionary recipe to $FILEHANDLE. Based on picardtools 2.5.0.
## Returns  : @commands
## Arguments: $referencefile_path     => Genome reference file
##          : $outfile_path           => Outfile path
##          : $stdoutfile_path        => Stdoutfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $FILEHANDLE             => Filehandle to write to
##          : $memory_allocation      => Memory allocation for java
##          : $temp_directory         => Redirect tmp files to java temp
##          : $java_use_large_pages   => Use java large pages
##          : $java_jar               => Java jar
##          : $create_index           => Create index

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $referencefile_path;
    my $outfile_path;
    my $stdoutfile_path;
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
        referencefile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$referencefile_path
        },
        outfile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfile_path
        },
        stdoutfile_path => {
            strict_type => 1,
            store       => \$stdoutfile_path,
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path,
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        memory_allocation    => { strict_type => 1, store => \$memory_allocation },
        temp_directory       => { strict_type => 1, store => \$temp_directory },
        java_jar             => { strict_type => 1, store => \$java_jar },
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

    ## Picardtools mergesamfiles
    push @commands, q{CreateSequenceDictionary};

    ## Picardtools base args
    @commands = picardtools_base(
        {
            commands_ref       => \@commands,
            referencefile_path => $referencefile_path,
            create_index       => $create_index,
        }
    );

    ## Output
    push @commands, q{OUTPUT=} . $outfile_path;

    push @commands,
      unix_standard_streams(
        {
            stdoutfile_path        => $stdoutfile_path,
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

sub picardtools_intervallisttools {

## Function : Perl wrapper for writing picardtools intervallisttools recipe to $FILEHANDLE. Based on picardtools 2.5.0.
##Returns  : @commands
## Arguments: $infile_paths_ref     => Infile paths {REF}
##          : $outfile_path         => Outfile path
##          : $referencefile_path   => Genome reference file
##          : $FILEHANDLE           => Sbatch filehandle to write to
##          : $stderrfile_path      => Stderrfile path
##          : $memory_allocation    => Memory allocation for java
##          : $temp_directory       => Redirect tmp files to java temp
##          : $java_use_large_pages => Use java large pages
##          : $java_jar             => Java jar
##          : $create_index         => Create index
##          : $padding              => The amount to pad each end of the intervals by before other operations are undertaken

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_paths_ref;
    my $outfile_path;
    my $referencefile_path;
    my $FILEHANDLE;
    my $stderrfile_path;
    my $memory_allocation;
    my $temp_directory;
    my $java_jar;

    ## Default(s)
    my $java_use_large_pages;
    my $create_index;
    my $padding;

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
        FILEHANDLE        => { store       => \$FILEHANDLE },
        stderrfile_path   => { strict_type => 1, store => \$stderrfile_path },
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
        padding => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$padding
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

    ## Picardtools intervallisttools
    push @commands, q{IntervalListTools};

    ## Picardtools base args
    @commands = picardtools_base(
        {
            commands_ref       => \@commands,
            referencefile_path => $referencefile_path,
            create_index       => $create_index,
        }
    );

    ##Options
    if ($padding) {

        push @commands, q{PADDING=} . $padding;
    }

    ## Infile
    push @commands, q{INPUT=} . join $SPACE . q{INPUT=}, @{$infile_paths_ref};

    ## Output
    push @commands, q{OUTPUT=} . $outfile_path;

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

sub picardtools_mergesamfiles {

## Function : Perl wrapper for writing picardtools mergesamfiles recipe to $FILEHANDLE. Based on picardtools 2.20.7.
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

## Function : Perl wrapper for writing picardtools markduplicates recipe to $FILEHANDLE. Based on picardtools 2.20.7.
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

## Function : Perl wrapper for writing picardtools gatherbamfiles recipe to $FILEHANDLE. Based on picardtools 2.20.7.
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

## Function : Perl wrapper for writing picardtools collectmultiplemetrics recipe to $FILEHANDLE. Based on picardtools 2.20.7.
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

## Function : Perl wrapper for writing picardtools collecthsmetrics recipe to $FILEHANDLE. Based on picardtools 2.20.7.
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

sub picardtools_sortvcf {

## Function : Perl wrapper for writing picardtools sortvcf recipe to $FILEHANDLE. Based on picardtools 2.14.1.
## Returns  : @commands
## Arguments: $create_index         => Create index
##          : $FILEHANDLE           => Sbatch filehandle to write to
##          : $java_jar             => Java jar
##          : $java_use_large_pages => Use java large pages
##          : $infile_paths_ref     => Infile paths {REF}
##          : $memory_allocation    => Memory allocation for java
##          : $outfile_path         => Outfile path
##          : $referencefile_path   => Genome reference file
##          : $sequence_dictionary  => Sequence dictionary
##          : $stderrfile_path      => Stderrfile path
##          : $temp_directory       => Redirect tmp files to java temp

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_paths_ref;
    my $java_jar;
    my $memory_allocation;
    my $outfile_path;
    my $referencefile_path;
    my $sequence_dictionary;
    my $stderrfile_path;
    my $temp_directory;

    ## Default(s)
    my $create_index;
    my $java_use_large_pages;

    my $tmpl = {
        create_index => {
            default     => q{false},
            allow       => [qw{ true false }],
            strict_type => 1,
            store       => \$create_index
        },
        FILEHANDLE       => { store => \$FILEHANDLE },
        infile_paths_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$infile_paths_ref
        },
        java_jar             => { strict_type => 1, store => \$java_jar },
        java_use_large_pages => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$java_use_large_pages
        },
        memory_allocation => { strict_type => 1, store => \$memory_allocation },
        outfile_path      => {
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
        sequence_dictionary => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sequence_dictionary
        },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        temp_directory  => { strict_type => 1, store => \$temp_directory },
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
    push @commands, q{SortVcf};

    ## Picardtools base args
    @commands = picardtools_base(
        {
            commands_ref       => \@commands,
            referencefile_path => $referencefile_path,
            create_index       => $create_index,
        }
    );

    push @commands, q{SEQUENCE_DICTIONARY=} . $sequence_dictionary;

    ## Infile
    push @commands, q{INPUT=} . join $SPACE . q{INPUT=}, @{$infile_paths_ref};

    ## Output
    push @commands, q{OUTPUT=} . $outfile_path;

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

sub sort_vcf {

## Function : Writes sbatch code to supplied filehandle to sort variants in vcf format.
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $FILEHANDLE            => SBATCH script FILEHANDLE to print to
##          : $infile_paths_ref      => Infiles to sort {REF}
##          : $outfile               => The sorted outfile
##          : $sequence_dict_file    => Human reference sequence dict file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $FILEHANDLE;
    my $infile_paths_ref;
    my $outfile;
    my $sequence_dict_file;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        FILEHANDLE       => { required => 1, defined => 1, store => \$FILEHANDLE, },
        infile_paths_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$infile_paths_ref
        },
        outfile => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfile,
        },
        sequence_dict_file => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sequence_dict_file
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    say {$FILEHANDLE} q{## Picard SortVcf};

    ## Writes java core commands to filehandle.
    picardtools_sortvcf(
        {
            FILEHANDLE           => $FILEHANDLE,
            infile_paths_ref     => \@{$infile_paths_ref},
            java_use_large_pages => $active_parameter_href->{java_use_large_pages},
            java_jar =>
              catfile( $active_parameter_href->{picardtools_path}, q{picard.jar} ),
            memory_allocation   => q{Xmx4g},
            outfile_path        => $outfile,
            referencefile_path  => $active_parameter_href->{human_genome_reference},
            sequence_dictionary => $sequence_dict_file,
            temp_directory      => $active_parameter_href->{temp_directory},
        }
    );
    return;
}

1;
