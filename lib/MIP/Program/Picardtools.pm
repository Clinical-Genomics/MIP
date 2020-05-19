package MIP::Program::Picardtools;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
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
use MIP::Language::Java qw{ java_core };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.07;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      picardtools_addorreplacereadgroups
      picardtools_base
      picardtools_createsequencedictionary
      picardtools_collecthsmetrics
      picardtools_collectmultiplemetrics
      picardtools_collectrnaseqmetrics
      picardtools_gatherbamfiles
      picardtools_intervallisttools
      picardtools_markduplicates
      picardtools_mergesamfiles
      picardtools_sortvcf
      sort_vcf
    };
}

sub picardtools_addorreplacereadgroups {

## Function : Perl wrapper for writing picardtools addorreplacereadgroups recipe to $filehandle. Based on picardtools 2.20.7.
## Returns  : @commands
## Arguments: $create_index            => Create index
##          : $filehandle              => Sbatch filehandle to write to
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
    my $filehandle;
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
        filehandle  => { store => \$filehandle },
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
                java_jar                  => $java_jar,
                java_use_large_pages      => $java_use_large_pages,
                memory_allocation         => $memory_allocation,
                picard_use_barclay_parser => 1,
                temp_directory            => $temp_directory,
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
    push @commands, q{-INPUT} . $SPACE . $infile_path;

    ## Output
    push @commands, q{-OUTPUT} . $SPACE . $outfile_path;

    ## Readgroup id
    push @commands, q{-RGID} . $SPACE . $readgroup_id;

    ## Readgroup library
    push @commands, q{-RGLB} . $SPACE . $readgroup_library;

    ## Readgroup platform
    push @commands, q{-RGPL} . $SPACE . $readgroup_platform;

    ## Readgroup platform unit
    push @commands, q{-RGPU} . $SPACE . $readgroup_platform_unit;

    ## Readgroup sample
    push @commands, q{-RGSM} . $SPACE . $readgroup_sample;

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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;

}

sub picardtools_base {

## Function : Perl wrapper for picardtools base. Based on Picardtools v2.20.7
## Returns  : @commands
## Arguments: $create_index       => Create a BAM index when writing a coordinate-sorted BAM file
##          : $commands_ref       => List of commands added earlier
##          : $filehandle         => Filehandle to write to
##          : $referencefile_path => Genome reference file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $commands_ref;
    my $filehandle;
    my $referencefile_path;

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
        filehandle         => { store => \$filehandle, },
        referencefile_path => {
            defined     => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = @{$commands_ref};

    unshift @commands, q{picard};

    if ( $create_index ne q{false} ) {

        push @commands, q{-CREATE_INDEX} . $SPACE . $create_index;
    }
    if ($referencefile_path) {

        push @commands, q{-R} . $SPACE . $referencefile_path;
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

sub picardtools_createsequencedictionary {

## Function : Perl wrapper for writing picardtools createsequencedictionary recipe to $filehandle. Based on picardtools 2.20.7.
## Returns  : @commands
## Arguments: $create_index           => Create index
##          : $filehandle             => Filehandle to write to
##          : $java_jar               => Java jar
##          : $java_use_large_pages   => Use java large pages
##          : $memory_allocation      => Memory allocation for java
##          : $outfile_path           => Outfile path
##          : $referencefile_path     => Genome reference file
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $temp_directory         => Redirect tmp files to java temp

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $java_jar;
    my $memory_allocation;
    my $outfile_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $temp_directory;

    ## Default(s)
    my $create_index;
    my $java_use_large_pages;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        create_index => {
            allow       => [qw{ true false }],
            default     => q{false},
            store       => \$create_index,
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
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        stderrfile_path_append => {
            store       => \$stderrfile_path_append,
            strict_type => 1,
        },
        stdoutfile_path => {
            store       => \$stdoutfile_path,
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
                java_jar                  => $java_jar,
                java_use_large_pages      => $java_use_large_pages,
                memory_allocation         => $memory_allocation,
                picard_use_barclay_parser => 1,
                temp_directory            => $temp_directory,
            }
        );
    }

    ## Picardtools mergesamfiles
    push @commands, q{CreateSequenceDictionary};

    ## Picardtools base args
    @commands = picardtools_base(
        {
            commands_ref       => \@commands,
            create_index       => $create_index,
            referencefile_path => $referencefile_path,
        }
    );

    ## Output
    push @commands, q{-OUTPUT} . $SPACE . $outfile_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdoutfile_path        => $stdoutfile_path,
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

sub picardtools_intervallisttools {

## Function : Perl wrapper for writing picardtools intervallisttools recipe to $filehandle. Based on picardtools 2.20.7.
##Returns  : @commands
## Arguments: $create_index         => Create index
##          : $filehandle           => Sbatch filehandle to write to
##          : $infile_paths_ref     => Infile paths {REF}
##          : $java_jar             => Java jar
##          : $java_use_large_pages => Use java large pages
##          : $memory_allocation    => Memory allocation for java
##          : $outfile_path         => Outfile path
##          : $padding              => The amount to pad each end of the intervals by before other operations are undertaken
##          : $referencefile_path   => Genome reference file
##          : $stderrfile_path      => Stderrfile path
##          : $temp_directory       => Redirect tmp files to java temp

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_paths_ref;
    my $java_jar;
    my $memory_allocation;
    my $outfile_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $temp_directory;

    ## Default(s)
    my $create_index;
    my $java_use_large_pages;
    my $padding;

    my $tmpl = {
        create_index => {
            allow       => [qw{ true false }],
            default     => q{false},
            store       => \$create_index,
            strict_type => 1,
        },
        filehandle       => { store => \$filehandle, },
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
        padding => {
            allow       => qr/ ^\d+$ /xsm,
            default     => 0,
            store       => \$padding,
            strict_type => 1,
        },
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

    # Stores commands depending on input parameters
    my @commands;

    ## Return java core commands
    if ($java_jar) {

        @commands = java_core(
            {
                java_jar                  => $java_jar,
                java_use_large_pages      => $java_use_large_pages,
                memory_allocation         => $memory_allocation,
                picard_use_barclay_parser => 1,
                temp_directory            => $temp_directory,
            }
        );
    }

    ## Picardtools intervallisttools
    push @commands, q{IntervalListTools};

    ## Picardtools base args
    @commands = picardtools_base(
        {
            commands_ref       => \@commands,
            create_index       => $create_index,
            referencefile_path => $referencefile_path,
        }
    );

    ##Options
    if ($padding) {

        push @commands, q{-PADDING} . $SPACE . $padding;
    }

    ## Infile
    push @commands, q{-INPUT} . $SPACE . join $SPACE . q{-INPUT} . $SPACE,
      @{$infile_paths_ref};

    ## Output
    push @commands, q{-OUTPUT} . $SPACE . $outfile_path;

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

sub picardtools_mergesamfiles {

## Function : Perl wrapper for writing picardtools mergesamfiles recipe to $filehandle. Based on picardtools 2.20.7.
## Returns  : @commands
## Arguments: $create_index           => Create index
##          : $filehandle             => Sbatch filehandle to write to
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
    my $filehandle;
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
        filehandle       => { store => \$filehandle },
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
                java_jar                  => $java_jar,
                java_use_large_pages      => $java_use_large_pages,
                memory_allocation         => $memory_allocation,
                picard_use_barclay_parser => 1,
                temp_directory            => $temp_directory,
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
        push @commands, q{-USE_THREADING} . $SPACE . $threading;
    }

    ## Infile
    push @commands, q{-INPUT} . $SPACE . join $SPACE . q{-INPUT} . $SPACE,
      @{$infile_paths_ref};

    ## Output
    push @commands, q{-OUTPUT} . $SPACE . $outfile_path;

    # Limit output to regions
    if ($regionsfile_path) {

        push @commands, q{-INTERVALS} . $SPACE . $regionsfile_path;
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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub picardtools_markduplicates {

## Function : Perl wrapper for writing picardtools markduplicates recipe to $filehandle. Based on picardtools 2.20.7.
## Returns  : @commands
## Arguments: $create_index               => Create index
##          : $filehandle                 => Sbatch filehandle to write to
##          : $infile_paths_ref           => Infile paths {REF}
##          : $java_jar                   => Java jar
##          : $java_use_large_pages       => Use java large pages
##          : $memory_allocation          => Memory allocation for java
##          : $metrics_file               => File to write duplication metrics to
##          : $optical_duplicate_distance => Max distance between optical duplicates
##          : $outfile_path               => Outfile path
##          : $referencefile_path         => Genome reference file
##          : $stderrfile_path            => Stderrfile path
##          : $stderrfile_path_append     => Append stderr info to file path
##          : $temp_directory             => Redirect tmp files to java temp

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $memory_allocation;
    my $metrics_file;
    my $infile_paths_ref;
    my $java_jar;
    my $optical_duplicate_distance;
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
        filehandle       => { store => \$filehandle, },
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
        optical_duplicate_distance => {
            allow       => [ undef, qr/ \A \d+ \z /xms ],
            store       => \$optical_duplicate_distance,
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
                java_use_large_pages      => $java_use_large_pages,
                java_jar                  => $java_jar,
                memory_allocation         => $memory_allocation,
                picard_use_barclay_parser => 1,
                temp_directory            => $temp_directory,
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

        push @commands, q{-METRICS_FILE} . $SPACE . $metrics_file;
    }
    if ($optical_duplicate_distance) {

        push @commands,
          q{-OPTICAL_DUPLICATE_PIXEL_DISTANCE} . $SPACE . $optical_duplicate_distance;
    }

    ## Infile
    push @commands, q{-INPUT} . $SPACE . join $SPACE . q{-INPUT} . $SPACE,
      @{$infile_paths_ref};

    ## Output
    push @commands, q{-OUTPUT} . $SPACE . $outfile_path;

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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub picardtools_gatherbamfiles {

## Function : Perl wrapper for writing picardtools gatherbamfiles recipe to $filehandle. Based on picardtools 2.20.7.
## Returns  : @commands
## Arguments: $create_index           => Create index
##          : $filehandle             => Sbatch filehandle to write to
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
    my $filehandle;
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
        filehandle       => { store => \$filehandle, },
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
                java_jar                  => $java_jar,
                java_use_large_pages      => $java_use_large_pages,
                memory_allocation         => $memory_allocation,
                picard_use_barclay_parser => 1,
                temp_directory            => $temp_directory,
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
    push @commands, q{-INPUT} . $SPACE . join $SPACE . q{-INPUT} . $SPACE,
      @{$infile_paths_ref};

    ## Output
    push @commands, q{-OUTPUT} . $SPACE . $outfile_path;

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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub picardtools_collectmultiplemetrics {

## Function : Perl wrapper for writing picardtools collectmultiplemetrics recipe to $filehandle. Based on picardtools 2.20.7.
## Returns  : @commands
## Arguments: $create_index           => Create index
##          : $filehandle             => Sbatch filehandle to write to
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
    my $filehandle;
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
        filehandle  => { store => \$filehandle, },
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
                java_jar                  => $java_jar,
                java_use_large_pages      => $java_use_large_pages,
                memory_allocation         => $memory_allocation,
                picard_use_barclay_parser => 1,
                temp_directory            => $temp_directory,
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
    push @commands, q{-INPUT} . $SPACE . $infile_path;

    ## Output
    push @commands, q{-OUTPUT} . $SPACE . $outfile_path;

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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub picardtools_collecthsmetrics {

## Function : Perl wrapper for writing picardtools collecthsmetrics recipe to $filehandle. Based on picardtools 2.20.7.
## Returns  : @commands
## Arguments: $bait_interval_file_paths_ref   => Interval list file(s) that contains the locations of the baits used {REF}
##          : $create_index                   => Create index
##          : $filehandle                     => Sbatch filehandle to write to
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
    my $filehandle;
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
        filehandle  => { store => \$filehandle, },
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
                java_use_large_pages      => $java_use_large_pages,
                java_jar                  => $java_jar,
                memory_allocation         => $memory_allocation,
                picard_use_barclay_parser => 1,
                temp_directory            => $temp_directory,
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
      q{-BAIT_INTERVALS} . $SPACE . join $SPACE . q{-BAIT_INTERVALS} . $SPACE,
      @{$bait_interval_file_paths_ref};

    ## Targets
    push @commands,
      q{-TARGET_INTERVALS} . $SPACE . join $SPACE . q{-TARGET_INTERVALS} . $SPACE,
      @{$target_interval_file_paths_ref};

    ## Infile
    push @commands, q{-INPUT} . $SPACE . $infile_path;

    ## Output
    push @commands, q{-OUTPUT} . $SPACE . $outfile_path;

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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub picardtools_sortvcf {

## Function : Perl wrapper for writing picardtools sortvcf recipe to $filehandle. Based on picardtools 2.20.7.
## Returns  : @commands
## Arguments: $create_index         => Create index
##          : $filehandle           => Sbatch filehandle to write to
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
    my $filehandle;
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
            allow       => [qw{ true false }],
            default     => q{false},
            store       => \$create_index,
            strict_type => 1,
        },
        filehandle       => { store => \$filehandle, },
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
        sequence_dictionary => {
            defined     => 1,
            required    => 1,
            store       => \$sequence_dictionary,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        temp_directory  => { store => \$temp_directory,  strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands;

    ## Return java core commands
    if ($java_jar) {

        @commands = java_core(
            {
                java_jar                  => $java_jar,
                java_use_large_pages      => $java_use_large_pages,
                memory_allocation         => $memory_allocation,
                picard_use_barclay_parser => 1,
                temp_directory            => $temp_directory,
            }
        );
    }

    ## Picardtools mergesamfiles
    push @commands, q{SortVcf};

    ## Picardtools base args
    @commands = picardtools_base(
        {
            commands_ref       => \@commands,
            create_index       => $create_index,
            referencefile_path => $referencefile_path,
        }
    );

    push @commands, q{-SEQUENCE_DICTIONARY} . $SPACE . $sequence_dictionary;

    ## Infile
    push @commands, q{-INPUT} . $SPACE . join $SPACE . q{-INPUT} . $SPACE,
      @{$infile_paths_ref};

    ## Output
    push @commands, q{-OUTPUT} . $SPACE . $outfile_path;

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

sub sort_vcf {

## Function : Writes sbatch code to supplied filehandle to sort variants in vcf format.
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $filehandle            => SBATCH script filehandle to print to
##          : $infile_paths_ref      => Infiles to sort {REF}
##          : $outfile               => The sorted outfile
##          : $sequence_dict_file    => Human reference sequence dict file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $filehandle;
    my $infile_paths_ref;
    my $outfile;
    my $sequence_dict_file;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        filehandle       => { defined => 1, required => 1, store => \$filehandle, },
        infile_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infile_paths_ref,
            strict_type => 1,
        },
        outfile => {
            defined     => 1,
            required    => 1,
            store       => \$outfile,
            strict_type => 1,
        },
        sequence_dict_file => {
            defined     => 1,
            required    => 1,
            store       => \$sequence_dict_file,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    say {$filehandle} q{## Picard SortVcf};

    ## Writes java core commands to filehandle.
    picardtools_sortvcf(
        {
            filehandle           => $filehandle,
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

sub picardtools_collectrnaseqmetrics {

## Function : Perl wrapper for writing picardtools CollectRnaSeqMetrics recipe to $filehandle. Based on picardtools 2.22.3.
## Returns  : @commands
## Arguments: $chart_outfile_path        => Interval list file(s) that contains the locations of the baits used {REF}
##          : $filehandle                => Sbatch filehandle to write to
##          : $gene_annotation_file_path => Gene annotations in refFlat format
##          : $infile_path               => Infile paths
##          : $java_jar                  => Java jar
##          : $java_use_large_pages      => Use java large pages
##          : $memory_allocation         => Memory allocation for java
##          : $outfile_path              => Outfile path
##          : $rrna_intervals_file_path  => rRNA intervals in interval_list format
##          : $stderrfile_path           => Stderrfile path
##          : $strand_specificity        => Strand specificity
##          : $temp_directory            => Redirect tmp files to java temp

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $chart_outfile_path;
    my $filehandle;
    my $gene_annotation_file_path;
    my $infile_path;
    my $java_jar;
    my $memory_allocation;
    my $outfile_path;
    my $rrna_intervals_file_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $strand_specificity;
    my $temp_directory;

    ## Default(s)
    my $java_use_large_pages;

    my $tmpl = {
        chart_outfile_path => {
            store       => \$chart_outfile_path,
            strict_type => 1,
        },
        filehandle                => { store => \$filehandle, },
        gene_annotation_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$gene_annotation_file_path,
            strict_type => 1,
        },
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
        rrna_intervals_file_path => {
            allow       => [ undef, qr{ [.]interval_list \z}xms ],
            store       => \$rrna_intervals_file_path,
            strict_type => 1,
        },
        stderrfile_path        => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append => {
            store       => \$stderrfile_path_append,
            strict_type => 1,
        },
        strand_specificity => {
            allow       => [qw{ unstranded forward_stranded reverse_stranded }],
            required    => 1,
            store       => \$strand_specificity,
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
                java_use_large_pages      => $java_use_large_pages,
                java_jar                  => $java_jar,
                memory_allocation         => $memory_allocation,
                picard_use_barclay_parser => 1,
                temp_directory            => $temp_directory,
            }
        );
    }

    ## Picardtools collecthsmetrics
    push @commands, q{CollectRnaSeqMetrics};

    ## Picardtools base args
    @commands = picardtools_base(
        {
            commands_ref => \@commands,
        }
    );

    if ($chart_outfile_path) {

        push @commands, q{-CHART_OUTPUT} . $SPACE . $chart_outfile_path;
    }

    push @commands, q{-REF_FLAT} . $SPACE . $gene_annotation_file_path;

    push @commands, q{-INPUT} . $SPACE . $infile_path;

    push @commands, q{-OUTPUT} . $SPACE . $outfile_path;

    if ($rrna_intervals_file_path) {

        push @commands, q{-RIBOSOMAL_INTERVALS} . $SPACE . $rrna_intervals_file_path;
    }

    if ( $strand_specificity eq q{unstranded} ) {

        push @commands, q{-STRAND_SPECIFICITY NONE};
    }
    elsif ( $strand_specificity eq q{forward_stranded} ) {

        push @commands, q{-STRAND_SPECIFICITY FIRST_READ_TRANSCRIPTION_STRAND};
    }
    elsif ( $strand_specificity eq q{reverse_stranded} ) {

        push @commands, q{-STRAND_SPECIFICITY SECOND_READ_TRANSCRIPTION_STRAND};
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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

1;
