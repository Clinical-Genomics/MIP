package MIP::Program::Alignment::Sambamba;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use Carp;
use utf8;    # Allow unicode characters in this script
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };

use FindBin qw{$Bin};    # Find directory of script
use File::Basename qw{dirname};
use File::Spec::Functions qw{catdir};

BEGIN {
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.00;

    # Inherit from Exporter to export functions and variables
    use base qw {Exporter};

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{sambamba_view sambamba_index sambamba_sort sambamba_markdup sambamba_flagstat sambamba_depth};
}

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use MIP::Unix::Standard_streams qw{unix_standard_streams};
use MIP::Unix::Write_to_file qw{unix_write_to_file};

use Params::Check qw{check allow last_error};
use Readonly;

## Constants
Readonly my $SPACE => q{ };

sub sambamba_view {

## sambamba_view

## Function : Perl wrapper for writing sambamba view recipe to $FILEHANDLE. Based on sambamba 0.6.5
## Returns  : "@commands"
## Arguments: $regions_ref, $FILEHANDLE, $infile_path, $outfile_path, $stderrfile_path, $with_header, $show_progress, $output_format, $referencefile_path
##          : $regions_ref        => The regions to process {REF}
##          : $FILEHANDLE         => Sbatch filehandle to write to
##          : $infile_path        => Infile path
##          : $outfile_path       => Outfile path
##          : $stderrfile_path    => Stderrfile path
##          : $referencefile_path => Reference for writing CRAM
##          : $with_header        => Include header
##          : $show_progress      => Show progress
##          : $output_format      => Output format

    my ($arg_href) = @_;

    ## Default(s)
    my $with_header;
    my $show_progress;
    my $output_format;

    ## Flatten argument(s)
    my $regions_ref;
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $referencefile_path;
    my $stderrfile_path_append;

    my $tmpl = {
        regions_ref =>
          { default => [], strict_type => 1, store => \$regions_ref },
        FILEHANDLE  => { required => 1, store => \$FILEHANDLE },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        outfile_path    => { strict_type => 1, store => \$outfile_path },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        referencefile_path =>
          { strict_type => 1, store => \$referencefile_path },
        with_header => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$with_header
        },
        show_progress => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$show_progress
        },
        output_format => {
            default     => qw{bam},
            allow       => [qw{sam bam cram json}],
            strict_type => 1,
            store       => \$output_format
        },
    };

    check( $tmpl, $arg_href, 1 )
      or croak qw{Could not parse arguments!};

    ## Array @commands stores commands depending on input parameters
    my @commands = qw{sambamba view};

    if ($with_header) {    # Include header

        push @commands, q{--with-header};
    }

    if ($output_format) {

        # Output format
        push @commands, q{--format} . $SPACE . $output_format;
    }

    if ($referencefile_path) {

        # Reference for writing CRAM
        push @commands, q{--ref-filename=} . $referencefile_path;
    }

    if ($show_progress) {

# Show progressbar in STDERR (works only for BAM files with no regions specified)
        push @commands, q{--show-progress};
    }

    if ($outfile_path) {

        # Specify output filename
        push @commands, q{--output-filename=} . $outfile_path;
    }

    ## Infile
    push @commands, $infile_path;

    if (@$regions_ref) {

        # Limit output to regions
        push @commands, join $SPACE, @{$regions_ref};
    }

    if ($stderrfile_path) {

        push @commands,
          unix_standard_streams(
            {
                stderrfile_path        => $stderrfile_path,
                stderrfile_path_append => $stderrfile_path_append,
            }
          );
    }
    if ($FILEHANDLE) {

        unix_write_to_file(
            {
                commands_ref => \@commands,
                separator    => $SPACE,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
    }
    return @commands;
}

sub sambamba_index {

## sambamba_index

## Function : Perl wrapper for writing sambamba index recipe to $FILEHANDLE. Based on sambamba 0.6.5
## Returns  : "@commands"
## Arguments: $FILEHANDLE, $infile_path, $stderrfile_path, $show_progress
##          : $FILEHANDLE      => Sbatch filehandle to write to
##          : $infile_path     => Infile path
##          : $stderrfile_path => Stderrfile path
##          : $show_progress   => Show progress

    my ($arg_href) = @_;

    ## Default(s)
    my $show_progress;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;

    my $tmpl = {
        FILEHANDLE  => { required => 1, store => \$FILEHANDLE },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        show_progress   => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$show_progress
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    ## @commands stores commands depending on input parameters
    my @commands = qw{sambamba index};

    if ($show_progress) {

# Show progressbar in STDERR (works only for BAM files with no regions specified)
        push @commands, q{--show-progress};
    }

    ## Infile
    push @commands, $infile_path;

    if ($stderrfile_path) {

        # Redirect stderr output to program specific stderr file
        push @commands,
          unix_standard_streams(
            {
                stderrfile_path        => $stderrfile_path,
                stderrfile_path_append => $stderrfile_path_append,
            }
          );
    }
    if ($FILEHANDLE) {

        unix_write_to_file(
            {
                commands_ref => \@commands,
                separator    => $SPACE,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
    }
    return @commands;
}

sub sambamba_sort {

## sambamba_sort

## Function : Perl wrapper for writing sambamba sort recipe to $FILEHANDLE. Based on sambamba 0.6.5
## Returns  : "@commands"
## Arguments: $FILEHANDLE, $infile_path, $outfile_path, $stderrfile_path, $show_progress, $memory_limit, $temp_directory
##          : $FILEHANDLE      => Sbatch filehandle to write to
##          : $infile_path     => Infile path
##          : $outfile_path    => Outfile path
##          : $stderrfile_path => Stderrfile path
##          : $show_progress   => Show progress
##          : $memory_limit    => Approximate total memory limit for all threads
##          : $temp_directory  => Directory for storing intermediate files; default is system directory for temporary files

    my ($arg_href) = @_;

    ## Default(s)
    my $show_progress;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $memory_limit;
    my $temp_directory;
    my $stderrfile_path_append;

    my $tmpl = {
        FILEHANDLE  => { required => 1, store => \$FILEHANDLE },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        outfile_path    => { strict_type => 1, store => \$outfile_path },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        memory_limit    => {
            allow       => qr/^\d+G$/,
            strict_type => 1,
            store       => \$memory_limit
        },
        temp_directory => { strict_type => 1, store => \$temp_directory },
        show_progress  => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$show_progress
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    ## @commands stores commands depending on input parameters
    my @commands = qw{sambamba sort};

    if ($show_progress) {

# Show progressbar in STDERR (works only for BAM files with no regions specified)
        push @commands, q{--show-progress};
    }

    if ($memory_limit) {

        # Approximate total memory limit for all threads
        push @commands, q{--memory-limit=} . $memory_limit;
    }

    if ($temp_directory) {

        # Directory for storing intermediate files
        push @commands, q{--tmpdir=} . $temp_directory;
    }

    ## Outfile
    if ($outfile_path) {

        push @commands, q{--out=} . $outfile_path;
    }

    ## Infile
    if ($infile_path) {

        push @commands, $infile_path;
    }
    if ($stderrfile_path) {

        # Redirect stderr output to program specific stderr file
        push @commands,
          unix_standard_streams(
            {
                stderrfile_path        => $stderrfile_path,
                stderrfile_path_append => $stderrfile_path_append,
            }
          );
    }
    if ($FILEHANDLE) {

        unix_write_to_file(
            {
                commands_ref => \@commands,
                separator    => $SPACE,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
    }
    return @commands;
}

sub sambamba_markdup {

## sambamba_markdup

## Function : Perl wrapper for writing sambamba markdup recipe to $FILEHANDLE. Based on sambamba 0.6.5
## Returns  : "@commands"
## Arguments: $infile_path, $outfile_path, $stderrfile_path, $FILEHANDLE, temp_directory, $show_progress, $hash_table_size, $overflow_list_size, $io_buffer_size
##          : $infile_path        => Infile path
##          : $outfile_path       => Outfile path
##          : $stderrfile_path    => Stderrfile path
##          : $FILEHANDLE         => Sbatch filehandle to write to
##          : $temp_directory     => Specify directory for temporary files
##          : $show_progress      => Show progress
##          : $hash_table_size    => Size of hash table for finding read pairs
##          : $overflow_list_size => Size of the overflow list where reads, thrown from the hash table, get a second chance to meet their pairs
##          : $io_buffer_size     => Two buffers of BUFFER_SIZE *megabytes* each are used for reading and writing BAM during the second pass

    my ($arg_href) = @_;

    ## Default(s)
    my $show_progress;

    ## Flatten argument(s)
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $temp_directory;
    my $hash_table_size;
    my $overflow_list_size;
    my $io_buffer_size;
    my $stderrfile_path_append;

    my $tmpl = {
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        outfile_path    => { strict_type => 1, store => \$outfile_path },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        FILEHANDLE      => { required    => 1, store => \$FILEHANDLE },
        temp_directory  => { strict_type => 1, store => \$temp_directory },
        hash_table_size => {
            allow       => qr/^\d+$/,
            strict_type => 1,
            store       => \$hash_table_size
        },
        overflow_list_size => {
            allow       => qr/^\d+$/,
            strict_type => 1,
            store       => \$overflow_list_size
        },
        io_buffer_size => {
            allow       => qr/^\d+$/,
            strict_type => 1,
            store       => \$io_buffer_size
        },
        show_progress => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$show_progress
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    ## @commands stores commands depending on input parameters
    my @commands = qw{sambamba markdup};

    if ($temp_directory) {

        # Specify directory for temporary files
        push @commands, q{--tmpdir=} . $temp_directory;
    }
    if ($hash_table_size) {

        # Size of hash table for finding read pairs
        push @commands, q{--hash-table-size=} . $hash_table_size;
    }
    if ($overflow_list_size) {

        # Size of the overflow list
        push @commands, q{--overflow-list-size=} . $overflow_list_size;
    }
    if ($io_buffer_size) {

        # Two buffers of BUFFER_SIZE *megabytes*
        push @commands, q{--io-buffer-size=} . $io_buffer_size;
    }
    if ($show_progress) {

# Show progressbar in STDERR (works only for BAM files with no regions specified)
        push @commands, q{--show-progress};
    }

    ## Infile
    push @commands, $infile_path;

    ## Outfile
    if ($outfile_path) {

        push( @commands, $outfile_path );
    }
    if ($stderrfile_path) {

        # Redirect stderr output to program specific stderr file
        push @commands,
          unix_standard_streams(
            {
                stderrfile_path        => $stderrfile_path,
                stderrfile_path_append => $stderrfile_path_append,
            }
          );
    }
    if ($FILEHANDLE) {

        unix_write_to_file(
            {
                commands_ref => \@commands,
                separator    => $SPACE,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
    }
    return @commands;
}

sub sambamba_flagstat {

## sambamba_flagstat

## Function : Perl wrapper for writing sambamba flagstat recipe to $FILEHANDLE. Based on sambamba 0.6.5
## Returns  : "@commands"
## Arguments: $FILEHANDLE, $infile_path, $outfile_path, $stderrfile_path, $show_progress, $memory_limit, $temp_directory
##          : $FILEHANDLE      => Sbatch filehandle to write to
##          : $infile_path     => Infile path
##          : $outfile_path    => Outfile path
##          : $stderrfile_path => Stderrfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $stderrfile_path_append;

    my $tmpl = {
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        outfile_path    => { strict_type => 1, store => \$outfile_path },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        FILEHANDLE => { store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    ##  @commands stores commands depending on input parameters
    my @commands = qw{sambamba flagstat};

    ## Infile
    if ($infile_path) {

        push @commands, $infile_path;
    }
    ## Outfile
    if ($outfile_path) {

        push @commands, q{>} . $SPACE . $outfile_path;
    }
    if ($stderrfile_path) {

        # Redirect stderr output to program specific stderr file
        push @commands,
          unix_standard_streams(
            {
                stderrfile_path        => $stderrfile_path,
                stderrfile_path_append => $stderrfile_path_append,
            }
          );
    }
    if ($FILEHANDLE) {

        unix_write_to_file(
            {
                commands_ref => \@commands,
                separator    => $SPACE,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
    }
    return @commands;
}

sub sambamba_depth {

## sambamba_depth

## Function : Perl wrapper for writing sambamba depth recipe to $FILEHANDLE. Based on sambamba 0.6.5
## Returns  : "@commands"
## Arguments: $depth_cutoffs_ref, $FILEHANDLE, $infile_path, $outfile_path, $stderrfile_path, $region, $filter, $min_base_quality, $mode, $fix_mate_overlap
##          : $depth_cutoffs_ref => Multiple thresholds can be provided, for each one an extra column will be added, the percentage of bases in the region where coverage is more than this value {REF}
##          : $FILEHANDLE        => Sbatch filehandle to write to
##          : $infile_path       => Infile path
##          : $outfile_path      => Outfile path
##          : $stderrfile_path   => Stderrfile path
##          : $region            => List or regions of interest or a single region in form chr:beg-end
##          : $filter            => Set custom filter for alignments
##          : $fix_mate_overlap  => Detect overlaps of mate reads and handle them on per-base basis
##          : $min_base_quality  => Don't count bases with lower base quality
##          : $mode              => Mode unit to print the statistics on

    my ($arg_href) = @_;

    ## Default(s)
    my $min_base_quality;
    my $mode;
    my $fix_mate_overlap;

    ## Flatten argument(s)
    my $depth_cutoffs_ref;
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $region;
    my $filter;
    my $stderrfile_path_append;

    my $tmpl = {
        depth_cutoffs_ref =>
          { default => [], strict_type => 1, store => \$depth_cutoffs_ref },
        FILEHANDLE  => { required => 1, store => \$FILEHANDLE },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        outfile_path     => { strict_type => 1, store => \$outfile_path },
        stderrfile_path  => { strict_type => 1, store => \$stderrfile_path },
        region           => { strict_type => 1, store => \$region },
        filter           => { strict_type => 1, store => \$filter },
        fix_mate_overlap => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$fix_mate_overlap
        },
        min_base_quality => {
            default     => 0,
            allow       => qr/^\d+$/,
            strict_type => 1,
            store       => \$min_base_quality
        },
        mode => {
            default     => qw{region},
            allow       => [qw{base region window}],
            strict_type => 1,
            store       => \$mode
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    ## @commands stores commands depending on input parameters
    my @commands = qw{sambamba depth};

    if ($mode) {

        push @commands, $mode;
    }
    if ($region) {    # Limit output to regions

        push @commands, q{--regions} . $SPACE . $region;
    }
    if ( @{ $depth_cutoffs_ref } ) {

        push @commands,
          q{--cov-threshold} . $SPACE . join q{ --cov-threshold} . $SPACE,
          @{ $depth_cutoffs_ref };
    }
    if ($min_base_quality) {

        push @commands, q{--min-base-quality} . $SPACE . $min_base_quality;
    }
    if ($fix_mate_overlap) {

        push @commands, q{--fix-mate-overlaps};
    }
    if ($filter) {

        push @commands, q{--filter} . $SPACE . $filter;
    }
    if ($outfile_path) {

        # Specify output filename
        push @commands, q{--output-filename=} . $outfile_path;
    }

    ## Infile
    push( @commands, $infile_path );

    if ($stderrfile_path) {

        # Redirect stderr output to program specific stderr file
        push @commands,
          unix_standard_streams(
            {
                stderrfile_path        => $stderrfile_path,
                stderrfile_path_append => $stderrfile_path_append,
            }
          );
    }
    if ($FILEHANDLE) {

        unix_write_to_file(
            {
                commands_ref => \@commands,
                separator    => $SPACE,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
    }
    return @commands;
}

1;
