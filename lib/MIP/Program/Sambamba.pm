package MIP::Program::Sambamba;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use File::Spec::Functions qw{ catdir catfile };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $DOT $NEWLINE $SEMICOLON $SPACE $UNDERSCORE };
use MIP::Environment::Executable qw{ get_executable_base_command };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.05;

    # Inherit from Exporter to export functions and variables
    use base qw {Exporter};

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ sambamba_depth split_and_index_aligment_file sambamba_flagstat sambamba_index sambamba_markdup sambamba_sort sambamba_view };
}

Readonly my $BASE_COMMAND => q{sambamba};

sub sambamba_view {

## Function : Perl wrapper for writing sambamba view recipe to $filehandle. Based on sambamba 0.6.5
## Returns  : @commands
## Arguments: $filehandle             => Sbatch filehandle to write to
##          : $infile_path            => Infile path
##          : $outfile_path           => Outfile path
##          : $output_format          => Output format
##          : $referencefile_path     => Reference for writing CRAM
##          : $regions_ref            => The regions to process {REF}
##          : $show_progress          => Show progress
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Stderrfile path append
##          : $with_header            => Include header

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $outfile_path;
    my $referencefile_path;
    my $regions_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;

    ## Default(s)
    my $output_format;
    my $show_progress;
    my $with_header;

    my $tmpl = {
        filehandle  => { required => 1, store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path  => { store => \$outfile_path, strict_type => 1, },
        output_format => {
            allow       => [qw{ sam bam cram json }],
            default     => q{ bam },
            store       => \$output_format,
            strict_type => 1,
        },
        referencefile_path => { store => \$referencefile_path, strict_type => 1, },
        regions_ref   => { default => [], store => \$regions_ref, strict_type => 1, },
        show_progress => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$show_progress,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        with_header => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$with_header,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands =
      ( get_executable_base_command( { base_command => $BASE_COMMAND, } ), qw{ view } );

    if ($with_header) {

        push @commands, q{--with-header};
    }

    if ($output_format) {

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

        push @commands, q{--output-filename=} . $outfile_path;
    }

    push @commands, $infile_path;

    if ( @{$regions_ref} ) {

        # Limit output to regions
        push @commands, join $SPACE, @{$regions_ref};
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

sub sambamba_index {

## Function : Perl wrapper for writing sambamba index recipe to $filehandle. Based on sambamba 0.6.5
## Returns  : @commands
## Arguments: $filehandle             => Sbatch filehandle to write to
##          : $infile_path            => Infile path
##          : $show_progress          => Show progress
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Stderrfile path append

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;

    ## Default(s)
    my $show_progress;

    my $tmpl = {
        filehandle  => { required => 1, store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        show_progress => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$show_progress,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands =
      ( get_executable_base_command( { base_command => $BASE_COMMAND, } ), qw{ index } );

    if ($show_progress) {

        # Show progressbar in STDERR (works only for BAM files with no regions specified)
        push @commands, q{--show-progress};
    }

    push @commands, $infile_path;

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

sub sambamba_sort {

## Function : Perl wrapper for writing sambamba sort recipe to $filehandle. Based on sambamba 0.6.5
## Returns  : @commands
## Arguments: $filehandle             => Sbatch filehandle to write to
##          : $infile_path            => Infile path
##          : $memory_limit           => Approximate total memory limit for all threads
##          : $outfile_path           => Outfile path
##          : $show_progress          => Show progress
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Stderrfile path append
##          : $temp_directory         => Directory for storing intermediate files; default is system directory for temporary files

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $memory_limit;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $temp_directory;

    ## Default(s)
    my $show_progress;

    my $tmpl = {
        filehandle  => { required => 1, store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        memory_limit => {
            allow       => qr/ ^\d+G$ /sxm,
            store       => \$memory_limit,
            strict_type => 1,
        },
        outfile_path  => { store => \$outfile_path, strict_type => 1, },
        show_progress => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$show_progress,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        temp_directory => { store => \$temp_directory, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands =
      ( get_executable_base_command( { base_command => $BASE_COMMAND, } ), qw{ sort } );

    if ($show_progress) {

        # Show progressbar in STDERR (works only for BAM files with no regions specified)
        push @commands, q{--show-progress};
    }

    if ($memory_limit) {

        # Approximate total memory limit for all threads
        push @commands, q{--memory-limit=} . $memory_limit;
    }

    if ($temp_directory) {

        push @commands, q{--tmpdir=} . $temp_directory;
    }

    if ($outfile_path) {

        push @commands, q{--out=} . $outfile_path;
    }

    if ($infile_path) {

        push @commands, $infile_path;
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

sub sambamba_markdup {

## Function : Perl wrapper for writing sambamba markdup recipe to $filehandle. Based on sambamba 0.6.5
## Returns  : @commands
## Arguments: $filehandle             => Sbatch filehandle to write to
##          : $hash_table_size        => Size of hash table for finding read pairs
##          : $infile_path            => Infile path
##          : $io_buffer_size         => Two buffers of BUFFER_SIZE *megabytes* each are used for reading and writing BAM during the second pass
##          : $overflow_list_size     => Size of the overflow list where reads, thrown from the hash table, get a second chance to meet their pairs
##          : $show_progress          => Show progress
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Stderrfile path append
##          : $stdoutfile_path        => Standard outfile path
##          : $temp_directory         => Specify directory for temporary files

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $hash_table_size;
    my $infile_path;
    my $io_buffer_size;
    my $overflow_list_size;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $temp_directory;

    ## Default(s)
    my $show_progress;

    my $tmpl = {
        filehandle      => { required => 1, store => \$filehandle, },
        hash_table_size => {
            allow       => qr/ ^\d+$ /sxm,
            store       => \$hash_table_size,
            strict_type => 1,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        io_buffer_size => {
            allow       => qr/ ^\d+$ /sxm,
            store       => \$io_buffer_size,
            strict_type => 1,
        },
        overflow_list_size => {
            allow       => qr/ ^\d+$ /sxm,
            store       => \$overflow_list_size,
            strict_type => 1,
        },
        show_progress => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$show_progress,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
        temp_directory  => { store => \$temp_directory,  strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = ( get_executable_base_command( { base_command => $BASE_COMMAND, } ),
        qw{ markdup } );

    if ($temp_directory) {

        push @commands, q{--tmpdir=} . $temp_directory;
    }

    if ($hash_table_size) {

        # Size of hash table for finding read pairs
        push @commands, q{--hash-table-size=} . $hash_table_size;
    }

    if ($overflow_list_size) {

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

    push @commands, $infile_path;

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

sub sambamba_flagstat {

## Function : Perl wrapper for writing sambamba flagstat recipe to $filehandle. Based on sambamba 0.6.5
## Returns  : @commands
## Arguments: $filehandle             => Sbatch filehandle to write to
##          : $infile_path            => Infile path
##          : $outfile_path           => Outfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Stderrfile path append

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;

    my $tmpl = {
        filehandle  => { store => \$filehandle },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path    => { store => \$outfile_path,    strict_type => 1, },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = (
        get_executable_base_command( { base_command => $BASE_COMMAND, } ),
        qw{ flagstat }
    );

    if ($infile_path) {

        push @commands, $infile_path;
    }

    if ($outfile_path) {

        push @commands, q{>} . $SPACE . $outfile_path;
    }

    # Redirect stderr output to program specific stderr file
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

sub sambamba_depth {

## Function : Perl wrapper for writing sambamba depth recipe to $filehandle. Based on sambamba 0.6.5
## Returns  : @commands
## Arguments: $depth_cutoffs_ref      => Multiple thresholds can be provided, for each one an extra column will be added, the percentage of bases in the region where coverage is more than this value {REF}
##          : $filehandle             => Sbatch filehandle to write to
##          : $filter                 => Set custom filter for alignments
##          : $fix_mate_overlap       => Detect overlaps of mate reads and handle them on per-base basis
##          : $infile_path            => Infile path
##          : $min_base_quality       => Don't count bases with lower base quality
##          : $mode                   => Mode unit to print the statistics on
##          : $outfile_path           => Outfile path
##          : $region                 => List or regions of interest or a single region in form chr:beg-end
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Stderrfile path append

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $depth_cutoffs_ref;
    my $filehandle;
    my $filter;
    my $infile_path;
    my $outfile_path;
    my $region;
    my $stderrfile_path;
    my $stderrfile_path_append;

    ## Default(s)
    my $fix_mate_overlap;
    my $min_base_quality;
    my $mode;

    my $tmpl = {
        depth_cutoffs_ref =>
          { default => [], store => \$depth_cutoffs_ref, strict_type => 1, },
        filehandle       => { required => 1,        store       => \$filehandle, },
        filter           => { store    => \$filter, strict_type => 1, },
        fix_mate_overlap => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$fix_mate_overlap,
            strict_type => 1,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        min_base_quality => {
            allow       => qr/ ^\d+$ /sxm,
            default     => 0,
            store       => \$min_base_quality,
            strict_type => 1,
        },
        mode => {
            allow       => [qw{ base region window }],
            default     => q{region},
            store       => \$mode,
            strict_type => 1,
        },
        outfile_path    => { store => \$outfile_path,    strict_type => 1, },
        region          => { store => \$region,          strict_type => 1, },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands =
      ( get_executable_base_command( { base_command => $BASE_COMMAND, } ), qw{ depth } );

    if ($mode) {

        push @commands, $mode;
    }

    if ($region) {

        push @commands, q{--regions} . $SPACE . $region;
    }

    if ( @{$depth_cutoffs_ref} ) {

        push @commands,
          q{--cov-threshold} . $SPACE . join q{ --cov-threshold} . $SPACE,
          @{$depth_cutoffs_ref};
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

        push @commands, q{--output-filename=} . $outfile_path;
    }

    push @commands, $infile_path;

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

sub split_and_index_aligment_file {

## Function : Split alignemnt file per contig and index new file. Creates the command line for xargs. Writes to sbatch filehandle and opens xargs filehandle
## Returns  : $xargs_file_counter
## Arguments: $contigs_ref         => Contigs to process {REF}
##          : $core_number         => Number of cores to use
##          : $filehandle          => Sbatch filehandle to write to
##          : $file_path           => File name - ususally sbatch
##          : $infile_path         => Infile path
##          : $outfile_path_prefix => Outfile path prefix
##          : $output_format       => Output format
##          : $recipe_info_path    => Program info path
##          : $temp_directory      => Temporary directory
##          : $xargsfilehandle     => XARGS filehandle to write to
##          : $xargs_file_counter  => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $contigs_ref;
    my $core_number;
    my $filehandle;
    my $file_path;
    my $infile_path;
    my $outfile_path_prefix;
    my $recipe_info_path;
    my $xargsfilehandle;

    ## Default(s)
    my $output_format;
    my $xargs_file_counter;

    my $tmpl = {
        contigs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$contigs_ref,
            strict_type => 1,
        },
        core_number => {
            defined     => 1,
            required    => 1,
            store       => \$core_number,
            strict_type => 1,
        },
        filehandle => { defined => 1, required => 1, store => \$filehandle, },
        file_path  => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        infile_path =>
          { defined => 1, required => 1, store => \$infile_path, strict_type => 1, },
        outfile_path_prefix => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path_prefix,
            strict_type => 1,
        },
        output_format => {
            allow       => [qw{ sam bam cram json }],
            default     => q{bam},
            store       => \$output_format,
            strict_type => 1,
        },
        recipe_info_path => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_info_path,
            strict_type => 1,
        },
        xargsfilehandle => { defined => 1, required => 1, store => \$xargsfilehandle, },
        xargs_file_counter => {
            allow       => qr/ ^\d+$ /xsm,
            default     => 0,
            store       => \$xargs_file_counter,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };

    my $xargs_file_path_prefix;

    my $file_suffix = $DOT . $output_format;

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number        => $core_number,
            filehandle         => $filehandle,
            file_path          => $file_path,
            recipe_info_path   => $recipe_info_path,
            xargsfilehandle    => $xargsfilehandle,
            xargs_file_counter => $xargs_file_counter,
        }
    );

    ## Split by contig
  CONTIG:
    foreach my $contig ( @{$contigs_ref} ) {

        ## Get parameters
        my $outfile_path =
          catfile( $outfile_path_prefix . $DOT . $contig . $file_suffix );
        my $stderrfile_path_view =
            $xargs_file_path_prefix
          . $UNDERSCORE . q{view}
          . $UNDERSCORE
          . $contig
          . $DOT
          . q{stderr.txt};

        sambamba_view(
            {
                filehandle      => $xargsfilehandle,
                infile_path     => $infile_path,
                outfile_path    => $outfile_path,
                output_format   => $output_format,
                regions_ref     => [$contig],
                show_progress   => 1,
                stderrfile_path => $stderrfile_path_view,
                with_header     => 1,
            }
        );

        # Seperate commands
        print {$xargsfilehandle} $SEMICOLON . $SPACE;

        my $stderrfile_path_index =
            $xargs_file_path_prefix
          . $UNDERSCORE
          . q{index}
          . $UNDERSCORE
          . $contig
          . $DOT
          . q{stderr.txt};
        sambamba_index(
            {
                filehandle      => $xargsfilehandle,
                infile_path     => $outfile_path,
                show_progress   => 1,
                stderrfile_path => $stderrfile_path_index,
            }
        );
        print {$xargsfilehandle} $NEWLINE;
    }
    return $xargs_file_counter;
}

1;
