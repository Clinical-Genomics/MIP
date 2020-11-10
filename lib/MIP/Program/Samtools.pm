package MIP::Program::Samtools;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $AMPERSAND $COMMA $SPACE };
use MIP::Environment::Executable qw{ get_executable_base_command };
use MIP::Processmanagement::Processes qw{ print_wait };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.12;

    # Inherit from Exporter to export functions and variables
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      samtools_base
      samtools_create_chromosome_files
      samtools_depth
      samtools_faidx
      samtools_flagstat
      samtools_idxstats
      samtools_index
      samtools_merge
      samtools_stats
      samtools_sort
      samtools_view
    };
}

Readonly my $BASE_COMMAND => q{samtools};

sub samtools_base {

## Function : Perl wrapper for samtools base. Based on Samtools 1.10
## Returns  : @commands
## Arguments: $commands_ref       => List of commands added earlier
##          : $filehandle         => Filehandle to write to
##          : $output_format      => Output format
##          : $referencefile_path => Reference file path (fasta)
##          : $thread_number      => Number of BAM/CRAM compression threads
##          : $write_index        => Write index while writing file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $commands_ref;
    my $filehandle;
    my $referencefile_path;

    ## Default(s)
    my $output_format;
    my $thread_number;
    my $write_index;

    my $tmpl = {
        commands_ref => {
            default     => [],
            store       => \$commands_ref,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        output_format => {
            allow       => [qw{ bam cram sam }],
            default     => q{bam},
            store       => \$output_format,
            strict_type => 1,
        },
        referencefile_path => {
            store       => \$referencefile_path,
            strict_type => 1,
        },
        thread_number => {
            allow       => [ undef, qr{ \A\d+\z }xsm, ],
            store       => \$thread_number,
            strict_type => 1,
        },
        write_index => {
            allow       => [ undef, 0, 1, ],
            default     => 0,
            store       => \$write_index,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = @{$commands_ref};

    if ($output_format) {

        push @commands, q{--output-fmt} . $SPACE . uc $output_format;
    }

    if ($referencefile_path) {

        push @commands, q{--reference} . $SPACE . $referencefile_path;
    }

    if ($thread_number) {

        push @commands, q{--threads} . $SPACE . $thread_number;
    }

    if ($write_index) {

        push @commands, q{--write-index};
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

sub samtools_create_chromosome_files {

## Function : Perl wrapper for writing chromosome files used by other scripts. Writes to filehandle.
## Returns  :
##          : $filehandle         => Sbatch filehandle to write to
##          : $infile_path        => Infile path
##          : $max_process_number => Max number of processeses
##          : $outfile_path       => Outfile path
##          : $regions_ref        => Regions to process {REF}
##          : $suffix             => Suffix to append to outfile_path
##          : $temp_directory     => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $max_process_number;
    my $regions_ref;
    my $suffix;
    my $temp_directory;

    my $tmpl = {
        filehandle  => { store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        max_process_number => {
            allow       => qr{ \A\d+\z }sxm,
            store       => \$max_process_number,
            strict_type => 1,
        },
        regions_ref => {
            default     => [],
            required    => 1,
            store       => \$regions_ref,
            strict_type => 1,
        },
        suffix         => { store => \$suffix, strict_type => 1, },
        temp_directory => {
            required    => 1,
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $process_batches_count = 1;

  CONTIG:
    while ( my ( $contig_index, $contig ) = each @{$regions_ref} ) {
        $process_batches_count = print_wait(
            {
                filehandle            => $filehandle,
                max_process_number    => $max_process_number,
                process_batches_count => $process_batches_count,
                process_counter       => $contig_index,
            }
        );

        samtools_faidx(
            {
                filehandle   => $filehandle,
                infile_path  => $infile_path,
                outfile_path => catfile( $temp_directory, $contig . $suffix ),
                regions_ref  => [$contig],
            }
        );
        say {$filehandle} $SPACE . $AMPERSAND;
    }
    return;
}

sub samtools_depth {

## Function : Perl wrapper for writing samtools depth recipe to $filehandle. Based on samtools 1.10 (using htslib 1.10).
## Returns  : @commands
## Arguments: $filehandle             => Filehandle to write to
##          : $infile_path            => Infile path
##          : $max_depth_treshold     => Set the depth value treshold, samtools depth defaults to 8000 if left unset
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $max_depth_treshold;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        max_depth_treshold => {
            allow   => qr{ \A\d+\z }xsm,
            defined => 1,
            store   => \$max_depth_treshold,
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands =
      ( get_executable_base_command( { base_command => $BASE_COMMAND, } ), qw{ depth } );

    if ($max_depth_treshold) {

        push @commands, q{-d} . $SPACE . $max_depth_treshold;
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

sub samtools_faidx {

## Function : Perl wrapper for writing samtools faidx recipe to $filehandle. Based on samtools 1.10 (using htslib 1.10).
## Returns  : @commands
##          : $filehandle             => Sbatch filehandle to write to
##          : $infile_path            => Infile path
##          : $outfile_path           => Outfile path
##          : $regions_ref            => The regions to process {REF}
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Stderrfile path append

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $outfile_path;
    my $regions_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;

    my $tmpl = {
        filehandle  => { store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path => { store => \$outfile_path, strict_type => 1, },
        regions_ref => { default => [], store => \$regions_ref, strict_type => 1, },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands =
      ( get_executable_base_command( { base_command => $BASE_COMMAND, } ), qw{ faidx } );

    push @commands, $infile_path;

    if ( @{$regions_ref} ) {

        push @commands, join $SPACE, @{$regions_ref};
    }

    if ($outfile_path) {

        push @commands, q{>} . $SPACE . $outfile_path;
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

sub samtools_flagstat {

## Function : Perl wrapper for writing samtools flag stat recipe to $filehandle. Based on samtools 1.10 (using htslib 1.10).
## Returns  : @commands
## Arguments: $filehandle             => Filehandle to write to
##          : $infile_path            => Infile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = (
        get_executable_base_command( { base_command => $BASE_COMMAND, } ),
        qw{ flagstat }
    );

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

sub samtools_idxstats {

## Function : Perl wrapper for writing samtools idxstats recipe to $filehandle. Based on samtools 1.10 (using htslib 1.10).
## Returns  : @commands
## Arguments: $filehandle             => Filehandle to write to
##          : $infile_path            => Infile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = (
        get_executable_base_command( { base_command => $BASE_COMMAND, } ),
        qw{ idxstats }
    );

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

sub samtools_index {

## Function : Perl wrapper for writing samtools index recipe to $filehandle. Based on samtools 1.10 (using htslib 1.10).
## Returns  : @commands
##          : $bai_format             => Generate BAI-format index for BAM files
##          : $filehandle             => Sbatch filehandle to write to
##          : $infile_path            => Infile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Stderrfile path append
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bai_format;
    my $filehandle;
    my $infile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        bai_format  => { store => \$bai_format, strict_type => 1, },
        filehandle  => { store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands =
      ( get_executable_base_command( { base_command => $BASE_COMMAND, } ), qw{ index } );

    if ($bai_format) {

        push @commands, q{-b};
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

sub samtools_merge {

## Function : Perl wrapper for writing samtools merge recipe to $filehandle. Based on samtools 1.10 (using htslib 1.10).
## Returns  : @commands
##          : $filehandle             => Sbatch filehandle to write to
##          : $force                  => Overwrite output file path
##          : $infile_paths_ref       => Infile paths {REF}
##          : $outfile_path           => Outfile path
##          : $output_format          => Output format
##          : $referencefile_path     => Reference file path (fasta)
##          : $region                 => Region to merge
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Stderrfile path append
##          : $stdoutfile_path        => Stdoutfile path
##          : $thread_number          => Number of BAM/CRAM compression threads
##          : $write_index            => Write index while writing file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_paths_ref;
    my $outfile_path;
    my $referencefile_path;
    my $region;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $thread_number;

    ## Default(s)
    my $force;
    my $output_format;
    my $write_index;

    my $tmpl = {
        filehandle => { store => \$filehandle, },
        force      => {
            allow       => [ undef, 0, 1, ],
            default     => 1,
            store       => \$force,
            strict_type => 1,
        },
        infile_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infile_paths_ref,
            strict_type => 1,
        },
        outfile_path  => { store => \$outfile_path, strict_type => 1, },
        output_format => {
            allow       => [qw{ bam cram sam }],
            default     => q{bam},
            store       => \$output_format,
            strict_type => 1,
        },
        referencefile_path => {
            store       => \$referencefile_path,
            strict_type => 1,
        },
        region => {
            store       => \$region,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => {
            strict_type => 1,
            store       => \$stdoutfile_path,
        },
        thread_number => {
            allow       => [ undef, qr{ \A\d+\z }xsm, ],
            store       => \$thread_number,
            strict_type => 1,
        },
        write_index => {
            allow       => [ undef, 0, 1, ],
            default     => 0,
            store       => \$write_index,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands =
      ( get_executable_base_command( { base_command => $BASE_COMMAND, } ), qw{ merge } );

    ## Samtools base args
    @commands = samtools_base(
        {
            commands_ref       => \@commands,
            output_format      => $output_format,
            referencefile_path => $referencefile_path,
            thread_number      => $thread_number,
            write_index        => $write_index,
        }
    );

    if ($force) {

        push @commands, q{-f};
    }

    if ($region) {

        push @commands, q{-R} . $SPACE . $region;
    }

    if ($outfile_path) {

        push @commands, $outfile_path;
    }

    push @commands, join $SPACE, @{$infile_paths_ref};

    # Redirect stderr output to program specific stderr file
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

sub samtools_sort {

## Function : Perl wrapper for writing samtools sort recipe to $filehandle. Based on samtools 1.10 (using htslib 1.10).
## Returns  : @commands
##          : $filehandle             => Sbatch filehandle to write to
##          : $infile_path            => Infile path
##          : $max_memory_per_thread  => Set maximum memory per thread; suffix K/M/G recognized
##          : $outfile_path           => Outfile path
##          : $output_format          => Output format
##          : $referencefile_path     => Reference file path (fasta)
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Stderrfile path append
##          : $stdoutfile_path        => Stdoutfile path
##          : $thread_number          => Number of BAM/CRAM compression threads
##          : $temp_file_path_prefix  => Write temporary files to path prefix
##          : $write_index            => Write index while writing file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $max_memory_per_thread;
    my $outfile_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $temp_file_path_prefix;
    my $thread_number;

    ## Default(s)
    my $output_format;
    my $write_index;

    my $tmpl = {
        filehandle  => { store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        max_memory_per_thread => {
            allow       => [ undef, qr{ \A \d+[G|K|M] \z }xsm, ],
            store       => \$max_memory_per_thread,
            strict_type => 1,
        },
        outfile_path  => { store => \$outfile_path, strict_type => 1, },
        output_format => {
            allow       => [qw{ bam cram sam }],
            default     => q{bam},
            store       => \$output_format,
            strict_type => 1,
        },
        referencefile_path => {
            store       => \$referencefile_path,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => {
            strict_type => 1,
            store       => \$stdoutfile_path,
        },
        temp_file_path_prefix => { store => \$temp_file_path_prefix, strict_type => 1, },
        thread_number         => {
            allow       => [ undef, qr{ \A\d+\z }xsm, ],
            store       => \$thread_number,
            strict_type => 1,
        },
        write_index => {
            allow       => [ undef, 0, 1, ],
            default     => 0,
            store       => \$write_index,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands =
      ( get_executable_base_command( { base_command => $BASE_COMMAND, } ), qw{ sort } );

    ## Samtools base args
    @commands = samtools_base(
        {
            commands_ref       => \@commands,
            output_format      => $output_format,
            referencefile_path => $referencefile_path,
            thread_number      => $thread_number,
            write_index        => $write_index,
        }
    );

    if ($max_memory_per_thread) {

        push @commands, q{-m} . $SPACE . $max_memory_per_thread;
    }
    if ($temp_file_path_prefix) {

        push @commands, q{-T} . $SPACE . $temp_file_path_prefix;
    }

    if ($outfile_path) {

        push @commands, q{-o} . $SPACE . $outfile_path;
    }

    push @commands, $infile_path;

    # Redirect stderr output to program specific stderr file
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

sub samtools_stats {

## Function : Perl wrapper for writing samtools stats recipe to $filehandle. Based on samtools 1.10 (using htslib 1.10).
## Returns  : @commands
##          : $auto_detect_input_format => Ignored (input format is auto-detected)
##          : $filehandle               => Sbatch filehandle to write to
##          : $infile_path              => Infile path
##          : $outfile_path             => Outfile path
##          : $regions_ref              => Regions to process {REF}
##          : $remove_overlap           => Remove overlaps of paired-end reads from coverage and base count computations
##          : $stderrfile_path          => Stderrfile path
##          : $stderrfile_path_append   => Stderrfile path append
##          : $stdoutfile_path          => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $outfile_path;
    my $regions_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $auto_detect_input_format;
    my $remove_overlap;

    my $tmpl = {
        auto_detect_input_format => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$auto_detect_input_format,
            strict_type => 1,
        },
        filehandle  => { store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path   => { store => \$outfile_path, strict_type => 1, },
        remove_overlap => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$remove_overlap,
            strict_type => 1,
        },
        regions_ref => { default => [], store => \$regions_ref, strict_type => 1, },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => {
            strict_type => 1,
            store       => \$stdoutfile_path,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands =
      ( get_executable_base_command( { base_command => $BASE_COMMAND, } ), qw{ stats } );

    if ($auto_detect_input_format) {

        push @commands, q{-s};
    }

    if ($remove_overlap) {

        push @commands, q{--remove-overlaps};
    }

    push @commands, $infile_path;

    if ( @{$regions_ref} ) {

        push @commands, join $SPACE, @{$regions_ref};
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

sub samtools_view {

## Function : Perl wrapper for writing samtools view recipe to $filehandle. Based on samtools 1.10 (using htslib 1.10).
## Returns  : "@commands"
##          : $auto_detect_input_format       => Ignored (input format is auto-detected)
##          : $exclude_reads_with_these_flags => Do not output alignments that match the bits set
##          : $filehandle                     => Sbatch filehandle to write to
##          : $fraction                       => Subsample the file to only a fraction of the alignments
##          : $infile_path                    => Infile path
##          : $outfile_path                   => Outfile path
##          : $output_format                  => Output format
##          : $referencefile_path             => Reference file path (fasta)
##          : $regions_ref                    => The regions to process {REF}
##          : $stderrfile_path                => Stderrfile path
##          : $stderrfile_path_append         => Stderrfile path append
##          : $thread_number                  => Number of BAM/CRAM compression threads
##          : $uncompressed_bam_output        => Uncompressed bam output
##          : $with_header                    => Include header
##          : $write_index                    => Write index while writing file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $exclude_reads_with_these_flags;
    my $filehandle;
    my $fraction;
    my $infile_path;
    my $outfile_path;
    my $referencefile_path;
    my $regions_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $thread_number;

    ## Default(s)
    my $auto_detect_input_format;
    my $output_format;
    my $uncompressed_bam_output;
    my $with_header;
    my $write_index;

    my $tmpl = {
        auto_detect_input_format => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$auto_detect_input_format,
            strict_type => 1,
        },
        exclude_reads_with_these_flags => {
            allow       => qr{ \A\d+\z }xsm,
            store       => \$exclude_reads_with_these_flags,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        fraction => {
            defined     => 1,
            store       => \$fraction,
            strict_type => 1,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path => {
            store       => \$outfile_path,
            strict_type => 1,
        },
        output_format => {
            allow       => [qw{ bam cram sam }],
            default     => q{bam},
            store       => \$output_format,
            strict_type => 1,
        },
        regions_ref => {
            default     => [],
            store       => \$regions_ref,
            strict_type => 1,
        },
        referencefile_path => {
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
        thread_number => {
            allow       => [ undef, qr{ \A\d+\z }xsm, ],
            store       => \$thread_number,
            strict_type => 1,
        },
        uncompressed_bam_output => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$uncompressed_bam_output,
            strict_type => 1,
        },
        with_header => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$with_header,
            strict_type => 1,
        },
        write_index => {
            allow       => [ undef, 0, 1, ],
            default     => 0,
            store       => \$write_index,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands =
      ( get_executable_base_command( { base_command => $BASE_COMMAND, } ), qw{ view } );

    ## Samtools base args
    @commands = samtools_base(
        {
            commands_ref       => \@commands,
            output_format      => $output_format,
            referencefile_path => $referencefile_path,
            thread_number      => $thread_number,
            write_index        => $write_index,
        }
    );

    if ($with_header) {

        push @commands, q{-h};
    }

    if ($auto_detect_input_format) {

        push @commands, q{-S};
    }

    if ($exclude_reads_with_these_flags) {
        push @commands, q{-F} . $SPACE . $exclude_reads_with_these_flags;
    }

    if ($outfile_path) {

        push @commands, q{-o} . $SPACE . $outfile_path;
    }

    if ($uncompressed_bam_output) {

        push @commands, q{-u};
    }

    if ($fraction) {

        push @commands, q{-s} . $SPACE . $fraction;
    }

    push @commands, $infile_path;

    if ( @{$regions_ref} ) {

        push @commands, join $SPACE, @{$regions_ref};
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

1;
