package MIP::Program::Alignment::Samtools;

use Carp;
use charnames qw{ :full :short };
use open qw{ :encoding(UTF-8) :std };
use strict;
use utf8;    #Allow unicode characters in this script
use warnings qw{ FATAL utf8 };
use warnings;

## CPANM
use File::Spec::Functions qw{ catdir catfile };
use Params::Check qw{ allow check last_error };
use Readonly;

## MIPs lib/
use MIP::Processmanagement::Processes qw{ print_wait };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.06;

    # Inherit from Exporter to export functions and variables
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      samtools_create_chromosome_files
      samtools_depth
      samtools_faidx
      samtools_idxstats
      samtools_index
      samtools_stats
      samtools_view };

}

## Constants
Readonly my $SPACE     => q{ };
Readonly my $COMMA     => q{,};
Readonly my $AMPERSAND => q{&};

sub samtools_view {

## Function : Perl wrapper for writing samtools view recipe to $filehandle. Based on samtools 1.9 (using htslib 1.9).
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
            allow       => [qw{ sam bam cram json }],
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
            allow       => qr{ \A\d+\z }xsm,
            store       => \$thread_number,
            strict_type => 1,
        },
        with_header => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$with_header,
            strict_type => 1,
        },
        uncompressed_bam_output => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$uncompressed_bam_output,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Array @commands stores commands depending on input parameters
    my @commands = qw{ samtools view };

    ## Options
    if ($thread_number) {

        #Number of threads
        push @commands, q{--threads} . $SPACE . $thread_number;
    }

    if ($with_header) {

        #Include header
        push @commands, q{-h};
    }

    if ($output_format) {

        #Output format
        push @commands, q{--output-fmt} . $SPACE . uc $output_format;
    }

    if ($auto_detect_input_format) {

        push @commands, q{-S};
    }

    if ($exclude_reads_with_these_flags) {
        push @commands, q{-F} . $SPACE . $exclude_reads_with_these_flags;
    }

    if ($referencefile_path) {

        push @commands, q{--reference} . $SPACE . $referencefile_path;
    }
    if ($outfile_path) {

        #Specify output filename
        push @commands, q{-o} . $SPACE . $outfile_path;
    }

    if ($uncompressed_bam_output) {

        push @commands, q{-u};
    }

    if ($fraction) {

        push @commands, q{-s} . $SPACE . $fraction;

    }

    ## Infile
    push @commands, $infile_path;

    if ( @{$regions_ref} ) {

        #Limit output to regions
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
            separator    => $SPACE,
            filehandle   => $filehandle,
        }
    );
    return @commands;
}

sub samtools_index {

## Function : Perl wrapper for writing samtools index recipe to $filehandle. Based on samtools 1.9 (using htslib 1.9).
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

    ## Array @commands stores commands depending on input parameters
    my @commands = qw{ samtools index };

    ## Options
    if ($bai_format) {

        # Generate BAI-format index for BAM files
        push @commands, q{-b};
    }

    ## Infile
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

## Function : Perl wrapper for writing samtools stats recipe to $filehandle. Based on samtools 1.9 (using htslib 1.9).
## Returns  : @commands
##          : $auto_detect_input_format => Ignored (input format is auto-detected)
##          : $filehandle               => Sbatch filehandle to write to
##          : $infile_path              => Infile path
##          : $outfile_path             => Outfile path
##          : $regions_ref              => Regions to process {REF}
##          : $remove_overlap           => Remove overlaps of paired-end reads from coverage and base count computations
##          : $stderrfile_path          => Stderrfile path
##          : $stderrfile_path_append   => Stderrfile path append
##          : $stdoutfile_path        => Stdoutfile path

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

    ## Array @commands stores commands depending on input parameters
    my @commands = qw{ samtools stats };

    if ($auto_detect_input_format) {

        push @commands, q{-s};
    }

    if ($remove_overlap) {

        push @commands, q{--remove-overlaps};
    }

    ## Infile
    push @commands, $infile_path;

    if ( @{$regions_ref} ) {

        # Limit output to regions
        push @commands, join $SPACE, @{$regions_ref};
    }

    if ($outfile_path) {

        # Specify output filename
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
            filehandle   => $filehandle,
            commands_ref => \@commands,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub samtools_faidx {

## Function : Perl wrapper for writing samtools faidx recipe to $filehandle. Based on samtools 1.9 (using htslib 1.9).
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

    ## Array @commands stores commands depending on input parameters
    my @commands = qw{ samtools faidx };

    ## Infile
    push @commands, $infile_path;

    if ( @{$regions_ref} ) {

        # Limit output to regions
        push @commands, join $SPACE, @{$regions_ref};
    }

    if ($outfile_path) {

        # Specify output filename
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

sub samtools_idxstats {

## Function : Perl wrapper for writing samtools idxstats recipe to $filehandle. Based on samtools 1.6 (using htslib 1.6).
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

    ## Stores commands depending on input parameters
    my @commands = qw{ samtools idxstats };

    ## Infile
    push @commands, $infile_path;

    ## Redirect stderr output to program specific stderr file
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

sub samtools_depth {

## Function : Perl wrapper for writing samtools depth recipe to $filehandle. Based on samtools 1.6 (using htslib 1.6).
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

    ## Stores commands depending on input parameters
    my @commands = qw{ samtools depth };

    ## Optionally set the read depth cutoff value
    if ($max_depth_treshold) {
        push @commands, q{-d} . $SPACE . $max_depth_treshold;
    }

    ## Infile
    push @commands, $infile_path;

    ## Redirect stderr output to program specific stderr file
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

1;
