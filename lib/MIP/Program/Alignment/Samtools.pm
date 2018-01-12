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
    our $VERSION = 1.03;

    # Inherit from Exporter to export functions and variables
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ samtools_view samtools_index samtools_stats samtools_mpileup samtools_faidx samtools_create_chromosome_files samtools_idxstats samtools_depth };

}

## Constants
Readonly my $SPACE     => q{ };
Readonly my $COMMA     => q{,};
Readonly my $AMPERSAND => q{&};

sub samtools_view {

## Function : Perl wrapper for writing samtools view recipe to $FILEHANDLE. Based on samtools 1.3.1 (using htslib 1.3.1).
## Returns  : "@commands"
##          : $auto_detect_input_format => Ignored (input format is auto-detected)
##          : $F_flag                   => Do not output alignments that match the bits set
##          : $FILEHANDLE               => Sbatch filehandle to write to
##          : $fraction                 => Subsample the file to only a fraction of the alignments
##          : $infile_path              => Infile path
##          : $outfile_path             => Outfile path
##          : $output_format            => Output format
##          : $regions_ref              => The regions to process {REF}
##          : $stderrfile_path          => Stderrfile path
##          : $stderrfile_path_append   => Stderrfile path append
##          : $thread_number            => Number of BAM/CRAM compression threads
##          : $uncompressed_bam_output  => Uncompressed bam output
##          : $with_header              => Include header

    my ($arg_href) = @_;

    ## Default(s)
    my $auto_detect_input_format;
    my $output_format;
    my $uncompressed_bam_output;
    my $with_header;

    ## Flatten argument(s)
    my $F_flag;
    my $FILEHANDLE;
    my $fraction;
    my $infile_path;
    my $outfile_path;
    my $regions_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $thread_number;

    my $tmpl = {
        auto_detect_input_format => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$auto_detect_input_format,
            strict_type => 1,
        },
        F_flag => {
            allow       => qr/^\d+$/,
            store       => \$F_flag,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
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
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        stderrfile_path_append => {
            store       => \$stderrfile_path_append,
            strict_type => 1,
        },
        thread_number => {
            allow       => qr/^\d+$/,
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
    my @commands = q{samtools view};

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

    if ($F_flag) {
        push @commands, q{-F} . $SPACE . $F_flag;
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
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    return @commands;
}

sub samtools_index {

## samtools_index

## Function : Perl wrapper for writing samtools index recipe to $FILEHANDLE. Based on samtools 1.3.1 (using htslib 1.3.1).
## Returns  : "@commands"
##          : $infile_path            => Infile path
##          : $stderrfile_path        => Stderrfile path
##          : $stdoutfile_path        => Stdoutfile path
##          : $FILEHANDLE             => Sbatch filehandle to write to
##          : $bai_format             => Generate BAI-format index for BAM files
##          : $stderrfile_path_append => Stderrfile path append

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $stdoutfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $bai_format;
    my $stderrfile_path_append;

    my $tmpl = {
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        FILEHANDLE => { store       => \$FILEHANDLE },
        bai_format => { strict_type => 1, store => \$bai_format },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Array @commands stores commands depending on input parameters
    my @commands = qw{ samtools index };

    ## Options
    if ($bai_format) {

        #Generate BAI-format index for BAM files
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    return @commands;
}

sub samtools_stats {

## samtools_stats

## Function : Perl wrapper for writing samtools stats recipe to $FILEHANDLE. Based on samtools 1.3.1 (using htslib 1.3.1).
## Returns  : "@commands"
##          : $regions_ref              => The regions to process {REF}
##          : $infile_path              => Infile path
##          : $outfile_path             => Outfile path
##          : $stderrfile_path          => Stderrfile path
##          : $FILEHANDLE               => Sbatch filehandle to write to
##          : $auto_detect_input_format => Ignored (input format is auto-detected)
##          : $stderrfile_path_append   => Stderrfile path append

    my ($arg_href) = @_;

    ## Default(s)
    my $auto_detect_input_format;

    ## Flatten argument(s)
    my $regions_ref;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $stderrfile_path_append;

    my $tmpl = {
        regions_ref =>
          { default => [], strict_type => 1, store => \$regions_ref },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        outfile_path    => { strict_type => 1, store => \$outfile_path },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        FILEHANDLE               => { store => \$FILEHANDLE },
        auto_detect_input_format => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$auto_detect_input_format
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Array @commands stores commands depending on input parameters
    my @commands = qw{ samtools stats };

    if ($auto_detect_input_format) {

        push @commands, q{-s};
    }

    ## Infile
    push @commands, $infile_path;

    if ( @{$regions_ref} ) {

        #Limit output to regions
        push @commands, join $SPACE, @{$regions_ref};
    }

    if ($outfile_path) {

        #Specify output filename
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    return @commands;
}

sub samtools_mpileup {

## samtools_mpileup

## Function : Perl wrapper for writing samtools mpileup recipe to $FILEHANDLE. Based on samtools 1.3.1 (using htslib 1.3.1).
## Returns  : "@commands"
##          : $infile_paths_ref                 => Infile paths {REF}
##          : $output_tags_ref                  => Optional tags to output {REF}
##          : $outfile_path                     => Outfile path
##          : $referencefile_path               => Reference sequence file
##          : $stderrfile_path                  => Stderrfile path
##          : $FILEHANDLE                       => Sbatch filehandle to write to
##          : $region                           => The regions to process {REF}
##          : $output_bcf                       => Generate genotype likelihoods in BCF format
##          : $per_sample_increased_sensitivity => Apply -m and -F per-sample for increased sensitivity
##          : $adjust_mq                        => Adjust mapping quality
##          : $stderrfile_path_append           => Stderrfile path append

    my ($arg_href) = @_;

    ## Default(s)
    my $per_sample_increased_sensitivity;
    my $adjust_mq;

    ## Flatten argument(s)
    my $infile_paths_ref;
    my $output_tags_ref;
    my $outfile_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $region;
    my $output_bcf;
    my $stderrfile_path_append;

    my $tmpl = {
        infile_paths_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$infile_paths_ref
        },
        output_tags_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$output_tags_ref
        },
        outfile_path       => { strict_type => 1, store => \$outfile_path },
        referencefile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$referencefile_path
        },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        FILEHANDLE => { store       => \$FILEHANDLE },
        region     => { strict_type => 1, store => \$region },
        output_bcf => { strict_type => 1, store => \$output_bcf },
        per_sample_increased_sensitivity => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$per_sample_increased_sensitivity
        },
        adjust_mq => {
            default     => 50,
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$adjust_mq
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Array @commands stores commands depending on input parameters
    my @commands = qw{ samtools mpileup };

    ## Options
    push @commands, q{--adjust-MQ} . $SPACE . $adjust_mq;

    if ($per_sample_increased_sensitivity) {

        push @commands, q{--per-sample-mF};
    }

    if ( @{$output_tags_ref} ) {

        push @commands, q{--output-tags} . $SPACE . join $COMMA,
          @{$output_tags_ref};
    }

    if ($region) {

        #Limit output to region
        push @commands, q{--region} . $SPACE . $region;
    }

    if ($referencefile_path) {

        #Reference sequence file
        push @commands, q{--fasta-ref} . $SPACE . $referencefile_path;
    }

    if ($output_bcf) {

        push @commands, q{--BCF};
    }

    if ($outfile_path) {

        #Specify output filename
        push @commands, q{--output} . $SPACE . $outfile_path;
    }

    ## Infile
    push @commands, join $SPACE, @{$infile_paths_ref};

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
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    return @commands;
}

sub samtools_faidx {

## samtools_faidx

## Function : Perl wrapper for writing samtools faidx recipe to $FILEHANDLE. Based on samtools 1.3.1 (using htslib 1.3.1).
## Returns  : "@commands"
##          : $regions_ref            => The regions to process {REF}
##          : $infile_path            => Infile path
##          : $outfile_path           => Outfile path
##          : $stderrfile_path        => Stderrfile path
##          : $FILEHANDLE             => Sbatch filehandle to write to
##          : $stderrfile_path_append => Stderrfile path append

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $regions_ref;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $stderrfile_path_append;

    my $tmpl = {
        regions_ref =>
          { default => [], strict_type => 1, store => \$regions_ref },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        outfile_path    => { strict_type => 1, store => \$outfile_path },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        FILEHANDLE => { store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Array @commands stores commands depending on input parameters
    my @commands = qw{ samtools faidx };

    ## Infile
    push @commands, $infile_path;

    if ( @{$regions_ref} ) {

        #Limit output to regions
        push @commands, join $SPACE, @{$regions_ref};
    }

    if ($outfile_path) {

        #Specify output filename
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    return @commands;
}

sub samtools_create_chromosome_files {

## Function : Perl wrapper for writing chromosome files used by other scripts. Writes to FILEHANDLE.
## Returns  :
##          : $regions_ref        => The regions to process {REF}
##          : $infile_path        => Infile path
##          : temp_directory      => Temporary directory
##          : $outfile_path       => Outfile path
##          : $suffix             => suffix to append to outfile_path
##          : $max_process_number => Max number of processeses
##          : $FILEHANDLE         => Sbatch filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $regions_ref;
    my $infile_path;
    my $temp_directory;
    my $suffix;
    my $max_process_number;
    my $FILEHANDLE;

    my $tmpl = {
        regions_ref => {
            required    => 1,
            default     => [],
            strict_type => 1,
            store       => \$regions_ref,
        },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path,
        },
        temp_directory => {
            required    => 1,
            strict_type => 1,
            store       => \$temp_directory,
        },
        suffix             => { strict_type => 1, store => \$suffix, },
        max_process_number => {
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$max_process_number,
        },
        FILEHANDLE => { store => \$FILEHANDLE, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $process_batches_count = 1;

    while ( my ( $contig_index, $contig ) = each @{$regions_ref} ) {
        $process_batches_count = print_wait(
            {
                process_counter       => $contig_index,
                max_process_number    => $max_process_number,
                process_batches_count => $process_batches_count,
                FILEHANDLE            => $FILEHANDLE,
            }
        );

        samtools_faidx(
            {
                regions_ref  => [$contig],
                infile_path  => $infile_path,
                outfile_path => catfile( $temp_directory, $contig . $suffix ),
                FILEHANDLE   => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $SPACE . $AMPERSAND;

    }
    return;
}

sub samtools_idxstats {

## Function : Perl wrapper for writing samtools idxstats recipe to $FILEHANDLE. Based on samtools 1.6 (using htslib 1.6).
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $infile_path            => Infile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
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
    my @commands = q{samtools idxstats};

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
            FILEHANDLE   => $FILEHANDLE,
            commands_ref => \@commands,
            separator    => $SPACE,

        }
    );
    return @commands;
}

sub samtools_depth {

## Function : Perl wrapper for writing samtools depth recipe to $FILEHANDLE. Based on samtools 1.6 (using htslib 1.6).
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $infile_path            => Infile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
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
    my @commands = q{samtools depth};

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
            FILEHANDLE   => $FILEHANDLE,
            commands_ref => \@commands,
            separator    => $SPACE,

        }
    );
    return @commands;
}

1;
