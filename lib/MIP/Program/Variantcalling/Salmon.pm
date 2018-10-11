package MIP::Program::Variantcalling::Salmon;

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
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ salmon_index salmon_quant };
}

## Constants
Readonly my $SPACE => q{ };

sub salmon_index {

## Function  : Perl wrapper for Salmon index, version 0.9.1.
## Returns   : @commands
## Arguments : $fasta_path             => Input reference fasta path, note salmon does not use the genome reference fasta, it uses a fasta file of transcripts
##           : $FILEHANDLE             => Filehandle to write to
##           : $outfile_path           => Outfile path
##           : $stderrfile_path        => Stderrfile path
##           : $stderrfile_path_append => Append stderr info to file path
##           : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $fasta_path;
    my $FILEHANDLE;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        fasta_path => {
            defined     => 1,
            required    => 1,
            store       => \$fasta_path,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
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
    my @commands = q{salmon index};

    # Options
    push @commands, q{--transcripts} . $SPACE . $fasta_path;

    push @commands, q{--index} . $SPACE . $outfile_path;

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

sub salmon_quant {

## Function  : Perl wrapper for Salmon quant, version 0.9.1.
## Returns   : @commands
## Arguments : $FILEHANDLE             => Filehandle to write to
##           : $gc_bias                => Correct for GC-bias
##           : $index_path             => Path to the index folder
##           : $lib                    => Library visit the salmon website for more  info
##           : $outfile_path           => The path of the  output directory
##           : $read_1_fastq_path      => Read 1 Fastq path
##           : $read_2_fastq_path      => Read 2 Fastq path
##           : $read_files_command     => command applied to the input FASTQ files
##           : $stderrfile_path        => Stderrfile path
##           : $stderrfile_path_append => Append stderr info to file path
##           : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $gc_bias;
    my $index_path;
    my $outfile_path;
    my $read_1_fastq_path;
    my $read_2_fastq_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $lib_type;
    my $read_files_command;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        gc_bias => {
            store       => \$gc_bias,
            strict_type => 1,
        },
        index_path => {
            defined     => 1,
            required    => 1,
            store       => \$index_path,
            strict_type => 1,
        },
        lib_type => {
            default     => q{ISF},
            store       => \$lib_type,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        read_1_fastq_path => {
            defined     => 1,
            required    => 1,
            store       => \$read_1_fastq_path,
            strict_type => 1,
        },
        read_2_fastq_path => {
            defined     => 1,
            store       => \$read_2_fastq_path,
            strict_type => 1,
        },
        read_files_command => {
            default     => q{pigz -dc},
            store       => \$read_files_command,
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
    my @commands = q{salmon quant};

    if ($gc_bias) {
        push @commands, q{--gcBias};
    }

    push @commands, q{--index} . $SPACE . $index_path;

# Library type, defines if the library is stranded or not, and the orientation of the reads, according to the documentation http://salmon.readthedocs.io/en/latest/library_type.html
    push @commands, q{--libType} . $SPACE . $lib_type;

    push @commands, q{--output} . $SPACE . $outfile_path;

# The input Fastq files, either single reads or paired. Salmon uses a bash command to stream the reads. Here, the default is <( pigz -dc file.fastq.gz )
    push @commands,
        q{-1}
      . $SPACE . q{<(}
      . $read_files_command
      . $SPACE
      . $read_1_fastq_path
      . $SPACE . q{)};

    if ($read_2_fastq_path) {
        push @commands,
            q{-2}
          . $SPACE . q{<(}
          . $read_files_command
          . $SPACE
          . $read_2_fastq_path
          . $SPACE . q{)};
    }

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
