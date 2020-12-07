package MIP::Program::Salmon;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $SPACE };
use MIP::Environment::Executable qw{ get_executable_base_command };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.07;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ salmon_index salmon_quant };
}

Readonly my $BASE_COMMAND => q{salmon};

sub salmon_index {

## Function  : Perl wrapper for Salmon index, version 0.9.1.
## Returns   : @commands
## Arguments : $fasta_path             => Input reference fasta path, note salmon does not use the genome reference fasta, it uses a fasta file of transcripts
##           : $filehandle             => Filehandle to write to
##           : $outfile_path           => Outfile path
##           : $stderrfile_path        => Stderrfile path
##           : $stderrfile_path_append => Append stderr info to file path
##           : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $fasta_path;
    my $filehandle;
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
        filehandle => {
            store => \$filehandle,
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

    my @commands =
      ( get_executable_base_command( { base_command => $BASE_COMMAND, } ), qw{ index } );

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
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,

        }
    );
    return @commands;
}

sub salmon_quant {

## Function  : Perl wrapper for Salmon quant, version 0.9.1.
## Returns   : @commands
## Arguments : $filehandle             => Filehandle to write to
##           : $gc_bias                => Correct for GC-bias
##           : $index_path             => Path to the index folder
##           : $libi_type              => Library visit the salmon website for more  info
##           : $outdir_path            => Path of the output directory
##           : $read_1_fastq_paths_ref => Read 1 Fastq paths
##           : $read_2_fastq_paths_ref => Read 2 Fastq paths
##           : $read_files_command     => command applied to the input FASTQ files
##           : $stderrfile_path        => Stderrfile path
##           : $stderrfile_path_append => Append stderr info to file path
##           : $stdoutfile_path        => Stdoutfile path
##           : $validate_mappings      => Validate mappings

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $gc_bias;
    my $index_path;
    my $outdir_path;
    my $read_1_fastq_paths_ref;
    my $read_2_fastq_paths_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $lib_type;
    my $read_files_command;
    my $validate_mappings;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
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
            allow       => [qw{ A ISF ISR MSF MSR OSR OSF }],
            default     => q{A},
            store       => \$lib_type,
            strict_type => 1,
        },
        outdir_path => {
            defined     => 1,
            required    => 1,
            store       => \$outdir_path,
            strict_type => 1,
        },
        read_1_fastq_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$read_1_fastq_paths_ref,
            strict_type => 1,
        },
        read_2_fastq_paths_ref => {
            default     => [],
            defined     => 1,
            store       => \$read_2_fastq_paths_ref,
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
        validate_mappings => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$validate_mappings,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands =
      ( get_executable_base_command( { base_command => $BASE_COMMAND, } ), qw{ quant } );

    if ($gc_bias) {
        push @commands, q{--gcBias};
    }

    push @commands, q{--index} . $SPACE . $index_path;

# Library type, defines if the library is stranded or not, and the orientation of the reads, according to the documentation http://salmon.readthedocs.io/en/latest/library_type.html
    push @commands, q{--libType} . $SPACE . $lib_type;

    push @commands, q{--output} . $SPACE . $outdir_path;

    if ($validate_mappings) {
        push @commands, q{--validateMappings};
    }

# The input Fastq files, either single reads or paired. Salmon uses a bash command to stream the reads. Here, the default is <( pigz -dc file.fastq.gz )
    push @commands,
        q{-1}
      . $SPACE . q{<(}
      . $read_files_command
      . $SPACE
      . join( $SPACE, @{$read_1_fastq_paths_ref} )
      . $SPACE . q{)};

    if ( @{$read_2_fastq_paths_ref} ) {
        push @commands,
            q{-2}
          . $SPACE . q{<(}
          . $read_files_command
          . $SPACE
          . join( $SPACE, @{$read_2_fastq_paths_ref} )
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
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,

        }
    );
    return @commands;
}

1;
