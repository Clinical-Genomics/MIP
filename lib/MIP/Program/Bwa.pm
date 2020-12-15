package MIP::Program::Bwa;

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
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $SPACE };
use MIP::Environment::Executable qw{ get_executable_base_command };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ bwa_index bwa_mem bwa_mem2_mem bwa_mem2_index run_bwamem };

}

Readonly my $BWA_BASE_COMMAND    => q{bwa};
Readonly my $BWA2_BASE_COMMAND   => q{bwa-mem2};
Readonly my $RUNBWA_BASE_COMMAND => q{run-bwamem};

sub bwa_index {

## Function : Perl wrapper for writing bwa mem recipe to $filehandle. Based on bwa 0.7.15-r1140.
## Returns  : @commands
## Arguments: $construction_algorithm => BWT construction algorithm: bwtsw, is or rb2 [auto]
##          : $filehandle             => Filehandle to write to
##          : $prefix                 => Prefix of the index [same as fasta name]
##          : $reference_genome       => Reference genome [.fasta]
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $construction_algorithm;
    my $filehandle;
    my $prefix;
    my $reference_genome;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)

    my $tmpl = {
        construction_algorithm => {
            allow       => [qw{ bwtsw rb2 }],
            store       => \$construction_algorithm,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        prefix => { defined => 1, required => 1, store => \$prefix, strict_type => 1, },
        reference_genome => {
            defined     => 1,
            required    => 1,
            store       => \$reference_genome,
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
        get_executable_base_command( { base_command => $BWA_BASE_COMMAND, } ),
        qw{ index }
    );

    push @commands, q{-p} . $SPACE . $prefix;

    if ($construction_algorithm) {

        push @commands, q{-a} . $SPACE . $construction_algorithm;
    }

    push @commands, $reference_genome;

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

sub bwa_mem {

## Function : Perl wrapper for writing bwa mem recipe to $filehandle. Based on bwa 0.7.15-r1140.
## Returns  : @commands
## Arguments: $filehandle              => Sbatch filehandle to write to
##          : $idxbase                 => Idxbase (human genome references and bwa mem idx files)
##          : $infile_path             => Infile path (read 1 or interleaved i.e. read 1 and 2)
##          : $interleaved_fastq_file  => Smart pairing
##          : $mark_split_as_secondary => Mark shorter split hits as secondary
##          : $read_group_header       => Read group header line, such as '@RG\tID:foo\tSM:bar'
##          : $second_infile_path      => Second infile path (read 2)
##          : $soft_clip_sup_align     => Use soft clipping for supplementary alignments
##          : $stderrfile_path         => Stderrfile path
##          : $stderrfile_path_append  => Stderrfile path append
##          : $stdoutfile_path         => Stdoutfile path
##          : $thread_number           => Number of threads

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $idxbase;
    my $infile_path;
    my $read_group_header;
    my $second_infile_path;
    my $soft_clip_sup_align;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $thread_number;

    ## Default(s)
    my $interleaved_fastq_file;
    my $mark_split_as_secondary;

    my $tmpl = {
        filehandle => { store => \$filehandle, },
        idxbase    => {
            defined     => 1,
            required    => 1,
            store       => \$idxbase,
            strict_type => 1,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        interleaved_fastq_file => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$interleaved_fastq_file,
            strict_type => 1,
        },
        mark_split_as_secondary => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$mark_split_as_secondary,
            strict_type => 1,
        },
        read_group_header   => { store => \$read_group_header,  strict_type => 1, },
        second_infile_path  => { store => \$second_infile_path, strict_type => 1, },
        soft_clip_sup_align => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$soft_clip_sup_align,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
        thread_number   => {
            allow       => qr/ ^\d+$ /xms,
            store       => \$thread_number,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands =
      ( get_executable_base_command( { base_command => $BWA_BASE_COMMAND, } ),
        qw{ mem } );

    if ($thread_number) {

        push @commands, q{-t} . $SPACE . $thread_number;
    }
    if ($interleaved_fastq_file) {

        push @commands, q{-p};
    }
    if ($mark_split_as_secondary) {

        push @commands, q{-M};
    }
    if ($soft_clip_sup_align) {

        push @commands, q{-Y};
    }
    if ($read_group_header) {

        push @commands, q{-R} . $SPACE . $read_group_header;
    }

    ## Human reference genome and bwa mem files
    push @commands, $idxbase;

    push @commands, $infile_path;

    ## Read 2
    if ($second_infile_path) {

        push @commands, $second_infile_path;
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

sub bwa_mem2_mem {

## Function : Perl wrapper for writing bwa mem 2 mem recipe to $filehandle. Based on bwa mem 2 2.0.
## Returns  : @commands
## Arguments: $filehandle              => Sbatch filehandle to write to
##          : $idxbase                 => Idxbase (human genome references and bwa mem idx files)
##          : $infile_path             => Infile path (read 1 or interleaved i.e. read 1 and 2)
##          : $interleaved_fastq_file  => Smart pairing
##          : $mark_split_as_secondary => Mark shorter split hits as secondary
##          : $read_group_header       => Read group header line, such as '@RG\tID:foo\tSM:bar'
##          : $second_infile_path      => Second infile path (read 2)
##          : $soft_clip_sup_align     => Use soft clipping for supplementary alignments
##          : $stderrfile_path         => Stderrfile path
##          : $stderrfile_path_append  => Stderrfile path append
##          : $stdoutfile_path         => Stdoutfile path
##          : $thread_number           => Number of threads

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $idxbase;
    my $infile_path;
    my $read_group_header;
    my $second_infile_path;
    my $soft_clip_sup_align;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $thread_number;

    ## Default(s)
    my $interleaved_fastq_file;
    my $mark_split_as_secondary;

    my $tmpl = {
        filehandle => { store => \$filehandle, },
        idxbase    => {
            defined     => 1,
            required    => 1,
            store       => \$idxbase,
            strict_type => 1,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        interleaved_fastq_file => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$interleaved_fastq_file,
            strict_type => 1,
        },
        mark_split_as_secondary => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$mark_split_as_secondary,
            strict_type => 1,
        },
        read_group_header   => { store => \$read_group_header,  strict_type => 1, },
        second_infile_path  => { store => \$second_infile_path, strict_type => 1, },
        soft_clip_sup_align => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$soft_clip_sup_align,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
        thread_number   => {
            allow       => qr/ ^\d+$ /xms,
            store       => \$thread_number,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = (
        get_executable_base_command( { base_command => $BWA2_BASE_COMMAND, } ),
        qw{ mem }
    );

    if ($thread_number) {

        push @commands, q{-t} . $SPACE . $thread_number;
    }
    if ($interleaved_fastq_file) {

        push @commands, q{-p};
    }
    if ($mark_split_as_secondary) {

        push @commands, q{-M};
    }
    if ($soft_clip_sup_align) {

        push @commands, q{-Y};
    }
    if ($read_group_header) {

        push @commands, q{-R} . $SPACE . $read_group_header;
    }

    ## Human reference genome and bwa mem files
    push @commands, $idxbase;

    push @commands, $infile_path;

    ## Read 2
    if ($second_infile_path) {

        push @commands, $second_infile_path;
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

sub bwa_mem2_index {

## Function : Perl wrapper for writing bwa mem 2 index recipe to $filehandle. Based on bwa mem 2 2.0.
## Returns  : @commands
## Arguments: $filehandle             => Filehandle to write to
##          : $prefix                 => Prefix of the index [same as fasta name]
##          : $reference_genome       => Reference genome [.fasta]
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $prefix;
    my $reference_genome;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        prefix => { defined => 1, required => 1, store => \$prefix, strict_type => 1, },
        reference_genome => {
            defined     => 1,
            required    => 1,
            store       => \$reference_genome,
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
        get_executable_base_command( { base_command => $BWA2_BASE_COMMAND, } ),
        qw{ index }
    );

    push @commands, q{-p} . $SPACE . $prefix;

    push @commands, $reference_genome;

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

sub run_bwamem {

## Function : Perl wrapper for writing run_bwamem recipe to $filehandle. Based on bwakit 0.7.12.
## Returns  : @commands
## Arguments: $filehandle             => Sbatch filehandle to write to
##          : $hla_typing             => Apply HLA typing
##          : $idxbase                => Idxbase (human genome references and bwa mem idx files)
##          : $infile_path            => Infile path (read 1 or interleaved i.e. read 1 and 2)
##          : $outfiles_prefix_path   => Prefix for output files
##          : $read_group_header      => Read group header line, such as '@RG\tID:foo\tSM:bar'
##          : $second_infile_path     => Second infile path (read 2)
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Stderrfile path append
##          : $thread_number          => Number of threads

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $idxbase;
    my $infile_path;
    my $outfiles_prefix_path;
    my $read_group_header;
    my $second_infile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $thread_number;

    ## Default(s)
    my $hla_typing;

    my $tmpl = {
        filehandle => { store => \$filehandle, },
        hla_typing => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$hla_typing,
            strict_type => 1,
        },
        idxbase => {
            defined     => 1,
            required    => 1,
            store       => \$idxbase,
            strict_type => 1,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfiles_prefix_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfiles_prefix_path,
            strict_type => 1,
        },
        read_group_header  => { store => \$read_group_header,  strict_type => 1, },
        second_infile_path => { store => \$second_infile_path, strict_type => 1, },
        stderrfile_path    => { store => \$stderrfile_path,    strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        thread_number => {
            allow       => qr/ ^\d+$ /xms,
            store       => \$thread_number,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands =
      ( get_executable_base_command( { base_command => $RUNBWA_BASE_COMMAND, } ), );

    if ($thread_number) {

        push @commands, q{-t} . $SPACE . $thread_number;
    }
    if ($hla_typing) {

        push @commands, q{-H};
    }
    if ($read_group_header) {

        push @commands, q{-R} . $SPACE . $read_group_header;
    }

    push @commands, q{-o} . $SPACE . $outfiles_prefix_path;

    ## Human reference genome and bwa mem files
    push @commands, $idxbase;

    push @commands, $infile_path;

    ## Read 2
    if ($second_infile_path) {

        push @commands, $second_infile_path;
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
