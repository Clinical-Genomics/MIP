package MIP::Program::Alignment::Bwa;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };

## CPANM
use Readonly;

## MIPs lib/
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ bwa_mem run_bwamem bwa_index };

}

## Constants
Readonly my $SPACE => q{ };

sub bwa_mem {

## Function : Perl wrapper for writing bwa mem recipe to $FILEHANDLE. Based on bwa 0.7.15-r1140.
## Returns  : @commands
## Arguments: $infile_path             => Infile path (read 1 or interleaved i.e. read 1 and 2)
##          : $idxbase                 => Idxbase (human genome references and bwa mem idx files)
##          : $second_infile_path      => Second infile path (read 2)
##          : $stdoutfile_path         => Stdoutfile path
##          : $stderrfile_path         => Stderrfile path
##          : $stderrfile_path_append  => Stderrfile path append
##          : $FILEHANDLE              => Sbatch filehandle to write to
##          : $thread_number           => Number of threads
##          : $read_group_header       => Read group header line, such as '@RG\tID:foo\tSM:bar'
##          : $mark_split_as_secondary => Mark shorter split hits as secondary
##          : $interleaved_fastq_file  => Smart pairing

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $idxbase;
    my $second_infile_path;
    my $stdoutfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;
    my $thread_number;
    my $read_group_header;

    ## Default(s)
    my $mark_split_as_secondary;
    my $interleaved_fastq_file;

    my $tmpl = {
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        idxbase =>
          { required => 1, defined => 1, strict_type => 1, store => \$idxbase },
        second_infile_path =>
          { strict_type => 1, store => \$second_infile_path },
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        FILEHANDLE    => { store => \$FILEHANDLE },
        thread_number => {
            allow       => qr/ ^\d+$ /xms,
            strict_type => 1,
            store       => \$thread_number
        },
        read_group_header => { strict_type => 1, store => \$read_group_header },
        mark_split_as_secondary => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$mark_split_as_secondary
        },
        interleaved_fastq_file => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$interleaved_fastq_file
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Bwa
    # Stores commands depending on input parameters
    my @commands = q{bwa mem};

    ##Options
    if ($thread_number) {

        # Number of threads
        push @commands, q{-t} . $SPACE . $thread_number;
    }
    if ($interleaved_fastq_file) {

        # Interleaved fastq file
        push @commands, q{-p};
    }
    if ($mark_split_as_secondary) {

        # Mark shorter split hits as secondary
        push @commands, q{-M};
    }
    if ($read_group_header) {

        # Read group header line
        push @commands, q{-R} . $SPACE . $read_group_header;
    }

    ## Human reference genome and bwa mem files
    push @commands, $idxbase;

    ## Infile
    push @commands, $infile_path;

    ## Read 2
    if ($second_infile_path) {

        push @commands, $second_infile_path;
    }

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

sub run_bwamem {

## Function : Perl wrapper for writing run_bwamem recipe to $FILEHANDLE. Based on bwakit 0.7.12.
## Returns  : @commands
## Arguments: $infile_path            => Infile path (read 1 or interleaved i.e. read 1 and 2)
##          : $idxbase                => Idxbase (human genome references and bwa mem idx files)
##          : $outfiles_prefix_path   => Prefix for output files
##          : $second_infile_path     => Second infile path (read 2)
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Stderrfile path append
##          : $FILEHANDLE             => Sbatch filehandle to write to
##          : $thread_number          => Number of threads
##          : $read_group_header      => Read group header line, such as '@RG\tID:foo\tSM:bar'
##          : $hla_typing             => Apply HLA typing

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $idxbase;
    my $outfiles_prefix_path;
    my $second_infile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;
    my $thread_number;
    my $read_group_header;

    ## Default(s)
    my $hla_typing;

    my $tmpl = {
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        idxbase =>
          { required => 1, defined => 1, strict_type => 1, store => \$idxbase },
        outfiles_prefix_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfiles_prefix_path
        },
        second_infile_path =>
          { strict_type => 1, store => \$second_infile_path },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        FILEHANDLE    => { store => \$FILEHANDLE },
        thread_number => {
            allow       => qr/ ^\d+$ /xms,
            strict_type => 1,
            store       => \$thread_number
        },
        read_group_header => { strict_type => 1, store => \$read_group_header },
        hla_typing        => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$hla_typing
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Run_bwamem
    # Stores commands depending on input parameters
    my @commands = q{run-bwamem};

    ## Options
    if ($thread_number) {

        # Number of threads
        push @commands, q{-t} . $SPACE . $thread_number;
    }
    if ($hla_typing) {

        # Apply HLA typing
        push @commands, q{-H};
    }
    if ($read_group_header) {

        # Read group header line
        push @commands, q{-R} . $SPACE . $read_group_header;
    }

    # Outfiles prefix
    push @commands, q{-o} . $SPACE . $outfiles_prefix_path;

    ## Human reference genome and bwa mem files
    push @commands, $idxbase;

    ## Infile
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;
}

sub bwa_index {

## Function : Perl wrapper for writing bwa mem recipe to $FILEHANDLE. Based on bwa 0.7.15-r1140.
## Returns  : @commands
## Arguments: $prefix                 => Prefix of the index [same as fasta name]
##          : $construction_algorithm => BWT construction algorithm: bwtsw, is or rb2 [auto]
##          : $reference_genome       => Reference genome [.fasta]
##          : $stdoutfile_path        => Stdoutfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $FILEHANDLE             => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $prefix;
    my $reference_genome;
    my $construction_algorithm;
    my $stdoutfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;

    ## Default(s)

    my $tmpl = {
        prefix =>
          { required => 1, defined => 1, strict_type => 1, store => \$prefix },
        reference_genome => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$reference_genome
        },
        construction_algorithm => {
            allow       => [qw{ bwtsw rb2 }],
            strict_type => 1,
            store       => \$construction_algorithm
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

    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = q{bwa index};

    push @commands, q{-p} . $SPACE . $prefix;

    if ($construction_algorithm) {

        push @commands, q{-a} . $SPACE . $construction_algorithm;
    }

    push @commands, $reference_genome;

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
1;
