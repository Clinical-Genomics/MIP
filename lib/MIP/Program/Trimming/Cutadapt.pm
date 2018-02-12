package MIP::Program::Trimming::Cutadapt;

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
use Readonly;

## MIPs lib/
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ cutadapt };
}

## Constants
Readonly my $SPACE => q{ };

sub cutadapt {

## Function : Perl wrapper for cutadapt.
## Returns  : @commands
## Arguments: $adapter_3_prime        => Sequence of an adapter ligated to the 3' end (paired data: of the first read)
##          : $adapter_5_prime        => Sequence of an adapter ligated to the 5' end (paired data: of the first read)
##          : $adapter_3_prime_second => 3' adapter to be removed from second read in a pair
##          : $adapter_5_prime_second => 5' adapter to be removed from second read in a pair
##          : $FILEHANDLE             => Filehandle to write to
##          : $infile_path            => Fastq file for read1
##          : $infile_path_second     => Fastq file for read1
##          : $outputfile_path        => Write trimmed reads to FILE
##          : $paired_filter          => Which of the reads in a paired-end read have to match the filtering criterion (any or both)
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $thread_numbter         => Number of CPU cores to use
##          : $min_trimmed_length     => Discard reads shorter than LENGTH

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $adapter_3_prime;
    my $adapter_5_prime;
    my $adapter_3_prime_second;
    my $adapter_5_prime_second;
    my $FILEHANDLE;
    my $infile_path;
    my $infile_path_second;
    my $min_trimmed_length;
    my $outputfile_path;
    my $paired_filter;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $thread_number;

    my $tmpl = {
        adapter_3_prime => {
            allow       => [ undef, qr/ ^\w+$ /xsm ],
            default     => undef,
            store       => \$adapter_3_prime,
            strict_type => 1,
        },
        adapter_5_prime => {
            allow       => [ undef, qr/ ^\w+$ /xsm ],
            default     => undef,
            store       => \$adapter_5_prime,
            strict_type => 1,
        },
        adapter_3_prime_second => {
            allow       => [ undef, qr/ ^\w+$ /xsm ],
            default     => undef,
            store       => \$adapter_3_prime_second,
            strict_type => 1,
        },
        adapter_5_prime_second => {
            allow       => [ undef, qr/ ^\w+$ /xsm ],
            default     => undef,
            store       => \$adapter_5_prime_second,
            strict_type => 1,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        infile_path_second => {
            store       => \$infile_path_second,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        min_trimmed_length => {
            allow       => qr/ ^\d+$ /xms,
            default     => 2,
            store       => \$min_trimmed_length,
            strict_type => 1,
        },
        outputfile_path => {
            store       => \$outputfile_path,
            strict_type => 1,
        },
        paired_filter => {
            allow       => [qw{ any both }],
            default     => qw{ both },
            store       => \$paired_filter,
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
        thread_number => {
            allow       => qr/ ^\d+$ /xms,
            default     => 2,
            store       => \$thread_number,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = q{cutadapt};

    push @commands, q{--cores=} . $thread_number;

    push @commands, q{--minimum-length=} . $min_trimmed_length;

    if ($adapter_3_prime) {

        push @commands, q{--adapter=} . $adapter_3_prime;
    }
    if ($adapter_5_prime) {

        push @commands, q{--front=} . $adapter_5_prime;
    }

    # conditions for paired-end reads
    if ($infile_path_second) {

        push @commands, q{--pair-filter=} . $paired_filter;
    }

    if ($infile_path_second && $adapter_3_prime_second ) {

        push @commands, q{-A} . $SPACE . $adapter_3_prime_second;
    }

    if ($infile_path_second && $adapter_5_prime_second ) {

        push @commands, q{-G} . $SPACE . $adapter_5_prime_second;
    }

    # input fastq file read1
    push @commands, $infile_path;

    # input fastq file read2
    if ($infile_path_second) {

        push @commands, $infile_path_second;
    }

    if ($outputfile_path) {

        push @commands, q{--output=} . $outputfile_path;
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
