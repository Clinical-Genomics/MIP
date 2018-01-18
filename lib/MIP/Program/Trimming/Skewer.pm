package MIP::Program::Trimming::Skewer;

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
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ skewer };
}

## Constants
Readonly my $SPACE => q{ };

sub skewer {

## Function : Perl wrapper for skewer, an adapter trimmer. Based on Version 0.2.2 (updated in April 4, 2016).
## Returns  : @commands
## Arguments:
##          : $adapter_sequence        => Adapter sequence pair/single-end reads
##          : $adapter_sequence_second => Adapter sequence for read2
##          : $compress_output         => Compress output in GZIP format
##          : $FILEHANDLE              => Filehandle to write to
##          : $infile_path             => Fastq file for read1
##          : $outsuffix               => Base name of output file
##          : $second_infile_path      => Fastq file for read2
##          : $stderrfile_path         => Stderrfile path
##          : $stderrfile_path_append  => Append stderr info to file path
##          : $stdoutfile_path         => Stdoutfile path
##          : $thread_number           => Number of threads
##          : $trim_mode               => Trimming mode: SE: head(5p), tail(3p), any. PE: pe(pairedend), mp(matepair), ap(amplicon).

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $adapter_sequence;
    my $adapter_sequence_second;
    my $compress_output;
    my $FILEHANDLE;
    my $infile_path;
    my $outsuffix;
    my $second_infile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $thread_number;
    my $trim_mode;

    my $tmpl = {
        adapter_sequence => {
            store       => \$adapter_sequence,
            strict_type => 1,
        },
        adapter_sequence_second => {
            store       => \$adapter_sequence_second,
            strict_type => 1,
        },
        compress_output => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$compress_output,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        infile_path => {
            required    => 1,
            defined     => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outsuffix => {
            required    => 1,
            store       => \$outsuffix,
            strict_type => 1,
        },
        second_infile_path => {
            store       => \$second_infile_path,
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
        trim_mode => {
            allow       => [qw{ any ap head mp pe tail }],
            default     => qw{ pe },
            store       => \$trim_mode,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = q{skewer};

    push @commands, q{-x} . $SPACE . $adapter_sequence;

    # If a second read is provided, trim adapter accordingly
    if ($second_infile_path) {
        push @commands, q{-y} . $SPACE . $adapter_sequence_second;
    }

    if ($compress_output) {

     # Compresses output with builtin --rsyncable. Not compatible with Mac OS X.
        push @commands, q{-z} . $SPACE;
    }

    push @commands, q{--mode} . $SPACE . $trim_mode;

    push @commands, q{--output} . $SPACE . $outsuffix;

    push @commands, q{--threads} . $SPACE . $thread_number;

    push @commands, $infile_path;

    if ( $second_infile_path && $adapter_sequence_second ) {
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
            FILEHANDLE   => $FILEHANDLE,
            commands_ref => \@commands,
            separator    => $SPACE,

        }
    );
    return @commands;
}

1;
