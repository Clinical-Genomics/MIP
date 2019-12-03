package MIP::Program::Gffcompare;

use 5.026;
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
use MIP::Constants qw{ $SPACE };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ gffcompare };
}

sub gffcompare {

## Function : Perl wrapper for GffCompare. Based on version 0.10.1.
## Returns  : @commands
## Arguments: $filehandle                 => Filehandle to write to
##          : $genome_sequence_path       => Genome sequence
##          : $gtf_reference_path         => Input GTF refrence file
##          : $ignore_non_overlapping_ref => Ignore reference transcripts that are not overlapped by any input gtf
##          : $infile_paths_ref           => Input bam file paths {REF}
##          : $outfile_path_prefix        => Path and prefix to prepend to outfiles
##          : $stderrfile_path            => Stderrfile path
##          : $stderrfile_path_append     => Append stderr info to file path
##          : $stdoutfile_path            => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $genome_sequence_path;
    my $gtf_reference_path;
    my $ignore_non_overlapping_ref;
    my $infile_paths_ref;
    my $outfile_path_prefix;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        genome_sequence_path => {
            store       => \$genome_sequence_path,
            strict_type => 1,
        },
        gtf_reference_path => {
            defined     => 1,
            required    => 1,
            store       => \$gtf_reference_path,
            strict_type => 1,
        },
        ignore_non_overlapping_ref => {
            allow       => [ undef, 0, 1 ],
            store       => \$ignore_non_overlapping_ref,
            strict_type => 1,
        },
        infile_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infile_paths_ref,
            strict_type => 1,
        },
        outfile_path_prefix => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path_prefix,
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

    my @commands = q{gffcompare};

    if ($genome_sequence_path) {

        push @commands, q{-s} . $SPACE . $genome_sequence_path;
    }

    push @commands, q{-r} . $SPACE . $gtf_reference_path;

    if ($ignore_non_overlapping_ref) {

        push @commands, q{-R};
    }

    push @commands, q{-o} . $SPACE . $outfile_path_prefix;

    push @commands, join $SPACE, @{$infile_paths_ref};

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
