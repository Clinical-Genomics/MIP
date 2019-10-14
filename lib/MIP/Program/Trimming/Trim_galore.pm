package MIP::Program::Trimming::Trim_galore;

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
    our @EXPORT_OK = qw{ trim_galore };
}

sub trim_galore {

## Function : Perl wrapper for Trim Galore. Based on version 0.6.4
## Returns  : @commands
## Arguments: $cores                  => Cores to be used by each process
##          : $fastqc                 => Run fastqc after trimming
##          : $FILEHANDLE             => Filehandle to write to
##          : $gzip_output            => Gzip output fastq file
##          : $infile_paths_ref       => Infile paths {REF}
##          : $outdir_path            => Outdirectory path
##          : $paired_reads           => Do paired end trimming
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdinfile_path         => Stdinfile path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $cores;
    my $FILEHANDLE;
    my $infile_paths_ref;
    my $outdir_path;
    my $paired_reads;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdinfile_path;
    my $stdoutfile_path;

    ## Default(s)
    my $fastqc;
    my $gzip_output;

    my $tmpl = {
        cores => {
            allow       => [ undef, qr/\A \d+ \z/xms ],
            store       => \$cores,
            strict_type => 1,
        },
        fastqc => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$fastqc,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        gzip_output => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$gzip_output,
            strict_type => 1,
        },
        infile_paths_ref => {
            default     => [],
            required    => 1,
            store       => \$infile_paths_ref,
            strict_type => 1,
        },
        outdir_path => {
            store       => \$outdir_path,
            strict_type => 1,
        },
        paired_reads => {
            allow       => [ undef, 0, 1 ],
            store       => \$paired_reads,
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
        stdinfile_path  => { store => \$stdinfile_path, strict_type => 1, },
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ trim_galore };

    if ($cores) {
        push @commands, q{--cores} . $SPACE . $cores;
    }

    if ($fastqc) {
        push @commands, q{--fastqc};
    }

    if ($gzip_output) {
        push @commands, q{--gzip};
    }

    if ($outdir_path) {
        push @commands, q{--output_dir} . $SPACE . $outdir_path;
    }

    if ($paired_reads) {
        push @commands, q{--paired};
    }

    push @commands, join $SPACE, @{$infile_paths_ref};

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdinfile_path         => $stdinfile_path,
            stdoutfile_path        => $stdoutfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,

        }
    );
    return @commands;
}

1;
