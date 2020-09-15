package MIP::Program::Telomerecat;

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

## MIPs lib/
use MIP::Constants qw{ $SPACE };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ telomerecat_bam2length };
}

sub telomerecat_bam2length {

## Function : Perl wrapper for telomerecat bam2length. Based on version 3.4.0
## Returns  : @commands
## Arguments: $filehandle             => Filehandle to write to
##          : $infile_paths_ref       => Bam files
##          : $outfile_path           => Outfile path
##          : $processes              => Number of processes
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdinfile_path         => Stdinfile path
##          : $stdoutfile_path        => Stdoutfile path
##          : $temp_directory         => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_paths_ref;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdinfile_path;
    my $stdoutfile_path;
    my $temp_directory;

    ## Default(s)
    my $processes;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        infile_paths_ref => {
            default     => [],
            required    => 1,
            store       => \$infile_paths_ref,
            strict_type => 1,
        },
        outfile_path => {
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        processes => {
            allow       => [ undef, qr/\A \d+ \z/xms ],
            default     => 1,
            store       => \$processes,
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
        temp_directory => { store => \$temp_directory, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ telomerecat bam2length };

    push @commands, q{--output} . $SPACE . $outfile_path;

    if ($processes) {

        push @commands, q{-p} . $SPACE . $processes;
    }

    if ($temp_directory) {

        push @commands, q{--temp_dir} . $SPACE . $temp_directory;
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
            filehandle   => $filehandle,
            separator    => $SPACE,

        }
    );
    return @commands;
}

1;
