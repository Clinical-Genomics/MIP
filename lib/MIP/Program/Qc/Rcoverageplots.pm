package MIP::Program::Qc::Rcoverageplots;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
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
    our @EXPORT_OK = qw{ covplots_genome };
}

## Constants
Readonly my $SPACE => q{ };

sub covplots_genome {

## covplots_genome

## Function : Perl wrapper for R script to generate coverage plots from bedtools
##            genomecov.
## Returns  : @commands

## Arguments: $infile_path, $outdirectory_path, $sample_id, $stdoutfile_path, $stderrfile_path, $stderrfile_path_append, $filehandle, $max_coverage_depth
##          : $infile_path            => Infile path
##          : $outdirectory_path      => Outfile path
##          : $sample_id              => Sample id
##          : $stdoutfile_path        => Stdoutfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $filehandle             => Filehandle to write to
##          : $max_coverage_depth     => Max coverage depth

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $outdirectory_path;
    my $sample_id;
    my $stdoutfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $filehandle;

    ## Default(s)
    my $max_coverage_depth;

    my $tmpl = {
        infile_path => {
            required    => 1,
            default     => 1,
            strict_type => 1,
            store       => \$infile_path,
        },
        outdirectory_path => {
            required    => 1,
            default     => 1,
            strict_type => 1,
            store       => \$outdirectory_path,
        },
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id
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
        filehandle => {
            store => \$filehandle,
        },
        max_coverage_depth => {
            default     => 30,
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$max_coverage_depth,
        },

    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = q{covplots_genome};

    ## Infile path
    push @commands, $infile_path;

    ## Sample id
    push @commands, $sample_id;

    ## X-axis max scale
    if ($max_coverage_depth) {

        push @commands, $max_coverage_depth;
    }

    ## Outfile path
    push @commands, $outdirectory_path;

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
            filehandle   => $filehandle,
        }
    );

    return @commands;
}

1;
