package MIP::Program::Chromograph;

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
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ chromograph };
}

sub chromograph {

## Function : Perl wrapper for chromograph.
## Returns  : @commands
## Arguments: $coverage_file_path     => Coverage data infile
##          : $FILEHANDLE             => Filehandle to write to
##          : $outdir_path            => Outdir path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdinfile_path         => Stdinfile path
##          : $stdoutfile_path        => Stdoutfile path
##          : $step                   => Bin size for wig file
##          : upd_regions_file_path   => UPD regions data file
##          : upd_sites_file_path     => UPD sites data file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $coverage_file_path;
    my $FILEHANDLE;
    my $outdir_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdinfile_path;
    my $stdoutfile_path;
    my $step;
    my $upd_regions_file_path;
    my $upd_sites_file_path;

    my $tmpl = {
        coverage_file_path => {
            store       => \$coverage_file_path,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        outdir_path => {
            defined     => 1,
            required    => 1,
            store       => \$outdir_path,
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
        step => {
            allow       => qr/\A \d+ \z/xms,
            store       => \$step,
            strict_type => 1,
        },
        upd_regions_file_path => {
            store       => \$upd_regions_file_path,
            strict_type => 1,
        },
        upd_sites_file_path => {
            store       => \$upd_sites_file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ chromograph };

    if ($coverage_file_path) {

        push @commands, q{--coverage} . $SPACE . $coverage_file_path;
    }

    push @commands, q{--outd} . $SPACE . $outdir_path;

    if ($step) {

        push @commands, q{--step} . $SPACE . $step;
    }

    if ($upd_regions_file_path) {

        push @commands, q{--regions} . $SPACE . $upd_regions_file_path;
    }

    if ($upd_sites_file_path) {

        push @commands, q{--upd} . $SPACE . $upd_sites_file_path;
    }

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
