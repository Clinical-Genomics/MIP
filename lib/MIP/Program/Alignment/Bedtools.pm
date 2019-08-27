package MIP::Program::Alignment::Bedtools;

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
use MIP::Constants qw{ $SPACE };
use MIP::Unix::Standard_streams qw{unix_standard_streams};
use MIP::Unix::Write_to_file qw{unix_write_to_file};

BEGIN {
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.03;

    # Inherit from Exporter to export functions and variables
    use base qw { Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ bedtools_genomecov bedtools_makewindows };

}

sub bedtools_genomecov {

## Function : Perl wrapper for writing bedtools genomecov recipe to $FILEHANDLE. Based on bedtools 2.26.0.
## Returns  : "@commands"
## Arguments: $bam_infile_path        => BAM infile path
##          : $depth_each_position    => Report the depth at each genome position (with zero-based coordinates)
##          : $FILEHANDLE             => Sbatch filehandle to write to
##          : $infile_path            => Infile paths
##          : $max_coverage           => Combine all positions with a depth >= max into a single bin in the histogram
##          : $referencefile_path     => Genome reference file
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Stderrfile path append
##          : $stdoutfile_path        => Outfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bam_infile_path;
    my $FILEHANDLE;
    my $infile_path;
    my $max_coverage;
    my $referencefile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $depth_each_position;

    my $tmpl = {
        bam_infile_path => {
            store       => \$bam_infile_path,
            strict_type => 1,
        },
        depth_each_position => {
            allow       => qr/\A \d+ \z/sxm,
            default     => 0,
            store       => \$depth_each_position,
            strict_type => 1,
        },
        FILEHANDLE  => { store => \$FILEHANDLE },
        infile_path => {
            strict_type => 1,
            store       => \$infile_path
        },
        max_coverage => {
            allow       => qr/\A \d+ \z/sxm,
            strict_type => 1,
            store       => \$max_coverage
        },
        referencefile_path => {
            strict_type => 1,
            store       => \$referencefile_path
        },
        stderrfile_path        => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append => { strict_type => 1, store => \$stderrfile_path_append },
        stdoutfile_path        => { strict_type => 1, store => \$stdoutfile_path },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    ## Array @commands stores commands depending on input parameters
    my @commands = qw{ bedtools genomecov };

    ## Options
    if ($bam_infile_path) {

        push @commands, q{-ibam} . $SPACE . $bam_infile_path;
    }
    if ($infile_path) {

        push @commands, q{-i} . $SPACE . $infile_path;
    }

    if ($max_coverage) {

        push @commands, q{-max} . $SPACE . $max_coverage;
    }
    if ($depth_each_position) {

        push @commands, q{-dz};
    }
    if ($referencefile_path) {

        push @commands, q{-g} . $SPACE . $referencefile_path;
    }

    # Redirect stderr output to program specific stderr file
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
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub bedtools_makewindows {

    ## Function : Perl wrapper for writing bedtools makewindows for bed files recipe to $FILEHANDLE. Based on bedtools 2.26.0.
    ## Returns  : "@commands"
    ## Arguments: $FILEHANDLE             => Sbatch filehandle to write to
    ##          : infile_bed_path         => Input BED file (with chrom,start,end fields).
    ##          : $stderrfile_path        => Stderrfile path
    ##          : $stderrfile_path_append => Stderrfile path append
    ##          : $stdoutfile_path        => Outfile path
    ##          : $step_size              => Step size (bp): i.e., how many base pairs to step before creating a new window.
    ##          : $window_size            => Divide each input interval (bp)

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_bed_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $step_size;
    my $window_size;

    my $tmpl = {
        FILEHANDLE      => { store => \$FILEHANDLE },
        infile_bed_path => {
            required    => 1,
            store       => \$infile_bed_path,
            strict_type => 1,
        },
        step_size => {
            allow       => qr/^\d+$/,
            store       => \$step_size,
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
        window_size => {
            allow       => qr/^\d+$/,
            required    => 1,
            store       => \$window_size,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    ## Array @commands stores commands depending on input parameters
    my @commands = qw{bedtools makewindow};

    ## Infile
    push @commands, q{-b} . $SPACE . $infile_bed_path;

    if ($step_size) {
        ## Step size
        push @commands, q{-s} . $SPACE . $step_size;
    }

    ## Window size
    push @commands, q{-w} . $SPACE . $window_size;

    # Redirect stderr output to program specific stderr file
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
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );
    return @commands;
}

1;
