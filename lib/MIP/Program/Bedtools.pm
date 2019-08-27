package MIP::Program::Bedtools;

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
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ bedtools_genomecov bedtools_intersectbed bedtools_makewindows };
}

sub bedtools_genomecov {

## Function : Perl wrapper for writing Bedtools genomecov recipe to already open $FILEHANDLE or return commands array. Based on Bedtools 2.26.0.
## Returns  : @commands
## Arguments: $bam_infile_path        => BAM infile path
##          : $depth_each_position    => Report the depth at each genome position (with zero-based coordinates)
##          : $FILEHANDLE             => Filehandle to write to
##          : $infile_path            => Infile path (bed/gff/vcf)
##          : $max_coverage           => Combine all positions with a depth >= max into a single bin in the histogram
##          : $referencefile_path     => Genome reference file
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

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
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        infile_path  => { store => \$infile_path, strict_type => 1, },
        max_coverage => {
            allow       => qr/\A \d+ \z/sxm,
            store       => \$max_coverage,
            strict_type => 1,
        },
        referencefile_path => {
            store       => \$referencefile_path,
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

    ## Stores commands depending on input parameters
    my @commands = qw{ bedtools genomecov };

    ## Options
    if ($bam_infile_path) {

        push @commands, q{-ibam} . $SPACE . $bam_infile_path;
    }

    ## Infile
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

sub bedtools_intersectbed {

## Function : Perl wrapper for writing Bedtools intersectbed recipe to already open $FILEHANDLE or return commands array. Based on Bedtools 2.26.0.
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $infile_path            => Infile path
##          : $intersectfile_path     => Intersect file (-b)
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $with_header            => Include header from infile in output

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $intersectfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $with_header;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        infile_path        => { store => \$infile_path, strict_type => 1, },
        intersectfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$intersectfile_path,
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
        with_header => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$with_header,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ intersectBed };

    ## Options
    if ($with_header) {

        push @commands, q{-header};
    }

    ## Infile
    if ($infile_path) {

        push @commands, q{-a} . $SPACE . $infile_path;
    }

    if ($intersectfile_path) {

        push @commands, q{-b} . $SPACE . $intersectfile_path;
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
##          : $infile_bed_path        => Input BED file (with chrom,start,end fields).
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
        FILEHANDLE      => { store => \$FILEHANDLE, },
        infile_bed_path => {
            required    => 1,
            store       => \$infile_bed_path,
            strict_type => 1,
        },
        step_size => {
            allow       => qr/\A \d+ \z/sxm,
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
            allow       => qr/\A \d+ \z/sxm,
            required    => 1,
            store       => \$window_size,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ bedtools makewindows };

    push @commands, q{-b} . $SPACE . $infile_bed_path;

    if ($step_size) {

        push @commands, q{-s} . $SPACE . $step_size;
    }

    push @commands, q{-w} . $SPACE . $window_size;

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
