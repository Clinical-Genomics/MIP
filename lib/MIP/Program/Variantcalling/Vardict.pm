package MIP::Program::Variantcalling::Vardict;

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
    our @EXPORT_OK = qw{ vardict };
}

## Constants
Readonly my $SPACE => q{ };

sub vardict {

    ## Function : Perl wrapper for Vardict, a variant caller.
    ## Returns  : @commands
    ## Arguments: $FILEHANDLE             => Filehandle to write to
    ##          : $referencefile_path     => Genome reference file
    ##          : $af_threshold           => The threshold for allele frequency
    ##          : $sample_name            => The sample name to be used directly, will overwrite -n option
    ##          : $infile_path_normal     => Infile path normal
    ##          : $infile_path_tumor      => Infile path tumor
    ##          : $out_chrom_start        => The column for chromosome
    ##          : $out_region_start       => The column for region start
    ##          : $out_region_end         => The column for region end
    ##          : $out_segment_annot      => The column for gene name, or segment annotation
    ##          : $infile_bed_region_info => Infile path for region info bed file
    ##          : $stderrfile_path        => Stderrfile path
    ##          : $stderrfile_path_append => Append stderr info to file path
    ##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)

    my $tmpl = {
        referencefile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$referencefile_path
        },
        af_threshold => {
            required => 1,
            defined  => 1,
            ## FIXME: fix the regex match for float point
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            strore      => \$af_thresold
        },
        sample_name => {
            required    => 1,
            allow       => qr/ ^\w+$ /xsm,
            strict_type => 1,
            strore      => \$sample_name
        },
        infile_path => {
            required    => 1,
            strict_type => 1,
            defined     => 1,
            store       => \$infile_path
        },
        out_chrom_start => {
            required    => 1,
            strict_type => 1,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$out_chrom_start
        },
        out_region_start => {
            required    => 1,
            strict_type => 1,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$out_region_start
        },
        out_region_end => {
            required    => 1,
            strict_type => 1,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$out_region_end
        },
        out_segment_annot => {
            required    => 1,
            strict_type => 1,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$out_segment_annot
        },
        infile_bed_region_info => {
            required    => 1,
            strict_type => 1,
            defined     => 1,
            store       => \$infile_bed_region_info
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path,
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append,
        },
        stdoutfile_path => {
            strict_type => 1,
            store       => \$stdoutfile_path,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = q{vardict};

    if ($referencefile_path) {
        push @commands, q{-G} . $SPACE . $referencefile_path;
    }
    if ($af_threshold) {
        push @commands, q{-f} . $SPACE . $af_threshold;
    }
    if ($sample_name) {
        push @commands, q{-N} . $SPACE . $sample_name;
    }
    ## FIXME: add two infile_path to concatenate as "tumor|normal" (with double quotation) for vardict input
    if ($infile_path) {
        push @commands, q{-b} . $SPACE . $infile_path;
    }
    if ($out_chrom_start) {
        push @commands, q{-c} . $SPACE . $out_chrom_start;
    }
    if ($out_region_start) {
        push @commands, q{-S} . $SPACE . $out_region_start;
    }
    if ($out_region_end) {
        push @commands, q{-E} . $SPACE . $out_region_end;
    }
    if ($out_segment_annot) {
        push @commands, q{-g} . $SPACE . $out_segment_annot;
    }
    if ($infile_bed_region_info) {
        push @commands, $SPACE . $infile_bed_region_info;
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
