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
Readonly my $SPACE  => q{ };
Readonly my $PIPE   => q{|};
Readonly my $DQUOTE => q{"};

sub vardict {

## Function : Perl wrapper for Vardict, a variant caller. Based on vardict 2017.09.24
## Returns  : @commands
## Arguments: $af_threshold           => Threshold for allele frequency
##          : $infile_bed_region_info => Infile path for region info bed file
##          : $infile_path_normal     => Infile path normal
##          : $infile_path_tumor      => Infile path tumor
##          : $FILEHANDLE             => Filehandle to write to
##          : $out_chrom_start        => Column for chromosome
##          : $out_region_start       => Column for region start
##          : $out_region_end         => Column for region end
##          : $out_segment_annotn     => Column for gene name, or segment annotation
##          : $referencefile_path     => Genome reference file
##          : $sample_name            => Sample name to be used directly, will overwrite -n option
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_bed_region_info;
    my $infile_paths_ref;
    my $out_chrom_start;
    my $out_region_start;
    my $out_region_end;
    my $out_segment_annotn;
    my $referencefile_path;
    my $sample_name;
    my $FILEHANDLE;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $af_threshold;

    my $tmpl = {
        af_threshold => {
            required => 1,
            defined  => 1,
            default  => 0.01,
            ## Exactly 2 decimal points after 0 or 1
            allow       => qr/ ^0.\d{1,2}$ | ^1$ /xsm,
            strict_type => 1,
            store       => \$af_threshold,
        },
        infile_bed_region_info => {
            required    => 1,
            strict_type => 1,
            defined     => 1,
            store       => \$infile_bed_region_info,
        },
        infile_paths_ref => {
            required    => 1,
            strict_type => 1,
            default     => [],
            defined     => 1,
            store       => \$infile_paths_ref,
        },
        out_chrom_start => {
            required    => 1,
            strict_type => 1,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$out_chrom_start,
        },
        out_region_start => {
            required    => 1,
            strict_type => 1,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$out_region_start,
        },
        out_region_end => {
            required    => 1,
            strict_type => 1,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$out_region_end,
        },
        out_segment_annotn => {
            required    => 1,
            strict_type => 1,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$out_segment_annotn,
        },
        referencefile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$referencefile_path,
        },
        sample_name => {
            required    => 1,
            allow       => qr/ ^\w+$ /xsm,
            strict_type => 1,
            store       => \$sample_name,
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

    push @commands, q{-f} . $SPACE . $af_threshold;

    ## Vardict requires input in form of "tumor_file_name|normal_file_name" or "tumor_file_name"
    push @commands,
      q{-b} . $SPACE . $DQUOTE . join( $PIPE, @{$infile_paths_ref} ) . $DQUOTE;

    push @commands, q{-c} . $SPACE . $out_chrom_start;

    push @commands, q{-S} . $SPACE . $out_region_start;

    push @commands, q{-E} . $SPACE . $out_region_end;

    push @commands, q{-g} . $SPACE . $out_segment_annotn;

    push @commands, q{-G} . $SPACE . $referencefile_path;

    push @commands, q{-N} . $SPACE . $sample_name;

    push @commands, $infile_bed_region_info;

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
