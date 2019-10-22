package MIP::Program::Qc::Rseqc;

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
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      rseqc_bam2wig
      rseqc_bam_stat
      rseqc_genebody_coverage2
      rseqc_infer_experiment
      rseqc_inner_distance
      rseqc_junction_annotation
      rseqc_junction_saturation
      rseqc_read_distribution
      rseqc_read_duplication
    };
}

## Constants
Readonly my $EQUAL           => q{=};
Readonly my $MIN_MAP_QUALITY => q{30};
Readonly my $SPACE           => q{ };

sub rseqc_bam2wig {

## Function : Perl wrapper for rseqc bam2wig.py. Version 3.0.0.
## Returns  : @commands
## Arguments: $chrom_size_file_path   => Tab separated file with chrom name and size
##          : $filehandle             => Filehandle to write to
##          : $infile_path            => Input file path
##          : $strand                 => Library strandedness
##          : $outfile_path_prefix    => Prefix for outfiles
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $chrom_size_file_path;
    my $filehandle;
    my $infile_path;
    my $strand;
    my $outfile_path_prefix;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        chrom_size_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$chrom_size_file_path,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        strand => {
            allow       => [ q{1++,1--,2+-,2-+}, q{1+-,1-+,2++,2--}, q{++,--}, q{+-,-+} ],
            store       => \$strand,
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

    ## Stores commands depending on input parameters
    my @commands = qw{ bam2wig.py };

    push @commands, q{--input-file} . $EQUAL . $infile_path;

    push @commands, q{--chromSize} . $EQUAL . $chrom_size_file_path;

    if ($strand) {
        push @commands, q{--strand} . $EQUAL . $strand;
    }

    push @commands, q{--out-prefix} . $EQUAL . $outfile_path_prefix;

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

sub rseqc_bam_stat {

## Function : Perl wrapper for rseqc bam_stat.py. Version 2.6.4.
## Returns  : @commands
## Arguments: $filehandle             => Filehandle to write to
##          : $infile_path            => Input file path
##          : $min_map_quality        => Minimum mapping quality
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $min_map_quality;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        min_map_quality => {
            allow       => qr/ ^\d+$ /sxm,
            default     => $MIN_MAP_QUALITY,
            store       => \$min_map_quality,
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
    my @commands = qw{ bam_stat.py };

    push @commands, q{--input-file} . $EQUAL . $infile_path;

    push @commands, q{--mapq} . $EQUAL . $min_map_quality;

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

sub rseqc_genebody_coverage2 {

## Function : Perl wrapper for rseqc geneBody_coverage.py. Version 3.0.0.
## Returns  : @commands
## Arguments: $bed_file_path          => Reference gene model
##          : $filehandle             => Filehandle to write to
##          : $infile_path            => Input BigWig file path
##          : $outfile_path_prefix    => Prefix for outfiles
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bed_file_path;
    my $filehandle;
    my $infile_path;
    my $outfile_path_prefix;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        bed_file_path => {
            defined     => 1,
            required    => 1,
            strict_type => 1,
            store       => \$bed_file_path,
        },
        filehandle => {
            store => \$filehandle,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path_prefix => {
            defined     => 1,
            required    => 1,
            strict_type => 1,
            store       => \$outfile_path_prefix,
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
    my @commands = qw{ geneBody_coverage2.py };

    push @commands, q{--input-file} . $EQUAL . $infile_path;

    push @commands, q{--refgene} . $EQUAL . $bed_file_path;

    push @commands, q{--out-prefix} . $EQUAL . $outfile_path_prefix;

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

sub rseqc_infer_experiment {

## Function : Perl wrapper for rseqc infer_experiment.py. Version 2.6.4.
## Returns  : @commands
## Arguments: $bed_file_path          => Transcript bed file path
##          : $infile_path            => Input file path
##          : $min_map_quality        => Minimum mapping quality
##          : $filehandle             => Filehandle to write to
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bed_file_path;
    my $filehandle;
    my $infile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $min_map_quality;

    my $tmpl = {
        bed_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$bed_file_path,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        min_map_quality => {
            allow       => qr/ ^\d+$ /sxm,
            default     => $MIN_MAP_QUALITY,
            strict_type => 1,
            store       => \$min_map_quality,
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
    my @commands = qw{ infer_experiment.py };

    push @commands, q{--input-file} . $EQUAL . $infile_path;

    push @commands, q{--refgene} . $EQUAL . $bed_file_path;

    push @commands, q{--mapq} . $EQUAL . $min_map_quality;

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

sub rseqc_inner_distance {

## Function : Perl wrapper for rseqc inner_distance.py. Version 2.6.4.
## Returns  : @commands
## Arguments: $bed_file_path          => Transcript bed file path
##          : $filehandle             => Filehandle to write to
##          : $infile_path            => Input file path
##          : $min_map_quality        => Minimum mapping quality
##          : $outfiles_path_prefix   => Outpath prefix for output files
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bed_file_path;
    my $filehandle;
    my $infile_path;
    my $outfiles_path_prefix;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $min_map_quality;

    my $tmpl = {
        bed_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$bed_file_path,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        min_map_quality => {
            allow       => qr/ ^\d+$ /sxm,
            default     => $MIN_MAP_QUALITY,
            strict_type => 1,
            store       => \$min_map_quality,
        },
        outfiles_path_prefix => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfiles_path_prefix,
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
    my @commands = qw{ inner_distance.py };

    push @commands, q{--input-file} . $EQUAL . $infile_path;

    push @commands, q{--refgene} . $EQUAL . $bed_file_path;

    push @commands, q{--mapq} . $EQUAL . $min_map_quality;

    push @commands, q{--out-prefix} . $EQUAL . $outfiles_path_prefix;

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

sub rseqc_junction_annotation {

## Function : Perl wrapper for rseqc junction_annotation.py. Version 2.6.4.
## Returns  : @commands
## Arguments: $bed_file_path          => Transcript bed file path
##          : $filehandle             => Filehandle to write to
##          : $infile_path            => Input file path
##          : $min_map_quality        => Minimum mapping quality
##          : $outfiles_path_prefix   => Outpath prefix for output files
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bed_file_path;
    my $filehandle;
    my $infile_path;
    my $outfiles_path_prefix;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $min_map_quality;

    my $tmpl = {
        bed_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$bed_file_path,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        min_map_quality => {
            allow       => qr/ ^\d+$ /sxm,
            default     => $MIN_MAP_QUALITY,
            strict_type => 1,
            store       => \$min_map_quality,
        },
        outfiles_path_prefix => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfiles_path_prefix,
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
    my @commands = qw{ junction_annotation.py };

    push @commands, q{--input-file} . $EQUAL . $infile_path;

    push @commands, q{--refgene} . $EQUAL . $bed_file_path;

    push @commands, q{--mapq} . $EQUAL . $min_map_quality;

    push @commands, q{--out-prefix} . $EQUAL . $outfiles_path_prefix;

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

sub rseqc_junction_saturation {

## Function : Perl wrapper for rseqc junction_saturation.py. Version 2.6.4.
## Returns  : @commands
## Arguments: $bed_file_path          => Transcript bed file path
##          : $filehandle             => Filehandle to write to
##          : $infile_path            => Input file path
##          : $min_map_quality        => Minimum mapping quality
##          : $outfiles_path_prefix   => Outpath prefix for output files
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bed_file_path;
    my $filehandle;
    my $infile_path;
    my $outfiles_path_prefix;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $min_map_quality;

    my $tmpl = {
        bed_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$bed_file_path,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        min_map_quality => {
            allow       => qr/ ^\d+$ /sxm,
            default     => $MIN_MAP_QUALITY,
            strict_type => 1,
            store       => \$min_map_quality,
        },
        outfiles_path_prefix => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfiles_path_prefix,
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
    my @commands = qw{ junction_saturation.py };

    push @commands, q{--input-file} . $EQUAL . $infile_path;

    push @commands, q{--refgene} . $EQUAL . $bed_file_path;

    push @commands, q{--mapq} . $EQUAL . $min_map_quality;

    push @commands, q{--out-prefix} . $EQUAL . $outfiles_path_prefix;

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

sub rseqc_read_distribution {

## Function : Perl wrapper for rseqc read_distribution.py. Version 2.6.4.
## Returns  : @commands
## Arguments: $bed_file_path          => Transcript bed file path
##          : $infile_path            => Input file path
##          : $filehandle             => Filehandle to write to
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bed_file_path;
    my $filehandle;
    my $infile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        bed_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$bed_file_path,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
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
    my @commands = qw{ read_distribution.py };

    push @commands, q{--input-file} . $EQUAL . $infile_path;

    push @commands, q{--refgene} . $EQUAL . $bed_file_path;

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

sub rseqc_read_duplication {

## Function : Perl wrapper for rseqc read_duplication.py. Version 2.6.4.
## Returns  : @commands
## Arguments: $filehandle             => Filehandle to write to
##          : $infile_path            => Input file path
##          : $min_map_quality        => Minimum mapping quality
##          : $outfiles_path_prefix   => Outpath prefix for output files
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $outfiles_path_prefix;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $min_map_quality;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        min_map_quality => {
            allow       => qr/ ^\d+$ /sxm,
            default     => $MIN_MAP_QUALITY,
            strict_type => 1,
            store       => \$min_map_quality,
        },
        outfiles_path_prefix => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfiles_path_prefix,
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
    my @commands = qw{ read_duplication.py };

    push @commands, q{--input-file} . $EQUAL . $infile_path;

    push @commands, q{--mapq} . $EQUAL . $min_map_quality;

    push @commands, q{--out-prefix} . $EQUAL . $outfiles_path_prefix;

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
