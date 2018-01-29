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
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ vardict vardict_var2vcf_paired vardict_var2vcf_single vardict_testsomatic vardict_teststrandbias };
}

## Constants
Readonly my $DQUOTE => q{"};
Readonly my $PIPE   => q{|};
Readonly my $SPACE  => q{ };

sub vardict {

## Function : Perl wrapper for Vardict, a variant caller. Based on vardict 2017.09.24
## Returns  : @commands
## Arguments: $af_threshold           => Threshold for allele frequency
##          : $FILEHANDLE             => Filehandle to write to
##          : $infile_bed_region_info => Infile path for region info bed file
##          : $infile_paths_ref       => Infile path array for normal and tumor bam file
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
    my $FILEHANDLE;
    my $infile_bed_region_info;
    my $infile_paths_ref;
    my $out_chrom_start;
    my $out_region_start;
    my $out_region_end;
    my $out_segment_annotn;
    my $referencefile_path;
    my $sample_name;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $af_threshold;

    my $tmpl = {
        af_threshold => {
            ## Exactly 2 decimal points after 0 or 1
            allow       => qr/ ^0.\d{1,2}$ | ^1$ /xsm,
            defined     => 1,
            default     => 0.01,
            store       => \$af_threshold,
            strict_type => 1,
        },
        infile_bed_region_info => {
            defined     => 1,
            required    => 1,
            store       => \$infile_bed_region_info,
            strict_type => 1,
        },
        infile_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infile_paths_ref,
            strict_type => 1,
        },
        out_chrom_start => {
            allow       => qr/ ^\d+$ /xsm,
            required    => 1,
            store       => \$out_chrom_start,
            strict_type => 1,
        },
        out_region_start => {
            allow       => qr/ ^\d+$ /xsm,
            required    => 1,
            store       => \$out_region_start,
            strict_type => 1,
        },
        out_region_end => {
            allow       => qr/ ^\d+$ /xsm,
            required    => 1,
            store       => \$out_region_end,
            strict_type => 1,
        },
        out_segment_annotn => {
            allow       => qr/ ^\d+$ /xsm,
            required    => 1,
            store       => \$out_segment_annotn,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
            required    => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        sample_name => {
            allow       => qr/ ^\w+$ /xsm,
            required    => 1,
            store       => \$sample_name,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
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
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,

        }
    );
    return @commands;
}

sub vardict_var2vcf_single {

## Function : Perl wrapper for var2vcf_singlelsls, converting variant output to vcf. Based on vardict 2017.09.24
## Returns  : @commands
## Arguments: $af_threshold           => Threshold for allele frequency
##          : $FILEHANDLE             => Filehandle to write to
##          : $sample_name            => Sample name to be used directly, will overwrite -n option
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $sample_name;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $af_threshold;

    my $tmpl = {
        af_threshold => {
            ## Exactly 2 decimal points after 0 or 1
            allow       => qr/ ^0.\d{1,2}$ | ^1$ /xsm,
            defined     => 1,
            default     => 0.01,
            required    => 1,
            store       => \$af_threshold,
            strict_type => 1,
        },
        sample_name => {
            allow       => qr/ ^\w+$ /xsm,
            required    => 1,
            store       => \$sample_name,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
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
    my @commands = q{var2vcf_valid.pl};

    push @commands, q{-f} . $SPACE . $af_threshold;

    push @commands, q{-N} . $SPACE . $sample_name;

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

sub vardict_var2vcf_paired {

    ## Function : Perl wrapper for var2vcf_paired, converting variant output to vcf. Based on vardict 2017.09.24
    ## Returns  : @commands
    ## Arguments: $af_threshold           => Threshold for allele frequency
    ##          : $FILEHANDLE             => Filehandle to write to
    ##          : $max_mm                 => The maximum mean mismatches allowed
    ##          : $max_pval               => The maximum p-valuem, set to 0 to keep all variants
    ##          : $sample_name            => Sample name to be used directly, will overwrite -n option
    ##          : $somatic_only           => Output only candidate somatic
    ##          : $stderrfile_path        => Stderrfile path
    ##          : $stderrfile_path_append => Append stderr info to file path
    ##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $sample_name;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $af_threshold;
    my $max_mm;
    my $max_pval;
    my $somatic_only;

    my $tmpl = {
        af_threshold => {
            ## Exactly 2 decimal points after 0 or 1
            allow       => qr/ ^0.\d{1,2}$ | ^1$ /xsm,
            defined     => 1,
            default     => 0.01,
            store       => \$af_threshold,
            strict_type => 1,
        },
        max_mm => {
            allow       => qr/ ^\d.\d*$ | ^\d*$ /xsm,
            defined     => 1,
            default     => 4.5,
            store       => \$max_mm,
            strict_type => 1,
        },
        max_pval => {
            ## Pvalue match must have a none zero digit as the last one
            allow       => qr/ ^0.\d*[1-9]$ /xsm,
            defined     => 1,
            default     => 0.9,
            store       => \$max_pval,
            strict_type => 1,
        },
        sample_name => {
            allow       => qr/ ^\w+$ /xsm,
            required    => 1,
            store       => \$sample_name,
            strict_type => 1,
        },
        somatic_only => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$somatic_only,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
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
    my @commands = q{var2vcf_paired.pl};

    if ($af_threshold) {

        # Allele frequency threshold
        push @commands, q{-f} . $SPACE . $af_threshold;
    }

    if ($max_mm) {

        # Maximum mean mismatches allowed
        push @commands, q{-m} . $SPACE . $max_mm;
    }

    if ($max_pval) {

        # Pvalue threshold for variants
        push @commands, q{-P} . $SPACE . $max_pval;
    }

    push @commands, q{-N} . $SPACE . $sample_name;

    if ($somatic_only) {

        # Output only somatic variants
        push @commands, q{-M} . $SPACE;
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

sub vardict_testsomatic {

## Function : Perl wrapper for testsomatic.R, A statistical test for AF of somatic variants. Based on vardict 2017.09.24
## Returns  : @commands
## Arguments:
##          : $FILEHANDLE             => Filehandle to write to
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
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
    my @commands = q{testsomatic.R};

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

sub vardict_teststrandbias {

## Function : Perl wrapper for testsomatic.R, A statistical test for AF for strand bias. Based on vardict 2017.09.24
## Returns  : @commands
## Arguments:
##          : $FILEHANDLE             => Filehandle to write to
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
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
    my @commands = q{teststrandbias.R};

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
