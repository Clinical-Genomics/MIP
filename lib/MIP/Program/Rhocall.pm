package MIP::Program::Rhocall;

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
    our @EXPORT_OK = qw{ rhocall_aggregate rhocall_annotate rhocall_viz };
}

sub rhocall_aggregate {

## Function : Perl wrapper for writing rhocall aggregate recipe to $FILEHANDLE or return commands array. Based on rhocall 0.3.
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $infile_path            => Infile path to read from
##          : $outfile_path           => Outfile path to write to
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path    => { store => \$outfile_path, strict_type => 1, },
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
    my @commands = qw{ rhocall aggregate };

    ## Options
    if ($outfile_path) {

        push @commands, q{--output} . $SPACE . $outfile_path;
    }

    push @commands, $infile_path;

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
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,

        }
    );
    return @commands;
}

sub rhocall_annotate {

## Function : Perl wrapper for writing rhocall annotate recipe to $FILEHANDLE or return commands array. Based on rhocall 0.3.
## Returns  : @commands
## Arguments: $bedfile_path           => BED file with AZ windows
##          : $FILEHANDLE             => Filehandle to write to
##          : $infile_path            => Infile path to read from
##          : $outfile_path           => Outfile path to write to
##          : $rohfile_path           => Rho style bcftools file
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $v14                    => Bcftools v1.4 or newer roh file including RG calls

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bedfile_path;
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $rohfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $v14;

    my $tmpl = {
        bedfile_path => { store => \$bedfile_path, strict_type => 1, },
        FILEHANDLE   => {
            store => \$FILEHANDLE,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path    => { store => \$outfile_path, strict_type => 1, },
        rohfile_path    => { store => \$rohfile_path, strict_type => 1, },
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
        v14 => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$v14,
            strict_type => 1,
        },

    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ rhocall annotate };

    ## Options
    if ($v14) {

        push @commands, q{--v14} . $SPACE;
    }

    if ($rohfile_path) {

        push @commands, q{-r} . $SPACE . $rohfile_path;
    }

    if ($bedfile_path) {

        push @commands, q{-b} . $SPACE . $bedfile_path;
    }

    if ($outfile_path) {

        push @commands, q{--output} . $SPACE . $outfile_path;
    }

    push @commands, $infile_path;

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
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,

        }
    );
    return @commands;
}

sub rhocall_viz {

## Function : Perl wrapper for writing rhocall viz recipe to $FILEHANDLE or return commands array. Based on rhocall 0.5.1.
## Returns  : @commands
## Arguments: $af_tag                 => Allele frequency to use
##          : $FILEHANDLE             => Filehandle to write to
##          : $infile_path            => Infile path to read from
##          : $outdir_path            => Directory path to write to
##          : $rohfile_path           => Rho style bcftools file
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $wig                    => Generate wig file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $af_tag;
    my $FILEHANDLE;
    my $infile_path;
    my $outdir_path;
    my $rohfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $wig;

    my $tmpl = {
        af_tag => {
            defined     => 1,
            store       => \$af_tag,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outdir_path => {
            defined     => 1,
            required    => 1,
            store       => \$outdir_path,
            strict_type => 1,
        },
        rohfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$rohfile_path,
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
        wig => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$wig,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ rhocall viz };

    if ($af_tag) {
        push @commands, q{--aftag} . $SPACE . $af_tag;
    }

    push @commands, q{--out_dir} . $SPACE . $outdir_path;

    push @commands, q{--roh} . $SPACE . $rohfile_path;

    if ($wig) {
        push @commands, q{--wig};
    }

    push @commands, $infile_path;

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
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,

        }
    );
    return @commands;
}
1;
