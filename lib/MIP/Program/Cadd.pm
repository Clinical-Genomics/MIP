package MIP::Program::Cadd;

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
use MIP::Constants qw { $SPACE };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ cadd cadd_install};
}

sub cadd {

## Function : Perl wrapper for dynamic annotation using CADD version 1.4.
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $genome_build           => Genome build
##          : $infile_path            => Infile path to read from
##          : $outfile_path           => Outfile path to write to
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdinfile_path         => Stdinfile path
##          : $stdoutfile_path        => Stdoutfile path
##          : $temp_dir_path          => Path to temp directory
##          : $version                => Version of CADD script

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $genome_build;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdinfile_path;
    my $stdoutfile_path;
    my $temp_dir_path;
    my $version;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        genome_build => {
            allow       => [qw{ GRCh37 GRCh38}],
            store       => \$genome_build,
            strict_type => 1,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
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
        temp_dir_path => {
            defined     => 1,
            required    => 1,
            store       => \$temp_dir_path,
            strict_type => 1,
        },
        version => {
            allow       => [qw{ v1.4 v1.5}],
            defined     => 1,
            required    => 1,
            store       => \$version,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ CADD.sh };

    ## Options
    if ($outfile_path) {

        push @commands, q{-o} . $SPACE . $outfile_path;
    }

    if ($version) {

        push @commands, q{-v} . $SPACE . $version;
    }
    if ($genome_build) {

        push @commands, q{-g} . $SPACE . $genome_build;
    }

    push @commands, q{-t} . $SPACE . $temp_dir_path;

    push @commands, $infile_path;

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

sub cadd_install {

## Function : Perl wrapper for installing/downloading CADD and its references. Based on CADD v1.5
## Returns  : @commands
## Arguments: $annotation_dir_path                 => Download path
##          : $download_annotations                => Download annotations
##          : $download_prescored_indel            => Download prescored indels
##          : $download_prescored_snv              => Download prescored snvs
##          : $download_prescored_with_annotations => Include annotations for prescored variants
##          : $FILEHANDLE                          => Filehandle to write to
##          : $stderrfile_path                     => Stderrfile path
##          : $stderrfile_path_append              => Append stderr info to file path
##          : $stdoutfile_path                     => Stdoutfile path
##          : $versions_ref                        => Version of CADD annotations {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $annotation_dir_path;
    my $download_annotations;
    my $download_prescored_indel;
    my $download_prescored_snv;
    my $download_prescored_with_annotations;
    my $FILEHANDLE;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $versions_ref;

    my $tmpl = {
        annotation_dir_path => {
            defined     => 1,
            required    => 1,
            store       => \$annotation_dir_path,
            strict_type => 1,
        },
        download_annotations => {
            allow       => [ undef, 0, 1 ],
            store       => \$download_annotations,
            strict_type => 1,
        },
        download_prescored_indel => {
            allow       => [ undef, 0, 1 ],
            store       => \$download_prescored_indel,
            strict_type => 1,
        },
        download_prescored_snv => {
            allow       => [ undef, 0, 1 ],
            store       => \$download_prescored_snv,
            strict_type => 1,
        },
        download_prescored_with_annotations => {
            allow       => [ undef, 0, 1 ],
            store       => \$download_prescored_with_annotations,
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
        versions_ref => {
            default     => [],
            required    => 1,
            store       => \$versions_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ bash install.sh };

    push @commands, q{-r} . $SPACE . $annotation_dir_path;

    if ($download_annotations) {

        push @commands, q{-a};
    }

    if ($download_prescored_indel) {

        push @commands, q{-i};
    }

    if ($download_prescored_snv) {

        push @commands, q{-n};
    }

    if ($download_prescored_with_annotations) {

        push @commands, q{-p};
    }

    push @commands, q{-v} . $SPACE . join q{ -v }, @{$versions_ref};

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
