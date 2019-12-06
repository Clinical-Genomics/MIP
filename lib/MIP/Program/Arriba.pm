package MIP::Program::Arriba;

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

    our $VERSION = 1.00;

    our @EXPORT_OK = qw{ arriba };
}

sub arriba {

## Function : Perl wrapper for Arriba commands module. Based on Arriba 1.1.0
## Returns  : @commands
## Arguments: $annotation_file_path       => Path to GTF file with annotations
##          : $blacklist_file_path        => Path to file with blacklist events
##          : $discarded_fusion_file_path => Path to write discarded fusion events to
##          : $filehandle                 => Filehandle to write to
##          : $genome_file_path           => Genome reference path
##          : $infile_path                => SAM/BAM file path
##          : $outfile_path               => Path to outfile
##          : $print_fusion_peptide       => Add fusion peptide sequence to outfile
##          : $print_fusion_transcript    => Add fusion transcript sequence to outfile
##          : $stderrfile_path            => Stderrfile path
##          : $stderrfile_path_append     => Append stderr info to file path
##          : $stdinfile_path             => Stdinfile path
##          : $stdoutfile_path            => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $annotation_file_path;
    my $blacklist_file_path;
    my $discarded_fusion_file_path;
    my $filehandle;
    my $genome_file_path;
    my $infile_path;
    my $outfile_path;
    my $print_fusion_peptide;
    my $print_fusion_transcript;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdinfile_path;
    my $stdoutfile_path;

    my $tmpl = {
        annotation_file_path => {
            required    => 1,
            store       => \$annotation_file_path,
            strict_type => 1,
        },
        blacklist_file_path => {
            store       => \$blacklist_file_path,
            strict_type => 1,
        },
        discarded_fusion_file_path => {
            store       => \$discarded_fusion_file_path,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        genome_file_path => {
            required    => 1,
            store       => \$genome_file_path,
            strict_type => 1,
        },
        infile_path => {
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path => {
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        print_fusion_peptide => {
            allow       => [ undef, 0, 1 ],
            store       => \$print_fusion_peptide,
            strict_type => 1,
        },
        print_fusion_transcript => {
            allow       => [ undef, 0, 1 ],
            store       => \$print_fusion_transcript,
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ arriba };

    push @commands, q{-g} . $SPACE . $annotation_file_path;

    if ($blacklist_file_path) {

        push @commands, q{-b} . $SPACE . $blacklist_file_path;
    }

    if ($discarded_fusion_file_path) {

        push @commands, q{-O} . $SPACE . $discarded_fusion_file_path;
    }

    push @commands, q{-a} . $SPACE . $genome_file_path;

    push @commands, q{-x} . $SPACE . $infile_path;

    push @commands, q{-o} . $SPACE . $outfile_path;

    if ($print_fusion_peptide) {

        push @commands, q{-P};
    }

    if ($print_fusion_transcript) {

        push @commands, q{-T};
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
            filehandle   => $filehandle,
            separator    => $SPACE,

        }
    );
    return @commands;
}

1;
