package MIP::Program::Utility::Fusion_filter;

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
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ fusion_filter_gtf_file_to_feature_seqs fusion_filter_prep_genome_lib };
}

## Constants
Readonly my $SPACE => q{ };

sub fusion_filter_gtf_file_to_feature_seqs {

## Function : Perl wrapper for fusion filter gtf_file_to_feature_seqs command to $FILEHANDLE or return commands array. Based on Fusion filter v0.5.0
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $gtf_path               => Input gtf path
##          : $referencefile_path     => Reference sequence file
##          : $seq_type               => Sequence type
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $gtf_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $seq_type;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        gtf_path => {
            defined     => 1,
            required    => 1,
            store       => \$gtf_path,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
            required    => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        seq_type => {
            allow       => [qw{ cDNA CDS prot }],
            default     => q{cDNA},
            store       => \$seq_type,
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
    my @commands = qw{ gtf_file_to_feature_seqs.pl };

    # Transcripts file
    push @commands, $gtf_path;

    # Reference sequence file
    push @commands, $referencefile_path;

    # Sequence type
    push @commands, $seq_type;

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

sub fusion_filter_prep_genome_lib {

## Function : Perl wrapper for fusion filter prep_genome_lib command to $FILEHANDLE or return commands array. Based on Fusion filter v0.5.0
## Returns  : @commands
## Arguments: $blast_pairs_file_path  => Blast pair file path
##          : $FILEHANDLE             => Filehandle to write to
##          : $gtf_path               => Input gtf path
##          : $output_dir_path        => Output directory path
##          : $referencefile_path     => Reference sequence file
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $thread_number          => Number of threads (CPUs)

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $blast_pairs_file_path;
    my $FILEHANDLE;
    my $gtf_path;
    my $output_dir_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $thread_number;

    ## Default(s)

    my $tmpl = {
        blast_pairs_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$blast_pairs_file_path,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        gtf_path => {
            defined     => 1,
            required    => 1,
            store       => \$gtf_path,
            strict_type => 1,
        },
        output_dir_path => {
            defined     => 1,
            required    => 1,
            store       => \$output_dir_path,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
            required    => 1,
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
        thread_number => {
            allow       => qr/ ^\d+$ /sxm,
            default     => 1,
            store       => \$thread_number,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ prep_genome_lib.pl };

    # Transcripts file
    push @commands, q{--gtf} . $SPACE . $gtf_path;

    # Reference sequence file
    push @commands, q{--genome_fa} . $SPACE . $referencefile_path;

    # Sequence type
    push @commands, q{--blast_pairs} . $SPACE . $blast_pairs_file_path;

    push @commands, q{--cpu} . $SPACE . $thread_number;

    push @commands, q{--output_dir} . $SPACE . $output_dir_path;

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
