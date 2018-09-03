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
    our @EXPORT_OK = qw{ fusion_filter_gtf_file_to_feature_seqs };
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
    my $seq_type;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)

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
            FILEHANDLE   => $FILEHANDLE,
            commands_ref => \@commands,
            separator    => $SPACE,

        }
    );
    return @commands;
}

1;
