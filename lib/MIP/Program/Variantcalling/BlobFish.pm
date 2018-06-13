package MIP::Program::Variantcalling::BlobFish;

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
    our @EXPORT_OK = qw{ blobfish_all_vs_all };
}

## Constants
Readonly my $SPACE => q{ };

sub blobfish_all_vs_all {

## Function : Perl wrapper for BlobFish --allvsall, version 0.1.0.
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $infile_paths_ref       => list of input salmon quant.sf files
##          : $in_samples_ref         => list of input sample ids, these should match the infile path list
##          : $outdir_path            => Output Directory path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $tx_to_gene_path            => the path to the transcripts to gene conversion file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_paths_ref;
    my $in_samples_ref;
    my $outdir_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $tx_to_gene_path;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        infile_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infile_paths_ref,
            strict_type => 1,
        },
        in_samples_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$in_samples_ref,
            strict_type => 1,
        },
        outdir_path => {
            defined     => 1,
            required    => 1,
            store       => \$outdir_path,
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
        tx_to_gene_path => {
            defined     => 1,
            required    => 1,
            store       => \$tx_to_gene_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = q{BlobFish.py --allvsall};

    # Options
    push @commands, q{--dir} . $SPACE . $outdir_path;
    push @commands, q{--paths} . $SPACE . join $SPACE, @{$infile_paths_ref};
    push @commands, q{--sample} . $SPACE . join $SPACE, @{$in_samples_ref};
    push @commands, q{--tx} . $SPACE . $tx_to_gene_path;

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
