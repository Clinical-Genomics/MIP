package MIP::Program::Variantcalling::Blobfish;

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
    our @EXPORT_OK = qw{ blobfish_allvsall };
}

## Constants
Readonly my $SPACE => q{ };

sub blobfish_allvsall {

## Function : Perl wrapper for BlobFish in allvsall mode (wrapper for DESeq2), based on version 0.0.2.
## Returns  : @commands
## Arguments: $conditions_ref         => Conditions for each indir path {REF}
##          : $filehandle             => Filehandle to write to
##          : $indir_paths_ref        => Path to salmon output directories {REF}
##          : $outdir_path            => Path to out directory
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $tx2gene_file_path      => Path to tx2gene file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conditions_ref;
    my $filehandle;
    my $indir_paths_ref;
    my $outdir_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $tx2gene_file_path;

    my $tmpl = {
        conditions_ref => {
            default     => [],
            required    => 1,
            store       => \$conditions_ref,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        indir_paths_ref => {
            default     => [],
            required    => 1,
            store       => \$indir_paths_ref,
            strict_type => 1,
        },
        outdir_path => {
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
        tx2gene_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$tx2gene_file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Each infile must have a condition
    croak q{Each indirectory must be matched with a condition}
      if ( scalar @{$indir_paths_ref} ne @{$conditions_ref} );

    ## Stores commands depending on input parameters
    my @commands = q{BlobFish.py --allvsall};

    ## Add infile paths
    push @commands, q{--paths} . $SPACE . join $SPACE, @{$indir_paths_ref};

    ## Add conditions
    push @commands, q{--conditions} . $SPACE . join $SPACE, @{$conditions_ref};

    ## Transcript to gene file
    push @commands, q{--tx} . $SPACE . $tx2gene_file_path;

    ## Outpath
    push @commands, q{--dir} . $SPACE . $outdir_path;

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
