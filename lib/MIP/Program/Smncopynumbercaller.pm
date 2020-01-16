package MIP::Program::Smncopynumbercaller;

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
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ smn_caller };
}

## Constants
Readonly my $GENOME_VERSION_19 => 19;
Readonly my $GENOME_VERSION_37 => 37;
Readonly my $GENOME_VERSION_38 => 38;

sub smn_caller {

## Function : Perl wrapper for sma caller. Version git commit: 4b2c1ad
## Returns  : @commands
## Arguments: $filehandle             => Filehandle to write to
##          : $genome_version         => Genome version
##          : $manifest_file_path     => Manifest file path containing absolute path to BAM files
##          : $outfile_prefix         => Outfile prefix
##          : $outdir_path            => Outdir path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdinfile_path         => Stdinfile path
##          : $stdoutfile_path        => Stdoutfile path
##          : $thread_number          => Number of threads to use in calling

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $genome_version;
    my $manifest_file_path;
    my $outdir_path;
    my $outfile_prefix;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdinfile_path;
    my $stdoutfile_path;
    my $thread_number;

    ## Default(s)

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        genome_version => {
            allow    => [ $GENOME_VERSION_19, $GENOME_VERSION_37, $GENOME_VERSION_38, ],
            defined  => 1,
            required => 1,
            store    => \$genome_version,
            strict_type => 1,
        },
        manifest_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$manifest_file_path,
            strict_type => 1,
        },
        outdir_path => {
            defined     => 1,
            required    => 1,
            store       => \$outdir_path,
            strict_type => 1,
        },
        outfile_prefix => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_prefix,
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
        thread_number => {
            allow       => qr{ \A \d+ \z }sxm,
            default     => 1,
            store       => \$thread_number,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ smn_caller.py };

    push @commands, q{--manifest} . $SPACE . $manifest_file_path;

    push @commands, q{--genome} . $SPACE . $genome_version;

    push @commands, q{--prefix} . $SPACE . $outfile_prefix;

    push @commands, q{--outDir} . $SPACE . $outdir_path;

    if ($thread_number) {

        push @commands, q{--threads} . $SPACE . $thread_number;
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
