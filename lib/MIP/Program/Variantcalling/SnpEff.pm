package MIP::Program::Variantcalling::Snpeff;

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
    our @EXPORT_OK = qw{ snpeff_ann };
}

## Constants
Readonly my $DASH  => q{-};
Readonly my $SPACE => q{ };

sub snpeff_ann {

## Function : Perl wrapper for writing snpeff ann recipe to already open $FILEHANDLE or return commands array. Based on SnpEff 4.2 (build 2015-12-05).
## Returns  : @commands
## Arguments: $config_file_path       => Config file path
##          : $FILEHANDLE             => Filehandle to write to
##          : $genome_build_version   => Genome build version
##          : $infile_path            => Infile path
##          : $outfile_path           => Outfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $verbosity            => Increase output verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $config_file_path;
    my $FILEHANDLE;
    my $genome_build_version;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $verbosity;

    my $tmpl = {
        config_file_path => { strict_type => 1, store => \$config_file_path },
        FILEHANDLE       => {
            store => \$FILEHANDLE,
        },
        genome_build_version => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$genome_build_version
        },
        infile_path     => { strict_type => 1, store => \$infile_path },
        outfile_path    => { strict_type => 1, store => \$outfile_path },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path,
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append,
        },
        stdoutfile_path => {
            strict_type => 1,
            store       => \$stdoutfile_path,
        },
        verbosity => {
            allow       => qr/^\w+$/,
            strict_type => 1,
            store       => \$verbosity
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = q{ann};

    ## Options
    if ($verbosity) {

        push @commands, $DASH . $verbosity;
    }

    if ($genome_build_version) {

        push @commands, $genome_build_version;
    }

    if ($config_file_path) {

        push @commands, q{-config} . $SPACE . $config_file_path;
    }

    ## Infile
    if ($infile_path) {

        push @commands, $infile_path;
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
            FILEHANDLE   => $FILEHANDLE,
            commands_ref => \@commands,
            separator    => $SPACE,

        }
    );
    return @commands;
}

1;
