package MIP::Program::Pigz;

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
    use base qw{Exporter};

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ pigz };
}

sub pigz {

## Function : Perl wrapper for writing pigz recipe to $filehandle or return commands array. Based on pigz 2.3.1.
## Returns  : @commands
## Arguments: $decompress             => Decompress
##          : $filehandle             => Filehandle to write to
##          : $infile_path            => Infile path
##          : $outfile_path           => Outfile path
##          : $processes              => Allow up to n compression threads
##          : $quiet                  => Suppress all warnings
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdout                 => Write on standard output, keep original files unchanged
##          : $verbose                => Verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $decompress;
    my $filehandle;
    my $infile_path;
    my $outfile_path;
    my $processes;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdout;

    ## Default(s)
    my $quiet;
    my $verbose;

    my $tmpl = {
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        stdout       => { store => \$stdout,       strict_type => 1, },
        decompress   => { store => \$decompress,   strict_type => 1, },
        outfile_path => { store => \$outfile_path, strict_type => 1, },
        processes    => {
            allow       => qr{ \A\d+\z }sxm,
            store       => \$processes,
            strict_type => 1,
        },
        filehandle      => { store => \$filehandle, },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        quiet => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$quiet,
            strict_type => 1,
        },
        verbose => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$verbose,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ pigz };

    if ($processes) {

        push @commands, q{--processes} . $SPACE . $processes;
    }

    if ($quiet) {

        push @commands, q{--quiet};
    }

    if ($verbose) {

        push @commands, q{--verbose};
    }

    if ($decompress) {

        push @commands, q{--decompress};
    }

    if ($stdout) {

        push @commands, q{--stdout};
    }

    push @commands, $infile_path;

    if ($outfile_path) {

        push @commands, q{>} . $SPACE . $outfile_path;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
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
