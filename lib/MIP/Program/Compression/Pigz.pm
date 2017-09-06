package MIP::Program::Compression::Pigz;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;    #Allow unicode characters in this script
use open qw{ :encoding(UTF-8) :std };
use charnames qw( :full :short );
use Carp;
use English qw{-no_match_vars};
use Params::Check qw{check allow last_error};

use Readonly;

use FindBin qw{$Bin};    #Find directory of script
use File::Basename qw{dirname};
use File::Spec::Functions qw{catdir};

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Unix::Standard_streams qw{unix_standard_streams};
use MIP::Unix::Write_to_file qw{unix_write_to_file};

BEGIN {
    require Exporter;
    use base qw{Exporter};

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{pigz};
}

## Constants
Readonly my $SPACE => q{ };

sub pigz {

## pigz

## Function : Perl wrapper for writing pigz recipe to $FILEHANDLE or return commands array. Based on pigz 2.3.1.
## Returns  : "@commands"
## Arguments: $FILEHANDLE, $stdout, $stderrfile_path, $stderrfile_path_append, $decompress, $infile_path, $outfile_path, $processes, $quiet, $verbose
##          : $FILEHANDLE             => Filehandle to write to
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdout                 => Write on standard output, keep original files unchanged
##          : $decompress             => Decompress
##          : $infile_path            => Infile path
##          : $outfile_path           => Outfile path
##          : $processes              => Allow up to n compression threads
##          : $quiet                  => Suppress all warnings
##          : $verbose                => Verbosity

    my ($arg_href) = @_;

    ## Default(s)
    my $quiet;
    my $verbose;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdout;
    my $decompress;
    my $infile_path;
    my $outfile_path;
    my $processes;

    my $tmpl = {
        FILEHANDLE      => { store       => \$FILEHANDLE },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        stdout      => { strict_type => 1, store => \$stdout },
        decompress  => { strict_type => 1, store => \$decompress },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        outfile_path => { strict_type => 1, store => \$outfile_path },
        processes    => {
            allow       => qr/^\d+$/,
            strict_type => 1,
            store       => \$processes
        },
        quiet => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$quiet
        },
        verbose => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$verbose
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## pigz
    my @commands = qw{pigz};    #Stores commands depending on input parameters

    ## Options
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

    if ($stdout) {    #Write to stdout stream

        push @commands, q{--stdout};
    }

    ## FILE
    push @commands, $infile_path;

    ## Outfile
    if ($outfile_path) {

        #Outfile
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;
}

1;
