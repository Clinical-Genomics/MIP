package MIP::Program::Compression::Bzip2;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;    #Allow unicode characters in this script
use open qw{  :encoding(UTF-8) :std};
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };

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
    our @EXPORT_OK = qw{ bzip2 };
}

## Constants
Readonly my $SPACE => q{ };

sub bzip2 {

## Function : Perl wrapper for writing bzip2 recipe to $filehandle or return commands array. Based on bzip2 v1.0.6
## Returns  : "@commands"
## Arguments: $infile_path            => Infile path
##          : $outfile_path           => Path to output file
##          : $stdout                 => Write on standard output
##          : $decompress             => Decompress bzip2 file
##          : $filehandle             => Filehandle to write to
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $force                  => Overwrite of output files
##          : $quiet                  => Suppress all warnings
##          : $verbose                => Verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $outfile_path;
    my $decompress;
    my $filehandle;
    my $stderrfile_path;
    my $stderrfile_path_append;

    ## Default(s)
    my $stdout;
    my $force;
    my $quiet;
    my $verbose;

    my $tmpl = {
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        outfile_path => {
            strict_type => 1,
            store       => \$outfile_path
        },
        decompress => {
            strict_type => 1,
            store       => \$decompress,
        },
        filehandle => {
            store => \$filehandle
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append
        },
        stdout => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$stdout
        },
        force => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$force,
        },
        quiet => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$quiet
        },
        verbose => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$verbose
        },

    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = q{bzip2};

    ## Options
    if ($quiet) {
        push @commands, q{--quiet};
    }

    if ($verbose) {
        push @commands, q{--verbose};
    }

    if ($force) {
        push @commands, q{--force};
    }

    if ($decompress) {
        push @commands, q{--decompress};
    }

    ## Write to stdout stream
    if ($stdout) {
        push @commands, q{--stdout};
    }

    ## Infile
    push @commands, $infile_path;

    ## Outfile
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
            separator    => $SPACE,
            filehandle   => $filehandle,
        }
    );

    return @commands;
}

1;
