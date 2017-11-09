package MIP::Program::Compression::Zip;

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
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ unzip };
}

## Constants
Readonly my $SPACE => q{ };

sub unzip {

## Function : Perl wrapper for writing unzip recipe to $FILEHANDLE or return commands array. Based on unzip v6.0
## Returns  : @commands
##          : $infile_path            => Infile path
##          : $outdir_path            => Path to output directory
##          : $stdout                 => Write on standard output
##          : $FILEHANDLE             => Filehandle to write to (scalar undefined)
##          : $stderrfile_path        => Stderrfile path (scalar )
##          : $stderrfile_path_append => Append stderr info to file path
##          : $quiet                  => Suppress all warnings
##          : $verbose                => Verbosity
##          : $force                  => Overwrite existing files

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $outdir_path;
    my $stdout;
    my $FILEHANDLE;
    my $stderrfile_path;
    my $stderrfile_path_append;

    ## Default(s)
    my $quiet;
    my $verbose;
    my $force;

    my $tmpl = {
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path,
        },
        outdir_path => {
            strict_type => 1,
            store       => \$outdir_path,
        },
        stdout => {
            strict_type => 1,
            store       => \$stdout,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path,
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append,
        },
        quiet => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$quiet,
        },
        verbose => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$verbose,
        },
        force => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$force,
        },

    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = q{unzip};

    ## Options
    if ($quiet) {
        push @commands, q{-q};
    }

    if ($verbose) {
        push @commands, q{-v};
    }

    if ($force) {
        push @commands, q{-o};
    }

    #Write to stdout stream
    if ($stdout) {
        push @commands, q{-p};
    }

    ## Infile
    push @commands, $infile_path;

    ## Outfile
    if ($outdir_path) {
        push @commands, q{-d} . $SPACE . $outdir_path;
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
