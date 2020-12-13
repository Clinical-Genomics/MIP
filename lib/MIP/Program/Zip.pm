package MIP::Program::Zip;

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

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ unzip };
}

sub unzip {

## Function : Perl wrapper for writing unzip recipe to $filehandle or return commands array. Based on unzip v6.0
## Returns  : @commands
##          : $filehandle             => Filehandle to write to (scalar undefined)
##          : $force                  => Overwrite existing files
##          : $infile_paths_ref       => Infile path
##          : $outdir_path            => Path to output directory
##          : $quiet                  => Suppress all warnings
##          : $stderrfile_path        => Stderrfile path (scalar )
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdout                 => Write on standard output
##          : $verbose                => Verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_paths_ref;
    my $outdir_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdout;

    ## Default(s)
    my $force;
    my $quiet;
    my $verbose;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        force => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$force,
            strict_type => 1,
        },
        infile_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infile_paths_ref,
            strict_type => 1,
        },
        outdir_path => {
            store       => \$outdir_path,
            strict_type => 1,
        },
        quiet => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$quiet,
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
        stdout => {
            store       => \$stdout,
            strict_type => 1,
        },
        verbose => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$verbose,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ unzip };

    if ($quiet) {

        push @commands, q{-q};
    }

    if ($verbose) {

        push @commands, q{-v};
    }

    if ($force) {

        push @commands, q{-o};
    }

    if ($stdout) {

        push @commands, q{-p};
    }

    push @commands, join $SPACE, @{$infile_paths_ref};

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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

1;
