package MIP::Program::Rsync;

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
use MIP::Constants qw{ $EQUALS $SPACE };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ rsync };
}

sub rsync {

## Function : Perl wrapper for rsync 3.1.2  protocol version 31.
## Returns  : @commands
## Arguments: $archive                => Archive mode
##          : $compress               => Compress file data during the transfer
##          : $copy_links             => Transform symlink into referent file/dir
##          : $destination            => Destination to rsync to
##          : $FILEHANDLE             => Filehandle to write to
##          : $source                 => Source to rsync from
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdinfile_path         => Stdinfile path
##          : $stdoutfile_path        => Stdoutfile path
##          : $temporary_dir          => Temporary dir

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $destination;
    my $FILEHANDLE;
    my $source;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdinfile_path;
    my $stdoutfile_path;
    my $temporary_dir;

    ## Default(s)
    my $archive;
    my $compress;
    my $copy_links;

    my $tmpl = {
        archive => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$archive,
            strict_type => 1,
        },
        compress => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$compress,
            strict_type => 1,
        },
        copy_links => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$copy_links,
            strict_type => 1,
        },
        destination => {
            defined     => 1,
            required    => 1,
            store       => \$destination,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        source => {
            defined     => 1,
            required    => 1,
            store       => \$source,
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
        temporary_dir => {
            store       => \$temporary_dir,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ rsync };

    if ($archive) {

        push @commands, q{--archive};
    }

    if ($compress) {

        push @commands, q{--compress};
    }

    if ($copy_links) {
        push @commands, q{--copy-links};
    }

    if ($temporary_dir) {

        push @commands, q{--temp-dir} . $EQUALS . $temporary_dir;
    }

    push @commands, $source;

    push @commands, $destination;

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
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,

        }
    );
    return @commands;
}

1;
