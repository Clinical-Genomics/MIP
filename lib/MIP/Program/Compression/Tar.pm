package MIP::Program::Compression::Tar;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
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
    our @EXPORT_OK = qw{ tar };
}

## Constants
Readonly my $SPACE => q{ };

sub tar {

## Function : Perl wrapper for writing tar command recipe to $FILEHANDLE or return commands array. Based on tar 1.23.
## Returns  : "@commands"

##Arguments : $extract                => Extract files from an archive
##          : $filter_gzip            => Filter the archive through gzip
##          : $file                   => Use archive file or device ARCHIVE
##          : $outdir_path            => Extract to other than current directory
##          : $stdoutfile_path        => Stdoutfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $FILEHANDLE             => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $extract;
    my $filter_gzip;
    my $file;
    my $outdirectory_path;
    my $stdoutfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;

    my $tmpl = {
        extract => {
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$extract
        },
        filter_gzip => {
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$filter_gzip
        },
        file => {
            strict_type => 1,
            store       => \$file
        },
        outdirectory_path => {
            strict_type => 1,
            store       => \$outdirectory_path,
        },
        stdoutfile_path => {
            strict_type => 1,
            store       => \$stdoutfile_path
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append
        },
        FILEHANDLE => {
            store => \$FILEHANDLE
        },

    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ tar };

    if ($extract) {

        push @commands, q{--extract};
    }
    if ($filter_gzip) {

        push @commands, q{-z};
    }
    if ($file) {

        push @commands, q{--file=} . $file;
    }
    if ($outdirectory_path) {
        push @commands, q{--directory=} . $outdirectory_path;
    }

    push @commands,
      unix_standard_streams(
        {
            stdoutfile_path        => $stdoutfile_path,
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
