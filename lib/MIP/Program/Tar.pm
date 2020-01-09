package MIP::Program::Tar;

use 5.026;
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
use MIP::Constants qw{ $SPACE };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.05;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ tar };
}

sub tar {

## Function: Perl wrapper for writing tar command recipe to $filehandle or return commands array. Based on tar 1.23.
## Returns : @commands
##Arguments: $create                 => Create tar archive
##         : $extract                => Extract files from an archive
##         : $file_path              => Use archive file or device ARCHIVE
##         : $filehandle             => Filehandle to write to
##         : $filter_gzip            => Filter the archive through gzip
##         : $in_paths_ref           => Files/directory to compress
##         : $outdir_path            => Extract to other than current directory
##         : $stderrfile_path        => Stderrfile path
##         : $stderrfile_path_append => Append stderr info to file path
##         : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $create;
    my $extract;
    my $file_path;
    my $filehandle;
    my $filter_gzip;
    my $in_paths_ref;
    my $outdirectory_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        create => {
            allow       => [ 0, 1 ],
            store       => \$create,
            strict_type => 1,
        },
        extract => {
            allow       => [ 0, 1 ],
            store       => \$extract,
            strict_type => 1,
        },
        file_path => {
            store       => \$file_path,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        filter_gzip => {
            allow       => [ 0, 1 ],
            store       => \$filter_gzip,
            strict_type => 1,
        },
        in_paths_ref => {
            default     => [],
            store       => \$in_paths_ref,
            strict_type => 1,
        },
        outdirectory_path => {
            store       => \$outdirectory_path,
            strict_type => 1,
        },
        stdoutfile_path => {
            store       => \$stdoutfile_path,
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ tar };

    if ($create) {

        push @commands, q{--create};
    }
    if ($extract) {

        push @commands, q{--extract};
    }
    if ($filter_gzip) {

        push @commands, q{-z};
    }
    if ($file_path) {

        push @commands, q{--file=} . $file_path;
    }
    if ($outdirectory_path) {

        push @commands, q{--directory=} . $outdirectory_path;
    }
    if ($in_paths_ref) {

        push @commands, @{$in_paths_ref};
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
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );

    return @commands;
}

1;
