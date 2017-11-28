package MIP::Program::Download::Wget_v5_10;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ wget };
}

sub wget {

## wget

## Function : Perl wrapper for writing wget recipe to a commands array. Based on GNU Wget 1.12, a non-interactive network retriever.
## Returns  : @commands
## Arguments: $outfile_path, $url, $stdoutfile_path, $stderrfile_path, stderrfile_path_append, $FILEHANDLE, $quiet, $verbose
##          : $outfile_path           => Outfile path. Write documents to FILE
##          : $url                    => Url to use for download
##          : $quiet                  => Quiet (no output)
##          : $verbose                => Verbosity

    my ($arg_href) = @_;

    ## Default(s)
    my $quiet;
    my $verbose;

    ## Flatten argument(s)
    my $outfile_path;
    my $url;

    ## Set active arguments
    $quiet        = $arg_href->{quiet};
    $verbose      = $arg_href->{verbose};
    $outfile_path = $arg_href->{outfile_path};
    $url          = $arg_href->{url};

    ## Wget
    # Stores commands depending on input parameters
    my @commands = qw{ wget };

    if ($quiet) {
        push @commands, q{--quiet};
    }

    if ($verbose) {
        push @commands, q{--verbose};
    }
    else {
        push @commands, q{--no-verbose};
    }

    ## URL
    push @commands, $url;

    ## Outfile
    if ($outfile_path) {
        push @commands, q{-O } . $outfile_path;
    }

    return @commands;
}

1;
