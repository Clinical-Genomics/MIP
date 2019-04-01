package MIP::IO::Files;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use FindBin qw{ $Bin };
use File::Basename qw{ basename dirname fileparse };
use File::Spec::Functions qw{ catdir catfile splitpath };
use strict;
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie;
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $AMPERSAND $ASTERISK $DOT $EMPTY_STR $NEWLINE $SPACE $UNDERSCORE };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.05;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ migrate_file };

}

sub migrate_file {

## Function : Copy file to from source ($infile_path) to destination ($outfile_path).
## Returns  : $infile_path_file_name
## Arguments: $FILEHANDLE      => Filehandle to write to
##          : $infile_path     => Infile path
##          : $outfile_path    => Outfile path
##          : $recursive       => Copy directories recursively
##          : $stderrfile_path => Stderrfile path
##          : $xargs           => Use xargs if defined {Optional}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $xargs;

    ## Default(s)
    my $recursive;

    my $tmpl = {
        FILEHANDLE  => { defined => 1, required => 1, store => \$FILEHANDLE, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        recursive => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$recursive,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        xargs           => { store => \$xargs,           strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Gnu::Coreutils qw{ gnu_cp };

    ## Split relative infile_path to file(s)
    my ( $infile_path_volume, $infile_path_directory, $infile_path_file_name ) =
      splitpath($infile_path);

    gnu_cp(
        {
            FILEHANDLE      => $FILEHANDLE,
            preserve        => 1,
            recursive       => $recursive,
            infile_path     => $infile_path,
            outfile_path    => $outfile_path,
            stderrfile_path => $stderrfile_path,
        }
    );

    # For print wait statement downstream
    if ( not defined $xargs ) {

        say {$FILEHANDLE} $AMPERSAND . $SPACE;
    }
    print {$FILEHANDLE} $NEWLINE;

    return $infile_path_file_name;
}

1;
