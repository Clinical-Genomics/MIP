package MIP::File::Decompression;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Basename qw{ fileparse };
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
use MIP::Constants qw{ $DOT $NEWLINE $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ decompress_file };
}

sub decompress_file {

## Function : Check if file needs to be decompress and write decompression if so
## Returns  :
## Arguments: $decompress_program => Decompress the downloaded file using program supplied
##          : $FILEHANDLE         => Filehandle to write to
##          : $outdir_path        => Outdirectory path
##          : $outfile_path       => Outfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $decompress_program;
    my $FILEHANDLE;
    my $outdir_path;
    my $outfile_path;

    my $tmpl = {
        decompress_program => { store => \$decompress_program, strict_type => 1, },
        FILEHANDLE  => { defined => 1, required => 1, store => \$FILEHANDLE, },
        outdir_path => {
            default     => cwd(),
            store       => \$outdir_path,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Parse::File qw{ parse_file_suffix };
    use MIP::Program::Compression::Gzip qw{ gzip };
    use MIP::Program::Compression::Tar qw{ tar };
    use MIP::Program::Compression::Zip qw{ unzip };

    return if ( not defined $decompress_program );

    if ( $decompress_program eq q{gzip} ) {

        ## Parse file suffix in filename.suffix(.gz).
        ## Removes suffix if matching else return undef
        my $outfile_path_no_suffix = parse_file_suffix(
            {
                file_name   => $outfile_path,
                file_suffix => $DOT . q{gz},
            }
        );

        gzip(
            {
                decompress   => 1,
                FILEHANDLE   => $FILEHANDLE,
                force        => 1,
                infile_path  => $outfile_path,
                outfile_path => $outfile_path_no_suffix,
                quiet        => 1,
                stdout       => 1,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    if ( $decompress_program eq q{unzip} ) {

        unzip(
            {
                FILEHANDLE  => $FILEHANDLE,
                force       => 1,
                infile_path => $outfile_path,
                outdir_path => $outdir_path,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    if ( $decompress_program eq q{tar} ) {

        tar(
            {
                extract           => 1,
                FILEHANDLE        => $FILEHANDLE,
                file_path         => $outfile_path,
                filter_gzip       => 1,
                outdirectory_path => $outdir_path,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }
    return 1;
}

1;
