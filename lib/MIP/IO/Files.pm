package MIP::IO::Files;

use strict;
use warnings;
use warnings qw( FATAL utf8 );
use utf8;    # Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use autodie;
use Params::Check qw[check allow last_error];
use File::Spec::Functions qw(catfile);

use FindBin qw($Bin);    # Find directory of script
use File::Basename qw(dirname);
use File::Spec::Functions qw(catdir);

##CPANM
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use MIP::Unix::Standard_streams qw(unix_standard_streams);
use MIP::Unix::Write_to_file qw(unix_write_to_file);

BEGIN {
    use base qw (Exporter);
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(migrate_file migrate_files);

}

##Constants
Readonly my $EMPTY_STR => q{};

sub migrate_files {

##migrate_files

##Function : Copies files from source to destination.
##Returns  : ""
##Arguments: $infiles_ref, $outfile_path, $FILEHANDLE, $indirectory, $core_number, $file_ending
##         : $infiles_ref  => Array of files to copy
##         : $outfile_path => Outfile path
##         : $FILEHANDLE   => Filehandle to write to
##         : $indirectory  => The directory for the files to be copied
##         : $core_number  => The number of cores that can be used
##         : $file_ending  => File ending for infiles. {Optional}

    my ($arg_href) = @_;

    ## Default(s)
    my $file_ending;

    ## Flatten argument(s)
    my $infiles_ref;
    my $outfile_path;
    my $FILEHANDLE;
    my $indirectory;
    my $core_number;

    my $tmpl = {
        infiles_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$infiles_ref
        },
        outfile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfile_path
        },
        FILEHANDLE  => { required => 1, defined => 1, store => \$FILEHANDLE },
        indirectory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$indirectory
        },
        core_number => {
            default     => 1,
            allow       => qr/^\d+$/,
            strict_type => 1,
            store       => \$core_number
        },
        file_ending => {
            default     => $EMPTY_STR,
            strict_type => 1,
            store       => \$file_ending
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    use MIP::Processmanagement::Processes qw(print_wait);

    my $process_batches_count = 1;

    say $FILEHANDLE q{## Copying file(s) to destination directory};

  INFILES:
    while ( my ( $file_index, $file ) = each $infiles_ref ) {

        $process_batches_count = print_wait(
            {
                process_counter       => $file_index,
                max_process_number    => $core_number,
                process_batches_count => $process_batches_count,
                FILEHANDLE            => $FILEHANDLE,
            }
        );

        ## Copies file to destination
        migrate_file(
            {
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => catfile( $indirectory, $file . $file_ending ),
                outfile_path => $outfile_path,
            }
        );
    }
    say $FILEHANDLE q{wait}, "\n";

    return;
}

sub migrate_file {

##migrate_file

##Function : Copy file to from source ($infile_path) to destination ($outfile_path).
##Returns  : "$infile_path_file_name"
##Arguments: $FILEHANDLE, $infile_path, $outfile_path, $stderrfile_path, $xargs
##         : $FILEHANDLE      => Filehandle to write to
##         : $infile_path     => Infile path
##         : $outfile_path    => Outfile path
##         : $stderrfile_path => Stderrfile path
##         : $xargs           => Use xargs if defined {Optional}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $xargs;

    my $tmpl = {
        FILEHANDLE  => { required => 1, defined => 1, store => \$FILEHANDLE },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        outfile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfile_path
        },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        xargs           => { strict_type => 1, store => \$xargs },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    use MIP::Gnu::Coreutils qw(gnu_cp);

    ## Split relative infile_path to file(s)
    my ( $infile_path_volume, $infile_path_directory, $infile_path_file_name )
      = File::Spec->splitpath($infile_path);

    gnu_cp(
        {
            FILEHANDLE      => $FILEHANDLE,
            preserve        => 1,
            infile_path     => $infile_path,
            outfile_path    => $outfile_path,
            stderrfile_path => $stderrfile_path,
        }
    );

    # For print wait statement downstream
    if ( !defined $xargs ) {

        say $FILEHANDLE q{& };
    }
    print $FILEHANDLE "\n";

    return $infile_path_file_name;
}

1;
