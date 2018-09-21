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
use lib catdir( dirname($Bin), q{ lib } );
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ migrate_file migrate_files xargs_migrate_contig_files };

}

## Constants
Readonly my $AMPERSAND  => q{&};
Readonly my $ASTERIX    => q{*};
Readonly my $DOT        => q{.};
Readonly my $EMPTY_STR  => q{};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub migrate_files {

## Function : Copies files from source to destination.
## Returns  :
## Arguments: $core_number  => The number of cores that can be used
##          : $FILEHANDLE   => Filehandle to write to
##          : $file_ending  => File ending for infiles
##          : $infiles_ref  => Array of files to copy {REF}
##          : $indirectory  => The directory for the files to be copied
##          : $outfile_path => Outfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $core_number;
    my $FILEHANDLE;
    my $indirectory;
    my $infiles_ref;
    my $outfile_path;

    ## Default(s)
    my $file_ending;

    my $tmpl = {
        core_number => {
            allow       => qr/ ^\d+$ /xsm,
            default     => 1,
            store       => \$core_number,
            strict_type => 1,
        },
        FILEHANDLE  => { defined => 1, required => 1, store => \$FILEHANDLE, },
        file_ending => {
            default     => $EMPTY_STR,
            store       => \$file_ending,
            strict_type => 1,
        },
        indirectory => {
            defined     => 1,
            required    => 1,
            store       => \$indirectory,
            strict_type => 1,
        },
        infiles_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infiles_ref,
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

    use MIP::Processmanagement::Processes qw{ print_wait };

    my $process_batches_count = 1;

    say {$FILEHANDLE} q{## Copying file(s) to destination directory};

  INFILES:
    while ( my ( $file_index, $file ) = each @{$infiles_ref} ) {

        $process_batches_count = print_wait(
            {
                FILEHANDLE            => $FILEHANDLE,
                max_process_number    => $core_number,
                process_batches_count => $process_batches_count,
                process_counter       => $file_index,
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
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    return;
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
    my ( $infile_path_volume, $infile_path_directory, $infile_path_file_name )
      = splitpath($infile_path);

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

sub xargs_migrate_contig_files {

## Function : Migrates file(s) to or from temporary directory (depending on supplied arguments) using xargs.
## Returns  : $xargs_file_counter
## Arguments: $contigs_ref        => Contigs to iterate over {REF}
##          : $core_number        => The number of cores to use
##          : $FILEHANDLE         => Sbatch filehandle to write to
##          : $file_ending        => File ending
##          : $file_path          => File name
##          : $first_command      => The inital command
##          : $indirectory        => In directory
##          : $infile             => Infile name without suffix attached
##          : $outdirectory       => Outdirectory
##          : $outfile            => Outfile name without suffix attached
##          : $program_info_path  => Program info path
##          : $temp_directory     => Temporary directory
##          : $XARGSFILEHANDLE    => XARGS filehandle to write to
##          : $xargs_file_counter => Xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $contigs_ref;
    my $FILEHANDLE;
    my $file_path;
    my $first_command;
    my $indirectory;
    my $infile;
    my $outdirectory;
    my $outfile;
    my $program_info_path;
    my $temp_directory;
    my $XARGSFILEHANDLE;

    ## Default(s)
    my $core_number;
    my $file_ending;
    my $xargs_file_counter;

    my $tmpl = {
        contigs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$contigs_ref,
            strict_type => 1,
        },
        core_number => {
            allow       => qr/ ^\d+$ /xsm,
            default     => 1,
            store       => \$core_number,
            strict_type => 1,
        },
        FILEHANDLE  => { defined => 1, required => 1, store => \$FILEHANDLE },
        file_ending => {
            default     => $DOT . q{vcf} . $ASTERIX,
            store       => \$file_ending,
            strict_type => 1,
        },
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        first_command => { store => \$first_command, strict_type => 1, },
        indirectory   => { store => \$indirectory,   strict_type => 1, },
        infile        => { store => \$infile,        strict_type => 1, },
        outdirectory  => { store => \$outdirectory,  strict_type => 1, },
        outfile       => { store => \$outfile,       strict_type => 1, },
        program_info_path =>
          { store => \$program_info_path, strict_type => 1, },
        temp_directory => {
            defined     => 1,
            required    => 1,
            store       => \$temp_directory,
            strict_type => 1,
        },
        XARGSFILEHANDLE =>
          { defined => 1, required => 1, store => \$XARGSFILEHANDLE, },
        xargs_file_counter => {
            allow       => qr/ ^\d+$ /xsm,
            default     => 0,
            store       => \$xargs_file_counter,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };

    ## Create file commands for xargs
    ( $xargs_file_counter, my $xargs_file_path_prefix ) = xargs_command(
        {
            core_number        => $core_number,
            FILEHANDLE         => $FILEHANDLE,
            file_path          => $file_path,
            first_command      => $first_command,
            program_info_path  => $program_info_path,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

  CONTIGS:
    foreach my $contig ( @{$contigs_ref} ) {

        my $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};

        if ( defined $infile ) {

            ## Get parameters
            my $infile_path =
              catfile( $indirectory, $infile . $DOT . $contig . $file_ending );

            ## Copy file(s) to temporary directory.
            migrate_file(
                {
                    FILEHANDLE      => $XARGSFILEHANDLE,
                    infile_path     => $infile_path,
                    outfile_path    => $temp_directory,
                    stderrfile_path => $stderrfile_path,
                    xargs           => q{xargs},
                }
            );
        }
        if ( defined $outfile && defined $outdirectory ) {

            ## Get parameters
            my $infile_path = catfile( $temp_directory,
                $outfile . $DOT . $contig . $file_ending );
            ## Copy file(s) from temporary directory.
            migrate_file(
                {
                    FILEHANDLE      => $XARGSFILEHANDLE,
                    infile_path     => $infile_path,
                    outfile_path    => $outdirectory,
                    stderrfile_path => $stderrfile_path,
                    xargs           => q{xargs},
                }
            );
        }
    }
    return $xargs_file_counter;
}

1;
