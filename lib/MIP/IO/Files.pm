package MIP::IO::Files;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use autodie;
use Params::Check qw{ check allow last_error };
use File::Spec::Functions qw{ catdir catfile splitpath };
use FindBin qw{ $Bin };
use File::Basename qw{ dirname };

##CPANM
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{ lib } );
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ migrate_file migrate_files xargs_migrate_contig_files };

}

##Constants
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
## Arguments: $infiles_ref  => Array of files to copy {REF}
##          : $outfile_path => Outfile path
##          : $FILEHANDLE   => Filehandle to write to
##          : $indirectory  => The directory for the files to be copied
##          : $core_number  => The number of cores that can be used
##          : $file_ending  => File ending for infiles

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infiles_ref;
    my $outfile_path;
    my $FILEHANDLE;
    my $indirectory;
    my $core_number;

    ## Default(s)
    my $file_ending;

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
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$core_number
        },
        file_ending => {
            default     => $EMPTY_STR,
            strict_type => 1,
            store       => \$file_ending
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
    say {$FILEHANDLE} q{wait}, $NEWLINE;

    return;
}

sub migrate_file {

## Function : Copy file to from source ($infile_path) to destination ($outfile_path).
## Returns  : $infile_path_file_name
## Arguments: $FILEHANDLE      => Filehandle to write to
##          : $infile_path     => Infile path
##          : $outfile_path    => Outfile path
##          : $stderrfile_path => Stderrfile path
##          : $xargs           => Use xargs if defined {Optional}

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

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Gnu::Coreutils qw{ gnu_cp };

    ## Split relative infile_path to file(s)
    my ( $infile_path_volume, $infile_path_directory, $infile_path_file_name )
      = splitpath($infile_path);

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
##          : $FILEHANDLE         => Sbatch filehandle to write to
##          : $XARGSFILEHANDLE    => XARGS filehandle to write to
##          : $file_path          => File name
##          : $temp_directory     => Temporary directory
##          : $program_info_path  => The program info path
##          : $infile             => Infile name without ending attached
##          : $indirectory        => In directory
##          : $outfile            => OutFile name without ending attached
##          : $outdirectory       => Out directory
##          : $core_number        => The number of cores to use
##          : $first_command      => The inital command
##          : $xargs_file_counter => The xargs file counter
##          : $file_ending        => File ending

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $contigs_ref;
    my $FILEHANDLE;
    my $XARGSFILEHANDLE;
    my $file_path;
    my $temp_directory;
    my $program_info_path;
    my $infile;
    my $indirectory;
    my $outfile;
    my $outdirectory;
    my $first_command;

    ## Default(s)
    my $file_ending;
    my $xargs_file_counter;
    my $core_number;

    my $tmpl = {
        contigs_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$contigs_ref
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
        XARGSFILEHANDLE =>
          { required => 1, defined => 1, store => \$XARGSFILEHANDLE },
        file_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_path
        },
        temp_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$temp_directory
        },
        program_info_path => { strict_type => 1, store => \$program_info_path },
        infile            => { strict_type => 1, store => \$infile },
        indirectory       => { strict_type => 1, store => \$indirectory },
        outfile           => { strict_type => 1, store => \$outfile },
        outdirectory      => { strict_type => 1, store => \$outdirectory },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter
        },
        core_number => {
            default     => 1,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$core_number
        },
        first_command => { strict_type => 1, store => \$first_command },
        file_ending   => {
            default     => $DOT . q{vcf} . $ASTERIX,
            strict_type => 1,
            store       => \$file_ending
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };

    ## Create file commands for xargs
    ( $xargs_file_counter, my $xargs_file_path_prefix ) = xargs_command(
        {
            FILEHANDLE         => $FILEHANDLE,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            file_path          => $file_path,
            program_info_path  => $program_info_path,
            core_number        => $core_number,
            xargs_file_counter => $xargs_file_counter,
            first_command      => $first_command,
        }
    );

  CONTIGS:
    foreach my $contig ( @{$contigs_ref} ) {

        my $stderrfile_path =
          $xargs_file_path_prefix . $DOT . $contig . $DOT . q{stderr.txt};

        if ( defined $infile ) {

            ## Get parameters
            my $infile_path = catfile( $indirectory,
                $infile . $UNDERSCORE . $contig . $file_ending );

            ## Copy file(s) to temporary directory.
            migrate_file(
                {
                    FILEHANDLE      => $XARGSFILEHANDLE,
                    infile_path     => $infile_path,
                    outfile_path    => $temp_directory,
                    xargs           => q{xargs},
                    stderrfile_path => $stderrfile_path,
                }
            );
        }
        if ( defined $outfile && defined $outdirectory ) {

            ## Get parameters
            my $infile_path = catfile( $temp_directory,
                $outfile . $UNDERSCORE . $contig . $file_ending );
            ## Copy file(s) from temporary directory.
            migrate_file(
                {
                    infile_path     => $infile_path,
                    outfile_path    => $outdirectory,
                    FILEHANDLE      => $XARGSFILEHANDLE,
                    xargs           => q{xargs},
                    stderrfile_path => $stderrfile_path,
                }
            );
        }
    }
    return $xargs_file_counter;
}

1;
