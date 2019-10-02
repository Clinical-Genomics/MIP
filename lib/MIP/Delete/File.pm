package MIP::Delete::File;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;

## CPANM
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $AMPERSAND $DOT $EMPTY_STR $NEWLINE $SPACE };
use MIP::Gnu::Coreutils qw{ gnu_rm };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ delete_contig_files delete_files };
}

sub delete_contig_files {

## Function : Delete contig files dictated by supplied array index and element.
## Returns  :
## Arguments: $core_number       => Number of cores to use
##          : $FILEHANDLE        => Sbatch filehandle to write to
##          : $file_elements_ref => Array to use for file iteration {REF}
##          : $file_ending       => File ending
##          : $file_name         => File name without ending attached
##          : $indirectory       => Indirectory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $core_number;
    my $FILEHANDLE;
    my $file_elements_ref;
    my $file_ending;
    my $file_name;
    my $indirectory;

    my $tmpl = {
        core_number => {
            defined     => 1,
            required    => 1,
            store       => \$core_number,
            strict_type => 1,
        },
        FILEHANDLE        => { defined => 1, required => 1, store => \$FILEHANDLE, },
        file_elements_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$file_elements_ref,
            strict_type => 1,
        },
        file_ending => {
            defined     => 1,
            required    => 1,
            store       => \$file_ending,
            strict_type => 1,
        },
        file_name => {
            defined     => 1,
            required    => 1,
            store       => \$file_name,
            strict_type => 1,
        },
        indirectory => {
            defined     => 1,
            required    => 1,
            store       => \$indirectory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Gnu::Coreutils qw{ gnu_rm };
    use MIP::Processmanagement::Processes qw{ print_wait };

    my $process_batches_count = 1;

    ## Remove infile at indirectory
    say {$FILEHANDLE} q{## Remove file at directory};

    while ( my ( $index, $element ) = each @{$file_elements_ref} ) {

        $process_batches_count = print_wait(
            {
                FILEHANDLE            => $FILEHANDLE,
                max_process_number    => $core_number,
                process_batches_count => $process_batches_count,
                process_counter       => $index,
            }
        );

        gnu_rm(
            {
                FILEHANDLE => $FILEHANDLE,
                force      => 1,
                infile_path =>
                  catfile( $indirectory, $file_name . $DOT . $element . $file_ending ),
            }
        );
        say {$FILEHANDLE} $AMPERSAND . $SPACE;
    }
    say {$FILEHANDLE} q{wait}, $NEWLINE;
    return;
}

sub delete_files {

## Function : Delete files
## Returns  :
## Arguments: $core_number => Number of cores that can be used
##          : $FILEHANDLE  => Filehandle to write to
##          : $file_ending => File ending for infiles. {Optional}
##          : $indirectory => The directory for the files to be copied
##          : $infiles_ref => The array of files to copy

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infiles_ref;
    my $FILEHANDLE;
    my $indirectory;
    my $core_number;

    ## Default(s)
    my $file_ending;

    my $tmpl = {
        core_number => {
            defined     => 1,
            required    => 1,
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Processmanagement::Processes qw{ print_wait };

    my $process_batches_count = 1;

    say {$FILEHANDLE} q{## Remove file(s)};
  FILE:
    while ( my ( $file_index, $file ) = each @{$infiles_ref} ) {

        $process_batches_count = print_wait(
            {
                FILEHANDLE            => $FILEHANDLE,
                max_process_number    => $core_number,
                process_batches_count => $process_batches_count,
                process_counter       => $file_index,
            }
        );

        ## Remove file
        gnu_rm(
            {
                FILEHANDLE  => $FILEHANDLE,
                force       => 1,
                infile_path => catfile( $indirectory, $file . $file_ending ),
            }
        );
        say {$FILEHANDLE} $AMPERSAND;
    }
    say {$FILEHANDLE} q{wait}, $NEWLINE;
    return;
}

1;
