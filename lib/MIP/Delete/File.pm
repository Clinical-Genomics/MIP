package MIP::Delete::File;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use File::Basename qw{ dirname };

## CPANM
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Gnu::Coreutils qw(gnu_rm);

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ delete_files };
}

## Constants
Readonly my $AMPERSAND => q{&};
Readonly my $EMPTY_STR => q{};
Readonly my $NEWLINE   => qq{\n};
Readonly my $SPACE     => q{ };

sub delete_files {

##delete_files

##Function : Delete files.
##Returns  : ""
##Arguments: $infiles_ref, $FILEHANDLE, $indirectory, $core_number, $file_ending
##         : $infiles_ref => The array of files to copy
##         : $FILEHANDLE  => Filehandle to write to
##         : $indirectory => The directory for the files to be copied
##         : $core_number => The number of cores that can be used
##         : $file_ending => File ending for infiles. {Optional}

    my ($arg_href) = @_;

    ## Default(s)
    my $file_ending;

    ## Flatten argument(s)
    my $infiles_ref;
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
        FILEHANDLE  => { required => 1, defined => 1, store => \$FILEHANDLE },
        indirectory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$indirectory
        },
        core_number => {
            required    => 1,
            defined     => 1,
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

    say {$FILEHANDLE} q{## Remove file(s)};
  FILE:
    while ( my ( $file_index, $file ) = each @{$infiles_ref} ) {

        $process_batches_count = print_wait(
            {
                process_counter       => $file_index,
                max_process_number    => $core_number,
                process_batches_count => $process_batches_count,
                FILEHANDLE            => $FILEHANDLE,
            }
        );

        ## Remove file
        gnu_rm(
            {
                infile_path => catfile( $indirectory, $file . $file_ending ),
                FILEHANDLE  => $FILEHANDLE,
                force       => 1,
            }
        );
        say {$FILEHANDLE} $AMPERSAND;
    }
    say {$FILEHANDLE} q{wait}, $NEWLINE;
    return;
}

1;
