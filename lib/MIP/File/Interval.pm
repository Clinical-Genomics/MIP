package MIP::File::Interval;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use File::Basename qw{ basename };
use File::Spec::Functions qw{ catfile };

## CPANM
use Readonly;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ generate_contig_interval_file };
}

## Constants
Readonly my $AMPERSAND  => q{&};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub generate_contig_interval_file {

## Function : Generate contig specific interval_list file
## Returns  :
## Arguments: $contigs_ref           => Contigs to split in file
##          : $exome_target_bed_file => Interval file to split
##          : $reference_dir         => MIP reference directory
##          : $outdirectory          => Outdirectory
##          : $file_ending           => File ending to add {Optional}
##          : $FILEHANDLE            => Filehandle to write to
##          : $max_cores_per_node    => Maximum core per node

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $contigs_ref;
    my $exome_target_bed_file;
    my $reference_dir;
    my $outdirectory;
    my $file_ending;
    my $FILEHANDLE;

    ## Default(s)
    my $max_cores_per_node;

    my $tmpl = {
        contigs_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$contigs_ref,
        },
        exome_target_bed_file => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$exome_target_bed_file,
        },
        reference_dir => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$reference_dir,
        },
        outdirectory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outdirectory,
        },
        file_ending => { strict_type => 1, store => \$file_ending },
        FILEHANDLE         => { store => \$FILEHANDLE, },
        max_cores_per_node => {
            default     => 1,
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$max_cores_per_node,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Processmanagement::Processes qw(print_wait);

    my $process_batches_count = 1;

    say {$FILEHANDLE} q{## Generate contig specific interval_list}, $NEWLINE;

  CONTIG:
    while ( my ( $contig_index, $contig ) = each @{$contigs_ref} ) {

        $process_batches_count = print_wait(
            {
                process_counter       => $contig_index,
                max_process_number    => $max_cores_per_node,
                process_batches_count => $process_batches_count,
                FILEHANDLE            => $FILEHANDLE,
            }
        );

        ## Splits a target file into new contig specific target file
        _split_interval_file_contigs(
            {
                FILEHANDLE   => $FILEHANDLE,
                indirectory  => $reference_dir,
                outdirectory => $outdirectory,
                infile       => basename($exome_target_bed_file),
                contig       => $contig,
                file_ending  => $file_ending,
            }
        );
    }
    say {$FILEHANDLE} q{wait}, $NEWLINE;
    return;
}

sub _split_interval_file_contigs {

## Function : Splits a target file into new contig specific target file
## Returns  :
## Arguments: $indirectory  => Indirectory
##          : $outdirectory => Outdirectory
##          : $infile       => Target file
##          : $contig       => Contig to extract
##          : $file_ending  => File ending to add
##          : $FILEHANDLE   => FILEHANDLE to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $indirectory;
    my $outdirectory;
    my $infile;
    my $contig;
    my $file_ending;
    my $FILEHANDLE;

    my $tmpl = {
        indirectory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$indirectory,
        },
        outdirectory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outdirectory,
        },
        infile => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile,
        },
        contig      => { store       => \$contig },
        file_ending => { strict_type => 1, store => \$file_ending },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## The contig to split
    if ( defined $contig ) {

        ### Select header and contig
        ## Execute perl
        print {$FILEHANDLE} q?perl -nae '?;

        ## If header line
        print {$FILEHANDLE} q?if($_=~/^\@/) { ?;

        ## Write to STDOUT
        print {$FILEHANDLE} q?print $_;} ?;

        ## If line begin with specific $contig
        print {$FILEHANDLE} q?elsif($_=~/^? . $contig . q?\s+/) { ?;

        ## Write to STDOUT
        print {$FILEHANDLE} q?print $_;}' ?;

        ## Infile
        print {$FILEHANDLE} catfile( $indirectory, $infile ) . $SPACE;

        my $outfile_path =
          catfile( $outdirectory, $contig . $UNDERSCORE . $infile );

        ## If file ending
        if ( defined $file_ending ) {

            ## Add supplied file ending
            $outfile_path .= $file_ending;
        }
        say {$FILEHANDLE} q{>} . $SPACE . $outfile_path . $SPACE . $AMPERSAND;
    }
    return;
}

1;
