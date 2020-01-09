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
    our $VERSION = 1.01;

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
## Returns  : %bed_file_path
## Arguments: $contigs_ref           => Contigs to split in file
##          : $exome_target_bed_file => Interval file to split
##          : $file_ending           => File ending to add {Optional}
##          : $filehandle            => Filehandle to write to
##          : $max_cores_per_node    => Maximum core per node
##          : $outdirectory          => Outdirectory
##          : $reference_dir         => MIP reference directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $contigs_ref;
    my $exome_target_bed_file;
    my $file_ending;
    my $filehandle;
    my $outdirectory;
    my $reference_dir;

    ## Default(s)
    my $max_cores_per_node;

    my $tmpl = {
        contigs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$contigs_ref,
            strict_type => 1,
        },
        exome_target_bed_file => {
            defined     => 1,
            required    => 1,
            store       => \$exome_target_bed_file,
            strict_type => 1,
        },
        file_ending => { strict_type => 1, store => \$file_ending },
        filehandle         => { store => \$filehandle, },
        max_cores_per_node => {
            allow       => qr/ ^\d+$ /sxm,
            default     => 1,
            store       => \$max_cores_per_node,
            strict_type => 1,
        },
        outdirectory => {
            defined     => 1,
            required    => 1,
            store       => \$outdirectory,
            strict_type => 1,
        },
        reference_dir => {
            defined     => 1,
            required    => 1,
            store       => \$reference_dir,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Processmanagement::Processes qw(print_wait);

    my %bed_file_path;
    my $process_batches_count = 1;

    say {$filehandle} q{## Generate contig specific interval_list}, $NEWLINE;

  CONTIG:
    while ( my ( $contig_index, $contig ) = each @{$contigs_ref} ) {

        $process_batches_count = print_wait(
            {
                process_counter       => $contig_index,
                max_process_number    => $max_cores_per_node,
                process_batches_count => $process_batches_count,
                filehandle            => $filehandle,
            }
        );

        ## Splits a target file into new contig specific target file
        my $contig_bed_file_path = _split_interval_file_contigs(
            {
                filehandle   => $filehandle,
                file_ending  => $file_ending,
                contig       => $contig,
                indirectory  => $reference_dir,
                infile       => basename($exome_target_bed_file),
                outdirectory => $outdirectory,
            }
        );
        $bed_file_path{$contig} = [$contig_bed_file_path];
    }
    say {$filehandle} q{wait}, $NEWLINE;

    return %bed_file_path;
}

sub _split_interval_file_contigs {

## Function : Splits a target file into new contig specific target file
## Returns  : $outfile_path
## Arguments: $indirectory  => Indirectory
##          : $outdirectory => Outdirectory
##          : $infile       => Target file
##          : $contig       => Contig to extract
##          : $file_ending  => File ending to add
##          : $filehandle   => filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $indirectory;
    my $outdirectory;
    my $infile;
    my $contig;
    my $file_ending;
    my $filehandle;

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
        filehandle => { required => 1, defined => 1, store => \$filehandle, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $outfile_path;

    ## The contig to split
    if ( defined $contig ) {

        ### Select header and contig
        ## Execute perl
        print {$filehandle} q?perl -nae '?;

        ## If header line
        print {$filehandle} q?if($_=~/^\@/) { ?;

        ## Write to STDOUT
        print {$filehandle} q?print $_;} ?;

        ## If line begin with specific $contig
        print {$filehandle} q?elsif($_=~/^? . $contig . q?\s+/) { ?;

        ## Write to STDOUT
        print {$filehandle} q?print $_;}' ?;

        ## Infile
        print {$filehandle} catfile( $indirectory, $infile ) . $SPACE;

        $outfile_path =
          catfile( $outdirectory, $contig . $UNDERSCORE . $infile );

        ## If file ending
        if ( defined $file_ending ) {

            ## Add supplied file ending
            $outfile_path .= $file_ending;
        }
        say {$filehandle} q{>} . $SPACE . $outfile_path . $SPACE . $AMPERSAND;
    }
    return $outfile_path;
}

1;
