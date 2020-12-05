package MIP::File::Interval;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

## MIPs lib/
use MIP::Constants qw{ $AMPERSAND $NEWLINE $SPACE $UNDERSCORE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ generate_contig_interval_file };
}

sub generate_contig_interval_file {

## Function : Generate contig specific interval list files
## Returns  : %bed_file_path
## Arguments: $contigs_ref           => Contigs to split in files
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
        file_ending        => { strict_type => 1, store => \$file_ending },
        filehandle         => { store       => \$filehandle, },
        max_cores_per_node => {
            allow       => qr/ \A \d+ \z /sxm,
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

    use MIP::Processmanagement::Processes qw{ print_wait };

    my %bed_file_path;
    my $process_batches_count = 1;

    say {$filehandle} q{## Generate contig specific interval_list}, $NEWLINE;

  CONTIG:
    while ( my ( $contig_index, $contig ) = each @{$contigs_ref} ) {

        $process_batches_count = print_wait(
            {
                filehandle            => $filehandle,
                max_process_number    => $max_cores_per_node,
                process_batches_count => $process_batches_count,
                process_counter       => $contig_index,
            }
        );

        ## Splits a target file into new contig specific target file
        my $contig_bed_file_path = _split_interval_file_contigs(
            {
                contig       => $contig,
                file_ending  => $file_ending,
                filehandle   => $filehandle,
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
## Arguments: $file_ending  => File ending to add
##          : $contig       => Contig to extract
##          : $filehandle   => filehandle to write to
##          : $indirectory  => Indirectory
##          : $infile       => Target file
##          : $outdirectory => Outdirectory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $contig;
    my $file_ending;
    my $filehandle;
    my $indirectory;
    my $infile;
    my $outdirectory;

    my $tmpl = {
        contig      => { store   => \$contig, },
        file_ending => { store   => \$file_ending, strict_type => 1, },
        filehandle  => { defined => 1, required => 1, store => \$filehandle, },
        indirectory => {
            defined     => 1,
            required    => 1,
            store       => \$indirectory,
            strict_type => 1,
        },
        infile => {
            defined     => 1,
            required    => 1,
            store       => \$infile,
            strict_type => 1,
        },
        outdirectory => {
            defined     => 1,
            required    => 1,
            store       => \$outdirectory,
            strict_type => 1,
        },
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

        $outfile_path = catfile( $outdirectory, $contig . $UNDERSCORE . $infile );

        if ( defined $file_ending ) {

            $outfile_path .= $file_ending;
        }
        say {$filehandle} q{>} . $SPACE . $outfile_path . $SPACE . $AMPERSAND;
    }
    return $outfile_path;
}

1;
