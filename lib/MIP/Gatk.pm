package MIP::Gatk;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib
use MIP::Constants qw{ $LOG_NAME $NEWLINE $TAB };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_gatk_sample_map_paths get_gatk_intervals };
}

sub check_gatk_sample_map_paths {

## Function : Check that the supplied gatk sample map file paths exists
## Returns  :
## Arguments: $sample_map_path => Sample map path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_map_path;

    my $tmpl = {
        sample_map_path => {
            defined     => 1,
            required    => 1,
            store       => \$sample_map_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Io::Read qw{ read_from_file };

    ## Constants
    Readonly my $FIELD_COUNTER => 2;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my @missing_files;

    my @lines = read_from_file(
        {
            chomp  => 1,
            format => q{line_by_line},
            path   => $sample_map_path,
        }
    );

  LINE:
    foreach my $line (@lines) {

        ## Get sample and file path (and check proper format)
        my ( $sample, $file_path, $unexpected_data ) =
          split $TAB, $line, $FIELD_COUNTER + 1;

        ## Make sure that we get what we expect
        if ( defined $unexpected_data ) {

            $log->logcarp( q{Unexpected trailing garbage at end of line '} . $line . q{':},
                $NEWLINE . $TAB . $unexpected_data . $NEWLINE );
        }

        ## Path exists
        next LINE if ( -e $file_path );

        my $error_msg =
            q{The supplied file path: }
          . $file_path
          . q{ from sample map file: }
          . $sample_map_path
          . q{ does not exist};
        push @missing_files, $error_msg;
    }
    return 1 if ( not @missing_files );

    ## Broadcast missing files
    foreach my $error_msg (@missing_files) {
        $log->fatal($error_msg);
    }
    exit 1;
}

sub get_gatk_intervals {

## Function : Generate and return interval hash
## Returns  : %gatk_intervals
## Arguments: $analysis_type         => Analysis type
##          : $contigs_ref           => Contigs to split in file
##          : $exome_target_bed_href => Exome target bed files lnked to sample ids
##          : $file_ending           => File ending to add {Optional}
##          : $filehandle            => Filehandle to write to
##          : $max_cores_per_node    => Maximum core per node
##          : $outdirectory          => Outdirectory
##          : $reference_dir         => MIP reference directory
##          : $sample_id             => Sample_id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $analysis_type;
    my $contigs_ref;
    my $exome_target_bed_href;
    my $file_ending;
    my $filehandle;
    my $outdirectory;
    my $reference_dir;
    my $sample_id;

    ## Default(s)
    my $max_cores_per_node;

    my $tmpl = {
        analysis_type => {
            defined     => 1,
            required    => 1,
            store       => \$analysis_type,
            strict_type => 1,
        },
        contigs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$contigs_ref,
            strict_type => 1,
        },
        exome_target_bed_href => {
            default     => {},
            store       => \$exome_target_bed_href,
            strict_type => 1,
        },
        file_ending => {
            store       => \$file_ending,
            strict_type => 1,
        },
        filehandle         => { store => \$filehandle, },
        max_cores_per_node => {
            allow       => qr/ \A \d+ \z /sxm,
            default     => 1,
            store       => \$max_cores_per_node,
            strict_type => 1,
        },
        outdirectory => {
            store       => \$outdirectory,
            strict_type => 1,
        },
        reference_dir => {
            store       => \$reference_dir,
            strict_type => 1,
        },
        sample_id => {
            store       => \$sample_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Active_parameter qw{ get_exome_target_bed_file };
    use MIP::Contigs qw{ generate_contig_interval_file };

    ## Store gatk interval as contigs
    my %gatk_intervals = map { $_ => [$_] } @{$contigs_ref};

    return %gatk_intervals if ( not $analysis_type eq q{wes} );

    my $exome_target_bed_file = get_exome_target_bed_file(
        {
            exome_target_bed_href => $exome_target_bed_href,
            file_ending           => $file_ending,
            sample_id             => $sample_id,
        }
    );

    ## Generate contig specific interval_list and return gatk interval as contig number with paths
    %gatk_intervals = generate_contig_interval_file(
        {
            contigs_ref           => $contigs_ref,
            exome_target_bed_file => $exome_target_bed_file,
            filehandle            => $filehandle,
            max_process_number    => $max_cores_per_node,
            outdirectory          => $outdirectory,
            reference_dir         => $reference_dir,
        }
    );
    return %gatk_intervals;
}

1;
