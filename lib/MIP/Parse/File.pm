package MIP::Parse::File;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ parse_fastq_infiles parse_fastq_infiles_format };

}

## Constants
Readonly my $EMPTY_STR => q{};

sub parse_fastq_infiles {

## Function : Parse fastq infiles for MIP processing.
## Returns  :
## Arguments: $active_parameter_href           => Active parameters for this analysis hash {REF}
##          : $file_info_href                  => File info hash {REF}
##          : $infile_both_strands_prefix_href => The infile(s) without the ".ending" and strand info {REF}
##          : $infile_lane_prefix_href         => Infile(s) without the ".ending" {REF}
##          : $log                             => Log object
##          : $sample_info_href                => Info on samples and family hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_both_strands_prefix_href;
    my $infile_lane_prefix_href;
    my $log;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        infile_both_strands_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_both_strands_prefix_href,
            strict_type => 1,
        },
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_lane_prefix_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::File qw{ check_interleaved };
    use MIP::Check::Parameter qw{ check_infile_contain_sample_id };
    use MIP::Get::File qw{ get_fastq_file_header_info get_read_length };
    use MIP::Set::File qw{ set_file_compression_features };
    use MIP::QC::Record qw{ add_infile_info };

  SAMPLE_ID:
    for my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        # Needed to be able to track when lanes are finished
        my $lane_tracker = 0;

        ## Unpack
        my $infiles_dir = $file_info_href->{$sample_id}{mip_infiles_dir};

      INFILE:
        while ( my ( $file_index, $file_name ) =
            each @{ $file_info_href->{$sample_id}{mip_infiles} } )
        {

            # Sequence read length
            my $read_length;

            # Is file interleaved
            my $is_interleaved;

            ## Set compression features
            my ( $is_file_compressed, $read_file_command ) =
              set_file_compression_features( { file_name => $file_name, } );

            ## If not compressed
            if ( not $is_file_compressed ) {

                ## Note: All files are rechecked downstream and uncompressed ones are gzipped automatically
                $file_info_href->{is_file_uncompressed}{$sample_id}++;
            }

            ## Parse infile according to filename convention
            my %infile_info =
              parse_fastq_infiles_format( { file_name => $file_name, } );

            ## Get read length and interleaved status from file
            if ( exists $infile_info{direction}
                and $infile_info{direction} == 1 )
            {

                ## Get sequence read length from file
                $read_length = get_read_length(
                    {
                        file_path => catfile( $infiles_dir, $file_name ),
                        read_file_command => $read_file_command,
                    }
                );

                ## Is file interleaved and have proper read direction
                $is_interleaved = check_interleaved(
                    {
                        file_path => catfile( $infiles_dir, $file_name ),
                        log       => $log,
                        read_file_command => $read_file_command,
                    }
                );
            }

            ## If filename convention is followed
            if ( keys %infile_info ) {

                ## Check that the sample_id provided and sample_id in infile name match.
                check_infile_contain_sample_id(
                    {
                        infile_name      => $file_name,
                        infile_sample_id => $infile_info{infile_sample_id},
                        log              => $log,
                        sample_id        => $sample_id,
                        sample_ids_ref =>
                          \@{ $active_parameter_href->{sample_ids} },
                    }
                );

                ## Adds information derived from infile name to hashes
                $lane_tracker = add_infile_info(
                    {
                        active_parameter_href => $active_parameter_href,
                        date                  => $infile_info{date},
                        direction             => $infile_info{direction},
                        file_info_href        => $file_info_href,
                        file_index            => $file_index,
                        flowcell              => $infile_info{flowcell},
                        index                 => $infile_info{index},
                        infile_both_strands_prefix_href =>
                          $infile_both_strands_prefix_href,
                        infile_lane_prefix_href => $infile_lane_prefix_href,
                        is_interleaved          => $is_interleaved,
                        lane                    => $infile_info{lane},
                        lane_tracker            => $lane_tracker,
                        read_length             => $read_length,
                        sample_id        => $infile_info{infile_sample_id},
                        sample_info_href => $sample_info_href,
                    }
                );
            }
            else {
                ## No regexp match i.e. file does not follow filename convention

                $log->warn(
                        q{Could not detect MIP file name convention for file: }
                      . $file_name
                      . q{.} );
                $log->warn(
                    q{Will try to find mandatory information from fastq header.}
                );

                ## Check that file name at least contains sample id
                if ( $file_name !~ /$sample_id/sxm ) {

                    $log->fatal(
q{Please check that the file name contains the sample_id.}
                    );
                    exit 1;
                }

                ## Get run info from fastq file header
                my %fastq_info_header = get_fastq_file_header_info(
                    {
                        file_path => catfile( $infiles_dir, $file_name ),
                        log       => $log,
                        read_file_command => $read_file_command,
                    }
                );

                if ( not exists $fastq_info_header{index} ) {

           # Special case since index is not present in fast headers casaava 1.4
                    $fastq_info_header{index} = $EMPTY_STR;
                }
                ## Adds information derived from infile name to hashes
                $lane_tracker = add_infile_info(
                    {
                        active_parameter_href => $active_parameter_href,
                        ## fastq format does not contain a date of the run,
                        ## so fake it with constant impossible date
                        date           => q{000101},
                        direction      => $fastq_info_header{direction},
                        file_index     => $file_index,
                        file_info_href => $file_info_href,
                        flowcell       => $fastq_info_header{flowcell},
                        index          => $fastq_info_header{index},
                        infile_both_strands_prefix_href =>
                          $infile_both_strands_prefix_href,
                        infile_lane_prefix_href => $infile_lane_prefix_href,
                        is_interleaved          => $is_interleaved,
                        lane                    => $fastq_info_header{lane},
                        lane_tracker            => $lane_tracker,
                        read_length             => $read_length,
                        sample_id               => $sample_id,
                        sample_info_href        => $sample_info_href,
                    }
                );

                $log->info(
                        q{Found following information from fastq header: lane=}
                      . $fastq_info_header{lane}
                      . q{ flow-cell=}
                      . $fastq_info_header{flowcell}
                      . q{ index=}
                      . $fastq_info_header{index}
                      . q{ direction=}
                      . $fastq_info_header{direction},
                );
                $log->warn(
q{Will add fake date '20010101' to follow file convention since this is not recorded in fastq header}
                );
            }
        }
    }
    return;
}

sub parse_fastq_infiles_format {

## Function : Parse infile according to MIP filename convention
## Returns  : %infile_info or undef
## Arguments: $file_name => File name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_name;

    my $tmpl = {
        file_name => {
            defined     => 1,
            required    => 1,
            store       => \$file_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Format::Mip qw{ fastq_file_name_regexp };

    # Store file name features
    my %file_info;

    ## Define MIP fastq file name formats matching regexp
    my %mip_file_name_regexp = fastq_file_name_regexp();

    # Parse file name
    my @file_features = $file_name =~ /$mip_file_name_regexp{regexp}/sxm;

  FEATURE:
    while ( my ( $index, $feature ) =
        each @{ $mip_file_name_regexp{features} } )
    {

        ## Return undef if not all expected features found
        return if ( not $file_features[$index] );

        # Store feature
        $file_info{$feature} = $file_features[$index];
    }
    return %file_info;
}

1;
