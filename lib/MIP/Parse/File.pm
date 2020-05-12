package MIP::Parse::File;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $DOT $EMPTY_STR $LOG_NAME $SPACE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.13;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ parse_fastq_infiles parse_file_suffix parse_io_outfiles };

}

sub parse_fastq_infiles {

## Function : Parse fastq infiles for MIP processing.
## Returns  :
## Arguments: $active_parameter_href           => Active parameters for this analysis hash {REF}
##          : $file_info_href                  => File info hash {REF}
##          : $infile_both_strands_prefix_href => The infile(s) without the ".ending" and strand info {REF}
##          : $infile_lane_prefix_href         => Infile(s) without the ".ending" {REF}
##          : $log                             => Log object
##          : $sample_info_href                => Info on samples and case hash {REF}

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

    use MIP::Fastq qw{ get_fastq_file_header_info };
    use MIP::File_info qw{
      get_sample_file_attribute
      parse_files_compression_status
      parse_sample_fastq_file_attributes
    };
    use MIP::Sample_info qw{ set_infile_info };

  SAMPLE_ID:
    for my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        # Needed to be able to track when lanes are finished
        my $lane_tracker = 0;
        my %fastq_info_header;

        my %file_info_sample = get_sample_file_attribute(
            {
                file_info_href => $file_info_href,
                sample_id      => $sample_id,
            }
        );
        ## Unpack
        my $infiles_dir = $file_info_sample{mip_infiles_dir};

      INFILE:
        while ( my ( $file_index, $file_name ) =
            each @{ $file_info_sample{mip_infiles} } )
        {

            my %infile_info = parse_sample_fastq_file_attributes(
                {
                    file_info_href => $file_info_href,
                    file_name      => $file_name,
                    infiles_dir    => $infiles_dir,
                    sample_id      => $sample_id,
                }
            );

            ## If filename convention is followed
            if ( not exists $infile_info{date} ) {

                ## No regexp match i.e. file does not follow filename convention
                $log->warn(q{Will try to find mandatory information from fastq header});

                ## Get run info from fastq file header
                %fastq_info_header = get_fastq_file_header_info(
                    {
                        file_path         => catfile( $infiles_dir, $file_name ),
                        read_file_command => $infile_info{read_file_command},
                    }
                );

                if ( not exists $fastq_info_header{index} ) {

                    # Special case since index is not present in fast headers casaava 1.4
                    $fastq_info_header{index} = $EMPTY_STR;
                }
            }
            if ( exists $infile_info{date} ) {

                ## Adds information derived from infile name to hashes
                $lane_tracker = set_infile_info(
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
                        is_interleaved          => $infile_info{is_interleaved},
                        lane                    => $infile_info{lane},
                        lane_tracker            => $lane_tracker,
                        read_length             => $infile_info{read_length},
                        sample_id               => $infile_info{infile_sample_id},
                        sample_info_href        => $sample_info_href,
                    }
                );
            }
            else {
                ## Adds information derived from infile name to hashes
                $lane_tracker = set_infile_info(
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
                        is_interleaved          => $infile_info{is_interleaved},
                        lane                    => $fastq_info_header{lane},
                        lane_tracker            => $lane_tracker,
                        read_length             => $infile_info{read_length},
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

            parse_files_compression_status(
                {
                    file_info_href => $file_info_href,
                    sample_id      => $sample_id,
                }
            );
        }
    }
    return;
}

sub parse_file_suffix {

## Function : Parse file suffix in filename.suffix(.gz). Removes suffix if matching else return undef
## Returns  : undef | $file_name
## Arguments: $file_name   => File name
##          : $file_suffix => File suffix to be removed

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_name;
    my $file_suffix;

    my $tmpl = {
        file_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_name
        },
        file_suffix => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_suffix
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my ($file_name_nosuffix) =
      $file_name =~ / (\S+)($file_suffix$ | $file_suffix.gz$) /xsm;

    return $file_name_nosuffix;
}

sub parse_io_outfiles {

## Function : Set and get the io files per chain, id and stream
## Returns  : %io
## Arguments: $chain_id               => Chain of recipe
##          : $file_info_href         => File info hash {REF}
##          : $file_name_prefixes     => Build outfile using file name prefix
##          : $file_name_prefixes_ref => Build outfile using file name prefixes {REF}
##          : $file_paths_ref         => File paths {REF}
##          : $id                     => Id (sample or case)
##          : $iterators_ref          => Build outfile using iterator (e.g contigs) {REF}
##          : $outdata_dir            => Outdata directory
##          : $parameter_href         => Parameter hash {REF}
##          : $recipe_name            => Recipe name
##          : $stream                 => Stream (out)
##          : $temp_directory         => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $chain_id;
    my $file_info_href;
    my $file_name_prefix;
    my $file_name_prefixes_ref;
    my $file_paths_ref;
    my $id;
    my $iterators_ref;
    my $outdata_dir;
    my $parameter_href;
    my $recipe_name;
    my $temp_directory;

    ## Default(s)
    my $stream;

    my $tmpl = {
        chain_id => {
            defined     => 1,
            required    => 1,
            store       => \$chain_id,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        file_name_prefix => {
            store       => \$file_name_prefix,
            strict_type => 1,
        },
        file_name_prefixes_ref => {
            default     => [],
            store       => \$file_name_prefixes_ref,
            strict_type => 1,
        },
        file_paths_ref => {
            default     => [],
            store       => \$file_paths_ref,
            strict_type => 1,
        },
        id => {
            defined     => 1,
            required    => 1,
            store       => \$id,
            strict_type => 1,
        },
        iterators_ref => {
            default     => [],
            store       => \$iterators_ref,
            strict_type => 1,
        },
        outdata_dir => {
            store       => \$outdata_dir,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        stream => {
            allow       => [qw{ out }],
            default     => q{out},
            store       => \$stream,
            strict_type => 1,
        },
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes };
    use MIP::Set::File qw{ set_io_files };

    my @file_paths = @{$file_paths_ref};

    ## Build default @file_paths
    if ( not @file_paths and $outdata_dir ) {

        my %rec_atr = get_recipe_attributes(
            {
                parameter_href => $parameter_href,
                recipe_name    => $recipe_name,
            }
        );
        my $outfile_tag    = $rec_atr{file_tag}       //= $EMPTY_STR;
        my $outfile_suffix = $rec_atr{outfile_suffix} //= $EMPTY_STR;
        my $directory = catdir( $outdata_dir, $id, $recipe_name );

        ## Default paths with iterators
        if ( @{$iterators_ref} and $file_name_prefix ) {

            ## Localize as we will mutate elements
            my @iterators = @{$iterators_ref};
            foreach my $iterator (@iterators) {

                ## Add "." if not empty string
                if ( $iterator ne $EMPTY_STR ) {

                    $iterator = $DOT . $iterator;
                }
            }
            @file_paths =
              map {
                catfile( $directory,
                    $file_name_prefix . $outfile_tag . $_ . $outfile_suffix )
              } @iterators;
        }
        ## Default paths without iterators
        else {

            ## $file_name_prefixes_needs to be set
            croak q{Missing argument!} if not @{$file_name_prefixes_ref};
            @file_paths =
              map { catfile( $directory, $_ . $outfile_tag . $outfile_suffix ) }
              @{$file_name_prefixes_ref};
        }
    }

    ## Set the io files per chain and stream
    set_io_files(
        {
            chain_id       => $chain_id,
            id             => $id,
            file_paths_ref => \@file_paths,
            file_info_href => $file_info_href,
            recipe_name    => $recipe_name,
            stream         => $stream,
            temp_directory => $temp_directory,
        }
    );

    my %io = get_io_files(
        {
            id             => $id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => $stream,
        }
    );

    return %io;
}

1;
