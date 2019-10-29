package MIP::Parse::File;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use FindBin qw{ $Bin };
use File::Spec::Functions qw{ catdir catfile splitdir };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use File::Path qw{ remove_tree };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants
  qw{ $BACKWARD_SLASH $COLON $DASH $DOT $DOUBLE_QUOTE $EMPTY_STR $PIPE $PLUS $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.06;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ parse_fastq_for_gender parse_fastq_infiles parse_fastq_infiles_format parse_file_suffix parse_io_outfiles };

}

sub get_number_of_male_reads {

## Function : Get the number of male reads by aligning fastq read chunk and counting "chrY" or "Y" aligned reads
## Returns  : $y_read_count
## Arguments: $commands_ref => Command array for cat {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $commands_ref;

    my $tmpl = {
        commands_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$commands_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use IPC::Cmd qw(run);

    ## Constants
    Readonly my $MAX_RANDOM_NUMBER => 10_000;

    ## Generate a random integer between 0-10,000.
    my $random_integer = int rand $MAX_RANDOM_NUMBER;

    ## Temporary bash file for commands
    my $bash_temp_file = catfile( $Bin,
        q{estimate_gender_from_reads} . $UNDERSCORE . $random_integer . q{.sh} );

    open my $filehandle, q{>}, $bash_temp_file
      or croak q{Cannot write to}
      . $SPACE
      . $bash_temp_file
      . $COLON
      . $SPACE
      . $OS_ERROR;

    ## Write to file
    say {$filehandle} join $SPACE, @{$commands_ref};

    my $cmds_ref = [ q{bash}, $bash_temp_file ];
    my ( $success, $error_message, $full_buf_ref, $stdout_buf_ref, $stderr_buf_ref ) =
      run( command => $cmds_ref, verbose => 0 );

    my $y_read_count = $stdout_buf_ref->[0];

    ## Clean-up
    close $filehandle;
    remove_tree($bash_temp_file);

    return $y_read_count;
}

sub build_stream_file_cmd {

## Function : Build command for streaming of chunk from fastq file(s)
## Returns  : @bwa_infiles
## Arguments: $fastq_files_ref => Fastq files {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $fastq_files_ref;

    my $tmpl = {
        fastq_files_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$fastq_files_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Gnu::Coreutils qw{ gnu_cat gnu_head gnu_tail };

    ## Constants
    Readonly my $BYTE_START_POS => 10_000;
    Readonly my $BYTE_STOP_POS  => $BYTE_START_POS * 100;

    my @bwa_infiles;

  FILE:
    foreach my $file_path ( @{$fastq_files_ref} ) {

        my @cmds_cat = gnu_cat( { infile_paths_ref => [$file_path], } );

        ## Check gzipped status of file path to choose correct cat binary (cat or zcat). Also prepend stream character.
        $cmds_cat[0] = _get_file_gzipped_status(
            {
                cmds_cat_ref => \@cmds_cat,
                file_path    => $file_path,
            }
        );

        push @cmds_cat, $PIPE;

        ## Start of byte chunk
        push @cmds_cat, gnu_tail( { number => $PLUS . $BYTE_START_POS, } );

        push @cmds_cat, $PIPE;

        ## End of byte chunk
        push @cmds_cat, gnu_head( { number => $BYTE_STOP_POS, } );

        $cmds_cat[-1] = $cmds_cat[-1] . q{)};

        ## Join for command line
        push @bwa_infiles, join $SPACE, @cmds_cat;
    }
    return @bwa_infiles;
}

sub get_sampling_fastq_files {

## Function : Get fastq files to sample reads from
## Returns  :
## Arguments: $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $infile_paths_ref        => Infile paths {REF}
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_lane_prefix_href;
    my $infile_paths_ref;
    my $sample_id;
    my $sample_info_href;

    my $tmpl = {
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_lane_prefix_href,
            strict_type => 1,
        },
        infile_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infile_paths_ref,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
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

    use MIP::Sample_info qw{ get_sequence_run_type get_sequence_run_type_is_interleaved };

    my @fastq_files;

    ## Perform per single-end or read pair
    my $paired_end_tracker = 0;
    my $is_interleaved_fastq;

  INFILE_PREFIX:
    while ( my ( $infile_index, $infile_prefix ) =
        each @{ $infile_lane_prefix_href->{$sample_id} } )
    {

        # Collect paired-end or single-end sequence run type
        my $sequence_run_type = get_sequence_run_type(
            {
                infile_lane_prefix => $infile_prefix,
                sample_id          => $sample_id,
                sample_info_href   => $sample_info_href,
            }
        );

        # Collect interleaved status for fastq file
        $is_interleaved_fastq = get_sequence_run_type_is_interleaved(
            {
                infile_lane_prefix => $infile_prefix,
                sample_id          => $sample_id,
                sample_info_href   => $sample_info_href,
            }
        );

        ## Infile(s)
        push @fastq_files, $infile_paths_ref->[$paired_end_tracker];

        # If second read direction is present
        if ( $sequence_run_type eq q{paired-end} ) {

            # Increment to collect correct read 2
            $paired_end_tracker = $paired_end_tracker + 1;
            push @fastq_files, $infile_paths_ref->[$paired_end_tracker];
        }
        ## Only perform once per sample and fastq file(s)
        last INFILE_PREFIX;
    }
    return $is_interleaved_fastq, @fastq_files;
}

sub update_gender_info {

## Function : Update gender info in active_parameter and update contigs depending on results.
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $file_info_href        => File info hash {REF}
##          : $log                   => Log object
##          : $sample_id             => Sample id
##          : $y_read_count          => Y read count

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $log;
    my $sample_id;
    my $y_read_count;

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
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
        y_read_count => {
            required => 1,
            store    => \$y_read_count,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Update::Contigs qw{ update_contigs_for_run };

    ## Constants
    Readonly my $MALE_THRESHOLD => 25;

    if ( $y_read_count > $MALE_THRESHOLD ) {

        $log->info(q{Found male according to fastq reads});
        ## Increment found_male
        $active_parameter_href->{found_male}++;
        $active_parameter_href->{gender_estimation}{$sample_id} = q{male};
    }
    else {
        $log->info(q{Found female according to fastq reads});

        ## Decrement found_male
        $active_parameter_href->{found_male}--;
        $active_parameter_href->{gender_estimation}{$sample_id} = q{female};

        ## Update contigs depending on settings in run (wes or if only male samples)
        update_contigs_for_run(
            {
                analysis_type_href  => \%{ $active_parameter_href->{analysis_type} },
                exclude_contigs_ref => \@{ $active_parameter_href->{exclude_contigs} },
                file_info_href      => $file_info_href,
                found_male          => $active_parameter_href->{found_male},
                log                 => $log,
            }
        );

    }
    return;
}

sub parse_fastq_for_gender {

## Function : Parse fastq infiles for gender. Update contigs depending on results.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $log                     => Log object
##          : $sample_info_href        => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
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

    use MIP::Gnu::Coreutils qw{ gnu_cut };
    use MIP::Gnu::Software::Gnu_grep qw{ gnu_grep };
    use MIP::Program::Alignment::Bwa qw{ bwa_mem };

    ## All sample ids have a gender - non need to continue
    return if ( not $active_parameter_href->{found_other} );

    ## Unpack
    my $referencefile_path = $active_parameter_href->{human_genome_reference};

  SAMPLE_ID:
    for my $sample_id ( @{ $active_parameter_href->{gender}{others} } ) {

        ## Only for sample with wgs analysis type
        next SAMPLE_ID
          if ( not $active_parameter_href->{analysis_type}{$sample_id} eq q{wgs} );

        $log->warn(qq{Detected gender "other/unknown" for sample_id: $sample_id});
        $log->warn(q{Sampling reads from fastq file to estimate gender});

        ### Estimate gender from reads

        ## Get infile directory
        my $infiles_dir = $file_info_href->{$sample_id}{mip_infiles_dir};

        ## Get fastq files to sample reads from
        my ( $is_interleaved_fastq, @fastq_files ) = get_sampling_fastq_files(
            {
                infile_lane_prefix_href => $infile_lane_prefix_href,
                infile_paths_ref        => $file_info_href->{$sample_id}{mip_infiles},
                sample_id               => $sample_id,
                sample_info_href        => $sample_info_href,
            }
        );

## Add infile dir to infiles
        @fastq_files = map { catfile( $infiles_dir, $_ ) } @fastq_files;

        ## Build command for streaming of chunk from fastq file(s)
        my @bwa_infiles = build_stream_file_cmd( { fastq_files_ref => \@fastq_files, } );

        ## Build bwa mem command
        my @commands = bwa_mem(
            {
                idxbase                 => $referencefile_path,
                infile_path             => $bwa_infiles[0],
                interleaved_fastq_file  => $is_interleaved_fastq,
                mark_split_as_secondary => 1,
                second_infile_path      => $bwa_infiles[1],
                thread_number           => 2,
            }
        );
        push @commands, $PIPE;

        ## Cut column with contig name
        push @commands,
          gnu_cut(
            {
                infile_path => $DASH,
                list        => 3,
            }
          );
        push @commands, $PIPE;

        ## Count contig name for male contig
        push @commands,
          gnu_grep(
            {
                pattern => $DOUBLE_QUOTE . q{chrY}
                  . $BACKWARD_SLASH
                  . $PIPE . q{Y}
                  . $DOUBLE_QUOTE,
                count => 1,
            }
          );

        ## Get the number of male reads by aligning fastq read chunk and counting "chrY" or "Y" aligned reads
        my $y_read_count = get_number_of_male_reads( { commands_ref => \@commands, } );

        ## Update gender info in active_parameter and update contigs depending on results
        update_gender_info(
            {
                active_parameter_href => $active_parameter_href,
                file_info_href        => $file_info_href,
                log                   => $log,
                sample_id             => $sample_id,
                y_read_count          => $y_read_count,
            }
        );

    }
    return;
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

    use MIP::Check::File qw{ check_interleaved };
    use MIP::Check::Parameter qw{ check_infile_contain_sample_id };
    use MIP::Get::File qw{ get_fastq_file_header_info get_read_length };
    use MIP::Set::File qw{ set_file_compression_features };
    use MIP::Sample_info qw{ set_infile_info };

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
            my %infile_info = parse_fastq_infiles_format( { file_name => $file_name, } );

            ## Get sequence read length from file
            $read_length = get_read_length(
                {
                    file_path         => catfile( $infiles_dir, $file_name ),
                    read_file_command => $read_file_command,
                }
            );

            ## Is file interleaved and have proper read direction
            $is_interleaved = check_interleaved(
                {
                    file_path         => catfile( $infiles_dir, $file_name ),
                    log               => $log,
                    read_file_command => $read_file_command,
                }
            );

            ## STAR does not support interleaved fastq files
            if ( ( $active_parameter_href->{analysis_type}{$sample_id} eq q{wts} )
                and $is_interleaved )
            {
                $log->fatal(q{MIP rd_rna does not support interleaved fastq files});
                $log->fatal(
                    q{Please deinterleave: } . catfile( $infiles_dir, $file_name ) );
                exit 1;
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
                        sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                    }
                );

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
                        is_interleaved          => $is_interleaved,
                        lane                    => $infile_info{lane},
                        lane_tracker            => $lane_tracker,
                        read_length             => $read_length,
                        sample_id               => $infile_info{infile_sample_id},
                        sample_info_href        => $sample_info_href,
                    }
                );
            }
            else {
                ## No regexp match i.e. file does not follow filename convention

                $log->warn( q{Could not detect MIP file name convention for file: }
                      . $file_name
                      . q{.} );
                $log->warn(q{Will try to find mandatory information from fastq header.});

                ## Check that file name at least contains sample id
                if ( $file_name !~ /$sample_id/sxm ) {

                    $log->fatal(
                        q{Please check that the file name contains the sample_id.});
                    exit 1;
                }

                ## Get run info from fastq file header
                my %fastq_info_header = get_fastq_file_header_info(
                    {
                        file_path         => catfile( $infiles_dir, $file_name ),
                        log               => $log,
                        read_file_command => $read_file_command,
                    }
                );

                if ( not exists $fastq_info_header{index} ) {

                    # Special case since index is not present in fast headers casaava 1.4
                    $fastq_info_header{index} = $EMPTY_STR;
                }
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
    while ( my ( $index, $feature ) = each @{ $mip_file_name_regexp{features} } ) {

        ## Return undef if not all expected features found
        return if ( not $file_features[$index] );

        # Store feature
        $file_info{$feature} = $file_features[$index];
    }
    return %file_info;
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

sub _get_file_gzipped_status {

## Function : Check gzipped status of file path to choose correct cat binary (cat or zcat). Also prepend stream character.
## Returns  : $cmd
## Arguments: $cmds_cat_ref => Command array for cat {REF}
##          : $file_path    => Fastq File path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path;
    my $cmds_cat_ref;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        cmds_cat_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$cmds_cat_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Parse::File qw{ parse_file_suffix };

    my $cmd;

    ## Parse file suffix in filename.suffix(.gz).
    ## Removes suffix if matching else return undef
    my $is_gzipped = parse_file_suffix(
        {
            file_name   => $file_path,
            file_suffix => $DOT . q{gz},
        }
    );
    if ($is_gzipped) {

        $cmd = q{<(z} . $cmds_cat_ref->[0];
    }
    else {

        $cmd = q{<(} . $cmds_cat_ref->[0];
    }
    return $cmd;
}

1;
