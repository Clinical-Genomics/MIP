package MIP::Parse::Gender;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use FindBin qw{ $Bin };
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
use MIP::Constants qw{ $BACKWARD_SLASH
  $CLOSE_PARENTHESIS
  $COLON
  $DASH
  $DOT
  $DOUBLE_QUOTE
  $LOG_NAME
  $OPEN_PARENTHESIS
  $PIPE
  $PLUS
  $SPACE
  $UNDERSCORE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ build_stream_file_cmd get_number_of_male_reads get_sampling_fastq_files parse_fastq_for_gender update_gender_info };
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
    Readonly my $BYTE_STOP_POS  => $BYTE_START_POS * 300;

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

        $cmds_cat[-1] = $cmds_cat[-1] . $CLOSE_PARENTHESIS;

        ## Join for command line
        push @bwa_infiles, join $SPACE, @cmds_cat;
    }
    return @bwa_infiles;
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
    use File::Path qw{ remove_tree };

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

sub get_sampling_fastq_files {

## Function : Get fastq files to sample reads from
## Returns  : $is_interleaved_fastq, @fastq_files
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

sub parse_fastq_for_gender {

## Function : Parse fastq infiles for gender. Update contigs depending on results.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $sample_info_href        => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
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
    use MIP::Program::Bwa qw{ bwa_mem };

    ## All sample ids have a gender - non need to continue
    return if ( not $active_parameter_href->{found_other} );

    ## Unpack
    my $log                = Log::Log4perl->get_logger($LOG_NAME);
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
                count   => 1,
                pattern => $DOUBLE_QUOTE . q{chrY}
                  . $BACKWARD_SLASH
                  . $PIPE . q{Y}
                  . $DOUBLE_QUOTE,
            }
          );

        ## Get the number of male reads by aligning fastq read chunk and counting "chrY" or "Y" aligned reads
        my $y_read_count = get_number_of_male_reads( { commands_ref => \@commands, } );

        ## Update gender info in active_parameter and update contigs depending on results
        update_gender_info(
            {
                active_parameter_href => $active_parameter_href,
                file_info_href        => $file_info_href,
                sample_id             => $sample_id,
                y_read_count          => $y_read_count,
            }
        );
    }
    return;
}

sub update_gender_info {

## Function : Update gender info in active_parameter and update contigs depending on results.
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $file_info_href        => File info hash {REF}
##          : $sample_id             => Sample id
##          : $y_read_count          => Y read count

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
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

    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Constants
    Readonly my $MALE_THRESHOLD => 36;

    if ( $y_read_count > $MALE_THRESHOLD ) {

        $log->info(q{Found male according to fastq reads});

        ## Increment found_male
        $active_parameter_href->{found_male}++;
        $active_parameter_href->{gender_estimation}{$sample_id} = q{male};
        return 1;
    }

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
    return 1;
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

        $cmd = q{<} . $OPEN_PARENTHESIS . q{z} . $cmds_cat_ref->[0];
    }
    else {

        $cmd = q{<} . $OPEN_PARENTHESIS . $cmds_cat_ref->[0];
    }
    return $cmd;
}

1;
