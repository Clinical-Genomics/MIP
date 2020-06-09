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
    our $VERSION = 1.09;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ build_stream_file_cmd
      get_number_of_male_reads
      get_sampling_fastq_files
      parse_fastq_for_gender
      update_gender_info
    };
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

    use MIP::Program::Gnu::Coreutils qw{ gnu_head gnu_tail };

    ## Constants
    Readonly my $BYTE_START_POS => 10_000;
    Readonly my $BYTE_STOP_POS  => $BYTE_START_POS * 300;

    my @bwa_infiles;

  FILE:
    foreach my $file_path ( @{$fastq_files_ref} ) {

        ## Check gzipped status of file path to choose correct cat binary (cat or gzip). Also prepend stream character.
        my @cmds_cat = _get_file_read_commands(
            {
                file_path => $file_path,
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

    use MIP::Environment::Child_process qw{ child_process };
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

    my $cmds_ref       = [ q{bash}, $bash_temp_file ];
    my %process_return = child_process(
        {
            commands_ref => $cmds_ref,
            process_type => q{ipc_cmd_run},
        }
    );

    my $y_read_count = $process_return{stdouts_ref}->[0];

    ## Clean-up
    close $filehandle;
    remove_tree($bash_temp_file);

    return $y_read_count;
}

sub get_sampling_fastq_files {

## Function : Get fastq files to sample reads from
## Returns  : $is_interleaved_fastq, @fastq_files
## Arguments: $file_info_sample_href => File info sample hash
##          : $infile_paths_ref      => Infile paths {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_sample_href;
    my $infile_paths_ref;

    my $tmpl = {
        file_info_sample_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_sample_href,
            strict_type => 1,
        },
        infile_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infile_paths_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @fastq_files;

    ## Perform per single-end or read pair
    my $paired_end_tracker = 0;

  INFILE_PREFIX:
    foreach
      my $infile_prefix ( @{ $file_info_sample_href->{no_direction_infile_prefixes} } )
    {

        push @fastq_files, $infile_paths_ref->[$paired_end_tracker];

        my $sequence_run_type =
          $file_info_sample_href->{$infile_prefix}{sequence_run_type};

        # If second read direction is present
        if ( $sequence_run_type eq q{paired-end} ) {

            # Increment to collect correct read 2
            $paired_end_tracker = $paired_end_tracker + 1;
            push @fastq_files, $infile_paths_ref->[$paired_end_tracker];
        }

        my $is_interleaved_fastq = $sequence_run_type eq q{interleaved} ? 1 : 0;

        ## Only perform once per sample and fastq file(s)
        return $is_interleaved_fastq, @fastq_files;
    }
    return;
}

sub parse_fastq_for_gender {

## Function : Parse fastq infiles for gender. Update contigs depending on results.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $consensus_analysis_type => Consensus analysis_type
##          : $file_info_href          => File info hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $consensus_analysis_type;
    my $file_info_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        consensus_analysis_type => {
            defined     => 1,
            required    => 1,
            store       => \$consensus_analysis_type,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Active_parameter qw{ get_active_parameter_attribute };
    use MIP::File_info qw{ get_sample_file_attribute };
    use MIP::Program::Gnu::Coreutils qw{ gnu_cut };
    use MIP::Program::Gnu::Software::Gnu_grep qw{ gnu_grep };
    use MIP::Program::Bwa qw{ bwa_mem };

    my $other_gender_count = get_active_parameter_attribute(
        {
            active_parameter_href => $active_parameter_href,
            attribute             => q{others},
            parameter_name        => q{gender},
        }
    );
    ## All sample ids have a gender - no need to continue
    return if ( not $other_gender_count );

    ## Unpack
    my $log                = Log::Log4perl->get_logger($LOG_NAME);
    my $referencefile_path = $active_parameter_href->{human_genome_reference};

  SAMPLE_ID:
    for my $sample_id ( @{ $active_parameter_href->{gender}{others} } ) {

        next SAMPLE_ID
          if ( not $active_parameter_href->{analysis_type}{$sample_id} eq q{wgs} );

        $log->warn(qq{Detected gender "other/unknown" for sample_id: $sample_id});
        $log->warn(q{Sampling reads from fastq file to estimate gender});

        ### Estimate gender from reads

        my %file_info_sample = get_sample_file_attribute(
            {
                file_info_href => $file_info_href,
                sample_id      => $sample_id,
            }
        );

        my $infiles_dir = $file_info_sample{mip_infiles_dir};

        my ( $is_interleaved_fastq, @fastq_files ) = get_sampling_fastq_files(
            {
                file_info_sample_href => \%file_info_sample,
                infile_paths_ref      => $file_info_sample{mip_infiles},
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
                active_parameter_href   => $active_parameter_href,
                consensus_analysis_type => $consensus_analysis_type,
                file_info_href          => $file_info_href,
                sample_id               => $sample_id,
                y_read_count            => $y_read_count,
            }
        );
    }
    return 1;
}

sub update_gender_info {

## Function : Update gender info in active_parameter and update contigs depending on results.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $consensus_analysis_type => Consensus analysis_type
##          : $file_info_href          => File info hash {REF}
##          : $sample_id               => Sample id
##          : $y_read_count            => Y read count

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $consensus_analysis_type;
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
        consensus_analysis_type => {
            defined     => 1,
            required    => 1,
            store       => \$consensus_analysis_type,
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

    use MIP::Active_parameter
      qw{ add_gender remove_sample_id_from_gender set_gender_estimation set_include_y };
    use MIP::Contigs qw{ update_contigs_for_run };

    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Constants
    Readonly my $MALE_THRESHOLD => 36;

    my $gender  = $y_read_count > $MALE_THRESHOLD ? q{male} : q{female};
    my $genders = $gender . q{s};
    $log->info(qq{Found $gender according to fastq reads});

    add_gender(
        {
            active_parameter_href => $active_parameter_href,
            sample_id             => $sample_id,
            gender                => $genders,
        }
    );

    ## For tracability
    set_gender_estimation(
        {
            active_parameter_href => $active_parameter_href,
            gender                => $gender,
            sample_id             => $sample_id,
        }
    );

    remove_sample_id_from_gender(
        {
            active_parameter_href => $active_parameter_href,
            gender                => q{others},
            sample_id             => $sample_id,
        }
    );

    set_include_y(
        {
            active_parameter_href => $active_parameter_href,
        }
    );

    ## Update contigs depending on settings in run (wes or if only male samples)
    update_contigs_for_run(
        {
            consensus_analysis_type => $consensus_analysis_type,
            exclude_contigs_ref     => \@{ $active_parameter_href->{exclude_contigs} },
            file_info_href          => $file_info_href,
            include_y               => $active_parameter_href->{include_y},
        }
    );
    return 1;
}

sub _get_file_read_commands {

## Function : Check gzipped status of file path to choose correct cat binary (cat or gzip). Also prepend stream character.
## Returns  : @read_cmds
## Arguments: $file_path => Fastq File path to check status for

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Path qw{ check_gzipped };
    use MIP::Program::Gnu::Coreutils qw{ gnu_cat };
    use MIP::Program::Gzip qw{ gzip };

    my @read_cmds;

    my $is_gzipped = check_gzipped(
        {
            file_name => $file_path,
        }
    );
    if ($is_gzipped) {

        @read_cmds = gzip(
            {
                decompress       => 1,
                infile_paths_ref => [$file_path],
                stdout           => 1,
            }
        );
    }
    else {

        @read_cmds = gnu_cat( { infile_paths_ref => [$file_path], } );
    }

    $read_cmds[0] = q{<} . $OPEN_PARENTHESIS . $read_cmds[0];

    return @read_cmds;
}

1;
