package MIP::Recipes::Analysis::Bwa_mem;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname fileparse };
use File::Spec::Functions qw{ catdir catfile devnull splitpath };
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
use MIP::Constants
  qw{ $ASTERISK $DOT $DOUBLE_QUOTE $EMPTY_STR $LOG_NAME $NEWLINE $PIPE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.18;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_bwa_mem analysis_run_bwa_mem };

}

sub analysis_bwa_mem {

## Function : Performs alignment of single and paired-end as well as interleaved fastq(.gz) files using the bwa mem binary
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_id;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $profile_base_command;
    my $temp_directory;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        case_id => {
            default     => $arg_href->{active_parameter_href}{case_id},
            store       => \$case_id,
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
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        profile_base_command => {
            default     => q{sbatch},
            store       => \$profile_base_command,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Bwa qw{ bwa_mem };
    use MIP::Program::Samtools qw{ samtools_stats samtools_view };
    use MIP::Program::Sambamba qw{ sambamba_sort };
    use MIP::Sample_info
      qw{ get_read_group get_sequence_run_type get_sequence_run_type_is_interleaved set_processing_metafile_in_sample_info set_recipe_metafile_in_sample_info set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Set MIP recipe name
    my $recipe_mode = $active_parameter_href->{$recipe_name};

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $sample_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
        }
    );
    my @infile_paths         = @{ $io{in}{file_paths} };
    my @infile_names         = @{ $io{in}{file_names} };
    my @infile_name_prefixes = @{ $io{in}{file_name_prefixes} };

    my $consensus_analysis_type = $parameter_href->{cache}{consensus_analysis_type};
    my $job_id_chain            = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            attribute      => q{chain},
        }
    );
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my %recipe_resource    = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $job_id_chain,
                id                     => $sample_id,
                file_info_href         => $file_info_href,
                file_name_prefixes_ref => \@{ $infile_lane_prefix_href->{$sample_id} },
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
            }
        )
    );

    my $outdir_path           = $io{out}{dir_path};
    my $outfile_suffix        = $io{out}{file_suffix};
    my @outfile_name_prefixes = @{ $io{out}{file_name_prefixes} };
    my @outfile_paths         = @{ $io{out}{file_paths} };
    my @outfile_path_prefixes = @{ $io{out}{file_path_prefixes} };

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Assign file tags
    my $outfile_tag =
      $file_info_href->{$sample_id}{$recipe_name}{file_tag};

    my $uncompressed_bam_output;
    if ( $outfile_suffix eq q{.bam} ) {

        # Used downstream in samtools view
        $uncompressed_bam_output = 1;
    }

    # Too avoid adjusting infile_index in submitting to jobs
    my $paired_end_tracker = 0;

    ## Perform per single-end or read pair
  INFILE_PREFIX:
    while ( my ( $infile_index, $infile_prefix ) =
        each @{ $infile_lane_prefix_href->{$sample_id} } )
    {

        ## Assign file features
        my $outfile_name_prefix = $outfile_name_prefixes[$infile_index];
        my $outfile_path        = $outfile_paths[$infile_index];
        my $outfile_path_prefix = $outfile_path_prefixes[$infile_index];

        # Collect paired-end or single-end sequence run type
        my $sequence_run_type = get_sequence_run_type(
            {
                infile_lane_prefix => $infile_prefix,
                sample_id          => $sample_id,
                sample_info_href   => $sample_info_href,
            }
        );

        # Collect interleaved status for fastq file
        my $is_interleaved_fastq = get_sequence_run_type_is_interleaved(
            {
                infile_lane_prefix => $infile_prefix,
                sample_id          => $sample_id,
                sample_info_href   => $sample_info_href,
            }
        );

        ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
        my ( $recipe_file_path, $recipe_info_path ) = setup_script(
            {
                active_parameter_href           => $active_parameter_href,
                core_number                     => $recipe_resource{core_number},
                directory_id                    => $sample_id,
                filehandle                      => $filehandle,
                job_id_href                     => $job_id_href,
                memory_allocation               => $recipe_resource{memory},
                log                             => $log,
                recipe_directory                => $recipe_name,
                recipe_name                     => $recipe_name,
                process_time                    => $recipe_resource{time},
                sleep                           => 1,
                source_environment_commands_ref => $recipe_resource{load_env_ref},
                temp_directory                  => $temp_directory,
            }
        );

        # Split to enable submission to %sample_info_qc later
        my ( $volume, $directory, $stderr_file ) =
          splitpath( $recipe_info_path . $DOT . q{stderr.txt} );

        ### SHELL:

        ### BWA MEM
        say {$filehandle} q{## Aligning reads with }
          . $recipe_name
          . q{ and sorting via Sambamba};

        ### Get parameters

        ## Infile(s)
        my $fastq_file_path = $infile_paths[$paired_end_tracker];
        my $second_fastq_file_path;

        # If second read direction is present
        if ( $sequence_run_type eq q{paired-end} ) {

            # Increment to collect correct read 2
            $paired_end_tracker     = $paired_end_tracker + 1;
            $second_fastq_file_path = $infile_paths[$paired_end_tracker];
        }

        ## Read group header line
        my %read_group = get_read_group(
            {
                infile_prefix    => $infile_prefix,
                platform         => $active_parameter_href->{platform},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        my @read_group_headers = (
            $DOUBLE_QUOTE . q{@RG} . q{\t},
            q{ID:} . $read_group{id} . q{\t},
            q{SM:} . $read_group{sm} . q{\t},
            q{PL:} . $read_group{pl} . q{\t},
            q{PU:} . $read_group{pu} . q{\t},
            q{LB:} . $read_group{lb} . $DOUBLE_QUOTE,
        );

        ## Prepare for downstream processing
        # Can be either infile or instream
        my $sambamba_sort_infile;

        # Prior to ALTs in reference genome
        bwa_mem(
            {
                filehandle              => $filehandle,
                idxbase                 => $referencefile_path,
                infile_path             => $fastq_file_path,
                interleaved_fastq_file  => $is_interleaved_fastq,
                mark_split_as_secondary => 1,
                read_group_header       => join( $EMPTY_STR, @read_group_headers ),
                soft_clip_sup_align => $active_parameter_href->{bwa_soft_clip_sup_align},
                second_infile_path  => $second_fastq_file_path,
                thread_number       => $recipe_resource{core_number},
            }
        );

        # Pipe SAM to BAM conversion of aligned reads
        print {$filehandle} $PIPE . $SPACE;

        samtools_view(
            {
                auto_detect_input_format => 1,
                filehandle               => $filehandle,
                infile_path              => q{-},
                thread_number            => $recipe_resource{core_number},
                uncompressed_bam_output  => $uncompressed_bam_output,
                with_header              => 1,
            }
        );
        print {$filehandle} $PIPE . $SPACE;

        ## Set sambamba sort input; Pipe from samtools view
        $sambamba_sort_infile =
          catfile( dirname( devnull() ), q{stdin} );

        ## Increment paired end tracker
        $paired_end_tracker++;

        ## Sort the output from bwa mem|run-bwamem
        sambamba_sort(
            {
                filehandle    => $filehandle,
                infile_path   => $sambamba_sort_infile,
                memory_limit  => $active_parameter_href->{bwa_sambamba_sort_memory_limit},
                outfile_path  => $outfile_path,
                show_progress => 1,
                temp_directory => $temp_directory,
            }
        );
        say {$filehandle} $NEWLINE;

        if ( $active_parameter_href->{bwa_mem_bamstats} ) {

            samtools_stats(
                {
                    auto_detect_input_format => 1,
                    filehandle               => $filehandle,
                    infile_path              => $outfile_path,
                    remove_overlap           => 1,
                }
            );
            print {$filehandle} $PIPE . $SPACE;

            ## Collect raw total sequences and reads mapped from samtools stats and calculate the percentage. Write it to stdout
            _add_percentage_mapped_reads_from_samtools(
                {
                    filehandle   => $filehandle,
                    outfile_path => $outfile_path_prefix . $DOT . q{stats},
                }
            );
            say {$filehandle} $NEWLINE;
        }

        if (    $active_parameter_href->{bwa_mem_cram}
            and $outfile_suffix ne q{.cram} )
        {

            say {$filehandle} q{## Create CRAM file from SAM|BAM};
            samtools_view(
                {
                    filehandle         => $filehandle,
                    infile_path        => $outfile_path,
                    outfile_path       => $outfile_path_prefix . $DOT . q{cram},
                    output_format      => q{cram},
                    referencefile_path => $referencefile_path,
                    with_header        => 1,
                }
            );
            say {$filehandle} $NEWLINE;
        }

        close $filehandle;

        if ( $recipe_mode == 1 ) {

            my $most_complete_format_key =
              q{most_complete} . $UNDERSCORE . substr $outfile_suffix, 1;
            set_processing_metafile_in_sample_info(
                {
                    metafile_tag     => $most_complete_format_key,
                    path             => $outfile_path,
                    sample_id        => $sample_id,
                    sample_info_href => $sample_info_href,
                }
            );

            if (    $active_parameter_href->{bwa_mem_cram}
                and $outfile_suffix ne q{.cram} )
            {

                # Required for analysisRunStatus check downstream
                my $qc_cram_path = $outfile_path_prefix . $DOT . q{cram};
                set_recipe_metafile_in_sample_info(
                    {
                        infile           => $outfile_name_prefix,
                        metafile_tag     => q{cram},
                        path             => $qc_cram_path,
                        recipe_name      => $recipe_name,
                        sample_id        => $sample_id,
                        sample_info_href => $sample_info_href,
                    }
                );

            }
            if ( $active_parameter_href->{bwa_mem_bamstats} ) {

                ## Collect QC metadata info for later use
                my $qc_stats_outfile = $outfile_path_prefix . $DOT . q{stats};
                set_recipe_outfile_in_sample_info(
                    {
                        infile           => $outfile_name_prefix,
                        path             => $qc_stats_outfile,
                        recipe_name      => q{bamstats},
                        sample_id        => $sample_id,
                        sample_info_href => $sample_info_href,
                    }
                );
            }

            set_recipe_outfile_in_sample_info(
                {
                    infile           => $outfile_name_prefix,
                    path             => $outfile_path,
                    recipe_name      => $recipe_name,
                    sample_id        => $sample_id,
                    sample_info_href => $sample_info_href,
                }
            );

            submit_recipe(
                {
                    base_command            => $profile_base_command,
                    case_id                 => $case_id,
                    dependency_method       => q{sample_to_sample_parallel},
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_chain            => $job_id_chain,
                    job_id_href             => $job_id_href,
                    log                     => $log,
                    recipe_file_path        => $recipe_file_path,
                    recipe_files_tracker    => $infile_index,
                    sample_id               => $sample_id,
                    submission_profile => $active_parameter_href->{submission_profile},
                }
            );
        }
    }
    return 1;
}

sub analysis_run_bwa_mem {

## Function : Performs alignment of single and paired-end as well as interleaved fastq(.gz) files using run-bwamem binary
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_id;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $profile_base_command;
    my $temp_directory;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        case_id => {
            default     => $arg_href->{active_parameter_href}{case_id},
            store       => \$case_id,
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
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        profile_base_command => {
            default     => q{sbatch},
            store       => \$profile_base_command,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Bwa qw{ bwa_mem run_bwamem };
    use MIP::Program::Samtools qw{ samtools_stats samtools_view };
    use MIP::Program::Sambamba qw{ sambamba_sort };
    use MIP::Sample_info
      qw{ get_read_group get_sequence_run_type set_processing_metafile_in_sample_info set_recipe_metafile_in_sample_info set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Set MIP recipe name
    my $recipe_mode = $active_parameter_href->{$recipe_name};

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $sample_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
            temp_directory => $temp_directory,
        }
    );
    my @infile_paths         = @{ $io{in}{file_paths} };
    my @infile_names         = @{ $io{in}{file_names} };
    my @infile_name_prefixes = @{ $io{in}{file_name_prefixes} };

    my $consensus_analysis_type = $parameter_href->{cache}{consensus_analysis_type};
    my $job_id_chain            = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            attribute      => q{chain},
        }
    );
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my %recipe_resource    = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $job_id_chain,
                id                     => $sample_id,
                file_info_href         => $file_info_href,
                file_name_prefixes_ref => \@{ $infile_lane_prefix_href->{$sample_id} },
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
                temp_directory         => $temp_directory,
            }
        )
    );

    my $outdir_path           = $io{out}{dir_path};
    my $outfile_suffix        = $io{out}{file_suffix};
    my @outfile_name_prefixes = @{ $io{out}{file_name_prefixes} };
    my @outfile_paths         = @{ $io{out}{file_paths} };
    my @outfile_path_prefixes = @{ $io{out}{file_path_prefixes} };

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Assign file tags
    my $outfile_tag =
      $file_info_href->{$sample_id}{$recipe_name}{file_tag};

    my $uncompressed_bam_output;
    if ( $outfile_suffix eq q{.bam} ) {

        # Used downstream in samtools view
        $uncompressed_bam_output = 1;
    }

    # Too avoid adjusting infile_index in submitting to jobs
    my $paired_end_tracker = 0;

    ## Perform per single-end or read pair
  INFILE_PREFIX:
    while ( my ( $infile_index, $infile_prefix ) =
        each @{ $infile_lane_prefix_href->{$sample_id} } )
    {

        ## Assign file features
        my $outfile_name_prefix = $outfile_name_prefixes[$infile_index];
        my $outfile_path        = $outfile_paths[$infile_index];
        my $outfile_path_prefix = $outfile_path_prefixes[$infile_index];

        # Collect paired-end or single-end sequence run type
        my $sequence_run_type = get_sequence_run_type(
            {
                infile_lane_prefix => $infile_prefix,
                sample_id          => $sample_id,
                sample_info_href   => $sample_info_href,
            }
        );

        ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
        my ( $recipe_file_path, $recipe_info_path ) = setup_script(
            {
                active_parameter_href           => $active_parameter_href,
                core_number                     => $recipe_resource{core_number},
                directory_id                    => $sample_id,
                filehandle                      => $filehandle,
                job_id_href                     => $job_id_href,
                memory_allocation               => $recipe_resource{memory},
                log                             => $log,
                recipe_directory                => $recipe_name,
                recipe_name                     => $recipe_name,
                process_time                    => $recipe_resource{time},
                sleep                           => 1,
                source_environment_commands_ref => $recipe_resource{load_env_ref},
                temp_directory                  => $temp_directory,
            }
        );

        # Split to enable submission to %sample_info_qc later
        my ( $volume, $directory, $stderr_file ) =
          splitpath( $recipe_info_path . $DOT . q{stderr.txt} );

        ### SHELL:

        ### BWA MEM
        say {$filehandle} q{## Aligning reads with }
          . $recipe_name
          . q{ and sorting via Sambamba};

        ### Get parameters

        ## Infile(s)
        my $fastq_file_path = $infile_paths[$paired_end_tracker];
        my $second_fastq_file_path;

        # If second read direction is present
        if ( $sequence_run_type eq q{paired-end} ) {

            # Increment to collect correct read 2
            $paired_end_tracker     = $paired_end_tracker + 1;
            $second_fastq_file_path = $infile_paths[$paired_end_tracker];
        }

        ## Read group header line
        my %read_group = get_read_group(
            {
                infile_prefix    => $infile_prefix,
                platform         => $active_parameter_href->{platform},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        my @read_group_headers = (
            $DOUBLE_QUOTE . q{@RG} . q{\t},
            q{ID:} . $read_group{id} . q{\t},
            q{SM:} . $read_group{sm} . q{\t},
            q{PL:} . $read_group{pl} . q{\t},
            q{PU:} . $read_group{pu} . q{\t},
            q{LB:} . $read_group{lb} . $DOUBLE_QUOTE,
        );

        ## Prepare for downstream processing
        # Can be either infile or instream
        my $sambamba_sort_infile;

        # If post to ALTs in reference genome
        run_bwamem(
            {
                filehandle           => $filehandle,
                hla_typing           => $active_parameter_href->{bwa_mem_hla},
                infile_path          => $fastq_file_path,
                idxbase              => $referencefile_path,
                outfiles_prefix_path => $outfile_path_prefix,
                read_group_header    => join( $EMPTY_STR, @read_group_headers ),
                second_infile_path   => $second_fastq_file_path,
                thread_number        => $recipe_resource{core_number},
            }
        );
        print {$filehandle} $PIPE . $SPACE;
        print {$filehandle} q{sh} . $SPACE;
        say   {$filehandle} $NEWLINE;

        ## Set sambamba sort input; Sort directly from run-bwakit
        $sambamba_sort_infile = $outfile_path_prefix . $DOT . q{aln} . $outfile_suffix;

        ## Increment paired end tracker
        $paired_end_tracker++;

        ## Sort the output from bwa mem|run-bwamem
        sambamba_sort(
            {
                filehandle    => $filehandle,
                infile_path   => $sambamba_sort_infile,
                memory_limit  => $active_parameter_href->{bwa_sambamba_sort_memory_limit},
                outfile_path  => $outfile_path,
                show_progress => 1,
                temp_directory => $temp_directory,
            }
        );
        say {$filehandle} $NEWLINE;

        if ( $active_parameter_href->{bwa_mem_bamstats} ) {

            samtools_stats(
                {
                    auto_detect_input_format => 1,
                    filehandle               => $filehandle,
                    infile_path              => $outfile_path,
                    remove_overlap           => 1,
                }
            );
            print {$filehandle} $PIPE . $SPACE;

            ## Collect raw total sequences and reads mapped from samtools stats and calculate the percentage. Write it to stdout
            _add_percentage_mapped_reads_from_samtools(
                {
                    filehandle   => $filehandle,
                    outfile_path => $outfile_path_prefix . $DOT . q{stats},
                }
            );
            say {$filehandle} $NEWLINE;
        }

        if (    $active_parameter_href->{bwa_mem_cram}
            and $outfile_suffix ne q{.cram} )
        {

            say {$filehandle} q{## Create CRAM file from SAM|BAM};
            samtools_view(
                {
                    filehandle         => $filehandle,
                    infile_path        => $outfile_path,
                    outfile_path       => $outfile_path_prefix . $DOT . q{cram},
                    output_format      => q{cram},
                    referencefile_path => $referencefile_path,
                    with_header        => 1,
                }
            );
            say {$filehandle} $NEWLINE;
        }

        close $filehandle;

        if ( $recipe_mode == 1 ) {

            my $most_complete_format_key =
              q{most_complete} . $UNDERSCORE . substr $outfile_suffix, 1;
            set_processing_metafile_in_sample_info(
                {
                    metafile_tag     => $most_complete_format_key,
                    path             => $outfile_path,
                    sample_id        => $sample_id,
                    sample_info_href => $sample_info_href,
                }
            );

            if (    $active_parameter_href->{bwa_mem_cram}
                and $outfile_suffix ne q{.cram} )
            {

                # Required for analysisRunStatus check downstream
                my $qc_cram_path = $outfile_path_prefix . $DOT . q{cram};
                set_recipe_metafile_in_sample_info(
                    {
                        infile           => $outfile_name_prefix,
                        metafile_tag     => q{cram},
                        path             => $qc_cram_path,
                        recipe_name      => $recipe_name,
                        sample_id        => $sample_id,
                        sample_info_href => $sample_info_href,
                    }
                );

            }
            if ( $active_parameter_href->{bwa_mem_bamstats} ) {

                ## Collect QC metadata info for later use
                my $qc_stats_outfile = $outfile_path_prefix . $DOT . q{stats};
                set_recipe_outfile_in_sample_info(
                    {
                        infile           => $outfile_name_prefix,
                        path             => $qc_stats_outfile,
                        recipe_name      => q{bamstats},
                        sample_id        => $sample_id,
                        sample_info_href => $sample_info_href,
                    }
                );
            }

            set_recipe_outfile_in_sample_info(
                {
                    infile           => $outfile_name_prefix,
                    path             => $outfile_path,
                    recipe_name      => $recipe_name,
                    sample_id        => $sample_id,
                    sample_info_href => $sample_info_href,
                }
            );

            submit_recipe(
                {
                    base_command            => $profile_base_command,
                    case_id                 => $case_id,
                    dependency_method       => q{sample_to_sample_parallel},
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_chain            => $job_id_chain,
                    job_id_href             => $job_id_href,
                    log                     => $log,
                    recipe_file_path        => $recipe_file_path,
                    recipe_files_tracker    => $infile_index,
                    sample_id               => $sample_id,
                    submission_profile => $active_parameter_href->{submission_profile},
                }
            );
        }
    }
    return 1;
}

sub _add_percentage_mapped_reads_from_samtools {

## Function : Collect raw total sequences and reads mapped from samtools stats and calculate the percentage. Write it to stdout.
## Returns  :
## Arguments: $filehandle   => Filehandle to write to
##          : $outfile_path => Outfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $outfile_path;

    my $tmpl = {
        filehandle   => { defined => 1, required => 1, store => \$filehandle, },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Add percentage mapped reads to samtools stats output
    # Execute perl
    print {$filehandle} q?perl -ne '?;

    # Initiate variables
    print {$filehandle} q?$raw; $map; ?;

    # Remove newline
    print {$filehandle} q?chomp $_; ?;

    # Always relay incoming line to stdout
    print {$filehandle} q?print $_, qq{\n}; ?;

    # Find raw total sequences
    print {$filehandle} q?if ($_=~/raw total sequences:\s+(\d+)/) { ?;

    # Assign raw total sequence
    print {$filehandle} q?$raw = $1; ?;

    # End if
    print {$filehandle} q?} ?;

    # Find reads mapped
    print {$filehandle} q?elsif ($_=~/reads mapped:\s+(\d+)/) { ?;

    # Assign reads mapped
    print {$filehandle} q?$map = $1; ?;

    # Calculate percentage
    print {$filehandle} q?my $percentage = ($map / $raw ) * 100; ?;

    # Write calculation to stdout
    print {$filehandle} q?print qq{percentage mapped reads:\t} . $percentage . qq{\n}?;

    # End elsif
    print {$filehandle} q?} ?;

    # End oneliner
    print {$filehandle} q?' ?;
    say   {$filehandle} q{>} . $SPACE . $outfile_path;

    return;
}

1;
