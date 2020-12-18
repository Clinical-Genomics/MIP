package MIP::Recipes::Analysis::Bwa_mem;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname fileparse };
use File::Spec::Functions qw{ catdir catfile devnull splitpath };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

## MIPs lib/
use MIP::Constants
  qw{ $AMPERSAND $ASTERISK $DOT $DOUBLE_QUOTE $LOG_NAME $NEWLINE $PIPE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_bwa_mem analysis_bwa_mem2 analysis_run_bwa_mem };

}

sub analysis_bwa_mem {

    ## Function : Performs alignment of single and paired-end as well as interleaved fastq(.gz) files using the bwa mem binary
    ## Returns  :
    ## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
    ##          : $case_id                 => Family id
    ##          : $file_info_href          => File info hash {REF}
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

    use MIP::File_info qw{get_sample_file_attribute};
    use MIP::Get::File qw{get_io_files};
    use MIP::Get::Parameter qw{get_recipe_attributes get_recipe_resources};
    use MIP::Parse::File qw{parse_io_outfiles};
    use MIP::Processmanagement::Processes qw{submit_recipe};
    use MIP::Program::Bwa qw{bwa_mem};
    use MIP::Program::Samtools qw{ samtools_index samtools_stats samtools_sort samtools_view};
    use MIP::Sample_info qw{
      get_rg_header_line
      set_recipe_metafile_in_sample_info
      set_recipe_outfile_in_sample_info};
    use MIP::Script::Setup_script qw{setup_script};

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
    my @infile_paths = @{ $io{in}{file_paths} };

    my $job_id_chain = get_recipe_attributes(
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

    my %file_info_sample = get_sample_file_attribute(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
    );

    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $job_id_chain,
                id                     => $sample_id,
                file_info_href         => $file_info_href,
                file_name_prefixes_ref => $file_info_sample{no_direction_infile_prefixes},
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
            }
        )
    );

    my $outfile_suffix        = $io{out}{file_suffix};
    my @outfile_name_prefixes = @{ $io{out}{file_name_prefixes} };
    my @outfile_paths         = @{ $io{out}{file_paths} };
    my @outfile_path_prefixes = @{ $io{out}{file_path_prefixes} };

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    my $output_format;
    my $uncompressed_bam_output;
    if ( $outfile_suffix eq q{.bam} ) {

        # Used downstream in samtools view
        $uncompressed_bam_output = 1;
        $output_format           = q{bam};
    }

    # Too avoid adjusting infile_index in submitting to jobs
    my $paired_end_tracker = 0;

    ## Perform per single-end or read pair
  INFILE_PREFIX:
    while ( my ( $infile_index, $infile_prefix ) =
        each @{ $file_info_sample{no_direction_infile_prefixes} } )
    {

        ## Assign file features
        my $outfile_name_prefix = $outfile_name_prefixes[$infile_index];
        my $outfile_path        = $outfile_paths[$infile_index];
        my $outfile_path_prefix = $outfile_path_prefixes[$infile_index];

        # Collect interleaved status for fastq file
        my $sequence_run_type = get_sample_file_attribute(
            {
                attribute      => q{sequence_run_type},
                file_info_href => $file_info_href,
                file_name      => $infile_prefix,
                sample_id      => $sample_id,
            }
        );
        my $is_interleaved_fastq = $sequence_run_type eq q{interleaved} ? 1 : 0;

        ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
        my ( $recipe_file_path, $recipe_info_path ) = setup_script(
            {
                active_parameter_href => $active_parameter_href,
                core_number           => $recipe_resource{core_number},
                directory_id          => $sample_id,
                filehandle            => $filehandle,
                job_id_href           => $job_id_href,
                memory_allocation     => $recipe_resource{memory},
                recipe_directory      => $recipe_name,
                recipe_name           => $recipe_name,
                process_time          => $recipe_resource{time},
                temp_directory        => $temp_directory,
            }
        );

        ### SHELL:

        ### BWA MEM
        say {$filehandle} q{## Aligning reads with } . $recipe_name . q{ and sorting via Samtools};

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

        ## Construct RG header line
        my $rg_header_line = get_rg_header_line(
            {
                infile_prefix    => $infile_prefix,
                platform         => $active_parameter_href->{platform},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
                separator        => q{\t},
            }
        );
        ## Add missing "@RG"
        $rg_header_line = $DOUBLE_QUOTE . q{@RG} . q{\t} . $rg_header_line . $DOUBLE_QUOTE;

        # Prior to ALTs in reference genome
        bwa_mem(
            {
                filehandle              => $filehandle,
                idxbase                 => $referencefile_path,
                infile_path             => $fastq_file_path,
                interleaved_fastq_file  => $is_interleaved_fastq,
                mark_split_as_secondary => 1,
                read_group_header       => $rg_header_line,
                soft_clip_sup_align     => $active_parameter_href->{bwa_soft_clip_sup_align},
                second_infile_path      => $second_fastq_file_path,
                thread_number           => $recipe_resource{core_number},
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

        ## Set samtools sort input; Pipe from samtools view
        my $samtools_sort_infile =
          catfile( dirname( devnull() ), q{stdin} );

        ## Increment paired end tracker
        $paired_end_tracker++;

        ## Sort the output from bwa mem
        samtools_sort(
            {
                filehandle            => $filehandle,
                infile_path           => $samtools_sort_infile,
                max_memory_per_thread => 2 . q{G},
                outfile_path          => $outfile_path,
                output_format         => $output_format,
                temp_file_path_prefix => catfile( $temp_directory, q{samtools_sort_temp} ),
                thread_number         => $recipe_resource{core_number},
            }
        );
        say {$filehandle} $NEWLINE;

        samtools_index(
            {
                filehandle  => $filehandle,
                infile_path => $outfile_path,
            }
        );
        say {$filehandle} $AMPERSAND . $NEWLINE;

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
            say {$filehandle} $AMPERSAND . $NEWLINE;
        }
        say {$filehandle} q{wait} . $NEWLINE;

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
                    base_command         => $profile_base_command,
                    case_id              => $case_id,
                    dependency_method    => q{sample_to_sample_parallel},
                    job_id_chain         => $job_id_chain,
                    job_id_href          => $job_id_href,
                    job_reservation_name => $active_parameter_href->{job_reservation_name},
                    log                  => $log,
                    max_parallel_processes_count_href =>
                      $file_info_href->{max_parallel_processes_count},
                    recipe_file_path     => $recipe_file_path,
                    recipe_files_tracker => $infile_index,
                    sample_id            => $sample_id,
                    submission_profile   => $active_parameter_href->{submission_profile},
                }
            );
        }
    }
    return 1;
}

sub analysis_bwa_mem2 {

## Function : Performs alignment of single and paired-end as well as interleaved fastq(.gz) files using the bwa mem 2 binary
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File info hash {REF}
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

    use MIP::File_info qw{ get_sample_file_attribute };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Bwa qw{ bwa_mem2_mem };
    use MIP::Program::Samtools qw{ samtools_index samtools_stats samtools_sort samtools_view };
    use MIP::Sample_info qw{
      get_rg_header_line
      set_recipe_metafile_in_sample_info
      set_recipe_outfile_in_sample_info };
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
    my @infile_paths = @{ $io{in}{file_paths} };

    my $job_id_chain = get_recipe_attributes(
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

    my %file_info_sample = get_sample_file_attribute(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
    );

    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $job_id_chain,
                id                     => $sample_id,
                file_info_href         => $file_info_href,
                file_name_prefixes_ref => $file_info_sample{no_direction_infile_prefixes},
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
            }
        )
    );

    my $outfile_suffix        = $io{out}{file_suffix};
    my @outfile_name_prefixes = @{ $io{out}{file_name_prefixes} };
    my @outfile_paths         = @{ $io{out}{file_paths} };
    my @outfile_path_prefixes = @{ $io{out}{file_path_prefixes} };

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    my $output_format;
    my $uncompressed_bam_output;
    if ( $outfile_suffix eq q{.bam} ) {

        # Used downstream in samtools view
        $uncompressed_bam_output = 1;
        $output_format           = q{bam};
    }

    # Too avoid adjusting infile_index in submitting to jobs
    my $paired_end_tracker = 0;

    ## Perform per single-end or read pair
  INFILE_PREFIX:
    while ( my ( $infile_index, $infile_prefix ) =
        each @{ $file_info_sample{no_direction_infile_prefixes} } )
    {

        ## Assign file features
        my $outfile_name_prefix = $outfile_name_prefixes[$infile_index];
        my $outfile_path        = $outfile_paths[$infile_index];
        my $outfile_path_prefix = $outfile_path_prefixes[$infile_index];

        # Collect interleaved status for fastq file
        my $sequence_run_type = get_sample_file_attribute(
            {
                attribute      => q{sequence_run_type},
                file_info_href => $file_info_href,
                file_name      => $infile_prefix,
                sample_id      => $sample_id,
            }
        );
        my $is_interleaved_fastq = $sequence_run_type eq q{interleaved} ? 1 : 0;

        ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
        my ( $recipe_file_path, $recipe_info_path ) = setup_script(
            {
                active_parameter_href => $active_parameter_href,
                core_number           => $recipe_resource{core_number},
                directory_id          => $sample_id,
                filehandle            => $filehandle,
                job_id_href           => $job_id_href,
                memory_allocation     => $recipe_resource{memory},
                recipe_directory      => $recipe_name,
                recipe_name           => $recipe_name,
                process_time          => $recipe_resource{time},
                temp_directory        => $temp_directory,
            }
        );

        ### SHELL:

        ### BWA MEM
        say {$filehandle} q{## Aligning reads with } . $recipe_name . q{ and sorting via Samtools};

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

        ## Construct RG header line
        my $rg_header_line = get_rg_header_line(
            {
                infile_prefix    => $infile_prefix,
                platform         => $active_parameter_href->{platform},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
                separator        => q{\t},
            }
        );
        ## Add missing "@RG"
        $rg_header_line = $DOUBLE_QUOTE . q{@RG} . q{\t} . $rg_header_line . $DOUBLE_QUOTE;
        my $samtools_view_outfile_path =
          $outfile_path_prefix . $UNDERSCORE . q{mem} . $outfile_suffix;

        # Prior to ALTs in reference genome
        bwa_mem2_mem(
            {
                filehandle              => $filehandle,
                idxbase                 => $referencefile_path,
                infile_path             => $fastq_file_path,
                interleaved_fastq_file  => $is_interleaved_fastq,
                mark_split_as_secondary => 1,
                read_group_header       => $rg_header_line,
                soft_clip_sup_align     => $active_parameter_href->{bwa_soft_clip_sup_align},
                second_infile_path      => $second_fastq_file_path,
                thread_number           => $recipe_resource{core_number},
            }
        );

        # Pipe SAM to BAM conversion of aligned reads
        print {$filehandle} $PIPE . $SPACE;

        samtools_view(
            {
                auto_detect_input_format => 1,
                filehandle               => $filehandle,
                infile_path              => q{-},
                outfile_path             => $samtools_view_outfile_path,
                thread_number            => $recipe_resource{core_number},
                uncompressed_bam_output  => $uncompressed_bam_output,
                with_header              => 1,
            }
        );
        say {$filehandle} $NEWLINE;

        ## Set samtools sort input;
        my $samtools_sort_infile = $samtools_view_outfile_path;

        ## Increment paired end tracker
        $paired_end_tracker++;

        ## Sort the output from bwa mem
        samtools_sort(
            {
                filehandle            => $filehandle,
                infile_path           => $samtools_sort_infile,
                max_memory_per_thread => 2 . q{G},
                outfile_path          => $outfile_path,
                output_format         => $output_format,
                temp_file_path_prefix => catfile( $temp_directory, q{samtools_sort_temp} ),
                thread_number         => $recipe_resource{core_number},
            }
        );
        say {$filehandle} $NEWLINE;

        samtools_index(
            {
                filehandle  => $filehandle,
                infile_path => $outfile_path,
            }
        );
        say {$filehandle} $AMPERSAND . $NEWLINE;

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
            say {$filehandle} $AMPERSAND . $NEWLINE;
        }
        say {$filehandle} q{wait} . $NEWLINE;

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
                    base_command         => $profile_base_command,
                    case_id              => $case_id,
                    dependency_method    => q{sample_to_sample_parallel},
                    job_id_chain         => $job_id_chain,
                    job_id_href          => $job_id_href,
                    job_reservation_name => $active_parameter_href->{job_reservation_name},
                    log                  => $log,
                    max_parallel_processes_count_href =>
                      $file_info_href->{max_parallel_processes_count},
                    recipe_file_path     => $recipe_file_path,
                    recipe_files_tracker => $infile_index,
                    sample_id            => $sample_id,
                    submission_profile   => $active_parameter_href->{submission_profile},
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

    use MIP::Environment::Executable qw{ get_executable_base_command };
    use MIP::File_info qw{get_sample_file_attribute};
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Bwa qw{ bwa_mem run_bwamem };
    use MIP::Program::Samtools qw{ samtools_index samtools_stats samtools_sort samtools_view };
    use MIP::Sample_info qw{
      get_rg_header_line
      set_recipe_metafile_in_sample_info
      set_recipe_outfile_in_sample_info };
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
    my @infile_paths = @{ $io{in}{file_paths} };

    my $job_id_chain = get_recipe_attributes(
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

    my %file_info_sample = get_sample_file_attribute(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
    );

    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $job_id_chain,
                id                     => $sample_id,
                file_info_href         => $file_info_href,
                file_name_prefixes_ref => $file_info_sample{no_direction_infile_prefixes},
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
                temp_directory         => $temp_directory,
            }
        )
    );

    my $outfile_suffix        = $io{out}{file_suffix};
    my @outfile_name_prefixes = @{ $io{out}{file_name_prefixes} };
    my @outfile_paths         = @{ $io{out}{file_paths} };
    my @outfile_path_prefixes = @{ $io{out}{file_path_prefixes} };

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    my $output_format;
    my $uncompressed_bam_output;
    if ( $outfile_suffix eq q{.bam} ) {

        # Used downstream in samtools view
        $uncompressed_bam_output = 1;
        $output_format           = q{bam};
    }

    # Too avoid adjusting infile_index in submitting to jobs
    my $paired_end_tracker = 0;

    ## Perform per single-end or read pair
  INFILE_PREFIX:
    while ( my ( $infile_index, $infile_prefix ) =
        each @{ $file_info_sample{no_direction_infile_prefixes} } )
    {

        ## Assign file features
        my $outfile_name_prefix = $outfile_name_prefixes[$infile_index];
        my $outfile_path        = $outfile_paths[$infile_index];
        my $outfile_path_prefix = $outfile_path_prefixes[$infile_index];

        my $sequence_run_type = get_sample_file_attribute(
            {
                attribute      => q{sequence_run_type},
                file_info_href => $file_info_href,
                file_name      => $infile_prefix,
                sample_id      => $sample_id,
            }
        );

        ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
        my ( $recipe_file_path, $recipe_info_path ) = setup_script(
            {
                active_parameter_href => $active_parameter_href,
                core_number           => $recipe_resource{core_number},
                directory_id          => $sample_id,
                filehandle            => $filehandle,
                job_id_href           => $job_id_href,
                memory_allocation     => $recipe_resource{memory},
                recipe_directory      => $recipe_name,
                recipe_name           => $recipe_name,
                process_time          => $recipe_resource{time},
                temp_directory        => $temp_directory,
            }
        );

        ### SHELL:

        ### BWA MEM
        say {$filehandle} q{## Aligning reads with } . $recipe_name . q{ and sorting via Samtools};

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

        ## Construct RG header line
        my $rg_header_line = get_rg_header_line(
            {
                infile_prefix    => $infile_prefix,
                platform         => $active_parameter_href->{platform},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
                separator        => q{\t},
            }
        );
        ## Add missing "@RG"
        $rg_header_line = $DOUBLE_QUOTE . q{@RG} . q{\t} . $rg_header_line . $DOUBLE_QUOTE;

        # If post to ALTs in reference genome
        run_bwamem(
            {
                filehandle           => $filehandle,
                hla_typing           => $active_parameter_href->{bwa_mem_hla},
                infile_path          => $fastq_file_path,
                idxbase              => $referencefile_path,
                outfiles_prefix_path => $outfile_path_prefix,
                read_group_header    => $rg_header_line,
                second_infile_path   => $second_fastq_file_path,
                thread_number        => $recipe_resource{core_number},
            }
        );
        print {$filehandle} $PIPE . $SPACE;
        print {$filehandle} get_executable_base_command( { base_command => q{bwakit}, } ) . q{ sh}
          . $SPACE;
        say {$filehandle} $NEWLINE;

        ## Set samtools sort input; Sort directly from run-bwakit
        my $samtools_sort_infile = $outfile_path_prefix . $DOT . q{aln} . $outfile_suffix;

        ## Increment paired end tracker
        $paired_end_tracker++;

        ## Sort the output from run-bwamem
        samtools_sort(
            {
                filehandle            => $filehandle,
                infile_path           => $samtools_sort_infile,
                max_memory_per_thread => 2 . q{G},
                outfile_path          => $outfile_path,
                output_format         => $output_format,
                temp_file_path_prefix => catfile( $temp_directory, q{samtools_sort_temp} ),
                thread_number         => $recipe_resource{core_number},
                write_index           => 1,
            }
        );
        say {$filehandle} $NEWLINE;

        samtools_index(
            {
                filehandle  => $filehandle,
                infile_path => $outfile_path,
            }
        );
        say {$filehandle} $AMPERSAND . $NEWLINE;

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
            say {$filehandle} $AMPERSAND . $NEWLINE;
        }
        say {$filehandle} q{wait} . $NEWLINE;

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
                    base_command         => $profile_base_command,
                    case_id              => $case_id,
                    dependency_method    => q{sample_to_sample_parallel},
                    job_id_chain         => $job_id_chain,
                    job_id_href          => $job_id_href,
                    job_reservation_name => $active_parameter_href->{job_reservation_name},
                    log                  => $log,
                    max_parallel_processes_count_href =>
                      $file_info_href->{max_parallel_processes_count},
                    recipe_file_path     => $recipe_file_path,
                    recipe_files_tracker => $infile_index,
                    sample_id            => $sample_id,
                    submission_profile   => $active_parameter_href->{submission_profile},
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

    use MIP::Environment::Executable qw{ get_executable_base_command };

    ## Add percentage mapped reads to samtools stats output

    my @commands = ( get_executable_base_command( { base_command => q{perl}, } ), );

    # Execute perl
    print {$filehandle} join $SPACE, @commands;
    print {$filehandle} q? -ne '?;

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
    print {$filehandle} q{>} . $SPACE . $outfile_path . $SPACE;

    return;
}

1;
