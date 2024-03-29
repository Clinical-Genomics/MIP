package MIP::Recipes::Analysis::Star_aln;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $ASTERISK $COMMA $DOT $EMPTY_STR $LOG_NAME $NEWLINE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_star_aln analysis_star_aln_mixed analysis_star_fusion_aln };

}

sub analysis_star_aln {

## Function : Alignment of fastq files using star aln
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
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
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File_info
      qw{ get_io_files get_sample_fastq_file_lanes get_sample_file_attribute parse_io_outfiles };
    use MIP::Program::Gnu::Coreutils qw{ gnu_mv gnu_rm };
    use MIP::Program::Samtools qw{ samtools_index };
    use MIP::Program::Star qw{ star_aln };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Sample_info qw{
      get_rg_header_line
      set_recipe_metafile_in_sample_info
      set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ## PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

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

## Build outfile_paths
    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Join infile lanes
    my $lanes_id = join $EMPTY_STR,
      get_sample_fastq_file_lanes(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
      );
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $recipe{job_id_chain},
                file_info_href         => $file_info_href,
                file_name_prefixes_ref =>
                  [ $sample_id . $UNDERSCORE . q{lanes} . $UNDERSCORE . $lanes_id ],
                id             => $sample_id,
                outdata_dir    => $active_parameter_href->{outdata_dir},
                parameter_href => $parameter_href,
                recipe_name    => $recipe_name,
            }
        )
    );
    my $outfile_name        = ${ $io{out}{file_names} }[0];
    my $outfile_path        = $io{out}{file_path};
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_suffix      = $io{out}{file_suffix};

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $recipe{core_number},
            directory_id          => $sample_id,
            filehandle            => $filehandle,
            job_id_href           => $job_id_href,
            memory_allocation     => $recipe{memory},
            process_time          => $recipe{time},
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
            temp_directory        => $temp_directory,
            ulimit_n              => $active_parameter_href->{star_ulimit_n},
        }
    );

    ### SHELL:
    say {$filehandle} q{## Performing fusion transcript detections using } . $recipe_name;

    ## Get infiles
    my $paired_end_tracker = 0;
    my @forward_files;
    my @reverse_files;
    my @read_groups;

    my %file_info_sample = get_sample_file_attribute(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
    );

    ## Perform per single-end or read pair
  INFILE_PREFIX:
    foreach my $infile_prefix ( @{ $file_info_sample{no_direction_infile_prefixes} } ) {

        my $sequence_run_type = get_sample_file_attribute(
            {
                attribute      => q{sequence_run_type},
                file_info_href => $file_info_href,
                file_name      => $infile_prefix,
                sample_id      => $sample_id,
            }
        );

        ## Add read one to file index array
        push @forward_files, $infile_paths[$paired_end_tracker];

        # If second read direction is present
        if ( $sequence_run_type eq q{paired-end} ) {

            # Increment to collect correct read 2 from infiles
            $paired_end_tracker++;

            ## Add read two file index array
            push @reverse_files, $infile_paths[$paired_end_tracker];
        }

        ## Increment paired end tracker
        $paired_end_tracker++;

        ## Construct RG header line
        my $rg_header_line = get_rg_header_line(
            {
                infile_prefix    => $infile_prefix,
                platform         => $active_parameter_href->{platform},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
                separator        => $SPACE,
            }
        );
        push @read_groups, $rg_header_line;
    }

    my @fastq_files =
      ( ( join $COMMA, @forward_files ), ( join $COMMA, @reverse_files ) );
    my $out_sam_attr_rgline = join $SPACE . $COMMA . $SPACE, @read_groups;

    my $referencefile_dir_path =
        $active_parameter_href->{star_aln_reference_genome}
      . $file_info_href->{star_aln_reference_genome}[0];

    my $out_sam_strand_field =
      $active_parameter_href->{library_type} eq q{unstranded} ? q{intronMotif} : undef;
    star_aln(
        {
            align_intron_max           => $active_parameter_href->{align_intron_max},
            align_mates_gap_max        => $active_parameter_href->{align_mates_gap_max},
            align_sjdb_overhang_min    => $active_parameter_href->{align_sjdb_overhang_min},
            chim_junction_overhang_min => $active_parameter_href->{chim_junction_overhang_min},
            chim_out_type              => $active_parameter_href->{chim_out_type},
            chim_segment_min           => $active_parameter_href->{chim_segment_min},
            filehandle                 => $filehandle,
            genome_dir_path            => $referencefile_dir_path,
            infile_paths_ref           => \@fastq_files,
            out_sam_attr_rgline        => $out_sam_attr_rgline,
            out_sam_strand_field       => $out_sam_strand_field,
            out_wig_norm               => q{None},
            out_wig_strand             => q{Unstranded},
            out_wig_type               => q{wiggle},
            outfile_name_prefix        => $outfile_path_prefix . $DOT,
            pe_overlap_nbases_min      => $active_parameter_href->{pe_overlap_nbases_min},
            quant_mode                 => q{GeneCounts},
            thread_number              => $recipe{core_number},
            two_pass_mode              => $active_parameter_href->{two_pass_mode},
        },
    );
    say {$filehandle} $NEWLINE;

    ## Rename bam file
    gnu_mv(
        {
            filehandle  => $filehandle,
            infile_path => $outfile_path_prefix
              . $DOT
              . q{Aligned.sortedByCoord.out}
              . $outfile_suffix,
            outfile_path => $outfile_path,
        }
    );
    say {$filehandle} $NEWLINE;

    samtools_index(
        {
            filehandle  => $filehandle,
            infile_path => $outfile_path,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Remove intermediary files
  FILE_TAG:
    foreach my $file_tag (qw{ _STARgenome _STARpass1 }) {

        gnu_rm(
            {
                filehandle  => $filehandle,
                force       => 1,
                infile_path => $outfile_path_prefix . $DOT . $file_tag,
                recursive   => 1,
            }
        );
        print {$filehandle} $NEWLINE;
    }

    ## Close filehandle
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe{mode} == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                infile           => $outfile_name,
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );
        my $qc_stats_outfile_path = $outfile_path_prefix . $DOT . q{Log.final.out};
        set_recipe_outfile_in_sample_info(
            {
                infile           => $outfile_name,
                path             => $qc_stats_outfile_path,
                recipe_name      => q{star_log},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command                      => $profile_base_command,
                case_id                           => $case_id,
                dependency_method                 => q{sample_to_sample},
                job_id_chain                      => $recipe{job_id_chain},
                job_id_href                       => $job_id_href,
                job_reservation_name              => $active_parameter_href->{job_reservation_name},
                log                               => $log,
                max_parallel_processes_count_href =>
                  $file_info_href->{max_parallel_processes_count},
                recipe_file_path   => $recipe_file_path,
                sample_id          => $sample_id,
                submission_profile => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

sub analysis_star_aln_mixed {

## Function : Alignment of mixed single and paired end fastq files using star aln
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
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
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File_info qw{ get_io_files get_sample_file_attribute parse_io_outfiles };
    use MIP::Program::Gnu::Coreutils qw{ gnu_mv gnu_rm };
    use MIP::Program::Samtools qw{ samtools_index };
    use MIP::Program::Star qw{ star_aln };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Sample_info qw{
      get_rg_header_line
      set_recipe_outfile_in_sample_info
      set_recipe_metafile_in_sample_info
      set_processing_metafile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ## PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            file_info_href => $file_info_href,
            id             => $sample_id,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
            temp_directory => $temp_directory,
        }
    );
    my @infile_paths = @{ $io{in}{file_paths} };
    my %recipe       = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
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
                chain_id               => $recipe{job_id_chain},
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

        # Collect paired-end or single-end sequence run mode
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
                core_number           => $recipe{core_number},
                directory_id          => $sample_id,
                filehandle            => $filehandle,
                job_id_href           => $job_id_href,
                memory_allocation     => $recipe{memory},
                recipe_directory      => $recipe_name,
                recipe_name           => $recipe_name,
                process_time          => $recipe{time},
                temp_directory        => $temp_directory,
                ulimit_n              => $active_parameter_href->{star_ulimit_n},
            }
        );

        ### SHELL

        ## Star aln
        say {$filehandle} q{## Aligning reads with } . $recipe_name;

        ### Get parameters
        ## Infile(s)
        my @fastq_files = $infile_paths[$paired_end_tracker];

        # If second read direction is present
        if ( $sequence_run_type eq q{paired-end} ) {

            # Increment to collect correct read 2 from %infile
            $paired_end_tracker++;
            push @fastq_files, $infile_paths[$paired_end_tracker];
        }
        my $referencefile_dir_path =
            $active_parameter_href->{star_aln_reference_genome}
          . $file_info_href->{star_aln_reference_genome}[0];

        ## Construct RG header line
        my $out_sam_attr_rgline = get_rg_header_line(
            {
                infile_prefix    => $infile_prefix,
                platform         => $active_parameter_href->{platform},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
                separator        => $SPACE,
            }
        );

        my $out_sam_strand_field =
          $active_parameter_href->{library_type} eq q{unstranded}
          ? q{intronMotif}
          : undef;
        star_aln(
            {
                filehandle                 => $filehandle,
                align_intron_max           => $active_parameter_href->{align_intron_max},
                align_mates_gap_max        => $active_parameter_href->{align_mates_gap_max},
                align_sjdb_overhang_min    => $active_parameter_href->{align_sjdb_overhang_min},
                chim_junction_overhang_min => $active_parameter_href->{chim_junction_overhang_min},
                chim_out_type              => $active_parameter_href->{chim_out_type},
                chim_segment_min           => $active_parameter_href->{chim_segment_min},
                genome_dir_path            => $referencefile_dir_path,
                infile_paths_ref           => \@fastq_files,
                out_sam_attr_rgline        => $out_sam_attr_rgline,
                out_sam_strand_field       => $out_sam_strand_field,
                outfile_name_prefix        => $outfile_path_prefix . $DOT,
                pe_overlap_nbases_min      => $active_parameter_href->{pe_overlap_nbases_min},
                thread_number              => $recipe{core_number},
                two_pass_mode              => $active_parameter_href->{two_pass_mode},
            },
        );
        say {$filehandle} $NEWLINE;

        ## Increment paired end tracker
        $paired_end_tracker++;

        ## Rename bam file
        gnu_mv(
            {
                filehandle  => $filehandle,
                infile_path => $outfile_path_prefix
                  . $DOT
                  . q{Aligned.sortedByCoord.out}
                  . $outfile_suffix,
                outfile_path => $outfile_path,
            }
        );
        say {$filehandle} $NEWLINE;

        samtools_index(
            {
                filehandle  => $filehandle,
                infile_path => $outfile_path,
            }
        );
        say {$filehandle} $NEWLINE;

        ## Remove intermediary files
      FILE_TAG:
        foreach my $file_tag (qw{ _STARgenome _STARpass1 }) {

            gnu_rm(
                {
                    filehandle  => $filehandle,
                    force       => 1,
                    infile_path => $outfile_path_prefix . $DOT . $file_tag,
                    recursive   => 1,
                }
            );
            print {$filehandle} $NEWLINE;
        }

        ## Close filehandles
        close $filehandle or $log->logcroak(q{Could not close filehandle});

        if ( $recipe{mode} == 1 ) {

            ## Collect QC metadata info for later use
            set_recipe_outfile_in_sample_info(
                {
                    infile           => $outfile_name_prefix,
                    path             => $outfile_path,
                    recipe_name      => $recipe_name,
                    sample_id        => $sample_id,
                    sample_info_href => $sample_info_href,
                }
            );

            my $star_aln_log = $outfile_path_prefix . $DOT . q{Log.final.out};
            set_recipe_metafile_in_sample_info(
                {
                    infile           => $outfile_name_prefix,
                    metafile_tag     => q{log},
                    path             => $star_aln_log,
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
                    job_id_chain         => $recipe{job_id_chain},
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

sub analysis_star_fusion_aln {

## Function : Alignment of fastq files using star aln
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $case_id               => Family id
##          : $file_info_href        => File_info hash {REF}
##          : $job_id_href           => Job id hash {REF}
##          : $parameter_href        => Parameter hash {REF}
##          : $profile_base_command  => Submission profile base command
##          : $recipe_name           => Program name
##          : $sample_id             => Sample id
##          : $sample_info_href      => Info on samples and case hash {REF}
##          : $temp_directory        => Temporary directory

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
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File_info
      qw{ get_io_files get_sample_fastq_file_lanes get_sample_file_attribute parse_io_outfiles };
    use MIP::Program::Gnu::Coreutils qw{ gnu_mv gnu_rm };
    use MIP::Program::Star qw{ star_aln };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Sample_info qw{
      get_rg_header_line
      set_recipe_metafile_in_sample_info
      set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ## PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

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

    ## Build outfile_paths
    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Join infile lanes
    my $lanes_id = join $EMPTY_STR,
      get_sample_fastq_file_lanes(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
      );
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $recipe{job_id_chain},
                file_info_href         => $file_info_href,
                file_name_prefixes_ref =>
                  [ $sample_id . $UNDERSCORE . q{lanes} . $UNDERSCORE . $lanes_id ],
                id             => $sample_id,
                outdata_dir    => $active_parameter_href->{outdata_dir},
                parameter_href => $parameter_href,
                recipe_name    => $recipe_name,
            }
        )
    );
    my $outfile_name        = ${ $io{out}{file_names} }[0];
    my $outfile_path        = $io{out}{file_path};
    my $outfile_path_prefix = $io{out}{file_path_prefix};

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $recipe{core_number},
            directory_id          => $sample_id,
            filehandle            => $filehandle,
            job_id_href           => $job_id_href,
            memory_allocation     => $recipe{memory},
            process_time          => $recipe{time},
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
            temp_directory        => $temp_directory,
            ulimit_n              => $active_parameter_href->{star_ulimit_n},
        }
    );

    ### SHELL:
    say {$filehandle} q{## Aligning reads for fusion transcript detections using } . $recipe_name;

    ## Get infiles
    my $paired_end_tracker = 0;
    my @forward_files;
    my @reverse_files;
    my @read_groups;

    my %file_info_sample = get_sample_file_attribute(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
    );

    ## Perform per single-end or read pair
  INFILE_PREFIX:
    foreach my $infile_prefix ( @{ $file_info_sample{no_direction_infile_prefixes} } ) {

        my $sequence_run_type = get_sample_file_attribute(
            {
                attribute      => q{sequence_run_type},
                file_info_href => $file_info_href,
                file_name      => $infile_prefix,
                sample_id      => $sample_id,
            }
        );

        ## Add read one to file index array
        push @forward_files, $infile_paths[$paired_end_tracker];

        # If second read direction is present
        if ( $sequence_run_type eq q{paired-end} ) {

            # Increment to collect correct read 2 from infiles
            $paired_end_tracker++;

            ## Add read two file index array
            push @reverse_files, $infile_paths[$paired_end_tracker];
        }

        ## Increment paired end tracker
        $paired_end_tracker++;

        ## Construct RG header line
        my $rg_header_line = get_rg_header_line(
            {
                infile_prefix    => $infile_prefix,
                platform         => $active_parameter_href->{platform},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
                separator        => $SPACE,
            }
        );
        push @read_groups, $rg_header_line;
    }

    my @fastq_files =
      ( ( join $COMMA, @forward_files ), ( join $COMMA, @reverse_files ) );
    my $out_sam_attr_rgline = join $SPACE . $COMMA . $SPACE, @read_groups;

    my $referencefile_dir_path =
      catdir( $active_parameter_href->{star_fusion_genome_lib_dir}, q{ref_genome.fa.star.idx} );

    my $out_sam_strand_field =
      $active_parameter_href->{library_type} eq q{unstranded} ? q{intronMotif} : undef;

    ## Parameters taken from https://github.com/STAR-Fusion/STAR-Fusion/wiki#KickstartMode
    ## Constants
    Readonly my $ALIGN_INTRON_MAX             => 100_000;
    Readonly my $ALIGN_MATES_GAP_MAX          => 100_000;
    Readonly my $ALIGN_SJDB_OVERHANG_MIN      => 10;
    Readonly my $ALIGN_SPLICED_MATEMAP_LMIN   => 30;
    Readonly my $CHIM_JUNCTION_OVERHANG_MIN   => 8;
    Readonly my $CHIM_MULTIMAP_NMAX           => 20;
    Readonly my $CHIM_MULTIMAP_SCORE_RANGE    => 3;
    Readonly my $CHIM_NONCHIM_SCORE_DROP_MIN  => 10;
    Readonly my $CHIM_SCORE_JUNCTION_NON_GTAG => -4;
    Readonly my $CHIM_SEGMENT_MIN             => 12;
    Readonly my $PE_OVERLAP_MMP               => 0.1;
    Readonly my $PE_OVERLAP_NBASES_MIN        => 12;

    star_aln(
        {
            align_insertion_flush                  => q{Right},
            align_intron_max                       => $ALIGN_INTRON_MAX,
            align_mates_gap_max                    => $ALIGN_MATES_GAP_MAX,
            align_sj_stitch_mismatch_nmax          => q{5 -1 5 5},
            align_sjdb_overhang_min                => $ALIGN_SJDB_OVERHANG_MIN,
            align_spliced_mate_map_lmin            => $ALIGN_SPLICED_MATEMAP_LMIN,
            align_spliced_mate_map_lmin_over_lmate => 0,
            chim_junction_overhang_min             => $CHIM_JUNCTION_OVERHANG_MIN,
            chim_multimap_nmax                     => $CHIM_MULTIMAP_NMAX,
            chim_multimap_score_range              => $CHIM_MULTIMAP_SCORE_RANGE,
            chim_nonchim_score_drop_min            => $CHIM_NONCHIM_SCORE_DROP_MIN,
            chim_out_junction_format               => 1,
            chim_score_junction_non_gtag           => $CHIM_SCORE_JUNCTION_NON_GTAG,
            chim_segment_min                       => $CHIM_SEGMENT_MIN,
            filehandle                             => $filehandle,
            genome_dir_path                        => $referencefile_dir_path,
            infile_paths_ref                       => \@fastq_files,
            out_sam_attr_rgline                    => $out_sam_attr_rgline,
            out_sam_strand_field                   => q{intronMotif},
            out_sam_unmapped                       => q{Within},
            outfile_name_prefix                    => $outfile_path_prefix . $DOT,
            pe_overlap_mmp                         => $PE_OVERLAP_MMP,
            pe_overlap_nbases_min                  => $PE_OVERLAP_NBASES_MIN,
            quant_mode                             => q{GeneCounts},
            thread_number                          => $recipe{core_number},
            two_pass_mode                          => q{None},
        },
    );
    say {$filehandle} $NEWLINE;

    ## Rename bam file
    gnu_mv(
        {
            filehandle   => $filehandle,
            infile_path  => $outfile_path_prefix . $DOT . q{Chimeric.out.junction},
            outfile_path => $outfile_path,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Close filehandle
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe{mode} == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                infile           => $outfile_name,
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command                      => $profile_base_command,
                case_id                           => $case_id,
                dependency_method                 => q{sample_to_sample},
                job_id_chain                      => $recipe{job_id_chain},
                job_id_href                       => $job_id_href,
                job_reservation_name              => $active_parameter_href->{job_reservation_name},
                log                               => $log,
                max_parallel_processes_count_href =>
                  $file_info_href->{max_parallel_processes_count},
                recipe_file_path   => $recipe_file_path,
                sample_id          => $sample_id,
                submission_profile => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}
1;
