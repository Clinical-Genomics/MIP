package MIP::Recipes::Analysis::Arriba;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile devnull };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $COMMA $DOT $EMPTY_STR $LOG_NAME $NEWLINE $PIPE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_arriba };

}

## Constants
Readonly my $CHIM_JUNCTION_OVERHANG_MIN => 10;
Readonly my $CHIM_MULTIMAP_NMAX         => 50;
Readonly my $CHIM_SCORE_DROP_MAX        => 30;
Readonly my $CHIM_SEGMENT_MIN           => 10;
Readonly my $CHIM_SEGMENT_READ_GAP_MAX  => 3;
Readonly my $OUT_FILTER_MISMATCH_NMAX   => 3;
Readonly my $OUT_FILTER_MULTIMAP_NMAX   => 50;

sub analysis_arriba {

## Function : Detect and visualize RNA fusions with Arriba
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Recipe name
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
    use MIP::Program::Gnu::Coreutils qw{ gnu_rm gnu_tee };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Arriba qw{ arriba };
    use MIP::Program::Samtools qw{ samtools_index samtools_sort };
    use MIP::Program::Star qw{ star_aln };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Sample_info qw{
      get_pedigree_sample_id_attributes
      get_rg_header_line
      set_file_path_to_store
      set_recipe_metafile_in_sample_info
      set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Unix::Write_to_file qw{ unix_write_to_file };

    Readonly my $SCALING_FACTOR => 0.75;

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            chain_id       => q{MAIN},
            id             => $sample_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
        }
    );
    my @infile_paths = @{ $io{in}{file_paths} };

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

    ### Get infiles
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

    my @arriba_commands = star_aln(
        {
            align_sj_stitch_mismatch_nmax          => q{5 -1 5 5},
            align_spliced_mate_map_lmin_over_lmate => q{0.5},
            chim_junction_overhang_min             => $CHIM_JUNCTION_OVERHANG_MIN,
            chim_multimap_nmax                     => $CHIM_MULTIMAP_NMAX,
            chim_out_type                          => q{WithinBAM HardClip},
            chim_score_drop_max                    => $CHIM_SCORE_DROP_MAX,
            chim_score_junction_non_gtag           => 0,
            chim_score_min                         => 1,
            chim_score_separation                  => 1,
            chim_segment_min                       => $CHIM_SEGMENT_MIN,
            chim_segment_read_gap_max              => $CHIM_SEGMENT_READ_GAP_MAX,
            genome_dir_path                        => $referencefile_dir_path,
            infile_paths_ref                       => \@fastq_files,
            out_bam_compression                    => 0,
            outfile_name_prefix                    => $outfile_path_prefix . $DOT,
            out_filter_mismatch_nmax               => $OUT_FILTER_MISMATCH_NMAX,
            out_filter_multimap_nmax               => $OUT_FILTER_MULTIMAP_NMAX,
            out_sam_attr_rgline                    => $out_sam_attr_rgline,
            out_sam_type                           => q{BAM Unsorted},
            out_sam_unmapped                       => q{Within},
            pe_overlap_nbases_min => $active_parameter_href->{pe_overlap_nbases_min},
            stdout_data_type      => q{BAM_Unsorted},
            thread_number         => $recipe{core_number},
            two_pass_mode         => q{None},
        },
    );
    push @arriba_commands, $PIPE;

    my $star_outfile_path = $outfile_path_prefix . $UNDERSCORE . q{unsorted} . $DOT . q{bam};
    push @arriba_commands,
      gnu_tee(
        {
            outfile_paths_ref => [$star_outfile_path],
        }
      );
    push @arriba_commands, $PIPE;

    push @arriba_commands,
      arriba(
        {
            annotation_file_path       => $active_parameter_href->{transcript_annotation},
            blacklist_file_path        => $active_parameter_href->{arriba_blacklist_path},
            discarded_fusion_file_path => $outfile_path_prefix
              . $UNDERSCORE
              . q{discarded}
              . $outfile_suffix,
            genome_file_path         => $active_parameter_href->{human_genome_reference},
            infile_path              => catfile( dirname( devnull() ), q{stdin} ),
            known_fusion_file_path   => $active_parameter_href->{arriba_known_fusion_path},
            outfile_path             => $outfile_path,
            protein_domain_file_path => $active_parameter_href->{fusion_protein_domain_path},
            tag_file_path            => $active_parameter_href->{arriba_known_fusion_path},
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@arriba_commands,
            filehandle   => $filehandle,
            separator    => $SPACE,

        }
    );
    say {$filehandle} $NEWLINE;

    ## Sort BAM before visualization
    my $sorted_bam_file = $outfile_path_prefix . $DOT . q{bam};
    my $thread_memory   = int( $recipe{memory} / $recipe{core_number} );
    ## Allow for process overhead
    $thread_memory = int( $thread_memory * $SCALING_FACTOR );
    samtools_sort(
        {
            filehandle            => $filehandle,
            infile_path           => $star_outfile_path,
            max_memory_per_thread => $thread_memory . q{G},
            outfile_path          => $sorted_bam_file,
            temp_file_path_prefix => $temp_directory,
            thread_number         => $recipe{core_number},
        }
    );
    say {$filehandle} $NEWLINE;

    samtools_index(
        {
            bai_format  => 1,
            filehandle  => $filehandle,
            infile_path => $sorted_bam_file,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Remove unsorted bam
    gnu_rm(
        {
            filehandle  => $filehandle,
            infile_path => $star_outfile_path,
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

        set_file_path_to_store(
            {
                format           => q{meta},
                id               => $sample_id,
                path             => $outfile_path,
                recipe_name      => $recipe_name,
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
