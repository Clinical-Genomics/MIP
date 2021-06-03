package MIP::Recipes::Analysis::Dragen_dna;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname basename };
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
use MIP::Constants qw{ $COMMA $DOT $LOG_NAME $NEWLINE $SPACE $TAB $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_dragen_dna_align_vc analysis_dragen_dna_joint_calling };

}

## Constants
Readonly my $INDENT_DEPTH_3 => 3;

sub analysis_dragen_dna_align_vc {

## Function : Dragen dna align and single sample variant calling analysis
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Format::Dragen qw{ create_dragen_fastq_list_sample_id };
    use MIP::File_info qw{ get_io_files get_sample_file_attribute parse_io_outfiles };
    use MIP::Program::Dragen qw{ dragen_dna_analysis };
    use MIP::Program::Ssh qw{ ssh };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Sample_info qw{
      get_read_group
      set_recipe_metafile_in_sample_info
      set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

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

    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id               => $recipe{job_id_chain},
                id                     => $sample_id,
                file_info_href         => $file_info_href,
                file_name_prefixes_ref => [$sample_id],
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
            }
        )
    );

    my $outdir_path         = $io{out}{dir_path};
    my $outfile_name_prefix = $io{out}{file_name_prefix};
    my $outfile_path        = $io{out}{file_path};

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ($recipe_file_path) = setup_script(
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
            set_errexit           => 0,
        }
    );

    ### SHELL:

    say {$filehandle} q{## } . $recipe_name;

    ## Get all sample fastq info for dragen as csv file
    my @dragen_fastq_list_lines;

    # Too avoid adjusting infile_index in submitting to jobs
    my $paired_end_tracker = 0;

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

        ## Read group header line
        my %read_group = get_read_group(
            {
                infile_prefix    => $infile_prefix,
                platform         => $active_parameter_href->{platform},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );
        ## Add read groups to line
        my @read_groups = qw{ id sm lb lane };

        push @dragen_fastq_list_lines, join $COMMA, @read_group{@read_groups};

        ## Infile(s)
        my $fastq_file_path = $infile_paths[$paired_end_tracker];

        ## Add file paths to line
        $dragen_fastq_list_lines[-1] .= $COMMA . $fastq_file_path;

        # If second read direction is present
        if ( $sequence_run_type eq q{paired-end} ) {

            # Increment to collect correct read 2
            $paired_end_tracker = $paired_end_tracker + 1;
            my $second_fastq_file_path = $infile_paths[$paired_end_tracker];

            ## Add file paths to line
            $dragen_fastq_list_lines[-1] .= $COMMA . $second_fastq_file_path;
        }
        ## Increment paired end tracker
        $paired_end_tracker++;
    }

    my $fastq_list_file_path =
      catfile( $outdir_path, $sample_id . $UNDERSCORE . q{fastq_list.csv} );
    create_dragen_fastq_list_sample_id(
        {
            fastq_list_file_path => $fastq_list_file_path,
            fastq_list_lines_ref => \@dragen_fastq_list_lines,
            log                  => $log,
        }
    );

    my @ssh_cmd = ssh(
        {
            filehandle       => $filehandle,
            user_at_hostname => $active_parameter_href->{dragen_user_at_hostname},
        }
    );
    print {$filehandle} $SPACE;
    my @cmds = dragen_dna_analysis(
        {
            alignment_output_format  => q{BAM},
            dbsnp_file_path          => $active_parameter_href->{dragen_dbsnp},
            dragen_hash_ref_dir_path => $active_parameter_href->{dragen_hash_ref_dir_path},
            enable_bam_indexing      => 1,
            enable_duplicate_marking => 1,
            enable_map_align         => 1,
            enable_map_align_output  => 1,
            enable_sort              => 1,
            enable_variant_caller    => 1,
            fastq_list_sample_id     => $sample_id,
            fastq_list_file_path     => $fastq_list_file_path,
            filehandle               => $filehandle,
            force                    => 1,
            outdirectory_path        => $outdir_path,
            outfile_prefix           => $outfile_name_prefix,
            sample_id                => $sample_id,
            vc_emit_ref_confidence   => q{GVCF},
        }
    );
    say {$filehandle} $NEWLINE;

    _dragen_wait_loop(
        {
            cmd        => join( $SPACE, @ssh_cmd, @cmds ),
            filehandle => $filehandle,
        }
    );

    ## Close filehandleS
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

sub analysis_dragen_dna_joint_calling {

## Function : Dragen dna joint calling analysis
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Recipe name
##          : $sample_info_href        => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $profile_base_command;

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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Pedigree qw{ create_fam_file };
    use MIP::File_info qw{ get_io_files parse_io_outfiles };
    use MIP::Program::Dragen qw{ dragen_dna_analysis };
    use MIP::Program::Ssh qw{ ssh };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Sample_info qw{ set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );
    my $temp_dir = $active_parameter_href->{temp_directory};

    ## Set and get the io files per chain, id and stream
    my %io = parse_io_outfiles(
        {
            chain_id               => $recipe{job_id_chain},
            id                     => $case_id,
            file_info_href         => $file_info_href,
            file_name_prefixes_ref => [$case_id],
            outdata_dir            => $active_parameter_href->{outdata_dir},
            parameter_href         => $parameter_href,
            recipe_name            => $recipe_name,
        }
    );

    my $outdir_path         = $io{out}{dir_path};
    my $outdir_path_prefix  = $io{out}{dir_path_prefix};
    my $outfile_name_prefix = $io{out}{file_name_prefix};
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_suffix      = $io{out}{file_suffix};
    my $outfile_path        = $outfile_path_prefix . $outfile_suffix;

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $recipe{core_number},
            directory_id          => $case_id,
            filehandle            => $filehandle,
            job_id_href           => $job_id_href,
            memory_allocation     => $recipe{memory},
            process_time          => $recipe{time},
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
            set_errexit           => 0,
        }
    );

    ## Collect infiles for all sample_ids
    my @dragen_infile_paths;
    while ( my ( $sample_id_index, $sample_id ) = each @{ $active_parameter_href->{sample_ids} } ) {

        ## Get the io infiles per chain and id
        my %sample_io = get_io_files(
            {
                id             => $sample_id,
                file_info_href => $file_info_href,
                parameter_href => $parameter_href,
                recipe_name    => $recipe_name,
                stream         => q{in},
            }
        );
        my $infile_path_prefix = $sample_io{in}{file_path_prefix};
        my $infile_suffix      = $sample_io{in}{file_suffix};
        my $infile_path        = $infile_path_prefix . $infile_suffix . q{.gz};

        ## Store infile path for each sample_id
        push @dragen_infile_paths, $infile_path;
    }

    ### SHELL:
    say {$filehandle} q{## } . $recipe_name;

    my $case_file_path;

    if ( $parameter_href->{cache}{trio} ) {

        $case_file_path = catfile( $outdir_path_prefix, $case_id . $DOT . q{fam} );

        ## Create .fam file to be used in variant calling analyses
        create_fam_file(
            {
                case_id          => $case_id,
                fam_file_path    => $case_file_path,
                filehandle       => $filehandle,
                sample_ids_ref   => $active_parameter_href->{sample_ids},
                sample_info_href => $sample_info_href,
            }
        );
    }

    ## DRAGEN combine gvcfs
    my @ssh_cmds = ssh(
        {
            filehandle       => $filehandle,
            user_at_hostname => $active_parameter_href->{dragen_user_at_hostname},
        }
    );
    print {$filehandle} $SPACE;
    my @combine_cmds = dragen_dna_analysis(
        {
            dbsnp_file_path            => $active_parameter_href->{dragen_dbsnp},
            dragen_hash_ref_dir_path   => $active_parameter_href->{dragen_hash_ref_dir_path},
            enable_combinegvcfs        => 1,
            filehandle                 => $filehandle,
            force                      => 1,
            outdirectory_path          => $outdir_path,
            outfile_prefix             => $outfile_name_prefix . $UNDERSCORE . q{comb},
            sample_gvcf_file_paths_ref => \@dragen_infile_paths,
        }
    );
    say {$filehandle} $NEWLINE;

    _dragen_wait_loop(
        {
            cmd        => join( $SPACE, @ssh_cmds, @combine_cmds ),
            filehandle => $filehandle,
        }
    );
    say {$filehandle} $NEWLINE;

    ## DRAGEN joint calling
    my $joint_calling_infile =
      catfile( $outdir_path, $outfile_name_prefix . $UNDERSCORE . q{comb.gvcf.gz} );

    ssh(
        {
            filehandle       => $filehandle,
            user_at_hostname => $active_parameter_href->{dragen_user_at_hostname},
        }
    );
    print {$filehandle} $SPACE;
    my @joint_call_cmds = dragen_dna_analysis(
        {
            dbsnp_file_path            => $active_parameter_href->{dragen_dbsnp},
            disable_vcf_compression    => 1,
            dragen_hash_ref_dir_path   => $active_parameter_href->{dragen_hash_ref_dir_path},
            enable_joint_genotyping    => 1,
            filehandle                 => $filehandle,
            force                      => 1,
            pedigree_file_path         => $case_file_path,
            outdirectory_path          => $outdir_path,
            outfile_prefix             => $outfile_name_prefix,
            sample_gvcf_file_paths_ref => [$joint_calling_infile],
        }
    );
    say {$filehandle} $NEWLINE;

    _dragen_wait_loop(
        {
            cmd        => join( $SPACE, @ssh_cmds, @joint_call_cmds ),
            filehandle => $filehandle,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Close filehandleS
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe{mode} == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command                      => $profile_base_command,
                case_id                           => $case_id,
                dependency_method                 => q{sample_to_case},
                job_id_chain                      => $recipe{job_id_chain},
                job_id_href                       => $job_id_href,
                log                               => $log,
                max_parallel_processes_count_href =>
                  $file_info_href->{max_parallel_processes_count},
                recipe_file_path   => $recipe_file_path,
                sample_ids_ref     => \@{ $active_parameter_href->{sample_ids} },
                submission_profile => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

sub _dragen_wait_loop {

## Function : Wait for dragen processor to become available
## Returns  :
## Arguments: $cmd           => Command to retry
##          : $filehandle    => Filehandle to write to
##          : $max_retries   => Maixum number of retries
##          : $time_to_sleep => Time to sleep before retrying (seconds)

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $cmd;
    my $filehandle;

    ## Default(s)
    my $max_retries;
    my $time_to_sleep;

    my $tmpl = {
        cmd => {
            defined     => 1,
            required    => 1,
            store       => \$cmd,
            strict_type => 1,
        },
        filehandle => {
            required => 1,
            store    => \$filehandle,
        },
        max_retries => {
            allow       => qr{ \A\d+\z }sxm,
            default     => 60,
            store       => \$max_retries,
            strict_type => 1,
        },
        time_to_sleep => {
            allow       => qr{ \A\d+\z }sxm,
            default     => 60,
            store       => \$time_to_sleep,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    say {$filehandle} q{status=$?};
    say {$filehandle} q{time_to_sleep=} . $time_to_sleep;

    say {$filehandle} q{if [ $status -eq 0 ]};
    say {$filehandle} q{then};
    say {$filehandle} $TAB, q{echo "$cmd command was successful"};

    say {$filehandle} q{else};
    say {$filehandle} q?for i in {1..? . $max_retries . q?}?;
    say {$filehandle} $TAB . q{do};
    say {$filehandle} $TAB x 2 . q{echo "$cmd failed $i times: Retrying in $time_to_sleep"};
    say {$filehandle} $TAB x 2, qq{sleep $time_to_sleep};
    say {$filehandle} $TAB x 2, $cmd;
    say {$filehandle} $TAB x 2, q{status=$?};
    say {$filehandle} $TAB x 2, q{if [ $status -eq 0 ]};
    say {$filehandle} $TAB x 2 . q{then};
    say {$filehandle} $TAB x $INDENT_DEPTH_3, q{echo "$cmd command was successful"};
    say {$filehandle} $TAB x $INDENT_DEPTH_3, q{break};
    say {$filehandle} $TAB x 2, q{fi};
    say {$filehandle} $TAB, q{done};
    say {$filehandle} q{fi};

    return;
}
1;
