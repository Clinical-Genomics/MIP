package MIP::Recipes::Analysis::Tiddit;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

## MIPs lib/
use MIP::Constants qw{ $AMPERSAND $ASTERISK $LOG_NAME $NEWLINE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_tiddit };

}

sub analysis_tiddit {

## Function : Call structural variants using Tiddit
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => The file info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $reference_dir           => MIP reference directory
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory

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
    my $reference_dir;
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
        reference_dir => {
            default     => $arg_href->{active_parameter_href}{reference_dir},
            store       => \$reference_dir,
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

    use MIP::Cluster qw{ get_core_number update_memory_allocation };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ print_wait submit_recipe };
    use MIP::Program::Svdb qw{ svdb_merge };
    use MIP::Program::Tiddit qw{ tiddit_sv };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Sample_info qw{ set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    my $max_cores_per_node = $active_parameter_href->{max_cores_per_node};
    my $modifier_core_number =
      scalar( @{ $active_parameter_href->{sample_ids} } );

    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

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

    my $outfile_suffix      = $io{out}{file_suffix};
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_path        = $io{out}{file_path};

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Update number of cores and memory depending on how many samples that are processed
    my $core_number = get_core_number(
        {
            max_cores_per_node   => $max_cores_per_node,
            modifier_core_number => $modifier_core_number,
            recipe_core_number   => $recipe{core_number},
        }
    );
    my $memory_allocation = update_memory_allocation(
        {
            node_ram_memory           => $active_parameter_href->{node_ram_memory},
            parallel_processes        => $core_number,
            process_memory_allocation => $recipe{memory},
        }
    );

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $core_number,
            directory_id          => $case_id,
            filehandle            => $filehandle,
            job_id_href           => $job_id_href,
            memory_allocation     => $memory_allocation,
            process_time          => $recipe{time},
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
            temp_directory        => $temp_directory,
        }
    );

    ### SHELL:

    my %tiddit_sample_file_info;

    ## Collect infiles for all sample_ids to enable migration to temporary directory
  SAMPLE_ID:
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
        my $infile_path        = $infile_path_prefix . $infile_suffix;

        $tiddit_sample_file_info{$sample_id}{in}  = $infile_path;
        $tiddit_sample_file_info{$sample_id}{out} = $outfile_path_prefix . $UNDERSCORE . $sample_id;
    }
    say {$filehandle} q{wait}, $NEWLINE;

    # Start counter
    my $process_batches_count = 1;

    ## Tiddit sv calling per sample id
  SAMPLE_ID:
    while ( my ( $sample_id_index, $sample_id ) = each @{ $active_parameter_href->{sample_ids} } ) {

        $process_batches_count = print_wait(
            {
                filehandle            => $filehandle,
                max_process_number    => $core_number,
                process_batches_count => $process_batches_count,
                process_counter       => $sample_id_index,
            }
        );

        ## Tiddit
        tiddit_sv(
            {
                filehandle                      => $filehandle,
                infile_path                     => $tiddit_sample_file_info{$sample_id}{in},
                minimum_number_supporting_pairs =>
                  $active_parameter_href->{tiddit_minimum_number_supporting_pairs},
                outfile_path_prefix => $tiddit_sample_file_info{$sample_id}{out},
                referencefile_path  => $active_parameter_href->{human_genome_reference},
            }
        );
        say {$filehandle} $AMPERSAND . $SPACE . $NEWLINE;
    }
    say {$filehandle} q{wait}, $NEWLINE;

    ## Get parameters
    ## Tiddit sample outfiles needs to be lexiographically sorted for svdb merge
    my @svdb_infile_paths =
      map { $tiddit_sample_file_info{$_}{out} . $outfile_suffix }
      @{ $active_parameter_href->{sample_ids} };

    svdb_merge(
        {
            filehandle       => $filehandle,
            infile_paths_ref => \@svdb_infile_paths,
            notag            => 1,
            stdoutfile_path  => $outfile_path,
        }
    );
    say {$filehandle} $NEWLINE;

    close $filehandle;

    if ( $recipe{mode} == 1 ) {

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
                job_reservation_name              => $active_parameter_href->{job_reservation_name},
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

1;
