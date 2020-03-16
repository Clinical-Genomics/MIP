package MIP::Recipes::Analysis::Fastqc;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw(fileparse);
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{:all};
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME $NEWLINE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{Exporter};

    # Set the version for version checking
    our $VERSION = 1.17;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_fastqc };

}

sub analysis_fastqc {

## Function : Raw sequence quality analysis using FASTQC.
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

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    use MIP::Cluster qw{ update_core_number_to_seq_mode update_memory_allocation };
    use MIP::Environment::Cluster qw{ check_max_core_number };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Program::Gnu::Coreutils qw{ gnu_mkdir };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ print_wait submit_recipe };
    use MIP::Program::Fastqc qw{ fastqc };
    use MIP::Sample_info qw{ set_recipe_outfile_in_sample_info };
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
            temp_directory => $temp_directory,
        }
    );
    my @infile_paths         = @{ $io{in}{file_paths} };
    my @infile_name_prefixes = @{ $io{in}{file_name_prefixes} };

    my $job_id_chain = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            attribute      => q{chain},
        }
    );
    my $recipe_mode     = $active_parameter_href->{$recipe_name};
    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );
    my $core_number = $recipe_resource{core_number};

    ## Outpaths
    my $outsample_directory =
      catdir( $active_parameter_href->{outdata_dir}, $sample_id, $recipe_name );
    my @outfile_paths =
      map {
        catdir( $outsample_directory, $_ . $UNDERSCORE . q{fastqc}, q{fastqc_data.txt} )
      } @infile_name_prefixes;

    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id       => $job_id_chain,
                id             => $sample_id,
                file_info_href => $file_info_href,
                file_paths_ref => \@outfile_paths,
                parameter_href => $parameter_href,
                recipe_name    => $recipe_name,
                temp_directory => $temp_directory,
            }
        )
    );

    my $outdir_path_prefix    = $io{out}{dir_path_prefix};
    my @outfile_name_prefixes = @{ $io{out}{file_name_prefixes} };
    @outfile_paths = @{ $io{out}{file_paths} };

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

  INFILE_LANE:
    foreach my $infile ( @{ $infile_lane_prefix_href->{$sample_id} } ) {

        ## Update the number of cores to be used in the analysis according to sequencing mode requirements
        $core_number = update_core_number_to_seq_mode(
            {
                core_number => $core_number,
                sequence_run_type =>
                  $sample_info_href->{sample}{$sample_id}{file}{$infile}
                  {sequence_run_type},
            }
        );
    }

    ## Limit number of cores requested to the maximum number of cores available per node
    $core_number = check_max_core_number(
        {
            core_number_requested => $core_number,
            max_cores_per_node    => $active_parameter_href->{max_cores_per_node},
        }
    );

    ## Update memory to match number of infiles
    my $memory_allocation = update_memory_allocation(
        {
            parallel_processes        => $core_number,
            process_memory_allocation => $recipe_resource{core_number},
            node_ram_memory           => $active_parameter_href->{node_ram_memory},
        }
    );

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ($recipe_file_path) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $sample_id,
            filehandle                      => $filehandle,
            job_id_href                     => $job_id_href,
            log                             => $log,
            memory_allocation               => $memory_allocation,
            process_time                    => $recipe_resource{time},
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            sleep                           => 1,
            source_environment_commands_ref => $recipe_resource{load_env_ref},
            temp_directory                  => $temp_directory,
        }
    );

    ### SHELL:

    say {$filehandle} q{## Create output dir};
    gnu_mkdir(
        {
            filehandle       => $filehandle,
            indirectory_path => $outsample_directory,
            parents          => 1,
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## } . $recipe_name;

    my $process_batches_count = 1;

  INFILE:
    while ( my ( $index, $infile_path ) = each @infile_paths ) {

        $process_batches_count = print_wait(
            {
                filehandle            => $filehandle,
                max_process_number    => $core_number,
                process_batches_count => $process_batches_count,
                process_counter       => $index,
            }
        );

        fastqc(
            {
                extract           => 1,
                filehandle        => $filehandle,
                infile_path       => $infile_path,
                outdirectory_path => $outsample_directory,
            }
        );
        say {$filehandle} q{&}, $NEWLINE;

        ## Collect QC metadata info for active recipe for later use
        if ( $recipe_mode == 1 ) {

            set_recipe_outfile_in_sample_info(
                {
                    infile           => $outfile_name_prefixes[$index],
                    path             => $outfile_paths[$index],
                    outdirectory     => $outsample_directory,
                    recipe_name      => $recipe_name,
                    sample_id        => $sample_id,
                    sample_info_href => $sample_info_href,
                }
            );
        }
    }
    say {$filehandle} q{wait}, $NEWLINE;

    close $filehandle;

    if ( $recipe_mode == 1 ) {

        submit_recipe(
            {
                base_command            => $profile_base_command,
                dependency_method       => q{sample_to_island},
                case_id                 => $case_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_chain            => $job_id_chain,
                job_id_href             => $job_id_href,
                job_reservation_name    => $active_parameter_href->{job_reservation_name},
                log                     => $log,
                recipe_file_path        => $recipe_file_path,
                sample_id               => $sample_id,
                submission_profile      => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

1;
