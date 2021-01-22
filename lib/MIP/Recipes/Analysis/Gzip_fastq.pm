package MIP::Recipes::Analysis::Gzip_fastq;

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
use MIP::Constants qw{ $DOT $LOG_NAME $NEWLINE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_gzip_fastq };

}

sub analysis_gzip_fastq {

## Function : Gzips fastq files
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
            defined     => 1,
            default     => {},
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

    use MIP::Cluster qw{ update_core_number_to_seq_mode };
    use MIP::Environment::Cluster qw{ check_max_core_number };
    use MIP::File_info
      qw{ get_io_files get_is_sample_files_compressed get_sample_file_attribute parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Gzip qw{ gzip };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Script::Setup_script qw{ setup_script };

    my $is_files_compressed = get_is_sample_files_compressed(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
    );
    ## No uncompressed fastq infiles for this sample_id
    return if ($is_files_compressed);

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
    my $indir_path_prefix    = $io{in}{dir_path_prefix};
    my @infile_names         = @{ $io{in}{file_names} };
    my @infile_name_prefixes = @{ $io{in}{file_name_prefixes} };
    my @infile_paths         = @{ $io{in}{file_paths} };

    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );
    my $core_number = $recipe{core_number};

    ## Outpaths
    my @outfile_paths =
      map { catdir( $indir_path_prefix, $_ . $DOT . q{fastq.gz} ) } @infile_name_prefixes;

    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id       => $recipe{job_id_chain},
                id             => $sample_id,
                file_info_href => $file_info_href,
                file_paths_ref => \@outfile_paths,
                parameter_href => $parameter_href,
                recipe_name    => $recipe_name,
            }
        )
    );

    ## Adjust according to number of infiles to process
    # One full lane on Hiseq takes approx. 2 h for gzip to process
    my $time = $recipe{time} * scalar @infile_names;

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    my %file_info_sample = get_sample_file_attribute(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
    );

  INFILE_LANE:
    foreach my $infile_prefix ( @{ $file_info_sample{no_direction_infile_prefixes} } ) {

        my $sequence_run_type = get_sample_file_attribute(
            {
                attribute      => q{sequence_run_type},
                file_info_href => $file_info_href,
                file_name      => $infile_prefix,
                sample_id      => $sample_id,
            }
        );
        ## Update the number of cores to be used in the analysis according to sequencing mode requirements
        $core_number = update_core_number_to_seq_mode(
            {
                core_number       => $core_number,
                sequence_run_type => $sequence_run_type,
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

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ($recipe_file_path) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $core_number,
            directory_id          => $sample_id,
            filehandle            => $filehandle,
            job_id_href           => $job_id_href,
            memory_allocation     => $recipe{memory},
            process_time          => $time,
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
        }
    );

    my $process_batches_count = 1;

# Used to print wait at the right times since infiles cannot be used (can be a mixture of .gz and .fast files)
    my $uncompressed_file_counter = 0;

    ## Gzip
    say {$filehandle} q{## } . $recipe_name;

  INFILE:
    while ( my ( $infile_index, $infile ) = each @infile_names ) {

        my $is_file_compressed = get_sample_file_attribute(
            {
                attribute      => q{is_file_compressed},
                file_info_href => $file_info_href,
                file_name      => $infile,
                sample_id      => $sample_id,
            }
        );
        ## For files ending with ".fastq" required since there can be a mixture (also .fastq.gz) within the sample dir
        if ( not $is_file_compressed ) {

            ## Using only $active_parameter{max_cores_per_node} cores
            if ( $uncompressed_file_counter ==
                $process_batches_count * $active_parameter_href->{max_cores_per_node} )
            {

                say {$filehandle} q{wait}, $NEWLINE;
                $process_batches_count = $process_batches_count + 1;
            }

            ## Perl wrapper for writing gzip recipe to $filehandle
            gzip(
                {
                    filehandle       => $filehandle,
                    infile_paths_ref => [ $infile_paths[$infile_index] ],
                }
            );
            say {$filehandle} q{&};
            $uncompressed_file_counter++;
        }
    }
    print {$filehandle} $NEWLINE;
    say {$filehandle} q{wait}, $NEWLINE;

    if ( $recipe{mode} == 1 ) {

        submit_recipe(
            {
                base_command         => $profile_base_command,
                case_id              => $case_id,
                dependency_method    => q{island_to_sample},
                job_id_chain         => $recipe{job_id_chain},
                job_id_href          => $job_id_href,
                job_reservation_name => $active_parameter_href->{job_reservation_name},
                log                  => $log,
                recipe_file_path     => $recipe_file_path,
                sample_id            => $sample_id,
                submission_profile   => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

1;
