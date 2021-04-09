package MIP::Recipes::Analysis::Star_fusion;

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
use MIP::Constants qw{ $ASTERISK $EMPTY_STR $LOG_NAME $NEWLINE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_star_fusion };

}

sub analysis_star_fusion {

## Function : Analysis recipe for star-fusion v1.8.0
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

    use MIP::File::Format::Star_fusion qw{ create_star_fusion_sample_file };
    use MIP::File_info qw{ get_io_files parse_io_outfiles get_sample_fastq_file_lanes };
    use MIP::Program::Gnu::Coreutils qw{ gnu_mv gnu_rm };
    use MIP::Program::Star_fusion qw{ star_fusion };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Sample_info qw{ set_file_path_to_store set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ## PREPROCESSING:

    ## Star fusion has a fixed sample_prefix
    Readonly my $STAR_FUSION_PREFIX => q{star-fusion.fusion_predictions};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            chain_id       => q{MAIN},
            file_info_href => $file_info_href,
            id             => $sample_id,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
            temp_directory => $temp_directory,
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
    my $outdir_path    = $io{out}{dir_path};
    my $outfile_path   = $io{out}{file_path};
    my $outfile_suffix = $io{out}{file_suffix};

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

# Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
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

    ## Star-fusion
    say {$filehandle} q{## Performing fusion transcript detections using } . $recipe_name;

    ## Create sample file
    my $sample_files_path = catfile( $outdir_path, $sample_id . q{_file.txt} );
    create_star_fusion_sample_file(
        {
            filehandle        => $filehandle,
            file_info_href    => $file_info_href,
            infile_paths_ref  => \@infile_paths,
            samples_file_path => $sample_files_path,
            sample_id         => $sample_id,
        }
    );
    say {$filehandle} $NEWLINE;

    star_fusion(
        {
            cpu                   => $recipe{core_number},
            examine_coding_effect => 1,
            filehandle            => $filehandle,
            fusion_inspector      => q{inspect},
            genome_lib_dir_path   => $active_parameter_href->{star_fusion_genome_lib_dir},
            min_junction_reads    => $active_parameter_href->{star_fusion_min_junction_reads},
            output_directory_path => $outdir_path,
            samples_file_path     => $sample_files_path,
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Rename outfile};
    my $star_fusion_outfile_path = catfile( $outdir_path, $STAR_FUSION_PREFIX . $outfile_suffix );
    gnu_mv(
        {
            filehandle   => $filehandle,
            infile_path  => $star_fusion_outfile_path,
            outfile_path => $outfile_path,
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Remove intermediary files};
  FILE_TAG:
    foreach my $file_tag (qw{ _STARgenome _STARpass1 _starF_checkpoints }) {

        gnu_rm(
            {
                filehandle  => $filehandle,
                force       => 1,
                infile_path => catdir( $outdir_path, $file_tag ),
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
