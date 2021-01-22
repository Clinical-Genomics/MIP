package MIP::Recipes::Analysis::Preseq;

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
use MIP::Constants qw{ $DOT $LOG_NAME $NEWLINE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_preseq };

}

sub analysis_preseq {

## Function : Recipe for running Preseq
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

    use MIP::File_info qw{ get_io_files parse_io_outfiles };
    use MIP::Program::Preseq qw{ preseq_lc_extrap };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
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
        }
    );

    my $infile_name_prefix = $io{in}{file_name_prefix};
    my $infile_path_prefix = $io{in}{file_path_prefix};
    my $infile_suffix      = $io{in}{file_suffix};
    my $infile_path        = $infile_path_prefix . $infile_suffix;

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
                file_name_prefixes_ref => [$infile_name_prefix],
                outdata_dir            => $active_parameter_href->{outdata_dir},
                parameter_href         => $parameter_href,
                recipe_name            => $recipe_name,
            }
        )
    );

    my $outfile_path = $io{out}{file_paths}->[0];
    my $outfile_name = $io{out}{file_names}->[0];

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
        }
    );

    ### SHELL:
    say {$filehandle} q{## } . $recipe_name;

    ## Run Preseq lc_extrap
    say {$filehandle} q{## Run lc_extrap};
    preseq_lc_extrap(
        {
            bam_infile_path => $infile_path,
            filehandle      => $filehandle,
            outfile_path    => $outfile_path,
            verbose         => 1,
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
                dependency_method                 => q{sample_to_island},
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
