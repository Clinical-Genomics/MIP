package MIP::Recipes::Analysis::Gens_generatedata;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use POSIX qw{ floor };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ %ANALYSIS $ASTERISK $DOT $LOG_NAME $NEWLINE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_gens_generatedata };

}

sub analysis_gens_generatedata {

## Function : Gens generatedata analysis recipe
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
        case_id_ref => {
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

    use MIP::File_info qw{ get_io_files parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Gens qw{ gens_generatedata };
    use MIP::Recipe qw{ parse_recipe_prerequisites };
    use MIP::Sample_info qw{ set_file_path_to_store set_recipe_outfile_in_sample_info };
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
    my $infile_name_prefix     = $io{in}{file_name_prefix};
    my $infile_tsv_path =$io{in}{file_path};

    my %vcf_io = get_io_files(
        {
            id             => $sample_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => q{deepvariant},
            stream         => q{out},
        }
    );
    my $infile_vcf_path = $vcf_io{out}{file_path};

    ## Get module parameters
    my %recipe = parse_recipe_prerequisites(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Outpaths
    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id         => $recipe{job_id_chain},
                id               => $sample_id,
                file_info_href   => $file_info_href,
                file_name_prefix => $infile_name_prefix,
                iterators_ref    => [qw{ baf cov }],
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                recipe_name      => $recipe_name,
            }
        )
    );
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my %outfile_path        = %{ $io{out}{file_path_href} };

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
        }
    );

    ### SHELL:

    ## Gens generatedata
    say {$filehandle} q{## Gens generatedata};

    ## generate_gens_data.pl
    gens_generatedata(
        {
            filehandle       => $filehandle,
            gnomad_positions => $active_parameter_href->{gens_gnomad_positions},
            infile_tsv_path  => $infile_tsv_path,
            infile_vcf_path  => $infile_vcf_path,
            outfile_prefix   => $outfile_path_prefix,
        }
    );

    close $filehandle;

    if ( $recipe{mode} == 1 ) {

        ## Collect QC metadata info for later use
        set_recipe_outfile_in_sample_info(
            {
                infile           => $infile_vcf_path,
                path             => $outfile_path{baf},
                recipe_name      => $recipe_name,
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

      FILE_TAG:
        foreach my $file_tag ( keys %outfile_path ) {

            set_file_path_to_store(
                {
                    format           => q{bed},
                    id               => $sample_id,
                    path             => $outfile_path{$file_tag},
                    path_index       => $outfile_path{$file_tag} . q{.tbi},
                    recipe_name      => $recipe_name,
                    sample_info_href => $sample_info_href,
                    tag              => $file_tag,
                }
            );

        }

        submit_recipe(
            {
                base_command         => $profile_base_command,
                case_id              => $case_id,
                dependency_method    => q{sample_to_island},
                job_id_chain         => $recipe{job_id_chain},
                job_id_href          => $job_id_href,
                job_reservation_name => $active_parameter_href->{job_reservation_name},
                log                  => $log,
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
