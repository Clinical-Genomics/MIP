package MIP::Recipes::Build::Star_fusion_prerequisites;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $DOT $NEWLINE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.05;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ build_star_fusion_prerequisites };

}

sub build_star_fusion_prerequisites {

## Function : Creates the STAR-Fusion prerequisites for the human genome and transcriptome
## Returns  :
## Arguments: $active_parameter_href        => Active parameters for this analysis hash {REF}
##          : $case_id                      => Case id
##          : $file_info_href               => File info hash {REF}
##          : $human_genome_reference       => Human genome reference
##          : $infile_lane_prefix_href      => Infile(s) without the ".ending" {REF}
##          : $job_id_href                  => Job id hash {REF}
##          : $log                          => Log object
##          : $parameter_href               => Parameter hash {REF}
##          : $recipe_name                  => Program name
##          : $parameter_build_suffixes_ref => The rtg reference associated directory suffixes {REF}
##          : $profile_base_command         => Submission profile base command
##          : $sample_info_href             => Info on samples and case hash {REF}
##          : $temp_directory               => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $log;
    my $parameter_href;
    my $recipe_name;
    my $parameter_build_suffixes_ref;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $human_genome_reference;
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
        human_genome_reference => {
            default     => $arg_href->{active_parameter_href}{human_genome_reference},
            store       => \$human_genome_reference,
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
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        parameter_build_suffixes_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$parameter_build_suffixes_ref,
            strict_type => 1,
        },
        profile_base_command => {
            default     => q{sbatch},
            store       => \$profile_base_command,
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

    use MIP::Get::Parameter qw{ get_recipe_resources };
    use MIP::Program::Gnu::Coreutils qw{ gnu_mkdir };
    use MIP::Language::Shell qw{ check_exist_and_move_file };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Star_fusion qw{ star_fusion_prep_genome_lib };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Constants
    Readonly my $MAX_RANDOM_NUMBER => 100_00;
    Readonly my $NUMBER_OF_CORES   => $active_parameter_href->{max_cores_per_node};
    Readonly my $PROCESSING_TIME   => 30;

    ## Unpack parameters
    my $job_id_chain    = $parameter_href->{$recipe_name}{chain};
    my $recipe_mode     = $active_parameter_href->{$recipe_name};
    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## filehandleS
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Generate a random integer between 0-10,000.
    my $random_integer = int rand $MAX_RANDOM_NUMBER;

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ($recipe_file_path) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $NUMBER_OF_CORES,
            directory_id                    => $case_id,
            filehandle                      => $filehandle,
            job_id_href                     => $job_id_href,
            log                             => $log,
            memory_allocation               => $recipe_resource{memory},
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            process_time                    => $PROCESSING_TIME,
            source_environment_commands_ref => $recipe_resource{load_env_ref},
        }
    );

    $log->warn( q{Will try to create required }
          . $human_genome_reference
          . q{ STAR-Fusion files before executing }
          . $recipe_name );

    say {$filehandle} q{## Building Star-Fusion dir files};
    ## Get parameters
    my $star_fusion_directory_tmp =
      catdir( $active_parameter_href->{star_fusion_reference_genome}
          . $UNDERSCORE
          . $random_integer );

    # Create temp dir
    gnu_mkdir(
        {
            filehandle       => $filehandle,
            indirectory_path => $star_fusion_directory_tmp,
            parents          => 1,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Build genome lib
    star_fusion_prep_genome_lib(
        {
            filehandle            => $filehandle,
            fusion_annot_lib_path => catfile(
                $active_parameter_href->{star_fusion_genome_lib_dir},
                q{fusion_annot_lib.gz}
            ),
            gtf_path           => $active_parameter_href->{transcript_annotation},
            output_dir_path    => $star_fusion_directory_tmp,
            pfam_db_path       => $active_parameter_href->{star_fusion_pfam_db},
            referencefile_path => $human_genome_reference,
            thread_number      => $NUMBER_OF_CORES,
        }
    );
    say {$filehandle} $NEWLINE;

  PREREQ:
    foreach my $suffix ( @{$parameter_build_suffixes_ref} ) {

        my $intended_file_path =
          $active_parameter_href->{star_fusion_reference_genome} . $suffix;

        ## Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
        check_exist_and_move_file(
            {
                filehandle          => $filehandle,
                intended_file_path  => $intended_file_path,
                temporary_file_path => $star_fusion_directory_tmp,
            }
        );
    }

    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe_mode == 1 ) {

        submit_recipe(
            {
                base_command       => $profile_base_command,
                dependency_method  => q{island_to_samples},
                case_id            => $case_id,
                job_id_href        => $job_id_href,
                log                => $log,
                job_id_chain       => $job_id_chain,
                recipe_file_path   => $recipe_file_path,
                sample_ids_ref     => \@{ $active_parameter_href->{sample_ids} },
                submission_profile => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

1;
