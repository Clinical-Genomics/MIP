package MIP::Recipes::Build::Star_prerequisites;

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

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ build_star_prerequisites };

}

## Constants
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub build_star_prerequisites {

## Function : Creates the Star prerequisites for the human genome
## Returns  :
## Arguments: $active_parameter_href        => Active parameters for this analysis hash {REF}
##          : $family_id                    => Family id
##          : $file_info_href               => File info hash {REF}
##          : $human_genome_reference       => Human genome reference
##          : $infile_lane_prefix_href      => Infile(s) without the ".ending" {REF}
##          : $job_id_href                  => Job id hash {REF}
##          : $log                          => Log object
##          : $parameter_href               => Parameter hash {REF}
##          : $program_name                 => Program name
##          : $parameter_build_suffixes_ref => The rtg reference associated directory suffixes {REF}
##          : $sample_info_href             => Info on samples and family hash {REF}
##          : $temp_directory               => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $log;
    my $parameter_href;
    my $program_name;
    my $parameter_build_suffixes_ref;
    my $sample_info_href;

    ## Default(s)
    my $family_id;
    my $human_genome_reference;
    my $temp_directory;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            store       => \$family_id,
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
        program_name => {
            defined     => 1,
            required    => 1,
            store       => \$program_name,
            strict_type => 1,
        },
        parameter_build_suffixes_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$parameter_build_suffixes_ref,
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

    use MIP::Get::Parameter qw{ get_module_parameters };
    use MIP::Gnu::Coreutils qw{ gnu_mkdir };
    use MIP::Language::Shell qw{ check_exist_and_move_file };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_no_dependency_add_to_samples };
    use MIP::Program::Alignment::Star qw{ star_genome_generate };
    use MIP::Recipes::Build::Human_genome_prerequisites
      qw{ build_human_genome_prerequisites };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Constants
    Readonly my $NUMBER_OF_CORES   => $active_parameter_href->{max_cores_per_node};
    Readonly my $MAX_RANDOM_NUMBER => 100_00;
    Readonly my $PROCESSING_TIME   => 3;
    Readonly my $READ_LENGTH       => 150;

    ## Set program mode
    my $program_mode = $active_parameter_href->{$program_name};

    ## Unpack parameters
    my $job_id_chain = $parameter_href->{$program_name}{chain};
    my ( $core_number, $time, @source_environment_cmds ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
    );

    ## FILEHANDLES
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Generate a random integer between 0-10,000.
    my $random_integer = int rand $MAX_RANDOM_NUMBER;

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($recipe_file_path) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            FILEHANDLE                      => $FILEHANDLE,
            directory_id                    => $family_id,
            core_number                     => $NUMBER_OF_CORES,
            job_id_href                     => $job_id_href,
            log                             => $log,
            program_directory               => $program_name,
            program_name                    => $program_name,
            process_time                    => $PROCESSING_TIME,
            source_environment_commands_ref => \@source_environment_cmds,
        }
    );

    build_human_genome_prerequisites(
        {
            active_parameter_href   => $active_parameter_href,
            FILEHANDLE              => $FILEHANDLE,
            file_info_href          => $file_info_href,
            infile_lane_prefix_href => $infile_lane_prefix_href,
            job_id_href             => $job_id_href,
            log                     => $log,
            parameter_build_suffixes_ref =>
              \@{ $file_info_href->{human_genome_reference_file_endings} },
            parameter_href   => $parameter_href,
            program_name     => $program_name,
            random_integer   => $random_integer,
            sample_info_href => $sample_info_href,
        }
    );

    if ( $parameter_href->{star_aln_reference_genome}{build_file} == 1 ) {

        $log->warn( q{Will try to create required }
              . $human_genome_reference
              . q{ star files before executing }
              . $program_name );

        say {$FILEHANDLE} q{## Building Star dir files};
        ## Get parameters
        my $star_directory_tmp =
            $active_parameter_href->{star_aln_reference_genome}
          . $UNDERSCORE
          . $random_integer;

        # Create temp dir
        gnu_mkdir(
            {
                FILEHANDLE       => $FILEHANDLE,
                indirectory_path => $star_directory_tmp,
                parents          => 1,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        star_genome_generate(
            {
                FILEHANDLE      => $FILEHANDLE,
                fasta_path      => $human_genome_reference,
                genome_dir_path => $star_directory_tmp,
                gtf_path        => $active_parameter_href->{transcript_annotation},
                read_length     => $READ_LENGTH,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

      PREREQ:
        foreach my $suffix ( @{$parameter_build_suffixes_ref} ) {

            my $intended_file_path =
              $active_parameter_href->{star_aln_reference_genome} . $suffix;

            ## Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
            check_exist_and_move_file(
                {
                    FILEHANDLE          => $FILEHANDLE,
                    intended_file_path  => $intended_file_path,
                    temporary_file_path => $star_directory_tmp,
                }
            );
        }

        ## Ensure that this subrutine is only executed once
        $parameter_href->{star_aln_reference_genome}{build_file} = 0;
    }

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

    if ( $program_mode == 1 ) {

        submit_recipe(
            {
                dependency_method  => q{island_to_samples},
                family_id          => $family_id,
                job_id_href        => $job_id_href,
                log                => $log,
                job_id_chain       => $job_id_chain,
                recipe_file_path   => $recipe_file_path,
                sample_ids_ref     => \@{ $active_parameter_href->{sample_ids} },
                submission_profile => $active_parameter_href->{submission_profile},
            }
        );
    }
    return;
}

1;
