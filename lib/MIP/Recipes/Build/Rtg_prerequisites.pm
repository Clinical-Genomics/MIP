package MIP::Recipes::Build::Rtg_prerequisites;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use autodie qw{ :all };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use File::Spec::Functions qw{ catdir catfile };

## CPANM
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ build_rtg_prerequisites };

}

##Constants
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub build_rtg_prerequisites {

## Function : Creates the Rtg prerequisites for the human genome
## Returns  :
## Arguments: $active_parameter_href      => Active parameters for this analysis hash {REF}
##          : $family_id                  => Family id
##          : $file_info_href             => File info hash {REF}
##          : $human_genome_reference     => Human genome reference
##          : $infile_lane_prefix_href    => Infile(s) without the ".ending" {REF}
##          : $job_id_href                => Job id hash {REF}
##          : $outaligner_dir             => Outaligner_dir used in the analysis
##          : $parameter_href             => Parameter hash {REF}
##          : $program_name               => Program name
##          : $rtg_directory_suffixes_ref => The rtg reference associated directory suffixes {REF}
##          : $sample_info_href           => Info on samples and family hash {REF}
##          : $temp_directory             => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $program_name;
    my $rtg_directory_suffixes_ref;
    my $sample_info_href;

    ## Default(s)
    my $family_id;
    my $human_genome_reference;
    my $outaligner_dir;
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
            default =>
              $arg_href->{active_parameter_href}{human_genome_reference},
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
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            store       => \$outaligner_dir,
            strict_type => 1,
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
        rtg_directory_suffixes_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$rtg_directory_suffixes_ref,
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

    use MIP::Language::Shell qw{ check_exist_and_move_file };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_no_dependency_add_to_samples };
    use MIP::Program::Qc::Rtg qw{ rtg_format };
    use MIP::Recipes::Build::Human_genome_prerequisites
      qw{ build_human_genome_prerequisites };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Constants
    Readonly my $MAX_RANDOM_NUMBER => 100_00;
    Readonly my $PROCESSING_TIME   => 3;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Unpack parameters
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};

    ## FILEHANDLES
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Generate a random integer between 0-10,000.
    my $random_integer = int rand $MAX_RANDOM_NUMBER;

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($file_path) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $family_id,
            job_id_href           => $job_id_href,
            program_directory     => $outaligner_dir,
            program_name          => $program_name,
            process_time          => $PROCESSING_TIME,
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
            parameter_href          => $parameter_href,
            program_name            => $program_name,
            random_integer          => $random_integer,
            sample_info_href        => $sample_info_href,
        }
    );

    if ( $parameter_href->{rtg_vcfeval_reference_genome}{build_file} == 1 ) {

        $log->warn( q{Will try to create required }
              . $human_genome_reference
              . q{ sdf files before executing }
              . $program_name );

        say {$FILEHANDLE} q{## Building SDF dir files};
        ## Get parameters
        my $sdf_directory_tmp =
            $active_parameter_href->{rtg_vcfeval_reference_genome}
          . $UNDERSCORE
          . $random_integer;
        rtg_format(
            {
                FILEHANDLE            => $FILEHANDLE,
                input_format          => q{fasta},
                reference_genome_path => $human_genome_reference,
                sdf_output_directory  => $sdf_directory_tmp,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

      PREREQ:
        foreach my $suffix ( @{$rtg_directory_suffixes_ref} ) {

            my $intended_file_path =
              $active_parameter_href->{rtg_vcfeval_reference_genome} . $suffix;

            ## Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
            check_exist_and_move_file(
                {
                    FILEHANDLE          => $FILEHANDLE,
                    intended_file_path  => $intended_file_path,
                    temporary_file_path => $sdf_directory_tmp,
                }
            );
        }

        ## Ensure that this subrutine is only executed once
        $parameter_href->{rtg_vcfeval_reference_genome}{build_file} = 0;
    }

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

    if ( $mip_program_mode == 1 ) {

        slurm_submit_job_no_dependency_add_to_samples(
            {
                job_id_href      => $job_id_href,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                family_id        => $family_id,
                path             => $job_id_chain,
                sbatch_file_name => $file_path,
                log              => $log,
            }
        );
    }
    return;
}

1;
