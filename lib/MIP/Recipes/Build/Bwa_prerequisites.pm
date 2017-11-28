package MIP::Recipes::Build::Bwa_prerequisites;

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
    our @EXPORT_OK = qw{ build_bwa_prerequisites };

}

##Constants
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub build_bwa_prerequisites {

## Function : Creates the Bwa prerequisites
## Returns  :
## Arguments: $parameter_href                       => Parameter hash {REF}
##          : $active_parameter_href                => Active parameters for this analysis hash {REF}
##          : $sample_info_href                     => Info on samples and family hash {REF}
##          : $file_info_href                       => File info hash {REF}
##          : $infile_lane_prefix_href              => Infile(s) without the ".ending" {REF}
##          : $job_id_href                          => Job id hash {REF}
##          : $bwa_build_reference_file_endings_ref => The bwa reference associated file endings {REF}
##          : $family_id                            => Family ID {REF}
##          : $outaligner_dir_ref                   => The outaligner_dir used in the analysis
##          : $program_name                         => Program name
##          : $family_id                            => Family id
##          : $temp_directory                       => Temporary directory
##          : $outaligner_dir                       => Outaligner_dir used in the analysis
##          : $human_genome_reference               => Human genome reference

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $bwa_build_reference_file_endings_ref;
    my $program_name;

    ## Default(s)
    my $family_id;
    my $temp_directory;
    my $outaligner_dir;
    my $human_genome_reference;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        bwa_build_reference_file_endings_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$bwa_build_reference_file_endings_ref
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory,
        },
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir,
        },
        human_genome_reference => {
            default =>
              $arg_href->{active_parameter_href}{human_genome_reference},
            strict_type => 1,
            store       => \$human_genome_reference
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Alignment::Bwa qw{ bwa_index };
    use MIP::Language::Shell qw{ check_exist_and_move_file };
    use MIP::Recipes::Build::Human_genome_prerequisites
      qw{ build_human_genome_prerequisites };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_no_dependency_add_to_samples };

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
            job_id_href           => $job_id_href,
            FILEHANDLE            => $FILEHANDLE,
            directory_id          => $family_id,
            program_name          => $program_name,
            program_directory     => $outaligner_dir,
            process_time          => $PROCESSING_TIME,
        }
    );

    build_human_genome_prerequisites(
        {
            parameter_href          => $parameter_href,
            active_parameter_href   => $active_parameter_href,
            sample_info_href        => $sample_info_href,
            file_info_href          => $file_info_href,
            infile_lane_prefix_href => $infile_lane_prefix_href,
            job_id_href             => $job_id_href,
            program_name            => $program_name,
            FILEHANDLE              => $FILEHANDLE,
            log                     => $log,
            random_integer          => $random_integer,
        }
    );

    if ( $parameter_href->{bwa_build_reference}{build_file} == 1 ) {

        $log->warn( q{Will try to create required }
              . $human_genome_reference
              . q{ index files before executing }
              . $program_name );

        say {$FILEHANDLE} q{## Building BWA index};
        ## Get parameters
        my $prefix = $human_genome_reference . $UNDERSCORE . $random_integer;
        bwa_index(
            {
                prefix                 => $prefix,
                construction_algorithm => q{bwtsw},
                reference_genome       => $human_genome_reference,
                FILEHANDLE             => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

      PREREQ_FILE:
        foreach my $file ( @{$bwa_build_reference_file_endings_ref} ) {

            my $intended_file_path = $human_genome_reference . $file;
            my $temporary_file_path =
              $human_genome_reference . $UNDERSCORE . $random_integer . $file;

            ## Checks if a file exists and moves the file in place if file is lacking or has a size of 0 bytes.
            check_exist_and_move_file(
                {
                    FILEHANDLE          => $FILEHANDLE,
                    intended_file_path  => $intended_file_path,
                    temporary_file_path => $temporary_file_path,
                }
            );
        }

        ## Ensure that this subrutine is only executed once
        $parameter_href->{bwa_build_reference}{build_file} = 0;
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
