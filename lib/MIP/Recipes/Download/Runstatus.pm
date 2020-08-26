package MIP::Recipes::Download::Runstatus;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME $NEWLINE $SPACE $TAB $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ download_runstatus };

}

sub download_runstatus {

## Function : Download run status
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this download hash {REF}
##          : $job_id_href           => The job_id hash {REF}
##          : $profile_base_command  => Submission profile base command
##          : $recipe_name           => Recipe name
##          : $temp_directory        => Temporary directory for recipe

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $job_id_href;
    my $recipe_name;

    ## Default(s)
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
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
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
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::Parameter qw{ get_recipe_resources };
    use MIP::Language::Shell qw{ check_mip_process_paths };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Processmanagement::Processes qw{ submit_recipe };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    my $reference_dir = $active_parameter_href->{reference_dir};

    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Set recipe mode
    my $recipe_mode = $active_parameter_href->{$recipe_name};

    ## Filehandle(s)
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $recipe_resource{core_number},
            directory_id                    => q{mip_download},
            filehandle                      => $filehandle,
            job_id_href                     => $job_id_href,
            memory_allocation               => $recipe_resource{memory},
            outdata_dir                     => $reference_dir,
            outscript_dir                   => $reference_dir,
            process_time                    => $recipe_resource{time},
            recipe_data_directory_path      => $active_parameter_href->{reference_dir},
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => $recipe_resource{load_env_ref},
        }
    );

    ### SHELL:

    say {$filehandle} q{## } . $recipe_name;

    ## Set status flagg
    say {$filehandle} q?STATUS="0"?;

    check_mip_process_paths(
        {
            filehandle => $filehandle,
            paths_ref  => $active_parameter_href->{runstatus_paths},
        }
    );

    say {$filehandle} q?if [ "$STATUS" -eq 1 ]; then?;
    say {$filehandle} $TAB . q?exit 1?;
    say {$filehandle} q?fi?, $NEWLINE;

    ## Close filehandleS
    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe_mode == 1 ) {

        submit_recipe(
            {
                base_command        => $profile_base_command,
                dependency_method   => q{add_to_all},
                job_dependency_type => q{afterany},
                job_id_chain        => q{PAN},
                job_id_href         => $job_id_href,
                log                 => $log,
                recipe_file_path    => $recipe_file_path,
                submission_profile  => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

1;
