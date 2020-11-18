package MIP::Recipes::Analysis::Mip_vercollect;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

## MIPs lib/
use MIP::Constants qw{ %CONTAINER_CMD $LOG_NAME $NEWLINE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.12;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_mip_vercollect };
}

sub analysis_mip_vercollect {

## Function : Collect executable versions for this analysis run.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $sample_info_href        => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
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
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Io::Write qw{ write_to_file };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Mip qw{ mip_vercollect };
    use MIP::Sample_info qw{ set_file_path_to_store set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    my $job_id_chain = get_recipe_attributes(
        {
            attribute      => q{chain},
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
        }
    );
    my $recipe_mode     = $active_parameter_href->{$recipe_name};
    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    my %io = parse_io_outfiles(
        {
            chain_id               => $job_id_chain,
            id                     => $case_id,
            file_info_href         => $file_info_href,
            file_name_prefixes_ref => [$case_id],
            outdata_dir            => $active_parameter_href->{outdata_dir},
            parameter_href         => $parameter_href,
            recipe_name            => $recipe_name,
        }
    );

    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_path        = $io{out}{file_path};

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ($recipe_file_path) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $recipe_resource{core_number},
            directory_id          => $case_id,
            filehandle            => $filehandle,
            job_id_href           => $job_id_href,
            memory_allocation     => $recipe_resource{memory},
            process_time          => $recipe_resource{time},
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
            set_pipefail          => 0,
        }
    );

    use MIP::Environment::Executable qw{ write_binaries_versions };
    use MIP::Program::Gnu::Coreutils qw{ gnu_rm };

    ### SHELL:

    say {$filehandle} q{## Remove potential previous mip version collect file};
    gnu_rm(
        {
            filehandle  => $filehandle,
            force       => 1,
            infile_path => $outfile_path,
        }
    );
    say {$filehandle} $NEWLINE;
    ## Write executable versions
    write_binaries_versions(
        {
            binary_info_href => \%CONTAINER_CMD,
            filehandle       => $filehandle,
            outfile_path     => $outfile_path,
        }
    );

#    my $binary_path_file_path = $outfile_path_prefix . $UNDERSCORE . q{binary_paths.yaml};

    ## Writes a YAML hash to file
    #    write_to_file(
    #        {
    #            data_href => \%CONTAINER_CMD,
    #            format    => q{yaml},
    #            path      => $binary_path_file_path,
    #        }
    #    );
    #    $log->info( q{Wrote: } . $binary_path_file_path );

    #    my $log_file_path = $outfile_path_prefix . $UNDERSCORE . q{vercollect.log};

    #    mip_vercollect(
    #       {
    #           filehandle    => $filehandle,
    #           infile_path   => $binary_path_file_path,
    #           log_file_path => $log_file_path,
    #           outfile_path  => $outfile_path,
    #       }
    #   );
    #   say {$filehandle} $NEWLINE;

    close $filehandle or $log->logcroak(q{Could not close filehandle});

    if ( $recipe_mode == 1 ) {

        set_recipe_outfile_in_sample_info(
            {
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        set_file_path_to_store(
            {
                format           => q{meta},
                id               => $case_id,
                path             => $outfile_path,
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command         => $profile_base_command,
                dependency_method    => q{add_to_all},
                job_dependency_type  => q{afterok},
                job_id_chain         => $job_id_chain,
                job_id_href          => $job_id_href,
                job_reservation_name => $active_parameter_href->{job_reservation_name},
                log                  => $log,
                recipe_file_path     => $recipe_file_path,
                submission_profile   => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

1;
