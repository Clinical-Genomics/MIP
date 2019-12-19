package MIP::Recipes::Pipeline::Install_rd_dna;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use List::Util qw{ any };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## MIPs lib/
use MIP::Constants qw{
  $COLON
  $DOT
  $LOG_NAME
  $NEWLINE
  $SINGLE_QUOTE
  $SPACE
};
use MIP::Script::Setup_script qw{ setup_install_script };
use MIP::Set::Parameter qw{ set_programs_for_installation };

## Recipes
use MIP::Recipes::Install::Conda qw{ install_conda_packages };
use MIP::Recipes::Install::Mip_scripts qw{ install_mip_scripts };
use MIP::Recipes::Install::Pip qw{ install_pip_packages };
use MIP::Recipes::Install::Post_installation qw{check_mip_installation };
use MIP::Recipes::Install::Singularity qw{ install_singularity_containers };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.17;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ pipeline_install_rd_dna };
}

sub pipeline_install_rd_dna {

## Function : Install recipes for rd_dna pipeline
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this download hash {REF}
##          : $quiet                 => Be quiet
##          : $verbose               => Verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;

    ## Default(s)
    my $quiet;
    my $verbose;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        quiet => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$quiet,
            strict_type => 1,
        },
        verbose => {
            default     => $arg_href->{active_parameter_href}{verbose},
            store       => \$verbose,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object now that log_file has been set
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Installation instruction file
    my $file_name_path = catfile( cwd(), q{mip.sh} );

    open my $filehandle, q{>}, $file_name_path
      or $log->logcroak( q{Cannot write to}
          . $SPACE
          . $SINGLE_QUOTE
          . $file_name_path
          . $SINGLE_QUOTE
          . $SPACE
          . $COLON
          . $OS_ERROR
          . $NEWLINE );

    ## Create bash file for writing install instructions
    setup_install_script(
        {
            active_parameter_href => $active_parameter_href,
            file_name             => $file_name_path,
            filehandle            => $filehandle,
            invoke_login_shell    => $active_parameter_href->{sbatch_mode},
            log                   => $log,
            remove_dir            => catfile( cwd(), $DOT . q{MIP} ),
            sbatch_mode           => $active_parameter_href->{sbatch_mode},
            set_errexit           => $active_parameter_href->{bash_set_errexit},
            set_nounset           => $active_parameter_href->{bash_set_nounset},
        }
    );

    $log->info( q{Writing install instructions to:} . $SPACE . $file_name_path );

    ## Source conda
    if ( not $active_parameter_href->{sbatch_mode} ) {
        say {$filehandle} q{## Source conda};
        say {$filehandle} q{source}
          . $SPACE
          . catfile( $active_parameter_href->{conda_path}, qw{ etc profile.d conda.sh } )
          . $NEWLINE;
    }

    my $env_name = $active_parameter_href->{environment_name};

    $log->info( q{Installing into environment: } . $env_name );

    ## Process input parameters to get a correct combination of programs that are to be installed
    set_programs_for_installation(
        {
            active_parameter_href => $active_parameter_href,
        }
    );

    ## Installing Conda packages
    install_conda_packages(
        {
            conda_env           => $active_parameter_href->{environment_name},
            conda_env_path      => $active_parameter_href->{conda_prefix_path},
            conda_no_update_dep => $active_parameter_href->{conda_no_update_dep},
            conda_packages_href => $active_parameter_href->{conda},
            filehandle          => $filehandle,
            quiet               => $active_parameter_href->{quiet},
            verbose             => $active_parameter_href->{verbose},
        }
    );

    ## Install PIP packages
    install_pip_packages(
        {
            conda_env         => $active_parameter_href->{environment_name},
            filehandle        => $filehandle,
            pip_packages_href => $active_parameter_href->{pip},
            quiet             => $active_parameter_href->{quiet},
        }
    );

    ## Pull and link containers
    install_singularity_containers(
        {
            active_parameter_href => $active_parameter_href,
            conda_env_path        => $active_parameter_href->{conda_prefix_path},
            container_href        => $active_parameter_href->{singularity},
            filehandle            => $filehandle,
            quiet                 => $active_parameter_href->{quiet},
            verbose               => $active_parameter_href->{verbose},
        }
    );

    ### Install shell programs
    ## Create dispatch table for shell installation subs
    my %shell_subs = ( mip_scripts => \&install_mip_scripts, );

    ## Launch shell installation subroutines
  SHELL_PROGRAM:
    for my $shell_program ( keys %{ $active_parameter_href->{shell} } ) {

        $shell_subs{$shell_program}->(
            {
                conda_environment => $active_parameter_href->{environment_name},
                conda_prefix_path => $active_parameter_href->{conda_prefix_path},
                filehandle        => $filehandle,
                program_parameters_href =>
                  \%{ $active_parameter_href->{shell}{$shell_program} },
                quiet   => $active_parameter_href->{quiet},
                verbose => $active_parameter_href->{verbose},
            }
        );
    }

    ## Write tests
    check_mip_installation(
        {
            active_parameter_href => $active_parameter_href,
            filehandle            => $filehandle,
        }
    );

    $log->info(q{Finished writing installation instructions for MIP});

    close $filehandle or $log->logcroak(q{Could not close filehandle});
    return;
}

1;
