package MIP::Recipes::Install::Docker;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile devnull };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants
  qw{ $AT $BACKWARD_SLASH $CLOSE_BRACE $DOLLAR_SIGN $DOUBLE_QUOTE $COLON $EMPTY_STR $LOG_NAME $OPEN_BRACE $NEWLINE $SINGLE_QUOTE $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_docker_containers };
}

sub install_docker_containers {

## Function : Pull container from singularity hub or docker hub
## Returns  :
## Arguments: $active_parameter_href => Active parameter hash {REF}
##          : $conda_env_path        => Path to conda environment
##          : $container_href        => Hash with container {REF}
##          : $filehandle            => Filehandle

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $conda_env_path;
    my $container_href;
    my $filehandle;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        conda_env_path => {
            required    => 1,
            store       => \$conda_env_path,
            strict_type => 1,
        },
        container_href => {
            default     => {},
            required    => 1,
            store       => \$container_href,
            strict_type => 1,
        },
        filehandle => {
            required => 1,
            store    => \$filehandle,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Environment::Container qw{ parse_container_uri };
    use MIP::Recipes::Install::Cadd qw{ install_cadd };
    use MIP::Recipes::Install::Vep qw{ install_vep };
    use MIP::Set::Parameter qw{ set_container_bind_paths };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    say {$filehandle} q{## Setup containers for docker};

    ## Containers requiring something extra
    my %finish_container_installation = (
        cadd => \&install_cadd,
        vep  => \&install_vep,
    );

    my $mip_bind_path =
      $BACKWARD_SLASH . $DOLLAR_SIGN . $OPEN_BRACE . q{MIP_BIND} . $CLOSE_BRACE;

  CONTAINER:
    while ( my ( $image_name, $image_href ) = each %{$container_href} ) {

        $log->info(
            q{Writing instructions for running image} . $COLON . $SPACE . $image_name );

        say {$filehandle} qq{## Setting up $image_name image};

        set_container_bind_paths(
            {
                bind_paths_ref => [$mip_bind_path],
                container_href => $image_href,
            }
        );

        ## Finishing touches for certain containers
        if ( $finish_container_installation{$image_name} ) {

            $finish_container_installation{$image_name}->(
                {
                    active_parameter_href => $active_parameter_href,
                    container_href        => $image_href,
                    container_path        => $image_href->{uri},
                    filehandle            => $filehandle,
                }
            );
        }

        ## Make available as exeuctable in bin
        setup_docker_executable(
            {
                conda_env_path         => $conda_env_path,
                executable_href        => $image_href->{executable},
                filehandle             => $filehandle,
                image_address          => $image_href->{uri},
                program_bind_paths_ref => $image_href->{program_bind_paths},
            }
        );
        print {$filehandle} $NEWLINE;

    }
    return 1;
}

sub setup_docker_executable {

## Function : Make docker executable available in conda bin
## Returns  :
## Arguments: $conda_env_path         => Path to conda environment
##          : $executable_href        => Hash with executables and their path in the container (if not in PATH) {REF}
##          : $filehandle             => Filehandle
##          : $image_address          => Address to docker image
##          : $program_bind_paths_ref => Extra static bind paths

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_env_path;
    my $executable_href;
    my $filehandle;
    my $image_address;
    my $program_bind_paths_ref;

    my $tmpl = {
        conda_env_path => {
            required    => 1,
            store       => \$conda_env_path,
            strict_type => 1,
        },
        executable_href => {
            default     => {},
            required    => 1,
            store       => \$executable_href,
            strict_type => 1,
        },
        filehandle => {
            required => 1,
            store    => \$filehandle,
        },
        image_address => {
            required    => 1,
            store       => \$image_address,
            strict_type => 1,
        },
        program_bind_paths_ref => {
            default     => $arg_href->{program_bind_paths_ref} ||= [],
            store       => \$program_bind_paths_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Language::Shell qw{ build_shebang };
    use MIP::Program::Docker qw{ docker_run };
    use MIP::Program::Gnu::Bash qw{ gnu_unset };
    use MIP::Program::Gnu::Coreutils qw{ gnu_chmod gnu_echo };

    my @shebang = build_shebang(
        {
            bash_bin_path => catfile( dirname( dirname( devnull() ) ), qw{ bin bash } ),
        }
    );

    my $bash_command =
        $BACKWARD_SLASH
      . $DOUBLE_QUOTE
      . $BACKWARD_SLASH
      . $DOLLAR_SIGN
      . $AT
      . $BACKWARD_SLASH
      . $DOUBLE_QUOTE;

  EXECUTABLE:
    foreach my $executable ( keys %{$executable_href} ) {

        my $container_executable = $executable_href->{$executable} ||= $executable;

        if ( $container_executable eq q{no_executable_in_image} ) {
            $container_executable = $EMPTY_STR;
        }
        my @docker_cmds = docker_run(
            {
                bind_paths_ref     => $program_bind_paths_ref,
                container_cmds_ref => [$container_executable],
                entrypoint         => $BACKWARD_SLASH
                  . $DOUBLE_QUOTE
                  . $BACKWARD_SLASH
                  . $DOUBLE_QUOTE,
                image  => $image_address,
                remove => 1,
            }
        );
        push @docker_cmds, $bash_command;

        my $proxy_executable_path = catfile( $conda_env_path, q{bin}, $executable );

        gnu_echo(
            {
                filehandle     => $filehandle,
                outfile_path   => $proxy_executable_path,
                string_wrapper => $SINGLE_QUOTE,
                strings_ref    => \@shebang,
            }
        );
        print {$filehandle} $NEWLINE;

        gnu_echo(
            {
                filehandle             => $filehandle,
                no_trailing_newline    => 1,
                stdoutfile_path_append => $proxy_executable_path,
                strings_ref            => [ join $SPACE, @docker_cmds ],
            }
        );
        print {$filehandle} $NEWLINE;
        gnu_chmod(
            {
                file_path  => $proxy_executable_path,
                filehandle => $filehandle,
                permission => q{a+x},
            }
        );
        print {$filehandle} $NEWLINE;

        ## Unset VEP_VERSION
        if ( $executable eq q{vep} ) {

            gnu_unset(
                {
                    bash_variable => q{VEP_VERSION},
                    filehandle    => $filehandle,
                }
            );
            print {$filehandle} $NEWLINE;
        }
    }
    return 1;
}

1;
