package MIP::Environment::Container;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

## MIPs lib/
use MIP::Constants
  qw{ $COLON $COMMA @CONTAINER_BIND_PATHS $CONTAINER_MANAGER $DOUBLE_QUOTE $EQUALS $SEMICOLON $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ get_recipe_executable_bind_path install_containers parse_containers parse_container_bind_paths parse_container_uri run_container set_executable_container_cmd };
}

sub get_recipe_executable_bind_path {

## Function : Get link between recipe and executables and set recipe_binds_path to executable
## Returns  : %recipe_executable_bind_path
## Arguments: $active_parameter_href => The active parameters for this analysis hash {REF}
##          : $parameter_href        => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Active_parameter qw{add_recipe_bind_paths};
    use MIP::Constants qw{@CONTAINER_BIND_PATHS};
    use MIP::Parameter qw{get_cache get_parameter_attribute};

    my %recipe_executable_bind_path;
    my @recipes = get_cache(
        {
            parameter_href => $parameter_href,
            parameter_name => q{recipe},
        }
    );

    foreach my $recipe_name (@recipes) {

        my @recipe_executables = get_parameter_attribute(
            {
                attribute      => q{program_executables},
                parameter_href => $parameter_href,
                parameter_name => $recipe_name,
            }
        );
        foreach my $executable (@recipe_executables) {

            my @export_bind_paths = @CONTAINER_BIND_PATHS;
            add_recipe_bind_paths(
                {
                    active_parameter_href => $active_parameter_href,
                    export_bind_paths_ref => \@export_bind_paths,
                    recipe_name           => $recipe_name,
                }
            );
            $recipe_executable_bind_path{$executable} = [@export_bind_paths];
        }
    }
    return %recipe_executable_bind_path;
}

sub install_containers {

## Function : Setup containers to use with docker or singularity
## Returns  :
## Arguments: $active_parameter_href => Active parameter hash {REF}
##          : $conda_env_path        => Path to conda environment
##          : $container_href        => Hash with container {REF}
##          : $container_manager     => Container manager
##          : $filehandle            => Filehandle

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $conda_env_path;
    my $container_href;
    my $container_manager;
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
        container_manager => {
            allow       => [qw{ docker singularity }],
            required    => 1,
            store       => \$container_manager,
            strict_type => 1,
        },
        filehandle => {
            required => 1,
            store    => \$filehandle,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Recipes::Install::Singularity qw{ install_singularity_containers };
    use MIP::Recipes::Install::Docker qw{ install_docker_containers };

    ## Return if no containers
    return if not keys %{$container_href};

    my %container_api = (
        docker => {
            arg_href => {
                active_parameter_href => $active_parameter_href,
                conda_env_path        => $conda_env_path,
                container_href        => $container_href,
                filehandle            => $filehandle,
            },
            method => \&install_docker_containers,
        },
        singularity => {
            arg_href => {
                active_parameter_href => $active_parameter_href,
                conda_env_path        => $conda_env_path,
                container_href        => $container_href,
                filehandle            => $filehandle,
            },
            method => \&install_singularity_containers,
        },
    );

    $container_api{$container_manager}{method}
      ->( { %{ $container_api{$container_manager}{arg_href} } } );

    return 1;
}

sub parse_containers {

## Function : Parse containers to set executable command based on current container manager
## Returns  :
## Arguments: $active_parameter_href => The active parameters for this analysis hash {REF}
##          : $parameter_href        => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Config qw{ get_install_containers };
    use MIP::Constants qw{ set_container_cmd };
    use MIP::Update::Parameters qw{ update_with_dynamic_config_parameters };

    %{ $active_parameter_href->{container} } =
      get_install_containers(
        { install_config_file => $active_parameter_href->{install_config_file}, } );

    my %dynamic_parameter = ( reference_dir => $active_parameter_href->{reference_dir}, );
    update_with_dynamic_config_parameters(
        {
            active_parameter_href  => $active_parameter_href->{container},
            dynamic_parameter_href => \%dynamic_parameter,
        }
    );

    my %container_cmd = set_executable_container_cmd(
        {
            active_parameter_href => $active_parameter_href,
            container_href        => $active_parameter_href->{container},
            container_manager     => $active_parameter_href->{container_manager},
            parameter_href        => $parameter_href,
        }
    );

    set_container_cmd( { container_cmd_href => \%container_cmd, } );

    return 1;
}

sub parse_container_bind_paths {

## Function : Parse container bind paths and add export command to array
## Returns  : $singularity_bind
## Arguments: $active_parameter_href       => The active parameters for this analysis hash {REF}
##          : $package_name                => Package name
##          : $source_environment_cmds_ref => Array with source environment commands {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $package_name;
    my $source_environment_cmds_ref;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        package_name => {
            defined     => 1,
            required    => 1,
            store       => \$package_name,
            strict_type => 1,
        },
        source_environment_cmds_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$source_environment_cmds_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Active_parameter qw{ add_recipe_bind_paths };
    use MIP::Environment::Path
      qw{ build_docker_bind_path_var build_singularity_bind_path_var };

    return if not $CONTAINER_MANAGER;

    my @export_bind_paths = @CONTAINER_BIND_PATHS;
    add_recipe_bind_paths(
        {
            active_parameter_href => $active_parameter_href,
            export_bind_paths_ref => \@export_bind_paths,
            recipe_name           => $package_name,
        }
    );

    my %container = (
        docker      => \&build_docker_bind_path_var,
        singularity => \&build_singularity_bind_path_var,
    );

    my $container_bind_var = $container{$CONTAINER_MANAGER}->(
        {
            bind_paths_ref => \@export_bind_paths,
        }
    );

    push @{$source_environment_cmds_ref}, $container_bind_var;

    return $container_bind_var;
}

sub parse_container_uri {

## Function : Parse container uri for selected container manager
## Returns  :
## Arguments: $container_manager => Container manager
##          : $uri_ref           => Container uri {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $container_manager;
    my $uri_ref;

    my $tmpl = {
        container_manager => {
            allow       => [qw{ docker singularity }],
            required    => 1,
            store       => \$container_manager,
            strict_type => 1,
        },
        uri_ref => {
            defined  => 1,
            required => 1,
            store    => \$uri_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return if $container_manager eq q{docker};

    if ( ${$uri_ref} =~ /\A quay|docker[.]io /xms ) {

        ${$uri_ref} = q{docker://} . ${$uri_ref};
    }

    return;
}

sub run_container {

## Function : Run a docker container or exec a singularity image
## Returns  : @commands
## Arguments: $bind_paths_ref         => Bind host directory to container {REF}
##          : $container_cmds_ref     => Cmds to be executed in container {REF}
##          : $container_manager      => Container manager
##          : $container_path         => Path to container
##          : $entrypoint             => Override container entrypoint
##          : $filehandle             => Filehandle to write to
##          : $image                  => Image to run
##          : $remove                 => Remove stopped container
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdinfile_path         => Stdinfile path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bind_paths_ref;
    my $container_cmds_ref;
    my $container_manager;
    my $container_path;
    my $filehandle;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdinfile_path;
    my $stdoutfile_path;

    ## Default(s)
    my $remove;

    my $tmpl = {
        bind_paths_ref => {
            default     => [],
            store       => \$bind_paths_ref,
            strict_type => 1,
        },
        container_cmds_ref => {
            default     => [],
            store       => \$container_cmds_ref,
            strict_type => 1,
        },
        container_manager => {
            allow       => [qw{ docker singularity }],
            required    => 1,
            store       => \$container_manager,
            strict_type => 1,
        },
        container_path => {
            defined     => 1,
            required    => 1,
            store       => \$container_path,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        remove => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$remove,
            strict_type => 1,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        stderrfile_path_append => {
            store       => \$stderrfile_path_append,
            strict_type => 1,
        },
        stdinfile_path  => { store => \$stdinfile_path, strict_type => 1, },
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Singularity qw{ singularity_exec };
    use MIP::Program::Docker qw{ docker_run };

  BIND_PATH:
    foreach my $bind_path ( @{$bind_paths_ref} ) {

        # Skip if mount path in image already is specified
        next BIND_PATH if $bind_path =~ m/[:]/xms;

        next BIND_PATH if $bind_path =~ m/\A \$/xms;

        $bind_path .= $COLON . $bind_path;
    }

    my %container_api = (
        docker => {
            arg_href => {
                bind_paths_ref         => $bind_paths_ref,
                container_cmds_ref     => $container_cmds_ref,
                filehandle             => $filehandle,
                image                  => $container_path,
                remove                 => $remove,
                stderrfile_path        => $stderrfile_path,
                stderrfile_path_append => $stderrfile_path_append,
                stdinfile_path         => $stdinfile_path,
                stdoutfile_path        => $stdoutfile_path,
            },
            method => \&docker_run,
        },
        singularity => {
            arg_href => {
                bind_paths_ref                 => $bind_paths_ref,
                filehandle                     => $filehandle,
                image                          => $container_path,
                singularity_container_cmds_ref => $container_cmds_ref,
                stderrfile_path                => $stdinfile_path,
                stderrfile_path_append         => $stderrfile_path_append,
                stdoutfile_path                => $stdoutfile_path,
            },
            method => \&singularity_exec,
        },
    );

    my @commands = $container_api{$container_manager}{method}
      ->( { %{ $container_api{$container_manager}{arg_href} } } );

    return @commands;
}

sub set_executable_container_cmd {

## Function : Set executable command depending on container manager
## Returns  :
## Arguments: $active_parameter_href => The active parameters for this analysis hash {REF}
##          : $container_href        => Containers hash {REF}
##          : $container_manager     => Container manager
##          : $parameter_href        => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $container_href;
    my $container_manager;
    my $parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        container_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$container_href,
            strict_type => 1,
        },
        container_manager => {
            allow       => [qw{ docker singularity }],
            required    => 1,
            store       => \$container_manager,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Constants qw{ @CONTAINER_BIND_PATHS };
    use Data::Diver qw{ Dive };
    use MIP::Program::Singularity qw{ singularity_exec };
    use MIP::Program::Docker qw{ docker_run };

    my %container_cmd;
    my @container_constant_bind_path = @CONTAINER_BIND_PATHS;

    my %recipe_executable_bind_path = get_recipe_executable_bind_path(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
        }
    );

  CONTAINER_NAME:
    foreach my $container_name ( keys %{$container_href} ) {

        parse_container_uri(
            {
                container_manager => $container_manager,
                uri_ref           => \$container_href->{$container_name}{uri},
            }
        );

      EXECUTABLE:
        while ( my ( $executable_name, $executable_path ) =
            each %{ $container_href->{$container_name}{executable} } )
        {

            ## Installation specific bind paths
            if (
                Dive(
                    $container_href, ( $container_name, q{bind_path}, $executable_name, )
                )
              )
            {
                push @{ $recipe_executable_bind_path{$executable_name} },
                  $container_href->{$container_name}{bind_path}{$executable_name};
            }
            my @bind_paths =
              exists $recipe_executable_bind_path{$executable_name}
              ? @{ $recipe_executable_bind_path{$executable_name} }
              : @container_constant_bind_path;

            my @cmds = run_container(
                {
                    bind_paths_ref    => \@bind_paths,
                    container_manager => $container_manager,
                    container_path    => $container_href->{$container_name}{uri},
                }
            );

            ## Do not add anything to @cmds
            if ( $executable_path and $executable_path eq q{no_executable_in_image} ) {
            }
            elsif ($executable_path) {

                push @cmds, $executable_path;
            }
            else {

                push @cmds, $executable_name;
            }
            $container_cmd{$executable_name} = join $SPACE, @cmds;
        }
    }
    return %container_cmd;
}

1;
