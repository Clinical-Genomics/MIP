package MIP::Environment::Container;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ fileparse };
use File::Path qw{ make_path };
use File::Spec::Functions qw{ catfile catdir };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

## MIPs lib/
use MIP::Constants
  qw{ $COLON $COMMA @CONTAINER_BIND_PATHS $CONTAINER_MANAGER $DOUBLE_QUOTE $EMPTY_STR $EQUALS $LOG_NAME $SEMICOLON $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      build_container_cmd
      check_installed_containers
      get_recipe_executable_bind_path
      parse_container_config
      parse_container_path
      parse_container_uri
      parse_containers
      pull_container
      run_container
      set_executable_container_cmd
    };
}

sub build_container_cmd {

    ## Function : Build executable command depending on container manager
    ## Returns  :
    ## Arguments: $active_parameter_href            => The active parameters for this analysis hash {REF}
    ##          : $container_href                   => Containers hash {REF}
    ##          : $container_manager                => Container manager
    ##          : $recipe_executable_bind_path_href => Recipe bind path hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $container_href;
    my $container_manager;
    my $recipe_executable_bind_path_href;

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
            allow       => [qw{docker singularity}],
            required    => 1,
            store       => \$container_manager,
            strict_type => 1,
        },
        recipe_executable_bind_path_href => {
            default     => {},
            defined     => 1,
            store       => \$recipe_executable_bind_path_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Data::Diver qw{ Dive };
    use List::MoreUtils qw { any };

    my @container_constant_bind_path = @CONTAINER_BIND_PATHS;
    my %container_cmd;

    my @gpu_executables =
      exists $active_parameter_href->{gpu_capable_executables}
      ? @{ $active_parameter_href->{gpu_capable_executables} }
      : [];

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
            if ( Dive( $container_href, ( $container_name, q{bind_path}, $executable_name, ) ) ) {
                push @{ $recipe_executable_bind_path_href->{$executable_name} },
                  $container_href->{$container_name}{bind_path}{$executable_name};
            }
            my @bind_paths =
              exists $recipe_executable_bind_path_href->{$executable_name}
              ? @{ $recipe_executable_bind_path_href->{$executable_name} }
              : @container_constant_bind_path;

            my $gpu_switch;
            if ( any { $_ eq $executable_name } @gpu_executables ) {
                $gpu_switch = $container_href->{$container_name}{gpu_support} ? 1 : 0;
            }

            ## Isolate singularity container
            my $no_home   = 0;
            my $clean_env = 0;
            if ( any { $_ eq $executable_name }
                @{ $active_parameter_href->{isolate_singularity_containers} } )
            {
                $no_home   = 1;
                $clean_env = 1;
            }

            my @cmds = run_container(
                {
                    bind_paths_ref    => \@bind_paths,
                    executable_name   => $executable_name,
                    container_manager => $container_manager,
                    container_path    => $container_href->{$container_name}{uri},
                    gpu_switch        => $gpu_switch,
                    no_home           => $no_home,
                    clean_env         => $clean_env,
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
    use MIP::Environment::Path qw{ reduce_dir_paths };
    use MIP::Language::Shell qw{ quote_bash_variable };
    use MIP::Parameter qw{get_cache get_parameter_attribute};

    my $temp_directory_quoted = quote_bash_variable(
        { string_with_variable_to_quote => $active_parameter_href->{temp_directory}, } );
    my $xdg_runtime_dir =
      $temp_directory_quoted . $COLON . catfile( $EMPTY_STR, qw{ run user }, q{$(id -u)} );

    my %recipe_executable_bind_path;
    my @recipes = get_cache(
        {
            parameter_href => $parameter_href,
            parameter_name => q{recipe},
        }
    );

  RECIPE:
    foreach my $recipe_name (@recipes) {

        my @recipe_executables = get_parameter_attribute(
            {
                attribute      => q{program_executables},
                parameter_href => $parameter_href,
                parameter_name => $recipe_name,
            }
        );
      RECIPE_EXECUTABLE:
        foreach my $executable (@recipe_executables) {

            my @export_bind_paths = @CONTAINER_BIND_PATHS;
            add_recipe_bind_paths(
                {
                    active_parameter_href => $active_parameter_href,
                    export_bind_paths_ref => \@export_bind_paths,
                    recipe_name           => $recipe_name,
                }
            );

            push @{ $recipe_executable_bind_path{$executable} }, @export_bind_paths;
        }
    }

  EXECUTABLE:
    foreach my $executable ( keys %recipe_executable_bind_path ) {

        ## Special case for xdg_runtime_dir, which always should be added
        push @{ $recipe_executable_bind_path{$executable} }, $xdg_runtime_dir;

        $recipe_executable_bind_path{$executable} = [
            reduce_dir_paths(
                {
                    dir_paths_ref => $recipe_executable_bind_path{$executable},
                }
            )
        ];
    }
    return %recipe_executable_bind_path;
}

sub check_installed_containers {

## Function : Parse containers to set executable command based on current container manager
## Returns  :
## Arguments: $container_href => Map of containers {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $container_href;

    my $tmpl = {
        container_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$container_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Path qw{ check_filesystem_objects_existance };
    use MIP::Validate::Data qw{ %constraint };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my @error_messages;

  CONTAINER:
    while ( my ( $container, $container_param_href ) = each %{$container_href} ) {

        ## Only run check for .sif file
        next CONTAINER if ( not $constraint{is_sif}->( $container_param_href->{uri} ) );

        my ( $is_ok, $error_message ) = check_filesystem_objects_existance(
            {
                object_name    => $container_param_href->{uri},
                object_type    => q{executable_file},
                parameter_name => $container,
            }
        );

        next CONTAINER if $is_ok;

        push @error_messages, $error_message;
    }

    return 1 if ( @error_messages == 0 );

  ERROR_MESSAGE:
    foreach my $error_message (@error_messages) {

        $log->fatal($error_message);
    }

    $log->fatal(q{Please install missing image files});
    exit 1;
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

    use MIP::Active_parameter qw{ update_with_dynamic_config_parameters };
    use MIP::Config qw{ get_install_containers };
    use MIP::Constants qw{ set_container_cmd };

    %{ $active_parameter_href->{container} } =
      get_install_containers(
        { container_config_file => $active_parameter_href->{container_config_file}, } );

    check_installed_containers(
        {
            container_href => $active_parameter_href->{container},
        }
    );

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

sub parse_container_config {

## Function : Parse container config and write to conda env
## Returns  :
## Arguments: $conda_environment_path => Conda environment path
##          : $container_href         => Container hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_environment_path;
    my $container_href;

    my $tmpl = {
        conda_environment_path => {
            store       => \$conda_environment_path,
            strict_type => 1,
        },
        container_href => {
            default     => {},
            store       => \$container_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Config qw{ get_install_containers };
    use MIP::Io::Write qw{ write_to_file };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Replace the uri with path
  IMAGE:
    foreach my $image ( keys %{$container_href} ) {

        if ( $container_href->{$image}{path} ) {

            $container_href->{$image}{uri} = delete $container_href->{$image}{path};
        }
    }

    ## Make directory if it doesn't exist
    my $container_config_dir_path = catdir( $conda_environment_path, qw{ bin templates } );

    make_path($container_config_dir_path);

    ## Fetch possible old config
    my $container_config_path = catfile( $container_config_dir_path, q{mip_container_config.yaml} );
    my %container_config      = ( not -e $container_config_path ) ? () : get_install_containers(
        {
            container_config_file => $container_config_path,
        }
    );

    ## Merge hashes
    %container_config = ( %container_config, %{$container_href} );

    $log->info( q{Writing container config to} . $COLON . $SPACE . $container_config_path );
    write_to_file(
        {
            data_href => \%container_config,
            format    => q{yaml},
            path      => $container_config_path,
        }
    );
    return;
}

sub parse_container_path {

## Function : Parse container download path
## Returns  :
## Arguments: $conda_environment_path   => Conda environment path
##          : $container_directory_path => Container directory path
##          : $container_href           => Container hash {REF}
##          : $container_manager        => Container manager
##          : $local_install            => Local install switch

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_environment_path;
    my $container_directory_path;
    my $container_href;
    my $container_manager;
    my $local_install;

    my $tmpl = {
        conda_environment_path => {
            store       => \$conda_environment_path,
            strict_type => 1,
        },
        container_directory_path => {
            store       => \$container_directory_path,
            strict_type => 1,
        },
        container_href => {
            default     => {},
            store       => \$container_href,
            strict_type => 1,
        },
        container_manager => {
            allow       => [qw{ docker singularity }],
            required    => 1,
            store       => \$container_manager,
            strict_type => 1,
        },
        local_install => {
            allow       => [ undef, 0, 1 ],
            store       => \$local_install,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Don't set variable for docker and unless local singularity installation
    return if $container_manager eq q{docker};
    return if not $local_install;

    ## Check for user input and default
    $container_directory_path = $container_directory_path
      // catdir( $conda_environment_path, q{bin} );

    make_path($container_directory_path);

    ## Set container outpath
    $container_href->{path} =
      catfile( $container_directory_path, fileparse( $container_href->{uri} ) . q{.sif} );

    $log->info( q{Saving image to} . $COLON . $SPACE . $container_href->{path} );

    return;
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

    return if ( ${$uri_ref} =~ m{ \A docker:[/]{2} }xms );

    if ( ${$uri_ref} =~ /\A ghcr|quay|docker[.]io /xms ) {

        ${$uri_ref} = q{docker://} . ${$uri_ref};
    }

    return;
}

sub pull_container {

## Function : Pull a docker or singularity container
## Returns  : @commands
## Arguments: $container_manager      => Container manager
##          : $container_path         => Path to container
##          : $filehandle             => Filehandle to write to
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdinfile_path         => Stdinfile path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $container_manager;
    my $container_uri;
    my $container_outpath;
    my $filehandle;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $force;

    my $tmpl = {
        container_manager => {
            allow       => [qw{ docker singularity }],
            required    => 1,
            store       => \$container_manager,
            strict_type => 1,
        },
        container_outpath => {
            store       => \$container_outpath,
            strict_type => 1,
        },
        container_uri => {
            defined     => 1,
            required    => 1,
            store       => \$container_uri,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        force => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$force,
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
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Singularity qw{ singularity_pull };
    use MIP::Program::Docker qw{ docker_pull };

    my %container_api = (
        docker => {
            arg_href => {
                filehandle             => $filehandle,
                image                  => $container_uri,
                stderrfile_path        => $stderrfile_path,
                stderrfile_path_append => $stderrfile_path_append,
                stdoutfile_path        => $stdoutfile_path,
            },
            method => \&docker_pull,
        },
        singularity => {
            arg_href => {
                container_uri          => $container_uri,
                filehandle             => $filehandle,
                force                  => $force,
                outfile_path           => $container_outpath,
                stderrfile_path        => $stderrfile_path,
                stderrfile_path_append => $stderrfile_path_append,
                stdoutfile_path        => $stdoutfile_path,
            },
            method => \&singularity_pull,
        },
    );

    my @commands = $container_api{$container_manager}{method}
      ->( { %{ $container_api{$container_manager}{arg_href} } } );

    return @commands;
}

sub run_container {

## Function : Run a docker container or exec a singularity image
## Returns  : @commands
## Arguments: $bind_paths_ref         => Bind host directory to container {REF}
##          : $clean_env              => Start with clean environment
##          : $container_cmds_ref     => Cmds to be executed in container {REF}
##          : $container_manager      => Container manager
##          : $container_path         => Path to container
##          : $entrypoint             => Override container entrypoint
##          : $executable_name        => Name of the executable
##          : $filehandle             => Filehandle to write to
##          : $gpu_switch             => Use nvidia experimental support
##          : $image                  => Image to run
##          : $no_home                => Don't mount home if it isn't the current working dir
##          : $remove                 => Remove stopped container
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdinfile_path         => Stdinfile path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bind_paths_ref;
    my $clean_env;
    my $container_cmds_ref;
    my $container_manager;
    my $container_path;
    my $executable_name;
    my $filehandle;
    my $gpu_switch;
    my $no_home;
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
        clean_env => {
            allow       => [ undef, 0, 1 ],
            store       => \$clean_env,
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
        executable_name => {
            store       => \$executable_name,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        gpu_switch => {
            allow       => [ undef, 0, 1 ],
            store       => \$gpu_switch,
            strict_type => 1,
        },
        no_home => {
            allow       => [ undef, 0, 1 ],
            store       => \$no_home,
            strict_type => 1,
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
                bind_paths_ref         => $bind_paths_ref,
                clean_env              => $clean_env,
                container_cmds_ref     => $container_cmds_ref,
                filehandle             => $filehandle,
                image                  => $container_path,
                gpu_switch             => $gpu_switch,
                no_home                => $no_home,
                stderrfile_path        => $stderrfile_path,
                stderrfile_path_append => $stderrfile_path_append,
                stdoutfile_path        => $stdoutfile_path,
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
            allow       => [qw{docker singularity}],
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

    my %recipe_executable_bind_path = get_recipe_executable_bind_path(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
        }
    );

    my %container_cmd = build_container_cmd(
        {
            active_parameter_href            => $active_parameter_href,
            container_href                   => $container_href,
            container_manager                => $container_manager,
            recipe_executable_bind_path_href => \%recipe_executable_bind_path,
        }
    );

    return %container_cmd;
}

1;
