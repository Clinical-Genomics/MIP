package MIP::Environment::Container;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
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

    # Set the version for version checking
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      build_container_cmd
      get_recipe_executable_bind_path
      parse_container_uri
      parse_containers
      run_container
      set_executable_container_cmd
    };
}

sub build_container_cmd {

    ## Function : Build executable command depending on container manager
    ## Returns  :
    ## Arguments: $active_parameter_href => The active parameters for this analysis hash {REF}
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

    my @container_constant_bind_path = @CONTAINER_BIND_PATHS;
    my %container_cmd;

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
                push @{ $recipe_executable_bind_path_href->{$executable_name} },
                  $container_href->{$container_name}{bind_path}{$executable_name};
            }
            my @bind_paths =
              exists $recipe_executable_bind_path_href->{$executable_name}
              ? @{ $recipe_executable_bind_path_href->{$executable_name} }
              : @container_constant_bind_path;

            my @cmds = run_container(
                {
                    active_parameter_href => $active_parameter_href,
                    bind_paths_ref        => \@bind_paths,
                    executable_name       => $executable_name,
                    container_manager     => $container_manager,
                    container_path        => $container_href->{$container_name}{uri},
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
    use MIP::Language::Shell qw{ quote_bash_variable };
    use MIP::Parameter qw{get_cache get_parameter_attribute};

    my $temp_directory_quoted = quote_bash_variable(
        { string_with_variable_to_quote => $active_parameter_href->{temp_directory}, } );
    my $xdg_runtime_dir =
        $temp_directory_quoted
      . $COLON
      . catfile( $EMPTY_STR, qw{ run user }, q{$(id -u)} );

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

            ## Special case for xdg_runtime_dir, which always should be added
            push @export_bind_paths, $xdg_runtime_dir;

            $recipe_executable_bind_path{$executable} = [@export_bind_paths];
        }
    }
    return %recipe_executable_bind_path;
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

    if ( ${$uri_ref} =~ /\A quay|docker[.]io /xms ) {

        ${$uri_ref} = q{docker://} . ${$uri_ref};
    }

    return;
}

sub run_container {

## Function : Run a docker container or exec a singularity image
## Returns  : @commands
## Arguments: $active_parameter_href  => The active parameters for this analysis hash {REF}
##            $bind_paths_ref         => Bind host directory to container {REF}
##          : $container_cmds_ref     => Cmds to be executed in container {REF}
##          : $container_manager      => Container manager
##          : $container_path         => Path to container
##          : $entrypoint             => Override container entrypoint
##          : $executable_name        => Name of the executable
##          : $filehandle             => Filehandle to write to
##          : $image                  => Image to run
##          : $remove                 => Remove stopped container
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdinfile_path         => Stdinfile path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $bind_paths_ref;
    my $container_cmds_ref;
    my $container_manager;
    my $container_path;
    my $executable_name;
    my $filehandle;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdinfile_path;
    my $stdoutfile_path;

    ## Default(s)
    my $remove;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
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
        executable_name => {
            store       => \$executable_name,
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
                active_parameter_href  => $active_parameter_href,
                bind_paths_ref         => $bind_paths_ref,
                executable_name        => $executable_name,
                filehandle             => $filehandle,
                image                  => $container_path,
                container_cmds_ref     => $container_cmds_ref,
                stderrfile_path        => $stdinfile_path,
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
