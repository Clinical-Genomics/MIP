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
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ install_containers parse_container_bind_paths parse_container_uri run_container };
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
                conda_env_path        => $active_parameter_href->{conda_prefix_path},
                container_href        => $active_parameter_href->{singularity},
                filehandle            => $filehandle,
            },
            method => \&install_docker_containers,
        },
        singularity => {
            arg_href => {
                active_parameter_href => $active_parameter_href,
                conda_env_path        => $active_parameter_href->{conda_prefix_path},
                container_href        => $active_parameter_href->{singularity},
                filehandle            => $filehandle,
            },
            method => \&install_singularity_containers,
        },
    );

    $container_api{$container_manager}{method}
      ->( { %{ $container_api{$container_manager}{arg_href} } } );

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

    use MIP::Environment::Path
      qw{ build_docker_bind_path_var build_singularity_bind_path_var reduce_dir_paths };

    ## For testing
    return if not $CONTAINER_MANAGER;

    my @export_bind_paths = @CONTAINER_BIND_PATHS;

    ## Look for extra bind paths
    if ( $active_parameter_href->{recipe_bind_path}{$package_name} ) {

        ## Add extra paths
        push @export_bind_paths,
          @{ $active_parameter_href->{recipe_bind_path}{$package_name} };

        ## Check for redundant paths
        @export_bind_paths = reduce_dir_paths( { dir_paths_ref => \@export_bind_paths } );
    }

    my %container = (
        docker      => \&_build_docker_bind_path_var,
        singularity => \&_build_singularity_bind_path_var,
    );

    my $container_bind_var = $container{$CONTAINER_MANAGER}->(
        {
            bind_paths_ref => \@export_bind_paths,
        }
    );

    push @{$source_environment_cmds_ref}, $container_bind_var;

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

    if ( ${$uri_ref} =~ /\A docker\.io /xms ) {

        ${$uri_ref} = q{docker://} . ${$uri_ref};
    }

    return;
}

sub run_container {

## Function : Run/Exec an image/container docker or singularity
## Returns  : @commands
## Arguments: $bind_paths_ref         => Bind host directory to container {REF}
##          : $container_cmds_ref     => Cmds to be executed in container {REF}
##          : $container_manager     => Container manager
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
                singularity_container          => $container_path,
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

sub _build_docker_bind_path_var {

## Function : Build bind path variable for use with docker
## Returns  : $mip_bind_var
## Arguments: $bind_paths_ref => Directories to be mounted {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bind_paths_ref;

    my $tmpl = {
        bind_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$bind_paths_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @mip_bind_paths = map { $_ . $COLON . $_ } @{$bind_paths_ref};

    my $mip_bind = join $SPACE . q{--volume} . $SPACE, @mip_bind_paths;

    my $mip_bind_var =
        q{export MIP_BIND}
      . $EQUALS
      . $DOUBLE_QUOTE
      . $mip_bind
      . $DOUBLE_QUOTE
      . $SEMICOLON;

    return $mip_bind_var;
}

sub _build_singularity_bind_path_var {

## Function : Build bind path variable for use with singularity
## Returns  : $singularity_bind_var
## Arguments: $bind_paths_ref => Directories to be mounted {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bind_paths_ref;

    my $tmpl = {
        bind_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$bind_paths_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $singularity_bind = join $COMMA, @{$bind_paths_ref};

    my $singularity_bind_var =
        q{export SINGULARITY_BIND}
      . $EQUALS
      . $DOUBLE_QUOTE
      . $singularity_bind
      . $DOUBLE_QUOTE
      . $SEMICOLON;

    return $singularity_bind_var;
}

1;
