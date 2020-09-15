package MIP::Set::Parameter;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use List::Util qw{ uniq };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.37;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      set_conda_path
      set_container_bind_paths
      set_programs_for_installation
    };
}

sub set_conda_path {

## Function : Set path to conda
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Environment::Path qw{ get_conda_path is_binary_in_path };

    ## Check if conda is in path
    is_binary_in_path(
        {
            binary => q{conda},
        }
    );

    ## Get path to conda
    my $conda_path = get_conda_path( {} );

    ## Set path to conda
    $active_parameter_href->{conda_path} = $conda_path;

    ## Set path to conda env
    my $environment_name = $active_parameter_href->{environment_name};
    $active_parameter_href->{conda_prefix_path} =
      catdir( $active_parameter_href->{conda_path}, q{envs}, $environment_name );

    return;
}

sub set_programs_for_installation {

## Function : Process the lists of programs that has been selected for installation
##          : and update the environment packages
## Returns  :
## Arguments: $active_parameter_href => The entire active parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Array::Utils qw{ array_minus };
    use Data::Diver qw{ Dive };
    use MIP::Get::Parameter qw{ get_programs_for_shell_installation };
    use MIP::Check::Installation qw{ check_and_add_dependencies };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Check that the options supplied are compatible with each other
    if (    ( scalar @{ $active_parameter_href->{skip_programs} } > 0 )
        and ( scalar @{ $active_parameter_href->{select_programs} } > 0 ) )
    {
        $log->fatal(
q{"--skip_programs" and "--select_programs" are mutually exclusive command line options}
        );
        exit 1;
    }

    ## Set programs to install depending on pipeline
    if ( scalar @{ $active_parameter_href->{select_programs} } == 0 ) {

      PIPELINE:
        foreach my $pipeline ( @{ $active_parameter_href->{pipelines} } ) {

            push @{ $active_parameter_href->{select_programs} },
              @{ $active_parameter_href->{$pipeline} };
        }
    }
    @{ $active_parameter_href->{select_programs} } =
      uniq @{ $active_parameter_href->{select_programs} };

    ## Get programs that are to be installed via shell
    my @shell_programs_to_install = get_programs_for_shell_installation(
        {
            conda_programs_href        => $active_parameter_href->{conda},
            log                        => $log,
            prefer_shell               => $active_parameter_href->{prefer_shell},
            shell_install_programs_ref => $active_parameter_href->{shell_install},
            shell_programs_href        => $active_parameter_href->{shell},
        }
    );

    ## Remove the conda packages that has been selected to be installed via SHELL
    delete @{ $active_parameter_href->{conda} }{@shell_programs_to_install};

    ## Delete shell programs that are to be installed via conda instead of shell
    my @shell_programs_to_delete = keys %{ $active_parameter_href->{shell} };
    @shell_programs_to_delete =
      array_minus( @shell_programs_to_delete, @shell_programs_to_install );
    delete @{ $active_parameter_href->{shell} }{@shell_programs_to_delete};

    ## Solve the installation when the skip_program or select_program parameter has been used
  INSTALL_MODE:
    foreach my $install_mode (qw{ conda pip shell container }) {

        ## Remove programs that are to be skipped
        delete @{ $active_parameter_href->{$install_mode} }
          { @{ $active_parameter_href->{skip_programs} } };

        ## Remove all non-selected programs
        if ( scalar @{ $active_parameter_href->{select_programs} } > 0 ) {
            my @non_selects = keys %{ $active_parameter_href->{$install_mode} };
            @non_selects =
              array_minus( @non_selects, @{ $active_parameter_href->{select_programs} } );
            delete @{ $active_parameter_href->{$install_mode} }{@non_selects};
        }
    }

    ## Check and add dependencies that are needed for shell programs if they are missing from the programs that are to be installed via conda.
  SHELL_PROGRAM:
    foreach my $shell_program ( keys %{ $active_parameter_href->{shell} } ) {
        my $dependency_href =
          Dive( $active_parameter_href->{shell}, $shell_program, q{conda_dependency} );

        next SHELL_PROGRAM if ( not defined $dependency_href );
        check_and_add_dependencies(
            {
                conda_program_href => $active_parameter_href->{conda},
                dependency_href    => $dependency_href,
                log                => $log,
                shell_program      => $shell_program,
            }
        );
    }
    return;
}

sub set_container_bind_paths {

## Function : Set/add bind paths to container hash
## Returns  :
## Arguments: $bind_paths_ref  => Active parameter hash {REF}
##          : $contaienr_href  => Container hah {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bind_paths_ref;
    my $container_href;

    my $tmpl = {
        bind_paths_ref => {
            default     => [],
            required    => 1,
            store       => \$bind_paths_ref,
            strict_type => 1,
        },
        container_href => {
            default     => {},
            required    => 1,
            store       => \$container_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    if ( $container_href->{program_bind_paths} ) {

        push @{ $container_href->{program_bind_paths} }, @{$bind_paths_ref};
    }
    else {
        $container_href->{program_bind_paths} = $bind_paths_ref;
    }

    return;
}

1;
