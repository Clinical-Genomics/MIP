package MIP::Set::Parameter;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use List::Util qw{ any };
use Readonly;

## MIPs lib/
use MIP::Constants
  qw{ $COLON $COMMA $CLOSE_BRACE $CLOSE_BRACKET $LOG_NAME $OPEN_BRACE $OPEN_BRACKET $SPACE $TAB };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.31;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      set_conda_path
      set_nist_file_name_path
      set_no_dry_run_parameters
      set_parameter_to_broadcast
      set_programs_for_installation
      set_recipe_mode
    };
}

## Constants
Readonly my $TWO         => 2;
Readonly my $ONE_HUNDRED => 100;

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

sub set_nist_file_name_path {

## Function : Set nist file name path by adding reference directory
## Returns  : 1
## Arguments: $active_parameter_href => Holds all set parameter for analysis
##          : $nist_parameters_ref   => Nist parameters to check

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $nist_parameters_ref;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        nist_parameters_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$nist_parameters_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Unpack
    my $reference_dir = $active_parameter_href->{reference_dir};

  NIST_PARAMETER:
    foreach my $nist_parameter ( @{$nist_parameters_ref} ) {

        # Alias
        my $nist_href = \%{ $active_parameter_href->{$nist_parameter} };

      NIST_VERSION:
        foreach my $nist_version ( keys %{$nist_href} ) {

          NIST_FILE:
            while ( my ( $nist_id, $file_name ) = each %{ $nist_href->{$nist_version} } )
            {

                ## Add reference directory to path
                $nist_href->{$nist_version}{$nist_id} =
                  catfile( $reference_dir, $file_name );
            }
        }
    }
    return 1;
}

sub set_no_dry_run_parameters {

## Function : Set parameters for true run i.e. not a dry run
## Returns  :
## Arguments: $analysis_date    => Analysis date
##          : $is_dry_run_all   => Dry run boolean
##          : $mip_version      => MIP version
##          : $sample_info_href => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $analysis_date;
    my $is_dry_run_all;
    my $mip_version;
    my $sample_info_href;

    my $tmpl = {
        analysis_date => {
            defined     => 1,
            required    => 1,
            store       => \$analysis_date,
            strict_type => 1,
        },
        is_dry_run_all => {
            allow       => [ 0, 1, undef ],
            required    => 1,
            store       => \$is_dry_run_all,
            strict_type => 1,
        },
        mip_version => {
            defined     => 1,
            required    => 1,
            store       => \$mip_version,
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

    return if ($is_dry_run_all);

    my %no_dry_run_info = (
        analysisrunstatus => q{not_finished},
        analysis_date     => $analysis_date,
        mip_version       => $mip_version,
    );

  KEY_VALUE_PAIR:
    while ( my ( $key, $value ) = each %no_dry_run_info ) {

        $sample_info_href->{$key} = $value;
    }

    return;
}

sub set_parameter_to_broadcast {

## Function : Set parameters to broadcast message
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $broadcasts_ref        => Holds the parameters info for broadcasting later {REF}
##          : $order_parameters_ref  => Order of parameters (for structured output) {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $broadcasts_ref;
    my $order_parameters_ref;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        broadcasts_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$broadcasts_ref,
            strict_type => 1,
        },
        order_parameters_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$order_parameters_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  PARAMETER:
    foreach my $parameter_name ( @{$order_parameters_ref} ) {

        next PARAMETER
          if ( not defined $active_parameter_href->{$parameter_name} );

        ## Hold parameters info
        my $info = q{Set } . $parameter_name . q{ to: };

        if ( ref $active_parameter_href->{$parameter_name} eq q{ARRAY} ) {

            $info = _parse_parameter_to_broadcast(
                {
                    info  => $info,
                    value => $active_parameter_href->{$parameter_name},
                }
            );

            ## Add info to broadcasts
            push @{$broadcasts_ref}, $info;
        }
        elsif ( ref $active_parameter_href->{$parameter_name} eq q{HASH} ) {

            $info = _parse_parameter_to_broadcast(
                {
                    info  => $info,
                    value => $active_parameter_href->{$parameter_name},
                }
            );

            ## Add info to broadcasts
            push @{$broadcasts_ref}, $info;
        }
        else {

            $info .= $active_parameter_href->{$parameter_name};

            ## Add info to broadcasts
            push @{$broadcasts_ref}, $info;
        }
    }
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
    foreach my $install_mode (qw{ conda pip shell singularity }) {

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

sub set_recipe_mode {

## Function : Set recipe mode
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis {REF}
##          : $mode                  => Mode to set
##          : $recipes_ref           => Recipes to set {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $mode;
    my $recipes_ref;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        mode => {
            allow       => [ 0, 1, $TWO ],
            defined     => 1,
            required    => 1,
            store       => \$mode,
            strict_type => 1,
        },
        recipes_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$recipes_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Set recipe mode
  RECIPE:
    foreach my $recipe ( @{$recipes_ref} ) {

        $active_parameter_href->{$recipe} = $mode;

        ## Broadcast
        $log->info(
            q{Set} . $SPACE . $recipe . $SPACE . q{to} . $COLON . $SPACE . $mode );
    }

    return;
}

sub _parse_parameter_to_broadcast {

## Function : Parse parameter to broadcast
## Returns  : $info
## Arguments: $info  => String to broadcast
##          : $value => Value to parse

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $info;
    my $value;

    my $tmpl = {
        info => {
            defined     => 1,
            required    => 1,
            store       => \$info,
            strict_type => 1,
        },
        value => {
            defined  => 1,
            required => 1,
            store    => \$value,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## HASH
    if ( ref $value eq q{HASH} ) {

        ## Start of hash
        $info .= $OPEN_BRACE;

      KEY:
        foreach my $key ( keys %{$value} ) {

            ## Key-value pairs
            $info .= $key . q{ => };

            if ( ref $value->{$key} eq q{HASH} ) {

                $info = _parse_parameter_to_broadcast(
                    {
                        info  => $info,
                        value => $value->{$key},
                    }
                );
                $info .= $COMMA . $SPACE;
                next KEY;
            }
            if ( ref $value->{$key} eq q{ARRAY} ) {

                $info .= $OPEN_BRACKET;

              ELEMENT:
                foreach my $element ( @{ $value->{$key} } ) {

                    $info = _parse_parameter_to_broadcast(
                        {
                            info  => $info,
                            value => $element,
                        }
                    );
                }
                ## Close array
                $info .= $CLOSE_BRACKET . $COMMA . $SPACE;
                next KEY;
            }
            if ( $value->{$key} ) {

                ## Scalar
                $info .= $value->{$key} . $COMMA . $SPACE;
            }
        }
        ## Close hash
        $info .= $CLOSE_BRACE;
        return $info;
    }
    ## ARRAY
    if ( ref $value eq q{ARRAY} ) {

        ## Open array
        $info .= $OPEN_BRACKET;

      ELEMENT:
        foreach my $element ( @{$value} ) {

            if ( ref $element eq q{HASH} ) {

                $info = _parse_parameter_to_broadcast(
                    {
                        info  => $info,
                        value => $element,
                    }
                );
                $info .= $COMMA . $SPACE;
                next ELEMENT;
            }
            if ( ref $element eq q{ARRAY} ) {

                $info .= $OPEN_BRACKET;

                foreach my $elements_ref ( @{$element} ) {

                    $info = _parse_parameter_to_broadcast(
                        {
                            info  => $info,
                            value => $elements_ref,
                        }
                    );
                }
                ## Close array
                $info .= $CLOSE_BRACKET . $COMMA . $SPACE;
                next ELEMENT;
            }
            if ($element) {

                ## Scalar
                $info .= $element . $COMMA . $SPACE;
            }
        }
        $info .= $CLOSE_BRACKET . $SPACE;
        return $info;
    }

    ## Scalar
    $info .= $value . $COMMA . $SPACE;
    return $info;
}

1;
