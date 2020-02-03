package MIP::Config;

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
use List::MoreUtils qw { any };
use Readonly;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_cmd_config_vs_definition_file
      parse_config
      parse_dynamic_config_parameters
      set_config_to_active_parameters
      set_default_config_dynamic_parameters
      write_mip_config };
}

sub check_cmd_config_vs_definition_file {

## Function : Compare keys from config and cmd with definitions file
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
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

    my @allowed_unique_keys =
      ( q{vcfparser_outfile_count}, $active_parameter_href->{case_id} );
    my @unique;

  ACTIVE_PARAMETER:
    foreach my $key ( keys %{$active_parameter_href} ) {

        ## Parameters from definitions file
        if ( not exists $parameter_href->{$key} ) {

            push @unique, $key;
        }
    }
  UNIQUE_KEYS:
    foreach my $unique_key (@unique) {

        ## Do not print if allowed_unique_keys that have been created dynamically from previous runs
        if ( not any { $_ eq $unique_key } @allowed_unique_keys ) {

            say {*STDERR} q{Found illegal key: }
              . $unique_key
              . q{ in config file or command line that is not defined in define_parameters.yaml};
            croak();
        }
    }
    return;
}

sub parse_config {

## Function : Parse config by loading from cli, removing prior dynamic entries,
##            setting to active_parameter and parsing config dynamic parameters
## Returns  : 1
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
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

    use MIP::Io::Read qw{ read_from_file };

    ## Loads a YAML file into an arbitrary hash and returns it.
    my %config_parameter = read_from_file(
        {
            format => q{yaml},
            path   => $active_parameter_href->{config_file},
        }
    );

    ## Remove previous analysis specific info not relevant for current run e.g. log file, which is read from pedigree or cmd
    my @remove_keys =
      qw{ found_female found_male found_other gender log_file dry_run_all };

    delete @config_parameter{@remove_keys};

## Set config parameters into %active_parameter unless $parameter
## has been supplied on the command line
    set_config_to_active_parameters(
        {
            active_parameter_href => $active_parameter_href,
            config_parameter_href => \%config_parameter,
        }
    );

    ## Compare keys from config and cmd (%active_parameter) with definitions file (%parameter)
    check_cmd_config_vs_definition_file(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
        }
    );

    my @config_dynamic_parameters = qw{ cluster_constant_path analysis_constant_path };

    ## Replace config parameter with cmd info for config dynamic parameter
    set_default_config_dynamic_parameters(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
            parameter_names_ref   => \@config_dynamic_parameters,
        }
    );

    ## Updates first the dynamic config parameters and then all other
    ## parameters to particular user/cluster following specifications
    parse_dynamic_config_parameters(
        {
            active_parameter_href         => $active_parameter_href,
            config_dynamic_parameters_ref => \@config_dynamic_parameters,
            parameter_href                => $parameter_href,
        }
    );
    return 1;
}

sub parse_dynamic_config_parameters {

## Function : Updates first the dynamic config parameters and then all other parameters to particular user/cluster following specifications
## Returns  :
## Arguments: $active_parameter_href         => Active parameters for this analysis hash {REF}
##          : $config_dynamic_parameters_ref => Config dynamic parameters
##          : $parameter_href                => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $config_dynamic_parameters_ref;
    my $parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        config_dynamic_parameters_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$config_dynamic_parameters_ref,
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

    use MIP::Update::Parameters qw{ update_dynamic_config_parameters };

    ## Loop through all config dynamic parameters and update value
  DYNAMIC_PARAM:
    foreach my $dynamic_param_name ( @{$config_dynamic_parameters_ref} ) {

        ## Updates the dynamic config parameters using supplied $case_id
        update_dynamic_config_parameters(
            {
                active_parameter_href => $active_parameter_href,
                dynamic_parameter_href =>
                  { case_id => $active_parameter_href->{case_id}, },
                parameter_name => $dynamic_param_name,
            }
        );
    }

    ## Map of dynamic parameters to update all other parameters
    my %dynamic_parameter = (
        cluster_constant_path  => $active_parameter_href->{cluster_constant_path},
        analysis_constant_path => $active_parameter_href->{analysis_constant_path},
        case_id                => $active_parameter_href->{case_id},
    );

    ## Loop through all parameters and update info
  PARAMETER:
    foreach my $parameter_name ( keys %{$parameter_href} ) {

        ## Updates the active parameters to particular user/cluster for dynamic config parameters following specifications. Leaves other entries untouched.
        update_dynamic_config_parameters(
            {
                active_parameter_href  => $active_parameter_href,
                dynamic_parameter_href => \%dynamic_parameter,
                parameter_name         => $parameter_name,
            }
        );
    }
    return 1;
}

sub set_config_to_active_parameters {

## Function : Add contig parameters to active_parameters if not already initilized from command line
## Returns  :
## Arguments: $active_parameter_href  => Active parameters for this analysis hash {REF}
##          : $config_parameter_href => Config parameters hash

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $config_parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        config_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$config_parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  PARAMETER:
    foreach my $parmeter_name ( keys %{$config_parameter_href} ) {

        ## Cmd initilized HASH
        next PARAMETER
          if ( ref $active_parameter_href->{$parmeter_name} eq qw{HASH}
            && keys %{ $active_parameter_href->{$parmeter_name} } );

        ## Cmd initilized ARRAY
        next PARAMETER
          if ( ref $active_parameter_href->{$parmeter_name} eq qw{ARRAY}
            && @{ $active_parameter_href->{$parmeter_name} } );

        ## Cmd initilized scalar
        next PARAMETER
          if ( defined $active_parameter_href->{$parmeter_name}
            and not ref $active_parameter_href->{$parmeter_name} );

        ### No input from cmd
        ## Add to active_parameter
        $active_parameter_href->{$parmeter_name} =
          $config_parameter_href->{$parmeter_name};
    }
    return;
}

sub set_default_config_dynamic_parameters {

## Function : Set default for config dynamic parameter using default definitions
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $parameter_href        => Parameter hash {REF}
##          : $parameter_names_ref   => MIP activate parameter names {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_href;
    my $parameter_names_ref;

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
        parameter_names_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$parameter_names_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  PARAMETER:
    foreach my $parameter_name ( @{$parameter_names_ref} ) {

        if ( exists $parameter_href->{$parameter_name}{default}
            and not defined $active_parameter_href->{$parameter_name} )
        {

            ## Transfer to active parameter
            $active_parameter_href->{$parameter_name} =
              $parameter_href->{$parameter_name}{default};
        }
    }
    return;
}

sub write_mip_config {

## Function : Write config file for analysis
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $log                   => Log object
##          : $remove_keys_ref       => Keys to remove before writing to file {REF}
##          : $sample_info_href      => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $log;
    my $remove_keys_ref;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        remove_keys_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$remove_keys_ref,
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

    use File::Basename qw{ dirname };
    use File::Path qw{ make_path };
    use MIP::File::Format::Yaml qw{ write_yaml };

    return if ( not $active_parameter_href->{config_file_analysis} );

    ## Create directory unless it already exists
    make_path( dirname( $active_parameter_href->{config_file_analysis} ) );

    ## Remove previous analysis specific info not relevant for current run e.g. log file, sample_ids which are read from pedigree or cmd
    delete @{$active_parameter_href}{ @{$remove_keys_ref} };

    ## Writes a YAML hash to file
    write_yaml(
        {
            yaml_href      => $active_parameter_href,
            yaml_file_path => $active_parameter_href->{config_file_analysis},
        }
    );
    $log->info( q{Wrote: } . $active_parameter_href->{config_file_analysis} );

    ## Add to sample_info for use downstream
    $sample_info_href->{config_file_analysis} =
      $active_parameter_href->{config_file_analysis};
    return;
}

1;
