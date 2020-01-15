package MIP::Parameter;

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
use List::MoreUtils qw { any };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $COLON $NEWLINE $SINGLE_QUOTE $SPACE $TAB };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.05;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      check_parameter_hash
      get_capture_kit
      get_order_of_parameters
      print_recipe
      set_cache
      set_cache_sample_id_parameter
    };
}

sub check_parameter_hash {

## Function : Evaluate parameters in parameters hash
## Returns  :
## Arguments: $file_path         => Path to yaml file
##          : $required_href     => Hash with required key {REF}
##          : $not_required_href => Hash with non required key {REF}
##          : $parameter_href    => Hash with parameters from yaml file {REF}

    my ($arg_href) = @_;

    ##Flatten argument(s)
    my $file_path;
    my $required_href;
    my $not_required_href;
    my $parameter_href;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        required_href => {
            default     => {},
            required    => 1,
            store       => \$required_href,
            strict_type => 1,
        },
        not_required_href => {
            default     => {},
            required    => 1,
            store       => \$not_required_href,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Check that required keys exists for each parameter
    _check_parameter_required_keys_exits(
        {
            file_path      => $file_path,
            required_href  => $required_href,
            parameter_href => $parameter_href,
        }
    );

    ## Test parameter for both required and not_required keys data type and values
    my @keys = ( \%{$required_href}, \%{$not_required_href} );

  KEY_HASH_REF:
    foreach my $key_href (@keys) {

        _check_parameter_keys(
            {
                file_path      => $file_path,
                key_href       => $key_href,
                parameter_href => $parameter_href,
            }
        );
    }
    return 1;
}

sub get_capture_kit {

## Function : Return a capture kit depending on user info. If $is_set_by_user is
##          : true, go a head and add capture kit no matter what the switch was.
## Returns  : $capture kit, "supported capture kit" or "undef"
## Arguments: $capture_kit                => Capture kit to add
##          : $is_set_by_user             => Has user supplied parameter {OPTIONAL}
##          : $supported_capture_kit_href => Supported capture kits hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $capture_kit;
    my $is_set_by_user;
    my $supported_capture_kit_href;

    my $tmpl = {
        capture_kit                => { store => \$capture_kit,    strict_type => 1, },
        is_set_by_user             => { store => \$is_set_by_user, strict_type => 1, },
        supported_capture_kit_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$supported_capture_kit_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Unpack
    my $supported_capture_kit = $supported_capture_kit_href->{$capture_kit};

    ## Set default or return supplied capture kit
    if ( not defined $is_set_by_user ) {

        return $supported_capture_kit if ( defined $supported_capture_kit );

        ## Return unchanged capture_kit string
        return $capture_kit;
    }
    ## Only add if user supplied no info on parameter
    if ( defined $is_set_by_user
        and not $is_set_by_user )
    {

        return $supported_capture_kit if ( defined $supported_capture_kit );

        ## Return unchanged capture_kit string
        return $capture_kit;
    }
    return;
}

sub get_order_of_parameters {

## Function : Get order of parameters as they appear in definition file(s)
## Returns  : @order_of_parameters
## Arguments: $define_parameters_files_ref => MIPs define parameters file(s)

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $define_parameters_files_ref;

    my $tmpl = {
        define_parameters_files_ref => {
            default     => [],
            store       => \$define_parameters_files_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Definition qw{ get_first_level_keys_order_from_definition_file };

    my @order_of_parameters;

  DEFINITION_FILE:
    foreach my $define_parameters_file ( @{$define_parameters_files_ref} ) {

        push @order_of_parameters,
          get_first_level_keys_order_from_definition_file(
            {
                file_path => $define_parameters_file,
            }
          );
    }
    return @order_of_parameters;
}

sub print_recipe {

## Function : Print all supported recipes in '-prm' mode if requested and then exit
## Returns  :
## Arguments: $order_parameters_ref => Order of addition to parameter array {REF}
##          : $parameter_href       => Parameter hash {REF}
##          : $print_recipe         => Print recipes switch
##          : $print_recipe_mode    => Mode to run recipes in

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $order_parameters_ref;
    my $parameter_href;
    my $print_recipe;

    ## Default(s)
    my $print_recipe_mode;

    my $tmpl = {
        order_parameters_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$order_parameters_ref,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        print_recipe => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$print_recipe,
            strict_type => 1,
        },
        print_recipe_mode => {
            allow       => [ undef, 0, 1, 2 ],
            default     => $arg_href->{print_recipe_mode} //= 2,
            store       => \$print_recipe_mode,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Definition qw{ get_first_level_keys_order_from_definition_file };
    use MIP::Parameter qw{ set_cache };

    ## Do not print
    return if ( not $print_recipe );

    set_cache(
        {
            aggregates_ref => [q{type:recipe}],
            parameter_href => $parameter_href,
        }
    );

  PARAMETER:
    foreach my $parameter ( @{$order_parameters_ref} ) {

        ## Only process recipes
        if (
            any { $_ eq $parameter }
            @{ $parameter_href->{cache}{recipe} }
          )
        {

            print {*STDOUT} q{--} . $parameter . $SPACE . $print_recipe_mode . $SPACE;

        }
    }
    print {*STDOUT} $NEWLINE;
    exit;
}

sub set_cache {

## Function : Set aggregate information from parameter hash to parameter cache
## Returns  :
## Arguments: $aggregates_ref => Data to aggregate {REF}
##          : $parameter_href => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $aggregates_ref;
    my $parameter_href;

    my $tmpl = {
        aggregates_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$aggregates_ref,
            strict_type => 1,
        },
        parameter_href => {
            default  => {},
            defined  => 1,
            required => 1,
            store    => \$parameter_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Constants
    Readonly my $FIELD_COUNTER => 2;

  PARAMETER:
    foreach my $parameter_name ( keys %{$parameter_href} ) {

      KEY_AND_STRING_TO_MATCH:
        foreach my $aggregate_element ( @{$aggregates_ref} ) {

            ## Split into key and string to match
            my ( $second_key, $string_to_match, $unexpected_data ) =
              split $COLON, $aggregate_element, $FIELD_COUNTER + 1;

            ## Make sure that we get what we expect
            if ( defined $unexpected_data ) {

                carp q{Unexpected trailing garbage at end of aggregate_element '}
                  . $aggregate_element
                  . q{':}, $NEWLINE . $TAB . $unexpected_data . $NEWLINE;
            }

            next KEY_AND_STRING_TO_MATCH
              if ( not defined $parameter_href->{$parameter_name}{$second_key} );

            next KEY_AND_STRING_TO_MATCH
              if ( $parameter_href->{$parameter_name}{$second_key} ne $string_to_match );

            push @{ $parameter_href->{cache}{$string_to_match} }, $parameter_name;
        }
    }
    return;
}

sub set_cache_sample_id_parameter {

## Function : Set parameter information to parameter cache at sample level
## Returns  :
## Arguments: $parameter_href  => Parameter hash {REF}
##          : $parameter_name  => Parameter to set
##          : $parameter_value => Parmeter value
##          : $sample_id       => Sample id to set parameter cache for

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $parameter_name;
    my $parameter_value;
    my $sample_id;

    my $tmpl = {
        parameter_href => {
            default  => {},
            defined  => 1,
            required => 1,
            store    => \$parameter_href,
        },
        parameter_name => {
            defined     => 1,
            required    => 1,
            store       => \$parameter_name,
            strict_type => 1,
        },
        parameter_value => {
            defined     => 1,
            required    => 1,
            store       => \$parameter_value,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    $parameter_href->{cache}{$sample_id}{$parameter_name} = $parameter_value;
    return;
}

sub _check_parameter_required_keys_exits {

## Function : Check that required keys exists
## Returns  :
## Arguments: $file_path      => Path to yaml file
##          : $required_href  => Hash with required key {REF}
##          : $parameter_href => Hash with parameters from yaml file {REF}

    my ($arg_href) = @_;

    ##Flatten argument(s)
    my $file_path;
    my $required_href;
    my $parameter_href;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        required_href => {
            default     => {},
            required    => 1,
            store       => \$required_href,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  PARAMETER:
    foreach my $parameter ( keys %{$parameter_href} ) {

      REQUIRED_KEY:
        foreach my $required_key ( keys %{$required_href} ) {

            next required_KEY
              if ( exists $parameter_href->{$parameter}{$required_key} );

            say {*STDERR} q{Missing required key: '}
              . $required_key
              . q{' for parameter: '}
              . $parameter
              . q{' in file: '}
              . $file_path
              . $SINGLE_QUOTE
              . $NEWLINE;
            croak();
        }
    }
    return;
}

sub _check_parameter_keys {

## Function : Evaluate parameter keys in hash
## Returns  :
## Arguments: $file_path      => Path to yaml file
##          : $key_href       => Hash with keys to check for parameter {REF}
##          : $parameter_href => Hash with parameters from yaml file {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path;
    my $key_href;
    my $parameter_href;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        key_href => {
            default     => {},
            required    => 1,
            store       => \$key_href,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  PARAMETER:
    foreach my $parameter ( keys %{$parameter_href} ) {

      KEY:
        foreach my $key ( keys %{$key_href} ) {

            next KEY if ( not exists $parameter_href->{$parameter}{$key} );

            ## Check key data type
            _check_parameter_data_type(
                {
                    file_path      => $file_path,
                    key            => $key,
                    key_href       => $key_href,
                    parameter      => $parameter,
                    parameter_href => $parameter_href,
                }
            );

            ## Check key values
            _check_parameter_values(
                {
                    file_path      => $file_path,
                    key            => $key,
                    key_href       => $key_href,
                    parameter      => $parameter,
                    parameter_href => $parameter_href,
                }
            );
        }
    }
    return;
}

sub _check_parameter_data_type {

## Function : Check key data type
## Returns  :
## Arguments: $file_path      => Path to yaml file
##          : $key            => Hash with non key
##          : $key_href       => Hash with key {REF}
##          : $parameter      => Parameter
##          : $parameter_href => Hash with parameters from yaml file {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path;
    my $key;
    my $key_href;
    my $parameter;
    my $parameter_href;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        key      => { defined => 1, required => 1, store => \$key, strict_type => 1, },
        key_href => {
            default     => {},
            required    => 1,
            store       => \$key_href,
            strict_type => 1,
        },
        parameter => {
            defined     => 1,
            required    => 1,
            store       => \$parameter,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Get actual parameter value data type
    my $data_type = ref $parameter_href->{$parameter}{$key};

    ## Array or hash
    if ($data_type) {

        ## Return if parameter value data type equals description
        return if ( $data_type eq $key_href->{$key}{key_data_type} );

        say {*STDERR} q{Found '}
          . $data_type
          . q{' but expected datatype '}
          . $key_href->{$key}{key_data_type}
          . q{' for parameter: '}
          . $parameter
          . q{' in key: '}
          . $key
          . q{' in file: '}
          . $file_path
          . $SINGLE_QUOTE
          . $NEWLINE;
        croak();
    }
    elsif ( $key_href->{$key}{key_data_type} ne q{SCALAR} ) {

        ## Wrong data_type
        say {*STDERR} q{Found 'SCALAR' but expected datatype '}
          . $key_href->{$key}{key_data_type}
          . q{' for parameter: '}
          . $parameter
          . q{' in key: '}
          . $key
          . q{' in file: '}
          . $file_path
          . $SINGLE_QUOTE
          . $NEWLINE;
        croak();
    }
    return;
}

sub _check_parameter_values {

## Function : Evaluate parameter key values
## Returns  :
## Arguments: $file_path      => Path to yaml file
##          : $key            => Hash with non  key
##          : $key_href       => Hash with  key {REF}
##          : $parameter      => Parameter
##          : $parameter_href => Hash with parameters from yaml file {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path;
    my $key;
    my $key_href;
    my $parameter;
    my $parameter_href;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        key      => { defined => 1, required => 1, store => \$key, strict_type => 1, },
        key_href => {
            default     => {},
            required    => 1,
            store       => \$key_href,
            strict_type => 1,
        },
        parameter => {
            defined     => 1,
            required    => 1,
            store       => \$parameter,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Check value(s)
    if ( $key_href->{$key}{values} ) {

        my $value = $parameter_href->{$parameter}{$key};

        if ( not( any { $_ eq $value } @{ $key_href->{$key}{values} } ) ) {

            say {*STDERR} q{Found illegal value '}
              . $value
              . q{' for parameter: '}
              . $parameter
              . q{' in key: '}
              . $key
              . q{' in file: '}
              . $file_path
              . $SINGLE_QUOTE
              . $NEWLINE
              . q{Allowed entries: '}
              . join( q{', '}, @{ $key_href->{$key}{values} } )
              . $SINGLE_QUOTE
              . $NEWLINE;
            croak();
        }
    }
    return;
}

1;
