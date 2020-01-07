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
use MIP::Constants qw{ $COLON $NEWLINE $SINGLE_QUOTE $TAB };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_parameter_hash get_order_of_parameters set_cache };
}

sub check_parameter_hash {

## Function : Evaluate parameters in parameters hash
## Returns  :
## Arguments: $file_path          => Path to yaml file
##          : $mandatory_href     => Hash with mandatory key {REF}
##          : $not_mandatory_href => Hash with non mandatory key {REF}
##          : $parameter_href     => Hash with parameters from yaml file {REF}

    my ($arg_href) = @_;

    ##Flatten argument(s)
    my $file_path;
    my $mandatory_href;
    my $not_mandatory_href;
    my $parameter_href;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        mandatory_href => {
            default     => {},
            required    => 1,
            store       => \$mandatory_href,
            strict_type => 1,
        },
        not_mandatory_href => {
            default     => {},
            required    => 1,
            store       => \$not_mandatory_href,
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

    ## Check that mandatory keys exists for each parameter
    _check_parameter_mandatory_keys_exits(
        {
            file_path      => $file_path,
            mandatory_href => $mandatory_href,
            parameter_href => $parameter_href,
        }
    );

    ## Test parameter for both mandatory and not_mandatory keys data type and values
    my @keys = ( \%{$mandatory_href}, \%{$not_mandatory_href} );

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

sub set_cache {

## Function : Set aggregate information from parameter hash to parameter cache
## Returns  :
## Arguments: $aggregates_ref => Data to aggregate {REF}
##          : $parameter_href => Parameter hash {REF}

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
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
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

sub _check_parameter_mandatory_keys_exits {

## Function : Check that mandatory keys exists
## Returns  :
## Arguments: $file_path      => Path to yaml file
##          : $mandatory_href => Hash with mandatory key {REF}
##          : $parameter_href => Hash with parameters from yaml file {REF}

    my ($arg_href) = @_;

    ##Flatten argument(s)
    my $file_path;
    my $mandatory_href;
    my $parameter_href;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        mandatory_href => {
            default     => {},
            required    => 1,
            store       => \$mandatory_href,
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

      MANDATORY_KEY:
        foreach my $mandatory_key ( keys %{$mandatory_href} ) {

            next MANDATORY_KEY
              if ( exists $parameter_href->{$parameter}{$mandatory_key} );

            say {*STDERR} q{Missing mandatory key: '}
              . $mandatory_key
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
