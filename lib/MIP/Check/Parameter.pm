package MIP::Check::Parameter;

use Carp;
use charnames qw{ :full :short };
use File::Spec::Functions qw{ catfile };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie;
use Email::Valid;
use Readonly;
use List::MoreUtils qw { any };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      check_allowed_array_values
      check_allowed_temp_directory
      check_cmd_config_vs_definition_file
      check_email_address check_parameter_hash
      check_gzipped
    };
}

## Constants
Readonly my $SINGLE_QUOTE => q{'};
Readonly my $NEWLINE      => qq{\n};

sub check_allowed_array_values {

## Function : Check that the array values are allowed
## Returns  :
## Arguments: $allowed_values_ref => Allowed values for parameter {REF}
##          : $values_ref         => Values for parameter {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $allowed_values_ref;
    my $values_ref;

    my $tmpl = {
        allowed_values_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$allowed_values_ref,
            strict_type => 1,
        },
        values_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$values_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %is_allowed;

    # Remake allowed values into keys in is_allowed hash
    map { $is_allowed{$_} = undef } @{$allowed_values_ref};

  VALUES:
    foreach my $value ( @{$values_ref} ) {

        # Test if value is allowed
        if ( not exists $is_allowed{$value} ) {

            return 0;
        }
    }

    # All ok
    return 1;
}

sub check_allowed_temp_directory {

## Function : Check that the temp directory value is allowed
## Returns  :
## Arguments: $temp_directory => Temp directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $temp_directory;

    my $tmpl = {
        temp_directory => {
            defined     => 1,
            required    => 1,
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %is_not_allowed = (
        q{/scratch}  => undef,
        q{/scratch/} => undef,
    );

    # Test if value is allowed
    if ( exists $is_not_allowed{$temp_directory} ) {

        return 0;
    }

    # All ok
    return 1;
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
      ( q{vcfparser_outfile_count}, $active_parameter_href->{family_id} );
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

            croak(  q{Found illegal key: }
                  . $unique_key
                  . q{ in config file or command line that is not defined in define_parameters.yaml}
            );
        }
    }
    return;
}

sub check_email_address {

## Function : Check the syntax of the email adress is valid and has a mail host.
## Returns  :
## Arguments: $email => The email adress
##          : $log   => Log object

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $email;
    my $log;

    my $tmpl = {
        email => {
            defined     => 1,
            required    => 1,
            store       => \$email,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Check syntax and mail host
    my $address = Email::Valid->address(
        -address => $email,
        -mxcheck => 1,
    );
    if ( not defined $address ) {

        $log->fatal( q{The supplied email: }
              . $email
              . q{ seem to be malformed according to }
              . Email::Valid->details() );
        exit 1;
    }
    return;
}

sub check_parameter_hash {

## Function : Evaluate parameters in parameters hash
## Returns  :
## Arguments: $file_path              => Path to yaml file
##          : $mandatory_key_href     => Hash with mandatory key {REF}
##          : $non_mandatory_key_href => Hash with non mandatory key {REF}
##          : $parameter_href         => Hash with parameters from yaml file {REF}

    my ($arg_href) = @_;

    ##Flatten argument(s)
    my $file_path;
    my $mandatory_key_href;
    my $non_mandatory_key_href;
    my $parameter_href;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        mandatory_key_href => {
            default     => {},
            required    => 1,
            store       => \$mandatory_key_href,
            strict_type => 1,
        },
        non_mandatory_key_href => {
            default     => {},
            required    => 1,
            store       => \$non_mandatory_key_href,
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
            file_path          => $file_path,
            mandatory_key_href => $mandatory_key_href,
            parameter_href     => $parameter_href,
        }
    );

    ## Test both mandatory and non_mandatory keys data type and values
    my @arguments = ( \%{$mandatory_key_href}, \%{$non_mandatory_key_href} );

  ARGUMENT_HASH_REF:
    foreach my $argument_href (@arguments) {

        ## Mandatory keys
        _check_parameter_keys(
            {
                file_path      => $file_path,
                key_href       => $argument_href,
                parameter_href => $parameter_href,
            }
        );
    }
    return;
}

sub _check_parameter_mandatory_keys_exits {

## Function : Check that mandatory keys exists
## Returns  :
## Arguments: $file_path          => Path to yaml file
##          : $mandatory_key_href => Hash with mandatory key {REF}
##          : $parameter_href     => Hash with parameters from yaml file {REF}

    my ($arg_href) = @_;

    ##Flatten argument(s)
    my $file_path;
    my $mandatory_key_href;
    my $parameter_href;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        mandatory_key_href => {
            default     => {},
            required    => 1,
            store       => \$mandatory_key_href,
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
        foreach my $mandatory_key ( keys %{$mandatory_key_href} ) {

            ## Mandatory key exists
            if ( not exists $parameter_href->{$parameter}{$mandatory_key} ) {

                croak(  q{Missing mandatory key: '}
                      . $mandatory_key
                      . q{' for parameter: '}
                      . $parameter
                      . q{' in file: '}
                      . $file_path
                      . $SINGLE_QUOTE
                      . $NEWLINE );
            }
        }
    }
    return;
}

sub _check_parameter_keys {

## Function : Evaluate parameter keys in hash
## Returns  :
## Arguments: $file_path      => Path to yaml file
##          : $key_href       => Hash with mandatory key {REF}
##          : $parameter_href => Hash with parameters from yaml file {REF}

    my ($arg_href) = @_;

    ##Flatten argument(s)
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

            ## Key exists in parameter
            if ( exists $parameter_href->{$parameter}{$key} ) {

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

                ## Evaluate key values
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
        key =>
          { defined => 1, required => 1, store => \$key, strict_type => 1, },
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

            croak(  q{Found illegal value '}
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
                  . $NEWLINE );
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
##          : $parameter_href => Hash with paremters from yaml file {REF}

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
        key =>
          { defined => 1, required => 1, store => \$key, strict_type => 1, },
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

    ## Check data_type
    my $data_type = ref( $parameter_href->{$parameter}{$key} );

    ## Array or hash
    if ($data_type) {

        ## Wrong data_type
        if ( not $data_type eq $key_href->{$key}{key_data_type} ) {

            croak(  q{Found '}
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
                  . $NEWLINE );
        }
    }
    elsif ( $key_href->{$key}{key_data_type} ne q{SCALAR} ) {

        ## Wrong data_type
        croak(  q{Found 'SCALAR' but expected datatype '}
              . $key_href->{$key}{key_data_type}
              . q{' for parameter: '}
              . $parameter
              . q{' in key: '}
              . $key
              . q{' in file: '}
              . $file_path
              . $SINGLE_QUOTE
              . $NEWLINE );
    }
    return;
}

sub check_gzipped {

## Function : Check if a file is gzipped.
## Returns  : "0 (=uncompressed)| 1 (=compressed)"
## Arguments: $file_name => File name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_name;

    my $tmpl = {
        file_name => {
            defined     => 1,
            required    => 1,
            store       => \$file_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $file_compression_status = 0;

    if ( $file_name =~ / .gz$ /xms ) {

        $file_compression_status = 1;
    }
    return $file_compression_status;
}

1;

