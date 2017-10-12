package MIP::Check::Parameter;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };

# Allow unicode characters in this script
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use autodie;
use Params::Check qw{ check allow last_error };
use File::Spec::Functions qw{ catfile };
use FindBin qw{ $Bin };

## CPANM
use Readonly;
use List::MoreUtils qw { any };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ check_parameter_hash check_allowed_array_values check_allowed_temp_directory};
}

##Constants
Readonly my $SINGLE_QUOTE => q{'};
Readonly my $NEWLINE => qq{\n};

sub check_parameter_hash {

## Function : Evaluate parameters in parameters hash
## Returns  : ""
## Arguments: $parameter_href         => Hash with parameters from yaml file {REF}
##          : $mandatory_key_href     => Hash with mandatory key {REF}
##          : $non_mandatory_key_href => Hash with non mandatory key {REF}
##          : $file_path              => Path to yaml file

    my ($arg_href) = @_;

    ##Flatten argument(s)
    my $parameter_href;
    my $mandatory_key_href;
    my $non_mandatory_key_href;
    my $file_path;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        mandatory_key_href => {
            required    => 1,
            default     => {},
            strict_type => 1,
            store       => \$mandatory_key_href,
        },
        non_mandatory_key_href => {
            required    => 1,
            default     => {},
            strict_type => 1,
            store       => \$non_mandatory_key_href,
        },
        file_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_path,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Check that mandatory keys exists for each parameter
    my $error_msg = _check_parameter_mandatory_keys_exits(
        {
            parameter_href     => $parameter_href,
            mandatory_key_href => $mandatory_key_href,
            file_path          => $file_path,
        }
    );

    if ($error_msg) {

        return $error_msg;
    }
    ## Test both mandatory and non_mandatory keys data type and values
    my @arguments = ( \%{$mandatory_key_href}, \%{$non_mandatory_key_href} );

  ARGUMENT_HASH_REF:
    foreach my $argument_href (@arguments) {

        ## Mandatory keys
        $error_msg = _check_parameter_keys(
            {
                parameter_href => $parameter_href,
                key_href       => $argument_href,
                file_path      => $file_path,
            }
        );
        if ($error_msg) {

            return $error_msg;
        }
    }
    return;
}

sub _check_parameter_mandatory_keys_exits {

## Function : Check that mandatory keys exists
## Returns  :
## Arguments: $parameter_href     => Hash with parameters from yaml file {REF}
##          : $mandatory_key_href => Hash with mandatory key {REF}
##          : $file_path          => Path to yaml file

    my ($arg_href) = @_;

    ##Flatten argument(s)
    my $parameter_href;
    my $mandatory_key_href;
    my $file_path;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        mandatory_key_href => {
            required    => 1,
            default     => {},
            strict_type => 1,
            store       => \$mandatory_key_href,
        },
        file_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_path,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  PARAMETER:
    foreach my $parameter ( keys %{$parameter_href} ) {

      MANDATORY_KEY:
        foreach my $mandatory_key ( keys %{$mandatory_key_href} ) {

            ## Mandatory key exists
            if ( not exists $parameter_href->{$parameter}{$mandatory_key} ) {

                my $error_msg =
                    q{Missing mandatory key: '}
                  . $mandatory_key
                  . q{' for parameter: '}
                  . $parameter
                  . q{' in file: '}
                  . $file_path . $SINGLE_QUOTE . $NEWLINE;
                return $error_msg;
            }
        }
    }
    return;
}

sub _check_parameter_keys {

## Function : Evaluate parameter keys in hash
## Returns  :
## Arguments: $parameter_href => Hash with parameters from yaml file {REF}
##          : $key_href       => Hash with mandatory key {REF}
##          : $file_path      => Path to yaml file

    my ($arg_href) = @_;

    ##Flatten argument(s)
    my $parameter_href;
    my $key_href;
    my $file_path;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        key_href => {
            required    => 1,
            default     => {},
            strict_type => 1,
            store       => \$key_href,
        },
        file_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_path,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $error_msg;

  PARAMETER:
    foreach my $parameter ( keys %{$parameter_href} ) {

      KEY:
        foreach my $key ( keys %{$key_href} ) {

            ## Key exists in parameter
            if ( exists $parameter_href->{$parameter}{$key} ) {

                ## Check key data type
                $error_msg = _check_parameter_data_type(
                    {
                        parameter_href => $parameter_href,
                        key_href       => $key_href,
                        parameter      => $parameter,
                        key            => $key,
                        file_path      => $file_path,
                    }
                );
                if ($error_msg) {

                    return $error_msg;
                }

                ## Evaluate key values
                $error_msg = _check_parameter_values(
                    {
                        parameter_href => $parameter_href,
                        key_href       => $key_href,
                        parameter      => $parameter,
                        key            => $key,
                        file_path      => $file_path,
                    }
                );
                if ($error_msg) {

                    return $error_msg;
                }
            }
        }
    }
    return;
}

sub _check_parameter_values {

## Function : Evaluate parameter key values
## Returns  :
## Arguments: $parameter_href => Hash with parameters from yaml file {REF}
##          : $key_href       => Hash with  key {REF}
##          : $key            => Hash with non  key
##          : $file_path      => Path to yaml file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $key_href;
    my $parameter;
    my $key;
    my $file_path;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        key_href => {
            required    => 1,
            default     => {},
            strict_type => 1,
            store       => \$key_href,
        },
        parameter => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$parameter,
        },
        key =>
          { required => 1, defined => 1, strict_type => 1, store => \$key, },
        file_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_path,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Check value(s)
    if ( $key_href->{$key}{values} ) {

        my $value = $parameter_href->{$parameter}{$key};

        if ( not( any { $_ eq $value } @{ $key_href->{$key}{values} } ) ) {

            my $error_msg =
                q{Found illegal value '}
              . $value
              . q{' for parameter: '}
              . $parameter
              . q{' in key: '}
              . $key
              . q{' in file: '}
              . $file_path . $SINGLE_QUOTE . $NEWLINE
              . q{Allowed entries: '}
              . join( q{', '}, @{ $key_href->{$key}{values} } ) . $SINGLE_QUOTE . $NEWLINE;
            return $error_msg;
        }
    }
    return;
}

sub _check_parameter_data_type {

## Function : Check key data type
## Returns  :
## Arguments: $parameter_href => Hash with paremters from yaml file {REF}
##          : $key_href       => Hash with key {REF}
##          : $key            => Hash with non key
##          : $file_path      => Path to yaml file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $key_href;
    my $parameter;
    my $key;
    my $file_path;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        key_href => {
            required    => 1,
            default     => {},
            strict_type => 1,
            store       => \$key_href
        },
        parameter => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$parameter
        },
        key =>
          { required => 1, defined => 1, strict_type => 1, store => \$key },
        file_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_path
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Check data_type
    my $data_type = ref( $parameter_href->{$parameter}{$key} );

    ## Array or hash
    if ($data_type) {

        ## Wrong data_type
        if ( not $data_type eq $key_href->{$key}{key_data_type} ) {

            my $error_msg =
                q{Found '}
              . $data_type
              . q{' but expected datatype '}
              . $key_href->{$key}{key_data_type}
              . q{' for parameter: '}
              . $parameter
              . q{' in key: '}
              . $key
              . q{' in file: '}
              . $file_path . $SINGLE_QUOTE . $NEWLINE;
            return $error_msg;
        }
    }
    elsif ( $key_href->{$key}{key_data_type} ne q{SCALAR} ) {

        ## Wrong data_type
        my $error_msg =
            q{Found 'SCALAR' but expected datatype '}
          . $key_href->{$key}{key_data_type}
          . q{' for parameter: '}
          . $parameter
          . q{' in key: '}
          . $key
          . q{' in file: '}
          . $file_path . $SINGLE_QUOTE . $NEWLINE;
        return $error_msg;
    }
    return;
}

sub check_allowed_array_values {

## Function : Check that the array values are allowed
## Returns  :
## Arguments: $allowed_values_ref, $values_ref
##          : $allowed_values_ref => Allowed values for parameter {REF}
##          : $values_ref         => Values for parameter {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $allowed_values_ref;
    my $values_ref;

    my $tmpl = {
        allowed_values_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$allowed_values_ref
        },
        values_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$values_ref
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
## Arguments: $temp_directory
##          : $temp_directory => Temp directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $temp_directory;

    my $tmpl = {
        temp_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$temp_directory
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

1;
