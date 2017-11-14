package MIP::Set::Parameter;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile splitpath };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use Readonly;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ set_config_to_active_parameters set_dynamic_parameter set_parameter_reference_dir_path set_parameter_to_broadcast };
}

## Constants
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };
Readonly my $TAB     => qq{\t};

sub set_config_to_active_parameters {

## Function : Add contig parameters to active_parameters if not already initilized from command line
## Returns  :
## Arguments: $config_parameter_href => Config parameters hash
##         : $active_parameter_href  => Active parameters for this analysis hash {REF}

    my ($arg_href) = @_;

## Flatten argument(s)
    my $config_parameter_href;
    my $active_parameter_href;

    my $tmpl = {
        config_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$config_parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
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
          if ( $active_parameter_href->{$parmeter_name}
            && ref( $active_parameter_href->{$parmeter_name} ) !~
            / HASH | ARRAY /sxm );

        ### No input from cmd
        ## Add to active_parameter
        $active_parameter_href->{$parmeter_name} =
          $config_parameter_href->{$parmeter_name};
    }
    return;
}

sub set_dynamic_parameter {

## Function : Sets dynamic aggregate information from definitions to parameter hash
## Returns  :
## Arguments: $parameter_href => Parameter hash {REF}
##          : $aggregates_ref => The data to aggregate and add to parameter hash{REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $aggregates_ref;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        aggregates_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$aggregates_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Constants
    Readonly my $RECORD_SEPARATOR => q{:};
    Readonly my $FIELD_COUNTER    => 2;

  PARAMETER:
    foreach my $parameter_name ( keys %{$parameter_href} ) {

      KEY_AND_STRING_TO_MATCH:
        foreach my $aggregate_element ( @{$aggregates_ref} ) {

            ## Split into key and string to match
            my ( $second_key, $string_to_match, $unexpected_data ) =
              split $RECORD_SEPARATOR, $aggregate_element, $FIELD_COUNTER + 1;

            ## Make sure that we get what we expect
            if ( defined $unexpected_data ) {

                carp
                  q{Unexpected trailing garbage at end of aggregate_element '}
                  . $aggregate_element
                  . q{':}, $NEWLINE . $TAB . $unexpected_data . $NEWLINE;
            }

            if ( defined $parameter_href->{$parameter_name}{$second_key}
                && $parameter_href->{$parameter_name}{$second_key} eq
                $string_to_match )
            {

                push @{ $parameter_href->{dynamic_parameter}{$string_to_match}
                  },
                  $parameter_name;
            }
        }
    }
    return;
}

sub set_parameter_reference_dir_path {

## Function : Set path for supplied reference(s) associated with parameter that should reside in the mip reference directory to full path.
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $parameter_name        => Parameter to update

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_name;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        parameter_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$parameter_name,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Unpack
    my $reference_dir = $active_parameter_href->{reference_dir};

    # $parameter can be array_ref, hash_ref, point to file or undef
    my $parameter = $active_parameter_href->{$parameter_name};

    return if ( not defined $parameter );

    if ( ref $parameter eq q{ARRAY} ) {

      FILE:
        foreach my $file ( @{$parameter} ) {

            ## Split to restate
            my ( $volume, $directory, $file_name ) = splitpath($file);

            ## Update original element - works since array_ref
            $file = catfile( $reference_dir, $file_name );
        }
        return;
    }
    elsif ( ref $parameter eq q{HASH} ) {

      FILE:
        foreach my $file ( keys %{$parameter} ) {

            ## Split to restate
            my ( $volume, $directory, $file_name ) = splitpath($file);

            ## Update original key with path and add potential annotation key
            ## by deleting original value (returns value deleted)
            $active_parameter_href->{$parameter_name}
              { catfile( $reference_dir, $file_name ) } =
              delete $active_parameter_href->{$parameter_name}{$file};
        }
        return;
    }
    else {

        ## File
        ## Split to restate
        my ( $volume, $directory, $file_name ) =
          splitpath( $active_parameter_href->{$parameter_name} );

        ## Restate to allow for changing mip reference directory between runs
        $active_parameter_href->{$parameter_name} = $file_name;

        ## Update original value
        my $path =
          catfile( $reference_dir, $active_parameter_href->{$parameter_name} );
        $active_parameter_href->{$parameter_name} = $path;

    }
    return;
}

sub set_parameter_to_broadcast {

## Function : Set parameters to broadcast message
## Returns  :
## Arguments: $parameter_href        => Holds all parameters
##          : $active_parameter_href => Active parameters for this analysis hash {REF}
    ##          : $order_parameters_ref  => Order of parameters (for structured output) {REF}
##          : $broadcasts_ref        => Holds the parameters info for broadcasting later {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $order_parameters_ref;
    my $broadcasts_ref;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        order_parameters_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$order_parameters_ref,
        },
        broadcasts_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$broadcasts_ref,
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

            ## Alias
            my $element_separator =
              $parameter_href->{$parameter_name}{element_separator};

            $info .= join $element_separator,
              @{ $active_parameter_href->{$parameter_name} };

            ## Add info to broadcasts
            push @{$broadcasts_ref}, $info;
        }
        elsif ( ref $active_parameter_href->{$parameter_name} eq q{HASH} ) {

            $info .= join q{,},
              map { qq{$_=$active_parameter_href->{$parameter_name}{$_}} }
              ( keys %{ $active_parameter_href->{$parameter_name} } );

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

1;
