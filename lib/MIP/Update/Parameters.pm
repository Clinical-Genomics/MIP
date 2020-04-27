package MIP::Update::Parameters;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
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
    our $VERSION = 1.09;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ update_dynamic_config_parameters update_with_dynamic_config_parameters };
}

sub update_dynamic_config_parameters {

## Function : Updates a scalar reference with dynamic parameter values. Leaves other entries untouched.
## Returns  :
## Arguments: $active_parameter_ref   => Active parameter for this analysis {REF}
##          : $dynamic_parameter_href => Map of dynamic parameters

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_ref;
    my $dynamic_parameter_href;

    my $tmpl = {
        active_parameter_ref => {
            required => 1,
            store    => \$active_parameter_ref,
        },
        dynamic_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$dynamic_parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  DYNAMIC_PARAMETER:
    while ( my ( $dynamic_parameter_name, $dynamic_parameter_value ) =
        each %{$dynamic_parameter_href} )
    {

        ${$active_parameter_ref} =~
          s/$dynamic_parameter_name!/$dynamic_parameter_value/xsmgi;
    }
    return;
}

sub update_with_dynamic_config_parameters {

## Function : Updates the config file to particular user/cluster with dynamic config parameters following specifications. Leaves other entries untouched.
## Returns  :
## Arguments: $active_parameter_href  => Active parameters for this analysis hash {REF}
##          : $dynamic_parameter_href => Map of dynamic parameters

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $dynamic_parameter_href;

    my $tmpl = {
        active_parameter_href => {
            required => 1,
            store    => \$active_parameter_href,
        },
        dynamic_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$dynamic_parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Return if variable isn't in use
    return if ( not defined $active_parameter_href );

    ## HASH
    if ( ref $active_parameter_href eq q{HASH} ) {

      KEY:
        foreach my $key ( keys %{$active_parameter_href} ) {

            next KEY if ( not defined $active_parameter_href->{$key} );

            if (   ref $active_parameter_href->{$key} eq q{HASH}
                or ref $active_parameter_href->{$key} eq q{ARRAY} )
            {

                update_with_dynamic_config_parameters(
                    {
                        active_parameter_href  => $active_parameter_href->{$key},
                        dynamic_parameter_href => $dynamic_parameter_href,
                    }
                );
            }

            ## Update parameter
            update_dynamic_config_parameters(
                {
                    active_parameter_ref   => \$active_parameter_href->{$key},
                    dynamic_parameter_href => $dynamic_parameter_href,
                }
            );
        }
    }
    ##  ARRAY
    elsif ( ref $active_parameter_href eq q{ARRAY} ) {

      ELEMENT:
        while ( my ( $element_index, $element ) = each @{$active_parameter_href} ) {

            next ELEMENT if ( not defined $element );

            if ( ref $element eq q{HASH} or ref $element eq q{ARRAY} ) {

                update_with_dynamic_config_parameters(
                    {
                        active_parameter_href  => $element,
                        dynamic_parameter_href => $dynamic_parameter_href,
                    }
                );
            }

            ## Update parameter
            update_dynamic_config_parameters(
                {
                    active_parameter_ref   => \$active_parameter_href->[$element_index],
                    dynamic_parameter_href => $dynamic_parameter_href,
                }
            );
        }
    }
    return;
}

1;
