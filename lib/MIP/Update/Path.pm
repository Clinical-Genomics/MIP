package MIP::Update::Path;

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

## MIPs lib/
use MIP::Constants qw{ $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ update_to_absolute_path };
}

sub update_to_absolute_path {

## Function : Change relative path to absolute path for parameters with update_path: aboslute path in definitions file
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $parameter_href        => Parameter hash {REF}

    my ($arg_href) = @_;

    ##Flatten argument(s)
    my $active_parameter_href;
    my $parameter_href;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Set::File qw{ set_absolute_path };
    use MIP::Parameter qw{ set_cache };

    ## Adds dynamic aggregate information from definitions to parameter hash
    # Collect all path that should be made absolute
    set_cache(
        {
            parameter_href => $parameter_href,
            aggregates_ref => [q{update_path:absolute_path}],
        }
    );

  DYNAMIC_PARAMETER:
    foreach my $parameter_name ( @{ $parameter_href->{cache}{absolute_path} } ) {

        ## If array
        if ( ref $active_parameter_href->{$parameter_name} eq q{ARRAY} ) {

            foreach my $parameter_value ( @{ $active_parameter_href->{$parameter_name} } )
            {

                $parameter_value = set_absolute_path(
                    {
                        path           => $parameter_value,
                        parameter_name => $parameter_name,
                    }
                );
            }
        }
        elsif ( ref $active_parameter_href->{$parameter_name} eq q{HASH} ) {
            ## Hash

            ## Alias
            my $parameter_value_ref = $active_parameter_href->{$parameter_name};

            ## Cannot use each since we are updating key
            foreach my $key ( keys %{$parameter_value_ref} ) {

                if ( defined $parameter_value_ref ) {

                    ## Find aboslute path for supplied key path or croaks and exists if path does not exists
                    my $updated_key = set_absolute_path(
                        {
                            path           => $key,
                            parameter_name => $parameter_name,
                        }
                    );

                    $parameter_value_ref->{$updated_key} =
                      delete( $parameter_value_ref->{$key} );
                }
            }
        }
        elsif ( exists $active_parameter_href->{$parameter_name}
            && defined $active_parameter_href->{$parameter_name} )
        {
            ## Scalar

            $active_parameter_href->{$parameter_name} = set_absolute_path(
                {
                    path           => $active_parameter_href->{$parameter_name},
                    parameter_name => $parameter_name,
                }
            );
        }
    }
    return;
}

1;
