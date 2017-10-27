package MIP::Set::Parameter;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };

## CPANM
use Readonly;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ set_dynamic_parameter };
}

## Constants
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };
Readonly my $TAB     => qq{\t};

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

1;
