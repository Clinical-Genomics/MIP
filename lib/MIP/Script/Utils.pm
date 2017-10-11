package MIP::Script::Utils;

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
use autodie;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ help set_default_array_parameters };
}

## Constants
Readonly my $DOT     => q{.};
Readonly my $NEWLINE => qq(\n);
Readonly my $SPACE   => q{ };

sub help {

## Function : Print help text and exit with supplied exit code
## Returns  :
## Arguments: $USAGE     => Help text
##          : $exit_code => Exit code

    my ($arg_href) = @_;

    ## Default(s)
    my $exit_code;

    ## Flatten argument(s)
    my $USAGE;

    my $tmpl = {
        USAGE =>
          { required => 1, defined => 1, strict_type => 1, store => \$USAGE },
        exit_code => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$exit_code,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    say {*STDOUT} $USAGE;
    exit $exit_code;
}

sub set_default_array_parameters {

## Function : Set default for array parameters unless parameter already exists in parameter hash
## Returns  :
## Arguments: $parameter_href       => Parameters hash {REF}
##          : $array_parameter_href => Hold the array parameter defaults as {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $array_parameter_href;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        array_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$array_parameter_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  PARAMETER:
    foreach my $parameter_name ( keys %{$array_parameter_href} ) {

        ## Unless parameter already exists
        if ( not @{ $parameter_href->{$parameter_name} } ) {

            $parameter_href->{$parameter_name} =
              $array_parameter_href->{$parameter_name}{default};
        }
    }
    return;
}

1;
