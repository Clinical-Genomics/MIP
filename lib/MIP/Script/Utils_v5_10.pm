package MIP::Script::Utils_v5_10;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };

BEGIN {
    use base qw{ Exporter };
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ help set_default_array_parameters };

}

sub help {

## help

## Function  : Print help text and exit with supplied exit code
## Returns   : ""
## Arguments : $USAGE, $exit_code
##           : $USAGE     => Help text
##           : $exit_code => Exit code

    my ($arg_href) = @_;

    my $USAGE     = $arg_href->{USAGE};
    my $exit_code = $arg_href->{exit_code};

    print STDOUT $USAGE, "\n";
    exit $exit_code;
}

sub set_default_array_parameters {

## set_default_array_parameters

## Function  : Set default for array parameters unless parameter already exists in parameter hash
## Returns   : ""
## Arguments : $parameter_href, $array_parameter_href
##           : $parameter_href       => Parameters hash {REF}
##           : $array_parameter_href => Hold the array parameter defaults as {REF}

    my ($arg_href) = @_;

    my $parameter_href       = $arg_href->{parameter_href};
    my $array_parameter_href = $arg_href->{array_parameter_href};

  PARAMETER_NAME:
    foreach my $parameter_name ( keys %{$array_parameter_href} ) {
        if ( not @{ $parameter_href->{$parameter_name} } ) {
            $parameter_href->{$parameter_name} =
              $array_parameter_href->{$parameter_name};
        }
    }

    return;
}

1;
