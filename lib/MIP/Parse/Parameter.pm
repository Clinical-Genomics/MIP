package MIP::Parse::Parameter;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile splitpath };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.18;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      parse_conda_env_name
      parse_download_reference_parameter
    };

}

sub parse_conda_env_name {

## Function : Build conda environment names depending on input parameters
## Returns  : $conda_environment_name
## Arguments: $base_name     => Degfault base environment name
##          : $date          => Date
##          : $environment   => Installation environment
##          : parameter_href => Parmeter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $base_name;
    my $date;
    my $environment;
    my $parameter_href;

    my $tmpl = {
        base_name => {
            defined     => 1,
            required    => 1,
            store       => \$base_name,
            strict_type => 1,
        },
        date => {
            defined     => 1,
            required    => 1,
            store       => \$date,
            strict_type => 1,
        },
        environment => {
            defined     => 1,
            required    => 1,
            store       => \$environment,
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

    my $environment_name = $parameter_href->{environment_name}{$environment};

    ## Give the env a default name if not given
    if ( not $environment_name ) {

        ## Strip the first character, i.e. 'e' from the environment string
        my $env_postfix = substr $environment, 1;
        $environment_name = $base_name . $UNDERSCORE . $env_postfix;
    }

    ## Prepend environemnt prefix
    if ( $parameter_href->{environment_prefix} ) {

        $environment_name =
          $parameter_href->{environment_prefix} . $UNDERSCORE . $environment_name;
    }

    ## Add environment date
    if ( $parameter_href->{add_environment_date} ) {

        $environment_name = $environment_name . $UNDERSCORE . $date;
    }

    ## Append environment suffix
    if ( $parameter_href->{environment_suffix} ) {

        $environment_name =
          $environment_name . $UNDERSCORE . $parameter_href->{environment_suffix};
    }

    return $environment_name;
}

sub parse_download_reference_parameter {

## Function : Remodel depending on if "--reference" was used or not as the user info is stored as a scalar per reference_id while yaml is stored as arrays per reference_id
## Returns  :
## Arguments: $reference_href => Reference hash {REF}

    my ($arg_href) = @_;

## Flatten argument(s)
    my $reference_href;

    my $tmpl = {
        reference_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$reference_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  VERSION_REF:
    foreach my $versions_ref ( values %{$reference_href} ) {

        if ( ref $versions_ref ne q{ARRAY} ) {

            ## Make scalar from CLI '--ref key=value' option into array
            $versions_ref = [$versions_ref];
        }
    }

    return;
}

1;
