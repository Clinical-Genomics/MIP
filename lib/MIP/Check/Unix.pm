package MIP::Check::Unix;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

BEGIN {

    use base qw{Exporter};
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.07;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_binary_in_path };
}

sub check_binary_in_path {

## Function : Scans through PATH for supplied binary
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis {REF}
##          : $binary                => Binary to search for
##          : $log                   => Log
##          : $program_name          => MIP program name (Analysis recipe switch)

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $binary;
    my $log;
    my $program_name;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        binary => {
            defined     => 1,
            required    => 1,
            store       => \$binary,
            strict_type => 1,
        },
        log => {
            store => \$log,
        },
        program_name => {
            store       => \$program_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Environment::Path qw{ is_binary_in_path };
    use MIP::Get::Parameter qw{ get_dynamic_conda_path };

    ## Search for binary in PATH in any MIP conda env defined by config
    ## or conda base
    my $env_binary_path = get_dynamic_conda_path(
        {
            active_parameter_href => $active_parameter_href,
            bin_file              => $binary,
            environment_key       => $program_name,
        }
    );

    ## Add binary for use downstream
    $active_parameter_href->{binary_path}{$binary} = $binary;

    ## Potential full path to binary
    my $binary_path = catfile( $env_binary_path, $binary );

    ## Search for binary in conda envs path
    if ( -e $binary_path ) {

        ## Update binary to full binary_path
        $active_parameter_href->{binary_path}{$binary} = $binary_path;

        $binary = $binary_path;
    }

    ## Test binary
    is_binary_in_path( { binary => $binary, } );

    return 1;
}

1;
