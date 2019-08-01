package MIP::Check::Unix;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use File::Spec::Functions qw{ catfile };
use IPC::Cmd qw{ can_run };
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
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_binary_in_path };
}

sub check_binary_in_path {

## Function  : Scans through PATH for supplied binary
## Returns   :
## Arguments : $active_parameter_href => Holds all set parameter for analysis {REF}
##           : $binary                => Binary to search for
##           : $log                   => Log
##           : $program_name          => MIP program name (Analysis recipe switch)

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

    use MIP::Get::Parameter qw{ get_dynamic_conda_path };

    ## Check if the binary has a singularity image
    if ( can_run( $active_parameter_href->{singularity_container}{$binary} ) ) {

        ## Broadcast successful scan through PATH for supplied binary
        _check_binary_broadcast_pass(
            {
                binary => catfile($binary),
                log    => $log,
            }
        );
        return 1;
    }

    ## Search for binary in PATH in any MIP conda env defined by config
    ## or conda base
    my $env_binary_path = get_dynamic_conda_path(
        {
            active_parameter_href => $active_parameter_href,
            bin_file              => $binary,
            environment_key       => $program_name,
        }
    );

    ## Search for binary in conda envs path
    if ( can_run( catfile( $env_binary_path, $binary ) ) ) {

        ## Broadcast successful scan through PATH for supplied binary
        _check_binary_broadcast_pass(
            {
                binary => catfile( $env_binary_path, $binary ),
                log    => $log,
            }
        );
        return 1;
    }

    # Search for binary in PATH
    if ( can_run($binary) ) {

        ## Broadcast successful scan through PATH for supplied binary
        _check_binary_broadcast_pass(
            {
                binary => $binary,
                log    => $log,
            }
        );
        return 1;
    }

    ## Broadcast scan through PATH for supplied binary when not found
    _check_binary_broadcast_fail(
        {
            binary => $binary,
            log    => $log,
        }
    );
    exit 1;
}

sub _check_binary_broadcast_fail {

## Function  : Broadcast scan through PATH for supplied binary when not found
## Returns   :
## Arguments : $binary => Binary to search for
##           : $log    => Log

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $binary;
    my $log;

    my $tmpl = {
        binary => {
            defined     => 1,
            required    => 1,
            store       => \$binary,
            strict_type => 1,
        },
        log => {
            store => \$log,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Broadcast binary not found
    if ($log) {

        $log->fatal( q{Could not detect } . $binary . q{ in PATH} );
    }
    else {

        say {*STDERR} q{Could not detect } . $binary . q{ in PATH};
    }
    return;
}

sub _check_binary_broadcast_pass {

## Function  : Broadcast successful scan through PATH for supplied binary
## Returns   :
## Arguments : $binary => Binary to search for
##           : $log    => Log

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $binary;
    my $log;

    my $tmpl = {
        binary => {
            defined     => 1,
            required    => 1,
            store       => \$binary,
            strict_type => 1,
        },
        log => {
            store => \$log,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Broadcast if found
    if ( defined $log ) {

        $log->info( q{Program check: } . $binary . q{ in PATH} );
    }
    else {

        say {*STDERR} q{Program check: } . $binary . q{ in PATH};
    }
    return 1;

}

1;
