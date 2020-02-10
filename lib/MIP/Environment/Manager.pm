package MIP::Environment::Manager;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ get_env_method_cmds };
}

sub get_env_method_cmds {

## Function : Get the standard load and unload env command for environment method
## Returns  : @env_method_cmds
## Arguments: $action     => What to do with the environment
##          : $env_method => Method used to load environment
##          : $env_name   => Name of environment

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $action;
    my $env_method;
    my $env_name;

    my $tmpl = {
        action => {
            allow       => [qw{ load unload }],
            defined     => 1,
            required    => 1,
            store       => \$action,
            strict_type => 1,
        },
        env_method => {
            allow       => [qw{ conda }],
            defined     => 1,
            required    => 1,
            store       => \$env_method,
            strict_type => 1,
        },
        env_name => {
            required    => 1,
            store       => \$env_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Package_manager::Conda qw{ conda_activate conda_deactivate };

    my %method_cmd = (
        conda => {
            load   => [ ( conda_activate(   { env_name => $env_name, } ), ) ],
            unload => [ ( conda_deactivate( {} ), ) ],
        },
    );
    return ( @{ $method_cmd{$env_method}{$action} } );
}
1;
