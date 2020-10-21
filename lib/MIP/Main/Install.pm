package MIP::Main::Install;

use 5.026;
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use Time::Piece;
use charnames qw{ :full :short };
use open qw{ :encoding(UTF-8) :std };
use utf8;
use warnings qw{ FATAL utf8 };
use warnings;

## MIPs lib/
use MIP::Active_parameter qw{ update_to_absolute_path };
use MIP::Config qw{ check_cmd_config_vs_definition_file set_config_to_active_parameters };
use MIP::Constants qw{ $COLON $MIP_VERSION $SPACE };
use MIP::Io::Read qw{ read_from_file };
use MIP::Log::MIP_log4perl qw{ get_log };
use MIP::Mip_install qw{ check_active_installation_parameters };
use MIP::Parameter qw{ set_default };
use MIP::Pipeline qw{ run_install_pipeline };
use MIP::Set::Parameter qw{ set_conda_path };

BEGIN {
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 2.18;

    # Functions and variables that can be optionally exported
    our @EXPORT_OK = qw{ mip_install };

}

sub mip_install {

## Function : Main script for generating MIP install scripts
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $parameter_href        => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
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

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments};

    ## Get local time
    my $date_time       = localtime;
    my $date_time_stamp = $date_time->datetime;
    my $date            = $date_time->ymd;

    ## Change relative path to absolute path for parameter with "update_path: absolute_path" in config
    update_to_absolute_path(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
        }
    );

    ### Config file
    ## Loads a YAML file into an arbitrary hash and returns it.
    my %config_parameter = read_from_file(
        {
            format => q{yaml},
            path   => $active_parameter_href->{config_file},
        }
    );

    ## Set config parameters into %active_parameter unless $parameter
    ## has been supplied on the command line
    set_config_to_active_parameters(
        {
            active_parameter_href => $active_parameter_href,
            config_parameter_href => \%config_parameter,
        }
    );

    ## Compare keys from config and cmd (%active_parameter) with definitions file (%parameter)
    check_cmd_config_vs_definition_file(
        {
            active_parameter_href => $active_parameter_href,
            parameter_href        => $parameter_href,
        }
    );

    ## Get log object and set log file in active parameters unless already set from cmd
    my $log = get_log(
        {
            active_parameter_href => $active_parameter_href,
            date                  => $date,
            date_time_stamp       => $date_time_stamp,
            log_name              => q{MIP_INSTALL},
            script                => q{mip_install},
        }
    );
    $log->info( q{MIP Version: } . $MIP_VERSION );
    $log->info( q{Writing log messages to}
          . $COLON
          . $SPACE
          . $active_parameter_href->{log_file} );

    ## Set default from parameter hash to active_parameter for uninitilized parameters
    set_default(
        {
            active_parameter_href => $active_parameter_href,
            custom_default_parameters_ref =>
              $parameter_href->{custom_default_parameters}{default},
            parameter_href => $parameter_href,
        }
    );

    ## Installation specific checks of active_parameter hash
    check_active_installation_parameters(
        {
            project_id  => $active_parameter_href->{project_id},
            sbatch_mode => $active_parameter_href->{sbatch_mode},
        }
    );

    ## Set path to conda
    set_conda_path(
        {
            active_parameter_href => $active_parameter_href,
        }
    );

    run_install_pipeline(
        {
            active_parameter_href => $active_parameter_href,
        }
    );

    return;
}

1;
