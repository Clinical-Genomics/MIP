package MIP::Main::Install;

use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Basename qw{ dirname basename fileparse };
use File::Spec::Functions qw{ catfile catdir devnull };
use Getopt::Long;
use IO::Handle;
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use Time::Piece;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use Readonly;

## MIPs lib/
use MIP::Active_parameter qw{ update_to_absolute_path };
use MIP::Check::Parameter qw{ check_active_installation_parameters };
use MIP::Config qw{ check_cmd_config_vs_definition_file set_config_to_active_parameters };
use MIP::Constants
  qw{ $COLON $COMMA $DOT $MIP_VERSION $NEWLINE $SINGLE_QUOTE $SPACE $UNDERSCORE };
use MIP::Io::Read qw{ read_from_file };
use MIP::Log::MIP_log4perl qw{ get_log };
use MIP::Parameter qw{ set_default };
use MIP::Set::Parameter qw{ set_conda_path };

## Recipes
use MIP::Recipes::Pipeline::Install_rd_dna qw{ pipeline_install_rd_dna };
use MIP::Recipes::Pipeline::Install_rd_rna qw{ pipeline_install_rd_rna };

## Constants
Readonly my $THREE     => 3;
Readonly my $MINUS_ONE => -1;

BEGIN {
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 2.12;

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

    ## Transfer to lexical variables
    # Parameters to include in each download run
    my %active_parameter = %{$active_parameter_href};

    # All parameters MIP install knows
    my %parameter = %{$parameter_href};

    ## Get local time
    my $date_time       = localtime;
    my $date_time_stamp = $date_time->datetime;
    my $date            = $date_time->ymd;

    # Catches name of current script
    my $script = _this_sub();

    ## Catches name of the calling module
    my $process = _parent_module();

    ## Build pipeline name
    my $pipeline = q{install} . $UNDERSCORE . lc $process;

    ## Change relative path to absolute path for parameter with "update_path: absolute_path" in config
    update_to_absolute_path(
        {
            active_parameter_href => \%active_parameter,
            parameter_href        => \%parameter,
        }
    );

    ### Config file
    ## Loads a YAML file into an arbitrary hash and returns it.
    my %config_parameter = read_from_file(
        {
            format => q{yaml},
            path   => $active_parameter{config_file},
        }
    );

    ## Set config parameters into %active_parameter unless $parameter
    ## has been supplied on the command line
    set_config_to_active_parameters(
        {
            active_parameter_href => \%active_parameter,
            config_parameter_href => \%config_parameter,
        }
    );

    ## Compare keys from config and cmd (%active_parameter) with definitions file (%parameter)
    check_cmd_config_vs_definition_file(
        {
            active_parameter_href => \%active_parameter,
            parameter_href        => \%parameter,
        }
    );

## Get log object and set log file in active parameters unless already set from cmd
    my $log = get_log(
        {
            active_parameter_href => \%active_parameter,
            date                  => $date,
            date_time_stamp       => $date_time_stamp,
            log_name              => uc $script,
            script                => $script,
        }
    );
    $log->info( q{MIP Version: } . $MIP_VERSION );
    $log->info(
        q{Writing log messages to} . $COLON . $SPACE . $active_parameter{log_file} );

    ## Set default from parameter hash to active_parameter for uninitilized parameters
    set_default(
        {
            active_parameter_href => \%active_parameter,
            custom_default_parameters_ref =>
              \@{ $parameter{custom_default_parameters}{default} },
            parameter_href => \%parameter,
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
            active_parameter_href => \%active_parameter,
        }
    );

    ## Store script, process and pipeline for broadcasting later
    $active_parameter{script}   = $script;
    $active_parameter{process}  = lc $process;
    $active_parameter{pipeline} = $pipeline;

    ## Create dispatch table of pipelines
    my %pipeline_table = (
        install_rd_dna => \&pipeline_install_rd_dna,
        install_rd_rna => \&pipeline_install_rd_rna,
    );

    $log->info( q{Pipeline type: } . $pipeline );
    $pipeline_table{$pipeline}->(
        {
            active_parameter_href => \%active_parameter,
            quiet                 => $active_parameter{quiet},
            verbose               => $active_parameter{verbose},
        }
    );
    return;
}

sub _this_sub {

## Function : Returns the name of the current subroutine
## Returns  : $this_sub
## Arguments:

    ## Get full path to subroutine
    my $this_sub = ( caller 1 )[$THREE];

    ## Isolate subroutine
    $this_sub = ( split /::/xms, $this_sub )[$MINUS_ONE];

    return $this_sub;
}

sub _parent_module {

## Function : Returns the name of the module that called this one
## Returns  : $parent_module
## Arguments:

    ## Get full path to module
    my $parent_module = ( caller 1 )[0];

    ## Isolate module
    $parent_module = ( split /::/xms, $parent_module )[$MINUS_ONE];

    return $parent_module;
}

1;
