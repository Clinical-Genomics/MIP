package MIP::Main::Download;

use 5.026;
use Carp;
use charnames qw( :full :short );
use Cwd qw{ abs_path };
use English qw{ -no_match_vars };
use File::Basename qw{ basename fileparse };
use File::Path qw{ make_path };
use warnings qw{ FATAL utf8 };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use Time::Piece;
use utf8;

## CPANM
use autodie qw{ open close :all };
use List::Util qw{ any };
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use MIP::Active_parameter qw{
  check_recipe_mode
  parse_recipe_resources
  set_load_env_environment
  update_recipe_mode_with_dry_run_all
  update_to_absolute_path
};
use MIP::Check::Download qw{ check_user_reference };
use MIP::Config qw{ check_cmd_config_vs_definition_file set_config_to_active_parameters };
use MIP::Constants
  qw{ $COLON $COMMA $DOT $LOG_NAME $MIP_VERSION $NEWLINE $SINGLE_QUOTE $SPACE $UNDERSCORE };
use MIP::Environment::Cluster qw{ check_max_core_number };
use MIP::Environment::User qw{ check_email_address };
use MIP::Io::Read qw{ read_from_file };
use MIP::Log::MIP_log4perl qw{ get_log };
use MIP::Parameter qw{
  parse_parameter_files
  set_cache
  set_default
};
use MIP::Parse::Parameter qw{ parse_download_reference_parameter };
use MIP::Pipeline qw{ run_download_pipeline };
use MIP::Recipes::Check qw{ check_recipe_exists_in_hash };
use MIP::Recipes::Parse qw{ parse_recipes };

BEGIN {
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.20;

    # Functions and variables that can be optionally exported
    our @EXPORT_OK = qw{ mip_download };

}

## Constants
Readonly my %RECIPE_PARAMETERS_TO_CHECK => (
    keys     => [qw{ recipe_core_number recipe_time }],
    elements => [qw{ associated_recipe  }],
);

sub mip_download {

## Function : Main script for generating MIP download scripts
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

    # All parameters MIP download knows
    my %parameter = %{$parameter_href};

    ## Get local time
    my $date_time       = localtime;
    my $date_time_stamp = $date_time->datetime;
    my $date            = $date_time->ymd;

    # Catches script name and removes ending
    my $script = fileparse( basename( $PROGRAM_NAME, $DOT . q{pl} ) );

    # Add specific MIP process
    $script .= $UNDERSCORE . q{download};

    ## Change relative path to absolute path for parameter with "update_path: absolute_path" in config
    update_to_absolute_path(
        {
            active_parameter_href => \%active_parameter,
            parameter_href        => \%parameter,
        }
    );

    ### Config file
    ## If config from cmd
    if ( exists $active_parameter{config_file}
        && defined $active_parameter{config_file} )
    {

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
    }

## Get log object and set log file in active parameters unless already set from cmd
    my $log = get_log(
        {
            active_parameter_href => \%active_parameter,
            date                  => $date,
            date_time_stamp       => $date_time_stamp,
            log_name              => uc q{mip_download},
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

    ## Make sure that we have lower case from user input
    @{ $active_parameter{reference_genome_versions} } =
      map { lc } @{ $active_parameter{reference_genome_versions} };

    ## Create reference dir if it does not exists
    make_path( $active_parameter{reference_dir} );

    ### Checks

    ## Parse existence of files and directories
    parse_parameter_files(
        {
            active_parameter_href => \%active_parameter,
            parameter_href        => \%parameter,
        }
    );

    ## Check email adress syntax and mail host
    check_email_address(
        {
            email => $active_parameter{email},
        }
    );

    ## Parameters that have keys or elements as MIP recipe names
    parse_recipes(
        {
            active_parameter_href   => \%active_parameter,
            parameter_href          => \%parameter,
            parameter_to_check_href => \%RECIPE_PARAMETERS_TO_CHECK,
        }
    );

    ## Check core number requested against environment provisioned
    parse_recipe_resources( { active_parameter_href => \%active_parameter, } );

    ## Adds dynamic aggregate information from definitions to parameter hash
    set_cache(
        {
            aggregates_ref => [
                ## Collects all recipes that MIP can handle
                q{type:recipe},
            ],
            parameter_href => \%parameter,
        }
    );

    ## Check correct value for recipe mode in MIP
    check_recipe_mode(
        {
            active_parameter_href => \%active_parameter,
            parameter_href        => \%parameter,
        }
    );

    ## Update recipe mode depending on dry_run_all flag
    update_recipe_mode_with_dry_run_all(
        {
            active_parameter_href => \%active_parameter,
            dry_run_all           => $active_parameter{dry_run_all},
            recipes_ref           => \@{ $parameter{cache}{recipe} },
        }
    );

    set_load_env_environment( { active_parameter_href => \%active_parameter, } );

    ## Remodel depending on if "--reference" was used or not as the user info is stored as a scalar per reference_id while yaml is stored as arrays per reference_id
    parse_download_reference_parameter(
        { reference_href => \%{ $active_parameter{reference} }, } );

    check_user_reference(
        {
            user_supplied_reference_ref => \%{ $active_parameter{reference} },
            reference_genome_versions_ref =>
              \@{ $active_parameter{reference_genome_versions} },
            reference_ref => \%{ $active_parameter{reference_feature} },
        }
    );

    $log->info(
q{Will write sbatch install instructions for references to individual sbatch scripts}
    );

    run_download_pipeline( { active_parameter_href => \%active_parameter, } );

    return;
}

## Investigate potential autodie error
if ( $EVAL_ERROR and $EVAL_ERROR->isa(q{autodie::exception}) ) {

    if ( $EVAL_ERROR->matches(q{default}) ) {

        say {*STDERR} q{Not an autodie error at all};
    }
    if ( $EVAL_ERROR->matches(q{open}) ) {

        say {*STDERR} q{Error from open};
    }
    if ( $EVAL_ERROR->matches(q{:io}) ) {

        say {*STDERR} q{Non-open, IO error.};
    }
}
elsif ($EVAL_ERROR) {

    say {*STDERR} q{A non-autodie exception.};
}

1;
