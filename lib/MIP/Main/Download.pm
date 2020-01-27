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
use MIP::Active_parameter qw{ update_to_absolute_path };
use MIP::Check::Download qw{ check_user_reference };
use MIP::Check::Parameter qw{ check_recipe_mode };
use MIP::Cluster qw{ check_max_core_number };
use MIP::Config qw{ check_cmd_config_vs_definition_file set_config_to_active_parameters };
use MIP::Constants
  qw{ $COLON $COMMA $DOT $MIP_VERSION $NEWLINE $SINGLE_QUOTE $SPACE $UNDERSCORE };
use MIP::File::Format::Yaml qw{ load_yaml };
use MIP::Log::MIP_log4perl qw{ get_log };
use MIP::Parameter qw{
  parse_parameter_files
  set_cache
  set_default
};
use MIP::Parse::Parameter qw{ parse_download_reference_parameter };
use MIP::Recipes::Check qw{ check_recipe_exists_in_hash };
use MIP::Recipes::Pipeline::Download_rd_dna qw{ pipeline_download_rd_dna };
use MIP::Recipes::Pipeline::Download_rd_rna qw{ pipeline_download_rd_rna };
use MIP::Update::Recipes qw{ update_recipe_mode_with_dry_run_all };
use MIP::User qw{ check_email_address };

BEGIN {
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.14;

    # Functions and variables that can be optionally exported
    our @EXPORT_OK = qw{ mip_download };

}

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
        my %config_parameter =
          load_yaml( { yaml_file => $active_parameter{config_file}, } );

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

    ## Parameters that have keys as MIP recipe names
    my @parameter_keys_to_check = (qw{ recipe_time recipe_core_number });
  PARAMETER_NAME:
    foreach my $parameter_name (@parameter_keys_to_check) {

        ## Test if key from query hash exists truth hash
        check_recipe_exists_in_hash(
            {
                parameter_name => $parameter_name,
                query_ref      => \%{ $active_parameter{$parameter_name} },
                truth_href     => \%parameter,
            }
        );
    }

    ## Parameters with key(s) that have elements as MIP recipe names
    my @parameter_element_to_check = qw{ associated_recipe };
  PARAMETER:
    foreach my $parameter ( keys %parameter ) {

      KEY:
        foreach my $parameter_name (@parameter_element_to_check) {

            next KEY if ( not exists $parameter{$parameter}{$parameter_name} );

            ## Test if element from query array exists truth hash
            check_recipe_exists_in_hash(
                {
                    parameter_name => $parameter_name,
                    query_ref      => \@{ $parameter{$parameter}{$parameter_name} },
                    truth_href     => \%parameter,
                }
            );
        }
    }

    ## Check that the module core number do not exceed the maximum per node
    foreach my $recipe_name ( keys %{ $active_parameter{recipe_core_number} } ) {

        ## Limit number of cores requested to the maximum number of cores available per node
        $active_parameter{recipe_core_number}{$recipe_name} = check_max_core_number(
            {
                max_cores_per_node => $active_parameter{max_cores_per_node},
                core_number_requested =>
                  $active_parameter{recipe_core_number}{$recipe_name},
            }
        );
    }

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
            log                   => $log,
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

    my $pipeline_type = $active_parameter{download_pipeline_type};

    ## Create dispatch table of pipelines
    my %pipeline = (
        rd_dna => \&pipeline_download_rd_dna,
        rd_rna => \&pipeline_download_rd_rna,
    );

    $log->info( q{Pipeline download type: } . $pipeline_type );
    $pipeline{$pipeline_type}->(
        {
            active_parameter_href => \%active_parameter,
            temp_directory        => $active_parameter{temp_directory},
        }
    );
    return;
}

##Investigate potential autodie error
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
