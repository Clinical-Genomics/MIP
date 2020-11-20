package MIP::Recipes::Install::Container;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ devnull };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

## MIPs lib/
use MIP::Constants qw{ $COLON $CONTAINER_MANAGER $LOG_NAME $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_containers };
}

sub install_containers {

## Function : Setup containers to use with docker or singularity
## Returns  :
## Arguments: $active_parameter_href => Active parameter hash {REF}
##          : $container_href        => Hash with container {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $container_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        container_href => {
            default     => {},
            required    => 1,
            store       => \$container_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Environment::Child_process qw{ child_process };
    use MIP::Environment::Container qw{ parse_container_uri run_container };
    use MIP::Recipes::Install::Vep qw{ install_vep };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Set verbosity
    my $stderr_file_path = $active_parameter_href->{quiet} ? devnull() : undef;

    ## Containers requiring something extra
    my %finish_container_installation = ( vep => \&install_vep, );

  CONTAINER:
    foreach my $container ( keys %{$container_href} ) {

        $log->info( q{Caching image} . $COLON . $SPACE . $container );

        parse_container_uri(
            {
                container_manager => $CONTAINER_MANAGER,
                uri_ref           => \$container_href->{$container}{uri},
            }
        );

        ## Get command for caching image
        my @container_commands = run_container(
            {
                container_cmds_ref => [qw{ ls }],
                container_path     => $container_href->{$container}{uri},
                container_manager  => $CONTAINER_MANAGER,
                stderrfile_path    => $stderr_file_path,
                stdoutfile_path    => devnull(),
            }
        );

        my %process_return = child_process(
            {
                commands_ref => [ join $SPACE, @container_commands ],
                process_type => q{ipc_cmd_run},
            }
        );

        if ( not $process_return{success} ) {

            $log->fatal(qq{$CONTAINER_MANAGER failed to cache $container});
            $log->logdie( $process_return{error_messages_ref}[0] );
        }

        ## Finishing touches for certain containers
        if ( $finish_container_installation{$container} ) {

            $finish_container_installation{$container}->(
                {
                    active_parameter_href => $active_parameter_href,
                    container_href        => $container_href->{$container},
                }
            );
        }
    }
    return 1;
}

1;
