package MIP::Recipes::Install::Post_installation;

use 5.026;
use Carp;
use Cwd qw{ abs_path };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile catdir };
use FindBin qw{ $Bin };
use Params::Check qw{ allow check last_error };
use charnames qw{ :full :short };
use open qw{ :encoding(UTF-8) :std };
use utf8;
use warnings qw{ FATAL utf8 };
use warnings;

## CPANM
use autodie qw{ :all };

## MIPs lib/
use MIP::Constants qw{ $COLON $CONTAINER_MANAGER $LOG_NAME $NEWLINE $SPACE $TAB };
use MIP::Environment::Child_process qw{ child_process };
use MIP::Environment::Container qw{ build_container_cmd };
use MIP::Environment::Executable
  qw{ get_binary_version_cmd get_binary_version get_executable };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.14;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      check_mip_installation
    };
}

sub check_mip_installation {

## Function : Write installation check oneliner to open filehandle
## Returns  :
## Arguments: $active_parameter_href => Active parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my %container_cmd = build_container_cmd(
        {
            active_parameter_href => $active_parameter_href,
            container_href        => $active_parameter_href->{container},
            container_manager     => $CONTAINER_MANAGER,
        }
    );

    my %capture_version_cmd = get_executable( {} );

  CONTAINER:
    while ( my ( $container, $container_href ) =
        each %{ $active_parameter_href->{container} } )
    {

        $log->info( q{Testing programs in image} . $COLON . $SPACE . $container );

      PROGRAM:
        foreach my $program ( keys %{ $container_href->{executable} } ) {

            if ( not $capture_version_cmd{$program} ) {

                $log->warn(
                    $TAB . q{No test available for} . $COLON . $SPACE . $program );
                next PROGRAM;
            }

            my %process_return = get_binary_version(
                {
                    binary                   => $program,
                    binary_path              => $container_cmd{$program},
                    capture_version_cmd_href => $capture_version_cmd{$program},
                    return_process_hash      => 1,
                }
            );

            if ( @{ $process_return{stdouts_ref} } ) {

                $log->info( $TAB
                      . q{Successfully launched}
                      . $COLON
                      . $SPACE
                      . $program
                      . $SPACE
                      . $process_return{stdouts_ref}[0] );
            }
            else {

                $log->warn( $TAB . q{Failed to launch} . $COLON . $SPACE . $program );
                $log->warn(
                    $TAB . q{Launch command returned} . $COLON . $NEWLINE . join $NEWLINE,
                    @{ $process_return{stderrs_ref} }
                );
            }
        }
    }
    return;
}

1;
