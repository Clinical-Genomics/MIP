package MIP::Recipes::Pipeline::Install;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use List::Util qw{ any };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## MIPs lib/
use MIP::Constants qw{
  $COLON
  $DOT
  $LOG_NAME
  $NEWLINE
  $SINGLE_QUOTE
  $SPACE
  set_container_constants
};
use MIP::Recipes::Install::Container qw{ install_containers };
use MIP::Recipes::Install::Mip_scripts qw{ install_mip_scripts };
use MIP::Recipes::Install::Post_installation qw{ check_mip_installation };
use MIP::Set::Parameter qw{ set_programs_for_installation };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ pipeline_install };
}

sub pipeline_install {

## Function : Install recipes for mip pipelines
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this download hash {REF}
##          : $quiet                 => Be quiet
##          : $verbose               => Verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;

    ## Default(s)
    my $quiet;
    my $verbose;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        quiet => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$quiet,
            strict_type => 1,
        },
        verbose => {
            default     => $arg_href->{active_parameter_href}{verbose},
            store       => \$verbose,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object now that log_file has been set
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    $log->info( q{Installing MIP into environment: }
          . $active_parameter_href->{environment_name} );

    set_container_constants( { active_parameter_href => $active_parameter_href, } );

    ## Process input parameters to get a correct combination of programs that are to be installed
    set_programs_for_installation(
        {
            active_parameter_href => $active_parameter_href,
        }
    );

    ## Cache containers
    install_containers(
        {
            active_parameter_href => $active_parameter_href,
            container_href        => $active_parameter_href->{container},
        }
    );

    ## Copy mip scripts
    install_mip_scripts(
        {
            active_parameter_href => $active_parameter_href,
        }
    );

    ## Check installation
    check_mip_installation(
        {
            active_parameter_href => $active_parameter_href,
        }
    );

    $log->info(q{Finished installing MIP});

    return;
}

1;
