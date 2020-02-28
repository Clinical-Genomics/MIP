package MIP::Recipes::Install::Cadd;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd q{abs_path};
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use List::MoreUtils qw{ any };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings qw{ FATAL utf8 };
use warnings;

## CPAN
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants
  qw{ $BACKWARD_SLASH $COLON $DOUBLE_QUOTE $FORWARD_SLASH $LOG_NAME $NEWLINE $SPACE };
use MIP::Set::Parameter qw{ set_container_bind_paths };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.05;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_cadd };
}

sub install_cadd {

## Function : Install/download references and annotations for CADD
## Returns  :
## Arguments: $active_parameter_href => Active parameter hash {REF}
##          : $contaienr_href        => Container hah {REF}
##          : $container_path        => Path to container
##          : $filehandle            => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $container_path;
    my $container_href;
    my $filehandle;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        container_href => {
            default     => {},
            required    => 1,
            store       => \$container_href,
            strict_type => 1,
        },
        container_path => {
            defined     => 1,
            required    => 1,
            store       => \$container_path,
            strict_type => 1,
        },
        filehandle => {
            defined => 1,
            store   => \$filehandle,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    my $reference_dir_path = $active_parameter_href->{reference_dir};

    if ( not $reference_dir_path ) {
        $log->warn(
q{Please supply a reference directory when installing CADD to use a static path}
        );
        $reference_dir_path =
            $BACKWARD_SLASH
          . $DOUBLE_QUOTE
          . $BACKWARD_SLASH
          . q{$MIP_BIND}
          . $BACKWARD_SLASH
          . $DOUBLE_QUOTE;
    }

    my $annotation_dir_path =
      catdir( $reference_dir_path, qw{ CADD-scripts data annotations } );

    ## Bind annotation dir path to pathh in container
    my $cadd_bind_path =
        $annotation_dir_path
      . $COLON
      . catdir( $FORWARD_SLASH, qw{ opt CADD-scripts data annotations } );

    ## Store annotation dir path for later
    set_container_bind_paths(
        {
            bind_paths_ref => [$cadd_bind_path],
            container_href => $container_href,
        }
    );
    return 1;
}

1;
