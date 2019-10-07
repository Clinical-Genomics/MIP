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
  use MIP::Constants qw{ $COLON $LOG_NAME $NEWLINE $SPACE };
use MIP::Gnu::Coreutils qw{ gnu_mkdir };
use MIP::Program::Cadd qw{ cadd_install };
use MIP::Program::Singularity qw{ singularity_exec };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_cadd };
}

sub install_cadd {

## Function : Install/download references and annotations for CADD
## Returns  :
## Arguments: $active_parameter_href => Active parameter hash {REF}
##          : $contaienr_href        => Container hah {REF}
##          : $container_path        => Path to VEP container
##          : $FILEHANDLE            => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $container_path;
    my $container_href;
    my $FILEHANDLE;

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
        FILEHANDLE => {
            defined => 1,
            store   => \$FILEHANDLE,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    my $reference_dir_path = $active_parameter_href->{reference_dir};

    if ( not $reference_dir_path ) {
        $log->fatal(q{Please supply a reference directory when installing CADD});
        exit 1;
    }

    my $annotation_dir_path =
      catdir( $reference_dir_path, qw{ CADD-scripts data annotations } );

    ## Bind annotation dir path to pathh in container
    my $cadd_bind_path =
      $annotation_dir_path . $COLON . catdir(qw{ / opt CADD-scripts data annotations });

    ## Store annotation dir path for later
    if ( $container_href->{program_bind_paths} ) {
        push @{ $container_href->{program_bind_paths} }, $cadd_bind_path;
    }
    else {
        $container_href->{program_bind_paths} = [$cadd_bind_path];
    }
    return;
}

1;
