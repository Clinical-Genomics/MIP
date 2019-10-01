package MIP::Recipes::Install::Cadd;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
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
##          : $conda_env             => Conda environment
##          : $conda_env_path        => Conda environment path
##          : $contaienr_href        => Container hah {REF}
##          : $container_path        => Path to VEP container
##          : $FILEHANDLE            => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $conda_env;
    my $conda_env_path;
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
        conda_env => {
            store       => \$conda_env,
            strict_type => 1,
        },
        conda_env_path => {
            defined     => 1,
            required    => 1,
            store       => \$conda_env_path,
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
            defined  => 1,
            required => 1,
            store    => \$FILEHANDLE,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    my @versions            = @{ $active_parameter_href->{cadd_versions} };
    my $annotation_dir_path = $active_parameter_href->{cadd_annotation_dir};

    ## Return if only API installation
    return if ( not $active_parameter_href->{cadd_download_annotations} );

    ## Download annotations
    say {$FILEHANDLE} q{## Download CADD annotations };
    $log->info(qq{Writing instructions for CADD annotation download});

    if ( not $annotation_dir_path ) {
        $annotation_dir_path =
          catdir( $conda_env_path, qw{ CADD-scripts data annotations } );
    }

    ## Store annotation dir path for later
    if ( not $container_href->{program_bind_paths} ) {
        $container_href->{program_bind_paths} = ();
    }
    my $cadd_bind_path =
      $annotation_dir_path . $COLON . catdir(qw{ opt CADD-scripts data annotations });
    push @{ $container_href->{program_bind_paths} }, $cadd_bind_path;

    ## Make sure that the cache directory exists
    if ( not -d $annotation_dir_path ) {
        say {$FILEHANDLE} q{## Create annotation directory};
        gnu_mkdir(
            {
                FILEHANDLE       => $FILEHANDLE,
                indirectory_path => $annotation_dir_path,
                parents          => 1,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    my @cadd_install_cmds = cadd_install(
        {
            annotation_dir_path  => $annotation_dir_path,
            download_annotations => 1,
            download_prescored_indel =>
              $active_parameter_href->{cadd_download_prescored_indel},
            download_prescored_snv =>
              $active_parameter_href->{cadd_download_prescored_snv},
            download_prescored_with_annotations =>
              $active_parameter_href->{cadd_download_prescored_with_annotations},
            versions_ref => $active_parameter_href->{cadd_versions},
        }
    );
    singularity_exec(
        {
            bind_paths_ref                 => [$annotation_dir_path],
            FILEHANDLE                     => $FILEHANDLE,
            singularity_container          => $container_path,
            singularity_container_cmds_ref => \@cadd_install_cmds,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    return;
}

1;
