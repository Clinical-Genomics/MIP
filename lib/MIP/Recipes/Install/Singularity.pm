package MIP::Recipes::Install::Singularity;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Check::Path qw{ check_future_filesystem_for_directory };
use MIP::Constants qw{ $COLON $LOG $NEWLINE $SPACE };
use MIP::Gnu::Coreutils qw{ gnu_mkdir };
use MIP::Log::MIP_log4perl qw{ retrieve_log };
use MIP::Package_manager::Singularity qw{ singularity_pull };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_singularity_containers };
}

sub install_singularity_containers {

## Function : Pull container from singularity hub or docker hub
## Returns  :
## Arguments: $conda_env          => Conda environemnt name
##          : $conda_env_path     => Path to conda environment
##          : $container_dir_path => Pull containers to this path
##          : $container_href     => Hash with container {REF}
##          : $FILEHANDLE         => Filehandle
##          : $quiet              => Optionally turn on quiet output
##          : $verbose            => Log debug messages

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_env;
    my $conda_env_path;
    my $container_dir_path;
    my $container_href;
    my $FILEHANDLE;
    my $quiet;
    my $verbose;

    my $tmpl = {
        conda_env => {
            store       => \$conda_env,
            strict_type => 1,
        },
        conda_env_path => {
            store       => \$conda_env_path,
            strict_type => 1,
        },
        container_dir_path => {
            store       => \$container_dir_path,
            strict_type => 1,
        },
        container_href => {
            default     => {},
            store       => \$container_href,
            strict_type => 1,
        },
        FILEHANDLE => {
            required => 1,
            store    => \$FILEHANDLE,
        },
        quiet => {
            allow       => [ undef, 0, 1 ],
            store       => \$quiet,
            strict_type => 1,
        },
        verbose => {
            allow       => [ undef, 0, 1 ],
            store       => \$verbose,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = retrieve_log(
        {
            log_name => $LOG,
            quiet    => $quiet,
            verbose  => $verbose,
        }
    );

    ## Return if no containers
    return if not keys %{$container_href};

    say {$FILEHANDLE} q{## Pull containers with Singularity};

    ## Set default path for containers
    if ( not $container_dir_path ) {
        $container_dir_path = catdir( $conda_env_path, qw{ share containers } );
    }

    ## Write check command to FILEHANDLE
    say {$FILEHANDLE} q{## Check for container path};
    my $dir_check = check_future_filesystem_for_directory(
        {
            directory_path => $container_dir_path,
        }
    );
    say {$FILEHANDLE} $dir_check . $NEWLINE;

  CONTAINER:
    foreach my $container ( keys %{$container_href} ) {

        $log->info( q{Writing instructions for pulling container}
              . $COLON
              . $SPACE
              . $container );

        say {$FILEHANDLE} q{## Setting up } . $container . q{ container};

        my $container_path = catfile( $container_dir_path, $container . q{.sif} );
        singularity_pull(
            {
                container_uri => $container_href->{$container}{uri},
                FILEHANDLE    => $FILEHANDLE,
                force         => 1,
                outfile_path  => $container_path,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }
    return 1;
}

1;
