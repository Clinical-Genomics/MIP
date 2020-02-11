package MIP::Recipes::Install::Conda;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings qw{ FATAL utf8 };
use warnings;

## CPAN
use Array::Utils qw{ intersect };
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants
  qw{ $BACKTICK $COLON $DOT $EQUALS $LOG_NAME $NEWLINE $PIPE $SPACE $UNDERSCORE };
use MIP::Gnu::Bash qw{ gnu_unset };
use MIP::Gnu::Coreutils qw{ gnu_ln };
use MIP::Log::MIP_log4perl qw{ retrieve_log };
use MIP::Program::Conda qw{ conda_create conda_install };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.21;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_conda_packages };

}

sub install_conda_packages {

## Function  : Install conda packages into a new or existing conda environment.
## Returns   :
## Arguments : $conda_env           => Name of conda environment
##           : $conda_env_path      => Path to conda environment (default: conda root)
##           : $conda_no_update_dep => Do not update dependencies
##           : $conda_packages_href => Hash holding conda packages and their version numbers {REF}
##           : $filehandle          => Filehandle to write to
##           : $quiet               => Log only warnings and above
##           : $verbose             => Log debug messages

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_env;
    my $conda_env_path;
    my $conda_packages_href;
    my $filehandle;
    my $quiet;
    my $verbose;

    ## Defaults
    my $conda_no_update_dep;

    my $tmpl = {
        conda_env => {
            required    => 1,
            store       => \$conda_env,
            strict_type => 1,
        },
        conda_env_path => {
            defined     => 1,
            required    => 1,
            store       => \$conda_env_path,
            strict_type => 1,
        },
        conda_no_update_dep => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$conda_no_update_dep,
            strict_type => 1,
        },
        conda_packages_href => {
            default     => {},
            required    => 1,
            store       => \$conda_packages_href,
            strict_type => 1,
        },
        filehandle => {
            required => 1,
            store    => \$filehandle,
        },
        quiet => {
            allow => [ undef, 0, 1 ],
            store => \$quiet,
        },
        verbose => {
            allow => [ undef, 0, 1 ],
            store => \$verbose,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments};

    ## Packages to be installed
    my @conda_packages = keys %{$conda_packages_href};

    ## Retrieve logger object
    my $log = retrieve_log(
        {
            log_name => $LOG_NAME,
            quiet    => $quiet,
            verbose  => $verbose,
        }
    );

    ## Create an array for conda packages that are to be installed from provided hash
    my @packages = _create_package_array(
        {
            package_href              => $conda_packages_href,
            package_version_separator => $EQUALS,
        }
    );

    if ( not -d $conda_env_path ) {

        ## Create conda environment and install packages
        $log->info( q{Writing installation instructions for environment: } . $conda_env );
        say {$filehandle} q{## Creating conda environment: }
          . $conda_env
          . q{ and install packages};
        conda_create(
            {
                conda_channels_ref => [qw{ bioconda conda-forge }],
                env_name           => $conda_env,
                filehandle         => $filehandle,
                packages_ref       => \@packages,
            }
        );
        say {$filehandle} $NEWLINE;
    }
    else {

        $log->warn( q{Conda environment: } . $conda_env . $SPACE . q{already exists} );
        $log->warn(q{Will try to install packages into existing environment});
        $log->info( q{Writing installation instructions for packages to environment: }
              . $conda_env );
        say {$filehandle} q{## Installing packages into existing environment};

        ## Return if no conda packages are to be installed
        return if ( scalar @packages == 0 );

        conda_install(
            {
                conda_channels_ref => [qw{ bioconda conda-forge }],
                no_update_dep      => $conda_no_update_dep,
                filehandle         => $filehandle,
                env_name           => $conda_env,
                packages_ref       => \@packages,
            }
        );
        say {$filehandle} $NEWLINE;
    }

    return;
}

sub _create_package_array {

## Function  : Takes a reference to hash of packages and creates an array with
##           : package and version joined with a supplied separator if value is defined.
##           : Also checks that the version number makes sense
## Returns   : "@packages"
## Arguments : $package_href              => Hash with packages {Hash}
##           : $package_version_separator => Scalar separating the package and the version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $package_href;
    my $package_version_separator;

    my $tmpl = {
        package_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$package_href,
            strict_type => 1,
        },
        package_version_separator => {
            defined     => 1,
            required    => 1,
            store       => \$package_version_separator,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @packages;

  PACKAGES:
    while ( my ( $package, $package_version ) = each %{$package_href} ) {

        if ( defined $package_version ) {

            push @packages, $package . $package_version_separator . $package_version;
        }
        else {
            push @packages, $package;
        }
    }
    return @packages;
}

1;
