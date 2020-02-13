package MIP::Recipes::Install::Pip;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings qw{ FATAL utf8 };
use warnings;

## CPAN
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME $NEWLINE $SPACE };
use MIP::Log::MIP_log4perl qw{ retrieve_log };
use MIP::Package_manager::Pip qw{ pip_install };
use MIP::Program::Conda qw{ conda_activate conda_deactivate };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.07;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_pip_packages };

}

sub install_pip_packages {

## Function : Recipe for install pip package(s).
## Returns  :
## Arguments: $conda_env         => Name of conda environment
##          : $filehandle        => Filehandle to write to
##          : $pip_packages_href => Hash with pip packages and their version numbers {REF}
##          : $quiet             => Optionally turn on quiet output
##          : $verbose           => Log debug messages

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_env;
    my $filehandle;
    my $pip_packages_href;
    my $quiet;
    my $verbose;

    my $tmpl = {
        conda_env => {
            defined     => 1,
            required    => 1,
            store       => \$conda_env,
            strict_type => 1,
        },
        filehandle => {
            required => 1,
            store    => \$filehandle,
        },
        pip_packages_href => {
            default => {},
            store   => \$pip_packages_href,
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

    ## Return if no packages are to be installed
    return if not keys %{$pip_packages_href};

    ## Retrieve logger object
    my $log = retrieve_log(
        {
            log_name => $LOG_NAME,
            quiet    => $quiet,
            verbose  => $verbose,
        }
    );

    $log->info(q{Writing install instructions for pip packages});
    say {$filehandle} q{### Install PIP packages};

    ## Install PIP packages in conda environment
    say {$filehandle} q{## Install PIP packages in conda environment:} . $SPACE
      . $conda_env;

    ## Activate conda environment
    conda_activate(
        {
            env_name   => $conda_env,
            filehandle => $filehandle,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Create an array for pip packages that are to be installed from provided hash
    my @packages = _create_package_array(
        {
            package_href              => $pip_packages_href,
            package_version_separator => q{==},
        }
    );

    ## Install PIP packages
    pip_install(
        {
            filehandle   => $filehandle,
            packages_ref => \@packages,
            quiet        => $quiet,
        }
    );
    say {$filehandle} $NEWLINE;

    say {$filehandle} q{## Deactivate conda environment};
    conda_deactivate(
        {
            filehandle => $filehandle,
        }
    );
    say {$filehandle} $NEWLINE;

    return;
}

sub _create_package_array {

##Function : Takes a reference to hash of packages and creates an array with
##         : package and version joined with a supplied separator if value is defined.
##         : Also checks that the version number makes sense
##Returns  : "@packages"
##Arguments: $package_href              => Hash with packages {REF}
##         : $package_version_separator => Scalar separating the package and the version

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

  PACKAGE:
    while ( my ( $package, $package_version ) = each %{$package_href} ) {

        if ( defined $package_version ) {

            # Check that the version number matches pattern
            if ( $package_version !~ qr/\d+.\d+ | \d+.\d+.\d+/xms ) {
                croak q{The version number does not match defiend pattern for }
                  . q{package: }
                  . $package
                  . $SPACE
                  . $package_version;
            }
            push @packages, $package . $package_version_separator . $package_version;
        }
        else {

            push @packages, $package;
        }
    }
    return @packages;
}

1;
