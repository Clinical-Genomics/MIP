package MIP::Recipes::Install::Conda;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };

#use IO::Handle;
#use File::Basename qw(dirname basename fileparse);
#use File::Spec::Functions qw(catfile catdir devnull);
use Readonly;

## MIPs lib/
use Program::Download::Wget qw(wget);
use MIP::Gnu::Bash qw(gnu_cd);
use MIP::Gnu::Coreutils qw(gnu_cp gnu_rm gnu_mv gnu_mkdir gnu_ln gnu_chmod );
use MIP::PacketManager::Conda
  qw{ conda_create conda_source_activate conda_source_deactivate conda_update conda_check conda_install };

use MIP::Check::Path qw{ check_dir_path_exist };


## Constants
Readonly my $SPACE   => q{ };
Readonly my $NEWLINE => qq{\n};


BEGIN {

    require Exporter;
    use base qw{Exporter};

    # Set the version for version checking
    our $VERSION = 1.0.0;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ setup_conda_env  install_bioconda_packages };

}


## Setup Conda env and install default packages

sub setup_conda_env {

## setup_conda_env
    
## Function  : Creates necessary conda environment and install package(s) from the default channel.
## Returns   :
## Arguments : $conda_packages_href, $conda_env, $conda_env_path, $FILEHANDLE, $conda_update
##           : $conda_packages_href => Ref to hash holding conda packages and their version numbers
##           : $conda_env           => Name of conda environment 
##           : $conda_env_path      => Path to conda environment (could bee path to root env)
##           : $FILEHANDLE          => Filehandle to wrire to
##           : $conda_update        => Update Conda if defined

    my ($arg_href) = @_;

    ## Defaults 
    
    ## Flatten argument(s)
    my $conda_packages_href;
    my $conda_env;          
    my $conda_env_path;     
    my $FILEHANDLE;
    my $conda_update;

    my $tmpl = {
        conda_packages_href => {
            required => 1,
            defined => 1,
            default => {},
            strict_type => 1,
            store => \$conda_packages_href,
        },
        conda_env => {
            strict_type => 1,
            store => \$conda_env,
        },
        conda_env_path => {
            required => 1,
            defined => 1,
            strict_type => 1,
            store => \$conda_env_path,
        },
        FILEHANDLE => {
            required => 1,
            store => \$FILEHANDLE,
        },
        conda_update => {
            store => \$conda_update,
        }
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Unix qw{ check_binary_in_path };
    use MIP::PacketManager::Conda qw{ conda_create conda_update conda_check conda_install };
    ## Scan the PATH for conda
    check_binary_in_path(
        {
            binary => q{conda},
        }
    );

    ## Check for active conda environment
    conda_check();
    
    ## Optionally update conda
    if ( $conda_update ) {
        say $FILEHANDLE q{## Updating Conda};
        conda_update(
            {
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }
    

    ## Setup the conda environment
    # Create arrray of packages that are to be installed
    my @packages = _create_package_array(
        {
            package_href              => $conda_packages_href,
            package_version_separator => q{=},
        }
    );

    if ( defined $conda_env ) {
        ## Check for existing conda environment
        if ( !-d $conda_env_path ) {
            ## Create conda environment and install packages
            say $FILEHANDLE q{## Creating conda environment: }
            . $conda_env . q{ and install packages};
            conda_create(
                {
                    env_name     => $conda_env,
                    packages_ref => \@packages,
                    FILEHANDLE   => $FILEHANDLE,
                }
            );
            say $FILEHANDLE $NEWLINE;
        }
        else {
            say STDERR q{Conda environment: } . $conda_env 
            . q{ already exists};
            say STDERR q{Will try to install packages into existing envronment};
            say {$FILEHANDLE} 
            q{## Installing conda packages into existing environment};
            conda_install(
                {
                    packages_ref => \@packages,
                    FILEHANDLE   => $FILEHANDLE,
                    env_name     => $conda_env,
                }
            );
            say {$FILEHANDLE} $NEWLINE;
        }
    }
    else {
        say {$FILEHANDLE}
        q{## Installing and/or updating python and packages in conda root};
        conda_install(
            {
                packages_ref => \@packages,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

}


sub install_bioconda_packages {


}

sub _create_package_array {

## create_package_list

##Function  : Takes a reference to hash of packages and creates an array with
##          : package and version joined with a supplied separator if value is defined.
##          : Also checks that the version number makes sense
##Returns   : "@packages"
##Arguments : $package_href, $package_version_separator
##          : $package_href              => Ref to hash with packages
##          : $package_version_separator => Scalar separating the package and the version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $package_href;
    my $package_version_separator;

    my $tmpl = {
        package_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$package_href
        },
        package_version_separator => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$package_version_separator
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @packages;

    PACKAGES:
    foreach ( keys %{$package_href} ) {
        if ( defined $package_href->{$_} ) {

            # Check that the version number matches pattern
            if ( $package_href->{$_} !~ qr/\d+.\d+ | \d+.\d+.\d+/xms ) {
                croak q{The version number does not match defiend pattern for }
                  . q{package: }
                  . $_;
            }
            push @packages, $_ . $package_version_separator
              . $package_href->{$_};
        }
        else {
            push @packages, $_;
        }
    }
    return @packages;
}


1;

## Install bioconda packages

#
#
## 
