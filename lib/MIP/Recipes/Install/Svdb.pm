package MIP::Recipes::Install::Svdb;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use Cwd;
use File::Spec::Functions qw{ catdir catfile };

## Cpanm
use Readonly;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_svdb };
}

## Constants
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub install_svdb {

## Function : Install SVDB
## Returns  : ""
## Arguments: $program_parameters_href => Hash with SVDB specific parameters {REF}
##          : $conda_prefix_path       => Conda prefix path
##          : $conda_environment       => Conda environment
##          : $noupdate                => Do not update
##          : $quiet                   => Be quiet
##          : $verbose                 => Set verbosity
##          : $FILEHANDLE              => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $svdb_parameters_href;
    my $conda_prefix_path;
    my $conda_environment;
    my $noupdate;
    my $quiet;
    my $verbose;
    my $FILEHANDLE;

    my $tmpl = {
        program_parameters_href => {
            required    => 1,
            default     => {},
            strict_type => 1,
            store       => \$svdb_parameters_href
        },
        conda_prefix_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$conda_prefix_path
        },
        conda_environment => {
            strict_type => 1,
            store       => \$conda_environment
        },
        noupdate => {
            strict_type => 1,
            store       => \$noupdate
        },
        quiet => {
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$quiet
        },
        verbose => {
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$verbose
        },
        FILEHANDLE => {
            required => 1,
            defined  => 1,
            store    => \$FILEHANDLE
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Modules
    use MIP::Gnu::Bash qw{ gnu_cd };
    use MIP::Gnu::Coreutils qw{ gnu_rm };
    use MIP::Log::MIP_log4perl qw{ retrieve_log };
    use MIP::Package_manager::Conda
      qw{ conda_source_activate conda_source_deactivate };
    use MIP::Package_manager::Pip qw{ check_pip_package pip_install };
    use MIP::Program::Download::Wget qw{ wget };
    use MIP::Program::Compression::Zip qw{ unzip };

    ## Unpack parameters
    my $svdb_version = $svdb_parameters_href->{version};

    ## Retrieve logger object
    my $log = retrieve_log(
        {
            log_name => q{mip_install::install_svdb},
            quiet    => $quiet,
            verbose  => $verbose,
        }
    );

    ## Store original working directory
    my $pwd = cwd();

    say {$FILEHANDLE} q{### Install SVDB};

    ## Check if svdb has been installed via pip
    my $svdb_status = check_pip_package(
        {
            package           => q{svdb},
            version           => $svdb_version,
            conda_environment => $conda_environment,
            conda_prefix_path => $conda_prefix_path,
        }
    );

    # Check if installation exists and is executable
    if ( -x catfile( $conda_prefix_path, qw{ bin SVDB } ) || $svdb_status ) {
        $log->info(
            q{SVDB is already installed in the specified conda environment.});
        if ($noupdate) {
            say {$FILEHANDLE}
              q{## SVDB is already installed, skippping reinstallation};
            say {$FILEHANDLE} $NEWLINE;
            return;
        }
        $log->warn(q{This will overwrite the current SVDB installation.});

    }

    $log->info(q{Writing instructions for SVDB installation via SHELL});

    ## Only activate conda environment if supplied by user
    if ($conda_environment) {
        ## Activate conda environment
        say $FILEHANDLE q{## Activate conda environment};
        conda_source_activate(
            {
                FILEHANDLE => $FILEHANDLE,
                env_name   => $conda_environment,
            }
        );
        say $FILEHANDLE $NEWLINE;
    }

    ## Move to miniconda environment
    say {$FILEHANDLE} q{## Move to conda environment};
    gnu_cd(
        {
            directory_path => $conda_prefix_path,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Download
    say {$FILEHANDLE} q{## Download SVDB};
    my $url =
        q{https://github.com/J35P312/SVDB/archive/SVDB-}
      . $svdb_version
      . $DOT . q{zip};
    my $svdb_zip_path = catfile( q{SVDB-} . $svdb_version . $DOT . q{zip} );
    wget(
        {
            url          => $url,
            FILEHANDLE   => $FILEHANDLE,
            quiet        => $quiet,
            verbose      => $verbose,
            outfile_path => $svdb_zip_path,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Extract
    say {$FILEHANDLE} q{## Extract};
    unzip(
        {
            infile_path => $svdb_zip_path,
            quiet       => $quiet,
            verbose     => $verbose,
            force       => 1,
            FILEHANDLE  => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Move to SVDB directory
    say {$FILEHANDLE} q{## Move to SVDB directory};
    gnu_cd(
        {
            directory_path => q{SVDB-SVDB-} . $svdb_version,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Pip install the downloaded SVDB package
    say {$FILEHANDLE} q{## Install};
    pip_install(
        {
            editable   => $DOT,
            quiet      => $quiet,
            verbose    => $verbose,
            FILEHANDLE => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Go back to subroutine origin
    say {$FILEHANDLE} q{## Moving back to original working directory};
    gnu_cd(
        {
            directory_path => $pwd,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Clean-up
    say {$FILEHANDLE} q{## Removing downloaded zip file};
    gnu_rm(
        {
            infile_path => catfile( $conda_prefix_path, $svdb_zip_path ),
            force       => 1,
            FILEHANDLE  => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Deactivate conda environment if environment exists
    if ($conda_environment) {
        say $FILEHANDLE q{## Deactivate conda environment};
        conda_source_deactivate(
            {
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say $FILEHANDLE $NEWLINE;
    }

    print {$FILEHANDLE} $NEWLINE;

    return;
}

1;
